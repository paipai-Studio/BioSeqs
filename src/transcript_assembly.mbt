///|
/// transcript_assembly.mbt — StringTie 风格 RNA-seq 转录本组装
/// (剪接图 + 贪心流路径求解, 本地确定性实现).
///
/// 数据流:
/// - SAM 读段 CIGAR N 操作提取剪接 junction(donor/acceptor + 计数);
/// - junction 位点切分参考轴成区段, 每碱基差分数组求区段平均深度;
/// - 剪接有向图: 区段为节点, junction 为边(容量 = 读段计数);
/// - 贪心流(简化 StringTie 网络流): 每轮枚举 source(入度 0)到
///   sink(出度 0)的全部路径, 取瓶颈(沿途区段剩余深度与边容量最小值)
///   最大者为一条 isoform, 沿途扣除流量; 瓶颈低于阈值停止;
/// - 剩余深度达阈值的孤立区段独立成单外显子 isoform;
/// - 输出合并外显子的转录本(coverage / read_count / FPKM)与 GTF 2.2.
///
/// 已知简化: 不做链感知(strand 由调用方指定或 "."),
/// 共享外显子的流量竞争由贪心顺序决定(先高瓶颈路径).
/// 数据不出本机.

///|
/// 一条剪接 junction.
pub(all) struct SpliceJunction {
  chrom : String
  /// 0-based 半开: 内含子 [donor, acceptor).
  donor : Int32
  acceptor : Int32
  /// 跨该 junction 的读段数.
  count : Int32
}

///|
/// 一条组装出的转录本 isoform.
pub(all) struct AssembledTranscript {
  id : String
  chrom : String
  strand : String
  /// 合并后的外显子区间(0-based 半开, 升序).
  exons : Array[(Int32, Int32)]
  /// 流求解得到的覆盖深度.
  coverage : Double
  /// 外显子总长(bp).
  effective_length : Int32
  /// 估计读段数(coverage * eff_len / avg_read_len).
  read_count : Double
  /// FPKM.
  fpkm : Double
}

///|
/// 组装汇总.
pub(all) struct AssemblyStats {
  n_reads : Int32
  n_junctions : Int32
  n_transcripts : Int32
  /// 全部 isoform FPKM 之和.
  total_fpkm : Double
}

// ---------------------------------------------------------------------------
// junction 提取
// ---------------------------------------------------------------------------

///|
/// 读段参考跨度终点(0-based 半开, CIGAR M/D/N/=/X 消耗参考).
fn ta_ref_end(rec : SamRecord) -> Int32 {
  let mut span = 0
  for elem in rec.cigar {
    let op = elem.op
    if op is Match || op is Delete || op is Skip || op is Equal || op is Diff {
      span = span + elem.length
    }
  }
  rec.pos - 1 + span
}

///|
/// 读段查询长度(M/I/S/=/X 消耗读段).
fn ta_query_len(rec : SamRecord) -> Int32 {
  let mut span = 0
  for elem in rec.cigar {
    let op = elem.op
    if op is Match ||
      op is Insert ||
      op is SoftClip ||
      op is Equal ||
      op is Diff {
      span = span + elem.length
    }
  }
  span
}

///|
/// 读段起点(0-based).
fn ta_start(rec : SamRecord) -> Int32 {
  rec.pos - 1
}

///|
/// 读段是否参与组装(mapped 主比对).
fn ta_eligible(rec : SamRecord) -> Bool {
  rec.rname != "*" &&
  rec.pos > 0 &&
  (rec.flag & 0x100) == 0 &&
  (rec.flag & 0x800) == 0
}

///|
/// 从 SAM 文本提取剪接 junction(CIGAR N 操作, 0-based 半开).
pub fn ta_junctions_from_sam(
  chrom : String,
  sam_text : String,
) -> Array[SpliceJunction] {
  let sam = parse_sam(sam_text)
  let counts : Map[String, Int32] = Map([], capacity=0)
  let ends : Map[String, (Int32, Int32)] = Map([], capacity=0)
  for rec in sam.records {
    if ta_eligible(rec) == false {
      continue
    }
    let mut ref_pos = ta_start(rec)
    for elem in rec.cigar {
      if elem.op is Skip {
        let donor = ref_pos
        let acceptor = ref_pos + elem.length
        let key = "\{donor}:\{acceptor}"
        let c = counts.get(key).unwrap_or(0)
        counts.set(key, c + 1)
        ends.set(key, (donor, acceptor))
        ref_pos = acceptor
      } else if elem.op is Match ||
        elem.op is Delete ||
        elem.op is Equal ||
        elem.op is Diff {
        ref_pos = ref_pos + elem.length
      }
    }
  }
  let out : Array[SpliceJunction] = Array(capacity=counts.length())
  for key in counts.keys() {
    let (donor, acceptor) = ends.get(key).unwrap()
    out.push(SpliceJunction::{
      chrom,
      donor,
      acceptor,
      count: counts.get(key).unwrap(),
    })
  }
  out.sort_by(fn(a, b) {
    if a.donor != b.donor {
      a.donor - b.donor
    } else {
      a.acceptor - b.acceptor
    }
  })
  out
}

// ---------------------------------------------------------------------------
// 组装
// ---------------------------------------------------------------------------

///|
/// 区段: 参考轴上被 junction 位点切分的区间.
priv struct TaSegment {
  start : Int32
  end_ : Int32
  /// 区段平均覆盖深度.
  depth : Double
}

///|
/// StringTie 风格转录本组装(贪心流路径求解).
///
/// min_flow 为路径瓶颈阈值(默认 2.0), max_transcripts 限制 isoform 数
/// (默认 10); strand 写入输出(不做链感知推断).
pub fn ta_assemble_transcripts(
  chrom : String,
  sam_text : String,
  strand~ : String = ".",
  min_flow~ : Double = 2.0,
  max_transcripts~ : Int32 = 10,
) -> Array[AssembledTranscript] {
  let sam = parse_sam(sam_text)
  let reads : Array[SamRecord] = Array(capacity=0)
  for rec in sam.records {
    if ta_eligible(rec) {
      reads.push(rec)
    }
  }
  if reads.length() == 0 {
    return []
  }
  // 1. junction 提取
  let junctions = ta_junctions_from_sam(chrom, sam_text)
  // 2. 边界点集合(junction 位点)
  let bounds : Array[Int32] = Array(capacity=0)
  for j in junctions {
    bounds.push(j.donor)
    bounds.push(j.acceptor)
  }
  bounds.sort_by(fn(a, b) { a - b })
  let uniq : Array[Int32] = Array(capacity=0)
  for b in bounds {
    if uniq.length() == 0 || uniq[uniq.length() - 1] != b {
      uniq.push(b)
    }
  }
  // 3. 覆盖范围与每碱基深度(差分)
  let mut min_start = 2147483647
  let mut max_end : Int32 = -2147483648
  let mut total_query_len = 0
  for rec in reads {
    let s = ta_start(rec)
    let e = ta_ref_end(rec)
    if s < min_start {
      min_start = s
    }
    if e > max_end {
      max_end = e
    }
    total_query_len = total_query_len + ta_query_len(rec)
  }
  let span_len : Int32 = max_end - min_start
  if span_len <= 0 {
    return []
  }
  let diff : Array[Int32] = Array::make(span_len + 1, 0)
  // 差分按 CIGAR M/=/X 块逐段累加(N/D 内含子不加深,
  // 否则跨 junction 读段会把内含子计入覆盖深度)
  for rec in reads {
    let mut ref_pos = ta_start(rec)
    for elem in rec.cigar {
      let op = elem.op
      if op is Match || op is Equal || op is Diff {
        let s = ref_pos - min_start
        let e = s + elem.length
        if e > s && s >= 0 && e <= span_len {
          diff[s] = diff[s] + 1
          diff[e] = diff[e] - 1
        }
        ref_pos = ref_pos + elem.length
      } else if op is Delete || op is Skip {
        ref_pos = ref_pos + elem.length
      }
    }
  }
  // 前缀和 -> 每碱基深度(按需积分)
  let depth_at = fn(p : Int32) -> Int32 {
    let idx = p - min_start
    if idx < 0 || idx >= span_len {
      return 0
    }
    let mut d = 0
    let mut i = 0
    while i <= idx {
      d = d + diff[i]
      i = i + 1
    }
    d
  }
  // 4. 区段(相邻边界间, 首尾裁剪到覆盖范围; 无 junction 时整体一段)
  let seg_bounds : Array[Int32] = Array(capacity=0)
  if uniq.length() == 0 {
    seg_bounds.push(min_start)
    seg_bounds.push(max_end)
  } else {
    seg_bounds.push(min_start)
    for b in uniq {
      if b > min_start && b < max_end {
        seg_bounds.push(b)
      }
    }
    seg_bounds.push(max_end)
  }
  let segments : Array[TaSegment] = Array(capacity=seg_bounds.length())
  for i = 0; i < seg_bounds.length() - 1; i = i + 1 {
    let s = seg_bounds[i]
    let e = seg_bounds[i + 1]
    if e <= s {
      continue
    }
    // 区段覆盖积分(逐碱基)
    let mut total = 0
    let mut p = s
    while p < e {
      total = total + depth_at(p)
      p = p + 1
    }
    // 区段实际覆盖范围(区间内有读段外显子块的子范围)
    let mut cov_s = -1
    let mut cov_e = -1
    let mut p2 = s
    while p2 < e {
      if depth_at(p2) > 0 {
        if cov_s < 0 {
          cov_s = p2
        }
        cov_e = p2 + 1
      }
      p2 = p2 + 1
    }
    if cov_s >= 0 {
      // 深度分母用覆盖子范围长度(边界区间内无覆盖的洞不稀释深度)
      let cov_len = cov_e - cov_s
      let avg = total.to_double() / cov_len.to_double()
      segments.push(TaSegment::{ start: cov_s, end_: cov_e, depth: avg, })
    }
  }
  let n_seg = segments.length()
  if n_seg == 0 {
    return []
  }
  // 5. 剪接图边: junction (donor, acceptor) -> 区段索引
  //    donor 所在区段的终点 == donor 的区段 i, acceptor 起点 == acceptor 的区段 j
  let seg_of_boundary = fn(b : Int32, upper : Bool) -> Int32? {
    // upper=true 找终点==b 的区段, false 找起点==b 的区段
    for i = 0; i < n_seg; i = i + 1 {
      if upper {
        if segments[i].end_ == b {
          return Some(i)
        }
      } else if segments[i].start == b {
        return Some(i)
      }
    }
    None
  }
  // 边容量: "i:j" -> 剩余容量
  let edge_cap : Map[String, Int32] = Map([], capacity=0)
  for j in junctions {
    match (seg_of_boundary(j.donor, true), seg_of_boundary(j.acceptor, false)) {
      (Some(i), Some(k)) =>
        if i != k {
          let key = "\{i}:\{k}"
          let c = edge_cap.get(key).unwrap_or(0)
          edge_cap.set(key, c + j.count)
        }
      _ => ()
    }
  }
  // 入度/出度
  let in_deg : Array[Int32] = Array::make(n_seg, 0)
  let out_deg : Array[Int32] = Array::make(n_seg, 0)
  for key in edge_cap.keys() {
    let parts = key.split(":").to_array()
    let i = bam_pp_parse_int_or(parts[0].to_owned(), -1)
    let k = bam_pp_parse_int_or(parts[1].to_owned(), -1)
    if i >= 0 && k >= 0 {
      out_deg[i] = out_deg[i] + 1
      in_deg[k] = in_deg[k] + 1
    }
  }
  // 邻接表
  let adj : Map[Int32, Array[Int32]] = Map([], capacity=0)
  for key in edge_cap.keys() {
    let parts = key.split(":").to_array()
    let i = bam_pp_parse_int_or(parts[0].to_owned(), -1)
    let k = bam_pp_parse_int_or(parts[1].to_owned(), -1)
    if i >= 0 && k >= 0 {
      match adj.get(i) {
        Some(arr) => arr.push(k)
        None => adj.set(i, [k])
      }
    }
  }
  // 6. 区段剩余深度
  let remain : Array[Double] = Array::make(n_seg, 0.0)
  for i = 0; i < n_seg; i = i + 1 {
    remain[i] = segments[i].depth
  }
  // 7. 贪心流(每轮迭代 DP 求最大瓶颈路径; 剪接图坐标单调,
  //    区段索引升序即拓扑序)
  let avg_read_len : Double = total_query_len.to_double() /
    reads.length().to_double()
  let n_reads = reads.length()
  let transcripts : Array[AssembledTranscript] = Array(capacity=0)
  let mut tid = 0
  while transcripts.length() < max_transcripts {
    // DP: bottleneck[j] = 到区段 j 的路径瓶颈最大值
    let bottleneck : Array[Double] = Array::make(n_seg, 0.0)
    let prev : Array[Int32] = Array::make(n_seg, -1)
    for i = 0; i < n_seg; i = i + 1 {
      if in_deg[i] == 0 && remain[i] > 0.0 {
        bottleneck[i] = remain[i]
      }
    }
    for i = 0; i < n_seg; i = i + 1 {
      if bottleneck[i] <= 0.0 {
        continue
      }
      // 计入本段容量
      if remain[i] < bottleneck[i] {
        bottleneck[i] = remain[i]
      }
      match adj.get(i) {
        Some(arr) =>
          for j in arr {
            let ekey = "\{i}:\{j}"
            let cap = edge_cap.get(ekey).unwrap_or(0).to_double()
            let cand = if cap < bottleneck[i] { cap } else { bottleneck[i] }
            if cand > bottleneck[j] {
              bottleneck[j] = cand
              prev[j] = i
            }
          }
        None => ()
      }
    }
    // 最优 sink(出度 0 且瓶颈最大)
    let mut best_sink = -1
    let mut best_flow = 0.0
    for j = 0; j < n_seg; j = j + 1 {
      let is_sink = match adj.get(j) {
        Some(arr) => arr.length() == 0
        None => true
      }
      if is_sink && bottleneck[j] > best_flow {
        best_flow = bottleneck[j]
        best_sink = j
      }
    }
    if best_sink < 0 || best_flow < min_flow {
      break
    }
    // 回溯路径
    let rev_path : Array[Int32] = Array(capacity=0)
    let mut cur = best_sink
    while cur >= 0 {
      rev_path.push(cur)
      cur = prev[cur]
    }
    if rev_path.length() < 2 {
      break
    }
    let path : Array[Int32] = Array(capacity=rev_path.length())
    let mut k = rev_path.length() - 1
    while k >= 0 {
      path.push(rev_path[k])
      k = k - 1
    }
    let bottleneck_val = best_flow
    // 扣除流量
    for node in path {
      remain[node] = remain[node] - bottleneck_val
    }
    for idx = 0; idx < path.length() - 1; idx = idx + 1 {
      let ekey = "\{path[idx]}:\{path[idx + 1]}"
      let cap = edge_cap.get(ekey).unwrap()
      edge_cap.set(ekey, (cap.to_double() - bottleneck_val).to_int())
    }
    // 外显子 = 路径区段(已由 junction 边分隔, 不再合并)
    let exons : Array[(Int32, Int32)] = Array(capacity=path.length())
    for node in path {
      exons.push((segments[node].start, segments[node].end_))
    }
    tid = tid + 1
    transcripts.push(
      ta_make_transcript(
        "TCONS_\{ta_pad8(tid)}",
        chrom,
        strand,
        exons,
        bottleneck_val,
        avg_read_len,
        n_reads,
      ),
    )
  }
  // 8. 剩余区段 -> 单外显子 isoform(无出边, 或出边全部不可用:
  //    边容量耗尽 / 下游区段剩余深度耗尽)
  for i = 0; i < n_seg && transcripts.length() < max_transcripts; i = i + 1 {
    if remain[i] >= min_flow {
      let drainable = match adj.get(i) {
        Some(arr) => {
          let mut any = false
          for j in arr {
            let cap = edge_cap.get("\{i}:\{j}").unwrap_or(0)
            if cap > 0 && remain[j] > 0.0 {
              any = true
              break
            }
          }
          any
        }
        None => false
      }
      if drainable == false {
        tid = tid + 1
        transcripts.push(
          ta_make_transcript(
            "TCONS_\{ta_pad8(tid)}",
            chrom,
            strand,
            [(segments[i].start, segments[i].end_)],
            remain[i],
            avg_read_len,
            n_reads,
          ),
        )
        remain[i] = 0.0
      }
    }
  }
  transcripts
}

///|
fn ta_pad8(n : Int32) -> String {
  let s = n.to_string()
  let mut out = ""
  let mut pad = 8 - s.length()
  while pad > 0 {
    out = out + "0"
    pad = pad - 1
  }
  out + s
}

///|
fn ta_make_transcript(
  id : String,
  chrom : String,
  strand : String,
  exons : Array[(Int32, Int32)],
  coverage : Double,
  avg_read_len : Double,
  n_reads : Int32,
) -> AssembledTranscript {
  let mut eff_len = 0
  for e in exons {
    eff_len = eff_len + (e.1 - e.0)
  }
  let read_count : Double = if avg_read_len > 0.0 {
    coverage * eff_len.to_double() / avg_read_len
  } else {
    0.0
  }
  let fpkm : Double = if n_reads > 0 && eff_len > 0 {
    read_count * 1.0e9 / (n_reads.to_double() * eff_len.to_double())
  } else {
    0.0
  }
  {
    id,
    chrom,
    strand,
    exons,
    coverage,
    effective_length: eff_len,
    read_count,
    fpkm,
  }
}

///|
/// 组装汇总.
pub fn ta_assembly_stats(
  chrom : String,
  sam_text : String,
  transcripts : Array[AssembledTranscript],
) -> AssemblyStats {
  let junctions = ta_junctions_from_sam(chrom, sam_text)
  let sam = parse_sam(sam_text)
  let mut n_reads = 0
  for rec in sam.records {
    if ta_eligible(rec) {
      n_reads = n_reads + 1
    }
  }
  let mut total_fpkm = 0.0
  for t in transcripts {
    total_fpkm = total_fpkm + t.fpkm
  }
  {
    n_reads,
    n_junctions: junctions.length(),
    n_transcripts: transcripts.length(),
    total_fpkm,
  }
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 转录本 -> GTF 2.2 文本(1-based 闭区间).
pub fn ta_transcripts_to_gtf(
  transcripts : Array[AssembledTranscript],
) -> String {
  let sb = StringBuilder()
  sb.write_string("#!genome-build. unknown\n")
  sb.write_string("##gff-version 2.2\n")
  for t in transcripts {
    let start = t.exons[0].0
    let end_ = t.exons[t.exons.length() - 1].1
    let cov = (t.coverage * 1000.0).to_int()
    sb.write_string(t.chrom)
    sb.write_string("\ttranscript_assembly\ttranscript\t")
    sb.write_string((start + 1).to_string())
    sb.write_char('\t')
    sb.write_string(end_.to_string())
    sb.write_char('\t')
    sb.write_string(cov.to_string())
    sb.write_char('\t')
    sb.write_string(if t.strand == "." { "+" } else { t.strand })
    sb.write_string("\t.\tgene_id \"TA_\{t.id}\"; transcript_id \"\{t.id}\";\n")
    let mut k = 1
    for e in t.exons {
      sb.write_string(t.chrom)
      sb.write_string("\ttranscript_assembly\texon\t")
      sb.write_string((e.0 + 1).to_string())
      sb.write_char('\t')
      sb.write_string(e.1.to_string())
      sb.write_string("\t.\t")
      sb.write_string(if t.strand == "." { "+" } else { t.strand })
      sb.write_string(
        "\t.\tgene_id \"TA_\{t.id}\"; transcript_id \"\{t.id}\"; exon_number \"\{k}\";\n",
      )
      k = k + 1
    }
  }
  sb.to_string()
}

///|
/// 转录本 -> JSON 值.
pub fn ta_transcripts_value(
  transcripts : Array[AssembledTranscript],
) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=transcripts.length())
  for t in transcripts {
    let exons : Array[JsonValue] = Array(capacity=t.exons.length())
    for e in t.exons {
      exons.push(JsonArr([JsonNum(e.0.to_double()), JsonNum(e.1.to_double())]))
    }
    arr.push(
      JsonObj(
        Map([
          ("id", JsonStr(t.id)),
          ("chrom", JsonStr(t.chrom)),
          ("strand", JsonStr(t.strand)),
          ("exons", JsonArr(exons)),
          ("coverage", JsonNum(t.coverage)),
          ("effective_length", JsonNum(t.effective_length.to_double())),
          ("read_count", JsonNum(t.read_count)),
          ("fpkm", JsonNum(t.fpkm)),
        ]),
      ),
    )
  }
  JsonArr(arr)
}

///|
/// 转录本 -> JSON 文本.
pub fn ta_transcripts_json(transcripts : Array[AssembledTranscript]) -> String {
  json_to_string(ta_transcripts_value(transcripts))
}

///|
/// 组装摘要行.
pub fn ta_assembly_stats_line(stats : AssemblyStats) -> String {
  "transcript-assembly: reads=\{stats.n_reads}, junctions=\{stats.n_junctions}, transcripts=\{stats.n_transcripts}, total_fpkm=\{stats.total_fpkm}"
}
