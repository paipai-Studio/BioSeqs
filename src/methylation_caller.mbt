///|
/// methylation_caller.mbt — Bismark 风格甲基化位点调用器
/// (methylation extractor, 本地确定性实现).
///
/// 语义(简化 Bismark methylation extractor):
/// - 输入 SAM 的读段为已比对到参考的转化读段(比对由外部完成,
///   如方向性文库的 3-letter 比对结果投影);
/// - 正链读段(flag 无 0x10, OT/CTOT 语义): C->T 转化 —
///   参考 C 位点上读段为 C 记 methylated, 为 T 记 unmethylated;
/// - 反链读段(flag 含 0x10, OB/CTOB 语义): G->A 转化 —
///   参考 G 位点上读段为 G 记 methylated, 为 A 记 unmethylated;
/// - 上下文按参考判定: CpG / CHG / CHH;
/// - 聚合逐 (pos, strand) 的 meth/unmeth 计数, 输出
///   Bismark TSV(与 bsseq_parse_bismark 闭环)与 BSRecord 数组
///   (直接进入 bsseq 平滑/DMR 统计链).
///
/// 已知简化: 不做 4-strand 比对选择(Bismark 由比对器决定方向),
/// 链语义由 SAM flag 推断; 不区分 OT/CTOT 的重复计数.
/// 数据不出本机.

///|
/// 单个甲基化位点的调用结果(聚合后).
pub(all) struct MethylationCall {
  chrom : String
  /// 1-based 位点坐标.
  pos : Int32
  /// OT/CpG 语境正链为 "+", OB 语境反链为 "-".
  strand : String
  /// CpG / CHG / CHH.
  context : String
  meth : Int32
  unmeth : Int32
}

///|
/// 调用汇总.
pub(all) struct MethylationCallSummary {
  n_reads : Int32
  /// 参与调用的读段数(mapped 主比对).
  n_cpg_sites : Int32
  n_chg_sites : Int32
  n_chh_sites : Int32
  /// CpG 全局甲基化率.
  cpg_methylation_rate : Double
}

// ---------------------------------------------------------------------------
// 提取
// ---------------------------------------------------------------------------

///|
fn mc_eligible(rec : SamRecord) -> Bool {
  rec.rname != "*" && rec.pos > 0 && (rec.flag & 0x100) == 0
}

///|
/// 参考位点的上下文判定.
/// 正链(C 位点): CpG = ref[p+1]=='G'; CHG = ref[p+2]=='G';
/// 反链(G 位点): CpG = ref[p-1]=='C'; CHG = ref[p-2]=='C'.
fn mc_context(reference : String, pos0 : Int32, reverse : Bool) -> String {
  let n = reference.length()
  let get = fn(i : Int32) -> Char {
    if i < 0 || i >= n {
      'N'
    } else {
      reference.unsafe_get(i).unsafe_to_char()
    }
  }
  if reverse == false {
    // C 位点: 看下游
    if get(pos0 + 1) == 'G' {
      return "CpG"
    }
    if get(pos0 + 2) == 'G' {
      return "CHG"
    }
    "CHH"
  } else {
    // G 位点: 看上游
    if get(pos0 - 1) == 'C' {
      return "CpG"
    }
    if get(pos0 - 2) == 'C' {
      return "CHG"
    }
    "CHH"
  }
}

///|
/// 从 SAM 文本提取逐位点甲基化计数.
///
/// 逐读段 CIGAR walk(M 块内参考-读段逐碱基比较),
/// 聚合键为 (pos, strand, context); 输出按 (pos, strand) 升序.
pub fn mc_extract(
  chrom : String,
  reference : String,
  sam_text : String,
) -> Array[MethylationCall] {
  let sam = parse_sam(sam_text)
  // 聚合键 "pos:strand" -> [context, meth, unmeth]
  let table : Map[String, (String, Int32, Int32)] = Map([], capacity=0)
  for rec in sam.records {
    if mc_eligible(rec) == false {
      continue
    }
    let reverse = (rec.flag & 0x10) != 0
    let strand = if reverse { "-" } else { "+" }
    let mut ref_pos = rec.pos - 1
    let read = rec.seq
    let mut qi = 0
    for elem in rec.cigar {
      let op = elem.op
      let len = elem.length
      if op is Match || op is Equal || op is Diff {
        for k = 0; k < len; k = k + 1 {
          if ref_pos >= 0 && ref_pos < reference.length() && qi < read.length() {
            let rb = reference.unsafe_get(ref_pos).unsafe_to_char()
            let qb = read.unsafe_get(qi).unsafe_to_char()
            let ctx : String? = if reverse == false {
              // OT: 参考 C, 读 C(meth)/T(unmeth)
              if rb == 'C' {
                if qb == 'C' {
                  Some("meth")
                } else if qb == 'T' {
                  Some("unmeth")
                } else {
                  None
                }
              } else {
                None
              }
              // OB: 参考 G, 读 G(meth)/A(unmeth)
            } else if rb == 'G' {
              if qb == 'G' {
                Some("meth")
              } else if qb == 'A' {
                Some("unmeth")
              } else {
                None
              }
            } else {
              None
            }
            match ctx {
              Some(kind) => {
                let context = mc_context(reference, ref_pos, reverse)
                let key = "\{ref_pos}:\{strand}"
                let (c0, m0, u0) = table.get(key).unwrap_or((context, 0, 0))
                if kind == "meth" {
                  table.set(key, (c0, m0 + 1, u0))
                } else {
                  table.set(key, (c0, m0, u0 + 1))
                }
              }
              None => ()
            }
          }
          ref_pos = ref_pos + 1
          qi = qi + 1
        }
      } else if op is Insert || op is SoftClip {
        qi = qi + len
      } else {
        // Delete / Skip / HardClip / Pad
        ref_pos = ref_pos + len
      }
    }
  }
  // 输出排序
  let keys : Array[String] = Array(capacity=table.length())
  for key in table.keys() {
    keys.push(key)
  }
  keys.sort_by(fn(a, b) { a.compare(b) })
  let out : Array[MethylationCall] = Array(capacity=keys.length())
  for key in keys {
    let (context, meth, unmeth) = table.get(key).unwrap()
    let colon = key.find(":").unwrap()
    let pos0 = parse_int(key[0:colon].to_owned())
    let strand = key[colon + 1:].to_owned()
    out.push(MethylationCall::{
      chrom,
      pos: pos0 + 1,
      strand,
      context,
      meth,
      unmeth,
    })
  }
  out
}

// ---------------------------------------------------------------------------
// 衔接 bsseq
// ---------------------------------------------------------------------------

///|
/// 调用结果 -> Bismark TSV 文本(bsseq_parse_bismark 闭环:
/// chrom pos strand meth unmeth context).
pub fn mc_bismark_tsv(calls : Array[MethylationCall]) -> String {
  let sb = StringBuilder()
  sb.write_string("Bismark methylation extractor output (mc)\n")
  for c in calls {
    sb.write_string(c.chrom)
    sb.write_char('\t')
    sb.write_string(c.pos.to_string())
    sb.write_char('\t')
    sb.write_string(c.strand)
    sb.write_char('\t')
    sb.write_string(c.meth.to_string())
    sb.write_char('\t')
    sb.write_string(c.unmeth.to_string())
    sb.write_char('\t')
    sb.write_string(c.context)
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 调用结果 -> BSRecord 数组(直接进入 bsseq 平滑/DMR 统计链).
pub fn mc_to_bsrecords(calls : Array[MethylationCall]) -> Array[BSRecord] {
  let out : Array[BSRecord] = Array(capacity=calls.length())
  for c in calls {
    out.push(
      BSRecord::new_with_context(
        c.chrom,
        c.pos,
        c.strand,
        c.meth,
        c.unmeth,
        c.context,
      ),
    )
  }
  out
}

// ---------------------------------------------------------------------------
// 汇总与导出
// ---------------------------------------------------------------------------

///|
/// 调用汇总(CpG 全局甲基化率等).
pub fn mc_summary(
  sam_text : String,
  calls : Array[MethylationCall],
) -> MethylationCallSummary {
  let sam = parse_sam(sam_text)
  let mut n_reads = 0
  for rec in sam.records {
    if mc_eligible(rec) {
      n_reads = n_reads + 1
    }
  }
  let mut n_cpg = 0
  let mut n_chg = 0
  let mut n_chh = 0
  let mut meth = 0
  let mut unmeth = 0
  for c in calls {
    match c.context {
      "CpG" => {
        n_cpg = n_cpg + 1
        meth = meth + c.meth
        unmeth = unmeth + c.unmeth
      }
      "CHG" => n_chg = n_chg + 1
      "CHH" => n_chh = n_chh + 1
      _ => ()
    }
  }
  let total = meth + unmeth
  let rate : Double = if total > 0 {
    meth.to_double() / total.to_double()
  } else {
    0.0
  }
  {
    n_reads,
    n_cpg_sites: n_cpg,
    n_chg_sites: n_chg,
    n_chh_sites: n_chh,
    cpg_methylation_rate: rate,
  }
}

///|
/// 调用结果 -> JSON 值.
pub fn mc_value(calls : Array[MethylationCall]) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=calls.length())
  for c in calls {
    arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(c.chrom)),
          ("pos", JsonNum(c.pos.to_double())),
          ("strand", JsonStr(c.strand)),
          ("context", JsonStr(c.context)),
          ("meth", JsonNum(c.meth.to_double())),
          ("unmeth", JsonNum(c.unmeth.to_double())),
        ]),
      ),
    )
  }
  JsonArr(arr)
}

///|
/// 调用结果 -> JSON 文本.
pub fn mc_json(calls : Array[MethylationCall]) -> String {
  json_to_string(mc_value(calls))
}

///|
/// 汇总摘要行.
pub fn mc_summary_line(summary : MethylationCallSummary) -> String {
  "methylation-caller: reads=\{summary.n_reads}, cpg=\{summary.n_cpg_sites}, chg=\{summary.n_chg_sites}, chh=\{summary.n_chh_sites}, cpg_rate=\{summary.cpg_methylation_rate}"
}
