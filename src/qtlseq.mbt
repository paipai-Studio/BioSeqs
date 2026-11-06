///|
/// qtlseq.mbt — QTL-seq / MutMap 风格 BSA 混池分离分析
/// (本地确定性实现).
///
/// 语义(QTL-seq / Takagi et al. 2013):
/// - 双混池(突变体池 mutant bulk + 野生型池 wild-type bulk)的
///   VCF 逐双等位 SNP 位点计算 SNP-index = alt AD / (ref AD + alt AD);
/// - ΔSNP-index = mutant index - wildtype index;
/// - 滑动窗口平均(窗口 bp + 步长 bp), |ΔSNP-index| 峰值区域即
///   与性状连锁的候选区间;
/// - MutMap(单池)模式: 野生型池缺省时 wt_index 恒 0,
///   delta 即 mutant SNP-index 自身;
/// - 候选区间: |mean delta| >= 阈值的连续窗口合并, BED 输出.
///
/// 输入约定: VCF FORMAT 含 AD 字段(逗号分隔 ref,alt 深度),
/// 样本名由调用方指定(mut_bulk / wt_bulk);
/// 深度 < min_depth 或 AD 缺失的位点跳过.
/// 数据不出本机.

///|
/// 分析参数.
pub(all) struct QtlseqParams {
  /// 单池最低深度.
  min_depth : Int32
  /// 滑窗大小(bp).
  window : Int32
  /// 步长(bp).
  step : Int32
  /// 候选区间 |ΔSNP-index| 阈值.
  delta_threshold : Double
}

///|
/// 默认参数(QTL-seq 常用量: 深度 >= 7, 1Mb 窗 100kb 步).
pub fn qtlseq_params() -> QtlseqParams {
  { min_depth: 7, window: 1000000, step: 100000, delta_threshold: 0.3, }
}

///|
/// 一个 SNP 位点的 SNP-index.
pub(all) struct QtlseqSite {
  chrom : String
  pos : Int32
  ref_allele : String
  alt : String
  /// 突变体池 alt 频率.
  mut_index : Double
  /// 野生型池 alt 频率(单池模式为 0).
  wt_index : Double
  mut_depth : Int32
  wt_depth : Int32
  /// ΔSNP-index.
  delta : Double
}

///|
/// 一个滑窗的聚合结果.
pub(all) struct QtlseqWindow {
  chrom : String
  /// 1-based 闭区间.
  start : Int32
  end_ : Int32
  n_snps : Int32
  mean_mut_index : Double
  mean_wt_index : Double
  mean_delta : Double
}

///|
/// 候选区间.
pub(all) struct QtlseqInterval {
  chrom : String
  /// 1-based 闭区间.
  start : Int32
  end_ : Int32
  n_snps : Int32
  /// 区间内最大 |ΔSNP-index|(窗口级).
  peak_delta : Double
  /// 峰值方向(mutant 富集为 +, 野生型富集为 -).
  direction : String
}

// ---------------------------------------------------------------------------
// SNP-index 计算
// ---------------------------------------------------------------------------

///|
/// 解析样本 AD 字段(ref,alt) -> (ref_depth, alt_depth).
/// 注: vcf.mbt 的 format 与 samples 数组均为整串按逗号拆分的产物
/// ("GT:AD" 单元素 / "0/1:10"+"10"), 此处先按冒号拆出 FORMAT 子键.
fn qtlseq_parse_ad(rec : VcfRecord, sample : String) -> (Int32, Int32)? {
  let format = match rec.format {
    Some(f) => f
    None => return None
  }
  if format.length() == 0 {
    return None
  }
  // format 键列表(整串 "GT:AD" 按冒号拆分)
  let keys = format[0].split(":").to_array()
  let mut ad_idx = -1
  for i = 0; i < keys.length(); i = i + 1 {
    if keys[i].to_owned() == "AD" {
      ad_idx = i
      break
    }
  }
  if ad_idx < 0 {
    return None
  }
  let fields = match rec.samples.get(sample) {
    Some(f) => f
    None => return None
  }
  // 还原样本整串(逗号拆分无损可逆)
  let raw_sb = StringBuilder()
  for i = 0; i < fields.length(); i = i + 1 {
    if i > 0 {
      raw_sb.write_char(',')
    }
    raw_sb.write_string(fields[i])
  }
  let raw = raw_sb.to_string()
  let sub = raw.split(":").to_array()
  if ad_idx >= sub.length() {
    return None
  }
  let ad = sub[ad_idx].to_owned()
  if ad == "." || ad.length() == 0 {
    return None
  }
  let parts = ad.split(",").to_array()
  if parts.length() < 2 {
    return None
  }
  let ref_d = parse_int(parts[0].to_owned())
  let alt_d = parse_int(parts[1].to_owned())
  Some((ref_d, alt_d))
}

///|
/// 逐位点 SNP-index 计算.
///
/// 双池模式: 两池 AD 均需达到 min_depth; 单池模式(wt_bulk 为
/// 空串): 仅突变体池, wt_index 恒 0(MutMap 语义).
pub fn qtlseq_snp_index(
  records : Array[VcfRecord],
  mut_bulk : String,
  wt_bulk~ : String = "",
  params~ : QtlseqParams = qtlseq_params(),
) -> Array[QtlseqSite] {
  let single = wt_bulk.length() == 0
  let out : Array[QtlseqSite] = Array(capacity=0)
  for rec in records {
    // 双等位 SNP
    if rec.is_snp() == false || rec.alt.length() != 1 {
      continue
    }
    let mut_ad = qtlseq_parse_ad(rec, mut_bulk)
    let (mut_ref, mut_alt) = match mut_ad {
      Some(v) => v
      None => continue
    }
    let mut_depth = mut_ref + mut_alt
    if mut_depth < params.min_depth {
      continue
    }
    let mut_index = mut_alt.to_double() / mut_depth.to_double()
    let (wt_index, wt_depth) : (Double, Int32) = if single {
      (0.0, 0)
    } else {
      let wt_ad = qtlseq_parse_ad(rec, wt_bulk)
      let (wt_ref, wt_alt) = match wt_ad {
        Some(v) => v
        None => continue
      }
      let d = wt_ref + wt_alt
      if d < params.min_depth {
        continue
      }
      (wt_alt.to_double() / d.to_double(), d)
    }
    out.push(QtlseqSite::{
      chrom: rec.chrom,
      pos: rec.pos,
      ref_allele: rec.ref_seq,
      alt: rec.alt[0],
      mut_index,
      wt_index,
      mut_depth,
      wt_depth,
      delta: mut_index - wt_index,
    })
  }
  out.sort_by(fn(a, b) { if a.pos != b.pos { a.pos - b.pos } else { 0 } })
  out
}

// ---------------------------------------------------------------------------
// 滑窗聚合
// ---------------------------------------------------------------------------

///|
/// 滑动窗口聚合(逐染色体, 窗口 1-based 闭区间, 起点固定网格).
pub fn qtlseq_windows(
  sites : Array[QtlseqSite],
  params~ : QtlseqParams = qtlseq_params(),
) -> Array[QtlseqWindow] {
  let out : Array[QtlseqWindow] = Array(capacity=0)
  if sites.length() == 0 || params.window <= 0 || params.step <= 0 {
    return out
  }
  // 按染色体分桶(保持位点顺序)
  let chrom_order : Array[String] = Array(capacity=0)
  let by_chrom : Map[String, Array[QtlseqSite]] = Map([], capacity=0)
  for s in sites {
    match by_chrom.get(s.chrom) {
      Some(arr) => arr.push(s)
      None => {
        chrom_order.push(s.chrom)
        by_chrom.set(s.chrom, [s])
      }
    }
  }
  for chrom in chrom_order {
    let ss = by_chrom.get(chrom).unwrap()
    let first_pos = ss[0].pos
    let last_pos = ss[ss.length() - 1].pos
    // 窗口网格: [start, start + window - 1] 闭区间
    let mut start = first_pos - (first_pos - 1) / params.step * params.step
    if start < 1 {
      start = 1
    }
    while start <= last_pos {
      let end_ = start + params.window - 1
      // 收集窗口内位点
      let mut n = 0
      let mut sum_mut = 0.0
      let mut sum_wt = 0.0
      for s in ss {
        if s.pos >= start && s.pos <= end_ {
          n = n + 1
          sum_mut = sum_mut + s.mut_index
          sum_wt = sum_wt + s.wt_index
        }
      }
      if n > 0 {
        out.push(QtlseqWindow::{
          chrom,
          start,
          end_,
          n_snps: n,
          mean_mut_index: sum_mut / n.to_double(),
          mean_wt_index: sum_wt / n.to_double(),
          mean_delta: (sum_mut - sum_wt) / n.to_double(),
        })
      }
      start = start + params.step
    }
  }
  out
}

// ---------------------------------------------------------------------------
// 候选区间
// ---------------------------------------------------------------------------

///|
/// 候选区间调用: |mean delta| >= 阈值的连续窗口(同染色体且相邻)
/// 合并为一个区间; 区间峰值为窗口级最大 |ΔSNP-index|.
pub fn qtlseq_intervals(
  windows : Array[QtlseqWindow],
  params~ : QtlseqParams = qtlseq_params(),
) -> Array[QtlseqInterval] {
  let out : Array[QtlseqInterval] = Array(capacity=0)
  let mut cur : Array[QtlseqWindow] = Array(capacity=0)
  let flush = fn(
    group : Array[QtlseqWindow],
    out : Array[QtlseqInterval],
  ) -> Unit {
    if group.length() == 0 {
      return
    }
    let chrom = group[0].chrom
    let start = group[0].start
    let end_ = group[group.length() - 1].end_
    let mut n = 0
    let mut peak = 0.0
    let mut direction = "+"
    for w in group {
      n = n + w.n_snps
      let abs_d = if w.mean_delta < 0.0 { -w.mean_delta } else { w.mean_delta }
      if abs_d > peak {
        peak = abs_d
        direction = if w.mean_delta >= 0.0 { "+" } else { "-" }
      }
    }
    out.push(QtlseqInterval::{
      chrom,
      start,
      end_,
      n_snps: n,
      peak_delta: peak,
      direction,
    })
  }
  for w in windows {
    let abs_d = if w.mean_delta < 0.0 { -w.mean_delta } else { w.mean_delta }
    let significant = abs_d >= params.delta_threshold
    let contiguous = cur.length() > 0 &&
      cur[cur.length() - 1].chrom == w.chrom &&
      w.start <= cur[cur.length() - 1].end_ + params.step
    if significant && (cur.length() == 0 || contiguous) {
      cur.push(w)
    } else {
      flush(cur, out)
      cur = Array(capacity=0)
      if significant {
        cur.push(w)
      }
    }
  }
  flush(cur, out)
  out
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 逐窗 TSV(chrom start end n_snps mean_mut mean_wt delta).
pub fn qtlseq_windows_tsv(windows : Array[QtlseqWindow]) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "chrom\tstart\tend\tn_snps\tmean_mut_index\tmean_wt_index\tdelta\n",
  )
  for w in windows {
    sb.write_string(w.chrom)
    sb.write_char('\t')
    sb.write_string(w.start.to_string())
    sb.write_char('\t')
    sb.write_string(w.end_.to_string())
    sb.write_char('\t')
    sb.write_string(w.n_snps.to_string())
    sb.write_char('\t')
    sb.write_string(w.mean_mut_index.to_string())
    sb.write_char('\t')
    sb.write_string(w.mean_wt_index.to_string())
    sb.write_char('\t')
    sb.write_string(w.mean_delta.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 候选区间 BED(0-based 半开).
pub fn qtlseq_intervals_bed(intervals : Array[QtlseqInterval]) -> String {
  let sb = StringBuilder()
  for iv in intervals {
    sb.write_string(iv.chrom)
    sb.write_char('\t')
    sb.write_string((iv.start - 1).to_string())
    sb.write_char('\t')
    sb.write_string(iv.end_.to_string())
    sb.write_char('\t')
    sb.write_string(
      "peak=\{iv.peak_delta};dir=\{iv.direction};snps=\{iv.n_snps}",
    )
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 位点/窗口/区间 -> JSON 值.
pub fn qtlseq_value(
  sites : Array[QtlseqSite],
  windows : Array[QtlseqWindow],
  intervals : Array[QtlseqInterval],
) -> JsonValue {
  let site_arr : Array[JsonValue] = Array(capacity=sites.length())
  for s in sites {
    site_arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(s.chrom)),
          ("pos", JsonNum(s.pos.to_double())),
          ("mut_index", JsonNum(s.mut_index)),
          ("wt_index", JsonNum(s.wt_index)),
          ("delta", JsonNum(s.delta)),
        ]),
      ),
    )
  }
  let win_arr : Array[JsonValue] = Array(capacity=windows.length())
  for w in windows {
    win_arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(w.chrom)),
          ("start", JsonNum(w.start.to_double())),
          ("end", JsonNum(w.end_.to_double())),
          ("n_snps", JsonNum(w.n_snps.to_double())),
          ("mean_delta", JsonNum(w.mean_delta)),
        ]),
      ),
    )
  }
  let iv_arr : Array[JsonValue] = Array(capacity=intervals.length())
  for iv in intervals {
    iv_arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(iv.chrom)),
          ("start", JsonNum(iv.start.to_double())),
          ("end", JsonNum(iv.end_.to_double())),
          ("peak_delta", JsonNum(iv.peak_delta)),
          ("direction", JsonStr(iv.direction)),
        ]),
      ),
    )
  }
  JsonObj(
    Map([
      ("n_sites", JsonNum(sites.length().to_double())),
      ("n_windows", JsonNum(windows.length().to_double())),
      ("n_intervals", JsonNum(intervals.length().to_double())),
      ("sites", JsonArr(site_arr)),
      ("windows", JsonArr(win_arr)),
      ("intervals", JsonArr(iv_arr)),
    ]),
  )
}

///|
/// JSON 文本.
pub fn qtlseq_json(
  sites : Array[QtlseqSite],
  windows : Array[QtlseqWindow],
  intervals : Array[QtlseqInterval],
) -> String {
  json_to_string(qtlseq_value(sites, windows, intervals))
}

///|
/// 摘要行.
pub fn qtlseq_summary_line(
  sites : Array[QtlseqSite],
  windows : Array[QtlseqWindow],
  intervals : Array[QtlseqInterval],
) -> String {
  let mut n_plus = 0
  let mut n_minus = 0
  for iv in intervals {
    if iv.direction == "+" {
      n_plus = n_plus + 1
    } else {
      n_minus = n_minus + 1
    }
  }
  "qtlseq: sites=\{sites.length()}, windows=\{windows.length()}, intervals=\{intervals.length()} (mutant-enriched=\{n_plus}, wt-enriched=\{n_minus})"
}
