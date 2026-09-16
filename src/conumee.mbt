///|
/// conumee：拷贝数分析统一接口。
///
/// 对应 Bioconductor conumee：包装 DNAcopy 的 CBS 分段、提供基因注释、
/// 样本间比较、质量控制（MAD/noise）、ASCII/SVG 可视化与 TSV 导出。

/// 拷贝数探针数据。
pub(all) struct CnmProbe {
  chr : String
  pos : Int
  logratio : Double
}

/// 分段结果。
pub(all) struct CnmSegment {
  chr : String
  start : Int
  end : Int
  n_probes : Int
  mean_logratio : Double
  copy_number : Double
  genes : Array[String]
}

/// 样本拷贝数结果。
pub(all) struct CnmSample {
  sample_id : String
  probes : Array[CnmProbe]
  segments : Array[CnmSegment]
  qc : CnmQC
}

/// 质量控制指标。
pub(all) struct CnmQC {
  n_probes : Int
  n_segments : Int
  median_logratio : Double
  mad_logratio : Double
  noise_score : Double
}

/// 样本间比较结果。
pub(all) struct CnmComparison {
  sample_a : String
  sample_b : String
  shared_segments : Array[CnmSegment]
  unique_a : Array[CnmSegment]
  unique_b : Array[CnmSegment]
  correlation : Double
}

/// 构造探针。
pub fn CnmProbe::new(chr : String, pos : Int, logratio : Double) -> CnmProbe {
  CnmProbe::{ chr: chr, pos: pos, logratio: logratio }
}

/// 构造样本。
pub fn CnmSample::new(sample_id : String, probes : Array[CnmProbe]) -> CnmSample {
  CnmSample::{
    sample_id: sample_id, probes: probes, segments: [],
    qc: CnmQC::{ n_probes: 0, n_segments: 0, median_logratio: 0.0, mad_logratio: 0.0, noise_score: 0.0 },
  }
}

/// 从三个平行数组构造样本（长度不一致返回 None）。
pub fn cnm_from_arrays(
  sample_id : String,
  chr : Array[String],
  pos : Array[Int],
  logratio : Array[Double],
) -> CnmSample? {
  let n = chr.length()
  if pos.length() != n || logratio.length() != n {
    return None
  }
  let probes : Array[CnmProbe] = Array::new()
  for i in 0..<n {
    probes.push(CnmProbe::new(chr[i], pos[i], logratio[i]))
  }
  Some(CnmSample::new(sample_id, probes))
}

/// 中位数。
pub fn cnm_median(vals : Array[Double]) -> Double {
  let n = vals.length()
  if n == 0 {
    return 0.0
  }
  let sorted = vals.copy()
  for i in 0..<n {
    let mut min_idx = i
    let j_start = i + 1
    for j in j_start..<n {
      if sorted[j] < sorted[min_idx] {
        min_idx = j
      }
    }
    if min_idx != i {
      let tmp = sorted[i]
      sorted[i] = sorted[min_idx]
      sorted[min_idx] = tmp
    }
  }
  if (n % 2) == 1 {
    sorted[(n - 1) / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

/// MAD（中位数绝对偏差 × 1.4826）。
pub fn cnm_mad(vals : Array[Double]) -> Double {
  let n = vals.length()
  if n == 0 {
    return 0.0
  }
  let med = cnm_median(vals)
  let devs : Array[Double] = Array::new()
  for v in vals {
    let d = v - med
    devs.push(if d < 0.0 { -d } else { d })
  }
  cnm_median(devs) * 1.4826
}

/// CBS 分段（复用 DNAcopy）。
pub fn cnm_segment(
  sample : CnmSample,
  alpha? : Double = 0.01,
  nperm? : Int = 10000,
  min_width? : Int = 2,
  undo_sd? : Double = 3.0,
) -> CnmSample {
  if sample.probes.length() == 0 {
    return sample
  }
  // 构建 Cna
  let n = sample.probes.length()
  let chrom_arr : Array[String] = Array::new()
  let pos_arr : Array[Int] = Array::new()
  let lr_arr : Array[Double] = Array::new()
  for p in sample.probes {
    chrom_arr.push(p.chr)
    pos_arr.push(p.pos)
    lr_arr.push(p.logratio)
  }
  let cna = Cna::new(chrom_arr, pos_arr, lr_arr) catch {
    SeqError(_) => return sample
  }
  let opts = DnaCopyOptions::of(
    alpha=alpha,
    nperm=nperm,
    min_width=min_width,
    undo_splits="sdundo",
    undo_sd=undo_sd,
    seed=42,
  )
  let result = segment(cna, options=opts) catch {
    SeqError(_) => return sample
  }
  // 转换分段
  let segments : Array[CnmSegment] = Array::new()
  for seg in result.segments {
    let cn = 2.0 * @math.pow(2.0, seg.seg_mean)
    segments.push(CnmSegment::{
      chr: seg.chrom, start: seg.start, end: seg.end,
      n_probes: seg.num_mark, mean_logratio: seg.seg_mean,
      copy_number: cn, genes: [],
    })
  }
  // QC
  let lr_vals : Array[Double] = Array::new()
  for p in sample.probes {
    lr_vals.push(p.logratio)
  }
  let mad = cnm_mad(lr_vals)
  let qc = CnmQC::{
    n_probes: n,
    n_segments: segments.length(),
    median_logratio: cnm_median(lr_vals),
    mad_logratio: mad,
    noise_score: mad,
  }
  CnmSample::{ sample_id: sample.sample_id, probes: sample.probes, segments: segments, qc: qc }
}

/// 为分段添加基因注释（基于重叠）。
pub fn cnm_annotate_genes(
  sample : CnmSample,
  gene_positions : Array[(String, Int, Int, String)],  // (chr, start, end, gene)
) -> CnmSample {
  let annotated : Array[CnmSegment] = Array::new()
  for seg in sample.segments {
    let genes : Array[String] = Array::new()
    for g in gene_positions {
      let (g_chr, g_start, g_end, g_symbol) = g
      if g_chr != seg.chr {
        continue
      }
      if g_start <= seg.end && g_end >= seg.start {
        genes.push(g_symbol)
      }
    }
    annotated.push(CnmSegment::{
      chr: seg.chr, start: seg.start, end: seg.end,
      n_probes: seg.n_probes, mean_logratio: seg.mean_logratio,
      copy_number: seg.copy_number, genes: genes,
    })
  }
  CnmSample::{ sample_id: sample.sample_id, probes: sample.probes, segments: annotated, qc: sample.qc }
}

/// 筛出拷贝数异常的分段（gain / loss）。
pub fn cnm_filter_segments(
  sample : CnmSample,
  gain_threshold? : Double = 2.5,
  loss_threshold? : Double = 1.5,
) -> Array[CnmSegment] {
  let gain = gain_threshold
  let loss = loss_threshold
  let filtered : Array[CnmSegment] = Array::new()
  for seg in sample.segments {
    if seg.copy_number >= gain || seg.copy_number <= loss {
      filtered.push(seg)
    }
  }
  filtered
}

/// 分段键（chr:start:end）。
fn cnm_seg_key(seg : CnmSegment) -> String {
  seg.chr + ":" + seg.start.to_string() + ":" + seg.end.to_string()
}

/// 两样本比较。
pub fn cnm_compare_samples(a : CnmSample, b : CnmSample) -> CnmComparison {
  // 探针水平相关性
  let corr_pairs : Array[(Double, Double)] = Array::new()
  for pa in a.probes {
    for pb in b.probes {
      if pa.chr == pb.chr && pa.pos == pb.pos {
        corr_pairs.push((pa.logratio, pb.logratio))
      }
    }
  }
  let mut correlation = 0.0
  if corr_pairs.length() > 1 {
    var sum_x = 0.0
    var sum_y = 0.0
    var sum_xy = 0.0
    var sum_x2 = 0.0
    var sum_y2 = 0.0
    for p in corr_pairs {
      let (x, y) = p
      sum_x = sum_x + x
      sum_y = sum_y + y
      sum_xy = sum_xy + x * y
      sum_x2 = sum_x2 + x * x
      sum_y2 = sum_y2 + y * y
    }
    let nn = corr_pairs.length().to_double()
    let numer = sum_xy - sum_x * sum_y / nn
    let denom = ((sum_x2 - sum_x * sum_x / nn) * (sum_y2 - sum_y * sum_y / nn)).sqrt()
    correlation = if denom > 0.0 { numer / denom } else { 0.0 }
  }
  // 分段比较
  let b_keys : Map[String, Int] = Map([], capacity=b.segments.length())
  for sb in b.segments {
    b_keys.set(cnm_seg_key(sb), 1)
  }
  let shared : Array[CnmSegment] = Array::new()
  let unique_a : Array[CnmSegment] = Array::new()
  for sa in a.segments {
    if b_keys.contains(cnm_seg_key(sa)) {
      shared.push(sa)
    } else {
      unique_a.push(sa)
    }
  }
  let a_keys : Map[String, Int] = Map([], capacity=a.segments.length())
  for sa in a.segments {
    a_keys.set(cnm_seg_key(sa), 1)
  }
  let unique_b : Array[CnmSegment] = Array::new()
  for sb in b.segments {
    if not(a_keys.contains(cnm_seg_key(sb))) {
      unique_b.push(sb)
    }
  }
  CnmComparison::{
    sample_a: a.sample_id, sample_b: b.sample_id,
    shared_segments: shared, unique_a: unique_a, unique_b: unique_b,
    correlation: correlation,
  }
}

/// 分段 TSV 导出。
pub fn cnm_segments_tsv(segments : Array[CnmSegment]) -> String {
  let sb = StringBuilder()
  sb.write_string("chr\tstart\tend\tn_probes\tmean_logratio\tcopy_number\tgenes\n")
  for seg in segments {
    var gene_list = ""
    for i in 0..<seg.genes.length() {
      if i > 0 { gene_list = gene_list + "," }
      gene_list = gene_list + seg.genes[i]
    }
    sb.write_string(
      seg.chr + "\t" +
      seg.start.to_string() + "\t" +
      seg.end.to_string() + "\t" +
      seg.n_probes.to_string() + "\t" +
      seg.mean_logratio.to_string() + "\t" +
      seg.copy_number.to_string() + "\t" +
      gene_list + "\n",
    )
  }
  sb.to_string()
}

/// ASCII 拷贝数图（+ gain / . normal / - loss）。
pub fn cnm_ascii_plot(
  segments : Array[CnmSegment],
  chr : String,
  width? : Int = 60,
) -> String {
  let width = width
  let chr_segs : Array[CnmSegment] = Array::new()
  for seg in segments {
    if seg.chr == chr {
      chr_segs.push(seg)
    }
  }
  if chr_segs.length() == 0 {
    return chr + ": no segments"
  }
  var min_pos = chr_segs[0].start
  var max_pos = chr_segs[0].end
  for seg in chr_segs {
    if seg.start < min_pos { min_pos = seg.start }
    if seg.end > max_pos { max_pos = seg.end }
  }
  let range = max_pos - min_pos + 1
  let scale = width.to_double() / range.to_double()
  let sb = StringBuilder()
  sb.write_string(chr + ": ")
  for i in 0..<width {
    let pos = min_pos + (i.to_double() / scale).to_int()
    let mut cn = 2.0
    for seg in chr_segs {
      if pos >= seg.start && pos <= seg.end {
        cn = seg.copy_number
        break
      }
    }
    let ch = if cn >= 3.0 { '+' } else if cn <= 1.0 { '-' } else { '.' }
    sb.write_char(ch)
  }
  sb.write_string("  (+ gain / - loss)")
  sb.to_string()
}
