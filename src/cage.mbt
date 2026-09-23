///|
/// CAGEr：CAGE（Cap Analysis of Gene Expression）数据转录起始位点（TSS）分析。
///
/// 对应 Bioconductor CAGEr：CTSS（TSS 起始位点，单碱基分辨率）聚类、
/// 归一化（simpleTPM / power-law）、差异 TSS 使用、TSS 注释（基因特征
/// 归类、伸长/宽峰形状分类）。

/// CTSS 位点（单碱基分辨率）。
pub(all) struct CageCtss {
  chr : String
  pos : Int            // 1-based 位置
  strand : String      // "+" / "-"
  count : Int          // 原始标签数
  normalized : Double  // 归一化后计数
}

/// 单个样本的 CTSS 集合。
pub(all) struct CageSample {
  sample_id : String
  ctss : Array[CageCtss]
  total_tags : Int
}

/// CAGE 数据集。
pub(all) struct CageSet {
  name : String
  samples : Array[CageSample]
  normalization : String
}

/// CTSS 聚类结果（TSS 簇）。
pub(all) struct CageCluster {
  chr : String
  strand : String
  start : Int          // 1-based 起始
  end : Int            // 1-based 结束（含）
  n_ctss : Int         // 簇内 CTSS 数
  total_tags : Int     // 簇内总标签数
  dominant_tss : Int   // 主导 TSS（最高标签）
  shape_score : Double // 形状评分（0=宽, 1=尖）
}

/// TSS 注释结果。
pub(all) struct CageAnnotation {
  tss_pos : Int
  strand : String
  gene_symbol : String
  feature : String     // promoter / proximal / distal / intergenic
  distance_to_tss : Int
  promoter_class : String  // sharp / broad
}

/// 构造 CTSS 位点。
pub fn CageCtss::new(
  chr : String,
  pos : Int,
  strand : String,
  count : Int,
) -> CageCtss {
  CageCtss::{ chr: chr, pos: pos, strand: strand, count: count, normalized: 0.0 }
}

/// 构造样本。
pub fn CageSample::new(sample_id : String, ctss : Array[CageCtss]) -> CageSample {
  let mut total = 0
  for c in ctss {
    total = total + c.count
  }
  CageSample::{ sample_id: sample_id, ctss: ctss, total_tags: total }
}

/// 创建空 CAGE 数据集。
pub fn cage_new(name : String) -> CageSet {
  CageSet::{ name: name, samples: [], normalization: "none" }
}

/// 添加样本（返回新数据集）。
pub fn cage_add_sample(cs : CageSet, sample : CageSample) -> CageSet {
  let samples = cs.samples.copy()
  samples.push(sample)
  CageSet::{ name: cs.name, samples: samples, normalization: cs.normalization }
}

/// 样本数量。
pub fn cage_n_samples(cs : CageSet) -> Int {
  cs.samples.length()
}

/// 全数据集总标签数。
pub fn cage_total_tags(cs : CageSet) -> Int {
  let mut total = 0
  for s in cs.samples {
    total = total + s.total_tags
  }
  total
}

/// simpleTPM 归一化（tags per million）。
pub fn cage_normalize_tpm(cs : CageSet) -> CageSet {
  let samples : Array[CageSample] = Array(capacity=0)
  for s in cs.samples {
    let nf = if s.total_tags > 0 {
      s.total_tags.to_double() / 1.0e6
    } else {
      1.0
    }
    let norm_ctss : Array[CageCtss] = Array(capacity=0)
    for c in s.ctss {
      let norm_val = if nf > 0.0 { c.count.to_double() / nf } else { 0.0 }
      norm_ctss.push(CageCtss::{
        chr: c.chr, pos: c.pos, strand: c.strand,
        count: c.count, normalized: norm_val,
      })
    }
    samples.push(CageSample::new(s.sample_id, norm_ctss))
  }
  CageSet::{ name: cs.name, samples: samples, normalization: "simpleTPM" }
}


/// power-law 归一化（CAGEr 风格：以占比次方缩放，压缩高表达位点）。
pub fn cage_normalize_powerlaw(cs : CageSet, alpha? : Double = 1.2) -> CageSet {
  let alpha = alpha
  let samples : Array[CageSample] = Array(capacity=0)
  for s in cs.samples {
    let norm_ctss : Array[CageCtss] = Array(capacity=0)
    for c in s.ctss {
      let frac = if s.total_tags > 0 {
        c.count.to_double() / s.total_tags.to_double()
      } else {
        0.0
      }
      let norm_val = if frac > 0.0 {
        @math.pow(frac, 1.0 / alpha)
      } else {
        0.0
      }
      norm_ctss.push(CageCtss::{
        chr: c.chr, pos: c.pos, strand: c.strand,
        count: c.count, normalized: norm_val,
      })
    }
    samples.push(CageSample::new(s.sample_id, norm_ctss))
  }
  CageSet::{ name: cs.name, samples: samples, normalization: "powerlaw" }
}

/// 合并所有样本的 CTSS（按 chr:pos:strand 聚合）。返回排好序的列表。
fn cage_merge_ctss_list(
  cs : CageSet,
) -> Array[(String, Int, String, Int)] {
  let merged : Map[String, (String, Int, String, Int)] = Map([], capacity=1024)
  for s in cs.samples {
    for c in s.ctss {
      let key = c.chr + ":" + c.pos.to_string() + ":" + c.strand
      let old = if merged.contains(key) {
        let (_, _, _, cnt) = merged.get(key).unwrap_or(("", 0, "+", 0))
        cnt
      } else {
        0
      }
      merged.set(key, (c.chr, c.pos, c.strand, old + c.count))
    }
  }
  let list : Array[(String, Int, String, Int)] = Array(capacity=0)
  for key in merged.keys() {
    list.push(merged.get(key).unwrap_or(("", 0, "+", 0)))
  }
  // 按 chr/strand/pos 升序
  let n = list.length()
  for i in 0..<n {
    let mut min_idx = i
    let j_start = i + 1
    for j in j_start..<n {
      let (ci, pi, si, _) = list[j]
      let (cm, pm, sm, _) = list[min_idx]
      if ci < cm || (ci == cm && si < sm) || (ci == cm && si == sm && pi < pm) {
        min_idx = j
      }
    }
    if min_idx != i {
      let tmp = list[i]
      list[i] = list[min_idx]
      list[min_idx] = tmp
    }
  }
  list
}

/// CTSS 聚类：在同一染色体/链上，相邻 CTSS 距离 ≤ max_dist 的聚成一簇。
pub fn cage_cluster(
  cs : CageSet,
  max_dist? : Int = 20,
  min_tags? : Int = 2,
) -> Array[CageCluster] {
  let max_dist = max_dist
  let min_tags = min_tags
  let list = cage_merge_ctss_list(cs)
  let n = list.length()
  let clusters : Array[CageCluster] = Array(capacity=0)
  let mut i = 0
  while i < n {
    let (chr_i, pos_i, strand_i, cnt_i) = list[i]
    if cnt_i < min_tags {
      i = i + 1
      continue
    }
    // 开启一个新簇
    let start_pos = pos_i
    let mut end_pos = pos_i
    let mut total = cnt_i
    let mut dominant_pos = pos_i
    let mut max_count = cnt_i
    let mut n_ctss = 1
    let mut j = i + 1
    while j < n {
      let (chr_j, pos_j, strand_j, cnt_j) = list[j]
      if chr_j != chr_i || strand_j != strand_i {
        break
      }
      if pos_j - end_pos > max_dist {
        break
      }
      if cnt_j >= min_tags {
        end_pos = pos_j
        total = total + cnt_j
        n_ctss = n_ctss + 1
        if cnt_j > max_count {
          max_count = cnt_j
          dominant_pos = pos_j
        }
      }
      j = j + 1
    }
    // 形状评分：主导 TSS 标签占比（高占比 → 尖峰 sharp）
    let shape_score = if total > 0 {
      max_count.to_double() / total.to_double()
    } else {
      0.0
    }
    clusters.push(CageCluster::{
      chr: chr_i, strand: strand_i,
      start: start_pos, end: end_pos,
      n_ctss: n_ctss, total_tags: total,
      dominant_tss: dominant_pos, shape_score: shape_score,
    })
    i = j
  }
  clusters
}

/// TSS 形状分类（sharp / broad）。
pub fn cage_classify_tss_shape(cluster : CageCluster, threshold? : Double = 0.5) -> String {
  let th = threshold
  if cluster.shape_score >= th { "sharp" } else { "broad" }
}

/// 两样本间差异 TSS 使用（log2 fold change，含伪计数）。
pub fn cage_diff_tss(
  sample_a : CageSample,
  sample_b : CageSample,
  pseudocount? : Double = 1.0,
) -> Array[(String, Int, String, Double, Double)] {
  let pseudo = pseudocount
  let nf_a = if sample_a.total_tags > 0 { sample_a.total_tags.to_double() / 1.0e6 } else { 1.0 }
  let nf_b = if sample_b.total_tags > 0 { sample_b.total_tags.to_double() / 1.0e6 } else { 1.0 }
  let merged : Map[String, (String, Int, String, Int, Int)] = Map([], capacity=1024)
  for c in sample_a.ctss {
    let key = c.chr + ":" + c.pos.to_string() + ":" + c.strand
    merged.set(key, (c.chr, c.pos, c.strand, c.count, 0))
  }
  for c in sample_b.ctss {
    let key = c.chr + ":" + c.pos.to_string() + ":" + c.strand
    let old = if merged.contains(key) {
      let (_, _, _, ca, _) = merged.get(key).unwrap_or(("", 0, "+", 0, 0))
      ca
    } else {
      0
    }
    merged.set(key, (c.chr, c.pos, c.strand, old, c.count))
  }
  let results : Array[(String, Int, String, Double, Double)] = Array(capacity=0)
  for key in merged.keys() {
    let (chr, pos, strand, ca, cb) = merged.get(key).unwrap_or(("", 0, "+", 0, 0))
    let norm_a = ca.to_double() / nf_a + pseudo
    let norm_b = cb.to_double() / nf_b + pseudo
    let log2fc = @math.ln(norm_b / norm_a) / @math.ln(2.0)
    results.push((chr, pos, strand, log2fc, norm_a + norm_b))
  }
  results
}

/// 注释 TSS 与基因特征的关系（基于距离）。
pub fn cage_annotate_tss(
  clusters : Array[CageCluster],
  gene_tss : Array[(String, Int, String, String)],  // (gene, tss_pos, strand, chr)
  promoter_upstream? : Int = 1000,
  promoter_downstream? : Int = 200,
) -> Array[CageAnnotation] {
  let upstream = promoter_upstream
  let downstream = promoter_downstream
  let annotations : Array[CageAnnotation] = Array(capacity=0)
  for cl in clusters {
    var min_dist = 1000000000
    var nearest_gene = ""
    var nearest_feature = "intergenic"
    for g in gene_tss {
      let (gene, tss_pos, g_strand, g_chr) = g
      if g_chr != cl.chr || g_strand != cl.strand {
        continue
      }
      let dist = cl.dominant_tss - tss_pos
      let abs_dist = if dist < 0 { -dist } else { dist }
      if abs_dist < min_dist {
        min_dist = abs_dist
        nearest_gene = gene
        if dist >= -upstream && dist <= downstream {
          nearest_feature = "promoter"
        } else if abs_dist <= 5000 {
          nearest_feature = "proximal"
        } else {
          nearest_feature = "distal"
        }
      }
    }
    let pclass = cage_classify_tss_shape(cl)
    annotations.push(CageAnnotation::{
      tss_pos: cl.dominant_tss,
      strand: cl.strand,
      gene_symbol: nearest_gene,
      feature: nearest_feature,
      distance_to_tss: min_dist,
      promoter_class: pclass,
    })
  }
  annotations
}

/// 簇统计摘要：(n, total_tags, avg_tags, avg_ctss, avg_shape)。
pub fn cage_cluster_summary(
  clusters : Array[CageCluster],
) -> (Int, Int, Double, Double, Double) {
  let n = clusters.length()
  if n == 0 {
    return (0, 0, 0.0, 0.0, 0.0)
  }
  var total_tags = 0
  var total_ctss = 0
  var total_shape = 0.0
  for cl in clusters {
    total_tags = total_tags + cl.total_tags
    total_ctss = total_ctss + cl.n_ctss
    total_shape = total_shape + cl.shape_score
  }
  let nd = n.to_double()
  (
    n,
    total_tags,
    total_tags.to_double() / nd,
    total_ctss.to_double() / nd,
    total_shape / nd,
  )
}
