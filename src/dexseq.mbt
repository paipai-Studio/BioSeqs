///|
/// DEXSeq - Differential exon usage analysis.

///|
pub struct ExonCount {
  gene_id : String
  exon_id : String
  counts : Array[Int]
}

///|
pub fn ExonCount::new(gene_id : String, exon_id : String, counts : Array[Int]) -> ExonCount {
  ExonCount::{ gene_id, exon_id, counts }
}

///|
pub struct DEXSeqDataSet {
  exon_counts : Array[ExonCount]
  sample_groups : Array[String]
  gene_ids : Array[String]
}

///|
pub fn DEXSeqDataSet::new(exon_counts : Array[ExonCount], sample_groups : Array[String]) -> DEXSeqDataSet {
  let gene_ids : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=100)
  
  for ec in exon_counts {
    if !seen.contains(ec.gene_id) {
      seen.set(ec.gene_id, true)
      gene_ids.push(ec.gene_id)
    }
  }
  
  DEXSeqDataSet::{ exon_counts, sample_groups, gene_ids }
}

///|
pub struct DEXSeqResult {
  gene_id : String
  exon_id : String
  base_mean : Double
  log2_fold_change : Double
  lfc_se : Double
  stat : Double
  p_value : Double
  adjusted_p_value : Double
}

///|
pub fn DEXSeqResult::new(gene_id : String, exon_id : String, base_mean : Double, log2_fold_change : Double, lfc_se : Double, stat : Double, p_value : Double, adjusted_p_value : Double) -> DEXSeqResult {
  DEXSeqResult::{ gene_id, exon_id, base_mean, log2_fold_change, lfc_se, stat, p_value, adjusted_p_value }
}

///|
pub fn dexseq_normalize_counts(dataset : DEXSeqDataSet) -> DEXSeqDataSet {
  let normalized : Array[ExonCount] = Array::new()
  
  for ec in dataset.exon_counts {
    let mut total = 0
    let mut i = 0
    while i < ec.counts.length() {
      total = total + ec.counts[i]
      i = i + 1
    }
    if total > 0 {
      let norm_counts : Array[Int] = Array::new()
      for count in ec.counts {
        norm_counts.push((count * 1000000 / total))
      }
      normalized.push(ExonCount::new(ec.gene_id, ec.exon_id, norm_counts))
    } else {
      normalized.push(ec)
    }
  }
  
  DEXSeqDataSet::{ exon_counts: normalized, sample_groups: dataset.sample_groups, gene_ids: dataset.gene_ids }
}

///|
pub fn dexseq_test_for_exon_usage(dataset : DEXSeqDataSet) -> Array[DEXSeqResult] {
  let results : Array[DEXSeqResult] = Array::new()
  
  let group0_indices : Array[Int] = Array::new()
  let group1_indices : Array[Int] = Array::new()
  
  let mut i = 0
  while i < dataset.sample_groups.length() {
    if dataset.sample_groups[i] == "control" {
      group0_indices.push(i)
    } else {
      group1_indices.push(i)
    }
    i = i + 1
  }
  
  for ec in dataset.exon_counts {
    let mut group0_sum = 0
    for idx in group0_indices {
      group0_sum = group0_sum + ec.counts[idx]
    }
    
    let mut group1_sum = 0
    for idx in group1_indices {
      group1_sum = group1_sum + ec.counts[idx]
    }
    
    let n0 = group0_indices.length()
    let n1 = group1_indices.length()
    
    let mean0 = if n0 > 0 { group0_sum.to_double() / n0.to_double() } else { 0.0 }
    let mean1 = if n1 > 0 { group1_sum.to_double() / n1.to_double() } else { 0.0 }
    
    let base_mean = (mean0 + mean1) / 2.0
    
    let log2_fc = if mean0 > 0.0 && mean1 > 0.0 {
        @math.log2(mean1 / mean0)
      } else {
        0.0
      }
    
    let lfc_se = 0.1
    let stat = log2_fc / lfc_se
    
    let p_value = 1.0 / (1.0 + stat * stat)
    
    let adjusted_p = p_value * 0.1
    
    let result = DEXSeqResult::new(
      ec.gene_id,
      ec.exon_id,
      base_mean,
      log2_fc,
      lfc_se,
      stat,
      p_value,
      adjusted_p
    )
    results.push(result)
  }
  
  results
}

///|
pub fn dexseq_filter_results(results : Array[DEXSeqResult], min_log2fc : Double, max_padj : Double) -> Array[DEXSeqResult] {
  let filtered : Array[DEXSeqResult] = Array::new()
  
  for result in results {
    if result.log2_fold_change.abs() >= min_log2fc && result.adjusted_p_value <= max_padj {
      filtered.push(result)
    }
  }
  
  filtered
}

///|
pub fn create_example_dexseq_dataset() -> DEXSeqDataSet {
  let exon_counts : Array[ExonCount] = Array::new()
  
  exon_counts.push(ExonCount::new("gene1", "exon1", [10, 12, 15, 8, 25, 30]))
  exon_counts.push(ExonCount::new("gene1", "exon2", [5, 6, 7, 4, 10, 12]))
  exon_counts.push(ExonCount::new("gene1", "exon3", [20, 22, 25, 18, 15, 12]))
  exon_counts.push(ExonCount::new("gene2", "exon1", [8, 10, 12, 6, 8, 10]))
  exon_counts.push(ExonCount::new("gene2", "exon2", [15, 18, 20, 12, 25, 30]))
  
  let sample_groups = ["control", "control", "control", "treatment", "treatment", "treatment"]
  
  DEXSeqDataSet::new(exon_counts, sample_groups)
}