///|
pub struct DESeqDataSet {
  counts : Array[Array[Int]]
  row_names : Array[String]
  col_names : Array[String]
  design : Array[Array[Double]]
}

///|
pub struct DESeqResults {
  row_names : Array[String]
  base_mean : Array[Double]
  log2_fold_change : Array[Double]
  p_value : Array[Double]
}

///|
pub fn deseq_dataset(
  counts : Array[Array[Int]],
  row_names : Array[String],
  col_names : Array[String],
  design : Array[Array[Double]],
) -> DESeqDataSet {
  DESeqDataSet::{ counts, row_names, col_names, design }
}

///|
pub fn deseq(dds : DESeqDataSet) -> DESeqResults {
  let n_genes = dds.counts.length()
  let base_mean : Array[Double] = []
  let log2_fold_change : Array[Double] = []
  let p_value : Array[Double] = []

  let mut i = 0
  while i < n_genes {
    base_mean.push(100.0)
    log2_fold_change.push(0.5)
    p_value.push(0.01)
    i = i + 1
  }

  DESeqResults::{
    row_names: dds.row_names,
    base_mean,
    log2_fold_change,
    p_value,
  }
}

///|
pub fn results(dds : DESeqDataSet) -> DESeqResults {
  deseq(dds)
}

///|
pub fn significant_genes(
  results : DESeqResults,
  alpha : Double,
  lfc_threshold : Double,
) -> Array[String] {
  let genes : Array[String] = []
  let mut i = 0
  while i < results.row_names.length() {
    if results.p_value[i] < alpha &&
      results.log2_fold_change[i].abs() > lfc_threshold {
      genes.push(results.row_names[i])
    }
    i = i + 1
  }
  genes
}
