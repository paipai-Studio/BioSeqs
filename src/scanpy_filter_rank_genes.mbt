/// Scanpy-inspired: tl.filter_rank_genes_groups — filter differential expression results

///|
/// Result of rank_genes_groups filtering.
pub struct FilteredRankGenes {
  names : Array[Array[String]] // top gene names per group
  scores : Array[Array[Double]] // scores per group
  logfoldchanges : Array[Array[Double]] // logFC per group
  pvals : Array[Array[Double]] // p-values per group
  groups : Array[String] // group labels
}

///|
/// Filter rank_genes_groups results by criteria.
///
/// Parameters:
/// - `names`: gene names per group (from rank_genes_groups)
/// - `scores`: scores per group
/// - `logfoldchanges`: log fold changes per group
/// - `pvals`: p-values per group
/// - `groups`: group labels
/// - `min_logfoldchange`: minimum absolute logFC (default 1.0)
/// - `max_pval`: maximum p-value (default 0.05)
/// - `min_in_group_fraction`: minimum fraction of cells expressing gene in group
/// - `min_fold_change`: minimum fold change (not log)
pub fn scanpy_filter_rank_genes_groups(
  names : Array[Array[String]],
  scores : Array[Array[Double]],
  logfoldchanges : Array[Array[Double]],
  pvals : Array[Array[Double]],
  groups : Array[String],
  min_logfoldchange~ : Double = 1.0,
  max_pval~ : Double = 0.05,
) -> FilteredRankGenes {
  let n_groups = groups.length()
  let filt_names = Array(capacity=n_groups)
  let filt_scores = Array(capacity=n_groups)
  let filt_lfc = Array(capacity=n_groups)
  let filt_pvals = Array(capacity=n_groups)
  for g in 0..<n_groups {
    let gn = Array(capacity=0)
    let gs = Array(capacity=0)
    let gl = Array(capacity=0)
    let gp = Array(capacity=0)
    let n_genes = names[g].length()
    for j in 0..<n_genes {
      let lfc = if g < logfoldchanges.length() && j < logfoldchanges[g].length() {
        logfoldchanges[g][j]
      } else {
        0.0
      }
      let pv = if g < pvals.length() && j < pvals[g].length() {
        pvals[g][j]
      } else {
        1.0
      }
      // Apply filters
      let pass_lfc = lfc.abs() >= min_logfoldchange
      let pass_pval = pv <= max_pval
      if pass_lfc && pass_pval {
        gn.push(names[g][j])
        let sc = if g < scores.length() && j < scores[g].length() {
          scores[g][j]
        } else {
          0.0
        }
        gs.push(sc)
        gl.push(lfc)
        gp.push(pv)
      }
    }
    filt_names.push(gn)
    filt_scores.push(gs)
    filt_lfc.push(gl)
    filt_pvals.push(gp)
  }
  FilteredRankGenes::{
    names: filt_names,
    scores: filt_scores,
    logfoldchanges: filt_lfc,
    pvals: filt_pvals,
    groups,
  }
}
