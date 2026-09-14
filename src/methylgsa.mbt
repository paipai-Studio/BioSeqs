///|
/// Gene set testing for methylation data, ported from Bioconductor's
/// `methylGSA` package.
///
/// The reference package tests whether gene sets are enriched among
/// CpG-associated genes.  Its `methylRRA` back-end aggregates per-gene
/// minimal CpG p-values two ways: OLS (mean of -log10 minimal p within
/// the set, scored against the gene universe) and Robust (Robust Rank
/// Aggregation: order statistics of rank-scaled p-values under
/// Beta(k, n-k+1)).  This port implements both aggregations plus
/// BH-FDR multiple testing correction.

///|
/// Per-CpG input: identifiers, p-values and a CpG-to-gene mapping.
pub struct MethylgsaInput {
  cpg_ids : Array[String]
  cpg_pvalues : Array[Double]
  cpg_to_gene : Map[String, String]
} derive(Eq)

///|
pub fn MethylgsaInput::new(
  cpg_ids~ : Array[String],
  cpg_pvalues~ : Array[Double],
  cpg_to_gene~ : Map[String, String],
) -> MethylgsaInput {
  MethylgsaInput::{ cpg_ids, cpg_pvalues, cpg_to_gene }
}

///|
/// Result for one gene set.
pub struct MethylgsaResult {
  set_name : String
  set_size : Int
  statistic : Double
  pvalue : Double
  qvalue : Double
} derive(Eq, Show, Debug)

///|
pub fn methylgsa_result_new(
  set_name~ : String,
  set_size~ : Int,
  statistic~ : Double,
  pvalue~ : Double,
  qvalue~ : Double,
) -> MethylgsaResult {
  MethylgsaResult::{ set_name, set_size, statistic, pvalue, qvalue }
}

///|
/// Per-gene minimal CpG p-value (first step of methylRRA).
pub fn methylgsa_gene_min_p(input : MethylgsaInput) -> Map[String, Double] {
  let gene_min : Map[String, Double] = Map::new()
  for i in 0..<input.cpg_ids.length() {
    let cpg = input.cpg_ids[i]
    let p = if i < input.cpg_pvalues.length() {
      input.cpg_pvalues[i]
    } else {
      1.0
    }
    match input.cpg_to_gene.get(cpg) {
      Some(g) => {
        match gene_min.get(g) {
          Some(cur) => {
            if p < cur {
              gene_min.set(g, p)
            }
          }
          None => gene_min.set(g, p)
        }
      }
      None => ()
    }
  }
  gene_min
}

///|
/// Mean of -log10(min p) over the genes present in `genes`, or 0.0
/// when none map.
pub fn methylgsa_ols_statistic(
  gene_min : Map[String, Double],
  genes : Array[String],
) -> (Double, Int) {
  let mut sum = 0.0
  let mut count = 0
  for g in genes {
    match gene_min.get(g) {
      Some(p) => {
        if p > 0.0 {
          let logp = methylgsa_log10(p)
          sum += 0.0 - logp
        } else {
          sum += 20.0
        }
        count += 1
      }
      None => ()
    }
  }
  if count == 0 {
    (0.0, 0)
  } else {
    (sum / count.to_double(), count)
  }
}

///|
/// log10 via bisection on 10^k.
pub fn methylgsa_log10(x : Double) -> Double {
  if x <= 0.0 {
    return -20.0
  }
  let mut lo = -20.0
  let mut hi = 20.0
  for _ in 0..<60 {
    let mid = (lo + hi) / 2.0
    let pow = methylgsa_pow10(mid)
    if pow < x {
      lo = mid
    } else {
      hi = mid
    }
  }
  (lo + hi) / 2.0
}

///|
fn methylgsa_pow10(x : Double) -> Double {
  let ln10 = 2.302585092994046
  methylgsa_exp(x * ln10)
}

///|
fn methylgsa_exp(x : Double) -> Double {
  if x > 40.0 {
    return methylgsa_exp(40.0) * methylgsa_exp(x - 40.0)
  }
  if x < -40.0 {
    return 0.0
  }
  let mut sum = 1.0
  let mut term = 1.0
  for i in 1..<30 {
    term = term * x / i.to_double()
    sum += term
    let abs_term = if term < 0.0 { 0.0 - term } else { term }
    if abs_term < 1.0e-12 {
      break
    }
  }
  sum
}

///|
/// Collect distinct gene set names.
fn methylgsa_collect_names(
  gene_sets : Map[String, Array[String]],
) -> Array[String] {
  let names : Array[String] = Array::new()
  for k in gene_sets.keys() {
    names.push(k)
  }
  names
}

///|
/// methylRRA-OLS gene set test.
///
/// Statistic: mean -log10(minimal CpG p) within each set.  Empirical
/// p-value: fraction of gene-universe background statistics that meet
/// or exceed the set statistic (one-sided, greater).
pub fn methylrra_ols(
  input : MethylgsaInput,
  gene_sets~ : Map[String, Array[String]],
) -> Array[MethylgsaResult] {
  let gene_min = methylgsa_gene_min_p(input)
  let universe : Array[String] = Array::new()
  for g in gene_min.keys() {
    universe.push(g)
  }
  let background : Array[Double] = Array::new()
  for g in universe {
    let singleton : Map[String, Array[String]] = Map::new()
    singleton.set("bg", [g])
    let (stat, _) = methylgsa_ols_statistic(gene_min, singleton["bg"])
    background.push(stat)
  }
  let results : Array[MethylgsaResult] = Array::new()
  for name in methylgsa_collect_names(gene_sets) {
    match gene_sets.get(name) {
      Some(genes) => {
        let (stat, size) = methylgsa_ols_statistic(gene_min, genes)
        let mut exceed = 0
        for b in background {
          if b >= stat {
            exceed += 1
          }
        }
        let p = if background.length() == 0 {
          1.0
        } else {
          exceed.to_double() / background.length().to_double()
        }
        results.push(
          MethylgsaResult::{
            set_name: name,
            set_size: size,
            statistic: stat,
            pvalue: p,
            qvalue: 1.0,
          },
        )
      }
      None => ()
    }
  }
  methylgsa_bh(results)
}

///|
/// methylRRA-Robust gene set test (Robust Rank Aggregation).
///
/// For each gene in a set, its minimal p-value is converted to a rank
/// fraction over the universe.  Within the set, fractions are sorted
/// ascending; each order statistic k of n under uniform yields a
/// Beta(k, n-k+1) CDF p-value.  The aggregation is min_k * n
/// (the classic RRA score), capped at 1.0.
pub fn methylrra_robust(
  input : MethylgsaInput,
  gene_sets~ : Map[String, Array[String]],
) -> Array[MethylgsaResult] {
  let gene_min = methylgsa_gene_min_p(input)
  let universe : Array[String] = Array::new()
  for g in gene_min.keys() {
    universe.push(g)
  }
  let n_universe = universe.length()
  let sorted_p : Array[Double] = Array::new()
  for g in universe {
    match gene_min.get(g) {
      Some(p) => sorted_p.push(p)
      None => ()
    }
  }
  let arr = sorted_p.copy()
  arr.sort_by(fn(a, b) {
    if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let results : Array[MethylgsaResult] = Array::new()
  for name in methylgsa_collect_names(gene_sets) {
    match gene_sets.get(name) {
      Some(genes) => {
        let fractions : Array[Double] = Array::new()
        let mut size = 0
        for g in genes {
          match gene_min.get(g) {
            Some(p) => {
              let mut rank = 0
              for q in arr {
                if q < p {
                  rank += 1
                }
              }
              let frac = if n_universe <= 1 {
                0.5
              } else {
                rank.to_double() / (n_universe - 1).to_double()
              }
              fractions.push(frac)
              size += 1
            }
            None => ()
          }
        }
        if size == 0 {
          results.push(
            MethylgsaResult::{
              set_name: name,
              set_size: 0,
              statistic: 0.0,
              pvalue: 1.0,
              qvalue: 1.0,
            },
          )
        } else {
          let fr = fractions.copy()
          fr.sort_by(fn(a, b) {
            if a < b {
              -1
            } else if a > b {
              1
            } else {
              0
            }
          })
          let mut best = 1.0
          for k in 0..<size {
            let beta_p = bsseq_betainc(
              fr[k],
              (k + 1).to_double(),
              (size - k).to_double(),
            )
            let score = beta_p * size.to_double()
            if score < best {
              best = score
            }
          }
          let agg = if best > 1.0 { 1.0 } else { best }
          results.push(
            MethylgsaResult::{
              set_name: name,
              set_size: size,
              statistic: agg,
              pvalue: agg,
              qvalue: 1.0,
            },
          )
        }
      }
      None => ()
    }
  }
  methylgsa_bh(results)
}

///|
/// Benjamini-Hochberg FDR correction over the results in place
/// order, filling `qvalue`.
pub fn methylgsa_bh(results : Array[MethylgsaResult]) -> Array[MethylgsaResult] {
  let n = results.length()
  if n == 0 {
    return results
  }
  let order : Array[Int] = Array::new()
  for i in 0..<n {
    order.push(i)
  }
  order.sort_by(fn(a, b) {
    let pa = results[a].pvalue
    let pb = results[b].pvalue
    if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  let mut prev_q = 1.0
  let mut idx = n - 1
  while idx >= 0 {
    let i = order[idx]
    let rank = (idx + 1).to_double()
    let raw = results[i].pvalue * n.to_double() / rank
    let q = if raw < prev_q { raw } else { prev_q }
    let qc = if q > 1.0 { 1.0 } else { q }
    results[i] = MethylgsaResult::{
      set_name: results[i].set_name,
      set_size: results[i].set_size,
      statistic: results[i].statistic,
      pvalue: results[i].pvalue,
      qvalue: qc,
    }
    prev_q = qc
    idx -= 1
  }
  results
}

///|
/// Sample methylation input: 8 CpGs across 5 genes.
pub fn methylgsa_sample_input() -> MethylgsaInput {
  MethylgsaInput::new(
    cpg_ids=["cg1", "cg2", "cg3", "cg4", "cg5", "cg6", "cg7", "cg8"],
    cpg_pvalues=[0.001, 0.010, 0.400, 0.900, 0.002, 0.700, 0.500, 0.030],
    cpg_to_gene={
      let m : Map[String, String] = Map::new()
      m.set("cg1", "GENE_A")
      m.set("cg2", "GENE_A")
      m.set("cg3", "GENE_B")
      m.set("cg4", "GENE_C")
      m.set("cg5", "GENE_D")
      m.set("cg6", "GENE_E")
      m.set("cg7", "GENE_A")
      m.set("cg8", "GENE_D")
      m
    },
  )
}

///|
/// Sample gene sets over the sample universe.
pub fn methylgsa_sample_gene_sets() -> Map[String, Array[String]] {
  let m : Map[String, Array[String]] = Map::new()
  m.set("SIGNIFICANT_SET", ["GENE_A", "GENE_D"])
  m.set("NEUTRAL_SET", ["GENE_C", "GENE_E"])
  m
}