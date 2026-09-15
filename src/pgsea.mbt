///|
/// Parametric Gene Set Enrichment Analysis (PGSEA).
///
/// Inspired by the R/Bioconductor package PGSEA (Kim & Volsky 2005), this
/// module scores gene set activity per sample using a parametric mean of
/// z-standardized gene expression, avoiding permutation tests used by GSEA.
///
/// Reference: Kim SY, Volsky DJ (2005) BMC Bioinformatics 6:144.

///|
/// Result of PGSEA for one gene set.
pub(all) struct PgseaResult {
  gene_set : String            // gene set name
  scores : Array[Double]       // score per sample (index-aligned with sample labels)
  mean_control : Double        // mean score in control samples
  mean_treatment : Double      // mean score in treatment samples
  t_stat : Double              // Welch t statistic (treatment vs control)
  p_value : Double             // two-sided p value
}

///|
/// Number of samples in the expression data.
fn pgsea_n_samples(expr : Array[Array[Double]]) -> Int {
  if expr.length() == 0 {
    0
  } else {
    expr[0].length()
  }
}

///|
/// Compute per-gene mean and standard deviation across all samples.
fn pgsea_gene_stats(expr : Array[Array[Double]]) -> (Array[Double], Array[Double]) {
  let n_genes = expr.length()
  let n_samples = pgsea_n_samples(expr)
  let means : Array[Double] = Array::make(n_genes, 0.0)
  let sds : Array[Double] = Array::make(n_genes, 0.0)
  for g in 0..<n_genes {
    let mut sum = 0.0
    for s in 0..<n_samples {
      sum = sum + expr[g][s]
    }
    let mean = if n_samples > 0 { sum / n_samples.to_double() } else { 0.0 }
    means[g] = mean
    let mut sq = 0.0
    for s in 0..<n_samples {
      let d = expr[g][s] - mean
      sq = sq + d * d
    }
    let sd = if n_samples > 1 { (sq / (n_samples - 1).to_double()).sqrt() } else { 0.0 }
    sds[g] = sd
  }
  (means, sds)
}

///|
/// Compute the parametric gene-set score for a single gene set and sample.
///
/// score = (1/sqrt(k)) * sum_{g in set} (x_gs - mu_g) / (sd_g + eps)
/// where k is the size of the gene set (in the data), x_gs is expression of
/// gene g in sample s, mu_g / sd_g are the gene's mean / sd across samples.
pub fn pgsea_score_one(
  expr : Array[Array[Double]],
  gene_names : Array[String],
  gene_set : Array[String],
  sample_idx : Int,
  means : Array[Double],
  sds : Array[Double],
) -> Double {
  let n_genes = expr.length()
  let mut accum = 0.0
  let mut k = 0
  for set_gene in gene_set {
    // locate the gene index
    let mut found = -1
    for g in 0..<n_genes {
      if gene_names[g] == set_gene {
        found = g
        break
      }
    }
    if found < 0 {
      continue
    }
    let g = found
    let sd = sds[g]
    let sd_eps = if sd > 1.0e-12 { sd } else { 1.0 }
    accum = accum + (expr[g][sample_idx] - means[g]) / sd_eps
    k = k + 1
  }
  if k == 0 {
    0.0
  } else {
    accum / k.to_double().sqrt()
  }
}

///|
/// Compute the full gene-set-by-sample score matrix.
/// Returns matrix with dimensions (n_gene_sets x n_samples).
pub fn pgsea_score_matrix(
  expr : Array[Array[Double]],
  gene_names : Array[String],
  gene_sets : Array[(String, Array[String])],
) -> (Array[Array[Double]], Array[String]) {
  let n_sets = gene_sets.length()
  let n_samples = pgsea_n_samples(expr)
  let (means, sds) = pgsea_gene_stats(expr)
  let matrix : Array[Array[Double]] = Array::make(n_sets, [])
  let set_names : Array[String] = Array::new()
  for i in 0..<n_sets {
    let (name, gs) = gene_sets[i]
    set_names.push(name)
    let row : Array[Double] = Array::make(n_samples, 0.0)
    for s in 0..<n_samples {
      row[s] = pgsea_score_one(expr, gene_names, gs, s, means, sds)
    }
    matrix[i] = row
  }
  (matrix, set_names)
}

///|
/// Run PGSEA for multiple gene sets, comparing two sample groups.
///
/// Parameters:
///   expr: gene-by-sample expression matrix
///   gene_names: gene symbols, index-aligned with expr rows
///   gene_sets: array of (set_name, member_genes)
///   sample_labels: per-sample group label (e.g., "control" / "treatment")
///   control_label: label designating the control group
///   treatment_label: label designating the treatment group
///
/// Returns one PgseaResult per gene set, with a Welch t-test between groups.
pub fn pgsea_run(
  expr : Array[Array[Double]],
  gene_names : Array[String],
  gene_sets : Array[(String, Array[String])],
  sample_labels : Array[String],
  control_label : String,
  treatment_label : String,
) -> Array[PgseaResult] {
  let (matrix, set_names) = pgsea_score_matrix(expr, gene_names, gene_sets)
  let n_samples = pgsea_n_samples(expr)
  let results : Array[PgseaResult] = Array::new()
  for i in 0..<set_names.length() {
    let scores = matrix[i]
    // Partition samples into control / treatment
    let control_scores : Array[Double] = Array::new()
    let treatment_scores : Array[Double] = Array::new()
    for s in 0..<n_samples {
      if sample_labels[s] == control_label {
        control_scores.push(scores[s])
      } else if sample_labels[s] == treatment_label {
        treatment_scores.push(scores[s])
      }
    }
    let mc = if control_scores.length() > 0 { pgsea_mean(control_scores) } else { 0.0 }
    let mt = if treatment_scores.length() > 0 { pgsea_mean(treatment_scores) } else { 0.0 }
    let n1 = control_scores.length()
    let n2 = treatment_scores.length()
    let v1 = if n1 > 0 { pgsea_variance(control_scores) } else { 0.0 }
    let v2 = if n2 > 0 { pgsea_variance(treatment_scores) } else { 0.0 }
    let t_stat = if n1 >= 2 && n2 >= 2 {
      let se_sq = v1 / n1.to_double() + v2 / n2.to_double()
      if se_sq < 1.0e-12 {
        0.0
      } else {
        (mt - mc) / se_sq.sqrt()
      }
    } else {
      0.0
    }
    let df = if n1 >= 2 && n2 >= 2 {
      let a = v1 / n1.to_double()
      let b = v2 / n2.to_double()
      let numer = (a + b) * (a + b)
      let denom = a * a / (n1 - 1).to_double() + b * b / (n2 - 1).to_double()
      if denom < 1.0e-12 { 1.0 } else { numer / denom }
    } else {
      1.0
    }
    let p_val = pgsea_t_pvalue(t_stat, df)
    results.push(PgseaResult::{
      gene_set: set_names[i],
      scores: scores,
      mean_control: mc,
      mean_treatment: mt,
      t_stat: t_stat,
      p_value: p_val,
    })
  }
  results
}

///|
/// Mean of a numeric array.
fn pgsea_mean(vals : Array[Double]) -> Double {
  if vals.length() == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in vals {
    s = s + v
  }
  s / vals.length().to_double()
}

///|
/// Sample variance of a numeric array.
fn pgsea_variance(vals : Array[Double]) -> Double {
  let n = vals.length()
  if n < 2 {
    return 0.0
  }
  let m = pgsea_mean(vals)
  let mut s = 0.0
  for v in vals {
    let d = v - m
    s = s + d * d
  }
  s / (n - 1).to_double()
}

///|
/// Two-sided p-value from a t statistic and degrees of freedom.
/// Uses the incomplete beta function (Student's t survival).
pub fn pgsea_t_pvalue(t : Double, df : Double) -> Double {
  let tt = t.abs()
  let x = df / (df + tt * tt)
  let ib = pgsea_incomplete_beta(df / 2.0, 0.5, x)
  ib
}

///|
/// Incomplete beta function via continued fraction (Numerical Recipes betacf).
fn pgsea_incomplete_beta(a : Double, b : Double, x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  let bt = @math.exp(pgsea_lgamma(a + b) - pgsea_lgamma(a) - pgsea_lgamma(b) +
    a * @math.ln(x) + b * @math.ln(1.0 - x))
  let use_sym = if x < (a + 1.0) / (a + b + 2.0) {
    bt * pgsea_betacf(a, b, x) / a
  } else {
    1.0 - bt * pgsea_betacf(b, a, 1.0 - x) / b
  }
  use_sym
}

///|
/// Continued fraction for the incomplete beta function.
fn pgsea_betacf(a : Double, b : Double, x : Double) -> Double {
  let maxit = 200
  let eps = 3.0e-10
  let fpmin = 1.0e-30
  let qab = a + b
  let qap = a + 1.0
  let qam = a - 1.0
  let mut c = 1.0
  let mut d = 1.0 - qab * x / qap
  if d.abs() < fpmin {
    d = fpmin
  }
  d = 1.0 / d
  let mut h = d
  let mut m = 1
  while m <= maxit {
    let m2 = 2 * m
    let mut aa = m.to_double() * (b - m.to_double()) * x /
      ((qam + m2.to_double()) * (a + m2.to_double()))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    h = h * d * c
    aa = -(a + m.to_double()) * (qab + m.to_double()) * x /
      ((a + m2.to_double()) * (qap + m2.to_double()))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    let del = d * c
    h = h * del
    if (del - 1.0).abs() < eps {
      break
    }
    m = m + 1
  }
  h
}

///|
/// Log-gamma function via Lanczos approximation.
fn pgsea_lgamma(x : Double) -> Double {
  let cof : Array[Double] = [
    76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5,
  ]
  let xx = x
  let mut y = x
  let mut tmp = xx + 5.5
  tmp = tmp - (xx + 0.5) * @math.ln(tmp)
  let mut ser = 1.000000000190015
  for j in 0..<6 {
    y = y + 1.0
    ser = ser + cof[j] / y
  }
  -tmp + @math.ln(2.5066282746310005 * ser / xx)
}

///|
/// Convenience: top gene sets ranked by significance (smallest p-value).
pub fn pgsea_top(results : Array[PgseaResult], k : Int) -> Array[PgseaResult] {
  let sorted : Array[PgseaResult] = Array::new()
  for r in results {
    sorted.push(r)
  }
  sorted.sort_by(fn(a, b) { a.p_value.compare(b.p_value) })
  let out : Array[PgseaResult] = Array::new()
  let limit = if k < sorted.length() { k } else { sorted.length() }
  for i in 0..<limit {
    out.push(sorted[i])
  }
  out
}

///|
/// Export PGSEA results as TSV.
pub fn pgsea_export(results : Array[PgseaResult]) -> String {
  let sb = StringBuilder()
  sb.write_string("gene_set\tmean_control\tmean_treatment\tt_stat\tp_value\n")
  for r in results {
    sb.write_string(
      r.gene_set + "\t" + r.mean_control.to_string() + "\t" +
      r.mean_treatment.to_string() + "\t" + r.t_stat.to_string() + "\t" +
      r.p_value.to_string() + "\n"
    )
  }
  sb.to_string()
}
