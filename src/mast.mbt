///|
pub struct MastData {
  expression : Array[Array[Double]]
  gene_names : Array[String]
  cell_names : Array[String]
  groups : Array[Int]
  cngeneson : Array[Double]
  n_genes : Int
  n_cells : Int
  n_groups : Int
}

pub fn MastData::new(
  expression : Array[Array[Double]],
  gene_names : Array[String],
  cell_names : Array[String],
  groups : Array[Int],
) -> MastData {
  let n_cells = expression.length()
  let n_genes = if n_cells > 0 { expression[0].length() } else { 0 }
  
  let cngeneson : Array[Double] = Array::make(n_cells, 0.0)
  let mut i = 0
  while i < n_cells {
    let mut detected = 0
    let mut j = 0
    while j < n_genes {
      if expression[i][j] > 0.0 {
        detected = detected + 1
      }
      j = j + 1
    }
    cngeneson[i] = if n_genes > 0 { detected.to_double() / n_genes.to_double() } else { 0.0 }
    i = i + 1
  }
  
  let mut max_group = 0
  let mut k = 0
  while k < groups.length() {
    if groups[k] > max_group {
      max_group = groups[k]
    }
    k = k + 1
  }
  
  MastData::{
    expression,
    gene_names,
    cell_names,
    groups,
    cngeneson,
    n_genes,
    n_cells,
    n_groups: max_group + 1,
  }
}

///|
pub struct HurdleResult {
  gene : Array[String]
  discrete_pval : Array[Double]
  continuous_pval : Array[Double]
  combined_pval : Array[Double]
  discrete_fdr : Array[Double]
  continuous_fdr : Array[Double]
  combined_fdr : Array[Double]
  log2_fc : Array[Double]
  detection_rate : Array[Double]
}

pub fn HurdleResult::new(n_genes : Int) -> HurdleResult {
  HurdleResult::{
    gene: Array::make(n_genes, ""),
    discrete_pval: Array::make(n_genes, 1.0),
    continuous_pval: Array::make(n_genes, 1.0),
    combined_pval: Array::make(n_genes, 1.0),
    discrete_fdr: Array::make(n_genes, 1.0),
    continuous_fdr: Array::make(n_genes, 1.0),
    combined_fdr: Array::make(n_genes, 1.0),
    log2_fc: Array::make(n_genes, 0.0),
    detection_rate: Array::make(n_genes, 0.0),
  }
}

pub fn HurdleResult::to_string(self : HurdleResult) -> String {
  let mut result = "HurdleResult\n"
  result = result + "=============\n"
  result = result + "Total genes: " + self.gene.length().to_string() + "\n\n"
  result = result + "Top 5 genes by combined p-value:\n"
  
  let indices : Array[Int] = Array::make(self.gene.length(), 0)
  let mut i = 0
  while i < self.gene.length() {
    indices[i] = i
    i = i + 1
  }
  
  let mut j = 0
  while j < mast_min_i(5, self.gene.length()) {
    let mut best_idx = j
    let mut best_pval = self.combined_pval[j]
    let mut k = j + 1
    while k < self.gene.length() {
      if self.combined_pval[k] < best_pval {
        best_pval = self.combined_pval[k]
        best_idx = k
      }
      k = k + 1
    }
    let temp = indices[j]
    indices[j] = indices[best_idx]
    indices[best_idx] = temp
    
    let idx = indices[j]
    result = result + "  " + self.gene[idx] + ": p=" + self.combined_pval[idx].to_string() + 
             " log2FC=" + self.log2_fc[idx].to_string() + "\n"
    j = j + 1
  }
  
  result
}

///|
fn mast_min_i(a : Int, b : Int) -> Int {
  if a < b { a } else { b }
}

///|
fn ln_mast(x : Double) -> Double {
  if x <= 0.0 {
    -1000.0
  } else {
    let mut result = 0.0
    let term = (x - 1.0) / (x + 1.0)
    let mut denom = 1.0
    let mut current = term
    let mut i = 0
    while i < 100 {
      result = result + current / denom
      current = current * term * term
      denom = denom + 2.0
      i = i + 1
    }
    2.0 * result
  }
}

///|
fn log2_mast(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    ln_mast(x) / ln_mast(2.0)
  }
}

///|
fn mast_exp(x : Double) -> Double {
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 20 {
    term = term * x / i.to_double()
    sum = sum + term
    i = i + 1
  }
  sum
}

///|
fn mean_d(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    0.0
  } else {
    let mut sum = 0.0
    let mut i = 0
    while i < n {
      sum = sum + arr[i]
      i = i + 1
    }
    sum / n.to_double()
  }
}

///|
fn var_d(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    0.0
  } else {
    let m = mean_d(arr)
    let mut sum = 0.0
    let mut i = 0
    while i < n {
      let diff = arr[i] - m
      sum = sum + diff * diff
      i = i + 1
    }
    sum / (n.to_double() - 1.0)
  }
}

///|
fn welch_t_test(x : Array[Double], y : Array[Double]) -> (Double, Double) {
  let n1 = x.length().to_double()
  let n2 = y.length().to_double()
  
  if n1 < 2.0 || n2 < 2.0 {
    return (0.0, 1.0)
  }
  
  let mean1 = mean_d(x)
  let mean2 = mean_d(y)
  let var1 = var_d(x)
  let var2 = var_d(y)
  
  if var1 == 0.0 && var2 == 0.0 {
    return (0.0, 1.0)
  }
  
  let se = sqrt_d_mast(var1 / n1 + var2 / n2)
  if se == 0.0 {
    return (0.0, 1.0)
  }
  
  let t_stat = (mean1 - mean2) / se
  
  let df_num = (var1 / n1 + var2 / n2) * (var1 / n1 + var2 / n2)
  let df_den = (var1 / n1) * (var1 / n1) / (n1 - 1.0) + (var2 / n2) * (var2 / n2) / (n2 - 1.0)
  let df = if df_den > 0.0 { df_num / df_den } else { n1 + n2 - 2.0 }
  
  let pval = t_test_pval(t_stat.abs(), df)
  
  (t_stat, pval)
}

///|
fn t_test_pval(t : Double, df : Double) -> Double {
  if df <= 0.0 {
    1.0
  } else {
    let x = df / (df + t * t)
    let p = beta_regularized(df / 2.0, 0.5, x)
    p
  }
}

///|
fn beta_regularized(a : Double, b : Double, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    let lbeta = ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
    let front = @math.pow(x, a) * @math.pow(1.0 - x, b) / mast_exp(lbeta) / a
    
    let mut result = front
    let mut term = front
    let mut n = 0
    while n < 200 {
      let numerator = (a + b + n.to_double()) * (a + n.to_double()) * x
      let denominator = (a + 2.0 * n.to_double()) * (a + 2.0 * n.to_double() + 1.0)
      term = term * numerator / denominator
      result = result + term
      
      if term.abs() < 0.0000000001 {
        break
      }
      
      n = n + 1
    }
    
    result
  }
}

///|
fn ln_gamma(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x < 1.0 {
    ln_gamma(x + 1.0) - ln_mast(x)
  } else {
    let z = x - 1.0
    let mut sum = 0.5 * ln_mast(2.0 * 3.141592653589793)
    sum = sum + (z + 0.5) * ln_mast(z + 1.0) - (z + 1.0)
    
    let coef = [
      1.0 / 12.0,
      -1.0 / 360.0,
      1.0 / 1260.0,
      -1.0 / 1680.0,
      1.0 / 1188.0,
    ]
    
    let mut i = 0
    while i < coef.length() {
      let divisor = @math.pow(z + 1.0, (2 * i + 1).to_double())
      sum = sum + coef[i] / divisor
      i = i + 1
    }
    
    sum
  }
}

///|
fn sqrt_d_mast(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    let mut z = x
    let mut prev = z
    let mut i = 0
    while i < 20 {
      z = (z + x / z) / 2.0
      if (z - prev).abs() < 0.0001 {
        break
      }
      prev = z
      i = i + 1
    }
    z
  }
}

///|
fn chisq_test_2x2(a : Int, b : Int, c : Int, d : Int) -> Double {
  let n = (a + b + c + d).to_double()
  if n == 0.0 {
    1.0
  } else {
    let e11 = (a + b).to_double() * (a + c).to_double() / n
    let e12 = (a + b).to_double() * (b + d).to_double() / n
    let e21 = (c + d).to_double() * (a + c).to_double() / n
    let e22 = (c + d).to_double() * (b + d).to_double() / n
    
    let chi2 = 
      (a.to_double() - e11) * (a.to_double() - e11) / e11 +
      (b.to_double() - e12) * (b.to_double() - e12) / e12 +
      (c.to_double() - e21) * (c.to_double() - e21) / e21 +
      (d.to_double() - e22) * (d.to_double() - e22) / e22
    
    chisq_pval(chi2, 1.0)
  }
}

///|
fn chisq_pval(chi2 : Double, df : Double) -> Double {
  if chi2 <= 0.0 || df <= 0.0 {
    1.0
  } else {
    let x = chi2 / 2.0
    let a = df / 2.0
    1.0 - gamma_regularized_p(a, x)
  }
}

///|
fn gamma_regularized_p(a : Double, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x < a + 1.0 {
    let mut sum = 1.0 / a
    let mut term = 1.0 / a
    let mut n = 1
    while n < 100 {
      term = term * x / (a + n.to_double())
      sum = sum + term
      if term.abs() < sum.abs() * 0.0000000001 {
        break
      }
      n = n + 1
    }
    sum * mast_exp(-x) * @math.pow(x, a) / mast_exp(ln_gamma(a))
  } else {
    1.0 - gamma_regularized_q(a, x)
  }
}

///|
fn gamma_regularized_q(a : Double, x : Double) -> Double {
  let mut w = 1.0
  let mut prev = 1.0
  let mut n = 0
  while n < 100 {
    let num = (n.to_double() - a) * ((n + 1).to_double() - a)
    let den = (2 * n + 1).to_double() + x - a
    w = den - num / w
    prev = prev * (den - num / w) / den
    n = n + 1
  }
  mast_exp(-x) * @math.pow(x, a) / (mast_exp(ln_gamma(a)) * w * x)
}

///|
fn fisher_exact_test(a : Int, b : Int, c : Int, d : Int) -> Double {
  let n = a + b + c + d
  let r1 = a + b
  let r2 = c + d
  let c1 = a + c
  
  let p_cutoff = hypergeometric_pmf(a, r1, r2, c1)
  
  let mut p_total = 0.0
  let min_a = mast_max_i(0, c1 - r2)
  let max_a = mast_min_i(r1, c1)
  
  let mut i = min_a
  while i <= max_a {
    let p = hypergeometric_pmf(i, r1, r2, c1)
    if p <= p_cutoff + 0.0000000001 {
      p_total = p_total + p
    }
    i = i + 1
  }
  
  if p_total > 1.0 { 1.0 } else { p_total }
}

///|
fn hypergeometric_pmf(k : Int, m : Int, n : Int, k_total : Int) -> Double {
  let log_m_choose_k = mast_log_choose(m, k)
  let log_n_choose_k = mast_log_choose(n, k_total - k)
  let log_total = mast_log_choose(m + n, k_total)
  mast_exp(log_m_choose_k + log_n_choose_k - log_total)
}

///|
fn mast_log_choose(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    -1000.0
  } else {
    ln_gamma((n + 1).to_double()) - ln_gamma((k + 1).to_double()) - ln_gamma((n - k + 1).to_double())
  }
}

///|
fn mast_max_i(a : Int, b : Int) -> Int {
  if a > b { a } else { b }
}

///|
pub fn zlm_hurdle(
  data : MastData,
  group1 : Int,
  group2 : Int,
) -> HurdleResult {
  let result = HurdleResult::new(data.n_genes)
  
  let mut g = 0
  while g < data.n_genes {
    result.gene[g] = data.gene_names[g]
    
    let group1_expr : Array[Double] = Array::new()
    let group2_expr : Array[Double] = Array::new()
    let group1_detected : Array[Bool] = Array::new()
    let group2_detected : Array[Bool] = Array::new()
    
    let mut i = 0
    while i < data.n_cells {
      if data.groups[i] == group1 {
        group1_expr.push(data.expression[i][g])
        group1_detected.push(data.expression[i][g] > 0.0)
      } else if data.groups[i] == group2 {
        group2_expr.push(data.expression[i][g])
        group2_detected.push(data.expression[i][g] > 0.0)
      }
      i = i + 1
    }
    
    let n1_detected = count_true(group1_detected)
    let n2_detected = count_true(group2_detected)
    let n1 = group1_detected.length()
    let n2 = group2_detected.length()
    
    let detection_rate = if n1 + n2 > 0 { 
      (n1_detected + n2_detected).to_double() / (n1 + n2).to_double() 
    } else { 
      0.0 
    }
    result.detection_rate[g] = detection_rate
    
    let discrete_pval = fisher_exact_test(
      n1_detected,
      n1 - n1_detected,
      n2_detected,
      n2 - n2_detected,
    )
    result.discrete_pval[g] = discrete_pval
    
    let group1_positive : Array[Double] = Array::new()
    let group2_positive : Array[Double] = Array::new()
    
    let mut j = 0
    while j < group1_expr.length() {
      if group1_expr[j] > 0.0 {
        group1_positive.push(group1_expr[j])
      }
      j = j + 1
    }
    
    let mut k = 0
    while k < group2_expr.length() {
      if group2_expr[k] > 0.0 {
        group2_positive.push(group2_expr[k])
      }
      k = k + 1
    }
    
    let continuous_pval = if group1_positive.length() >= 2 && group2_positive.length() >= 2 {
      let (_, pval) = welch_t_test(group1_positive, group2_positive)
      pval
    } else {
      1.0
    }
    result.continuous_pval[g] = continuous_pval
    
    let combined_stat = -2.0 * (ln_mast(discrete_pval + 0.000000000000000000000000000001) + ln_mast(continuous_pval + 0.000000000000000000000000000001))
    let combined_p = chisq_pval(combined_stat, 4.0)
    result.combined_pval[g] = combined_p
    
    let mean1 = mean_d(group1_expr)
    let mean2 = mean_d(group2_expr)
    let log2_fc = if mean2 > 0.0 { log2_mast((mean1 + 1.0) / (mean2 + 1.0)) } else { 0.0 }
    result.log2_fc[g] = log2_fc
    
    g = g + 1
  }
  
  let discrete_fdr = bh_fdr(result.discrete_pval)
  let continuous_fdr = bh_fdr(result.continuous_pval)
  let combined_fdr = bh_fdr(result.combined_pval)
  
  let mut m = 0
  while m < data.n_genes {
    result.discrete_fdr[m] = discrete_fdr[m]
    result.continuous_fdr[m] = continuous_fdr[m]
    result.combined_fdr[m] = combined_fdr[m]
    m = m + 1
  }
  
  result
}

///|
fn count_true(arr : Array[Bool]) -> Int {
  let mut count = 0
  let mut i = 0
  while i < arr.length() {
    if arr[i] {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
fn bh_fdr(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  let fdr : Array[Double] = Array::make(n, 1.0)
  
  let indices : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }
  
  let mut j = 0
  while j < n - 1 {
    let mut min_idx = j
    let mut min_p = pvalues[j]
    let mut k = j + 1
    while k < n {
      if pvalues[k] < min_p {
        min_p = pvalues[k]
        min_idx = k
      }
      k = k + 1
    }
    let temp = indices[j]
    indices[j] = indices[min_idx]
    indices[min_idx] = temp
    j = j + 1
  }
  
  let mut prev_fdr = 1.0
  let mut m = n - 1
  while m >= 0 {
    let p = pvalues[indices[m]]
    let fdr_val = p * n.to_double() / (m + 1).to_double()
    let corrected = if fdr_val < prev_fdr { fdr_val } else { prev_fdr }
    fdr[indices[m]] = if corrected > 1.0 { 1.0 } else { corrected }
    prev_fdr = fdr[indices[m]]
    m = m - 1
  }
  
  fdr
}

///|
pub struct MastSummary {
  n_de_combined : Int
  n_de_discrete : Int
  n_de_continuous : Int
  top_genes : Array[String]
  threshold : Double
}

pub fn MastSummary::new(result : HurdleResult, threshold : Double) -> MastSummary {
  let mut n_combined = 0
  let mut n_discrete = 0
  let mut n_continuous = 0
  
  let top : Array[String] = Array::new()
  let mut i = 0
  while i < result.gene.length() {
    if result.combined_fdr[i] < threshold {
      n_combined = n_combined + 1
    }
    if result.discrete_fdr[i] < threshold {
      n_discrete = n_discrete + 1
    }
    if result.continuous_fdr[i] < threshold {
      n_continuous = n_continuous + 1
    }
    
    if i < 10 {
      top.push(result.gene[i])
    }
    
    i = i + 1
  }
  
  MastSummary::{
    n_de_combined: n_combined,
    n_de_discrete: n_discrete,
    n_de_continuous: n_continuous,
    top_genes: top,
    threshold,
  }
}

pub fn MastSummary::to_string(self : MastSummary) -> String {
  let mut result = "MAST Summary\n"
  result = result + "============\n"
  result = result + "FDR threshold: " + self.threshold.to_string() + "\n\n"
  result = result + "DE genes (combined): " + self.n_de_combined.to_string() + "\n"
  result = result + "DE genes (discrete): " + self.n_de_discrete.to_string() + "\n"
  result = result + "DE genes (continuous): " + self.n_de_continuous.to_string() + "\n\n"
  result = result + "Top genes:\n"
  let mut i = 0
  while i < self.top_genes.length() && i < 5 {
    result = result + "  " + self.top_genes[i] + "\n"
    i = i + 1
  }
  result
}
