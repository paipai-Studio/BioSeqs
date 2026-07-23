///|
/// genefilter module, inspired by Bioconductor genefilter package.
/// Provides gene filtering functions with statistical tests.

///|
fn gf_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = if x > 1.0 { x } else { 1.0 }
  let mut i = 0
  while i < 60 {
    let next = (guess + x / guess) / 2.0
    let diff = if next > guess { next - guess } else { guess - next }
    guess = next
    if diff < 0.0000000000001 {
      return guess
    }
    i = i + 1
  }
  guess
}

///|
fn gf_sort_double(arr : Array[Double]) -> Array[Double] {
  let n = arr.length()
  let result = Array::make(n, 0.0)
  
  let mut i = 0
  while i < n {
    result[i] = arr[i]
    i = i + 1
  }
  
  i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      if result[j] < result[i] {
        let tmp = result[i]
        result[i] = result[j]
        result[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  
  result
}

///|
pub struct FilterResult {
  passing_genes : Array[String]
  statistics : Array[Double]
  p_values : Array[Double]
}

///|
pub fn FilterResult::new(passing_genes : Array[String], statistics : Array[Double], p_values : Array[Double]) -> FilterResult {
  FilterResult::{ passing_genes, statistics, p_values }
}

///|
pub struct GeneExpression {
  gene_ids : Array[String]
  expression : Array[Array[Double]]
  groups : Array[String]
}

///|
pub fn GeneExpression::new(gene_ids : Array[String], expression : Array[Array[Double]], groups : Array[String]) -> GeneExpression {
  GeneExpression::{ gene_ids, expression, groups }
}

///|
pub fn row_ttest(expression : GeneExpression, group1 : String, group2 : String) -> FilterResult {
  let passing_genes : Array[String] = Array::new()
  let statistics : Array[Double] = Array::new()
  let p_values : Array[Double] = Array::new()
  
  let mut i = 0
  while i < expression.gene_ids.length() {
    let expr = expression.expression[i]
    
    let group1_vals : Array[Double] = Array::new()
    let group2_vals : Array[Double] = Array::new()
    
    let mut j = 0
    while j < expression.groups.length() {
      if expression.groups[j] == group1 {
        group1_vals.push(expr[j])
      } else if expression.groups[j] == group2 {
        group2_vals.push(expr[j])
      }
      j = j + 1
    }
    
    if group1_vals.length() >= 2 && group2_vals.length() >= 2 {
      let n1 = group1_vals.length().to_double()
      let n2 = group2_vals.length().to_double()
      
      let mean1 = calc_mean(group1_vals)
      let mean2 = calc_mean(group2_vals)
      
      let var1 = calc_variance(group1_vals)
      let var2 = calc_variance(group2_vals)
      
      let se = gf_sqrt(var1 / n1 + var2 / n2)
      
      let t_stat = if se == 0.0 { 0.0 } else { (mean1 - mean2) / se }
      
      statistics.push(t_stat)
      p_values.push(0.0)
      passing_genes.push(expression.gene_ids[i])
    }
    
    i = i + 1
  }
  
  FilterResult::new(passing_genes, statistics, p_values)
}

///|
pub fn row_wilcoxon(expression : GeneExpression, group1 : String, group2 : String) -> FilterResult {
  let passing_genes : Array[String] = Array::new()
  let statistics : Array[Double] = Array::new()
  let p_values : Array[Double] = Array::new()
  
  let mut i = 0
  while i < expression.gene_ids.length() {
    let expr = expression.expression[i]
    
    let group1_vals : Array[Double] = Array::new()
    let group2_vals : Array[Double] = Array::new()
    
    let mut j = 0
    while j < expression.groups.length() {
      if expression.groups[j] == group1 {
        group1_vals.push(expr[j])
      } else if expression.groups[j] == group2 {
        group2_vals.push(expr[j])
      }
      j = j + 1
    }
    
    if group1_vals.length() >= 2 && group2_vals.length() >= 2 {
      let w_stat = wilcoxon_statistic(group1_vals, group2_vals)
      
      statistics.push(w_stat)
      p_values.push(0.0)
      passing_genes.push(expression.gene_ids[i])
    }
    
    i = i + 1
  }
  
  FilterResult::new(passing_genes, statistics, p_values)
}

///|
fn wilcoxon_statistic(group1 : Array[Double], group2 : Array[Double]) -> Double {
  let n1 = group1.length()
  let n2 = group2.length()
  
  let mut rank_sum = 0.0
  let mut i = 0
  while i < n1 {
    let x = group1[i]
    let mut rank = 1.0
    
    let mut j = 0
    while j < n1 {
      if j != i && group1[j] < x {
        rank = rank + 1.0
      } else if j != i && group1[j] == x {
        rank = rank + 0.5
      }
      j = j + 1
    }
    
    j = 0
    while j < n2 {
      if group2[j] < x {
        rank = rank + 1.0
      } else if group2[j] == x {
        rank = rank + 0.5
      }
      j = j + 1
    }
    
    rank_sum = rank_sum + rank
    i = i + 1
  }
  
  rank_sum
}

///|
pub fn variance_filter(expression : GeneExpression, threshold : Double) -> Array[String] {
  let passing_genes : Array[String] = Array::new()
  
  let mut i = 0
  while i < expression.gene_ids.length() {
    let expr = expression.expression[i]
    let var = calc_variance(expr)
    
    if var >= threshold {
      passing_genes.push(expression.gene_ids[i])
    }
    
    i = i + 1
  }
  
  passing_genes
}

///|
pub fn cv_filter(expression : GeneExpression, threshold : Double) -> Array[String] {
  let passing_genes : Array[String] = Array::new()
  
  let mut i = 0
  while i < expression.gene_ids.length() {
    let expr = expression.expression[i]
    let m = calc_mean(expr)
    let cv_val = if m == 0.0 { 0.0 } else { gf_sqrt(calc_variance(expr)) / m }
    
    if cv_val >= threshold {
      passing_genes.push(expression.gene_ids[i])
    }
    
    i = i + 1
  }
  
  passing_genes
}

///|
pub fn row_quantile_filter(expression : GeneExpression, q_low : Double, q_high : Double) -> Array[String] {
  let passing_genes : Array[String] = Array::new()
  
  let all_vals : Array[Double] = Array::new()
  
  let mut i = 0
  while i < expression.expression.length() {
    let mut j = 0
    while j < expression.expression[i].length() {
      all_vals.push(expression.expression[i][j])
      j = j + 1
    }
    i = i + 1
  }
  
  let sorted_vals = gf_sort_double(all_vals)
  
  let low_idx = (all_vals.length().to_double() * q_low).to_int()
  let high_idx = (all_vals.length().to_double() * q_high).to_int()
  
  let low_threshold = if low_idx < sorted_vals.length() { sorted_vals[low_idx] } else { 0.0 }
  let high_threshold = if high_idx < sorted_vals.length() { sorted_vals[high_idx] } else { 0.0 }
  
  let mut i = 0
  while i < expression.gene_ids.length() {
    let expr = expression.expression[i]
    let gene_min = calc_min(expr)
    let gene_max = calc_max(expr)
    
    if gene_min >= low_threshold && gene_max <= high_threshold {
      passing_genes.push(expression.gene_ids[i])
    }
    
    i = i + 1
  }
  
  passing_genes
}

///|
fn calc_mean(vals : Array[Double]) -> Double {
  if vals.length() == 0 {
    return 0.0
  }
  
  let mut sum = 0.0
  let mut i = 0
  while i < vals.length() {
    sum = sum + vals[i]
    i = i + 1
  }
  
  sum / vals.length().to_double()
}

///|
fn calc_variance(vals : Array[Double]) -> Double {
  if vals.length() < 2 {
    return 0.0
  }
  
  let m = calc_mean(vals)
  let mut sum_sq = 0.0
  
  let mut i = 0
  while i < vals.length() {
    sum_sq = sum_sq + @math.pow(vals[i] - m, 2.0)
    i = i + 1
  }
  
  sum_sq / (vals.length() - 1).to_double()
}

///|
fn calc_min(vals : Array[Double]) -> Double {
  if vals.length() == 0 {
    return 0.0
  }
  
  let mut m = vals[0]
  let mut i = 1
  while i < vals.length() {
    if vals[i] < m {
      m = vals[i]
    }
    i = i + 1
  }
  
  m
}

///|
fn calc_max(vals : Array[Double]) -> Double {
  if vals.length() == 0 {
    return 0.0
  }
  
  let mut m = vals[0]
  let mut i = 1
  while i < vals.length() {
    if vals[i] > m {
      m = vals[i]
    }
    i = i + 1
  }
  
  m
}

///|
pub fn create_example_gene_expression() -> GeneExpression {
  let gene_ids = ["gene1", "gene2", "gene3", "gene4", "gene5"]
  let expression = [
    [10.5, 12.3, 11.1, 25.0, 27.5, 26.2],
    [5.8, 6.2, 5.5, 7.0, 6.8, 7.2],
    [100.2, 98.5, 102.1, 50.0, 48.5, 52.0],
    [15.7, 16.2, 15.9, 14.8, 15.2, 14.5],
    [3.2, 3.5, 3.1, 3.0, 3.3, 2.9],
  ]
  let groups = ["control", "control", "control", "treatment", "treatment", "treatment"]
  
  GeneExpression::new(gene_ids, expression, groups)
}