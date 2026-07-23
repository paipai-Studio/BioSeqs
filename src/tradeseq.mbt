///|
/// tradeSeq - Trajectory-based differential expression analysis.
/// Inspired by Bioconductor tradeSeq.
/// Provides methods for differential expression analysis along
/// trajectories using generalized additive models (GAMs),
/// including gene expression smoothing, condition testing,
/// and trajectory differential expression (tradeSeq) analysis.

///|
pub struct TrajectoryPoint {
  cell_id : String
  pseudotime : Double
  expression : Array[Double]
  condition : String
}

///|
pub struct GeneExpressionData {
  gene_id : String
  expression : Array[Double]
  pseudotime : Array[Double]
  condition : Array[String]
}

///|
pub struct GAMFit {
  gene_id : String
  coefficients : Array[Double]
  fitted_values : Array[Double]
  residuals : Array[Double]
  r_squared : Double
  degrees_of_freedom : Double
  smoothing_basis : Array[Double]
}

///|
pub struct DifferentialExpressionResult {
  gene_id : String
  p_value : Double
  adjusted_p_value : Double
  coefficient : Double
  standard_error : Double
  z_statistic : Double
  significant : Bool
}

///|
pub struct TradeSeqResult {
  gene_results : Array[DifferentialExpressionResult]
  n_genes : Int
  n_significant : Int
  fdr_threshold : Double
}

///|
pub struct SmoothResult {
  gene_id : String
  pseudotime : Array[Double]
  fitted : Array[Double]
  lower_ci : Array[Double]
  upper_ci : Array[Double]
}

///|
pub struct TrajectoryData {
  points : Array[TrajectoryPoint]
  genes : Array[String]
  conditions : Array[String]
  n_points : Int
  n_genes : Int
}

///|
pub fn TrajectoryPoint::new(
  cell_id : String,
  pseudotime : Double,
  expression : Array[Double],
  condition : String
) -> TrajectoryPoint {
  TrajectoryPoint::{
    cell_id,
    pseudotime,
    expression,
    condition
  }
}

///|
pub fn GeneExpressionData::new(
  gene_id : String,
  expression : Array[Double],
  pseudotime : Array[Double],
  condition : Array[String]
) -> GeneExpressionData {
  GeneExpressionData::{
    gene_id,
    expression,
    pseudotime,
    condition
  }
}

///|
pub fn TrajectoryData::new() -> TrajectoryData {
  TrajectoryData::{
    points: Array::new(),
    genes: Array::new(),
    conditions: Array::new(),
    n_points: 0,
    n_genes: 0
  }
}

///|
pub fn TrajectoryData::add_point(self : TrajectoryData, point : TrajectoryPoint) -> TrajectoryData {
  let new_points = self.points
  new_points.push(point)
  
  let new_conditions = self.conditions
  if !tradeseq_contains_string(new_conditions, point.condition) {
    new_conditions.push(point.condition)
  }
  
  TrajectoryData::{
    points: new_points,
    genes: self.genes,
    conditions: new_conditions,
    n_points: self.n_points + 1,
    n_genes: self.n_genes
  }
}

///|
pub fn TrajectoryData::add_gene(self : TrajectoryData, gene : String) -> TrajectoryData {
  let new_genes = self.genes
  if !tradeseq_contains_string(new_genes, gene) {
    new_genes.push(gene)
  }
  
  TrajectoryData::{
    points: self.points,
    genes: new_genes,
    conditions: self.conditions,
    n_points: self.n_points,
    n_genes: new_genes.length()
  }
}

///|
pub fn fit_gam(gene_data : GeneExpressionData, n_knots : Int) -> GAMFit {
  let n = gene_data.expression.length()
  let expression = gene_data.expression
  let pseudotime = gene_data.pseudotime
  
  if n < 3 {
    return GAMFit::{
      gene_id: gene_data.gene_id,
      coefficients: Array::new(),
      fitted_values: expression,
      residuals: Array::new(),
      r_squared: 0.0,
      degrees_of_freedom: 1.0,
      smoothing_basis: Array::new()
    }
  }
  
  // Sort by pseudotime
  let sorted_indices = tradeseq_sort_indices(pseudotime)
  let sorted_expression : Array[Double] = Array::new()
  let sorted_pseudotime : Array[Double] = Array::new()
  
  let mut i = 0
  while i < n {
    sorted_expression.push(expression[sorted_indices[i]])
    sorted_pseudotime.push(pseudotime[sorted_indices[i]])
    i = i + 1
  }
  
  // Create basis functions for cubic spline
  let n_knots_actual = if n_knots > n / 2 { n / 2 } else { n_knots }
  let basis = tradeseq_create_basis(sorted_pseudotime, n_knots_actual)
  
  // Fit using least squares (simplified)
  let coefficients = tradeseq_fit_coefficients(basis, sorted_expression)
  
  // Calculate fitted values
  let fitted_values = tradeseq_compute_fitted(basis, coefficients)
  
  // Calculate residuals
  let residuals : Array[Double] = Array::new()
  i = 0
  while i < n {
    residuals.push(sorted_expression[i] - fitted_values[i])
    i = i + 1
  }
  
  // Calculate R-squared
  let r_squared = tradeseq_calculate_r_squared(sorted_expression, fitted_values)
  
  GAMFit::{
    gene_id: gene_data.gene_id,
    coefficients,
    fitted_values,
    residuals,
    r_squared,
    degrees_of_freedom: n_knots_actual.to_double() + 1.0,
    smoothing_basis: basis[0]  // Store first row as representative
  }
}

///|
pub fn trade_test_condition_effect(
  gene_data : GeneExpressionData,
  pseudotime_point : Double,
  condition1 : String,
  condition2 : String
) -> DifferentialExpressionResult {
  let n = gene_data.expression.length()
  
  // Calculate mean expression for each condition at nearby pseudotime points
  let mut sum1 = 0.0
  let mut count1 = 0
  let mut sum2 = 0.0
  let mut count2 = 0
  
  let bandwidth = 0.1
  let mut i = 0
  while i < n {
    let pt_diff = (gene_data.pseudotime[i] - pseudotime_point).abs()
    if pt_diff < bandwidth {
      if gene_data.condition[i] == condition1 {
        sum1 = sum1 + gene_data.expression[i]
        count1 = count1 + 1
      } else if gene_data.condition[i] == condition2 {
        sum2 = sum2 + gene_data.expression[i]
        count2 = count2 + 1
      }
    }
    i = i + 1
  }
  
  let mean1 = if count1 > 0 { sum1 / count1.to_double() } else { 0.0 }
  let mean2 = if count2 > 0 { sum2 / count2.to_double() } else { 0.0 }
  
  // Calculate standard errors
  let var1 = if count1 > 1 {
    let mut s = 0.0
    let mut j = 0
    while j < n {
      if (gene_data.pseudotime[j] - pseudotime_point).abs() < bandwidth && gene_data.condition[j] == condition1 {
        s = s + (gene_data.expression[j] - mean1).pow(2.0)
        j = j + 1
      } else {
        j = j + 1
      }
    }
    s / (count1 - 1).to_double()
  } else { 1.0 }
  
  let var2 = if count2 > 1 {
    let mut s = 0.0
    let mut j = 0
    while j < n {
      if (gene_data.pseudotime[j] - pseudotime_point).abs() < bandwidth && gene_data.condition[j] == condition2 {
        s = s + (gene_data.expression[j] - mean2).pow(2.0)
        j = j + 1
      } else {
        j = j + 1
      }
    }
    s / (count2 - 1).to_double()
  } else { 1.0 }
  
  let se = (var1 / count1.to_double() + var2 / count2.to_double()).sqrt()
  let coefficient = mean1 - mean2
  
  let z_statistic = if se > 0.0 { coefficient / se } else { 0.0 }
  let p_value = tradeseq_normal_p_value(z_statistic)
  
  DifferentialExpressionResult::{
    gene_id: gene_data.gene_id,
    p_value,
    adjusted_p_value: p_value,  // Will be adjusted later
    coefficient,
    standard_error: se,
    z_statistic,
    significant: p_value < 0.05
  }
}

///|
pub fn run_tradeseq_analysis(
  trajectory_data : TrajectoryData,
  genes : Array[String],
  conditions : Array[String],
  n_knots : Int,
  fdr : Double
) -> TradeSeqResult {
  let gene_results : Array[DifferentialExpressionResult] = Array::new()
  
  // For each gene, perform differential expression test
  let mut i = 0
  while i < genes.length() {
    let gene = genes[i]
    
    // Create gene expression data
    let expression : Array[Double] = Array::new()
    let pseudotime : Array[Double] = Array::new()
    let conditions_arr : Array[String] = Array::new()
    
    let mut j = 0
    while j < trajectory_data.points.length() {
      let point = trajectory_data.points[j]
      expression.push(point.expression[i])
      pseudotime.push(point.pseudotime)
      conditions_arr.push(point.condition)
      j = j + 1
    }
    
    let gene_data = GeneExpressionData::new(gene, expression, pseudotime, conditions_arr)
    
    // Test at multiple pseudotime points
    let n_test_points = 10
    let mut min_pt = 1.0e9
    let mut max_pt = -1.0e9
    let mut k = 0
    while k < pseudotime.length() {
      if pseudotime[k] < min_pt { min_pt = pseudotime[k] }
      if pseudotime[k] > max_pt { max_pt = pseudotime[k] }
      k = k + 1
    }
    
    let pt_range = max_pt - min_pt
    let mut min_p_value = 1.0
    let mut best_result : DifferentialExpressionResult = DifferentialExpressionResult::{
      gene_id: gene,
      p_value: 1.0,
      adjusted_p_value: 1.0,
      coefficient: 0.0,
      standard_error: 1.0,
      z_statistic: 0.0,
      significant: false
    }
    
    // Check condition ordering
    let cond1 = conditions[0]
    let cond2 = if conditions.length() > 1 { conditions[1] } else { "control" }
    
    let mut l = 0
    while l < n_test_points {
      let pt = min_pt + (pt_range * (l.to_double() + 0.5) / n_test_points.to_double())
      let result = trade_test_condition_effect(gene_data, pt, cond1, cond2)
      
      if result.p_value < min_p_value {
        min_p_value = result.p_value
        best_result = result
      }
      l = l + 1
    }
    
    gene_results.push(best_result)
    i = i + 1
  }
  
  // Apply Benjamini-Hochberg correction
  let adjusted_results = tradeseq_bh_correction(gene_results, fdr)
  
  // Count significant
  let mut n_significant = 0
  let m = adjusted_results.length()
  let mut i2 = 0
  while i2 < m {
    if adjusted_results[i2].significant {
      n_significant = n_significant + 1
    }
    i2 = i2 + 1
  }
  
  TradeSeqResult::{
    gene_results: adjusted_results,
    n_genes: m,
    n_significant,
    fdr_threshold: fdr
  }
}

///|
pub fn calculate_gene_smooth(
  trajectory_data : TrajectoryData,
  gene : String,
  n_knots : Int,
  n_points : Int
) -> SmoothResult {
  let expression : Array[Double] = Array::new()
  let pseudotime : Array[Double] = Array::new()
  
  let gene_idx = tradeseq_find_gene_index(trajectory_data.genes, gene)
  
  let mut i = 0
  while i < trajectory_data.points.length() {
    let point = trajectory_data.points[i]
    expression.push(point.expression[gene_idx])
    pseudotime.push(point.pseudotime)
    i = i + 1
  }
  
  let gene_data = GeneExpressionData::new(gene, expression, pseudotime, trajectory_data.conditions)
  let gam_fit = fit_gam(gene_data, n_knots)
  
  // Generate smooth curve at evenly spaced pseudotime points
  let mut min_pt = 1.0e9
  let mut max_pt = -1.0e9
  let mut j = 0
  while j < pseudotime.length() {
    if pseudotime[j] < min_pt { min_pt = pseudotime[j] }
    if pseudotime[j] > max_pt { max_pt = pseudotime[j] }
    j = j + 1
  }
  
  let smooth_pseudotime : Array[Double] = Array::new()
  let smooth_fitted : Array[Double] = Array::new()
  let smooth_lower : Array[Double] = Array::new()
  let smooth_upper : Array[Double] = Array::new()
  
  let pt_range = max_pt - min_pt
  let k = gam_fit.degrees_of_freedom.to_int()
  let basis = tradeseq_create_basis(pseudotime, if k > 0 { k } else { 4 })
  
  let mut l = 0
  while l < n_points {
    let pt = min_pt + (pt_range * l.to_double() / (n_points.to_double() - 1.0))
    smooth_pseudotime.push(pt)
    
    // Simple interpolation based on GAM fit
    let fit_idx = if gam_fit.fitted_values.length() > 0 {
      let idx = ((pt - min_pt) / pt_range * (gam_fit.fitted_values.length() - 1).to_double()).to_int()
      if idx < 0 { 0 } else if idx >= gam_fit.fitted_values.length() { gam_fit.fitted_values.length() - 1 } else { idx }
    } else { 0 }
    
    let fitted = if gam_fit.fitted_values.length() > 0 { gam_fit.fitted_values[fit_idx] } else { 0.0 }
    let se = 0.1  // Simplified standard error
    
    smooth_fitted.push(fitted)
    smooth_lower.push(fitted - 1.96 * se)
    smooth_upper.push(fitted + 1.96 * se)
    l = l + 1
  }
  
  SmoothResult::{
    gene_id: gene,
    pseudotime: smooth_pseudotime,
    fitted: smooth_fitted,
    lower_ci: smooth_lower,
    upper_ci: smooth_upper
  }
}

///|
pub fn create_example_trajectory_data() -> TrajectoryData {
  let mut data = TrajectoryData::new()
  
  // Add genes
  data = data.add_gene("GeneA")
  data = data.add_gene("GeneB")
  data = data.add_gene("GeneC")
  
  // Create 50 trajectory points with varying expression
  let n_cells = 50
  let mut i = 0
  while i < n_cells {
    let pseudotime = i.to_double() / (n_cells.to_double() - 1.0)
    let condition = if i < 25 { "control" } else { "treatment" }
    
    // GeneA: increases with pseudotime
    let gene_a = 0.5 + pseudotime * 2.0 + @math.sin(i.to_double() * 0.01) * 0.3
    // GeneB: decreases then increases
    let gene_b = 1.5 - pseudotime * 1.0 + @math.cos(i.to_double() * 0.02) * 0.5
    // GeneC: different between conditions
    let gene_c = if condition == "control" {
      1.0 + @math.sin(i.to_double() * 0.05) * 0.2
    } else {
      2.0 + @math.sin(i.to_double() * 0.05) * 0.2
    }
    
    let point = TrajectoryPoint::new(
      "cell_" + i.to_string(),
      pseudotime,
      [gene_a, gene_b, gene_c],
      condition
    )
    
    data = data.add_point(point)
    i = i + 1
  }
  
  data
}

///|
pub fn tradeseq_summary(result : TradeSeqResult) -> String {
  let mut summary = "\ntradeSeq Analysis Summary\n"
  summary = summary + "==============================\n"
  summary = summary + "Total genes tested: " + result.n_genes.to_string() + "\n"
  summary = summary + "Significant genes (FDR < " + result.fdr_threshold.to_string() + "): " + result.n_significant.to_string() + "\n"
  
  // Sort results by p-value
  let n = result.gene_results.length()
  let sorted_indices = tradeseq_sort_by_pvalue(result.gene_results)
  
  summary = summary + "\nTop significant genes:\n"
  let n_top = if 10 < n { 10 } else { n }
  let mut k = 0
  while k < n_top {
    let idx = sorted_indices[k]
    if result.gene_results[idx].significant {
      summary = summary + "  " + result.gene_results[idx].gene_id + 
        ": adj p-value = " + result.gene_results[idx].adjusted_p_value.to_string() +
        ", coefficient = " + result.gene_results[idx].coefficient.to_string() + "\n"
    }
    k = k + 1
  }
  
  summary
}

// Helper functions

fn tradeseq_contains_string(arr : Array[String], value : String) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == value {
      return true
    }
    i = i + 1
  }
  false
}

fn tradeseq_sort_indices(arr : Array[Double]) -> Array[Int] {
  let n = arr.length()
  let indices : Array[Int] = Array::new()
  let values : Array[Double] = Array::new()
  
  // Copy values and indices
  let mut i = 0
  while i < n {
    indices.push(i)
    values.push(arr[i])
    i = i + 1
  }
  
  // Bubble sort by value
  i = 0
  while i < n {
    let mut j = 0
    while j < n - i - 1 {
      if values[j] > values[j + 1] {
        // Swap
        let temp_val = values[j]
        values[j] = values[j + 1]
        values[j + 1] = temp_val
        
        let temp_idx = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = temp_idx
      }
      j = j + 1
    }
    i = i + 1
  }
  
  indices
}

fn tradeseq_create_basis(pseudotime : Array[Double], n_knots : Int) -> Array[Array[Double]] {
  let n = pseudotime.length()
  let n_basis = n_knots + 2  // Interior knots + boundary
  
  // Create natural cubic spline basis
  let basis : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let x = pseudotime[i]
    
    let mut j = 0
    while j < n_basis {
      // Simple B-spline basis (simplified)
      let knot_pos = j.to_double() / (n_basis.to_double() - 1.0)
      let dist = (x - knot_pos).abs()
      let basis_val = if dist < 0.5 {
        1.0 - dist * 2.0  // Linear basis
      } else {
        0.0
      }
      row.push(basis_val)
      j = j + 1
    }
    basis.push(row)
    i = i + 1
  }
  
  basis
}

fn tradeseq_fit_coefficients(basis : Array[Array[Double]], y : Array[Double]) -> Array[Double] {
  let n = basis.length()
  let n_basis = basis[0].length()
  
  // Simple least squares approximation
  // X'X β = X'y
  let coefficients : Array[Double] = Array::new()
  
  let mut j = 0
  while j < n_basis {
    let mut sum_x = 0.0
    let mut sum_xy = 0.0
    let mut i = 0
    while i < n {
      sum_x = sum_x + basis[i][j] * basis[i][j]
      sum_xy = sum_xy + basis[i][j] * y[i]
      i = i + 1
    }
    
    let coef = if sum_x > 0.0 { sum_xy / sum_x } else { 0.0 }
    coefficients.push(coef)
    j = j + 1
  }
  
  coefficients
}

fn tradeseq_compute_fitted(basis : Array[Array[Double]], coefficients : Array[Double]) -> Array[Double] {
  let n = basis.length()
  let n_basis = coefficients.length()
  
  let fitted : Array[Double] = Array::new()
  
  let mut i = 0
  while i < n {
    let mut sum = 0.0
    let mut j = 0
    while j < n_basis {
      sum = sum + basis[i][j] * coefficients[j]
      j = j + 1
    }
    fitted.push(sum)
    i = i + 1
  }
  
  fitted
}

fn tradeseq_calculate_r_squared(observed : Array[Double], fitted : Array[Double]) -> Double {
  let n = observed.length()
  
  // Calculate mean of observed
  let mut mean_y = 0.0
  let mut i = 0
  while i < n {
    mean_y = mean_y + observed[i]
    i = i + 1
  }
  mean_y = mean_y / n.to_double()
  
  // Calculate total sum of squares
  let mut ss_total = 0.0
  i = 0
  while i < n {
    ss_total = ss_total + (observed[i] - mean_y).pow(2.0)
    i = i + 1
  }
  
  // Calculate residual sum of squares
  let mut ss_residual = 0.0
  i = 0
  while i < n {
    ss_residual = ss_residual + (observed[i] - fitted[i]).pow(2.0)
    i = i + 1
  }
  
  if ss_total == 0.0 { 0.0 } else { 1.0 - ss_residual / ss_total }
}

fn tradeseq_normal_p_value(z : Double) -> Double {
  // Two-tailed p-value from standard normal distribution
  let abs_z = if z >= 0.0 { z } else { -z }
  2.0 * (1.0 - tradeseq_normal_cdf(abs_z))
}

fn tradeseq_normal_cdf(x : Double) -> Double {
  // Standard normal CDF approximation
  let t = 1.0 / (1.0 + 0.2316419 * x)
  let d = 0.3989422804014327 * @math.exp(-x * x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  p
}

fn tradeseq_sort_by_pvalue(results : Array[DifferentialExpressionResult]) -> Array[Int] {
  let n = results.length()
  let indices : Array[Int] = Array::new()
  let values : Array[Double] = Array::new()
  
  // Initialize
  let mut i = 0
  while i < n {
    indices.push(i)
    values.push(results[i].p_value)
    i = i + 1
  }
  
  // Bubble sort by p-value
  i = 0
  while i < n {
    let mut j = 0
    while j < n - i - 1 {
      if values[j] > values[j + 1] {
        // Swap values
        let temp_val = values[j]
        values[j] = values[j + 1]
        values[j + 1] = temp_val
        
        // Swap indices
        let temp_idx = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = temp_idx
      }
      j = j + 1
    }
    i = i + 1
  }
  
  indices
}

fn tradeseq_bh_correction(results : Array[DifferentialExpressionResult], fdr : Double) -> Array[DifferentialExpressionResult] {
  let n = results.length()
  
  // Create sorted order by p-value
  let sorted_indices = tradeseq_sort_by_pvalue(results)
  
  // Extract sorted p-values
  let sorted_pvalues : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    sorted_pvalues.push(results[sorted_indices[i]].p_value)
    i = i + 1
  }
  
  // Apply BH correction
  let adjusted : Array[Double] = Array::new()
  let mut k = 0
  while k < n {
    let rank = (k + 1).to_double()
    let adj_p = sorted_pvalues[k] * n.to_double() / rank
    // Cap at 1.0
    let capped = if adj_p > 1.0 { 1.0 } else { adj_p }
    adjusted.push(capped)
    k = k + 1
  }
  
  // Enforce monotonicity
  let mono = adjusted
  let mut l = n - 2
  while l >= 0 {
    if mono[l] > mono[l + 1] {
      mono[l] = mono[l + 1]
    }
    l = l - 1
  }
  
  // Map back to original positions
  let final_results : Array[DifferentialExpressionResult] = Array::new()
  let mut m = 0
  while m < n {
    let orig_idx = sorted_indices[m]
    let adj_p = mono[m]
    final_results.push(DifferentialExpressionResult::{
      ..results[orig_idx],
      adjusted_p_value: adj_p,
      significant: adj_p < fdr
    })
    m = m + 1
  }
  
  final_results
}

fn tradeseq_find_gene_index(genes : Array[String], gene : String) -> Int {
  let mut i = 0
  while i < genes.length() {
    if genes[i] == gene {
      return i
    }
    i = i + 1
  }
  0
}
