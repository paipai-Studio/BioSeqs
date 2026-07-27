///|
/// Bio.PopGen Advanced - Advanced population genetics statistics
/// and neutrality tests.
///
/// This module implements classic neutrality tests including:
/// - Tajima's D (Tajima, 1989)
/// - Fu & Li's D and F (Fu & Li, 1993)
/// - McDonald-Kreitman test (McDonald & Kreitman, 1991)
/// - Watterson's theta (1975)
/// - Various allele frequency spectrum statistics
///
/// These tests are fundamental for detecting selection and
/// demographic effects in population genetic data.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// Represents a single nucleotide polymorphism (SNP).
pub struct PolymorphicSite {
  /// Position in the sequence (1-based)
  position : Int
  /// Derived allele frequency (0.0 to 1.0)
  derived_freq : Double
  /// Number of derived alleles observed
  derived_count : Int
  /// Total sample size (number of chromosomes)
  sample_size : Int
  /// Whether this is a singleton (derived count = 1)
  is_singleton : Bool
  /// Whether this is a replacement (nonsynonymous) site
  is_replacement : Bool
}

///|
/// Configuration for a neutrality test.
pub struct NeutralityTestConfig {
  /// Sample size (number of sequences/individuals)
  sample_size : Int
  /// Number of sequence sites (length of alignment)
  sequence_length : Int
  /// Mutation rate per site
  mutation_rate : Double
  /// Whether to use folded frequency spectrum
  use_folded : Bool
}

///|
/// Default neutrality test configuration.
pub fn NeutralityTestConfig::default(sample_size : Int, sequence_length : Int) -> NeutralityTestConfig {
  NeutralityTestConfig::{
    sample_size,
    sequence_length,
    mutation_rate: 0.0,
    use_folded: false,
  }
}

///|
/// Result of a neutrality test.
pub struct NeutralityTestResult {
  /// Test name
  test_name : String
  /// Test statistic value
  statistic : Double
  /// p-value (if calculable)
  p_value : Double?
  /// Whether the test is significant (p < 0.05)
  is_significant : Bool
  /// Interpretation
  interpretation : String
  /// Standard error of the statistic
  standard_error : Double
  /// Sample size used
  sample_size : Int
}

///|
/// Result of the McDonald-Kreitman test.
pub struct MKTestResult {
  /// Number of nonsynonymous polymorphisms (within species)
  p_nonsyn : Int
  /// Number of synonymous polymorphisms (within species)
  p_syn : Int
  /// Number of nonsynonymous substitutions (between species)
  d_nonsyn : Int
  /// Number of synonymous substitutions (between species)
  d_syn : Int
  /// Neutrality index (NI = (p_nonsyn/p_syn) / (d_nonsyn/d_syn))
  neutrality_index : Double
  /// p-value from Fisher's exact test
  p_value : Double
  /// Interpretation
  interpretation : String
  /// Whether selection is detected
  has_selection : Bool
}

///|
/// Full allele frequency spectrum results.
pub struct AlleleFrequencySpectrum {
  /// Sample size
  sample_size : Int
  /// Number of singletons (i = 1)
  singletons : Int
  /// Number of doubletons (i = 2)
  doubletons : Int
  /// Total number of segregating sites (S)
  segregating_sites : Int
  /// Total number of derived mutations (eta)
  total_mutations : Int
  /// The folded frequency spectrum
  folded_spectrum : Array[Int]
  /// The unfolded frequency spectrum
  unfolded_spectrum : Array[Int]
}

// ============================================================================
// PolymorphicSite Implementation
// ============================================================================

///|
/// Create a new PolymorphicSite.
pub fn new_polymorphic_site(
  position : Int,
  derived_count : Int,
  sample_size : Int,
  is_replacement : Bool,
) -> PolymorphicSite {
  let derived_freq = if sample_size > 0 {
    derived_count.to_double() / sample_size.to_double()
  } else {
    0.0
  }
  PolymorphicSite::{
    position,
    derived_freq,
    derived_count,
    sample_size,
    is_singleton: derived_count == 1,
    is_replacement,
  }
}

///|
/// Create a singleton site.
pub fn new_singleton(position : Int, sample_size : Int, is_replacement : Bool) -> PolymorphicSite {
  new_polymorphic_site(position, 1, sample_size, is_replacement)
}

///|
/// Create a doubleton site.
pub fn new_doubleton(position : Int, sample_size : Int, is_replacement : Bool) -> PolymorphicSite {
  new_polymorphic_site(position, 2, sample_size, is_replacement)
}

// ============================================================================
// Watterson's Theta
// ============================================================================

///|
/// Calculate Watterson's theta (theta_W = S / a1).
/// This is the standard neutral estimator of theta
/// based on the number of segregating sites.
///
/// Reference: Watterson, G. A. (1975).
/// "On the number of segregating sites in genetical models
/// without recombination." Theor. Popul. Biol. 7:256-276.
pub fn popgen_watterson_theta(segregating_sites : Int, sample_size : Int) -> Double {
  if sample_size <= 1 || segregating_sites == 0 {
    return 0.0
  }
  let a1 = watterson_a1(sample_size)
  if a1 == 0.0 {
    return 0.0
  }
  segregating_sites.to_double() / a1
}

///|
/// Calculate the harmonic sum a1 = sum(1/i) for i = 1 to n-1.
/// Used in Watterson's theta calculation.
pub fn watterson_a1(sample_size : Int) -> Double {
  if sample_size <= 1 {
    return 0.0
  }
  let mut sum = 0.0
  for i = 1; i < sample_size; i = i + 1 {
    sum += 1.0 / i.to_double()
  }
  sum
}

///|
/// Calculate theta pi (average number of pairwise differences).
/// This is Tajima's estimator of theta.
pub fn tajima_theta_pi(pairwise_differences : Double) -> Double {
  pairwise_differences
}

// ============================================================================
// Tajima's D
// ============================================================================

///|
/// Calculate Tajima's D statistic.
///
/// D = (theta_pi - theta_W) / standard_error
///
/// A negative D indicates an excess of low-frequency variants
/// (suggesting population expansion or positive selection).
/// A positive D indicates a deficiency of low-frequency variants
/// (suggesting population contraction or balancing selection).
///
/// Reference: Tajima, F. (1989).
/// "Statistical method for testing the neutral mutation hypothesis
/// by DNA polymorphism." Genetics 123:585-595.
pub fn popgen_tajima_d(
  theta_pi : Double,
  segregating_sites : Int,
  sample_size : Int,
  sequence_length : Int,
) -> NeutralityTestResult {
  let theta_w = popgen_watterson_theta(segregating_sites, sample_size)
  let se = tajima_d_standard_error(sample_size, segregating_sites, theta_w)
  
  let d_value = if se > 0.0 {
    (theta_pi - theta_w) / se
  } else {
    0.0
  }
  
  // Approximate p-value using beta distribution
  let p_val = approximate_tajima_d_pvalue(d_value, sample_size)
  
  let interpretation = if d_value < -1.0 {
    "Negative D: excess of low-frequency variants, suggests population expansion or positive selection"
  } else if d_value > 1.0 {
    "Positive D: deficiency of low-frequency variants, suggests population contraction or balancing selection"
  } else {
    "Neutral pattern: consistent with equilibrium neutrality"
  }
  
  NeutralityTestResult::{
    test_name: "Tajima's D",
    statistic: d_value,
    p_value: p_val,
    is_significant: if p_val != None { p_val.unwrap() < 0.05 } else { false },
    interpretation,
    standard_error: se,
    sample_size,
  }
}

///|
/// Calculate the standard error of Tajima's D.
pub fn tajima_d_standard_error(n : Int, s : Int, theta_w : Double) -> Double {
  if n <= 1 || s <= 0 || theta_w <= 0.0 {
    return 0.0
  }
  
  let a1 = watterson_a1(n)
  let a2 = watterson_a2(n)
  let b1 = tajima_b1(n)
  let b2 = tajima_b2(n)
  
  // Var(theta_pi - theta_W) ≈ c1*S + c2*S*(S-1)
  let theta_var = (b1 / a1) * theta_w + (b2 / (a1 * a2)) * theta_w * theta_w
  theta_var.sqrt()
}

///|
/// Calculate a2 = sum(1/i^2) for i = 1 to n-1.
pub fn watterson_a2(sample_size : Int) -> Double {
  if sample_size <= 1 {
    return 0.0
  }
  let mut sum = 0.0
  for i = 1; i < sample_size; i = i + 1 {
    sum += 1.0 / (i.to_double() * i.to_double())
  }
  sum
}

///|
/// Calculate b1 for Tajima's D variance.
pub fn tajima_b1(n : Int) -> Double {
  if n <= 1 {
    return 0.0
  }
  // b1 = (n+1)/(3*(n-1))
  (n.to_double() + 1.0) / (3.0 * (n.to_double() - 1.0))
}

///|
/// Calculate b2 for Tajima's D variance.
pub fn tajima_b2(n : Int) -> Double {
  if n <= 1 {
    return 0.0
  }
  let a1 = watterson_a1(n)
  let a2 = watterson_a2(n)
  // b2 = (2*(n^2 + n + 3)) / (9*n*(n-1)) - (a2/(a1^2))
  let term1 = 2.0 * (n.to_double() * n.to_double() + n.to_double() + 3.0) / (9.0 * n.to_double() * (n.to_double() - 1.0))
  let term2 = a2 / (a1 * a1)
  term1 - term2
}

///|
/// Approximate p-value for Tajima's D.
/// Uses the beta distribution approximation.
pub fn approximate_tajima_d_pvalue(d : Double, n : Int) -> Double? {
  if n <= 1 {
    return None
  }
  
  // Approximation: D approximately follows a scaled beta distribution
  // Mean ≈ 0, Var ≈ 1
  // Transform to approximate normal
  let abs_d = d.abs()
  
  // Simple approximation for p-value
  // Two-sided test
  let p = 2.0 * (1.0 - popgen_normal_cdf(abs_d))
  
  if p < 0.0001 {
    Some(0.0001)
  } else if p > 1.0 {
    Some(1.0)
  } else {
    Some(p)
  }
}

///|
/// Approximate standard normal CDF.
pub fn popgen_normal_cdf(x : Double) -> Double {
  // Abramowitz and Stegun approximation
  let t = 1.0 / (1.0 + 0.2316419 * x.abs())
  let d = 0.3989422804014327 * @math.exp(-x * x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if x > 0.0 {
    1.0 - p
  } else {
    p
  }
}

// ============================================================================
// Fu & Li's D Test
// ============================================================================

///|
/// Calculate Fu & Li's D statistic.
///
/// D = (n * singleton - S) / standard_error
///
/// Compares the number of singleton (derived) mutations
/// to the total number of segregating sites.
///
/// Reference: Fu, Y.-X. & Li, W.-H. (1993).
/// "Statistical tests of neutrality of mutations."
/// Genetics 133:693-709.
pub fn fu_li_d(
  singletons : Int,
  segregating_sites : Int,
  sample_size : Int,
) -> NeutralityTestResult {
  if sample_size <= 1 {
    return NeutralityTestResult::{
      test_name: "Fu & Li's D",
      statistic: 0.0,
      p_value: None,
      is_significant: false,
      interpretation: "Insufficient sample size",
      standard_error: 0.0,
      sample_size,
    }
  }
  
  let n = sample_size.to_double()
  let a1 = watterson_a1(sample_size)
  let theta = if a1 > 0.0 { segregating_sites.to_double() / a1 } else { 0.0 }
  
  // D = (singletons/n - theta) / SE
  let observed = singletons.to_double() / n
  let expected = theta
  let se = fu_li_d_standard_error(sample_size, segregating_sites, theta)
  
  let d_value = if se > 0.0 {
    (observed - expected) / se
  } else {
    0.0
  }
  
  let p_val = approximate_tajima_d_pvalue(d_value, sample_size)
  
  let interpretation = if d_value < -1.0 {
    "Negative D: excess of singletons, suggests population expansion or positive selection"
  } else if d_value > 1.0 {
    "Positive D: deficiency of singletons, suggests balancing selection or population contraction"
  } else {
    "Neutral pattern"
  }
  
  NeutralityTestResult::{
    test_name: "Fu & Li's D",
    statistic: d_value,
    p_value: p_val,
    is_significant: if p_val != None { p_val.unwrap() < 0.05 } else { false },
    interpretation,
    standard_error: se,
    sample_size,
  }
}

///|
/// Standard error for Fu & Li's D.
pub fn fu_li_d_standard_error(n : Int, s : Int, theta : Double) -> Double {
  if n <= 1 || s <= 0 || theta <= 0.0 {
    return 0.0
  }
  
  let a1 = watterson_a1(n)
  let a2 = watterson_a2(n)
  
  // Simplified variance calculation
  let theta_var = theta * (1.0 / a1 + a2 / (a1 * a1 * n.to_double()))
  theta_var.sqrt()
}

///|
/// Calculate Fu & Li's F statistic.
///
/// F = (singletons - S * a1) / standard_error
///
/// This uses a different weighting of the frequency spectrum.
pub fn fu_li_f(
  singletons : Int,
  segregating_sites : Int,
  sample_size : Int,
) -> NeutralityTestResult {
  if sample_size <= 1 {
    return NeutralityTestResult::{
      test_name: "Fu & Li's F",
      statistic: 0.0,
      p_value: None,
      is_significant: false,
      interpretation: "Insufficient sample size",
      standard_error: 0.0,
      sample_size,
    }
  }
  
  let a1 = watterson_a1(sample_size)
  let a2 = watterson_a2(sample_size)
  
  // Expected singletons under neutrality
  let expected = segregating_sites.to_double() * a1 / (a1 + a2)
  let theta = if a1 > 0.0 { segregating_sites.to_double() / a1 } else { 0.0 }
  
  // Standard error
  let se = (theta * (a1 * a2) / (a1 + a2)).sqrt()
  
  let f_value = if se > 0.0 {
    (singletons.to_double() - expected) / se
  } else {
    0.0
  }
  
  let p_val = approximate_tajima_d_pvalue(f_value, sample_size)
  
  let interpretation = if f_value < -1.0 {
    "Negative F: excess of singletons relative to overall diversity"
  } else if f_value > 1.0 {
    "Positive F: deficiency of singletons relative to overall diversity"
  } else {
    "Neutral pattern"
  }
  
  NeutralityTestResult::{
    test_name: "Fu & Li's F",
    statistic: f_value,
    p_value: p_val,
    is_significant: if p_val != None { p_val.unwrap() < 0.05 } else { false },
    interpretation,
    standard_error: se,
    sample_size,
  }
}

// ============================================================================
// McDonald-Kreitman Test
// ============================================================================

///|
/// Perform the McDonald-Kreitman test.
///
/// Compares nonsynonymous and synonymous:
/// - Polymorphisms (within species): p_nonsyn, p_syn
/// - Fixations (between species): d_nonsyn, d_syn
///
/// Under neutrality: p_nonsyn/p_syn ≈ d_nonsyn/d_syn
///
/// If p_nonsyn/p_syn > d_nonsyn/d_syn:
///   - Positive selection (excess of replacement polymorphisms)
///
/// If p_nonsyn/p_syn < d_nonsyn/d_syn:
///   - Purifying selection (excess of replacement fixations)
///
/// Reference: McDonald, J. H. & Kreitman, M. (1991).
/// "Adaptive protein evolution at the Adh locus in Drosophila."
/// Nature 351:652-654.
pub fn mcdonald_kreitman_test(
  p_nonsyn : Int,
  p_syn : Int,
  d_nonsyn : Int,
  d_syn : Int,
) -> MKTestResult {
  // Calculate ratios
  let poly_ratio = if p_syn > 0 {
    p_nonsyn.to_double() / p_syn.to_double()
  } else {
    0.0
  }
  
  let div_ratio = if d_syn > 0 {
    d_nonsyn.to_double() / d_syn.to_double()
  } else {
    0.0
  }
  
  // Neutrality index
  let ni = if div_ratio > 0.0 {
    poly_ratio / div_ratio
  } else {
    1.0
  }
  
  // Fisher's exact test p-value
  let p_val = fisher_exact_test_2x2(p_nonsyn, p_syn, d_nonsyn, d_syn)
  
  // Interpretation
  let (interpretation, has_selection) = if p_val < 0.05 {
    if poly_ratio > div_ratio * 1.5 {
      ("Significant positive selection: excess of replacement polymorphisms over fixations", true)
    } else if div_ratio > poly_ratio * 1.5 {
      ("Significant purifying selection: excess of replacement fixations over polymorphisms", true)
    } else {
      ("Significant departure from neutrality (non-neutral pattern)", true)
    }
  } else {
    ("No significant departure from neutrality", false)
  }
  
  MKTestResult::{
    p_nonsyn,
    p_syn,
    d_nonsyn,
    d_syn,
    neutrality_index: ni,
    p_value: p_val,
    interpretation,
    has_selection,
  }
}

///|
/// Simple Fisher's exact test for 2x2 contingency table.
pub fn fisher_exact_test_2x2(a : Int, b : Int, c : Int, d : Int) -> Double {
  // Use chi-square approximation for simplicity
  // Better would be to use the exact hypergeometric distribution
  let n = (a + b + c + d).to_double()
  if n <= 0.0 {
    return 1.0
  }
  
  // Expected values
  let row1 = (a + b).to_double()
  let row2 = (c + d).to_double()
  let col1 = (a + c).to_double()
  let col2 = (b + d).to_double()
  
  let e_a = row1 * col1 / n
  let e_b = row1 * col2 / n
  let e_c = row2 * col1 / n
  let e_d = row2 * col2 / n
  
  // Chi-square statistic
  let chi2_a = if e_a > 0.0 {
    let diff = (a.to_double() - e_a).abs() - 0.5
    diff * diff / e_a
  } else { 0.0 }
  let chi2_b = if e_b > 0.0 {
    let diff = (b.to_double() - e_b).abs() - 0.5
    diff * diff / e_b
  } else { 0.0 }
  let chi2_c = if e_c > 0.0 {
    let diff = (c.to_double() - e_c).abs() - 0.5
    diff * diff / e_c
  } else { 0.0 }
  let chi2_d = if e_d > 0.0 {
    let diff = (d.to_double() - e_d).abs() - 0.5
    diff * diff / e_d
  } else { 0.0 }
  let chi2 = chi2_a + chi2_b + chi2_c + chi2_d
  
  // Approximate p-value from chi-square with df=1
  // Using the relationship between chi-square and normal
  let p = @math.exp(-chi2 / 2.0)
  
  if p < 0.0001 { 0.0001 } else if p > 1.0 { 1.0 } else { p }
}

///|
/// McDonald-Kreitman test from site arrays.
pub fn mktest_from_sites(
  poly_sites : Array[PolymorphicSite],
  fixed_sites : Array[PolymorphicSite],
) -> MKTestResult {
  let mut p_nonsyn = 0
  let mut p_syn = 0
  let mut d_nonsyn = 0
  let mut d_syn = 0
  
  for site in poly_sites {
    if site.is_replacement {
      p_nonsyn += 1
    } else {
      p_syn += 1
    }
  }
  
  for site in fixed_sites {
    if site.is_replacement {
      d_nonsyn += 1
    } else {
      d_syn += 1
    }
  }
  
  mcdonald_kreitman_test(p_nonsyn, p_syn, d_nonsyn, d_syn)
}

// ============================================================================
// Allele Frequency Spectrum
// ============================================================================

///|
/// Calculate the allele frequency spectrum from polymorphic sites.
pub fn calculate_afs(sites : Array[PolymorphicSite], sample_size : Int) -> AlleleFrequencySpectrum {
  let mut singletons = 0
  let mut doubletons = 0
  let mut segregating_sites = 0
  let mut total_mutations = 0
  
  let unfolded = Array::make(sample_size, 0)
  let folded = Array::make((sample_size + 1) / 2, 0)
  
  for site in sites {
    segregating_sites += 1
    total_mutations += site.derived_count
    
    if site.derived_count == 1 {
      singletons += 1
    } else if site.derived_count == 2 {
      doubletons += 1
    }
    
    // Unfolded spectrum
    if site.derived_count < sample_size {
      unfolded[site.derived_count] += 1
    }
    
    // Folded spectrum (use minimum of derived and ancestral)
    let minor = if site.derived_count < sample_size - site.derived_count {
      site.derived_count
    } else {
      sample_size - site.derived_count
    }
    if minor > 0 && minor < (sample_size + 1) / 2 {
      folded[minor] += 1
    }
  }
  
  AlleleFrequencySpectrum::{
    sample_size,
    singletons,
    doubletons,
    segregating_sites,
    total_mutations,
    folded_spectrum: folded,
    unfolded_spectrum: unfolded,
  }
}

///|
/// Calculate Tajima's D from allele frequency spectrum.
pub fn tajima_d_from_afs(afs : AlleleFrequencySpectrum, sequence_length : Int) -> NeutralityTestResult {
  let n = afs.sample_size
  if n <= 1 || afs.segregating_sites <= 0 {
    return NeutralityTestResult::{
      test_name: "Tajima's D",
      statistic: 0.0,
      p_value: None,
      is_significant: false,
      interpretation: "Insufficient data",
      standard_error: 0.0,
      sample_size: n,
    }
  }
  
  // Calculate theta_pi (average number of pairwise differences per site)
  let mut theta_pi = 0.0
  let n_float = n.to_double()
  
  for i = 1; i < n; i = i + 1 {
    let count = if i < afs.unfolded_spectrum.length() { afs.unfolded_spectrum[i] } else { 0 }
    let contrib = (2.0 * i.to_double() * (n_float - i.to_double())) / (n_float * (n_float - 1.0))
    theta_pi += contrib * count.to_double()
  }
  
  // Normalize by sequence length
  if sequence_length > 0 {
    theta_pi = theta_pi / sequence_length.to_double()
  }
  
  popgen_tajima_d(theta_pi, afs.segregating_sites, n, sequence_length)
}

///|
/// Calculate Fu & Li's D from allele frequency spectrum.
pub fn fu_li_d_from_afs(afs : AlleleFrequencySpectrum) -> NeutralityTestResult {
  fu_li_d(afs.singletons, afs.segregating_sites, afs.sample_size)
}

// ============================================================================
// Combined Neutrality Analysis
// ============================================================================

///|
/// Run a comprehensive set of neutrality tests.
pub struct NeutralityAnalysisResult {
  tajima_d : NeutralityTestResult
  fu_li_d : NeutralityTestResult
  fu_li_f : NeutralityTestResult
  watterson_theta : Double
  tajima_theta_pi : Double
  segregating_sites : Int
  sample_size : Int
  is_neutral : Bool
}

///|
/// Run comprehensive neutrality analysis.
pub fn run_neutrality_analysis(
  sites : Array[PolymorphicSite],
  sample_size : Int,
  sequence_length : Int,
) -> NeutralityAnalysisResult {
  let afs = calculate_afs(sites, sample_size)
  
  let theta_w = popgen_watterson_theta(afs.segregating_sites, sample_size)
  let result_d = tajima_d_from_afs(afs, sequence_length)
  let result_fu_d = fu_li_d_from_afs(afs)
  let result_fu_f = fu_li_f(afs.singletons, afs.segregating_sites, sample_size)
  
  // Calculate theta_pi
  let mut theta_pi = 0.0
  let n_float = sample_size.to_double()
  
  for i = 1; i < sample_size; i = i + 1 {
    let count = if i < afs.unfolded_spectrum.length() { afs.unfolded_spectrum[i] } else { 0 }
    let contrib = (2.0 * i.to_double() * (n_float - i.to_double())) / (n_float * (n_float - 1.0))
    theta_pi += contrib * count.to_double()
  }
  if sequence_length > 0 {
    theta_pi = theta_pi / sequence_length.to_double()
  }
  
  // Determine if neutral
  let is_neutral = !result_d.is_significant && 
    !result_fu_d.is_significant && 
    !result_fu_f.is_significant
  
  NeutralityAnalysisResult::{
    tajima_d: result_d,
    fu_li_d: result_fu_d,
    fu_li_f: result_fu_f,
    watterson_theta: theta_w,
    tajima_theta_pi: theta_pi,
    segregating_sites: afs.segregating_sites,
    sample_size,
    is_neutral,
  }
}

///|
/// Format a neutrality analysis result as a string.
pub fn neutrality_analysis_to_string(result : NeutralityAnalysisResult) -> String {
  let mut s = ""
  s = s + "=== Neutrality Analysis Results ===\n\n"
  s = s + "Sample size: "
  s = s + result.sample_size.to_string()
  s = s + "\n"
  s = s + "Segregating sites: "
  s = s + result.segregating_sites.to_string()
  s = s + "\n\n"
  s = s + "Watterson's theta (theta_W): "
  s = s + result.watterson_theta.to_string()
  s = s + "\n"
  s = s + "Tajima's theta (theta_pi): "
  s = s + result.tajima_theta_pi.to_string()
  s = s + "\n\n"
  
  // Tajima's D
  s = s + "Tajima's D: "
  s = s + result.tajima_d.statistic.to_string()
  if result.tajima_d.p_value != None {
    s = s + " (p = "
    s = s + result.tajima_d.p_value.unwrap().to_string()
    s = s + ")"
  }
  if result.tajima_d.is_significant {
    s = s + " *** SIGNIFICANT ***"
  }
  s = s + "\n  "
  s = s + result.tajima_d.interpretation
  s = s + "\n\n"
  
  // Fu & Li's D
  s = s + "Fu & Li's D: "
  s = s + result.fu_li_d.statistic.to_string()
  if result.fu_li_d.p_value != None {
    s = s + " (p = "
    s = s + result.fu_li_d.p_value.unwrap().to_string()
    s = s + ")"
  }
  if result.fu_li_d.is_significant {
    s = s + " *** SIGNIFICANT ***"
  }
  s = s + "\n  "
  s = s + result.fu_li_d.interpretation
  s = s + "\n\n"
  
  // Fu & Li's F
  s = s + "Fu & Li's F: "
  s = s + result.fu_li_f.statistic.to_string()
  if result.fu_li_f.p_value != None {
    s = s + " (p = "
    s = s + result.fu_li_f.p_value.unwrap().to_string()
    s = s + ")"
  }
  if result.fu_li_f.is_significant {
    s = s + " *** SIGNIFICANT ***"
  }
  s = s + "\n  "
  s = s + result.fu_li_f.interpretation
  s = s + "\n\n"
  
  s = s + "Overall: "
  if result.is_neutral {
    s = s + "CONSISTENT WITH NEUTRALITY"
  } else {
    s = s + "SIGNIFICANT DEPARTURE FROM NEUTRALITY DETECTED"
  }
  s = s + "\n"
  
  s
}