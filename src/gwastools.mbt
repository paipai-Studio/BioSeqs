///|
/// GWASTools (Genome-Wide Association Studies) quality-control and association
/// core: minor allele frequency, missing call rate, Hardy-Weinberg equilibrium
/// exact test, and the allele-based chi-square association test.

///|
pub struct GwasAssocResult {
  chi2 : Double
  p_value : Double
}

///|
/// Minor allele frequency from genotype counts (hom_ref/het/hom_alt). Returns
/// the smaller of the two allele frequencies (NaN for a zero-total locus).
pub fn gwas_minor_allele_frequency(hom_ref : Int, het : Int, hom_alt : Int) -> Double {
  let n = hom_ref + het + hom_alt
  if n == 0 {
    return 0.0 / 0.0
  }
  let p_ref = (2.0 * hom_ref.to_double() + het.to_double()) / (2.0 * n.to_double())
  let p_alt = 1.0 - p_ref
  if p_ref < p_alt {
    p_ref
  } else {
    p_alt
  }
}

///|
/// Missing call rate given the number of called genotypes vs total:
/// `(total - called) / total` (NaN for a zero-total input).
pub fn gwas_missing_call_rate(called : Int, total : Int) -> Double {
  if total <= 0 {
    return 0.0 / 0.0
  }
  (total - called).to_double() / total.to_double()
}

///|
/// Natural log of n! (n >= 0), exact for the moderate n seen in exact HWE
/// testing; used to keep the HWE test numerically stable.
fn ln_factorial(n : Int) -> Double {
  let mut s = 0.0
  let mut i = 2
  while i <= n {
    s = s + @math.ln(i.to_double())
    i = i + 1
  }
  s
}

///|
/// Log-probability of a genotype configuration under HWE given allele counts:
/// `ln[ n! (n_A!)(n_B!) 2^het / ((aa!)(het!)(bb!)(2n)!) ]`.
fn hwe_lnp(n : Int, n_a : Int, n_b : Int, aa : Int, het : Int, bb : Int, ln2 : Double) -> Double {
  ln_factorial(n) + ln_factorial(n_a) + ln_factorial(n_b) +
  het.to_double() * ln2 -
  ln_factorial(aa) -
  ln_factorial(het) -
  ln_factorial(bb) -
  ln_factorial(2 * n)
}

///|
/// Hardy-Weinberg equilibrium exact test (GWASTools `exactHWE`): conditional on
/// the observed allele counts, sum the probabilities of every genotype
/// configuration no more likely than the observed one. Returns the two-sided
/// exact p-value.
pub fn gwas_hwe_exact(hom_ref : Int, het : Int, hom_alt : Int) -> Double {
  let n = hom_ref + het + hom_alt
  if n == 0 {
    return 0.0 / 0.0
  }
  let n_a = 2 * hom_ref + het
  let n_b = 2 * n - n_a
  let ln2 = @math.ln(2.0)
  let obs_lnp = hwe_lnp(n, n_a, n_b, hom_ref, het, hom_alt, ln2)
  let min_aa = if n_a - n < 0 { 0 } else { n_a - n }
  let max_aa = n_a / 2
  let mut p = 0.0
  let mut aa = min_aa
  while aa <= max_aa {
    let hh = n_a - 2 * aa
    let bb = n - aa - hh
    if hh >= 0 && bb >= 0 {
      let lnp = hwe_lnp(n, n_a, n_b, aa, hh, bb, ln2)
      if lnp <= obs_lnp + 1.0e-12 {
        p = p + @math.exp(lnp)
      }
    }
    aa = aa + 1
  }
  if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Allele-based case/control association chi-square test (GWASTools
/// `assocTest`): a 2×2 allele contingency table. Returns the 1-df chi-square
/// statistic and its two-sided asymptotic p-value (normal approximation at
/// `sqrt(chi2)`).
pub fn gwas_chisq_assoc(
  case_ref : Int,
  case_alt : Int,
  ctrl_ref : Int,
  ctrl_alt : Int,
) -> GwasAssocResult {
  let a = case_ref.to_double()
  let c = case_alt.to_double()
  let b = ctrl_ref.to_double()
  let d = ctrl_alt.to_double()
  let n = a + b + c + d
  if n == 0.0 {
    return GwasAssocResult::{ chi2: 0.0, p_value: 1.0 }
  }
  let denom = (a + b) * (c + d) * (a + c) * (b + d)
  let chi2 = if denom == 0.0 {
    0.0
  } else {
    n * (a * d - b * c) * (a * d - b * c) / denom
  }
  GwasAssocResult::{ chi2, p_value: stat_pvalue_zscore(chi2.sqrt()) }
}