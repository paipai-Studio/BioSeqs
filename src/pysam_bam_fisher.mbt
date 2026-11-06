/// pysam_bam_fisher.mbt — Fisher's exact test for strand bias
/// Inspired by pysam/samtools fisher

///|
/// Fisher's exact test for 2x2 contingency table.
/// Computes the probability of the observed arrangement.
pub fn pysam_fisher_exact(a : Int32, b : Int32, c : Int32, d : Int32) -> Double {
  // 2x2 table:
  // | a  b |
  // | c  d |
  let n = a + b + c + d
  if n == 0 {
    return 1.0
  }
  // Log factorial approach to avoid overflow
  let log_p = pysam_lnfact(a + b) +
    pysam_lnfact(c + d) +
    pysam_lnfact(a + c) +
    pysam_lnfact(b + d) -
    pysam_lnfact(a) -
    pysam_lnfact(b) -
    pysam_lnfact(c) -
    pysam_lnfact(d) -
    pysam_lnfact(n)
  @math.exp(log_p)
}

///|
/// Log factorial using Stirling's approximation for large values.
fn pysam_lnfact(n : Int32) -> Double {
  if n <= 1 {
    return 0.0
  }
  if n <= 100 {
    // Exact for small values
    let mut result = 0.0
    for i = 2; i <= n; i = i + 1 {
      result = result + @math.ln(i.to_double())
    }
    return result
  }
  // Stirling: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2*pi*n)
  let nd = n.to_double()
  nd * @math.ln(nd) - nd + 0.5 * @math.ln(6.283185307179586 * nd)
}

///|
/// Strand bias test.
/// Counts: forward ref, reverse ref, forward alt, reverse alt.
/// Returns p-value for strand bias.
pub fn pysam_strand_bias_test(
  fwd_ref : Int32,
  rev_ref : Int32,
  fwd_alt : Int32,
  rev_alt : Int32,
) -> Double {
  // 2x2 table:
  // | fwd_ref  rev_ref |
  // | fwd_alt  rev_alt |
  pysam_fisher_exact(fwd_ref, rev_ref, fwd_alt, rev_alt)
}

///|
/// Read pair orientation bias test.
pub fn pysam_orientation_bias(
  f1f2 : Int32, // read1 forward, read2 forward
  f1r2 : Int32, // read1 forward, read2 reverse
  r1f2 : Int32, // read1 reverse, read2 forward
  r1r2 : Int32, // read1 reverse, read2 reverse
) -> Double {
  pysam_fisher_exact(f1f2, f1r2, r1f2, r1r2)
}

///|
/// Phred-scaled p-value.
pub fn pysam_phred_p(p : Double) -> Double {
  if p <= 0.0 {
    999.0
  } else {
    -10.0 * @math.log10(p)
  }
}
