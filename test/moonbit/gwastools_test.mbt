///|
/// Tests for GWASTools QC + association core.

///|
fn gwas_dclose(a : Double, b : Double) -> Bool {
  let d = a - b
  let ad = if d < 0.0 { 0.0 - d } else { d }
  ad < 1.0e-4
}

///|
test "gwas_minor_allele_frequency" {
  assert_eq(@src.gwas_minor_allele_frequency(10, 0, 0), 0.0)
  assert_eq(@src.gwas_minor_allele_frequency(0, 10, 0), 0.5)
  assert_eq(@src.gwas_minor_allele_frequency(5, 5, 0), 0.25)
}

///|
test "gwas_missing_call_rate" {
  assert_eq(@src.gwas_missing_call_rate(150, 200), 0.25)
  assert_eq(@src.gwas_missing_call_rate(100, 100), 0.0)
}

///|
test "gwas_hwe_exact_symmetric_heterozygotes" {
  // n=3, all heterozygotes: conditional exact p = 0.4
  let p = @src.gwas_hwe_exact(0, 3, 0)
  assert_eq(gwas_dclose(p, 0.4), true)
}

///|
test "gwas_hwe_exact_most_likely_configuration" {
  // p(1,2,1) is the most likely configuration -> exact p = 1.0
  let p = @src.gwas_hwe_exact(1, 2, 1)
  assert_eq(gwas_dclose(p, 1.0), true)
}

///|
test "gwas_chisq_assoc" {
  // 2x2 allele table: case ref/alt = 100/50, ctrl ref/alt = 80/70
  // chi2 = 50/9 ~ 5.5556, p ~ 0.0184
  let r = @src.gwas_chisq_assoc(100, 50, 80, 70)
  assert_eq(gwas_dclose(r.chi2, 50.0 / 9.0), true)
  assert_eq(gwas_dclose(r.p_value, 0.01842), true)
}