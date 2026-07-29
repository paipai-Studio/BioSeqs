///|
/// Tests for Bio.PopGen Advanced - Advanced population genetics statistics.

///|
/// Test PolymorphicSite creation.
test "poly_site_create" {
  let site = @src.new_polymorphic_site(100, 3, 10, false)
  assert_eq(site.position, 100)
  assert_eq(site.derived_count, 3)
  assert_eq(site.sample_size, 10)
  assert_eq(site.is_replacement, false)
  assert_eq(site.is_singleton, false)
}

///|
/// Test singleton site.
test "poly_site_singleton" {
  let site = @src.new_singleton(100, 20, true)
  assert_eq(site.is_singleton, true)
  assert_eq(site.derived_count, 1)
  assert_eq(site.is_replacement, true)
  assert_eq(site.derived_freq, 0.05)
}

///|
/// Test doubleton site.
test "poly_site_doubleton" {
  let site = @src.new_doubleton(200, 20, false)
  assert_eq(site.is_singleton, false)
  assert_eq(site.derived_count, 2)
  assert_eq(site.derived_freq, 0.1)
}

///|
/// Test Watterson's theta calculation.
test "watterson_theta_basic" {
  let theta = @src.popgen_watterson_theta(5, 10)
  // With 5 segregating sites and n=10, theta should be positive
  assert_true(theta > 0.0)
}

///|
/// Test Watterson's theta with zero sites.
test "watterson_theta_zero" {
  let theta = @src.popgen_watterson_theta(0, 10)
  assert_eq(theta, 0.0)
}

///|
/// Test Watterson's theta with small sample.
test "watterson_theta_small_sample" {
  let theta = @src.popgen_watterson_theta(5, 1)
  assert_eq(theta, 0.0)
}

///|
/// Test watterson_a1 calculation.
test "watterson_a1" {
  let a1 = @src.watterson_a1(10)
  // a1 = sum(1/i) for i=1..9 = 1 + 1/2 + 1/3 + ... + 1/9
  let expected = 1.0 + 0.5 + 0.3333333 + 0.25 + 0.2 + 0.1666667 + 0.1428571 + 0.125 + 0.1111111
  assert_true((a1 - expected).abs() < 0.01)
}

///|
/// Test watterson_a1 with n=1.
test "watterson_a1_one" {
  let a1 = @src.watterson_a1(1)
  assert_eq(a1, 0.0)
}

///|
/// Test watterson_a2 calculation.
test "watterson_a2" {
  let a2 = @src.watterson_a2(10)
  // a2 = sum(1/i^2) for i=1..9
  assert_true(a2 > 0.0)
}

///|
/// Test Tajima's D calculation.
test "tajima_d_basic" {
  let result = @src.popgen_tajima_d(0.01, 10, 10, 1000)
  assert_eq(result.test_name, "Tajima's D")
  assert_true(!result.statistic.is_nan())
}

///|
/// Test Tajima's D negative (excess low-frequency variants).
test "tajima_d_negative" {
  // Create a scenario where theta_pi < theta_W
  // This happens when there are many low-frequency variants
  let result = @src.popgen_tajima_d(0.005, 20, 20, 1000)
  // With many singletons (high S relative to pi), D should be negative or close to 0
  assert_true(!result.statistic.is_nan())
}

///|
/// Test Tajima's D p-value.
test "tajima_d_pvalue" {
  let result = @src.popgen_tajima_d(0.01, 10, 10, 1000)
  if result.p_value != None {
    assert_true(result.p_value.unwrap() >= 0.0)
    assert_true(result.p_value.unwrap() <= 1.0)
  }
}

///|
/// Test Fu & Li's D calculation.
test "fu_li_d_basic" {
  let result = @src.fu_li_d(5, 10, 10)
  assert_eq(result.test_name, "Fu & Li's D")
  assert_true(!result.statistic.is_nan())
}

///|
/// Test Fu & Li's D with insufficient data.
test "fu_li_d_small_sample" {
  let result = @src.fu_li_d(1, 1, 1)
  assert_eq(result.statistic, 0.0)
}

///|
/// Test Fu & Li's F calculation.
test "fu_li_f_basic" {
  let result = @src.fu_li_f(5, 10, 10)
  assert_eq(result.test_name, "Fu & Li's F")
  assert_true(!result.statistic.is_nan())
}

///|
/// Test Fu & Li's F with insufficient data.
test "fu_li_f_small_sample" {
  let result = @src.fu_li_f(1, 1, 1)
  assert_eq(result.statistic, 0.0)
}

///|
/// Test McDonald-Kreitman test.
test "mktest_basic" {
  let result = @src.mcdonald_kreitman_test(20, 30, 10, 40)
  assert_eq(result.p_nonsyn, 20)
  assert_eq(result.p_syn, 30)
  assert_eq(result.d_nonsyn, 10)
  assert_eq(result.d_syn, 40)
  assert_true(result.neutrality_index > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

///|
/// Test MK test with purifying selection.
test "mktest_purifying_selection" {
  // More nonsynonymous fixations than polymorphisms suggests purifying selection
  let result = @src.mcdonald_kreitman_test(5, 40, 20, 30)
  assert_true(result.neutrality_index < 1.0)
}

///|
/// Test MK test with positive selection.
test "mktest_positive_selection" {
  // More nonsynonymous polymorphisms than fixations suggests positive selection
  let result = @src.mcdonald_kreitman_test(30, 20, 5, 40)
  assert_true(result.neutrality_index > 1.0)
}

///|
/// Test MK test with zero values.
test "mktest_zero_values" {
  let result = @src.mcdonald_kreitman_test(0, 0, 0, 0)
  assert_eq(result.neutrality_index, 1.0)
  assert_eq(result.p_value, 1.0)
}

///|
/// Test MK test from site arrays.
test "mktest_from_sites" {
  let poly_sites : Array[@src.PolymorphicSite] = Array::new()
  let fixed_sites : Array[@src.PolymorphicSite] = Array::new()
  
  // Add polymorphic sites
  poly_sites.push(@src.new_polymorphic_site(100, 2, 20, true))  // nonsynonymous
  poly_sites.push(@src.new_polymorphic_site(200, 3, 20, false)) // synonymous
  poly_sites.push(@src.new_polymorphic_site(300, 1, 20, true))  // nonsynonymous
  
  // Add fixed sites
  fixed_sites.push(@src.new_polymorphic_site(400, 20, 20, true))  // nonsynonymous
  fixed_sites.push(@src.new_polymorphic_site(500, 20, 20, false)) // synonymous
  
  let result = @src.mktest_from_sites(poly_sites, fixed_sites)
  assert_eq(result.p_nonsyn, 2)
  assert_eq(result.p_syn, 1)
  assert_eq(result.d_nonsyn, 1)
  assert_eq(result.d_syn, 1)
}

///|
/// Test allele frequency spectrum calculation.
test "afs_calculate" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  sites.push(@src.new_singleton(100, 20, false))
  sites.push(@src.new_singleton(200, 20, true))
  sites.push(@src.new_doubleton(300, 20, false))
  sites.push(@src.new_polymorphic_site(400, 5, 20, false))
  
  let afs = @src.calculate_afs(sites, 20)
  assert_eq(afs.singletons, 2)
  assert_eq(afs.doubletons, 1)
  assert_eq(afs.segregating_sites, 4)
}

///|
/// Test allele frequency spectrum with empty sites.
test "afs_empty" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  let afs = @src.calculate_afs(sites, 10)
  assert_eq(afs.segregating_sites, 0)
  assert_eq(afs.singletons, 0)
}

///|
/// Test Tajima's D from AFS.
test "tajima_d_from_afs" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  sites.push(@src.new_singleton(100, 20, false))
  sites.push(@src.new_singleton(200, 20, false))
  sites.push(@src.new_doubleton(300, 20, false))
  
  let afs = @src.calculate_afs(sites, 20)
  let result = @src.tajima_d_from_afs(afs, 1000)
  assert_eq(result.test_name, "Tajima's D")
  assert_true(!result.statistic.is_nan())
}

///|
/// Test Fu & Li's D from AFS.
test "fu_li_d_from_afs" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  sites.push(@src.new_singleton(100, 20, false))
  sites.push(@src.new_doubleton(200, 20, false))
  
  let afs = @src.calculate_afs(sites, 20)
  let result = @src.fu_li_d_from_afs(afs)
  assert_eq(result.test_name, "Fu & Li's D")
}

///|
/// Test normal CDF approximation.
test "normal_cdf" {
  let p0 = @src.popgen_normal_cdf(0.0)
  assert_true((p0 - 0.5).abs() < 0.01)
  
  let p1 = @src.popgen_normal_cdf(1.96)
  assert_true(p1 > 0.96 && p1 < 0.98)
  
  let p_neg1 = @src.popgen_normal_cdf(-1.96)
  assert_true(p_neg1 > 0.02 && p_neg1 < 0.04)
}

///|
/// Test comprehensive neutrality analysis.
test "neutrality_analysis" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  // Create a mix of sites
  for i = 1; i < 10; i = i + 1 {
    sites.push(@src.new_singleton(i * 100, 30, i % 2 == 0))
  }
  for i = 1; i < 5; i = i + 1 {
    sites.push(@src.new_doubleton(500 + i * 100, 30, false))
  }
  for i = 1; i < 3; i = i + 1 {
    sites.push(@src.new_polymorphic_site(1000 + i * 100, 5, 30, i == 1))
  }
  
  let result = @src.run_neutrality_analysis(sites, 30, 10000)
  assert_true(!result.tajima_d.statistic.is_nan())
  assert_true(!result.fu_li_d.statistic.is_nan())
  assert_true(!result.fu_li_f.statistic.is_nan())
  assert_true(result.watterson_theta >= 0.0)
  assert_true(result.tajima_theta_pi >= 0.0)
  assert_eq(result.segregating_sites, 15)
  assert_eq(result.sample_size, 30)
}

///|
/// Test neutrality analysis with empty data.
test "neutrality_analysis_empty" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  let result = @src.run_neutrality_analysis(sites, 10, 1000)
  assert_eq(result.segregating_sites, 0)
}

///|
/// Test neutrality analysis output formatting.
test "neutrality_analysis_to_string" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  sites.push(@src.new_singleton(100, 20, false))
  sites.push(@src.new_singleton(200, 20, true))
  sites.push(@src.new_doubleton(300, 20, false))
  
  let result = @src.run_neutrality_analysis(sites, 20, 5000)
  let output = @src.neutrality_analysis_to_string(result)
  assert_true(output.contains("Neutrality Analysis Results"))
  assert_true(output.contains("Sample size:"))
  assert_true(output.contains("Segregating sites:"))
  assert_true(output.contains("Tajima's D:"))
  assert_true(output.contains("Fu & Li's D:"))
  assert_true(output.contains("Fu & Li's F:"))
}

///|
/// Test fisher exact test.
test "fisher_exact_test" {
  let p = @src.fisher_exact_test_2x2(10, 20, 5, 30)
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
/// Test fisher exact test with empty table.
test "fisher_exact_test_empty" {
  let p = @src.fisher_exact_test_2x2(0, 0, 0, 0)
  assert_eq(p, 1.0)
}

///|
/// Test NeutralityTestConfig.
test "neutrality_config" {
  let config = @src.NeutralityTestConfig::default(20, 10000)
  assert_eq(config.sample_size, 20)
  assert_eq(config.sequence_length, 10000)
  assert_eq(config.mutation_rate, 0.0)
  assert_false(config.use_folded)
}

///|
/// Test Tajima's B1 calculation.
test "tajima_b1" {
  let b1 = @src.tajima_b1(10)
  // b1 = (10+1)/(3*(10-1)) = 11/27 ≈ 0.4074
  assert_true((b1 - 0.4074).abs() < 0.01)
}

///|
/// Test Tajima's B2 calculation.
test "tajima_b2" {
  let b2 = @src.tajima_b2(10)
  assert_true(!b2.is_nan())
}

///|
/// Test multiple singletons detection.
test "multiple_singletons" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  for i = 1; i < 15; i = i + 1 {
    sites.push(@src.new_singleton(i * 50, 25, i % 3 == 0))
  }
  
  let afs = @src.calculate_afs(sites, 25)
  assert_eq(afs.singletons, 14)
  assert_eq(afs.segregating_sites, 14)
}

///|
/// Test high frequency variants.
test "high_freq_variants" {
  let sites : Array[@src.PolymorphicSite] = Array::new()
  // Create high-frequency variants (derived count close to sample size)
  sites.push(@src.new_polymorphic_site(100, 24, 25, false))
  sites.push(@src.new_polymorphic_site(200, 23, 25, true))
  sites.push(@src.new_polymorphic_site(300, 22, 25, false))
  
  let afs = @src.calculate_afs(sites, 25)
  assert_eq(afs.singletons, 0) // No singletons
}