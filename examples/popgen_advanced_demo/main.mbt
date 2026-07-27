///|
/// Advanced Population Genetics Demo.
/// Demonstrates neutrality tests, allele frequency spectrum analysis,
/// and the McDonald-Kreitman test.

fn main {
  println("=== Advanced Population Genetics Demo ===\n")

  println("1. Creating polymorphic sites...")
  let sites : Array[@src.PolymorphicSite] = Array::new()
  let sample_size = 30
  let mut i = 0
  while i < 10 {
    sites.push(@src.new_polymorphic_site(i * 100 + 1, 1, sample_size, i % 3 == 0))
    i = i + 1
  }
  let mut j = 0
  while j < 5 {
    sites.push(@src.new_polymorphic_site(1000 + j * 100, 2, sample_size, false))
    j = j + 1
  }
  let mut k = 0
  while k < 3 {
    sites.push(@src.new_polymorphic_site(2000 + k * 100, 5, sample_size, k == 1))
    k = k + 1
  }
  println("   Created \{sites.length()} polymorphic sites (sample size = \{sample_size})")

  println("\n2. Computing Watterson's theta...")
  let theta_w = @src.popgen_watterson_theta(sites.length(), sample_size)
  println("   Theta_W = S / a1 = \{theta_w}")
  let a1 = @src.watterson_a1(sample_size)
  println("   a1 = sum(1/i) for i=1..n-1 = \{a1}")

  println("\n3. Computing Tajima's D...")
  let sequence_length = 5000
  let result_d = @src.popgen_tajima_d(0.02, sites.length(), sample_size, sequence_length)
  println("   Tajima's D = \{result_d.statistic}")
  println("   Standard error = \{result_d.standard_error}")
  if result_d.p_value != None {
    println("   p-value = \{result_d.p_value.unwrap()}")
  }
  println("   Significant (p < 0.05): \{result_d.is_significant}")
  println("   Interpretation: \{result_d.interpretation}")

  println("\n4. Computing Fu & Li's D...")
  let singletons = 10
  let result_fu_d = @src.fu_li_d(singletons, sites.length(), sample_size)
  println("   Fu & Li's D = \{result_fu_d.statistic}")
  println("   Standard error = \{result_fu_d.standard_error}")
  if result_fu_d.p_value != None {
    println("   p-value = \{result_fu_d.p_value.unwrap()}")
  }
  println("   Significant (p < 0.05): \{result_fu_d.is_significant}")
  println("   Interpretation: \{result_fu_d.interpretation}")

  println("\n5. Computing Fu & Li's F...")
  let result_fu_f = @src.fu_li_f(singletons, sites.length(), sample_size)
  println("   Fu & Li's F = \{result_fu_f.statistic}")
  println("   Standard error = \{result_fu_f.standard_error}")
  if result_fu_f.p_value != None {
    println("   p-value = \{result_fu_f.p_value.unwrap()}")
  }
  println("   Significant (p < 0.05): \{result_fu_f.is_significant}")
  println("   Interpretation: \{result_fu_f.interpretation}")

  println("\n6. McDonald-Kreitman test...")
  let mk_result = @src.mcdonald_kreitman_test(25, 35, 8, 42)
  println("   Polymorphisms: p_nonsyn = \{mk_result.p_nonsyn}, p_syn = \{mk_result.p_syn}")
  println("   Fixations: d_nonsyn = \{mk_result.d_nonsyn}, d_syn = \{mk_result.d_syn}")
  println("   Neutrality index (NI) = \{mk_result.neutrality_index}")
  println("   p-value = \{mk_result.p_value}")
  println("   Has selection: \{mk_result.has_selection}")
  println("   Interpretation: \{mk_result.interpretation}")

  println("\n7. Allele frequency spectrum...")
  let afs = @src.calculate_afs(sites, sample_size)
  println("   Sample size: \{afs.sample_size}")
  println("   Segregating sites: \{afs.segregating_sites}")
  println("   Singletons: \{afs.singletons}")
  println("   Doubletons: \{afs.doubletons}")
  println("   Total mutations: \{afs.total_mutations}")

  println("\n8. Running comprehensive neutrality analysis...")
  let result = @src.run_neutrality_analysis(sites, sample_size, sequence_length)
  let report = @src.neutrality_analysis_to_string(result)
  println("\{report}")

  println("=== Advanced Population Genetics Demo Complete ===")
}