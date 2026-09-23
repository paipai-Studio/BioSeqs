///|
/// Tests for methylGSA — methylation gene set testing.

///|
test "methylgsa gene min p" {
  let input = @src.methylgsa_sample_input()
  let min_p = @src.methylgsa_gene_min_p(input)
  match min_p.get("GENE_A") {
    Some(p) => assert_true((p - 0.001).abs() < 1.0e-12)
    None => assert_true(false)
  }
  match min_p.get("GENE_D") {
    Some(p) => assert_true((p - 0.002).abs() < 1.0e-12)
    None => assert_true(false)
  }
  match min_p.get("GENE_E") {
    Some(p) => assert_true((p - 0.7).abs() < 1.0e-12)
    None => assert_true(false)
  }
  assert_true(min_p.get("GENE_Z") is None)
}

///|
test "methylrra_ols ranks significant set first" {
  let input = @src.methylgsa_sample_input()
  let results = @src.methylrra_ols(
    input,
    gene_sets=@src.methylgsa_sample_gene_sets(),
  )
  assert_eq(results.length(), 2)
  assert_eq(results[0].set_name, "SIGNIFICANT_SET")
  assert_eq(results[0].set_size, 2)
  assert_true(results[0].pvalue <= results[1].pvalue)
}

///|
test "methylrra_ols qvalues within unit range" {
  let input = @src.methylgsa_sample_input()
  let results = @src.methylrra_ols(
    input,
    gene_sets=@src.methylgsa_sample_gene_sets(),
  )
  for r in results {
    assert_true(r.qvalue >= 0.0)
    assert_true(r.qvalue <= 1.0)
  }
}

///|
test "methylrra_robust ranks significant set first" {
  let input = @src.methylgsa_sample_input()
  let results = @src.methylrra_robust(
    input,
    gene_sets=@src.methylgsa_sample_gene_sets(),
  )
  assert_eq(results.length(), 2)
  let sig = results[0]
  let neu = results[1]
  assert_eq(sig.set_name, "SIGNIFICANT_SET")
  assert_true(sig.pvalue <= neu.pvalue)
  assert_true(sig.statistic <= 1.0)
}

///|
test "methylrra_robust empty set" {
  let input = @src.methylgsa_sample_input()
  let sets : Map[String, Array[String]] = Map([], capacity=2)
  sets.set("EMPTY", ["NO_SUCH_GENE"])
  let results = @src.methylrra_robust(input, gene_sets=sets)
  assert_eq(results.length(), 1)
  assert_eq(results[0].set_size, 0)
  assert_eq(results[0].pvalue, 1.0)
}

///|
test "methylgsa bh correction monotone" {
  let results : Array[@src.MethylgsaResult] = Array(capacity=0)
  results.push(
    @src.methylgsa_result_new(
      set_name="a",
      set_size=1,
      statistic=1.0,
      pvalue=0.01,
      qvalue=1.0,
    ),
  )
  results.push(
    @src.methylgsa_result_new(
      set_name="b",
      set_size=1,
      statistic=0.5,
      pvalue=0.04,
      qvalue=1.0,
    ),
  )
  results.push(
    @src.methylgsa_result_new(
      set_name="c",
      set_size=1,
      statistic=0.2,
      pvalue=0.03,
      qvalue=1.0,
    ),
  )
  let adjusted = @src.methylgsa_bh(results)
  for r in adjusted {
    assert_true(r.qvalue >= r.pvalue)
    assert_true(r.qvalue <= 1.0)
  }
  assert_true(adjusted[0].qvalue <= adjusted[1].qvalue)
}

///|
test "methylgsa log10 sanity" {
  assert_true((@src.methylgsa_log10(10.0) - 1.0).abs() < 1.0e-4)
  assert_true((@src.methylgsa_log10(100.0) - 2.0).abs() < 1.0e-4)
  assert_true((@src.methylgsa_log10(0.001) + 3.0).abs() < 1.0e-4)
}

///|
test "methylgsa ols statistic of empty set is zero" {
  let input = @src.methylgsa_sample_input()
  let min_p = @src.methylgsa_gene_min_p(input)
  let (stat, size) = @src.methylgsa_ols_statistic(min_p, ["MISSING"])
  assert_eq(stat, 0.0)
  assert_eq(size, 0)
}