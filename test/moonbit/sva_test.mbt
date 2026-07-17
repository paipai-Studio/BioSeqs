///|
test "sva_combat_result_new" {
  let adjusted = [[1.0, 2.0], [3.0, 4.0]]
  let gamma_hat = [[0.5, 0.6], [0.7, 0.8]]
  let delta_hat = [[1.0, 1.1], [1.2, 1.3]]
  let result = @src.ComBatResult::new(adjusted, gamma_hat, delta_hat)

  assert_eq(result.adjusted.length(), 2)
  assert_eq(result.gamma_hat.length(), 2)
  assert_eq(result.delta_hat.length(), 2)
}

///|
test "sva_sva_result_new" {
  let sv = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let pprob = [0.5, 0.8, 0.3]
  let result = @src.SvaResult::new(sv, pprob, 2)

  assert_eq(result.sv.length(), 2)
  assert_eq(result.pprob_gam.length(), 3)
  assert_eq(result.n_sv, 2)
}

///|
test "sva_sample_data" {
  let (dat, genes, samples, batch) = @src.sample_expression_matrix()

  assert_eq(genes.length(), 20)
  assert_eq(samples.length(), 12)
  assert_eq(dat.length(), 20)
  assert_eq(dat[0].length(), 12)
  assert_eq(batch.length(), 12)
}

///|
test "sva_combat" {
  let (dat, _genes, _samples, batch) = @src.sample_expression_matrix()

  let result = @src.combat(dat, batch)

  assert_eq(result.adjusted.length(), dat.length())
  assert_eq(result.adjusted[0].length(), dat[0].length())
  assert_eq(result.gamma_hat.length(), dat.length())
  assert_eq(result.delta_hat.length(), dat.length())
}

///|
test "sva_combat_simple" {
  let (dat, _genes, _samples, batch) = @src.sample_expression_matrix()

  let adjusted = @src.combat_simple(dat, batch)

  assert_eq(adjusted.length(), dat.length())
  assert_eq(adjusted[0].length(), dat[0].length())
}

///|
test "sva_sva" {
  let (dat, _genes, _samples, _batch) = @src.sample_expression_matrix()
  let mod = [[1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0]]

  let result = @src.sva(dat, mod, n_sv = 2)

  assert_eq(result.n_sv, 2)
  assert_eq(result.sv.length(), 2)
  assert_eq(result.pprob_gam.length(), dat.length())
}
