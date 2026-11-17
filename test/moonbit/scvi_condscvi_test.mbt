///|
/// Tests for CondSCVI (conditional VAE) and the cells2spatial DestVI bridge.

///|
test "scvi_condscvi_train_and_profiles" {
  // two cell types with complementary gene programs
  let sc = [
    [50.0, 50.0, 50.0, 1.0, 1.0, 1.0],
    [55.0, 45.0, 52.0, 0.0, 1.0, 2.0],
    [48.0, 52.0, 49.0, 1.0, 0.0, 1.0],
    [1.0, 1.0, 2.0, 45.0, 50.0, 48.0],
    [2.0, 0.0, 1.0, 52.0, 48.0, 50.0],
    [0.0, 2.0, 1.0, 49.0, 51.0, 47.0],
  ]
  let labels = [0, 0, 0, 1, 1, 1]
  let model = @src.scvi_condscvi_init_model(6, 2)
  let trained = @src.scvi_condscvi_train(model, sc, labels)
  assert_eq(trained.trained, true)
  // loss is finite
  assert_true(trained.final_loss == trained.final_loss)
  // latent: n_cells × n_latent
  let lat = @src.scvi_condscvi_latent(trained, sc, labels)
  assert_eq(lat.length(), 6)
  assert_eq(lat[0].length(), 5)
  // denoised expression keeps the library size
  let den = @src.scvi_condscvi_normalized_expression(trained, sc, labels)
  assert_eq(den.length(), 6)
  let lib0 = den[0][0] +
    den[0][1] +
    den[0][2] +
    den[0][3] +
    den[0][4] +
    den[0][5]
  assert_true(
    @moonbit.assert_num_close(
      lib0,
      sc[0][0] + sc[0][1] + sc[0][2] + sc[0][3] + sc[0][4] + sc[0][5],
      1.0e-6,
    ),
  )
  // reference profiles: gene × cell_type, type A dominated by genes 0-2
  assert_eq(trained.profiles.length(), 6)
  assert_eq(trained.profiles[0].length(), 2)
  let a_first = (
      trained.profiles[0][0] + trained.profiles[1][0] + trained.profiles[2][0]
    ) /
    3.0
  let a_last = (
      trained.profiles[3][0] + trained.profiles[4][0] + trained.profiles[5][0]
    ) /
    3.0
  assert_true(a_first > a_last)
  let b_first = (
      trained.profiles[0][1] + trained.profiles[1][1] + trained.profiles[2][1]
    ) /
    3.0
  let b_last = (
      trained.profiles[3][1] + trained.profiles[4][1] + trained.profiles[5][1]
    ) /
    3.0
  assert_true(b_last > b_first)
}

///|
test "scvi_cells2spatial" {
  let sc = [
    [50.0, 50.0, 50.0, 1.0, 1.0, 1.0],
    [55.0, 45.0, 52.0, 0.0, 1.0, 2.0],
    [1.0, 1.0, 2.0, 45.0, 50.0, 48.0],
    [2.0, 0.0, 1.0, 52.0, 48.0, 50.0],
  ]
  let labels = [0, 0, 1, 1]
  let model = @src.scvi_condscvi_init_model(
    6,
    2,
    config=@src.CondScviConfig::default(),
  )
  let trained = @src.scvi_condscvi_train(model, sc, labels)
  // spots: pure A, pure B, 50/50 mix
  let spatial = [
    [100.0, 100.0, 100.0, 2.0, 2.0, 2.0],
    [2.0, 2.0, 4.0, 90.0, 100.0, 96.0],
    [52.0, 50.0, 52.0, 46.0, 48.0, 50.0],
  ]
  let props = @src.scvi_cells2spatial(trained, spatial)
  assert_eq(props.length(), 3)
  assert_eq(props[0].length(), 2)
  // proportions roughly sum to 1
  for s in 0..<3 {
    let total = props[s][0] + props[s][1]
    assert_true(total > 0.5 && total < 1.5)
  }
  // pure A spot is A-dominant, pure B spot is B-dominant
  assert_true(props[0][0] > props[0][1])
  assert_true(props[1][1] > props[1][0])
  // the mixed spot is balanced
  assert_true(@moonbit.assert_num_close(props[2][0], 0.5, 0.35))
}
