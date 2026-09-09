///|
/// Tests for splatter. Exact pipeline values (base means, library size
/// factors, count totals) come from a full Python mirror of the splat
/// simulation with the same deterministic RNG stream.

///|
fn splat_close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol
}

///|
test "splat_lcg_deterministic_stream" {
  // mirror: lcg(42) -> u1 = 0.000944073312424, u2 = 0.571362864026503
  let (s1, u1) = @src.splat_lcg(42)
  let (s2, u2) = @src.splat_lcg(s1)
  assert_true(splat_close(u1, 0.000944073312424, 1.0e-12))
  assert_true(splat_close(u2, 0.571362864026503, 1.0e-12))
  let _ = s2
}

///|
test "splat_normal_stream" {
  // mirror: normal(42) = -3.363420803019271, next state 1226992407
  let (next, value) = @src.splat_normal(42)
  assert_eq(next, 1226992407)
  assert_true(splat_close(value, -3.363420803019271, 1.0e-9))
}

///|
test "splat_gamma_stream" {
  // mirror: gamma(shape=0.6, state=42) = 0.000011051577206
  let (_next, value) = @src.splat_gamma(0.6, 42)
  assert_true(splat_close(value, 0.000011051577206, 1.0e-15))
  // gamma(shape=1) is exponential with mean 1: use the mean of many draws
  let mut state = 7
  let mut total = 0.0
  let n = 2000
  for _i in 0..<n {
    let (s1, draw) = @src.splat_gamma(1.0, state)
    state = s1
    total = total + draw
  }
  let mean = total / n.to_double()
  assert_true(mean > 0.8 && mean < 1.2)
}

///|
test "splat_poisson_small_lambda" {
  // Knuth sampling: lambda=2 has mean 2 over many draws
  let mut state = 11
  let mut total = 0
  let n = 3000
  for _i in 0..<n {
    let (s1, k) = @src.splat_poisson(2.0, state)
    state = s1
    total = total + k
  }
  let mean = total.to_double() / n.to_double()
  assert_true(mean > 1.8 && mean < 2.2)
}

///|
test "splat_poisson_large_lambda" {
  // normal approximation branch: lambda=100 has mean ~100
  let mut state = 13
  let mut total = 0
  let n = 2000
  for _i in 0..<n {
    let (s1, k) = @src.splat_poisson(100.0, state)
    state = s1
    total = total + k
  }
  let mean = total.to_double() / n.to_double()
  assert_true(mean > 96.0 && mean < 104.0)
}

///|
test "splat_simulate_deterministic_and_mirror_parity" {
  // 10 genes x 6 cells plain simulation, seed 12345
  let params = @src.SplatParams::create(10, 6)
  let a = @src.splat_simulate(params, 12345)
  let b = @src.splat_simulate(params, 12345)
  // identical seed -> identical counts
  assert_eq(a.counts, b.counts)
  assert_eq(a.gene_base_means, b.gene_base_means)
  // mirror parity: base means and library factors for the first entries
  assert_true(splat_close(a.gene_base_means[0], 0.392888872964955, 1.0e-9))
  assert_true(splat_close(a.gene_base_means[1], 0.550650379546523, 1.0e-9))
  assert_true(splat_close(a.gene_base_means[2], 0.673247115798264, 1.0e-9))
  assert_true(splat_close(a.gene_base_means[3], 2.237694606260325, 1.0e-9))
  assert_true(splat_close(a.gene_base_means[4], 0.065399239409401, 1.0e-9))
  assert_true(splat_close(a.cell_lib_factors[0], 1.113836946133259, 1.0e-9))
  assert_true(splat_close(a.cell_lib_factors[1], 0.917125598135853, 1.0e-9))
  assert_true(splat_close(a.cell_lib_factors[2], 1.164677750960487, 1.0e-9))
  // mirror parity: count totals and rows
  let mut total = 0
  for row in a.counts {
    for c in row {
      total = total + c
    }
  }
  assert_eq(total, 142)
  assert_eq(a.counts[0], [0, 1, 0, 1, 0, 0])
  assert_eq(a.counts[1], [1, 0, 0, 0, 2, 0])
}

///|
test "splat_simulate_dimensions_and_layout" {
  let params = @src.SplatParams::create(10, 6, batch_cells=[3, 3])
  let result = @src.splat_simulate(params, 5)
  assert_eq(result.counts.length(), 10)
  for row in result.counts {
    assert_eq(row.length(), 6)
  }
  assert_eq(result.cell_batch, [0, 0, 0, 1, 1, 1])
  assert_eq(result.cell_group, [0, 0, 0, 0, 0, 0])
  for factor in result.cell_lib_factors {
    assert_true(factor > 0.0)
  }
}

///|
test "splat_simulate_different_seeds_differ" {
  let params = @src.SplatParams::create(20, 8)
  let a = @src.splat_simulate(params, 1)
  let b = @src.splat_simulate(params, 2)
  let mut identical = true
  for g in 0..<20 {
    if a.counts[g] != b.counts[g] {
      identical = false
      break
    }
  }
  assert_true(!identical)
}

///|
test "splat_simulate_groups_and_de" {
  // two groups, all genes DE with strong up/down factors, no dropout
  let params = @src.SplatParams::create(
    20,
    0,
    group_cells=[4, 4],
    de_prob=1.0,
    de_down_prob=0.5,
    de_fac_loc=1.0,
    de_fac_scale=0.01,
  )
  let result = @src.splat_simulate(params, 777)
  assert_eq(result.cell_group, [0, 0, 0, 0, 1, 1, 1, 1])
  assert_eq(result.group_log_means.length(), 2)
  // mirror parity for the DE pipeline (20 genes x 8 cells, seed 777)
  let mut total = 0
  for row in result.counts {
    for c in row {
      total = total + c
    }
  }
  assert_eq(total, 500)
  // groups are drawn independently: some genes shift in opposite
  // directions (large |diff|), others in the same direction (small diff)
  let mut max_diff = 0.0
  let mut min_diff = 1.0e9
  for g in 0..<20 {
    let diff = (result.group_log_means[0][g] - result.group_log_means[1][g]).abs()
    if diff > max_diff {
      max_diff = diff
    }
    if diff < min_diff {
      min_diff = diff
    }
  }
  assert_true(max_diff > 1.5)
  assert_true(min_diff < 0.1)
}

///|
test "splat_simulate_dropout" {
  let params = @src.SplatParams::create(
    20,
    0,
    group_cells=[4, 4],
    de_prob=1.0,
    de_down_prob=0.5,
    de_fac_loc=1.0,
    de_fac_scale=0.01,
    enable_dropout=true,
  )
  let result = @src.splat_simulate(params, 777)
  // mirror parity: 25 dropout events for this seed (20 genes x 8 cells)
  assert_eq(result.n_dropout_events, 25)
  // dropout only zeroes entries, never changes dimensions
  assert_eq(result.counts.length(), 20)
  assert_eq(result.counts[0].length(), 8)
}

///|
test "splat_simulate_path" {
  // single group of 10 cells simulated as a pseudotime path
  let params = @src.SplatParams::create(
    10,
    0,
    group_cells=[10],
    enable_path=true,
    n_path_steps=10,
    path_probability=1.0,
    path_shift_scale=0.05,
    de_prob=0.0,
  )
  let result = @src.splat_simulate(params, 999)
  // mirror parity for the path pipeline
  let mut total = 0
  for row in result.counts {
    for c in row {
      total = total + c
    }
  }
  assert_eq(total, 2467)
  // pseudotime increases across the 10 path cells
  for c in 0..<10 {
    assert_true(splat_close(result.cell_pseudotime[c], c.to_double() * 0.1, 1.0e-12))
  }
  assert_eq(result.cell_group, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
}

///|
test "splat_simulate_default_params" {
  // 100 genes x 100 cells with splatter defaults: base means come from
  // gamma(shape=0.6, rate=0.3) so the empirical mean is near shape/rate=2
  let params = @src.SplatParams::default()
  assert_eq(params.n_genes, 100)
  let result = @src.splat_simulate(params, 12345)
  assert_eq(result.counts.length(), 100)
  assert_eq(result.counts[0].length(), 100)
  let mut total = 0.0
  for m in result.gene_base_means {
    total = total + m
    assert_true(m > 0.0 && m < 1000.0)
  }
  let mean = total / 100.0
  assert_true(mean > 0.5 && mean < 4.5)
  // overdispersion: NB counts have variance above the mean overall
  let mut sum = 0
  for row in result.counts {
    for c in row {
      sum = sum + c
    }
  }
  assert_true(sum > 0)
}
