///|
/// Benchmark infrastructure for SFS and coalescent operations.
/// Verifies that SfsStats (one pass) produces the same results as
/// separate function calls (multiple passes), and scales to large inputs.

// =========================================================================
// Benchmark helpers
// =========================================================================

///|
/// Generate n sequences of length L with segregating sites.
/// seq[0] is all-A (ancestral). Each subsequent sequence introduces
/// mutations at distinct positions.
fn bench_make_sequences(n : Int, seq_len : Int) -> Array[String] {
  let seqs : Array[String] = Array(capacity=n)
  // Ancestral: all A
  let ancestral = "A".repeat(seq_len)
  seqs.push(ancestral)
  for i in 1..<n {
    let chars : Array[Int] = Array(capacity=seq_len)
    for j in 0..<seq_len {
      chars.push(65) // 'A'
    }
    // Introduce mutations at positions unique to this sequence
    let mut pos = i
    while pos < seq_len {
      chars[pos] = 71 // 'G'
      pos = pos + n
    }
    let mut s = ""
    for c in chars {
      s = s + Char::to_string(c.unsafe_to_char())
    }
    seqs.push(s)
  }
  seqs
}

///|
/// Generate a random-ish alignment with controlled diversity.
fn bench_make_diverse_sequences(n : Int, seq_len : Int) -> Array[String] {
  let seqs : Array[String] = Array(capacity=n)
  for i in 0..<n {
    let mut s = ""
    for j in 0..<seq_len {
      // Pseudo-random base selection based on position
      let base = (i * 7 + j * 3) % 4
      let c = match base {
        0 => "A"
        1 => "C"
        2 => "G"
        _ => "T"
      }
      s = s + c
    }
    seqs.push(s)
  }
  seqs
}

// =========================================================================
// Scale benchmarks
// =========================================================================

///|
test "bench: SfsStats 10 seqs x 100 sites" {
  let seqs = bench_make_sequences(10, 100)
  let stats = @src.SfsStats::from_sequences(seqs)
  assert_true(stats.segregating_sites() >= 0)
  assert_true(stats.nucleotide_diversity() >= 0.0)
}

///|
test "bench: SfsStats 20 seqs x 200 sites" {
  let seqs = bench_make_sequences(20, 200)
  let stats = @src.SfsStats::from_sequences(seqs)
  // All statistics from one pass
  let _ = stats.tajima_d()
  let _ = stats.fu_li_d()
  let _ = stats.fu_li_f()
  let _ = stats.fay_wu_h()
  let _ = stats.zeng_e()
  assert_true(not(stats.tajima_d().is_nan()))
}

///|
test "bench: SfsStats 50 seqs x 500 sites" {
  let seqs = bench_make_sequences(50, 500)
  let stats = @src.SfsStats::from_sequences(seqs)
  let _ = stats.tajima_d()
  let _ = stats.fu_li_f()
  let _ = stats.fay_wu_h()
  assert_true(not(stats.tajima_d().is_nan()))
}

///|
test "bench: SfsStats diverse 20 seqs x 100 sites" {
  let seqs = bench_make_diverse_sequences(20, 100)
  let stats = @src.SfsStats::from_sequences(seqs)
  assert_true(stats.segregating_sites() >= 0)
  let _ = stats.tajima_d()
  let _ = stats.fu_li_d()
}

// =========================================================================
// Correctness: cached vs uncached must match
// =========================================================================

///|
test "bench: SfsStats matches separate calls" {
  let seqs = bench_make_sequences(10, 50)
  let stats = @src.SfsStats::from_sequences(seqs)
  // Compare each statistic
  let s_sep = @src.phylo_segregating_sites(seqs)
  let pi_sep = @src.phylo_pi(seqs)
  let eta_sep = @src.phylo_count_singletons(seqs)
  let sfs_sep = @src.phylo_observed_sfs(seqs)
  assert_eq(stats.segregating_sites(), s_sep)
  assert_true((stats.nucleotide_diversity() - pi_sep).abs() < 1.0e-10)
  assert_eq(stats.singleton_count(), eta_sep)
  let sfs_c = stats.sfs()
  for i in 0..<sfs_c.length() {
    assert_eq(sfs_c[i], sfs_sep[i])
  }
}

///|
test "bench: all test statistics match separate" {
  let seqs = bench_make_sequences(8, 40)
  let stats = @src.SfsStats::from_sequences(seqs)
  // Tajima D
  let d_c = stats.tajima_d()
  let d_s = @src.phylo_tajima_d(
    stats.segregating_sites(), stats.n_sequences(),
    stats.nucleotide_diversity()
  )
  assert_true((d_c - d_s).abs() < 1.0e-10)
  // Fu Li D
  let fd_c = stats.fu_li_d()
  let fd_s = @src.phylo_fu_li_d(
    stats.segregating_sites(), stats.n_sequences(),
    stats.singleton_count()
  )
  assert_true((fd_c - fd_s).abs() < 1.0e-10)
  // Fay Wu H
  let h_c = stats.fay_wu_h()
  let h_s = @src.phylo_fay_wu_h(
    stats.sfs(), stats.n_sequences(), stats.nucleotide_diversity()
  )
  assert_true((h_c - h_s).abs() < 1.0e-10)
}

// =========================================================================
// Coalescent simulation benchmarks
// =========================================================================

///|
test "bench: coalescent simulation 20 lineages" {
  let rng = @src.Rng::new(42)
  let result = @src.phylo_simulate_coalescent(20, rng)
  assert_eq(@src.phylo_count_leaves(result.tree), 20)
}

///|
test "bench: coalescent simulation 50 lineages" {
  let rng = @src.Rng::new(123)
  let result = @src.phylo_simulate_coalescent(50, rng)
  assert_eq(@src.phylo_count_leaves(result.tree), 50)
}

///|
test "bench: coalescent mean depth 10 reps" {
  let rng = @src.Rng::new(99)
  let depth = @src.phylo_coalescent_mean_depth(10, 1000.0, 10, rng)
  assert_true(depth > 0.0)
}

///|
test "bench: expected SFS computation" {
  let sfs = @src.phylo_expected_sfs(20, 1.0)
  assert_eq(sfs.length(), 19)
  // E[ξ_1] = θ/1 = 1.0
  assert_true((sfs[0] - 1.0).abs() < 1.0e-10)
}

///|
test "bench: watterson theta at scale" {
  // θ_W = S / H_{n-1}; for n=100, S=50: H_99 ≈ 5.177
  // θ_W ≈ 50 / 5.177 ≈ 9.66
  let theta = @src.phylo_watterson_theta(50, 100)
  assert_true(theta > 8.0)
  assert_true(theta < 11.0)
}
