///|
/// Tests for Site Frequency Spectrum (SFS).

test "sfs: expected unfolded basic" {
  // E[ξ_i] = θ / i
  let sfs = @src.phylo_expected_sfs(5, 1.0)
  assert_eq(sfs.length(), 4)
  assert_true((sfs[0] - 1.0).abs() < 0.001)  // θ/1 = 1
  assert_true((sfs[1] - 0.5).abs() < 0.001)  // θ/2 = 0.5
  assert_true((sfs[2] - (1.0 / 3.0)).abs() < 0.001)  // θ/3
  assert_true((sfs[3] - 0.25).abs() < 0.001)  // θ/4
}

///|
test "sfs: expected unfolded with theta" {
  let sfs = @src.phylo_expected_sfs(4, 10.0)
  assert_eq(sfs.length(), 3)
  assert_true((sfs[0] - 10.0).abs() < 0.001)  // 10/1
  assert_true((sfs[1] - 5.0).abs() < 0.001)   // 10/2
  assert_true((sfs[2] - (10.0 / 3.0)).abs() < 0.001)  // 10/3
}

///|
test "sfs: expected folded" {
  // n=5: folded classes i=1,2
  // η_1 = θ/1 + θ/4 = 1 + 0.25 = 1.25
  // η_2 = θ/2 + θ/3 = 0.5 + 0.333 = 0.833
  let sfs = @src.phylo_expected_sfs_folded(5, 1.0)
  assert_eq(sfs.length(), 2)
  assert_true((sfs[0] - 1.25).abs() < 0.001)
  assert_true((sfs[1] - (5.0 / 6.0)).abs() < 0.001)
}

///|
test "sfs: expected folded even n" {
  // n=4: folded classes i=1,2
  // η_1 = θ/1 + θ/3
  // η_2 = θ/2 (i = n/2, no folding)
  let sfs = @src.phylo_expected_sfs_folded(4, 1.0)
  assert_eq(sfs.length(), 2)
  assert_true((sfs[0] - (1.0 + 1.0 / 3.0)).abs() < 0.001)
  assert_true((sfs[1] - 0.5).abs() < 0.001)
}

///|
test "sfs: expected sfs n=1 empty" {
  let sfs = @src.phylo_expected_sfs(1, 1.0)
  assert_eq(sfs.length(), 0)
}

///|
test "sfs: observed sfs basic" {
  // 4 sequences, 5 sites
  // Seq0: AAAAA (ancestral)
  // Seq1: AAATA (1 derived at pos 3)
  // Seq2: ATAAA (1 derived at pos 2) → singleton
  // Seq3: AATTA (1 derived at pos 2, 1 at pos 3) → doubleton
  let seqs = [
    "AAAAA",
    "AAATA",
    "ATAAA",
    "AATTA",
  ]
  let sfs = @src.phylo_observed_sfs(seqs)
  // n=4 → SFS has 3 entries (i=1,2,3)
  assert_eq(sfs.length(), 3)
  // pos 0: all A → not segregating
  // pos 1: A,A,T,T → derived count=2 → doubleton
  // pos 2: A,T,A,T → derived count=2 → doubleton
  // pos 3: A,A,A,A → not segregating
  // Wait, let me recalculate:
  // pos 0: A,A,A,A → not seg
  // pos 1: A,A,T,T → derived=2 → SFS[1]++
  // pos 2: A,T,A,T → derived=2 → SFS[1]++
  // pos 3: A,A,A,A → not seg
  // pos 4: A,A,A,A → not seg
  // pos 1: A,A,T,A → singleton
  // pos 2: A,A,A,T → singleton
  // pos 3: A,T,A,T → doubleton
  // SFS = [2, 1, 0]
  assert_eq(sfs[0], 2)  // singletons
  assert_eq(sfs[1], 1)  // doubletons
  assert_eq(sfs[2], 0)  // tripletons
}

///|
test "sfs: observed sfs with singletons" {
  let seqs = [
    "AAAA",
    "ATAA",  // singleton at pos 1
    "AATA",  // singleton at pos 2
    "AAAT",  // singleton at pos 3
  ]
  let sfs = @src.phylo_observed_sfs(seqs)
  assert_eq(sfs[0], 3)  // 3 singletons
  assert_eq(sfs[1], 0)
  assert_eq(sfs[2], 0)
}

///|
test "sfs: observed folded sfs" {
  let seqs = [
    "AAAA",
    "ATAA",
    "AATA",
    "AAAT",
  ]
  let sfs = @src.phylo_observed_sfs_folded(seqs)
  // n=4, half=2 → 2 classes
  assert_eq(sfs.length(), 2)
  // 3 singletons (minor=1)
  assert_eq(sfs[0], 3)
  assert_eq(sfs[1], 0)
}

///|
test "sfs: segregating sites" {
  let seqs = [
    "AAAA",
    "ATAA",
    "AATA",
    "AAAT",
  ]
  let s = @src.phylo_segregating_sites(seqs)
  assert_eq(s, 3)
}

///|
test "sfs: count singletons" {
  let seqs = [
    "AAAA",
    "ATAA",
    "AATA",
    "AAAT",
  ]
  let s = @src.phylo_count_singletons(seqs)
  assert_eq(s, 3)
}

///|
test "sfs: watterson theta" {
  // θ_W = S / H_{n-1}
  // n=4, S=3: H_3 = 1 + 1/2 + 1/3 = 1.833
  // θ_W = 3 / 1.833 = 1.636
  let theta = @src.phylo_watterson_theta(3, 4)
  assert_true((theta - 3.0 / (1.0 + 0.5 + 1.0 / 3.0)).abs() < 0.001)
}

///|
test "sfs: watterson theta zero" {
  let theta = @src.phylo_watterson_theta(0, 4)
  assert_eq(theta, 0.0)
}

///|
test "sfs: nucleotide diversity" {
  // 3 sequences, 3 sites
  // AAA, ATA, AAT
  // pairs: (0,1): 1 diff, (0,2): 1 diff, (1,2): 2 diffs
  // π = (1+1+2) / 3 = 4/3
  let seqs = ["AAA", "ATA", "AAT"]
  let pi = @src.phylo_pi(seqs)
  assert_true((pi - 4.0 / 3.0).abs() < 0.001)
}

///|
test "sfs: pi identical sequences" {
  let seqs = ["AAA", "AAA", "AAA"]
  let pi = @src.phylo_pi(seqs)
  assert_eq(pi, 0.0)
}

///|
test "sfs: tajima d basic" {
  // Simple test: n=10, S=5, π=0.5
  let d = @src.phylo_tajima_d(5, 10, 0.5)
  // Just check it produces a finite value
  assert_true(d == d)  // not NaN
}

///|
test "sfs: tajima d zero segregating" {
  let d = @src.phylo_tajima_d(0, 10, 0.0)
  assert_eq(d, 0.0)
}

///|
test "sfs: fu li d basic" {
  let d = @src.phylo_fu_li_d(5, 10, 2)
  assert_true(d == d)  // not NaN
}

///|
test "sfs: fu li f basic" {
  let f = @src.phylo_fu_li_f(5, 10, 0.5, 2)
  assert_true(f == f)  // not NaN
}

///|
test "sfs: fay wu h basic" {
  // SFS = [3, 1, 0] (3 singletons, 1 doubleton, 0 tripleton)
  // θ_H = 2/(n(n-1)) * (1²*3 + 2²*1 + 3²*0) = 2/12 * 7 = 7/6
  // H = π - θ_H
  let sfs = [3, 1, 0]
  let h = @src.phylo_fay_wu_h(sfs, 4, 2.0)
  let theta_h = 2.0 / (4.0 * 3.0) * (1.0 * 3.0 + 4.0 * 1.0 + 9.0 * 0.0)
  assert_true((h - (2.0 - theta_h)).abs() < 0.001)
}

///|
test "sfs: zeng e basic" {
  // θ_L = 1/(n-1) * Σ i*ξ_i
  // SFS = [3, 1, 0], n=4
  // θ_L = 1/3 * (1*3 + 2*1 + 3*0) = 5/3
  // E = π - θ_L
  let sfs = [3, 1, 0]
  let e = @src.phylo_zeng_e(sfs, 4, 2.0)
  let theta_l = (1.0 * 3.0 + 2.0 * 1.0 + 3.0 * 0.0) / 3.0
  assert_true((e - (2.0 - theta_l)).abs() < 0.001)
}

///|
test "sfs: num bins" {
  // Sturge's rule: 1 + log2(n)
  // n=10: 1 + log2(10) ≈ 1 + 3.32 = 4.32 → 4
  let bins = @src.phylo_sfs_num_bins(10)
  assert_true(bins >= 3 && bins <= 5)
}

///|
test "sfs: sfs to string" {
  let sfs = @src.phylo_expected_sfs(4, 1.0)
  let s = @src.phylo_sfs_to_string(sfs)
  assert_true(s.contains("Site Frequency Spectrum"))
  assert_true(s.contains("E[ξ_i]"))
}

///|
test "sfs: expected sfs sums to watterson" {
  // Σ E[ξ_i] = θ * H_{n-1} = E[S]
  let n = 10
  let theta = 2.0
  let sfs = @src.phylo_expected_sfs(n, theta)
  let mut total = 0.0
  for v in sfs {
    total = total + v
  }
  let expected_s = @src.phylo_watterson_expected(n, theta)
  assert_true((total - expected_s).abs() < 0.001)
}

///|
test "sfs: observed sfs single sequence" {
  let seqs = ["AAAA"]
  let sfs = @src.phylo_observed_sfs(seqs)
  assert_eq(sfs.length(), 0)
}

///|
test "sfs: tajima d negative with excess singletons" {
  // Excess rare variants → negative Tajima D
  // 10 sequences with many singletons
  let seqs = [
    "AAAAAAAAAA",
    "TAAAAAAAAA",
    "ATAAAAAAAA",
    "AATAAAAAAA",
    "AAATAAAAAA",
    "AAAATAAAAA",
    "AAAAATAAAA",
    "AAAAAATAAA",
    "AAAAAAATAA",
    "AAAAAAAATA",
  ]
  let s = @src.phylo_segregating_sites(seqs)
  let pi = @src.phylo_pi(seqs)
  let d = @src.phylo_tajima_d(s, 10, pi)
  // With excess singletons, D should be negative
  assert_true(d < 0.0)
}

///|
test "sfs: all sequences identical" {
  let seqs = ["AAAA", "AAAA", "AAAA", "AAAA"]
  let s = @src.phylo_segregating_sites(seqs)
  assert_eq(s, 0)
  let pi = @src.phylo_pi(seqs)
  assert_eq(pi, 0.0)
  let singletons = @src.phylo_count_singletons(seqs)
  assert_eq(singletons, 0)
}

///|
test "sfs: fu li d known value" {
  // Reference: n=10, S=15, eta=3
  // a_n = H_9 = 2.828968, b_n = sum 1/i^2 = 1.539768
  // v_D = 1 + a_n^2/(b_n*(n-2)) - 2*a_n/((n-1)*b_n) = 1.241415
  // u_D = a_n - 1 - v_D = 0.587553
  // D = (15 - 2.828968*3) / sqrt(0.587553*15 + 1.241415*225) = 0.383700
  let d = @src.phylo_fu_li_d(15, 10, 3)
  assert_true((d - 0.383700).abs() < 0.001)
}

///|
test "sfs: fu li f known value" {
  // Reference: n=10, S=15, eta=3, pi=3.5
  // F = (pi - eta) / sqrt(u_F*S + v_F*S^2) = 0.037150
  let f = @src.phylo_fu_li_f(15, 10, 3.5, 3)
  assert_true((f - 0.037150).abs() < 0.001)
}

///|
test "sfs: fay wu h known value" {
  // SFS = [3, 1, 0], n=4, pi=2.0
  // theta_H = 2/(4*3) * (1*3 + 4*1 + 0) = 2/12 * 7 = 7/6 = 1.166667
  // H = 2.0 - 1.166667 = 0.833333
  let sfs = [3, 1, 0]
  let h = @src.phylo_fay_wu_h(sfs, 4, 2.0)
  assert_true((h - 0.833333).abs() < 0.001)
}

///|
test "sfs: fu li d n too small" {
  // n=2: Fu & Li D requires n >= 3 (n-2 in denominator)
  // Should return 0.0 for n <= 2
  let d = @src.phylo_fu_li_d(5, 2, 1)
  assert_eq(d, 0.0)
  // n=3 should work
  let d3 = @src.phylo_fu_li_d(5, 3, 1)
  assert_true(not(d3.is_nan()))
}

///|
// ===== SFS Stats Context Tests =====

///|
test "sfs stats: basic from sequences" {
  // seq[0] = ancestral
  // pos 0: A/A/A/A → not segregating
  // pos 1: A/A/G/G → ξ_2 = 1
  // pos 2: A/G/G/G → ξ_3 = 1
  // pos 3: A/G/G/G → ξ_3 = 1
  let seqs = ["AAAA", "AAGG", "AGGG", "AGGG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  assert_eq(stats.n_sequences(), 4)
  assert_eq(stats.segregating_sites(), 3)
  // singletons: no position has derived_count == 1
  assert_eq(stats.singleton_count(), 0)
}

///|
test "sfs stats: singletons detected" {
  // pos 0: A/G/A/A → ξ_1 = 1 (singleton)
  // pos 1: A/A/A/G → ξ_1 = 1 (singleton)
  let seqs = ["AA", "GA", "AA", "AG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  assert_eq(stats.segregating_sites(), 2)
  assert_eq(stats.singleton_count(), 2)
}

///|
test "sfs stats: pi matches separate call" {
  let seqs = ["AAAA", "AAGG", "AGGG", "AGGG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let pi_sep = @src.phylo_pi(seqs)
  assert_true((stats.nucleotide_diversity() - pi_sep).abs() < 1.0e-10)
}

///|
test "sfs stats: s matches separate call" {
  let seqs = ["AAAA", "AAGG", "AGGG", "AGGG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let s_sep = @src.phylo_segregating_sites(seqs)
  assert_eq(stats.segregating_sites(), s_sep)
}

///|
test "sfs stats: sfs matches separate call" {
  let seqs = ["AAAA", "AAGG", "AGGG", "AGGG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let sfs_sep = @src.phylo_observed_sfs(seqs)
  let sfs_ctx = stats.sfs()
  assert_eq(sfs_ctx.length(), sfs_sep.length())
  for i in 0..<sfs_ctx.length() {
    assert_eq(sfs_ctx[i], sfs_sep[i])
  }
}

///|
test "sfs stats: tajima d matches" {
  let seqs = ["AAAA", "AAGG", "AGGG", "AGGG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let d_ctx = stats.tajima_d()
  let d_sep = @src.phylo_tajima_d(
    stats.segregating_sites(), stats.n_sequences(),
    stats.nucleotide_diversity()
  )
  assert_true((d_ctx - d_sep).abs() < 1.0e-10)
}

///|
test "sfs stats: fu li d matches" {
  let seqs = ["AA", "GA", "AA", "AG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let d_ctx = stats.fu_li_d()
  let d_sep = @src.phylo_fu_li_d(
    stats.segregating_sites(), stats.n_sequences(),
    stats.singleton_count()
  )
  assert_true((d_ctx - d_sep).abs() < 1.0e-10)
}

///|
test "sfs stats: fu li f matches" {
  let seqs = ["AA", "GA", "AA", "AG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let f_ctx = stats.fu_li_f()
  let f_sep = @src.phylo_fu_li_f(
    stats.segregating_sites(), stats.n_sequences(),
    stats.nucleotide_diversity(), stats.singleton_count()
  )
  assert_true((f_ctx - f_sep).abs() < 1.0e-10)
}

///|
test "sfs stats: fay wu h matches" {
  let seqs = ["AAAA", "AAGG", "AGGG", "AGGG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let h_ctx = stats.fay_wu_h()
  let h_sep = @src.phylo_fay_wu_h(
    stats.sfs(), stats.n_sequences(), stats.nucleotide_diversity()
  )
  assert_true((h_ctx - h_sep).abs() < 1.0e-10)
}

///|
test "sfs stats: zeng e matches" {
  let seqs = ["AAAA", "AAGG", "AGGG", "AGGG"]
  let stats = @src.SfsStats::from_sequences(seqs)
  let e_ctx = stats.zeng_e()
  let e_sep = @src.phylo_zeng_e(
    stats.sfs(), stats.n_sequences(), stats.nucleotide_diversity()
  )
  assert_true((e_ctx - e_sep).abs() < 1.0e-10)
}

///|
test "sfs stats: all statistics from one pass" {
  // The key benefit: all stats computed from one SfsStats
  let seqs = ["ACGTACGT", "ACGTACGA", "ACGTACAA", "ACGTAAAA"]
  let stats = @src.SfsStats::from_sequences(seqs)
  // All methods reuse cached S, π, η, SFS — no rescanning
  let _ = stats.tajima_d()
  let _ = stats.fu_li_d()
  let _ = stats.fu_li_f()
  let _ = stats.fay_wu_h()
  let _ = stats.zeng_e()
  let _ = stats.watterson_theta()
  // Just verify no crash and values are finite
  assert_true(not(stats.tajima_d().is_nan()))
  assert_true(not(stats.fu_li_d().is_nan()))
  assert_true(not(stats.fu_li_f().is_nan()))
  assert_true(not(stats.fay_wu_h().is_nan()))
  assert_true(not(stats.zeng_e().is_nan()))
}

///|
test "sfs stats: identical sequences" {
  let seqs = ["AAAA", "AAAA", "AAAA"]
  let stats = @src.SfsStats::from_sequences(seqs)
  assert_eq(stats.segregating_sites(), 0)
  assert_eq(stats.singleton_count(), 0)
  assert_true(stats.nucleotide_diversity() < 1.0e-10)
}

///|
test "sfs stats: single sequence" {
  let stats = @src.SfsStats::from_sequences(["AAAA"])
  assert_eq(stats.n_sequences(), 1)
  assert_eq(stats.segregating_sites(), 0)
}
