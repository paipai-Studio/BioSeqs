///|
/// Tests for bioframe_stats.mbt (bf_coverage_per_interval) +
/// bioframe_symdiff.mbt (bf_setdiff / bf_symmetric_diff / bf_jaccard).

// ================================================================
// bf_coverage_per_interval
// ================================================================

test "bf_coverage_per_interval: simple partial overlap" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 50, 150)]
  let cov = @src.bf_coverage_per_interval(a, b)
  assert_eq(cov.length(), 1)
  assert_eq(cov[0], 50)
}

test "bf_coverage_per_interval: merged B covers two A" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 200, 300),
  ]
  let b = [@src.BedInterval::new("chr1", 50, 250)]
  let cov = @src.bf_coverage_per_interval(a, b)
  assert_eq(cov.length(), 2)
  assert_eq(cov[0], 50)
  assert_eq(cov[1], 50)
}

test "bf_coverage_per_interval: empty A → empty output" {
  let cov = @src.bf_coverage_per_interval([], [@src.BedInterval::new("chr1", 0, 100)])
  assert_eq(cov.length(), 0)
}

test "bf_coverage_per_interval: empty B → all zeros" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 200, 300),
  ]
  let cov = @src.bf_coverage_per_interval(a, [])
  assert_eq(cov.length(), 2)
  assert_eq(cov[0], 0)
  assert_eq(cov[1], 0)
}

test "bf_coverage_per_interval: full coverage" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 0, 200)]
  let cov = @src.bf_coverage_per_interval(a, b)
  assert_eq(cov[0], 100)
}

test "bf_coverage_per_interval: B on different chrom → 0" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr2", 0, 200)]
  let cov = @src.bf_coverage_per_interval(a, b)
  assert_eq(cov[0], 0)
}

// ---- bf_coverage_per_interval_fraction ------------------------------------

test "bf_coverage_per_interval_fraction: basic" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 50, 100)]
  let f = @src.bf_coverage_per_interval_fraction(a, b)
  assert_eq(f.length(), 1)
  assert_true((f[0] - 0.5).abs() < 0.0001)
}

test "bf_coverage_per_interval_fraction: zero-length A → 0" {
  let a = [@src.BedInterval::new("chr1", 0, 0)]
  let b = [@src.BedInterval::new("chr1", 0, 100)]
  let f = @src.bf_coverage_per_interval_fraction(a, b)
  assert_eq(f[0], 0.0)
}

// ================================================================
// bf_setdiff  (A - B)
// ================================================================

test "bf_setdiff: splits one interval in two" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 25, 50)]
  let out = @src.bf_setdiff(a, b)
  assert_eq(out.length(), 2)
  assert_eq(out[0].start, 0)
  assert_eq(out[0].end, 25)
  assert_eq(out[1].start, 50)
  assert_eq(out[1].end, 100)
}

test "bf_setdiff: no overlap → A unchanged" {
  let a = [@src.BedInterval::new("chr1", 0, 50)]
  let b = [@src.BedInterval::new("chr1", 100, 200)]
  let out = @src.bf_setdiff(a, b)
  assert_eq(out.length(), 1)
  assert_eq(out[0].start, 0)
  assert_eq(out[0].end, 50)
}

test "bf_setdiff: full coverage → empty" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 0, 200)]
  let out = @src.bf_setdiff(a, b)
  assert_eq(out.length(), 0)
}

test "bf_setdiff: B in middle of A → two fragments" {
  let a = [@src.BedInterval::new("chr1", 0, 200)]
  let b = [@src.BedInterval::new("chr1", 50, 150)]
  let out = @src.bf_setdiff(a, b)
  assert_eq(out.length(), 2)
  assert_eq(out[0].start, 0)
  assert_eq(out[0].end, 50)
  assert_eq(out[1].start, 150)
  assert_eq(out[1].end, 200)
}

test "bf_setdiff: empty B → A unchanged" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 200, 300),
  ]
  let out = @src.bf_setdiff(a, [])
  assert_eq(out.length(), 2)
}

test "bf_setdiff: B on different chrom → no effect" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr2", 50, 150)]
  let out = @src.bf_setdiff(a, b)
  assert_eq(out.length(), 1)
  assert_eq(out[0].start, 0)
  assert_eq(out[0].end, 100)
}

// ================================================================
// bf_symmetric_diff
// ================================================================

test "bf_symmetric_diff: A minus B plus B minus A" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 50, 150)]
  let out = @src.bf_symmetric_diff(a, b)
  assert_eq(out.length(), 2)
  assert_eq(out[0].start, 0)
  assert_eq(out[0].end, 50)
  assert_eq(out[1].start, 100)
  assert_eq(out[1].end, 150)
}

// ================================================================
// bf_jaccard
// ================================================================

test "bf_jaccard: identical sets → 1.0" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 0, 100)]
  let j = @src.bf_jaccard(a, b)
  assert_true((j - 1.0).abs() < 0.0001)
}

test "bf_jaccard: no overlap → 0.0" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 200, 300)]
  let j = @src.bf_jaccard(a, b)
  assert_true((j - 0.0).abs() < 0.0001)
}

// ================================================================
// bf_sets_equal
// ================================================================

test "bf_sets_equal: identical merged sets → true" {
  let a = [
    @src.BedInterval::new("chr1", 0, 50),
    @src.BedInterval::new("chr1", 40, 100),
  ]
  let b = [@src.BedInterval::new("chr1", 0, 100)]
  assert_true(@src.bf_sets_equal(a, b))
}

test "bf_sets_equal: different sets → false" {
  let a = [@src.BedInterval::new("chr1", 0, 100)]
  let b = [@src.BedInterval::new("chr1", 0, 50)]
  assert_true(!@src.bf_sets_equal(a, b))
}
