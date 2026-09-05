///|
/// Tests for the regioneR-style region permutation testing module, validated
/// against Bioconductor regioneR::permTest / randomizeRegions behavior.

///|
/// Build a single-chromosome GRanges with star strands.
fn rr_make(seqname : String, ranges : Array[(Int, Int)]) -> @src.GRanges {
  let n = ranges.length()
  @src.granges(Array::make(n, seqname), ranges, Array::make(n, @src.strand_star()))
}

///|
/// Small two-chromosome genome used across the tests: chr1 = 10kb, chr2 = 5kb.
fn rr_genome() -> @src.GenomeLayout raise {
  @src.GenomeLayout::new([("chr1", 10000), ("chr2", 5000)])
}

///|
test "regioner_genome_layout_validation" {
  let result = try {
    ignore(@src.GenomeLayout::new([]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
  let result2 = try {
    ignore(@src.GenomeLayout::new([("chr1", 0)]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result2)
  let result3 = try {
    ignore(@src.GenomeLayout::new([("chr1", 100), ("chr1", 200)]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result3)
  let genome = rr_genome()
  assert_eq(genome.chroms.length(), 2)
}

///|
test "regioner_randomize_preserves_widths_and_chromosomes" {
  let a = rr_make("chr1", [(100, 200), (500, 600), (2000, 2100)])
  let genome = rr_genome()
  let r1 = @src.randomize_regions(a, genome, options=@src.RegionPermOptions::of(seed=42))
  assert_eq(@src.granges_length(r1), 3)
  let mut i = 0
  while i < 3 {
    // per-chromosome default: same chromosome, same width, inside bounds
    assert_eq(r1.seqnames[i], "chr1")
    assert_eq(r1.widths[i], a.widths[i])
    assert_true(r1.starts[i] >= 1)
    assert_true(r1.ends[i] <= 10000)
    assert_eq(r1.ends[i] - r1.starts[i] + 1, a.widths[i])
    i = i + 1
  }
  // reproducible with the same seed
  let r2 = @src.randomize_regions(a, genome, options=@src.RegionPermOptions::of(seed=42))
  assert_eq(r2.starts, r1.starts)
  // a different seed produces a different placement
  let r3 = @src.randomize_regions(a, genome, options=@src.RegionPermOptions::of(seed=43))
  assert_false(r3.starts == r1.starts)
}

///|
test "regioner_randomize_respects_mask" {
  let genome = rr_genome()
  // mask the first 9000 bases of chr1: only 9001..10000 remain
  let mask = rr_make("chr1", [(1, 9000)])
  let a = rr_make("chr1", [(100, 199), (200, 299), (300, 399)])
  let mut trial = 0
  while trial < 10 {
    let r = @src.randomize_regions(
      a,
      genome,
      mask=mask,
      options=@src.RegionPermOptions::of(seed=trial + 1),
    )
    let mut i = 0
    while i < 3 {
      assert_true(r.starts[i] >= 9001)
      assert_true(r.ends[i] <= 10000)
      i = i + 1
    }
    trial = trial + 1
  }
  // a mask leaving no room raises
  let full_mask = rr_make("chr1", [(1, 10000)])
  let result = try {
    ignore(@src.randomize_regions(a, genome, mask=full_mask))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
}

///|
test "regioner_randomize_out_of_bounds_raises" {
  let genome = rr_genome()
  let result = try {
    ignore(@src.randomize_regions(rr_make("chr3", [(1, 10)]), genome))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
  let result2 = try {
    ignore(@src.randomize_regions(rr_make("chr1", [(9995, 10010)]), genome))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result2)
}

///|
test "regioner_randomize_global_pool" {
  // per_chromosome = false lets chr1 regions move to chr2
  let a = rr_make("chr1", [(100, 199), (300, 399), (500, 599), (700, 799)])
  let genome = rr_genome()
  let mut any_on_chr2 = 0
  let mut trial = 0
  while trial < 10 {
    let r = @src.randomize_regions(
      a,
      genome,
      options=@src.RegionPermOptions::of(seed=trial + 1, per_chromosome=false),
    )
    assert_eq(@src.granges_length(r), 4)
    let mut i = 0
    while i < 4 {
      assert_true(r.seqnames[i] == "chr1" || r.seqnames[i] == "chr2")
      assert_eq(r.widths[i], a.widths[i])
      if r.seqnames[i] == "chr2" {
        any_on_chr2 = any_on_chr2 + 1
        assert_true(r.ends[i] <= 5000)
      } else {
        assert_true(r.ends[i] <= 10000)
      }
      i = i + 1
    }
    trial = trial + 1
  }
  // with a 1/3 chance per region, 40 draws must reach chr2 at least once
  assert_true(any_on_chr2 > 0)
}

///|
test "regioner_randomize_multi_chromosome" {
  let genome = rr_genome()
  // two identical-width regions on different chromosomes must stay on their
  // own chromosomes (guards the per-chromosome interval cache)
  let a = @src.granges(
    ["chr1", "chr2"],
    [(100, 199), (200, 299)],
    [@src.strand_star(), @src.strand_star()],
  )
  let r = @src.randomize_regions(a, genome, options=@src.RegionPermOptions::of(seed=9))
  assert_eq(r.seqnames[0], "chr1")
  assert_eq(r.seqnames[1], "chr2")
  assert_eq(r.widths[0], 100)
  assert_eq(r.widths[1], 100)
  assert_true(r.ends[0] <= 10000)
  assert_true(r.ends[1] <= 5000)
}

///|
test "regioner_circular_randomize" {
  let a = rr_make("chr1", [(100, 199), (5000, 5099), (9000, 9099)])
  let genome = rr_genome()
  let r = @src.circular_randomize_regions(a, genome, options=@src.RegionPermOptions::of(seed=5))
  assert_eq(@src.granges_length(r), 3)
  let mut i = 0
  while i < 3 {
    assert_eq(r.seqnames[i], "chr1")
    assert_eq(r.widths[i], a.widths[i])
    assert_true(r.starts[i] >= 1)
    assert_true(r.ends[i] <= 10000)
    i = i + 1
  }
  // deterministic and different from the input placement
  let r2 = @src.circular_randomize_regions(a, genome, options=@src.RegionPermOptions::of(seed=5))
  assert_eq(r2.starts, r.starts)
  assert_false(r.starts == a.starts)
}

///|
test "regioner_stat_num_overlaps" {
  let a = rr_make("chr1", [(100, 200), (300, 400), (500, 600)])
  let b = rr_make("chr1", [(150, 250), (350, 450), (600, 700), (150, 160)])
  // every pairwise overlap counted
  assert_eq(@src.stat_num_overlaps(a, b), 4)
  // each region of a counted once
  assert_eq(@src.stat_num_overlaps(a, b, count_once=true), 3)
  // cross-chromosome overlaps do not count
  let c = rr_make("chr2", [(150, 250)])
  assert_eq(@src.stat_num_overlaps(a, c), 0)
}

///|
test "regioner_stat_mean_distance" {
  let b = rr_make("chr1", [(300, 400)])
  // gap between 200 and 300 is 99 bases
  assert_eq(@src.stat_mean_distance(rr_make("chr1", [(100, 200)]), b), 99.0)
  // overlapping regions are at distance 0
  assert_eq(@src.stat_mean_distance(rr_make("chr1", [(350, 380)]), b), 0.0)
  // regions on chromosomes absent from b are skipped
  let two = @src.granges(
    ["chr1", "chr2"],
    [(100, 200), (100, 200)],
    [@src.strand_star(), @src.strand_star()],
  )
  assert_eq(@src.stat_mean_distance(two, b), 99.0)
  // no distances at all -> NaN
  let empty = rr_make("chr2", [])
  assert_true(@src.stat_mean_distance(rr_make("chr1", [(100, 200)]), empty).is_nan())
}

///|
test "regioner_overlap_perm_test_associated" {
  let genome = rr_genome()
  // b: 20 probes of width 51 spaced 450 bases apart on chr1
  let ranges : Array[(Int, Int)] = []
  let mut k = 0
  while k < 20 {
    ranges.push((250 + k * 450, 300 + k * 450))
    k = k + 1
  }
  let b = rr_make("chr1", ranges)
  // a: the first 5 probes of b, so a is strongly associated with b
  let a_ranges : Array[(Int, Int)] = []
  let mut j = 0
  while j < 5 {
    a_ranges.push(ranges[j])
    j = j + 1
  }
  let a = rr_make("chr1", a_ranges)
  let res = @src.overlap_perm_test(
    a,
    b,
    genome,
    options=@src.RegionPermOptions::of(ntimes=200, seed=1),
  )
  // observed: each of the 5 regions overlaps exactly one probe
  assert_eq(res.statistic, 5.0)
  assert_eq(res.permuted.length(), 200)
  assert_true(res.p_value <= 0.05)
  assert_true(res.z_score > 1.5)
  // permutation distribution is mostly below the observed value
  let mut ge = 0
  for p in res.permuted {
    if p >= res.statistic {
      ge = ge + 1
    }
  }
  assert_true(ge <= 10)
}

///|
test "regioner_perm_test_null_and_alternatives" {
  let genome = rr_genome()
  // a: 5 scattered probes, b: a single probe in between
  let a = rr_make("chr1", [(100, 149), (2000, 2049), (4000, 4049), (6000, 6049), (8000, 8049)])
  let b = rr_make("chr1", [(4900, 4949)])
  // with alternative = "less", a low observed overlap is significant
  let res_less = @src.perm_test(
    a,
    b,
    genome,
    fn(x, y) { @src.stat_num_overlaps(x, y, count_once=true).to_double() },
    options=@src.RegionPermOptions::of(ntimes=100, seed=3, alternative="less"),
  )
  assert_eq(res_less.statistic, 0.0)
  // the observed overlap (0, the minimum) is never significant under "less"
  assert_true(res_less.p_value >= 0.5)
  // the observed overlap lies at or below the permutation mean
  assert_true(res_less.z_score <= 0.0)
  // two.sided p-value stays within [0, 1]
  let res_two = @src.perm_test(
    a,
    b,
    genome,
    fn(x, y) { @src.stat_num_overlaps(x, y, count_once=true).to_double() },
    options=@src.RegionPermOptions::of(ntimes=50, seed=3, alternative="two.sided"),
  )
  assert_true(res_two.p_value >= 0.0 && res_two.p_value <= 1.0)
  // invalid alternative and ntimes raise
  let result = try {
    ignore(
      @src.perm_test(
        a,
        b,
        genome,
        fn(x, y) { @src.stat_num_overlaps(x, y).to_double() },
        options=@src.RegionPermOptions::of(alternative="bogus"),
      ),
    )
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
  let result2 = try {
    ignore(
      @src.perm_test(
        a,
        b,
        genome,
        fn(x, y) { @src.stat_num_overlaps(x, y).to_double() },
        options=@src.RegionPermOptions::of(ntimes=0),
      ),
    )
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result2)
}

///|
test "regioner_perm_test_reproducible" {
  let genome = rr_genome()
  let a = rr_make("chr1", [(100, 199), (300, 399)])
  let b = rr_make("chr1", [(2000, 2049)])
  let run = fn() -> @src.PermTestResult raise {
    @src.overlap_perm_test(a, b, genome, options=@src.RegionPermOptions::of(ntimes=50, seed=11))
  }
  let r1 = run()
  let r2 = run()
  assert_eq(r1.permuted, r2.permuted)
  assert_eq(r1.p_value, r2.p_value)
  assert_eq(r1.z_score, r2.z_score)
  assert_eq(r1.alternative, "greater")
  assert_eq(r1.ntimes, 50)
}

///|
test "regioner_perm_test_with_mask" {
  let genome = rr_genome()
  // mask forces randomized regions into 9001..10000, far away from b
  let mask = rr_make("chr1", [(1, 9000)])
  let a = rr_make("chr1", [(2000, 2049), (2100, 2149)])
  let b = rr_make("chr1", [(9500, 9549)])
  let res = @src.overlap_perm_test(
    a,
    b,
    genome,
    mask=mask,
    options=@src.RegionPermOptions::of(ntimes=100, seed=1),
  )
  // randomized regions can never reach b, so the observed overlap of 0 is
  // not surprising
  assert_eq(res.statistic, 0.0)
  assert_eq(res.p_value, 1.0)
}
