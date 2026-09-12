///|
/// Tests for the DNAcopy CBS segmentation module, validated against
/// Bioconductor DNAcopy::segment() behavior.

///|
/// Build a Cna with one probe per 1kb step starting at 1000.
fn dc_make_cna(chrom : Array[String], vals : Array[Double]) -> @src.Cna raise {
  let maploc : Array[Int] = []
  let mut i = 0
  while i < vals.length() {
    maploc.push((i + 1) * 1000)
    i = i + 1
  }
  @src.Cna::new(chrom, maploc, vals)
}

///|
test "dnacopy_segment_four_levels" {
  // 60 probes: normal(0) x15, gain(+1) x15, loss(-1) x15, normal x15
  let vals : Array[Double] = []
  let mut i = 0
  while i < 60 {
    vals.push(
      if i < 15 {
        0.0
      } else if i < 30 {
        1.0
      } else if i < 45 {
        -1.0
      } else {
        0.0
      },
    )
    i = i + 1
  }
  let chrom : Array[String] = Array::make(60, "chr1")
  let cna = dc_make_cna(chrom, vals)
  let res = @src.segment(cna, options=@src.DnaCopyOptions::of(nperm=500))
  assert_eq(res.segments.length(), 4)
  assert_eq(res.segments[0].start, 1000)
  assert_eq(res.segments[0].end, 15000)
  assert_eq(res.segments[0].num_mark, 15)
  assert_eq(res.segments[0].seg_mean, 0.0)
  assert_eq(res.segments[1].start, 16000)
  assert_eq(res.segments[1].end, 30000)
  assert_eq(res.segments[1].seg_mean, 1.0)
  assert_eq(res.segments[2].start, 31000)
  assert_eq(res.segments[2].end, 45000)
  assert_eq(res.segments[2].seg_mean, -1.0)
  assert_eq(res.segments[3].start, 46000)
  assert_eq(res.segments[3].end, 60000)
  assert_eq(res.segments[3].seg_mean, 0.0)
}

///|
test "dnacopy_segment_constant_data" {
  let vals : Array[Double] = Array::make(30, 0.0)
  let chrom : Array[String] = Array::make(30, "chr1")
  let cna = dc_make_cna(chrom, vals)
  let res = @src.segment(cna, options=@src.DnaCopyOptions::of(nperm=200))
  assert_eq(res.segments.length(), 1)
  assert_eq(res.segments[0].num_mark, 30)
  assert_eq(res.segments[0].seg_mean, 0.0)
  assert_eq(res.segments[0].start, 1000)
  assert_eq(res.segments[0].end, 30000)
}

///|
test "dnacopy_segment_single_probe" {
  // fewer probes than 2*min_width: single segment, no changepoints
  let cna = dc_make_cna(["chr1"], [0.5])
  let res = @src.segment(cna, options=@src.DnaCopyOptions::of(nperm=100))
  assert_eq(res.segments.length(), 1)
  assert_eq(res.segments[0].num_mark, 1)
  assert_eq(res.segments[0].seg_mean, 0.5)
}

///|
test "dnacopy_segment_per_chromosome" {
  // chr1: normal(0) x10 + gain(1) x10; chr2: flat(0) x10; chr3: 1 probe
  let chrom : Array[String] = []
  let vals : Array[Double] = []
  let mut i = 0
  while i < 20 {
    chrom.push("chr1")
    vals.push(if i < 10 { 0.0 } else { 1.0 })
    i = i + 1
  }
  i = 0
  while i < 10 {
    chrom.push("chr2")
    vals.push(0.0)
    i = i + 1
  }
  chrom.push("chr3")
  vals.push(-0.5)
  let cna = dc_make_cna(chrom, vals)
  let res = @src.segment(cna, options=@src.DnaCopyOptions::of(nperm=500))
  // chromosomes keep first-appearance order; per-chromosome segmentation
  assert_eq(res.segments.length(), 4)
  assert_eq(res.segments[0].chrom, "chr1")
  assert_eq(res.segments[0].num_mark, 10)
  assert_eq(res.segments[0].seg_mean, 0.0)
  assert_eq(res.segments[1].chrom, "chr1")
  assert_eq(res.segments[1].num_mark, 10)
  assert_eq(res.segments[1].seg_mean, 1.0)
  assert_eq(res.segments[2].chrom, "chr2")
  assert_eq(res.segments[2].num_mark, 10)
  assert_eq(res.segments[3].chrom, "chr3")
  assert_eq(res.segments[3].num_mark, 1)
  assert_eq(res.segments[3].seg_mean, -0.5)
}

///|
test "dnacopy_sdundo_merges_when_sd_large" {
  // four levels with a small alternating noise pattern; a huge undo_sd
  // must merge all adjacent segments into one
  let vals : Array[Double] = []
  let mut i = 0
  while i < 60 {
    let level = if i < 15 {
      0.0
    } else if i < 30 {
      1.0
    } else if i < 45 {
      -1.0
    } else {
      0.0
    }
    vals.push(level + (if i % 2 == 0 { 0.01 } else { -0.01 }))
    i = i + 1
  }
  let chrom : Array[String] = Array::make(60, "chr1")
  let cna = dc_make_cna(chrom, vals)
  let res = @src.segment(
    cna,
    options=@src.DnaCopyOptions::of(
      nperm=500,
      undo_splits="sdundo",
      undo_sd=1000.0,
    ),
  )
  assert_eq(res.segments.length(), 1)
  assert_eq(res.segments[0].num_mark, 60)
}

///|
test "dnacopy_sdundo_keeps_when_sd_tiny" {
  let vals : Array[Double] = []
  let mut i = 0
  while i < 60 {
    let level = if i < 15 {
      0.0
    } else if i < 30 {
      1.0
    } else if i < 45 {
      -1.0
    } else {
      0.0
    }
    vals.push(level + (if i % 2 == 0 { 0.01 } else { -0.01 }))
    i = i + 1
  }
  let chrom : Array[String] = Array::make(60, "chr1")
  let cna = dc_make_cna(chrom, vals)
  let res = @src.segment(
    cna,
    options=@src.DnaCopyOptions::of(
      nperm=500,
      undo_splits="sdundo",
      undo_sd=0.001,
    ),
  )
  assert_eq(res.segments.length(), 4)
}

///|
test "dnacopy_invalid_min_width_raises" {
  let cna = dc_make_cna(Array::make(3, "chr1"), [0.0, 1.0, 0.0])
  let result = try {
    ignore(@src.segment(cna, options=@src.DnaCopyOptions::of(min_width=1)))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
  let result2 = try {
    ignore(@src.segment(cna, options=@src.DnaCopyOptions::of(min_width=6)))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result2)
}

///|
test "dnacopy_invalid_undo_splits_raises" {
  let cna = dc_make_cna(Array::make(4, "chr1"), [0.0, 1.0, 0.0, 0.0])
  let result = try {
    ignore(
      @src.segment(cna, options=@src.DnaCopyOptions::of(undo_splits="prune")),
    )
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
}

///|
test "dnacopy_cna_length_mismatch_raises" {
  let result = try {
    ignore(@src.Cna::new(["chr1", "chr1"], [1, 2], [0.0]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
}

///|
test "dnacopy_segment_mean_rounding" {
  // mean of 10 probes alternating +/- 0.00003 is 0; use 1/3-like mean
  let vals : Array[Double] = Array::make(9, 1.0)
  vals.push(0.0)
  let chrom : Array[String] = Array::make(10, "chr1")
  let cna = dc_make_cna(chrom, vals)
  let res = @src.segment(cna, options=@src.DnaCopyOptions::of(nperm=200))
  // 9/10 = 0.9 exactly representable at 4 decimals
  assert_eq(res.segments.length(), 1)
  assert_eq(res.segments[0].seg_mean, 0.9)
}
