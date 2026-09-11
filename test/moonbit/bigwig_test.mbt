///|
fn bigwig_target(name : String, length : Int) -> @src.BigBedTarget {
  @src.BigBedTarget::create(name, length) catch {
    _ => abort("test chromosome should be valid")
  }
}

///|
fn bigwig_approx(
  left : Double,
  right : Double,
  tolerance? : Double = 1.0e-9,
) -> Bool {
  (left - right).abs() <= tolerance
}

///|
fn bigwig_lcg_values(count : Int) -> Array[Double] {
  let values : Array[Double] = []
  let mut seed = 1
  for _ in 0..<count {
    seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF
    values.push(seed.to_double() / 1073741824.0 - 1.0)
  }
  values
}

///|
fn bigwig_raises(operation : () -> Unit raise @src.BigWigError) -> Bool {
  try {
    operation()
    false
  } catch {
    _ => true
  }
}

///|
fn bigwig_parsed(bytes : Array[Int]) -> @src.BigWigFile {
  @src.bigwig_parse(bytes) catch {
    @src.BigWigError(msg) => abort("BigWig fixture should parse: " + msg)
  }
}

///|
fn bigwig_query(
  file : @src.BigWigFile,
  chrom : String,
  start? : Int = 0,
  end_? : Int? = None,
) -> Array[@src.BigWigInterval] {
  file.intervals(chrom, start~, end_~) catch {
    _ => abort("BigWig interval query should not fail")
  }
}

///|
fn bigwig_per_base(
  file : @src.BigWigFile,
  chrom : String,
  start : Int,
  end_ : Int,
) -> Array[Double] {
  file.values(chrom, start=start, end_=Some(end_)) catch {
    _ => abort("BigWig values query should not fail")
  }
}

///|
fn bigwig_binned(
  file : @src.BigWigFile,
  chrom : String,
  start? : Int = 0,
  end_? : Int? = None,
  n_bins? : Int = 1,
  stat? : String = "mean",
) -> Array[Double?] {
  file.stats(chrom, start~, end_~, n_bins~, stat~) catch {
    _ => abort("BigWig stats query should not fail")
  }
}

///|
fn bigwig_present(value : Double?) -> Double {
  match value {
    Some(number) => number
    None => abort("expected a covered BigWig bin")
  }
}

///|
fn bigwig_roundtrip_file(compress : Bool) -> @src.BigWigFile {
  let targets = [
    bigwig_target("chr1", 1000),
    bigwig_target("chr2", 500),
    bigwig_target("chr3", 32000),
  ]
  let intervals : Array[@src.BigWigInterval] = []
  intervals.push(@src.BigWigInterval::create("chr1", 0, 10, 1.0))
  intervals.push(@src.BigWigInterval::create("chr1", 10, 20, 2.5))
  intervals.push(@src.BigWigInterval::create("chr1", 20, 30, -3.0))
  intervals.push(@src.BigWigInterval::create("chr2", 0, 100, 7.5))
  let chr3_values = bigwig_lcg_values(3000)
  for i in 0..<3000 {
    intervals.push(
      @src.BigWigInterval::create("chr3", i * 10, i * 10 + 10, chr3_values[i]),
    )
  }
  let config = @src.BigWigWriteConfig::create(compress~) catch {
    _ => abort("BigWig config should be valid")
  }
  let bytes = @src.bigwig_write(targets, intervals, config~) catch {
    _ => abort("BigWig roundtrip should serialize")
  }
  bigwig_parsed(bytes)
}

///|
test "bigwig roundtrip chromosomes and header" {
  let file = bigwig_roundtrip_file(true)
  inspect(file.chroms(), content="[(chr1, 1000), (chr2, 500), (chr3, 32000)]")
  let header = file.header()
  inspect(header.version, content="4")
  inspect(header.zoom_levels, content="0")
  inspect(file.section_count(), content="4")
  let summary = file.total_summary()
  inspect(summary.valid_count, content="30130")
  if !bigwig_approx(summary.min_val, -3.0) {
    abort("summary min mismatch")
  }
  if !bigwig_approx(summary.max_val, 7.5) {
    abort("summary max mismatch")
  }
  let expected_values = bigwig_lcg_values(3000)
  let mut expected_sum = 1.0 * 10.0 + 2.5 * 10.0 + -3.0 * 10.0 + 7.5 * 100.0
  let mut expected_sumsq =
    1.0 * 1.0 * 10.0 +
    2.5 * 2.5 * 10.0 +
    3.0 * 3.0 * 10.0 +
    7.5 * 7.5 * 100.0
  for i in 0..<3000 {
    expected_sum = expected_sum + expected_values[i] * 10.0
    expected_sumsq = expected_sumsq +
      expected_values[i] * expected_values[i] * 10.0
  }
  if !bigwig_approx(summary.sum_data, expected_sum, tolerance=1.0e-3) {
    abort("summary sum mismatch")
  }
  if !bigwig_approx(summary.sum_squares, expected_sumsq, tolerance=1.0e-3) {
    abort("summary sum of squares mismatch")
  }
}

///|
test "bigwig roundtrip intervals across sections and compression" {
  let file = bigwig_roundtrip_file(true)
  let raw_file = bigwig_roundtrip_file(false)
  for candidate in [file, raw_file] {
    let chr1 = bigwig_query(candidate, "chr1")
    inspect(chr1.length(), content="3")
    inspect(chr1[0].start, content="0")
    inspect(chr1[0].end, content="10")
    if !bigwig_approx(chr1[0].value, 1.0) {
      abort("chr1[0] value mismatch")
    }
    inspect(chr1[2].end, content="30")
    let chr2 = bigwig_query(candidate, "chr2")
    inspect(chr2.length(), content="1")
    if !bigwig_approx(chr2[0].value, 7.5) {
      abort("chr2 value mismatch")
    }
    let chr3 = bigwig_query(candidate, "chr3")
    inspect(chr3.length(), content="3000")
    let expected = bigwig_lcg_values(3000)
    inspect(chr3[0].start, content="0")
    inspect(chr3[2727].end, content="27280")
    inspect(chr3[2999].start, content="29990")
    if !bigwig_approx(chr3[2999].value, expected[2999], tolerance=1.0e-6) {
      abort("multi-section value mismatch")
    }
    // Overlapping query returns whole entries unclipped (libBigWig semantics).
    let partial = bigwig_query(candidate, "chr1", start=5, end_=Some(25))
    inspect(partial.length(), content="3")
  }
}

///|
test "bigwig values with NaN gaps" {
  for candidate in [bigwig_roundtrip_file(true), bigwig_roundtrip_file(false)] {
    let values = bigwig_per_base(candidate, "chr1", 0, 35)
    inspect(values.length(), content="35")
    if !bigwig_approx(values[0], 1.0) || !bigwig_approx(values[15], 2.5) ||
      !bigwig_approx(values[25], -3.0) {
      abort("covered values mismatch")
    }
    if !values[30].is_nan() || !values[34].is_nan() {
      abort("uncovered bases must be NaN")
    }
  }
}

///|
test "bigwig raw statistics" {
  let file = bigwig_roundtrip_file(true)
  if !bigwig_approx(bigwig_present(bigwig_binned(file, "chr1")[0]), 5.0 / 30.0) {
    abort("chr1 mean mismatch")
  }
  if !bigwig_approx(bigwig_present(bigwig_binned(file, "chr1", stat="cov")[0]), 0.03) {
    abort("chr1 cov mismatch")
  }
  if !bigwig_approx(bigwig_present(bigwig_binned(file, "chr1", stat="sum")[0]), 5.0) {
    abort("chr1 sum mismatch")
  }
  if !bigwig_approx(bigwig_present(bigwig_binned(file, "chr1", stat="min")[0]), -3.0) {
    abort("chr1 min mismatch")
  }
  if !bigwig_approx(bigwig_present(bigwig_binned(file, "chr1", stat="max")[0]), 2.5) {
    abort("chr1 max mismatch")
  }
  let expected_std = ((162.5 - 25.0 / 30.0) / 29.0).sqrt()
  if !bigwig_approx(
    bigwig_present(bigwig_binned(file, "chr1", stat="std")[0]),
    expected_std,
  ) {
    abort("chr1 std mismatch")
  }
  // Four bins over [0,100): [0.8, -3, gap, gap].
  let four = bigwig_binned(file, "chr1", end_=Some(100), n_bins=4)
  inspect(four.length(), content="4")
  if !bigwig_approx(bigwig_present(four[0]), 0.8) {
    abort("bin 0 mean mismatch")
  }
  if !bigwig_approx(bigwig_present(four[1]), -3.0) {
    abort("bin 1 mean mismatch")
  }
  if !(four[2] is None) || !(four[3] is None) {
    abort("gap bins must be None")
  }
  // chr3 aggregate across both sections: exact full coverage.
  let chr3_sum = bigwig_present(bigwig_binned(file, "chr3", stat="sum")[0])
  let expected_values = bigwig_lcg_values(3000)
  let mut expected = 0.0
  for value in expected_values {
    expected = expected + value * 10.0
  }
  if !bigwig_approx(chr3_sum, expected, tolerance=1.0e-3) {
    abort("chr3 sum mismatch")
  }
  // Cumulative-floor bins can differ in width by one base; coverage is the
  // covered-base fraction of each bin (authoritative pyBigWig value: the last
  // bin is 8667/10667, not an even 1.0).
  let cov = bigwig_binned(file, "chr3", n_bins=3, stat="cov")
  let chr3_length = 32000
  let chr3_covered = 30000
  for i in 0..<3 {
    let bin_start = (i.to_double() * chr3_length.to_double() / 3.0).floor()
    let bin_end = ((i + 1).to_double() * chr3_length.to_double() / 3.0).floor()
    let covered = (chr3_covered.to_double().min(bin_end) - bin_start).max(0.0)
    let expected_cov = covered / (bin_end - bin_start)
    if !bigwig_approx(bigwig_present(cov[i]), expected_cov) {
      abort("chr3 coverage mismatch")
    }
  }
}

///|
test "bigwig parses pyBigWig nozoom fixture" {
  let file = bigwig_parsed(bigwig_nozoom_bytes())
  inspect(file.chroms(), content="[(chr1, 1000), (chr2, 500)]")
  inspect(file.section_count(), content="2")
  let summary = file.total_summary()
  inspect(summary.valid_count, content="130")
  if !bigwig_approx(summary.min_val, -3.0) ||
    !bigwig_approx(summary.max_val, 7.5) ||
    !bigwig_approx(summary.sum_data, 755.0) ||
    !bigwig_approx(summary.sum_squares, 5787.5) {
    abort("nozoom summary mismatch")
  }
  let chr1 = bigwig_query(file, "chr1")
  inspect(chr1.length(), content="3")
  if !bigwig_approx(chr1[2].value, -3.0) {
    abort("nozoom chr1 value mismatch")
  }
  let chr2 = bigwig_query(file, "chr2")
  inspect(chr2.length(), content="1")
  if !bigwig_approx(chr2[0].value, 7.5) {
    abort("chr2 value mismatch")
  }
  let values = bigwig_per_base(file, "chr2", 0, 101)
  inspect(values.length(), content="101")
  if !bigwig_approx(values[99], 7.5) || !values[100].is_nan() {
    abort("chr2 per-base values mismatch")
  }
}

///|
test "bigwig uses zoom data for wide bins" {
  let file = bigwig_parsed(bigwig_small_bytes())
  inspect(file.zoom_levels().length(), content="1")
  inspect(file.zoom_levels()[0].reduction_level, content="512")
  // chr1 bin width 1000 >= reduction 512: zoom aggregation.
  let expectations : Array[(String, Double)] = [
    ("mean", 5.0 / 30.0),
    ("min", -3.0),
    ("max", 2.5),
    ("cov", 0.03),
    ("sum", 5.0),
    ("std", 2.3610829387546217),
  ]
  for i in 0..<expectations.length() {
    let (stat, expected) = expectations[i]
    let value = bigwig_present(bigwig_binned(file, "chr1", stat=stat)[0])
    if !bigwig_approx(value, expected, tolerance=1.0e-5) {
      abort("zoom " + stat + " mismatch: " + value.to_string())
    }
  }
  // chr2 bin width 500 < reduction 512: exact raw fallback.
  if !bigwig_approx(bigwig_present(bigwig_binned(file, "chr2")[0]), 7.5) {
    abort("chr2 raw fallback mean mismatch")
  }
  if !bigwig_approx(bigwig_present(bigwig_binned(file, "chr2", stat="cov")[0]), 0.2) {
    abort("chr2 raw fallback cov mismatch")
  }
  let five = bigwig_binned(file, "chr2", n_bins=5)
  if !bigwig_approx(bigwig_present(five[0]), 7.5) {
    abort("chr2 first bin mismatch")
  }
  if !(five[1] is None) || !(five[4] is None) {
    abort("chr2 trailing bins must be None")
  }
  // Narrow queries stay exact even with a zoom level present.
  let narrow = bigwig_binned(file, "chr1", end_=Some(100), n_bins=4)
  if !bigwig_approx(bigwig_present(narrow[0]), 0.8) {
    abort("narrow bin 0 mismatch")
  }
  if !bigwig_approx(bigwig_present(narrow[1]), -3.0) {
    abort("narrow bin 1 mismatch")
  }
  let three = bigwig_binned(file, "chr1", end_=Some(30), n_bins=3)
  if !bigwig_approx(bigwig_present(three[0]), 1.0) ||
    !bigwig_approx(bigwig_present(three[1]), 2.5) ||
    !bigwig_approx(bigwig_present(three[2]), -3.0) {
    abort("narrow three-bin stats mismatch")
  }
}

///|
test "bigwig decodes variableStep fixedStep and bedGraph sections" {
  let file = bigwig_parsed(bigwig_custom_bytes())
  inspect(file.chroms(), content="[(chr1, 500), (chr2, 400), (chr3, 300)]")
  let chr1 = bigwig_query(file, "chr1")
  inspect(chr1.length(), content="3")
  inspect(chr1[0].start, content="10")
  inspect(chr1[0].end, content="15")
  if !bigwig_approx(chr1[0].value, 1.5) {
    abort("variableStep first value mismatch")
  }
  inspect(chr1[1].start, content="25")
  inspect(chr1[1].end, content="30")
  inspect(chr1[2].start, content="100")
  inspect(chr1[2].end, content="105")
  if !bigwig_approx(chr1[2].value, 3.75) {
    abort("variableStep value mismatch")
  }
  let chr2 = bigwig_query(file, "chr2")
  inspect(chr2.length(), content="4")
  inspect(chr2[0].start, content="20")
  inspect(chr2[0].end, content="24")
  if !bigwig_approx(chr2[0].value, 0.5) {
    abort("fixedStep first value mismatch")
  }
  inspect(chr2[3].start, content="50")
  inspect(chr2[3].end, content="54")
  if !bigwig_approx(chr2[3].value, 4.0) {
    abort("fixedStep value mismatch")
  }
  let chr3 = bigwig_query(file, "chr3")
  inspect(chr3.length(), content="2")
  inspect(chr3[0].start, content="0")
  inspect(chr3[0].end, content="7")
  if !bigwig_approx(chr3[0].value, 9.5) {
    abort("bedGraph first value mismatch")
  }
  inspect(chr3[1].end, content="200")
  let summary = file.total_summary()
  inspect(summary.valid_count, content="138")
  if !bigwig_approx(summary.min_val, -2.25) ||
    !bigwig_approx(summary.max_val, 9.5) ||
    !bigwig_approx(summary.sum_data, 4.5) ||
    !bigwig_approx(summary.sum_squares, 912.875, tolerance=1.0e-6) {
    abort("custom summary mismatch")
  }
  // bin0 [0,100): covered 16 bases, sum 23 -> mean 23/16.
  let bins = bigwig_binned(file, "chr2", end_=Some(400), n_bins=4)
  if !bigwig_approx(bigwig_present(bins[0]), 1.4375) {
    abort("custom chr2 mean mismatch")
  }
  if !(bins[1] is None) || !(bins[3] is None) {
    abort("custom trailing bins must be None")
  }
}

///|
test "bigwig writer validates input" {
  let targets = [
    bigwig_target("chr1", 100),
    bigwig_target("chr2", 100),
  ]
  let good = [
    @src.BigWigInterval::create("chr1", 0, 10, 1.0),
    @src.BigWigInterval::create("chr2", 0, 5, 2.0),
  ]
  let _ = @src.bigwig_write(targets, good) catch {
    _ => abort("valid input should write")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let bad = [@src.BigWigInterval::create("chr3", 0, 10, 1.0)]
    let _ = @src.bigwig_write(targets, bad)
  }) {
    abort("unknown chromosome should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let bad = [@src.BigWigInterval::create("chr1", 90, 101, 1.0)]
    let _ = @src.bigwig_write(targets, bad)
  }) {
    abort("out-of-range interval should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let bad = [@src.BigWigInterval::create("chr1", 20, 20, 1.0)]
    let _ = @src.bigwig_write(targets, bad)
  }) {
    abort("empty interval should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let bad = [
      @src.BigWigInterval::create("chr1", 0, 10, 1.0),
      @src.BigWigInterval::create("chr1", 5, 15, 2.0),
    ]
    let _ = @src.bigwig_write(targets, bad)
  }) {
    abort("overlapping intervals should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let bad = [@src.BigWigInterval::create("chr1", 0, 10, 1.0e39)]
    let _ = @src.bigwig_write(targets, bad)
  }) {
    abort("f32 overflow should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let bad = [@src.BigWigInterval::create("chr1", 0, 10, 0.0 / 0.0)]
    let _ = @src.bigwig_write(targets, bad)
  }) {
    abort("NaN value should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = @src.BigWigWriteConfig::create(block_size=0)
  }) {
    abort("invalid block size should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let dup = [bigwig_target("chr1", 100), bigwig_target("chr1", 100)]
    let _ = @src.bigwig_write(dup, [])
  }) {
    abort("duplicate chromosome names should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = @src.bigwig_write([], [])
  }) {
    abort("missing chromosomes should fail")
  }
}

///|
test "bigwig empty track roundtrips" {
  let targets = [bigwig_target("chr1", 100)]
  let bytes = @src.bigwig_write(targets, []) catch {
    _ => abort("empty track should write")
  }
  let file = bigwig_parsed(bytes)
  inspect(file.section_count(), content="0")
  inspect(file.total_summary().valid_count, content="0")
  inspect(bigwig_query(file, "chr1").length(), content="0")
  let values = bigwig_per_base(file, "chr1", 0, 10)
  inspect(values.length(), content="10")
  if !values[0].is_nan() {
    abort("empty track values must be NaN")
  }
  if !(bigwig_binned(file, "chr1")[0] is None) {
    abort("empty track stats must be None")
  }
  let cov = bigwig_binned(file, "chr1", stat="cov")
  if !bigwig_approx(bigwig_present(cov[0]), 0.0) {
    abort("empty coverage must be zero")
  }
}

///|
test "bigwig reader rejects malformed input" {
  let good = bigwig_nozoom_bytes()
  let bad_magic = good.copy()
  bad_magic[0] = 0
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = @src.bigwig_parse(bad_magic)
  }) {
    abort("bad magic should fail")
  }
  let truncated = good.copy()
  truncated.truncate(good.length() - 4)
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = @src.bigwig_parse(truncated)
  }) {
    abort("missing trailing magic should fail")
  }
  let bad_tail = good.copy()
  bad_tail[bad_tail.length() - 1] = bad_tail[bad_tail.length() - 1] ^ 0xFF
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = @src.bigwig_parse(bad_tail)
  }) {
    abort("corrupt trailing magic should fail")
  }
  let file = bigwig_parsed(good)
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = file.intervals("nope")
  }) {
    abort("unknown chromosome query should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = file.intervals("chr1", start=10, end_=Some(5))
  }) {
    abort("reversed range should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = file.intervals("chr1", start=0, end_=Some(1001))
  }) {
    abort("range past chromosome end should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = file.stats("chr1", n_bins=0)
  }) {
    abort("zero bins should fail")
  }
  if !bigwig_raises(fn() raise @src.BigWigError {
    let _ = file.stats("chr1", stat="median")
  }) {
    abort("unknown statistic should fail")
  }
}
