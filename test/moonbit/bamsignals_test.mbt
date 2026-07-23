///|
/// bamsignals module tests - BAM signal extraction and analysis

///|
test "bamsig_create_example_regions" {
  let regions = @src.bamsig_create_example_regions()
  assert_eq(regions.length(), 4)
}

///|
test "bamsig_create_example_records" {
  let records = @src.bamsig_create_example_records()
  assert_eq(records.length(), 6)
}

///|
test "bamsig_params_new" {
  let params = @src.BamsigParams::new()
  assert_true(params.filter_dup == false)
  assert_eq(params.min_mapq, 10)
  assert_true(params.single_end == false)
  assert_eq(params.extend_len, 150)
}

///|
test "bamsig_region_methods" {
  let region = @src.BamsigRegion::new("chr1", 100, 200, "TestRegion")
  assert_eq(region.bamsig_chrom(), "chr1")
  assert_eq(region.bamsig_start(), 100)
  assert_eq(region.bamsig_end(), 200)
  assert_eq(region.bamsig_width(), 100)
  assert_eq(region.bamsig_id(), "TestRegion")
}

///|
test "bamsig_record_filtering" {
  let records = @src.bamsig_create_example_records()
  let params = @src.BamsigParams::new()

  let mut passing = 0
  for r in records {
    if r.bamsig_passes_filter(params) {
      passing = passing + 1
    }
  }
  assert_true(passing >= 1)

  let strict_params = @src.BamsigParams::new().bamsig_set_min_mapq(30)
  let mut strict_passing = 0
  for r in records {
    if r.bamsig_passes_filter(strict_params) {
      strict_passing = strict_passing + 1
    }
  }
  assert_true(strict_passing <= passing)

  let dup_params = @src.BamsigParams::new().bamsig_set_filter_dup(true)
  let mut dup_passing = 0
  for r in records {
    if r.bamsig_passes_filter(dup_params) {
      dup_passing = dup_passing + 1
    }
  }
  assert_true(dup_passing <= passing)
}

///|
test "bamsig_count_signals" {
  let regions = @src.bamsig_create_example_regions()
  let records = @src.bamsig_create_example_records()
  let params = @src.BamsigParams::new()

  let chr1_records : Array[@src.BamsigRecord] = Array::new()
  let chr2_records : Array[@src.BamsigRecord] = Array::new()
  for r in records {
    if r.chrom == "chr1" {
      chr1_records.push(r)
    } else {
      chr2_records.push(r)
    }
  }

  let records_by_chrom : Map[String, Array[@src.BamsigRecord]] = Map(
    [("chr1", chr1_records), ("chr2", chr2_records)],
    capacity=2,
  )

  let counts = @src.bamsig_count_signals(regions, records_by_chrom, params)
  assert_eq(counts.length(), 4)
  for i = 0; i < counts.length(); i = i + 1 {
    assert_eq(counts[i].length(), 1)
    assert_true(counts[i][0] >= 0.0)
  }
}

///|
test "bamsig_create_example_signal" {
  let signal = @src.bamsig_create_example_signal()
  assert_true(signal.bamsig_n_regions() > 0)
  assert_true(signal.bamsig_n_samples() > 0)

  let n_regions = signal.bamsig_n_regions()
  for i = 0; i < n_regions; i = i + 1 {
    let region_id = signal.bamsig_get_region_id(i)
    assert_true(region_id.length() > 0)
    let count = signal.bamsig_get_count(i, 0)
    assert_true(count >= 0.0)
    let norm_count = signal.bamsig_get_norm_count(i, 0)
    assert_true(norm_count >= 0.0)
  }

  let row_sum = signal.bamsig_row_sum(0)
  assert_true(row_sum >= 0.0)
  let col_sum = signal.bamsig_col_sum(0)
  assert_true(col_sum >= 0.0)
}

///|
test "bamsig_signal_filter" {
  let signal = @src.bamsig_create_example_signal()
  let filtered = signal.bamsig_filter(30.0, 1)
  assert_true(filtered.bamsig_n_regions() <= signal.bamsig_n_regions())
}

///|
test "bamsig_chromatin_state" {
  let signal = @src.bamsig_create_example_signal()
  let promoter_regions = [0, 1]
  let enhancer_regions = [2]
  let repressed_regions = [3]

  let chrom_state = @src.bamsig_analyze_chromatin_states(
    signal, promoter_regions, enhancer_regions, repressed_regions
  )

  assert_true(chrom_state.bamsig_n_states() > 0)
  assert_true(chrom_state.bamsig_active_promoter() >= 0.0)
  assert_true(chrom_state.bamsig_active_enhancer() >= 0.0)
  assert_true(chrom_state.bamsig_repressed() >= 0.0)

  for i = 0; i < chrom_state.bamsig_n_states(); i = i + 1 {
    let label = chrom_state.bamsig_get_label(i)
    assert_true(label.length() > 0)
    let freq = chrom_state.bamsig_get_frequency(i)
    assert_true(!freq.is_nan())
  }
}

///|
test "bamsig_params_setters" {
  let params = @src.BamsigParams::new()
  let params3 = params.bamsig_set_min_mapq(20)
  assert_eq(params3.min_mapq, 20)

  let params4 = params.bamsig_set_extend(200)
  assert_eq(params4.extend_len, 200)
}