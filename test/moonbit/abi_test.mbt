///|
/// ABI trace file format parsing module tests.

///|
test "abi_trace_new" {
  let trace = @src.AbiTrace::new()
  assert_eq(trace.sample_name, "")
  assert_eq(trace.base_order, "GATC")
  assert_eq(trace.bases, "")
  assert_eq(trace.num_bases, 0)
  assert_eq(trace.num_points, 0)
  assert_eq(trace.traces.length(), 0)
}

///|
test "abi_date_creation" {
  let date = @src.AbiDate::new(2024, 1, 15)
  assert_eq(date.year, 2024)
  assert_eq(date.month, 1)
  assert_eq(date.day, 15)
}

///|
test "abi_time_creation" {
  let time = @src.AbiTime::new(10, 30, 45)
  assert_eq(time.hours, 10)
  assert_eq(time.minutes, 30)
  assert_eq(time.seconds, 45)
}

///|
test "abi_date_to_string" {
  let date = @src.AbiDate::new(2024, 3, 5)
  let s = date.to_string()
  assert_eq(s, "2024-03-05")
}

///|
test "abi_time_to_string" {
  let time = @src.AbiTime::new(8, 5, 3)
  let s = time.to_string()
  assert_eq(s, "08:05:03")
}

///|
test "abi_date_to_string_no_padding" {
  let date = @src.AbiDate::new(2024, 12, 25)
  let s = date.to_string()
  assert_eq(s, "2024-12-25")
}

///|
test "abi_time_to_string_no_padding" {
  let time = @src.AbiTime::new(14, 30, 0)
  let s = time.to_string()
  assert_eq(s, "14:30:00")
}

///|
test "abi_dir_entry_new" {
  let entry = @src.AbiDirEntry::new("DATA", 1)
  assert_eq(entry.tag_name, "DATA")
  assert_eq(entry.tag_number, 1)
  assert_eq(entry.element_type, 0)
  assert_eq(entry.element_size, 0)
  assert_eq(entry.num_elements, 0)
  assert_eq(entry.data_size, 0)
}

///|
test "abi_sample_trace" {
  let trace = @src.abi_sample_trace()
  assert_eq(trace.sample_name, "Sample001")
  assert_eq(trace.capillary_id, "CAP-1")
  assert_eq(trace.base_order, "GATC")
  assert_eq(trace.bases, "ACGTACGTACGT")
  assert_eq(trace.num_bases, 12)
  assert_eq(trace.base_positions.length(), 12)
  assert_eq(trace.qualities.length(), 12)
  assert_eq(trace.traces.length(), 4) // 4 channels
  assert_eq(trace.trace_names.length(), 4)
  assert_eq(trace.num_points, 200)
}

///|
test "abi_sample_trace_quality_at" {
  let trace = @src.abi_sample_trace()
  assert_eq(trace.quality_at(0), 20)
  assert_eq(trace.quality_at(2), 40)
  assert_eq(trace.quality_at(-1), -1)
  assert_eq(trace.quality_at(100), -1)
}

///|
test "abi_sample_trace_peak_at" {
  let trace = @src.abi_sample_trace()
  assert_eq(trace.peak_at(0), 10)
  assert_eq(trace.peak_at(1), 25)
  assert_eq(trace.peak_at(-1), -1)
  assert_eq(trace.peak_at(100), -1)
}

///|
test "abi_sample_trace_seq_length" {
  let trace = @src.abi_sample_trace()
  assert_eq(trace.seq_length(), 12)
}

///|
test "abi_sample_trace_trace_for_base" {
  let trace = @src.abi_sample_trace()
  let g_trace = trace.trace_for_base("G")
  assert_eq(g_trace.length(), 200)
  // Peak at position 40 (G peak)
  assert_eq(g_trace[40], 800)
  let a_trace = trace.trace_for_base("A")
  assert_eq(a_trace[10], 800)
}

///|
test "abi_sample_trace_run_date" {
  let trace = @src.abi_sample_trace()
  assert_eq(trace.run_date.year, 2024)
  assert_eq(trace.run_date.month, 1)
  assert_eq(trace.run_date.day, 15)
}

///|
test "abi_sample_trace_run_time" {
  let trace = @src.abi_sample_trace()
  assert_eq(trace.run_time.hours, 10)
  assert_eq(trace.run_time.minutes, 30)
  assert_eq(trace.run_time.seconds, 0)
}

///|
test "abi_sample_bytes_magic" {
  let bytes = @src.abi_sample_bytes()
  // First 4 bytes should be "ABIF" = 0x41, 0x42, 0x49, 0x46
  assert_eq(bytes[0], 0x41)
  assert_eq(bytes[1], 0x42)
  assert_eq(bytes[2], 0x49)
  assert_eq(bytes[3], 0x46)
}

///|
test "abi_parse_sample_bytes" {
  let bytes = @src.abi_sample_bytes()
  let trace = @src.abi_parse(bytes)
  // Should parse the PBAS tag with "ACGT"
  assert_eq(trace.bases, "ACGT")
  assert_eq(trace.num_bases, 4)
  // Should parse FWO_ tag with "GATC"
  assert_eq(trace.base_order, "GATC")
}

///|
test "abi_parse_empty_bytes" {
  let bytes : Array[Int] = Array::new()
  let trace = @src.abi_parse(bytes)
  // Empty input should return an empty trace
  assert_eq(trace.bases, "")
  assert_eq(trace.num_bases, 0)
}

///|
test "abi_parse_bad_magic" {
  let bytes : Array[Int] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  let trace = @src.abi_parse(bytes)
  // Bad magic should return an empty trace
  assert_eq(trace.bases, "")
}

///|
test "abi_encode_basic" {
  let trace = @src.abi_sample_trace()
  let bytes = @src.abi_encode(trace)
  // Should start with ABIF magic
  assert_eq(bytes[0], 0x41)
  assert_eq(bytes[1], 0x42)
  assert_eq(bytes[2], 0x49)
  assert_eq(bytes[3], 0x46)
}
