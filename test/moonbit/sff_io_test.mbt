///|
/// Bio.SeqIO.SffIO module tests - Standard Flowgram Format binary parsing

///|
test "sff_header_new" {
  let header = @src.SffHeader::new(2, 4, "TACG", "TCAG")
  assert_eq(header.n_reads, 2)
  assert_eq(header.n_flows, 4)
  assert_eq(header.flow_chars, "TACG")
  assert_eq(header.key_sequence, "TCAG")
  assert_eq(header.key_length, 4)
  // Header length should be padded to multiple of 8
  assert_eq(header.header_length % 8, 0)
}

///|
test "sff_header_methods" {
  let header = @src.SffHeader::new(3, 8, "TACGTACG", "TCAG")
  assert_eq(header.sff_n_reads(), 3)
  assert_eq(header.sff_n_flows(), 8)
  assert_eq(header.sff_flow_chars(), "TACGTACG")
  assert_eq(header.sff_key_sequence(), "TCAG")
}

///|
test "sff_header_validate" {
  let header = @src.SffHeader::new(1, 4, "TACG", "TCAG")
  let errors = header.sff_validate()
  assert_eq(errors.length(), 0)
}

///|
test "sff_read_new" {
  let qualities = [10, 20, 30, 25]
  let flowgram = [1.0, 0.0, 1.0, 0.0]
  let flow_index = [1, 3, 1, 3]
  let read = @src.SffRead::new("READ001", "ACGT", qualities, flowgram, flow_index, 1, 4, 0, 0)
  assert_eq(read.name, "READ001")
  assert_eq(read.bases, "ACGT")
  assert_eq(read.n_bases, 4)
  assert_eq(read.name_length, 7)
  // Read header length should be padded to multiple of 8
  assert_eq(read.read_header_length % 8, 0)
}

///|
test "sff_read_methods" {
  let read = @src.sff_create_example_read()
  assert_eq(read.sff_name(), "TEST001")
  assert_eq(read.sff_bases(), "ACGTACGT")
  assert_eq(read.sff_n_bases(), 8)
  let quals = read.sff_qualities()
  assert_eq(quals.length(), 8)
  assert_eq(quals[0], 20)
  let fg = read.sff_flowgram()
  assert_eq(fg.length(), 4)
  let fi = read.sff_flow_index()
  assert_eq(fi.length(), 8)
}

///|
test "sff_read_clipped" {
  let read = @src.sff_create_example_read()
  // clip_qual_left=1, clip_qual_right=8, so clipped is the full sequence
  let clipped = read.sff_clipped_bases()
  assert_eq(clipped, "ACGTACGT")
  let clipped_quals = read.sff_clipped_qualities()
  assert_eq(clipped_quals.length(), 8)
}

///|
test "sff_read_mean_quality" {
  let read = @src.sff_create_example_read()
  let mean_q = read.sff_mean_quality()
  // Qualities: [20, 25, 30, 28, 26, 24, 22, 20], mean = 195/8 = 24.375
  assert_true((mean_q - 24.375).abs() < 1.0e-9)
}

///|
test "sff_file_new" {
  let file = @src.sff_create_example_file()
  assert_eq(file.sff_n_reads(), 2)
  assert_eq(file.sff_total_bases(), 20) // 10 + 10
}

///|
test "sff_file_methods" {
  let file = @src.sff_create_example_file()
  let header = file.sff_header()
  assert_eq(header.sff_n_reads(), 2)
  let reads = file.sff_reads()
  assert_eq(reads.length(), 2)
  let r0 = file.sff_get_read(0)
  assert_eq(r0.sff_name(), "READ0001")
  let r1 = file.sff_get_read(1)
  assert_eq(r1.sff_name(), "READ0002")
}

///|
test "sff_file_find_read" {
  let file = @src.sff_create_example_file()
  let found = file.sff_find_read("READ0001")
  assert_true(found is Some(_))
  assert_eq(found.unwrap().sff_name(), "READ0001")
  let found2 = file.sff_find_read("READ0002")
  assert_true(found2 is Some(_))
  let missing = file.sff_find_read("MISSING")
  assert_true(missing is None)
}

///|
test "sff_encode_parse_roundtrip" {
  let file = @src.sff_create_example_file()
  let bytes = @src.sff_encode(file)
  // Bytes should be non-empty
  assert_true(bytes.length() > 100)
  // Parse it back
  let parsed_opt = @src.sff_parse(bytes)
  assert_true(parsed_opt is Some(_))
  let parsed = parsed_opt.unwrap()
  // Header should match
  assert_eq(parsed.sff_n_reads(), 2)
  let header = parsed.sff_header()
  assert_eq(header.flow_chars, "TACG")
  assert_eq(header.key_sequence, "TCAG")
  assert_eq(header.n_flows, 4)
  // Reads should match
  let r0 = parsed.sff_get_read(0)
  assert_eq(r0.sff_name(), "READ0001")
  assert_eq(r0.sff_bases(), "ACGTACGTAC")
  assert_eq(r0.sff_n_bases(), 10)
  let r1 = parsed.sff_get_read(1)
  assert_eq(r1.sff_name(), "READ0002")
  assert_eq(r1.sff_bases(), "TGCATGCATG")
}

///|
test "sff_parse_invalid_magic" {
  // Create bytes with wrong magic number
  let bytes = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 40, 0, 4, 0, 4, 0, 1]
  let result = @src.sff_parse(bytes)
  assert_true(result is None)
}

///|
test "sff_encode_u32" {
  let file = @src.sff_create_example_file()
  let bytes = @src.sff_encode(file)
  // First 4 bytes should be the magic number 0x2E736666
  // 0x2E = 46, 0x73 = 115, 0x66 = 102, 0x66 = 102
  assert_eq(bytes[0], 46)
  assert_eq(bytes[1], 115)
  assert_eq(bytes[2], 102)
  assert_eq(bytes[3], 102)
}

///|
test "sff_encode_parse_single_read" {
  let header = @src.SffHeader::new(1, 4, "TACG", "TCAG")
  let read = @src.sff_create_example_read()
  let file = @src.SffFile::new(header, [read])
  let bytes = @src.sff_encode(file)
  let parsed_opt = @src.sff_parse(bytes)
  assert_true(parsed_opt is Some(_))
  let parsed = parsed_opt.unwrap()
  assert_eq(parsed.sff_n_reads(), 1)
  let r0 = parsed.sff_get_read(0)
  assert_eq(r0.sff_name(), "TEST001")
  assert_eq(r0.sff_bases(), "ACGTACGT")
  // Qualities should round-trip
  let quals = r0.sff_qualities()
  assert_eq(quals[0], 20)
  assert_eq(quals[7], 20)
}

///|
test "sff_create_example_header" {
  let header = @src.sff_create_example_header()
  assert_eq(header.sff_n_reads(), 1)
  assert_eq(header.sff_n_flows(), 4)
  assert_eq(header.sff_flow_chars(), "TACG")
}

///|
test "sff_pad_to_8" {
  // Test via header_length: 31 + 4 + 4 = 39, padded to 40
  let header = @src.SffHeader::new(1, 4, "TACG", "TCAG")
  assert_eq(header.header_length, 40)
  // 31 + 8 + 4 = 43, padded to 48
  let header2 = @src.SffHeader::new(1, 8, "TACGTACG", "TCAG")
  assert_eq(header2.header_length, 48)
}
