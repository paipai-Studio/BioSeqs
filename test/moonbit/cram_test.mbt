///|
/// Tests for CRAM ITF-8 encoding, binary write round-trip, and SAM conversion.

// ============================================================
// ITF-8 encode/decode tests
// ============================================================

test "itf8 single byte" {
  // 0..127 → 1 byte
  let buf : Array[Int] = []
  @src.itf8_encode(0, buf)
  @src.itf8_encode(127, buf)
  assert_eq(buf.length(), 2)
  assert_eq(buf[0], 0)
  assert_eq(buf[1], 127)

  // decode back
  
  let s2 = StringBuilder()
  s2.write_char(buf[0].to_uint16().unsafe_to_char())
  s2.write_char(buf[1].to_uint16().unsafe_to_char())
  let data = s2.to_string()
  let (v0, p0) = @src.itf8_decode(data, 0)
  assert_eq(v0, 0)
  assert_eq(p0, 1)
  let (v1, p1) = @src.itf8_decode(data, 1)
  assert_eq(v1, 127)
  assert_eq(p1, 2)
}

test "itf8 two bytes" {
  // 128..16383 → 2 bytes
  let buf : Array[Int] = []
  @src.itf8_encode(128, buf)
  @src.itf8_encode(16383, buf)
  assert_eq(buf.length(), 4) // 2 + 2

  let sb = StringBuilder()
  for i = 0; i < buf.length(); i = i + 1 {
    sb.write_char(buf[i].to_uint16().unsafe_to_char())
  }
  let data = sb.to_string()
  let (v0, p0) = @src.itf8_decode(data, 0)
  assert_eq(v0, 128)
  let (v1, p1) = @src.itf8_decode(data, p0)
  assert_eq(v1, 16383)
}

test "itf8 three bytes" {
  // 16384..2097151 → 3 bytes
  let buf : Array[Int] = []
  @src.itf8_encode(16384, buf)
  @src.itf8_encode(0x1FFFF, buf)
  assert_eq(buf.length(), 6) // 3 + 3

  let sb = StringBuilder()
  for i = 0; i < buf.length(); i = i + 1 {
    sb.write_char(buf[i].to_uint16().unsafe_to_char())
  }
  let data = sb.to_string()
  let (v0, p0) = @src.itf8_decode(data, 0)
  assert_eq(v0, 16384)
  let (v1, _) = @src.itf8_decode(data, p0)
  assert_eq(v1, 0x1FFFF)
}

test "itf8 four bytes" {
  // 2097152..268435455 → 4 bytes
  let buf : Array[Int] = []
  @src.itf8_encode(0x1000000, buf)
  assert_eq(buf.length(), 4)

  let sb = StringBuilder()
  for i = 0; i < buf.length(); i = i + 1 {
    sb.write_char(buf[i].to_uint16().unsafe_to_char())
  }
  let data = sb.to_string()
  let (v, _) = @src.itf8_decode(data, 0)
  assert_eq(v, 0x1000000)
}

test "itf8 five bytes" {
  // > 268435455 → 5 bytes
  let buf : Array[Int] = []
  @src.itf8_encode(0x10000000, buf)
  assert_eq(buf.length(), 5)

  let sb = StringBuilder()
  for i = 0; i < buf.length(); i = i + 1 {
    sb.write_char(buf[i].to_uint16().unsafe_to_char())
  }
  let data = sb.to_string()
  let (v, _) = @src.itf8_decode(data, 0)
  assert_eq(v, 0x10000000)
}

// ============================================================
// CRAM round-trip: parse → write → parse
// ============================================================

test "cram roundtrip sample" {
  // Parse the sample CRAM content
  let original = @src.sample_cram_content()
  let file1 = @src.parse_cram(original)
  // Write it back
  let bytes = file1.to_cram_bytes()
  // Re-parse
  let file2 = @src.parse_cram(bytes)

  // Compare header
  assert_eq(file2.header.version.major, file1.header.version.major)
  assert_eq(file2.header.version.minor, file1.header.version.minor)
  assert_eq(file2.header.references.length(), file1.header.references.length())
  assert_eq(file2.header.references[0].name, "chr1")
  assert_eq(file2.header.references[0].length, 100)
  assert_eq(file2.header.sam_header, file1.header.sam_header)

  // Compare records
  assert_eq(file2.records.length(), file1.records.length())
  assert_eq(file2.records.length(), 1)
  let r1 = file1.records[0]
  let r2 = file2.records[0]
  assert_eq(r2.ref_id, r1.ref_id)
  assert_eq(r2.pos, r1.pos)
  assert_eq(r2.read_features.read_name, "read1")
  assert_eq(r2.read_features.mapping_quality, 60)
  assert_eq(r2.read_features.sequence, "ATGCATGC")
  assert_eq(r2.read_features.quality, "!!!!!!!!")
  assert_eq(r2.read_features.cigar.length(), 1)
  assert_eq(r2.read_features.cigar[0].op is @src.BamCigarOp::Match, true)
  assert_eq(r2.read_features.cigar[0].length, 100)
}

// ============================================================
// CRAM → SAM conversion
// ============================================================

test "cram to_sam" {
  let file = @src.sample_cram_file().unwrap()
  let sam = file.to_sam()

  // Should contain SAM header
  assert_true(sam.contains("@SQ"))
  assert_true(sam.contains("chr1"))

  // Should contain the read
  assert_true(sam.contains("read1"))
  // 1-based position: pos=1 → 2
  assert_true(sam.contains("\t2\t"))
  // mapq=60
  assert_true(sam.contains("\t60\t"))
  // cigar 100M
  assert_true(sam.contains("100M"))
  // sequence
  assert_true(sam.contains("ATGCATGC"))
}

test "cram record flag helpers" {
  let file = @src.sample_cram_file().unwrap()
  let rec = file.records[0]
  // flags=0 → not paired, not unmapped, not reverse
  assert_false(rec.is_paired())
  assert_false(rec.is_unmapped())
  assert_false(rec.is_reverse())
}