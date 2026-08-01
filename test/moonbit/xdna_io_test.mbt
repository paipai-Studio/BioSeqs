///|
/// Tests for the XDNA I/O module (Bio.SeqIO.XdnaIO equivalent).
///
/// These tests exercise the public API exposed in `src/xdna_io.mbt`:
/// the XdnaSeqType enum, XdnaRecord / XdnaFile structs, the 2-bit DNA
/// packing helpers, the big-endian u32 helpers, the byte/hex serializers
/// (`xdna_to_bytes` / `xdna_from_bytes` / `xdna_write` / `xdna_read`),
/// the SeqRecord integration (`xdna_to_seq_records` / `xdna_from_seq_records`),
/// and the sample data (`xdna_sample_file`).

// ============================================================================
// 1. XdnaSeqType - to_int / from_int / to_string / description
// ============================================================================

///|
test "xdna_seq_type_to_int" {
  // XdnaSeqType constructors cannot be constructed directly from outside the
  // package; use from_int to obtain a value, then verify its to_int round-trip.
  assert_eq(@src.XdnaSeqType::from_int(0).to_int(), 0)
  assert_eq(@src.XdnaSeqType::from_int(1).to_int(), 1)
  assert_eq(@src.XdnaSeqType::from_int(2).to_int(), 2)
  assert_eq(@src.XdnaSeqType::from_int(3).to_int(), 3)
}

///|
test "xdna_seq_type_from_int" {
  assert_true(@src.XdnaSeqType::from_int(0) is @src.XdnaSeqType::DnaType)
  assert_true(@src.XdnaSeqType::from_int(1) is @src.XdnaSeqType::RnaType)
  assert_true(@src.XdnaSeqType::from_int(2) is @src.XdnaSeqType::ProteinType)
  assert_true(@src.XdnaSeqType::from_int(3) is @src.XdnaSeqType::UnknownType)
}

///|
test "xdna_seq_type_from_int_round_trip" {
  let dna = @src.XdnaSeqType::from_int(0)
  let rna = @src.XdnaSeqType::from_int(1)
  let protein = @src.XdnaSeqType::from_int(2)
  let unknown = @src.XdnaSeqType::from_int(3)
  assert_true(@src.XdnaSeqType::from_int(dna.to_int()) is @src.XdnaSeqType::DnaType)
  assert_true(@src.XdnaSeqType::from_int(rna.to_int()) is @src.XdnaSeqType::RnaType)
  assert_true(
    @src.XdnaSeqType::from_int(protein.to_int()) is @src.XdnaSeqType::ProteinType,
  )
  assert_true(
    @src.XdnaSeqType::from_int(unknown.to_int()) is @src.XdnaSeqType::UnknownType,
  )
}

///|
test "xdna_seq_type_from_int_invalid_maps_to_unknown" {
  // Any value outside 0-3 maps to UnknownType.
  assert_true(@src.XdnaSeqType::from_int(4) is @src.XdnaSeqType::UnknownType)
  assert_true(@src.XdnaSeqType::from_int(-1) is @src.XdnaSeqType::UnknownType)
  assert_true(@src.XdnaSeqType::from_int(99) is @src.XdnaSeqType::UnknownType)
}

///|
test "xdna_seq_type_to_string" {
  assert_eq(@src.XdnaSeqType::from_int(0).to_string(), "DnaType")
  assert_eq(@src.XdnaSeqType::from_int(1).to_string(), "RnaType")
  assert_eq(@src.XdnaSeqType::from_int(2).to_string(), "ProteinType")
  assert_eq(@src.XdnaSeqType::from_int(3).to_string(), "UnknownType")
}

///|
test "xdna_seq_type_description" {
  assert_eq(@src.XdnaSeqType::from_int(0).description(), "DNA")
  assert_eq(@src.XdnaSeqType::from_int(1).description(), "RNA")
  assert_eq(@src.XdnaSeqType::from_int(2).description(), "Protein")
  assert_eq(@src.XdnaSeqType::from_int(3).description(), "Unknown")
}

// ============================================================================
// 2. XdnaRecord - creation, accessors, to_string
// ============================================================================

///|
test "xdna_record_new_and_accessors" {
  let rec = @src.XdnaRecord::new(
    "seq1",
    "ATGCATGC",
    @src.XdnaSeqType::from_int(0),
    "annotation text",
  )
  assert_eq(rec.name(), "seq1")
  assert_eq(rec.sequence(), "ATGCATGC")
  assert_eq(rec.seq_type().to_int(), 0)
  assert_eq(rec.annotations(), "annotation text")
  // Default checksum is 0.
  assert_eq(rec.checksum(), 0)
}

///|
test "xdna_record_new_with_checksum" {
  let rec = @src.XdnaRecord::new(
    "seq2",
    "GATTACA",
    @src.XdnaSeqType::from_int(0),
    "",
    checksum=12345,
  )
  assert_eq(rec.name(), "seq2")
  assert_eq(rec.sequence(), "GATTACA")
  assert_eq(rec.checksum(), 12345)
}

///|
test "xdna_record_new_rna_type" {
  let rec = @src.XdnaRecord::new(
    "rna1",
    "AUGCAUGC",
    @src.XdnaSeqType::from_int(1),
    "",
  )
  assert_true(rec.seq_type() is @src.XdnaSeqType::RnaType)
  assert_eq(rec.sequence(), "AUGCAUGC")
}

///|
test "xdna_record_new_protein_type" {
  let rec = @src.XdnaRecord::new(
    "prot1",
    "MKLVGV",
    @src.XdnaSeqType::from_int(2),
    "a protein",
  )
  assert_true(rec.seq_type() is @src.XdnaSeqType::ProteinType)
  assert_eq(rec.annotations(), "a protein")
}

///|
test "xdna_record_to_string" {
  let rec = @src.XdnaRecord::new(
    "seq1",
    "ATGCATGC",
    @src.XdnaSeqType::from_int(0),
    "note",
    checksum=42,
  )
  let s = rec.to_string()
  assert_true(s.contains("XdnaRecord(name=seq1"))
  assert_true(s.contains("seq_type=DnaType"))
  assert_true(s.contains("length=8"))
  assert_true(s.contains("annotations_length=4"))
  assert_true(s.contains("checksum=42"))
}

// ============================================================================
// 3. XdnaFile - creation, records, version, n_records
// ============================================================================

///|
test "xdna_file_new_is_empty" {
  let file = @src.XdnaFile::new()
  assert_eq(file.n_records(), 0)
  assert_eq(file.records().length(), 0)
  assert_eq(file.version(), 1)
}

///|
test "xdna_file_add_record_and_n_records" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("a", "ACGT", @src.XdnaSeqType::from_int(0), ""))
  assert_eq(file.n_records(), 1)
  file.add_record(@src.XdnaRecord::new("b", "TTTT", @src.XdnaSeqType::from_int(0), ""))
  assert_eq(file.n_records(), 2)
}

///|
test "xdna_file_records_returns_copy" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("a", "ACGT", @src.XdnaSeqType::from_int(0), ""))
  let recs = file.records()
  // Mutating the returned array should not affect the file.
  recs.push(@src.XdnaRecord::new("x", "GGGG", @src.XdnaSeqType::from_int(0), ""))
  assert_eq(file.n_records(), 1)
  assert_eq(recs.length(), 2)
}

///|
test "xdna_file_records_access" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("seq1", "ATGC", @src.XdnaSeqType::from_int(0), ""))
  file.add_record(@src.XdnaRecord::new("seq2", "GATTACA", @src.XdnaSeqType::from_int(0), ""))
  let recs = file.records()
  assert_eq(recs[0].name(), "seq1")
  assert_eq(recs[0].sequence(), "ATGC")
  assert_eq(recs[1].name(), "seq2")
  assert_eq(recs[1].sequence(), "GATTACA")
}

// ============================================================================
// 4. Base encoding - xdna_base_to_code / xdna_code_to_base
// ============================================================================

///|
test "xdna_base_to_code_valid" {
  assert_eq(@src.xdna_base_to_code('A'), 0)
  assert_eq(@src.xdna_base_to_code('C'), 1)
  assert_eq(@src.xdna_base_to_code('G'), 2)
  assert_eq(@src.xdna_base_to_code('T'), 3)
}

///|
test "xdna_base_to_code_invalid" {
  assert_eq(@src.xdna_base_to_code('N'), -1)
  assert_eq(@src.xdna_base_to_code('X'), -1)
  assert_eq(@src.xdna_base_to_code('Z'), -1)
  assert_eq(@src.xdna_base_to_code('5'), -1)
  // Lowercase is not handled.
  assert_eq(@src.xdna_base_to_code('a'), -1)
}

///|
test "xdna_code_to_base_valid" {
  assert_eq(@src.xdna_code_to_base(0), 'A')
  assert_eq(@src.xdna_code_to_base(1), 'C')
  assert_eq(@src.xdna_code_to_base(2), 'G')
  assert_eq(@src.xdna_code_to_base(3), 'T')
}

///|
test "xdna_code_to_base_invalid_returns_N" {
  assert_eq(@src.xdna_code_to_base(4), 'N')
  assert_eq(@src.xdna_code_to_base(-1), 'N')
  assert_eq(@src.xdna_code_to_base(99), 'N')
}

// ============================================================================
// 5. Pack / unpack sequences
// ============================================================================

///|
test "xdna_pack_sequence_simple" {
  // ATGC: A=0, T=3, G=2, C=1 -> 0b00_11_10_01 = 57
  assert_eq(@src.xdna_pack_sequence("ATGC"), [57])
}

///|
test "xdna_pack_sequence_odd_length" {
  // ATGCA (5 bases): 2 bytes. Byte 0 = 57 (ATGC), byte 1 = 0 (A + padding).
  assert_eq(@src.xdna_pack_sequence("ATGCA"), [57, 0])
}

///|
test "xdna_pack_sequence_empty" {
  assert_eq(@src.xdna_pack_sequence("").length(), 0)
}

///|
test "xdna_unpack_sequence_simple" {
  assert_eq(@src.xdna_unpack_sequence([57], 4), "ATGC")
}

///|
test "xdna_unpack_sequence_odd_length" {
  assert_eq(@src.xdna_unpack_sequence([57, 0], 5), "ATGCA")
}

///|
test "xdna_unpack_sequence_empty" {
  assert_eq(@src.xdna_unpack_sequence([], 0), "")
}

///|
test "xdna_pack_unpack_round_trip_simple" {
  let seq = "ATGC"
  let packed = @src.xdna_pack_sequence(seq)
  let unpacked = @src.xdna_unpack_sequence(packed, seq.length())
  assert_eq(unpacked, seq)
}

///|
test "xdna_pack_unpack_round_trip_odd" {
  let seq = "ATGCA"
  let packed = @src.xdna_pack_sequence(seq)
  let unpacked = @src.xdna_unpack_sequence(packed, seq.length())
  assert_eq(unpacked, seq)
}

///|
test "xdna_pack_unpack_round_trip_long" {
  let seq = "ATGCATGCATGCATGCATGC"
  let packed = @src.xdna_pack_sequence(seq)
  assert_eq(packed.length(), 5)
  let unpacked = @src.xdna_unpack_sequence(packed, seq.length())
  assert_eq(unpacked, seq)
}

///|
test "xdna_pack_unpack_round_trip_single_base" {
  let seq = "A"
  let packed = @src.xdna_pack_sequence(seq)
  assert_eq(packed.length(), 1)
  let unpacked = @src.xdna_unpack_sequence(packed, seq.length())
  assert_eq(unpacked, seq)
}

///|
test "xdna_pack_unpack_all_bases" {
  // Verify each of A/C/G/T packs and unpacks correctly.
  for base in ["A", "C", "G", "T"] {
    let packed = @src.xdna_pack_sequence(base)
    let unpacked = @src.xdna_unpack_sequence(packed, 1)
    assert_eq(unpacked, base)
  }
}

// ============================================================================
// 6. Big-endian u32 helpers
// ============================================================================

///|
test "xdna_write_u32_be_known_value" {
  // 0x12345678 -> [0x12, 0x34, 0x56, 0x78]
  assert_eq(@src.xdna_write_u32_be(0x12345678), [0x12, 0x34, 0x56, 0x78])
}

///|
test "xdna_write_u32_be_zero" {
  assert_eq(@src.xdna_write_u32_be(0), [0, 0, 0, 0])
}

///|
test "xdna_read_u32_be_known_value" {
  assert_eq(@src.xdna_read_u32_be([0x12, 0x34, 0x56, 0x78], 0), 0x12345678)
}

///|
test "xdna_read_write_u32_be_round_trip" {
  let values = [0, 1, 255, 256, 65535, 65536, 0x12345678, 0xFFFFFFFF]
  for v in values {
    let bytes = @src.xdna_write_u32_be(v)
    let read_back = @src.xdna_read_u32_be(bytes, 0)
    assert_eq(read_back, v)
  }
}

///|
test "xdna_read_u32_be_at_offset" {
  // Reading at a non-zero offset should work.
  let data = [0, 0, 0, 0, 0x12, 0x34, 0x56, 0x78]
  assert_eq(@src.xdna_read_u32_be(data, 4), 0x12345678)
}

// ============================================================================
// 7. Serialization round-trip: xdna_to_bytes -> xdna_from_bytes
// ============================================================================

///|
test "xdna_to_bytes_from_bytes_round_trip_dna" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("seq1", "ATGCATGC", @src.XdnaSeqType::from_int(0), ""))
  file.add_record(@src.XdnaRecord::new("seq2", "GATTACA", @src.XdnaSeqType::from_int(0), ""))
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 2)
  let recs = file2.records()
  assert_eq(recs[0].name(), "seq1")
  assert_eq(recs[0].sequence(), "ATGCATGC")
  assert_true(recs[0].seq_type() is @src.XdnaSeqType::DnaType)
  assert_eq(recs[1].name(), "seq2")
  assert_eq(recs[1].sequence(), "GATTACA")
}

///|
test "xdna_to_bytes_from_bytes_round_trip_with_annotations" {
  let file = @src.XdnaFile::new()
  file.add_record(
    @src.XdnaRecord::new("annot1", "ATGC", @src.XdnaSeqType::from_int(0), "some annotation"),
  )
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 1)
  let rec = file2.records()[0]
  assert_eq(rec.name(), "annot1")
  assert_eq(rec.sequence(), "ATGC")
  assert_eq(rec.annotations(), "some annotation")
}

///|
test "xdna_to_bytes_from_bytes_round_trip_rna" {
  // RNA sequences are stored as raw ASCII (not 2-bit packed).
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("rna1", "AUGCAUGC", @src.XdnaSeqType::from_int(1), ""))
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 1)
  let rec = file2.records()[0]
  assert_eq(rec.name(), "rna1")
  assert_eq(rec.sequence(), "AUGCAUGC")
  assert_true(rec.seq_type() is @src.XdnaSeqType::RnaType)
}

///|
test "xdna_to_bytes_from_bytes_round_trip_protein" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("prot1", "MKLVGV", @src.XdnaSeqType::from_int(2), ""))
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 1)
  let rec = file2.records()[0]
  assert_eq(rec.name(), "prot1")
  assert_eq(rec.sequence(), "MKLVGV")
  assert_true(rec.seq_type() is @src.XdnaSeqType::ProteinType)
}

///|
test "xdna_to_bytes_includes_checksum" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("seq1", "ATGC", @src.XdnaSeqType::from_int(0), ""))
  let bytes = @src.xdna_to_bytes(file)
  // The last 4 bytes are the checksum (big-endian u32).
  let checksum = @src.xdna_read_u32_be(bytes, bytes.length() - 4)
  // The checksum is the sum of all preceding bytes; for a non-empty file it
  // should be strictly positive.
  assert_true(checksum > 0)
  // After round-trip, each record's checksum field is set to the stored value.
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.records()[0].checksum(), checksum)
}

///|
test "xdna_from_bytes_preserves_version" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("seq1", "ATGC", @src.XdnaSeqType::from_int(0), ""))
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.version(), 1)
}

// ============================================================================
// 8. Hex string I/O round-trip: xdna_write -> xdna_read
// ============================================================================

///|
test "xdna_write_produces_lowercase_hex" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("A", "ATGC", @src.XdnaSeqType::from_int(0), ""))
  let hex = @src.xdna_write(file)
  // Every character must be a lowercase hex digit.
  for c in hex {
    let code = c.to_int()
    let is_digit = code >= '0'.to_int() && code <= '9'.to_int()
    let is_lower = code >= 'a'.to_int() && code <= 'f'.to_int()
    assert_true(is_digit || is_lower)
  }
}

///|
test "xdna_write_read_round_trip" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("seq1", "ATGCATGC", @src.XdnaSeqType::from_int(0), ""))
  file.add_record(@src.XdnaRecord::new("seq2", "GATTACA", @src.XdnaSeqType::from_int(0), ""))
  let hex = @src.xdna_write(file)
  let file2 = @src.xdna_read(hex)
  assert_eq(file2.n_records(), 2)
  let recs = file2.records()
  assert_eq(recs[0].name(), "seq1")
  assert_eq(recs[0].sequence(), "ATGCATGC")
  assert_eq(recs[1].name(), "seq2")
  assert_eq(recs[1].sequence(), "GATTACA")
}

///|
test "xdna_write_read_round_trip_with_annotations" {
  let file = @src.XdnaFile::new()
  file.add_record(
    @src.XdnaRecord::new("a1", "GGGG", @src.XdnaSeqType::from_int(0), "note here"),
  )
  let hex = @src.xdna_write(file)
  let file2 = @src.xdna_read(hex)
  let rec = file2.records()[0]
  assert_eq(rec.name(), "a1")
  assert_eq(rec.sequence(), "GGGG")
  assert_eq(rec.annotations(), "note here")
}

// ============================================================================
// 9. xdna_get_sequence_string
// ============================================================================

///|
test "xdna_get_sequence_string" {
  let rec = @src.XdnaRecord::new("seq1", "ATGCATGC", @src.XdnaSeqType::from_int(0), "")
  assert_eq(@src.xdna_get_sequence_string(rec), "ATGCATGC")
}

///|
test "xdna_get_sequence_string_empty" {
  let rec = @src.XdnaRecord::new("empty", "", @src.XdnaSeqType::from_int(0), "")
  assert_eq(@src.xdna_get_sequence_string(rec), "")
}

// ============================================================================
// 10. SeqRecord conversion - xdna_to_seq_records / xdna_from_seq_records
// ============================================================================

///|
test "xdna_to_seq_records" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("seq1", "ATGCATGC", @src.XdnaSeqType::from_int(0), ""))
  file.add_record(@src.XdnaRecord::new("seq2", "GATTACA", @src.XdnaSeqType::from_int(0), ""))
  let records = @src.xdna_to_seq_records(file)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].name, "seq1")
  assert_eq(records[0].description, "seq1")
  assert_eq(records[0].seq.to_string(), "ATGCATGC")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "GATTACA")
}

///|
test "xdna_from_seq_records" {
  let records = [
    @src.SeqRecord::new(@src.Seq::new("ATGCATGC"), id="seq1"),
    @src.SeqRecord::new(@src.Seq::new("GATTACA"), id="seq2"),
  ]
  let file = @src.xdna_from_seq_records(records, @src.XdnaSeqType::from_int(0))
  assert_eq(file.n_records(), 2)
  let recs = file.records()
  assert_eq(recs[0].name(), "seq1")
  assert_eq(recs[0].sequence(), "ATGCATGC")
  assert_true(recs[0].seq_type() is @src.XdnaSeqType::DnaType)
  assert_eq(recs[1].name(), "seq2")
  assert_eq(recs[1].sequence(), "GATTACA")
}

///|
test "xdna_seq_records_round_trip" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("seq1", "ATGC", @src.XdnaSeqType::from_int(0), ""))
  let records = @src.xdna_to_seq_records(file)
  let file2 = @src.xdna_from_seq_records(records, @src.XdnaSeqType::from_int(0))
  assert_eq(file2.n_records(), 1)
  assert_eq(file2.records()[0].name(), "seq1")
  assert_eq(file2.records()[0].sequence(), "ATGC")
}

// ============================================================================
// 11. Sample file validation
// ============================================================================

///|
test "xdna_sample_file_structure" {
  let file = @src.xdna_sample_file()
  assert_eq(file.n_records(), 2)
  assert_eq(file.version(), 1)
}

///|
test "xdna_sample_file_records" {
  let file = @src.xdna_sample_file()
  let recs = file.records()
  assert_eq(recs[0].name(), "seq1")
  assert_eq(recs[0].sequence(), "ATGCATGC")
  assert_true(recs[0].seq_type() is @src.XdnaSeqType::DnaType)
  assert_eq(recs[1].name(), "seq2")
  assert_eq(recs[1].sequence(), "GATTACA")
  assert_true(recs[1].seq_type() is @src.XdnaSeqType::DnaType)
}

///|
test "xdna_sample_file_round_trip_bytes" {
  let file = @src.xdna_sample_file()
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 2)
  let recs = file2.records()
  assert_eq(recs[0].name(), "seq1")
  assert_eq(recs[0].sequence(), "ATGCATGC")
  assert_eq(recs[1].name(), "seq2")
  assert_eq(recs[1].sequence(), "GATTACA")
}

///|
test "xdna_sample_file_round_trip_hex" {
  let file = @src.xdna_sample_file()
  let hex = @src.xdna_write(file)
  let file2 = @src.xdna_read(hex)
  assert_eq(file2.n_records(), 2)
  let recs = file2.records()
  assert_eq(recs[0].name(), "seq1")
  assert_eq(recs[0].sequence(), "ATGCATGC")
  assert_eq(recs[1].name(), "seq2")
  assert_eq(recs[1].sequence(), "GATTACA")
}

// ============================================================================
// 12. Edge cases - empty file, single record, single base
// ============================================================================

///|
test "edge_case_empty_file_to_bytes" {
  let file = @src.XdnaFile::new()
  let bytes = @src.xdna_to_bytes(file)
  // An empty file still has the 4-byte trailing checksum (which is 0).
  assert_eq(bytes.length(), 4)
  assert_eq(@src.xdna_read_u32_be(bytes, 0), 0)
}

///|
test "edge_case_empty_file_round_trip" {
  let file = @src.XdnaFile::new()
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 0)
}

///|
test "edge_case_empty_file_hex_round_trip" {
  let file = @src.XdnaFile::new()
  let hex = @src.xdna_write(file)
  let file2 = @src.xdna_read(hex)
  assert_eq(file2.n_records(), 0)
}

///|
test "edge_case_single_record" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("only", "ACGTACGT", @src.XdnaSeqType::from_int(0), ""))
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 1)
  let rec = file2.records()[0]
  assert_eq(rec.name(), "only")
  assert_eq(rec.sequence(), "ACGTACGT")
}

///|
test "edge_case_single_base" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("one", "A", @src.XdnaSeqType::from_int(0), ""))
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 1)
  let rec = file2.records()[0]
  assert_eq(rec.name(), "one")
  assert_eq(rec.sequence(), "A")
}

///|
test "edge_case_single_base_hex_round_trip" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("one", "G", @src.XdnaSeqType::from_int(0), ""))
  let hex = @src.xdna_write(file)
  let file2 = @src.xdna_read(hex)
  assert_eq(file2.n_records(), 1)
  assert_eq(file2.records()[0].sequence(), "G")
}

///|
test "edge_case_empty_sequence" {
  let file = @src.XdnaFile::new()
  file.add_record(@src.XdnaRecord::new("empty", "", @src.XdnaSeqType::from_int(0), ""))
  let bytes = @src.xdna_to_bytes(file)
  let file2 = @src.xdna_from_bytes(bytes)
  assert_eq(file2.n_records(), 1)
  let rec = file2.records()[0]
  assert_eq(rec.name(), "empty")
  assert_eq(rec.sequence(), "")
}

///|
test "edge_case_to_seq_records_empty_file" {
  let file = @src.XdnaFile::new()
  let records = @src.xdna_to_seq_records(file)
  assert_eq(records.length(), 0)
}

///|
test "edge_case_from_seq_records_empty" {
  let file = @src.xdna_from_seq_records([], @src.XdnaSeqType::from_int(0))
  assert_eq(file.n_records(), 0)
}
