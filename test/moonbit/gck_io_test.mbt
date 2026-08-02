///|
/// Tests for the GCK I/O module (Bio.SeqIO.GckIO equivalent).
///
/// These tests exercise the public API exposed in `src/gck_io.mbt`:
/// the GckSeqType enum, GckFeature / GckFile structs, the hex/byte
/// serializers (`gck_parse` / `gck_write` / `gck_from_bytes` /
/// `gck_to_bytes`), the SeqRecord integration (`gck_to_seqrecord` /
/// `gck_from_seqrecord`), and utility functions (`gck_count_features_by_type`,
/// `gck_seq_length`, `gck_summary`, `gck_sample`).

// ============================================================================
// 1. GckSeqType - to_int / from_int / to_string
// ============================================================================

///|
test "gck_seq_type_to_int" {
  assert_eq(@src.GckSeqType::from_int(0).to_int(), 0)
  assert_eq(@src.GckSeqType::from_int(1).to_int(), 1)
  assert_eq(@src.GckSeqType::from_int(2).to_int(), 2)
  assert_eq(@src.GckSeqType::from_int(3).to_int(), 3)
}

///|
test "gck_seq_type_from_int" {
  assert_true(@src.GckSeqType::from_int(0) is @src.GckSeqType::GckDna)
  assert_true(@src.GckSeqType::from_int(1) is @src.GckSeqType::GckRna)
  assert_true(@src.GckSeqType::from_int(2) is @src.GckSeqType::GckProtein)
  assert_true(@src.GckSeqType::from_int(3) is @src.GckSeqType::GckUnknown)
}

///|
test "gck_seq_type_from_int_round_trip" {
  let dna = @src.GckSeqType::from_int(0)
  let rna = @src.GckSeqType::from_int(1)
  let protein = @src.GckSeqType::from_int(2)
  assert_true(@src.GckSeqType::from_int(dna.to_int()) is @src.GckSeqType::GckDna)
  assert_true(@src.GckSeqType::from_int(rna.to_int()) is @src.GckSeqType::GckRna)
  assert_true(
    @src.GckSeqType::from_int(protein.to_int()) is @src.GckSeqType::GckProtein,
  )
}

///|
test "gck_seq_type_from_int_invalid_maps_to_unknown" {
  assert_true(@src.GckSeqType::from_int(99) is @src.GckSeqType::GckUnknown)
  assert_true(@src.GckSeqType::from_int(-1) is @src.GckSeqType::GckUnknown)
}

///|
test "gck_seq_type_to_string" {
  assert_eq(@src.GckSeqType::from_int(0).to_string(), "DNA")
  assert_eq(@src.GckSeqType::from_int(1).to_string(), "RNA")
  assert_eq(@src.GckSeqType::from_int(2).to_string(), "Protein")
  assert_eq(@src.GckSeqType::from_int(3).to_string(), "Unknown")
}

// ============================================================================
// 2. GckFile - construction and defaults
// ============================================================================

///|
test "gck_file_new_defaults" {
  let file = @src.GckFile::new()
  assert_eq(@src.gck_seq_length(file), 0)
  assert_eq(file.features.length(), 0)
  assert_eq(file.sequence, "")
  assert_false(file.is_circular)
  assert_eq(file.version, 1)
}

// ============================================================================
// 3. GckFeature - construction
// ============================================================================

///|
test "gck_feature_new" {
  let f = @src.GckFeature::new(
    name="AmpR",
    type_="CDS",
    direction=1,
    segments=[(0, 19)],
  )
  assert_eq(f.name, "AmpR")
  assert_eq(f.type_, "CDS")
  assert_eq(f.direction, 1)
  assert_eq(f.segments.length(), 1)
  let (s, e) = f.segments[0]
  assert_eq(s, 0)
  assert_eq(e, 19)
}

///|
test "gck_feature_new_multiple_segments" {
  let f = @src.GckFeature::new(
    name="gene1",
    type_="CDS",
    direction=2,
    segments=[(0, 10), (20, 30), (40, 50)],
  )
  assert_eq(f.name, "gene1")
  assert_eq(f.direction, 2)
  assert_eq(f.segments.length(), 3)
}

// ============================================================================
// 4. gck_sample - sample data
// ============================================================================

///|
test "gck_sample_structure" {
  let file = @src.gck_sample()
  assert_true(@src.gck_seq_length(file) > 0)
  assert_true(file.features.length() >= 2)
  assert_true(file.is_circular)
}

///|
test "gck_sample_feature_types" {
  let file = @src.gck_sample()
  assert_true(@src.gck_count_features_by_type(file, "CDS") >= 1)
}

// ============================================================================
// 5. gck_count_features_by_type
// ============================================================================

///|
test "gck_count_features_by_type_empty" {
  let file = @src.GckFile::new()
  assert_eq(@src.gck_count_features_by_type(file, "CDS"), 0)
}

///|
test "gck_count_features_by_type_multiple" {
  let file = @src.GckFile::new()
  file.add_feature(
    @src.GckFeature::new(name="a", type_="CDS", direction=1, segments=[(0, 10)]),
  )
  file.add_feature(
    @src.GckFeature::new(name="b", type_="CDS", direction=1, segments=[(20, 30)]),
  )
  file.add_feature(
    @src.GckFeature::new(
      name="c",
      type_="promoter",
      direction=1,
      segments=[(40, 50)],
    ),
  )
  assert_eq(@src.gck_count_features_by_type(file, "CDS"), 2)
  assert_eq(@src.gck_count_features_by_type(file, "promoter"), 1)
  assert_eq(@src.gck_count_features_by_type(file, "terminator"), 0)
}

// ============================================================================
// 6. gck_seq_length
// ============================================================================

///|
test "gck_seq_length_empty" {
  let file = @src.GckFile::new()
  assert_eq(@src.gck_seq_length(file), 0)
}

///|
test "gck_seq_length_nonempty" {
  let file = @src.GckFile::new()
  file.set_sequence("ATGGCTAGCAA")
  assert_eq(@src.gck_seq_length(file), 11)
}

// ============================================================================
// 7. gck_to_bytes / gck_from_bytes round trip
// ============================================================================

///|
test "gck_to_bytes_from_bytes_round_trip" {
  let file = @src.gck_sample()
  let bytes = @src.gck_to_bytes(file)
  let parsed = @src.gck_from_bytes(bytes)
  assert_eq(parsed.sequence, file.sequence)
  assert_eq(parsed.is_circular, file.is_circular)
  assert_eq(parsed.features.length(), file.features.length())
  assert_eq(parsed.version, file.version)
}

///|
test "gck_write_parse_round_trip" {
  let file = @src.gck_sample()
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.sequence, file.sequence)
  assert_eq(parsed.is_circular, file.is_circular)
  assert_eq(parsed.features.length(), file.features.length())
}

///|
test "gck_to_bytes_starts_with_magic" {
  let file = @src.GckFile::new()
  file.set_sequence("ATGC")
  let bytes = @src.gck_to_bytes(file)
  // Magic "GCK_FILE" = 8 bytes
  assert_true(bytes.length() >= 8)
  // First byte should be 'G' (71)
  assert_eq(bytes[0], 71)
  // Second byte should be 'C' (67)
  assert_eq(bytes[1], 67)
  // Third byte should be 'K' (75)
  assert_eq(bytes[2], 75)
}

///|
test "gck_parse_invalid_magic_returns_empty" {
  // Hex for "XXXXXXXX" + minimal data — should not crash, returns empty file
  let bad_hex = "5858585858585858"
  let parsed = @src.gck_parse(bad_hex)
  assert_eq(parsed.sequence, "")
  assert_eq(parsed.features.length(), 0)
}

///|
test "gck_parse_empty_hex" {
  let parsed = @src.gck_parse("")
  assert_eq(parsed.sequence, "")
}

///|
test "gck_parse_preserves_features" {
  let file = @src.gck_sample()
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.features.length(), file.features.length())
  // Feature names and types should round-trip
  assert_eq(parsed.features[0].name, file.features[0].name)
  assert_eq(parsed.features[0].type_, file.features[0].type_)
}

///|
test "gck_parse_preserves_seq_type" {
  let file = @src.gck_sample()
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.seq_type.to_int(), file.seq_type.to_int())
}

///|
test "gck_parse_preserves_topology" {
  let file = @src.GckFile::new()
  file.set_sequence("ATGCATGC")
  file.set_circular(true)
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_true(parsed.is_circular)
}

///|
test "gck_parse_preserves_linear_topology" {
  let file = @src.GckFile::new()
  file.set_sequence("ATGCATGC")
  file.set_circular(false)
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_false(parsed.is_circular)
}

// ============================================================================
// 8. gck_to_seqrecord / gck_from_seqrecord
// ============================================================================

///|
test "gck_to_seqrecord" {
  let file = @src.gck_sample()
  let rec = @src.gck_to_seqrecord(file)
  assert_eq(rec.seq.data, file.sequence)
}

///|
test "gck_from_seqrecord" {
  let rec = @src.SeqRecord::new(@src.Seq::new("ATGGCTAGCAA"), id="test")
  let file = @src.gck_from_seqrecord(rec)
  assert_eq(file.sequence, "ATGGCTAGCAA")
}

///|
test "gck_seqrecord_round_trip" {
  let file = @src.gck_sample()
  let rec = @src.gck_to_seqrecord(file)
  let file2 = @src.gck_from_seqrecord(rec)
  assert_eq(file2.sequence, file.sequence)
}

// ============================================================================
// 9. gck_summary
// ============================================================================

///|
test "gck_summary_contains_key_info" {
  let file = @src.gck_sample()
  let summary = @src.gck_summary(file)
  assert_true(summary.contains("GCK File Summary:"))
  assert_true(summary.contains("Sequence Length:"))
  assert_true(summary.contains("Features:"))
}

///|
test "gck_summary_empty_file" {
  let file = @src.GckFile::new()
  let summary = @src.gck_summary(file)
  assert_true(summary.contains("GCK File Summary:"))
  assert_true(summary.contains("Features: 0"))
}

// ============================================================================
// 10. Custom file construction and serialization
// ============================================================================

///|
test "gck_custom_file_round_trip" {
  let file = @src.GckFile::new()
  file.set_sequence("ATGAAATAG")
  file.set_circular(false)
  file.add_feature(
    @src.GckFeature::new(name="gene1", type_="CDS", direction=1, segments=[(0, 8)]),
  )
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.sequence, "ATGAAATAG")
  assert_false(parsed.is_circular)
  assert_eq(parsed.features.length(), 1)
  assert_eq(parsed.features[0].name, "gene1")
}

///|
test "gck_multiple_features_round_trip" {
  let file = @src.GckFile::new()
  file.set_sequence("ATGCATGCATGCATGC")
  file.set_circular(true)
  file.add_feature(
    @src.GckFeature::new(name="f1", type_="CDS", direction=1, segments=[(0, 3)]),
  )
  file.add_feature(
    @src.GckFeature::new(
      name="f2",
      type_="promoter",
      direction=2,
      segments=[(4, 7)],
    ),
  )
  file.add_feature(
    @src.GckFeature::new(
      name="f3",
      type_="terminator",
      direction=0,
      segments=[(8, 11)],
    ),
  )
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.features.length(), 3)
  assert_eq(parsed.features[0].name, "f1")
  assert_eq(parsed.features[1].name, "f2")
  assert_eq(parsed.features[2].name, "f3")
  assert_eq(parsed.features[0].direction, 1)
  assert_eq(parsed.features[1].direction, 2)
  assert_eq(parsed.features[2].direction, 0)
}

///|
test "gck_multi_segment_feature_round_trip" {
  let file = @src.GckFile::new()
  file.set_sequence("ATGCATGCATGC")
  file.add_feature(
    @src.GckFeature::new(
      name="gene",
      type_="CDS",
      direction=1,
      segments=[(0, 2), (4, 6), (8, 10)],
    ),
  )
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.features.length(), 1)
  assert_eq(parsed.features[0].segments.length(), 3)
  let (s0, e0) = parsed.features[0].segments[0]
  let (s1, e1) = parsed.features[0].segments[1]
  let (s2, e2) = parsed.features[0].segments[2]
  assert_eq(s0, 0)
  assert_eq(e0, 2)
  assert_eq(s1, 4)
  assert_eq(e1, 6)
  assert_eq(s2, 8)
  assert_eq(e2, 10)
}

// ============================================================================
// 11. Edge cases
// ============================================================================

///|
test "gck_feature_direction_values" {
  let f1 = @src.GckFeature::new(
    name="f1",
    type_="CDS",
    direction=0,
    segments=[(0, 10)],
  )
  let f2 = @src.GckFeature::new(
    name="f2",
    type_="CDS",
    direction=1,
    segments=[(0, 10)],
  )
  let f3 = @src.GckFeature::new(
    name="f3",
    type_="CDS",
    direction=2,
    segments=[(0, 10)],
  )
  assert_eq(f1.direction, 0)
  assert_eq(f2.direction, 1)
  assert_eq(f3.direction, 2)
}

///|
test "gck_seq_type_protein_round_trip" {
  let file = @src.GckFile::new()
  file.set_seq_type(@src.GckSeqType::from_int(2)) // Protein
  file.set_sequence("MKTAYIAKQR")
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.seq_type.to_int(), 2)
  assert_eq(parsed.seq_type.to_string(), "Protein")
  assert_eq(parsed.sequence, "MKTAYIAKQR")
}

///|
test "gck_version_round_trip" {
  let file = @src.GckFile::new()
  file.set_version(2)
  file.set_sequence("ATGC")
  let hex = @src.gck_write(file)
  let parsed = @src.gck_parse(hex)
  assert_eq(parsed.version, 2)
}
