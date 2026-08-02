///|
/// Tests for the SnapGene I/O module (Bio.SeqIO.SnapGeneIO equivalent).
///
/// These tests exercise the public API exposed in `src/snapgene_io.mbt`:
/// the SnapgeneSeqType enum, SnapgeneFeature / SnapgenePrimer / SnapgeneFile
/// structs, the hex/byte serializers (`snapgene_parse` / `snapgene_write` /
/// `snapgene_from_bytes` / `snapgene_to_bytes`), the SeqRecord integration
/// (`snapgene_to_seqrecord` / `snapgene_from_seqrecord`), and utility
/// functions (`snapgene_count_features_by_type`, `snapgene_seq_length`,
/// `snapgene_summary`, `snapgene_sample`).

// ============================================================================
// 1. SnapgeneSeqType - to_string / from_string
// ============================================================================

///|
test "snapgene_seq_type_to_string" {
  assert_eq(@src.SnapgeneSeqType::from_string("DNA").to_string(), "DNA")
  assert_eq(@src.SnapgeneSeqType::from_string("RNA").to_string(), "RNA")
  assert_eq(@src.SnapgeneSeqType::from_string("protein").to_string(), "protein")
  assert_eq(@src.SnapgeneSeqType::from_string("unknown").to_string(), "unknown")
}

///|
test "snapgene_seq_type_from_string" {
  assert_true(@src.SnapgeneSeqType::from_string("DNA") is @src.SnapgeneSeqType::SnapgeneDna)
  assert_true(@src.SnapgeneSeqType::from_string("RNA") is @src.SnapgeneSeqType::SnapgeneRna)
  assert_true(
    @src.SnapgeneSeqType::from_string("protein") is @src.SnapgeneSeqType::SnapgeneProtein,
  )
  assert_true(
    @src.SnapgeneSeqType::from_string("xxx") is @src.SnapgeneSeqType::SnapgeneUnknown,
  )
}

///|
test "snapgene_seq_type_from_string_round_trip" {
  let dna = @src.SnapgeneSeqType::from_string("DNA")
  let rna = @src.SnapgeneSeqType::from_string("RNA")
  let protein = @src.SnapgeneSeqType::from_string("protein")
  assert_eq(@src.SnapgeneSeqType::from_string(dna.to_string()).to_string(), "DNA")
  assert_eq(@src.SnapgeneSeqType::from_string(rna.to_string()).to_string(), "RNA")
  assert_eq(
    @src.SnapgeneSeqType::from_string(protein.to_string()).to_string(),
    "protein",
  )
}

// ============================================================================
// 2. SnapgeneFile - construction and defaults
// ============================================================================

///|
test "snapgene_file_new_defaults" {
  let file = @src.SnapgeneFile::new()
  assert_eq(@src.snapgene_seq_length(file), 0)
  assert_eq(file.features.length(), 0)
  assert_eq(file.primers.length(), 0)
  assert_eq(file.sequence, "")
  assert_eq(file.comments, "")
  assert_eq(file.notes, "")
  assert_eq(file.export_version, "")
  assert_false(file.is_circular)
}

// ============================================================================
// 3. SnapgeneFeature - construction
// ============================================================================

///|
test "snapgene_feature_new" {
  let qualifiers : Map[String, String] = Map([], capacity=0)
  qualifiers["gene"] = "lacZ"
  let f = @src.SnapgeneFeature::new(
    name="orfA",
    type_="CDS",
    direction=1,
    segments=[(0, 99)],
    qualifiers=qualifiers,
  )
  assert_eq(f.name, "orfA")
  assert_eq(f.type_, "CDS")
  assert_eq(f.direction, 1)
  assert_eq(f.segments.length(), 1)
  // 0-based inclusive segments
  let (s, e) = f.segments[0]
  assert_eq(s, 0)
  assert_eq(e, 99)
}

///|
test "snapgene_feature_new_empty_qualifiers" {
  let empty_q : Map[String, String] = Map([], capacity=0)
  let f = @src.SnapgeneFeature::new(
    name="promoter1",
    type_="promoter",
    direction=-1,
    segments=[(10, 30), (40, 60)],
    qualifiers=empty_q,
  )
  assert_eq(f.name, "promoter1")
  assert_eq(f.type_, "promoter")
  assert_eq(f.direction, -1)
  assert_eq(f.segments.length(), 2)
}

// ============================================================================
// 4. SnapgenePrimer - construction
// ============================================================================

///|
test "snapgene_primer_new" {
  let p = @src.SnapgenePrimer::new(
    name="Fwd",
    sequence="ATGGCTAGCAA",
    description="forward primer",
    bind_to=1,
    binds_at_position=0,
  )
  assert_eq(p.name, "Fwd")
  assert_eq(p.sequence, "ATGGCTAGCAA")
  assert_eq(p.description, "forward primer")
  assert_eq(p.bind_to, 1)
  assert_eq(p.binds_at_position, 0)
}

// ============================================================================
// 5. snapgene_sample - sample data
// ============================================================================

///|
test "snapgene_sample_structure" {
  let file = @src.snapgene_sample()
  assert_true(@src.snapgene_seq_length(file) > 0)
  assert_eq(file.features.length(), 2)
  assert_eq(file.primers.length(), 2)
  assert_true(file.is_circular)
  assert_eq(file.comments, "Sample plasmid")
  assert_eq(file.export_version, "1.0")
}

///|
test "snapgene_sample_feature_types" {
  let file = @src.snapgene_sample()
  assert_eq(@src.snapgene_count_features_by_type(file, "CDS"), 1)
  assert_eq(@src.snapgene_count_features_by_type(file, "promoter"), 1)
  assert_eq(@src.snapgene_count_features_by_type(file, "terminator"), 0)
}

// ============================================================================
// 6. snapgene_count_features_by_type
// ============================================================================

///|
test "snapgene_count_features_by_type_empty" {
  let file = @src.SnapgeneFile::new()
  assert_eq(@src.snapgene_count_features_by_type(file, "CDS"), 0)
}

///|
test "snapgene_count_features_by_type_multiple" {
  let file = @src.SnapgeneFile::new()
  let q : Map[String, String] = Map([], capacity=0)
  file.add_feature(
    @src.SnapgeneFeature::new(
      name="a",
      type_="CDS",
      direction=1,
      segments=[(0, 10)],
      qualifiers=q,
    ),
  )
  file.add_feature(
    @src.SnapgeneFeature::new(
      name="b",
      type_="CDS",
      direction=1,
      segments=[(20, 30)],
      qualifiers=q,
    ),
  )
  file.add_feature(
    @src.SnapgeneFeature::new(
      name="c",
      type_="promoter",
      direction=1,
      segments=[(40, 50)],
      qualifiers=q,
    ),
  )
  assert_eq(@src.snapgene_count_features_by_type(file, "CDS"), 2)
  assert_eq(@src.snapgene_count_features_by_type(file, "promoter"), 1)
}

// ============================================================================
// 7. snapgene_seq_length
// ============================================================================

///|
test "snapgene_seq_length_empty" {
  let file = @src.SnapgeneFile::new()
  assert_eq(@src.snapgene_seq_length(file), 0)
}

///|
test "snapgene_seq_length_nonempty" {
  let file = @src.SnapgeneFile::new()
  file.set_sequence("ATGGCTAGCAA")
  assert_eq(@src.snapgene_seq_length(file), 11)
}

// ============================================================================
// 8. snapgene_to_bytes / snapgene_from_bytes round trip
// ============================================================================

///|
test "snapgene_to_bytes_from_bytes_round_trip" {
  let file = @src.snapgene_sample()
  let bytes = @src.snapgene_to_bytes(file)
  let parsed = @src.snapgene_from_bytes(bytes)
  assert_eq(parsed.sequence, file.sequence)
  assert_eq(parsed.is_circular, file.is_circular)
  assert_eq(parsed.features.length(), file.features.length())
  assert_eq(parsed.primers.length(), file.primers.length())
}

///|
test "snapgene_write_parse_round_trip" {
  let file = @src.snapgene_sample()
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  assert_eq(parsed.sequence, file.sequence)
  assert_eq(parsed.is_circular, file.is_circular)
  assert_eq(parsed.features.length(), file.features.length())
  assert_eq(parsed.primers.length(), file.primers.length())
}

///|
test "snapgene_to_bytes_empty_file" {
  let file = @src.SnapgeneFile::new()
  let bytes = @src.snapgene_to_bytes(file)
  // Should at least contain a comments packet header (5 bytes)
  assert_true(bytes.length() >= 5)
  let parsed = @src.snapgene_from_bytes(bytes)
  assert_eq(parsed.sequence, "")
}

///|
test "snapgene_parse_preserves_comments" {
  let file = @src.snapgene_sample()
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  assert_eq(parsed.comments, file.comments)
  assert_eq(parsed.export_version, file.export_version)
}

///|
test "snapgene_parse_preserves_features" {
  let file = @src.snapgene_sample()
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  assert_eq(parsed.features.length(), 2)
  // Features should preserve name and type
  assert_eq(parsed.features[0].name, "ORF1")
  assert_eq(parsed.features[0].type_, "CDS")
  assert_eq(parsed.features[1].name, "Promoter1")
  assert_eq(parsed.features[1].type_, "promoter")
}

///|
test "snapgene_parse_preserves_primers" {
  let file = @src.snapgene_sample()
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  assert_eq(parsed.primers.length(), 2)
  assert_eq(parsed.primers[0].name, "Fwd_Primer")
  assert_eq(parsed.primers[0].sequence, "ATGGCTAGCAAT")
  assert_eq(parsed.primers[1].name, "Rev_Primer")
}

// ============================================================================
// 9. snapgene_to_seqrecord / snapgene_from_seqrecord
// ============================================================================

///|
test "snapgene_to_seqrecord" {
  let file = @src.snapgene_sample()
  let rec = @src.snapgene_to_seqrecord(file)
  assert_eq(rec.seq.data, file.sequence)
}

///|
test "snapgene_from_seqrecord" {
  let rec = @src.SeqRecord::new(@src.Seq::new("ATGGCTAGCAA"), id="test")
  let file = @src.snapgene_from_seqrecord(rec)
  assert_eq(file.sequence, "ATGGCTAGCAA")
}

///|
test "snapgene_seqrecord_round_trip" {
  let file = @src.snapgene_sample()
  let rec = @src.snapgene_to_seqrecord(file)
  let file2 = @src.snapgene_from_seqrecord(rec)
  assert_eq(file2.sequence, file.sequence)
}

// ============================================================================
// 10. snapgene_summary
// ============================================================================

///|
test "snapgene_summary_contains_key_info" {
  let file = @src.snapgene_sample()
  let summary = @src.snapgene_summary(file)
  assert_true(summary.contains("SnapGene File Summary:"))
  assert_true(summary.contains("Features: 2"))
  assert_true(summary.contains("Primers: 2"))
}

///|
test "snapgene_summary_empty_file" {
  let file = @src.SnapgeneFile::new()
  let summary = @src.snapgene_summary(file)
  assert_true(summary.contains("SnapGene File Summary:"))
  assert_true(summary.contains("Features: 0"))
  assert_true(summary.contains("Primers: 0"))
}

// ============================================================================
// 11. Custom file construction and serialization
// ============================================================================

///|
test "snapgene_custom_file_round_trip" {
  let file = @src.SnapgeneFile::new()
  file.set_sequence("ATGAAATAG")
  file.set_circular(false)
  file.set_comments("custom")
  let q : Map[String, String] = Map([], capacity=0)
  file.add_feature(
    @src.SnapgeneFeature::new(
      name="gene1",
      type_="CDS",
      direction=1,
      segments=[(0, 8)],
      qualifiers=q,
    ),
  )
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  assert_eq(parsed.sequence, "ATGAAATAG")
  assert_false(parsed.is_circular)
  assert_eq(parsed.features.length(), 1)
  assert_eq(parsed.features[0].name, "gene1")
}

///|
test "snapgene_rna_type_round_trip" {
  let file = @src.SnapgeneFile::new()
  file.set_seq_type(@src.SnapgeneSeqType::from_string("RNA"))
  file.set_sequence("AUGAAAUAG")
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  assert_eq(parsed.seq_type.to_string(), "RNA")
  assert_eq(parsed.sequence, "AUGAAAUAG")
}

///|
test "snapgene_notes_round_trip" {
  let file = @src.SnapgeneFile::new()
  file.set_sequence("ATGC")
  file.set_notes("These are notes")
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  assert_eq(parsed.notes, "These are notes")
}

///|
test "snapgene_export_version_default" {
  // When export_version is empty, the writer should still produce valid output
  let file = @src.SnapgeneFile::new()
  file.set_sequence("ATGC")
  file.set_export_version("")
  let hex = @src.snapgene_write(file)
  let parsed = @src.snapgene_parse(hex)
  // Default version "1.0" should be applied
  assert_eq(parsed.export_version, "1.0")
}

// ============================================================================
// 12. Edge cases
// ============================================================================

///|
test "snapgene_parse_empty_hex" {
  let parsed = @src.snapgene_parse("")
  assert_eq(parsed.sequence, "")
  assert_eq(parsed.features.length(), 0)
}

///|
test "snapgene_to_bytes_starts_with_comments_packet" {
  // First packet must be Comments (type 0)
  let file = @src.SnapgeneFile::new()
  file.set_sequence("ATGC")
  let bytes = @src.snapgene_to_bytes(file)
  assert_eq(bytes[0], 0) // Comments packet type
}

///|
test "snapgene_feature_direction_values" {
  let q : Map[String, String] = Map([], capacity=0)
  let f1 = @src.SnapgeneFeature::new(
    name="f1",
    type_="CDS",
    direction=1,
    segments=[(0, 10)],
    qualifiers=q,
  )
  let f2 = @src.SnapgeneFeature::new(
    name="f2",
    type_="CDS",
    direction=-1,
    segments=[(0, 10)],
    qualifiers=q,
  )
  let f3 = @src.SnapgeneFeature::new(
    name="f3",
    type_="CDS",
    direction=0,
    segments=[(0, 10)],
    qualifiers=q,
  )
  assert_eq(f1.direction, 1)
  assert_eq(f2.direction, -1)
  assert_eq(f3.direction, 0)
}
