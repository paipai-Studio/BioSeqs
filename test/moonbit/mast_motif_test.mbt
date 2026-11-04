///|
/// MAST motif search output parsing module tests.

///|
test "mast_motif_match_new" {
  let m = @src.MastMotifMatch::new("motif1", 50, 0.001)
  assert_eq(m.motif_id, "motif1")
  assert_eq(m.position, 50)
  assert_eq(m.pvalue, 0.001)
  assert_eq(m.length, 0)
  assert_eq(m.score, 0.0)
  assert_eq(m.strand, "+")
  assert_eq(m.sequence, "")
}

///|
test "mast_motif_seq_new" {
  let seq = @src.MastMotifSeq::new("seq1")
  assert_eq(seq.name, "seq1")
  assert_eq(seq.description, "")
  assert_eq(seq.length, 0)
  assert_eq(seq.combined_pvalue, 1.0)
  assert_eq(seq.combined_evalue, 0.0)
  assert_eq(seq.matches.length(), 0)
  assert_eq(seq.diagram, "")
}

///|
test "mast_motif_seq_add_match" {
  let seq = @src.MastMotifSeq::new("seq1")
  let m = @src.MastMotifMatch::new("motif1", 50, 0.001)
  seq.add_match(m)
  assert_eq(seq.num_matches(), 1)
}

///|
test "mast_motif_seq_best_match" {
  let seq = @src.MastMotifSeq::new("seq1")
  seq.add_match(@src.MastMotifMatch::new("m1", 10, 0.01))
  seq.add_match(@src.MastMotifMatch::new("m2", 20, 0.001))
  seq.add_match(@src.MastMotifMatch::new("m3", 30, 0.1))
  match seq.best_match() {
    Some(best) => assert_eq(best.motif_id, "m2")
    None => assert_true(false)
  }
}

///|
test "mast_motif_seq_best_match_empty" {
  let seq = @src.MastMotifSeq::new("seq1")
  match seq.best_match() {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "mast_motif_info_new" {
  let info = @src.MastMotifInfo::new("1", 10)
  assert_eq(info.id, "1")
  assert_eq(info.width, 10)
  assert_eq(info.alphabet, "ACGT")
}

///|
test "mast_motif_record_new" {
  let record = @src.MastMotifRecord::new()
  assert_eq(record.version, "")
  assert_eq(record.database, "")
  assert_eq(record.motifs.length(), 0)
  assert_eq(record.sequences.length(), 0)
  assert_eq(record.motif_count, 0)
  assert_eq(record.hit_count, 0)
}

///|
test "mast_motif_record_num_hits" {
  let record = @src.MastMotifRecord::new()
  assert_eq(record.num_hits(), 0)
  let seq = @src.MastMotifSeq::new("seq1")
  record.sequences.push(seq)
  assert_eq(record.num_hits(), 1)
}

///|
test "mast_motif_record_find_sequence" {
  let record = @src.MastMotifRecord::new()
  record.sequences.push(@src.MastMotifSeq::new("seq1"))
  record.sequences.push(@src.MastMotifSeq::new("seq2"))
  match record.find_sequence("seq1") {
    Some(seq) => assert_eq(seq.name, "seq1")
    None => assert_true(false)
  }
  match record.find_sequence("nonexistent") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "mast_motif_record_motif_matches" {
  let record = @src.MastMotifRecord::new()
  let seq1 = @src.MastMotifSeq::new("seq1")
  seq1.add_match(@src.MastMotifMatch::new("m1", 10, 0.01))
  seq1.add_match(@src.MastMotifMatch::new("m2", 20, 0.02))
  record.sequences.push(seq1)
  let seq2 = @src.MastMotifSeq::new("seq2")
  seq2.add_match(@src.MastMotifMatch::new("m1", 30, 0.03))
  record.sequences.push(seq2)
  let m1_matches = record.motif_matches("m1")
  assert_eq(m1_matches.length(), 2)
  let m2_matches = record.motif_matches("m2")
  assert_eq(m2_matches.length(), 1)
}

///|
test "mast_motif_record_best_sequence" {
  let record = @src.mast_motif_sample_record()
  // seq1 has evalue 3.4e-03, seq2 has evalue 1.1e-01
  // seq1 has lower (better) evalue
  match record.best_sequence() {
    Some(best) => assert_eq(best.name, "seq1")
    None => assert_true(false)
  }
}

///|
test "mast_motif_record_best_sequence_empty" {
  let record = @src.MastMotifRecord::new()
  match record.best_sequence() {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "mast_motif_parse_version" {
  let text = "MAST version 5.0.3\n"
  let record = @src.mast_motif_parse(text)
  assert_eq(record.version, "5.0.3")
}

///|
test "mast_motif_parse_empty" {
  let record = @src.mast_motif_parse("")
  assert_eq(record.hit_count, 0)
}

///|
test "mast_motif_parse_sample" {
  let record = @src.mast_motif_parse(@src.mast_motif_sample_text())
  // Should parse version and database
  assert_true(record.version.length() > 0)
  assert_true(record.database.length() > 0)
}

///|
test "mast_motif_sample_record" {
  let record = @src.mast_motif_sample_record()
  assert_eq(record.version, "5.0.3")
  assert_eq(record.database, "sample_db.fasta")
  assert_eq(record.motifs.length(), 3)
  assert_eq(record.motifs[0].id, "1")
  assert_eq(record.motifs[0].width, 10)
  assert_eq(record.sequences.length(), 2)
  assert_eq(record.sequences[0].name, "seq1")
  assert_eq(record.sequences[0].length, 500)
  assert_eq(record.sequences[0].matches.length(), 2)
}

///|
test "mast_motif_match_to_string" {
  let m = @src.MastMotifMatch::new("motif1", 50, 0.001)
  let s = @src.mast_motif_match_to_string(m)
  assert_true(s.contains("motif1"))
  assert_true(s.contains("50"))
  assert_true(s.contains("0.001"))
}

///|
test "mast_motif_seq_to_string" {
  let record = @src.mast_motif_sample_record()
  let seq = record.sequences[0]
  let s = @src.mast_motif_seq_to_string(seq)
  assert_true(s.contains("seq1"))
  assert_true(s.contains("test sequence 1"))
  assert_true(s.contains("500"))
}

///|
test "mast_motif_record_to_string" {
  let record = @src.mast_motif_sample_record()
  let s = @src.mast_motif_record_to_string(record)
  assert_true(s.contains("MAST Record"))
  assert_true(s.contains("5.0.3"))
  assert_true(s.contains("sample_db.fasta"))
}

///|
test "mast_motif_parse_with_database" {
  let text = "MAST version 5.0.3\nDATABASE= test_db.fasta\n"
  let record = @src.mast_motif_parse(text)
  assert_eq(record.version, "5.0.3")
  assert_eq(record.database, "test_db.fasta")
}
