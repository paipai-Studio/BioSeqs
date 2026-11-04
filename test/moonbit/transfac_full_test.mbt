///|
/// TRANSFAC full motif format parsing module tests.

///|
test "transfac_reference_new" {
  let r = @src.TransfacReference::new()
  assert_eq(r.number, "")
  assert_eq(r.authors, "")
  assert_eq(r.title, "")
  assert_eq(r.journal, "")
  assert_eq(r.pubmed_id, "")
}

///|
test "transfac_motif_full_new" {
  let m = @src.TransfacMotifFull::new(accession="T00001")
  assert_eq(m.accession, "T00001")
  assert_eq(m.identifier, "")
  assert_eq(m.description, "")
  assert_eq(m.binding_factors.length(), 0)
  assert_eq(m.comments, "")
  assert_eq(m.name, "")
  assert_eq(m.species, "")
  assert_eq(m.pfm.length(), 0)
  assert_eq(m.alphabet.length(), 4) // default A,C,G,T
  assert_eq(m.width, 0)
  assert_eq(m.references.length(), 0)
}

///|
test "transfac_motif_full_alphabet_index" {
  let m = @src.TransfacMotifFull::new(accession="T00001")
  assert_eq(m.alphabet_index("A"), 0)
  assert_eq(m.alphabet_index("C"), 1)
  assert_eq(m.alphabet_index("G"), 2)
  assert_eq(m.alphabet_index("T"), 3)
  assert_eq(m.alphabet_index("N"), -1)
}

///|
test "transfac_record_full_new" {
  let record = @src.TransfacRecordFull::new()
  assert_eq(record.version, "")
  assert_eq(record.motifs.length(), 0)
  assert_eq(record.num_motifs(), 0)
}

///|
test "transfac_record_full_find" {
  let record = @src.TransfacRecordFull::new()
  let m = @src.TransfacMotifFull::new(accession="T00001")
  record.motifs.push(m)
  match record.find("T00001") {
    Some(found) => assert_eq(found.accession, "T00001")
    None => assert_true(false)
  }
  match record.find("nonexistent") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "transfac_parse_full_version" {
  let text = "VV  January 2024\nXX\n//\n"
  let record = @src.transfac_parse_full(text)
  assert_eq(record.version, "January 2024")
}

///|
test "transfac_parse_full_empty" {
  let record = @src.transfac_parse_full("")
  assert_eq(record.num_motifs(), 0)
}

///|
test "transfac_parse_full_single_motif" {
  let text = "AC  T00001\nXX\nID  AP-1\nXX\nDE  activator protein 1\nXX\nNA  AP-1\nXX\nOS  Homo sapiens\nXX\nBF  Jun\nXX\nBF  Fos\nXX\nP0  A C G T\n01  4 0 0 0\n02  0 0 4 0\n03  0 4 0 0\n04  0 0 0 4\nXX\n//\n"
  let record = @src.transfac_parse_full(text)
  assert_eq(record.num_motifs(), 1)
  let m = record.motifs[0]
  assert_eq(m.accession, "T00001")
  assert_eq(m.identifier, "AP-1")
  assert_eq(m.description, "activator protein 1")
  assert_eq(m.name, "AP-1")
  assert_eq(m.species, "Homo sapiens")
  assert_eq(m.binding_factors.length(), 2)
  assert_eq(m.binding_factors[0], "Jun")
  assert_eq(m.binding_factors[1], "Fos")
  assert_eq(m.width, 4)
  assert_eq(m.pfm.length(), 4) // 4 rows for A, C, G, T
}

///|
test "transfac_parse_full_sample_text" {
  let record = @src.transfac_parse_full(@src.transfac_sample_text_full())
  assert_eq(record.num_motifs(), 2)
  assert_eq(record.motifs[0].accession, "T00001")
  assert_eq(record.motifs[0].identifier, "AP-1")
  assert_eq(record.motifs[1].accession, "T00002")
  assert_eq(record.motifs[1].identifier, "NF-kappaB")
}

///|
test "transfac_parse_full_pfm" {
  let record = @src.transfac_parse_full(@src.transfac_sample_text_full())
  let m = record.motifs[0]
  // PFM should have 4 rows (A, C, G, T) and 5 columns
  assert_eq(m.pfm.length(), 4)
  assert_eq(m.width, 5)
  // First column: A=4, C=0, G=0, T=0
  assert_eq(m.pfm[0][0], 4)
  assert_eq(m.pfm[1][0], 0)
  assert_eq(m.pfm[2][0], 0)
  assert_eq(m.pfm[3][0], 0)
}

///|
test "transfac_motif_full_consensus" {
  let m = @src.transfac_sample_motif_full()
  let cons = m.consensus()
  // PFM: A=[4,0,0,0,3], C=[0,0,4,0,1], G=[0,4,0,0,0], T=[0,0,0,4,0]
  // Consensus should be AGCTA
  assert_eq(cons, "AGCTA")
}

///|
test "transfac_motif_full_total_at" {
  let m = @src.transfac_sample_motif_full()
  // Position 0: A=4, C=0, G=0, T=0 => total=4
  assert_eq(m.total_at(0), 4)
  // Position 4: A=3, C=1, G=0, T=0 => total=4
  assert_eq(m.total_at(4), 4)
  // Out of range
  assert_eq(m.total_at(10), 0)
}

///|
test "transfac_motif_full_frequency" {
  let m = @src.transfac_sample_motif_full()
  // Position 0: A=4/4=1.0
  assert_eq(m.frequency("A", 0), 1.0)
  // Position 4: A=3/4=0.75
  assert_eq(m.frequency("A", 4), 0.75)
  assert_eq(m.frequency("C", 4), 0.25)
}

///|
test "transfac_motif_full_to_string" {
  let m = @src.transfac_sample_motif_full()
  let s = @src.transfac_motif_full_to_string(m)
  assert_true(s.contains("AC  T00001"))
  assert_true(s.contains("ID  AP-1"))
  assert_true(s.contains("P0"))
  assert_true(s.contains("//"))
}

///|
test "transfac_record_full_to_string" {
  let record = @src.transfac_parse_full(@src.transfac_sample_text_full())
  let s = @src.transfac_record_full_to_string(record)
  assert_true(s.contains("VV"))
  assert_true(s.contains("AC  T00001"))
  assert_true(s.contains("ID  AP-1"))
}

///|
test "transfac_parse_full_multiple_motifs" {
  let record = @src.transfac_parse_full(@src.transfac_sample_text_full())
  // First motif
  assert_eq(record.motifs[0].binding_factors.length(), 2)
  assert_eq(record.motifs[0].species, "Homo sapiens")
  // Second motif
  assert_eq(record.motifs[1].species, "Mus musculus")
  assert_eq(record.motifs[1].width, 4)
}

///|
test "transfac_parse_full_comments" {
  let record = @src.transfac_parse_full(@src.transfac_sample_text_full())
  let m = record.motifs[0]
  assert_true(m.comments.length() > 0)
  assert_true(m.comments.contains("AP-1"))
}

///|
test "transfac_parse_full_second_motif_pfm" {
  let record = @src.transfac_parse_full(@src.transfac_sample_text_full())
  let m = record.motifs[1]
  // PFM for NF-kappaB: 4 positions
  assert_eq(m.pfm.length(), 4)
  assert_eq(m.width, 4)
  // Position 0: A=0, C=0, G=8, T=0
  assert_eq(m.pfm[0][0], 0)
  assert_eq(m.pfm[1][0], 0)
  assert_eq(m.pfm[2][0], 8)
  assert_eq(m.pfm[3][0], 0)
}
