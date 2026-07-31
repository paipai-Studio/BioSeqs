// Whitebox tests for new SeqIO functionality:
//   - EMBL format parser
//   - PIR/NBRF format parser
//   - Tab format parser + writer
//   - GenBank format writer
//   - seqio_convert helper
// Run `moon test` to execute.

// --- EMBL Parser Tests ---

///|
test "parse_embl_basic_single" {
  let embl_text =
    "ID   HSBGLOD; SV 1; linear; genomic DNA; STD; HUM; 500 BP.\n" +
    "XX\n" +
    "AC   M12345;\n" +
    "XX\n" +
    "DE   Human beta-globin gene region.\n" +
    "XX\n" +
    "SQ   Sequence 500 BP; 125 A; 125 C; 125 G; 125 T; 0 other;\n" +
    "     atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc\n" +
    "    61 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc\n" +
    "//\n"
  let records = @src.parse_embl(embl_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "M12345")
  assert_eq(records[0].name, "HSBGLOD")
  assert_eq(records[0].description, "Human beta-globin gene region")
  assert_eq(records[0].seq.to_string(), "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC")
}

///|
test "parse_embl_multiple_records" {
  let embl_text =
    "ID   SEQ1; SV 1; linear; DNA; STD; 10 BP.\n" +
    "XX\n" +
    "AC   ACC001;\n" +
    "XX\n" +
    "DE   First sequence.\n" +
    "XX\n" +
    "SQ   Sequence 10 BP; 2 A; 3 C; 2 G; 3 T; 0 other;\n" +
    "     aaccggtttt\n" +
    "//\n" +
    "ID   SEQ2; SV 1; linear; DNA; STD; 8 BP.\n" +
    "XX\n" +
    "AC   ACC002;\n" +
    "XX\n" +
    "DE   Second sequence description.\n" +
    "XX\n" +
    "SQ   Sequence 8 BP; 2 A; 2 C; 2 G; 2 T; 0 other;\n" +
    "     acgtacgt\n" +
    "//\n"
  let records = @src.parse_embl(embl_text)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "ACC001")
  assert_eq(records[0].name, "SEQ1")
  assert_eq(records[0].seq.to_string(), "AACCGGTTTT")
  assert_eq(records[1].id, "ACC002")
  assert_eq(records[1].name, "SEQ2")
  assert_eq(records[1].description, "Second sequence description")
  assert_eq(records[1].seq.to_string(), "ACGTACGT")
}

///|
test "parse_embl_no_ac_uses_id" {
  let embl_text =
    "ID   MYSEQ; SV 1; linear; DNA; STD; 4 BP.\n" +
    "XX\n" +
    "DE   Sequence without accession.\n" +
    "XX\n" +
    "SQ   Sequence 4 BP; 1 A; 1 C; 1 G; 1 T; 0 other;\n" +
    "     acgt\n" +
    "//\n"
  let records = @src.parse_embl(embl_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "MYSEQ")
  assert_eq(records[0].name, "MYSEQ")
  assert_eq(records[0].seq.to_string(), "ACGT")
}

// --- PIR/NBRF Parser Tests ---

///|
test "parse_pir_basic" {
  let pir_text =
    ">P1;ALB_HUMAN\n" +
    "Serum albumin precursor (Human).\n" +
    "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAK*"
  let records = @src.parse_pir(pir_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "ALB_HUMAN")
  assert_eq(records[0].name, "ALB_HUMAN")
  assert_eq(records[0].description, "Serum albumin precursor (Human).")
  assert_true(records[0].seq.to_string().has_prefix("MKWVTFISLL"))
  assert_true(records[0].seq.to_string().has_suffix("FAK"))
}

///|
test "parse_pir_multiple" {
  let pir_text =
    ">P1;PROT1\n" +
    "Protein one.\n" +
    "MVLSQDEVVCF* \n" +
    ">P1;PROT2\n" +
    "Protein two.\n" +
    "MEEPQSDPSVEPPLS*"
  let records = @src.parse_pir(pir_text)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "PROT1")
  assert_eq(records[0].seq.to_string(), "MVLSQDEVVCF")
  assert_eq(records[1].id, "PROT2")
  assert_eq(records[1].seq.to_string(), "MEEPQSDPSVEPPLS")
}

///|
test "parse_pir_sequence_spaces_ignored" {
  let pir_text =
    ">DL;DNA1\n" +
    "Linear DNA fragment.\n" +
    "atgc atgc atgc *"
  let records = @src.parse_pir(pir_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "DNA1")
  assert_eq(records[0].seq.to_string(), "atgcatgcatgc")
}

// --- Tab Format Tests ---

///|
test "parse_tab_basic" {
  let tab_text = "ID001\tACGTACGT\nID002\tTTTTCCCC"
  let records = @src.parse_tab(tab_text)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "ID001")
  assert_eq(records[0].seq.to_string(), "ACGTACGT")
  assert_eq(records[1].id, "ID002")
  assert_eq(records[1].seq.to_string(), "TTTTCCCC")
}

///|
test "parse_tab_blank_lines_ignored" {
  let tab_text = "\n\nID001\tACGT\n\nID002\tTTTT\n\n"
  let records = @src.parse_tab(tab_text)
  assert_eq(records.length(), 2)
}

///|
test "parse_tab_sequence_spaces_ignored" {
  let tab_text = "ID001\tACGT ACGT ACGT"
  let records = @src.parse_tab(tab_text)
  assert_eq(records[0].seq.to_string(), "ACGTACGTACGT")
}

///|
test "parse_tab_two_tabs_error" {
  let tab_text = "ID001\ttab1\ttab2"
  try {
    let _ = @src.parse_tab(tab_text)
    fail("Expected SeqIOError")
  } catch {
    _ => assert_true(true)
  }
}

///|
test "parse_tab_no_tab_error" {
  let tab_text = "ID001 ACGT"
  try {
    let _ = @src.parse_tab(tab_text)
    fail("Expected SeqIOError")
  } catch {
    _ => assert_true(true)
  }
}

///|
test "write_tab_basic" {
  let records = [
    @src.SeqRecord::new(@src.Seq::new("ACGT"), id="r1", description="r1"),
    @src.SeqRecord::new(@src.Seq::new("TTTT"), id="r2", description="r2"),
  ]
  let out = @src.write_tab(records)
  assert_eq(out, "r1\tACGT\nr2\tTTTT\n")
}

///|
test "write_tab_roundtrip" {
  let tab_text = "s1\tACGTACGT\ns2\tGGGGCCCC\ns3\tAAAA\n"
  let records = @src.parse_tab(tab_text)
  let out = @src.write_tab(records)
  assert_eq(out, tab_text)
}

// --- GenBank Writer Tests ---

///|
test "write_genbank_basic_structure" {
  let rec = @src.SeqRecord::new(
    @src.Seq::new("atgcatgcatgcatgcatgc"),
    id="NM_001",
    name="TestGene",
    description="Test gene for unit testing",
  )
  let out = @src.write_genbank([rec])
  // Check that required sections are present
  assert_true(out.contains("LOCUS"))
  assert_true(out.contains("TestGene"))
  assert_true(out.contains("DEFINITION"))
  assert_true(out.contains("Test gene for unit testing."))
  assert_true(out.contains("ACCESSION"))
  assert_true(out.contains("NM_001"))
  assert_true(out.contains("VERSION"))
  assert_true(out.contains("ORIGIN"))
  assert_true(out.contains("//"))
}

///|
test "write_genbank_sequence_format" {
  let seq = "atgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgc"
  let rec = @src.SeqRecord::new(@src.Seq::new(seq), id="X1", name="LOC1")
  let out = @src.write_genbank([rec])
  // ORIGIN lines should be present with the sequence
  assert_true(out.contains("atgcatgcat"))
  // Should contain position label
  assert_true(out.contains("      1 "))
}

///|
test "write_genbank_multiple_records" {
  let r1 = @src.SeqRecord::new(@src.Seq::new("ACGT"), id="A1", name="L1", description="First")
  let r2 = @src.SeqRecord::new(@src.Seq::new("TGCA"), id="B2", name="L2", description="Second")
  let out = @src.write_genbank([r1, r2])
  // Should have two LOCUS and two terminators
  assert_eq(out.split("LOCUS").length(), 3) // includes 1 before first match
  assert_eq(out.split("//").length(), 3)
}

// --- seqio_parse/write unified interface Tests ---

///|
test "seqio_parse_embl_through_unified" {
  let embl =
    "ID   S1; SV 1; linear; DNA; STD; 4 BP.\n" +
    "XX\n" +
    "AC   A001;\n" +
    "XX\n" +
    "DE   A test.\n" +
    "XX\n" +
    "SQ   Sequence 4 BP; 1 A; 1 C; 1 G; 1 T; 0 other;\n" +
    "     acgt\n" +
    "//\n"
  try {
    let records = @src.seqio_parse(embl, "embl")
    assert_eq(records.length(), 1)
    assert_eq(records[0].id, "A001")
  } catch {
    e => fail("seqio_parse embl raised error: " + e.to_string())
  }
}

///|
test "seqio_parse_pir_through_unified" {
  let pir = ">P1;P1\nA protein.\nMACE*"
  try {
    let records = @src.seqio_parse(pir, "pir")
    assert_eq(records.length(), 1)
    assert_eq(records[0].id, "P1")
    assert_eq(records[0].seq.to_string(), "MACE")
    // Also try "nbrf" alias
    let records2 = @src.seqio_parse(pir, "nbrf")
    assert_eq(records2.length(), 1)
  } catch {
    e => fail("seqio_parse pir raised error: " + e.to_string())
  }
}

///|
test "seqio_parse_tab_through_unified" {
  let tab = "s1\tACGT\ns2\tTGCA"
  try {
    let records = @src.seqio_parse(tab, "tab")
    assert_eq(records.length(), 2)
  } catch {
    e => fail("seqio_parse tab raised error: " + e.to_string())
  }
}

///|
test "seqio_write_genbank_tab_through_unified" {
  let records = [
    @src.SeqRecord::new(@src.Seq::new("ACGT"), id="r1"),
  ]
  try {
    let gb_text = @src.seqio_write(records, "genbank")
    assert_true(gb_text.contains("LOCUS"))
    let gb_text2 = @src.seqio_write(records, "gb")
    assert_eq(gb_text, gb_text2)
    let tab_text = @src.seqio_write(records, "tab")
    assert_true(tab_text.contains("\t"))
  } catch {
    e => fail("seqio_write raised error: " + e.to_string())
  }
}

// --- seqio_convert Tests ---

///|
test "seqio_convert_fasta_to_tab" {
  let fasta = ">r1 First record\nACGTACGT\n>r2 Second\nTTTTCCCC"
  try {
    let tab = @src.seqio_convert(fasta, "fasta", "tab")
    assert_eq(tab, "r1\tACGTACGT\nr2\tTTTTCCCC\n")
  } catch {
    e => fail("seqio_convert fasta->tab raised error: " + e.to_string())
  }
}

///|
test "seqio_convert_tab_to_fasta" {
  let tab = "x1\tAAAA\nx2\tCCCCGGGG"
  try {
    let fasta = @src.seqio_convert(tab, "tab", "fasta")
    let parsed_back = @src.parse_fasta(fasta)
    assert_eq(parsed_back.length(), 2)
    assert_eq(parsed_back[0].id, "x1")
    assert_eq(parsed_back[0].seq.to_string(), "AAAA")
    assert_eq(parsed_back[1].id, "x2")
    assert_eq(parsed_back[1].seq.to_string(), "CCCCGGGG")
  } catch {
    e => fail("seqio_convert tab->fasta raised error: " + e.to_string())
  }
}

///|
test "seqio_convert_genbank_write_to_genbank_roundtrip" {
  // Build a GenBank file from scratch, parse back, check consistency
  let orig_seq = "ATGCATGCATGC"
  let rec = @src.SeqRecord::new(
    @src.Seq::new(orig_seq),
    id="ACC999",
    name="LOC999",
    description="Round trip test gene",
  )
  try {
    let gb_text = @src.seqio_write([rec], "genbank")
    let parsed = @src.seqio_parse(gb_text, "genbank")
    assert_eq(parsed.length(), 1)
    assert_eq(parsed[0].name, "LOC999")
    // GenBank parser uppercases the sequence; the ORIGIN section writer uses lowercase,
    // but parse_genbank normalizes to uppercase.
    assert_eq(parsed[0].seq.to_string(), "ATGCATGCATGC")
  } catch {
    e => fail("gb roundtrip error: " + e.to_string())
  }
}

///|
test "seqio_convert_unknown_format_errors" {
  try {
    let _ = @src.seqio_convert("s1\tACGT", "tab", "no_such_format_xyz")
    fail("Expected SeqIOError for unknown output format")
  } catch {
    _ => assert_true(true)
  }
  try {
    let _ = @src.seqio_convert("not valid", "no_such_format_in", "tab")
    fail("Expected SeqIOError for unknown input format")
  } catch {
    _ => assert_true(true)
  }
}

// --- seqio_read + seqio_to_dict Tests ---

///|
test "seqio_read_single_record" {
  let tab = "only\tACGT"
  try {
    let rec = @src.seqio_read(tab, "tab")
    assert_eq(rec.id, "only")
    assert_eq(rec.seq.to_string(), "ACGT")
  } catch {
    e => fail("seqio_read error: " + e.to_string())
  }
}

///|
test "seqio_read_too_many_raises" {
  let tab = "s1\tACGT\ns2\tTGCA"
  try {
    let _ = @src.seqio_read(tab, "tab")
    fail("Expected SeqIOError for multiple records")
  } catch {
    _ => assert_true(true)
  }
}

///|
test "seqio_to_dict_basic" {
  let records = [
    @src.SeqRecord::new(@src.Seq::new("A"), id="g1"),
    @src.SeqRecord::new(@src.Seq::new("C"), id="g2"),
  ]
  let d = @src.seqio_to_dict(records)
  assert_eq(d.length(), 2)
  match d.get("g1") {
    Some(r) => assert_eq(r.seq.to_string(), "A")
    None => fail("g1 missing")
  }
  match d.get("g2") {
    Some(r) => assert_eq(r.seq.to_string(), "C")
    None => fail("g2 missing")
  }
}
