///|
/// Tests for Bio.SearchIO.InterproscanIO module.

// ---------------------------------------------------------------------------
// Record construction and accessors
// ---------------------------------------------------------------------------

test "ips_record_construction" {
  let go : Array[String] = Array::new()
  go.push("GO:0003674")
  go.push("GO:0008150")
  let r = @src.InterproScanRecord::new(
    protein_id="sp|P12345|PROT_HUMAN",
    md5="abc123def456",
    seq_length=500,
    analysis="Pfam",
    signature_acc="PF00001",
    signature_desc="p53 DNA-binding domain",
    start_pos=1,
    stop_pos=100,
    score=150.5,
    status="T",
    date="01-Jan-2024",
    ipr_acc="IPR000001",
    ipr_desc="p53, DNA-binding domain",
    go_terms=go,
  )
  assert_eq(r.protein_id(), "sp|P12345|PROT_HUMAN")
  assert_eq(r.md5(), "abc123def456")
  assert_eq(r.seq_length(), 500)
  assert_eq(r.analysis(), "Pfam")
  assert_eq(r.signature_acc(), "PF00001")
  assert_eq(r.signature_desc(), "p53 DNA-binding domain")
  assert_eq(r.start_pos(), 1)
  assert_eq(r.stop_pos(), 100)
  assert_eq(r.score(), 150.5)
  assert_eq(r.status(), "T")
  assert_eq(r.date(), "01-Jan-2024")
  assert_eq(r.ipr_acc(), "IPR000001")
  assert_eq(r.ipr_desc(), "p53, DNA-binding domain")
  assert_eq(r.go_terms().length(), 2)
  assert_eq(r.go_terms()[0], "GO:0003674")
  assert_eq(r.go_terms()[1], "GO:0008150")
}

test "ips_record_accessors_empty_go" {
  let empty_go : Array[String] = Array::new()
  let r = @src.InterproScanRecord::new(
    protein_id="sp|Q9Y6K9|XYZ_HUMAN",
    md5="def789",
    seq_length=300,
    analysis="SMART",
    signature_acc="SM00100",
    signature_desc="domain",
    start_pos=10,
    stop_pos=200,
    score=80.0,
    status="T",
    date="02-Jan-2024",
    ipr_acc="IPR000002",
    ipr_desc="domain desc",
    go_terms=empty_go,
  )
  assert_eq(r.go_terms().length(), 0)
}

test "ips_record_to_string" {
  let go : Array[String] = Array::new()
  go.push("GO:0003674")
  let r = @src.InterproScanRecord::new(
    protein_id="sp|P12345|PROT_HUMAN",
    md5="abc123",
    seq_length=500,
    analysis="Pfam",
    signature_acc="PF00001",
    signature_desc="p53",
    start_pos=1,
    stop_pos=100,
    score=150.5,
    status="T",
    date="01-Jan-2024",
    ipr_acc="IPR000001",
    ipr_desc="p53 domain",
    go_terms=go,
  )
  let s = r.to_string()
  assert_true(s.contains("protein_id=sp|P12345|PROT_HUMAN"))
  assert_true(s.contains("analysis=Pfam"))
  assert_true(s.contains("signature=PF00001"))
  assert_true(s.contains("score=150.5"))
  assert_true(s.contains("ipr=IPR000001"))
  assert_true(s.contains("go_terms=GO:0003674"))
}

// ---------------------------------------------------------------------------
// TSV parsing
// ---------------------------------------------------------------------------

test "ips_parse_single_record" {
  let t = "\t"
  let data = "sp|P12345|PROT_HUMAN" + t + "abc123" + t + "500" + t + "Pfam" + t + "PF00001" + t + "p53" + t + "1" + t + "100" + t + "150.5" + t + "T" + t + "01-Jan-2024" + t + "IPR000001" + t + "p53 domain" + t + "GO:0003674\n"
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 1)
  assert_eq(records[0].protein_id(), "sp|P12345|PROT_HUMAN")
}

test "ips_parse_multiple_records" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 8)
}

test "ips_parse_all_fields_correct" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let r = records[0]
  assert_eq(r.protein_id(), "sp|P12345|PROT_HUMAN")
  assert_eq(r.md5(), "abc123def456")
  assert_eq(r.seq_length(), 500)
  assert_eq(r.analysis(), "Pfam")
  assert_eq(r.signature_acc(), "PF00001")
  assert_eq(r.signature_desc(), "p53 DNA-binding domain")
  assert_eq(r.start_pos(), 1)
  assert_eq(r.stop_pos(), 100)
  assert_eq(r.score(), 150.5)
  assert_eq(r.status(), "T")
  assert_eq(r.date(), "01-Jan-2024")
  assert_eq(r.ipr_acc(), "IPR000001")
  assert_eq(r.ipr_desc(), "p53, DNA-binding domain")
  assert_eq(r.go_terms().length(), 2)
  assert_eq(r.go_terms()[0], "GO:0003674")
  assert_eq(r.go_terms()[1], "GO:0008150")
}

test "ips_parse_multiple_go_terms" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  // Record 4 (index 3) has 3 GO terms.
  let r = records[3]
  assert_eq(r.go_terms().length(), 3)
  assert_eq(r.go_terms()[0], "GO:0005524")
  assert_eq(r.go_terms()[1], "GO:0006280")
  assert_eq(r.go_terms()[2], "GO:0003676")
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "ips_parse_empty_input" {
  let records = @src.parse_interproscan("")
  assert_eq(records.length(), 0)
}

test "ips_parse_only_comments" {
  let data = "# InterProScan output\n# version 5.0\n# another comment\n"
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 0)
}

test "ips_parse_blank_lines" {
  let data = "\n\n  \n\n"
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 0)
}

test "ips_parse_no_hits" {
  let t = "\t"
  let data = "sp|P99999|NOPROT_HUMAN" + t + "xyz789" + t + "100" + t + "No hits\n" +
    "sp|P12345|PROT_HUMAN" + t + "abc123" + t + "500" + t + "Pfam" + t + "PF00001" + t + "p53" + t + "1" + t + "100" + t + "150.5" + t + "T" + t + "01-Jan-2024" + t + "IPR000001" + t + "p53 domain" + t + "GO:0003674\n"
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 1)
  assert_eq(records[0].protein_id(), "sp|P12345|PROT_HUMAN")
}

test "ips_parse_only_no_hits" {
  let t = "\t"
  let data = "sp|P99999|NOPROT_HUMAN" + t + "xyz789" + t + "100" + t + "No hits\n"
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 0)
}

test "ips_parse_missing_ipr" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  // Record 3 (index 2) has "-" for IPR fields.
  let r = records[2]
  assert_eq(r.ipr_acc(), "")
  assert_eq(r.ipr_desc(), "")
}

test "ips_parse_missing_go_terms" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  // Record 3 (index 2) has "-" for GO terms.
  let r = records[2]
  assert_eq(r.go_terms().length(), 0)
}

test "ips_parse_missing_go_terms_dash" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  // Record 7 (index 6) has "-" for GO terms.
  let r = records[6]
  assert_eq(r.go_terms().length(), 0)
}

test "ips_parse_score_edge_cases" {
  let t = "\t"
  // Score is "-" → 0.0
  let data1 = "prot1" + t + "md5a" + t + "100" + t + "Pfam" + t + "PF01" + t + "desc" + t + "1" + t + "50" + t + "-" + t + "T" + t + "01-Jan-2024" + t + "-" + t + "-" + t + "-\n"
  let r1 = @src.parse_interproscan(data1)
  assert_eq(r1.length(), 1)
  assert_eq(r1[0].score(), 0.0)
  // Score is "NaN" → 0.0
  let data2 = "prot2" + t + "md5b" + t + "100" + t + "Pfam" + t + "PF02" + t + "desc" + t + "1" + t + "50" + t + "NaN" + t + "T" + t + "01-Jan-2024" + t + "-" + t + "-" + t + "-\n"
  let r2 = @src.parse_interproscan(data2)
  assert_eq(r2.length(), 1)
  assert_eq(r2[0].score(), 0.0)
  // Score is empty → 0.0
  let data3 = "prot3" + t + "md5c" + t + "100" + t + "Pfam" + t + "PF03" + t + "desc" + t + "1" + t + "50" + t + "" + t + "T" + t + "01-Jan-2024" + t + "-" + t + "-" + t + "-\n"
  let r3 = @src.parse_interproscan(data3)
  assert_eq(r3.length(), 1)
  assert_eq(r3[0].score(), 0.0)
}

test "ips_parse_fewer_columns" {
  let t = "\t"
  // Line with only 10 columns should be padded to 14.
  let data = "prot1" + t + "md5a" + t + "100" + t + "Pfam" + t + "PF01" + t + "desc" + t + "1" + t + "50" + t + "100.0" + t + "T\n"
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 1)
  assert_eq(records[0].date(), "")
  assert_eq(records[0].ipr_acc(), "")
  assert_eq(records[0].ipr_desc(), "")
  assert_eq(records[0].go_terms().length(), 0)
}

test "ips_parse_tsv_alias" {
  // parse_interproscan and parse_interproscan_tsv should produce the same results.
  let data = @src.interproscan_sample_data()
  let r1 = @src.parse_interproscan(data)
  let r2 = @src.parse_interproscan_tsv(data)
  assert_eq(r1.length(), r2.length())
  for i in 0..<r1.length() {
    assert_eq(r1[i].protein_id(), r2[i].protein_id())
    assert_eq(r1[i].analysis(), r2[i].analysis())
    assert_eq(r1[i].signature_acc(), r2[i].signature_acc())
  }
}

test "ips_parse_skips_comment_lines_in_sample" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  // Sample has 2 comment lines and 8 data lines.
  assert_eq(records.length(), 8)
  // None of the records should have analysis == "#".
  for r in records {
    assert_true(r.analysis() != "# InterProScan output (sample)")
  }
}

// ---------------------------------------------------------------------------
// Sample data validation
// ---------------------------------------------------------------------------

test "ips_sample_data_has_8_records" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  assert_eq(records.length(), 8)
}

test "ips_sample_data_proteins" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let proteins = @src.interproscan_get_unique_proteins(records)
  assert_eq(proteins.length(), 3)
  assert_eq(proteins[0], "sp|P12345|PROT_HUMAN")
  assert_eq(proteins[1], "sp|Q9Y6K9|XYZ_HUMAN")
  assert_eq(proteins[2], "sp|O15143|ABC_HUMAN")
}

test "ips_sample_data_signatures" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let sigs = @src.interproscan_get_unique_signatures(records)
  assert_eq(sigs.length(), 8)
  assert_eq(sigs[0], "PF00001")
  assert_eq(sigs[1], "SM00100")
  assert_eq(sigs[2], "PS00001")
}

// ---------------------------------------------------------------------------
// Filtering
// ---------------------------------------------------------------------------

test "ips_filter_by_analysis_pfam" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let pfam = @src.interproscan_filter_by_analysis(records, "Pfam")
  // Pfam records: PF00001, PF00002, PF00003 = 3 records.
  assert_eq(pfam.length(), 3)
  for r in pfam {
    assert_eq(r.analysis(), "Pfam")
  }
}

test "ips_filter_by_analysis_smart" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let smart = @src.interproscan_filter_by_analysis(records, "SMART")
  assert_eq(smart.length(), 2)
  for r in smart {
    assert_eq(r.analysis(), "SMART")
  }
}

test "ips_filter_by_analysis_no_match" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let none = @src.interproscan_filter_by_analysis(records, "NonExistent")
  assert_eq(none.length(), 0)
}

test "ips_filter_by_protein" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let prot1 = @src.interproscan_filter_by_protein(records, "sp|P12345|PROT_HUMAN")
  // P12345 has 3 records: Pfam, SMART, PROSITEPATTERNS.
  assert_eq(prot1.length(), 3)
  for r in prot1 {
    assert_eq(r.protein_id(), "sp|P12345|PROT_HUMAN")
  }
}

test "ips_filter_by_protein_multiple" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let prot3 = @src.interproscan_filter_by_protein(records, "sp|O15143|ABC_HUMAN")
  // O15143 has 3 records: Pfam, HAMMER, SMART.
  assert_eq(prot3.length(), 3)
}

test "ips_filter_by_protein_no_match" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let none = @src.interproscan_filter_by_protein(records, "sp|FAKE|FAKE_HUMAN")
  assert_eq(none.length(), 0)
}

test "ips_filter_on_empty_records" {
  let empty : Array[@src.InterproScanRecord] = Array::new()
  assert_eq(@src.interproscan_filter_by_analysis(empty, "Pfam").length(), 0)
  assert_eq(@src.interproscan_filter_by_protein(empty, "prot1").length(), 0)
}

// ---------------------------------------------------------------------------
// Unique value extraction
// ---------------------------------------------------------------------------

test "ips_get_unique_proteins" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let proteins = @src.interproscan_get_unique_proteins(records)
  assert_eq(proteins.length(), 3)
  assert_true(proteins.contains("sp|P12345|PROT_HUMAN"))
  assert_true(proteins.contains("sp|Q9Y6K9|XYZ_HUMAN"))
  assert_true(proteins.contains("sp|O15143|ABC_HUMAN"))
}

test "ips_get_unique_signatures" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let sigs = @src.interproscan_get_unique_signatures(records)
  assert_eq(sigs.length(), 8)
  assert_true(sigs.contains("PF00001"))
  assert_true(sigs.contains("SM00100"))
  assert_true(sigs.contains("PS00001"))
  assert_true(sigs.contains("PS50122"))
  assert_true(sigs.contains("PF00002"))
  assert_true(sigs.contains("PF00003"))
  assert_true(sigs.contains("HM00001"))
  assert_true(sigs.contains("SM00150"))
}

test "ips_get_unique_signatures_empty" {
  let empty : Array[@src.InterproScanRecord] = Array::new()
  let sigs = @src.interproscan_get_unique_signatures(empty)
  assert_eq(sigs.length(), 0)
}

test "ips_get_go_terms" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let go = @src.interproscan_get_go_terms(records)
  // GO:0003674 appears in records 1 and 2 (deduplicated).
  // Unique GO terms: GO:0003674, GO:0008150, GO:0005524, GO:0006280,
  //   GO:0003676, GO:0004672, GO:0016887, GO:0044183 = 8 unique.
  assert_eq(go.length(), 8)
  assert_true(go.contains("GO:0003674"))
  assert_true(go.contains("GO:0008150"))
  assert_true(go.contains("GO:0005524"))
  assert_true(go.contains("GO:0006280"))
  assert_true(go.contains("GO:0003676"))
  assert_true(go.contains("GO:0004672"))
  assert_true(go.contains("GO:0016887"))
  assert_true(go.contains("GO:0044183"))
}

test "ips_get_go_terms_empty" {
  let empty : Array[@src.InterproScanRecord] = Array::new()
  let go = @src.interproscan_get_go_terms(empty)
  assert_eq(go.length(), 0)
}

test "ips_get_go_terms_dedup" {
  // Two records with the same GO term should produce only one entry.
  let t = "\t"
  let data = "prot1" + t + "md5a" + t + "100" + t + "Pfam" + t + "PF01" + t + "desc" + t + "1" + t + "50" + t + "100.0" + t + "T" + t + "01-Jan-2024" + t + "IPR01" + t + "d1" + t + "GO:0003674\n" +
    "prot2" + t + "md5b" + t + "200" + t + "SMART" + t + "SM01" + t + "desc" + t + "1" + t + "50" + t + "100.0" + t + "T" + t + "01-Jan-2024" + t + "IPR02" + t + "d2" + t + "GO:0003674|GO:0008150\n"
  let records = @src.parse_interproscan(data)
  let go = @src.interproscan_get_go_terms(records)
  assert_eq(go.length(), 2)
  assert_eq(go[0], "GO:0003674")
  assert_eq(go[1], "GO:0008150")
}

// ---------------------------------------------------------------------------
// Grouping
// ---------------------------------------------------------------------------

test "ips_group_by_protein" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let groups = @src.interproscan_group_by_protein(records)
  assert_eq(groups.length(), 3)
  // First group: P12345 with 3 records.
  let (pid0, recs0) = groups[0]
  assert_eq(pid0, "sp|P12345|PROT_HUMAN")
  assert_eq(recs0.length(), 3)
  // Second group: Q9Y6K9 with 2 records.
  let (pid1, recs1) = groups[1]
  assert_eq(pid1, "sp|Q9Y6K9|XYZ_HUMAN")
  assert_eq(recs1.length(), 2)
  // Third group: O15143 with 3 records.
  let (pid2, recs2) = groups[2]
  assert_eq(pid2, "sp|O15143|ABC_HUMAN")
  assert_eq(recs2.length(), 3)
}

test "ips_group_by_protein_preserves_order" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let groups = @src.interproscan_group_by_protein(records)
  // First-seen order: P12345, Q9Y6K9, O15143.
  assert_eq(groups[0].0, "sp|P12345|PROT_HUMAN")
  assert_eq(groups[1].0, "sp|Q9Y6K9|XYZ_HUMAN")
  assert_eq(groups[2].0, "sp|O15143|ABC_HUMAN")
}

test "ips_group_by_protein_empty" {
  let empty : Array[@src.InterproScanRecord] = Array::new()
  let groups = @src.interproscan_group_by_protein(empty)
  assert_eq(groups.length(), 0)
}

test "ips_group_by_protein_single" {
  let t = "\t"
  let data = "prot1" + t + "md5a" + t + "100" + t + "Pfam" + t + "PF01" + t + "desc" + t + "1" + t + "50" + t + "100.0" + t + "T" + t + "01-Jan-2024" + t + "IPR01" + t + "d1" + t + "GO:0003674\n"
  let records = @src.parse_interproscan(data)
  let groups = @src.interproscan_group_by_protein(records)
  assert_eq(groups.length(), 1)
  assert_eq(groups[0].0, "prot1")
  assert_eq(groups[0].1.length(), 1)
}

// ---------------------------------------------------------------------------
// Summary
// ---------------------------------------------------------------------------

test "ips_summary_basic" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let s = @src.interproscan_summary(records)
  assert_true(s.contains("InterProScan Summary"))
  assert_true(s.contains("Total records: 8"))
  assert_true(s.contains("Unique proteins: 3"))
  assert_true(s.contains("Unique signatures: 8"))
  assert_true(s.contains("Unique analyses: 5"))
}

test "ips_summary_lists_analyses" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let s = @src.interproscan_summary(records)
  assert_true(s.contains("Pfam"))
  assert_true(s.contains("SMART"))
  assert_true(s.contains("PROSITEPATTERNS"))
  assert_true(s.contains("Profile"))
  assert_true(s.contains("HAMMER"))
}

test "ips_summary_empty" {
  let empty : Array[@src.InterproScanRecord] = Array::new()
  let s = @src.interproscan_summary(empty)
  assert_true(s.contains("Total records: 0"))
  assert_true(s.contains("Unique proteins: 0"))
  assert_true(s.contains("Unique signatures: 0"))
  assert_true(s.contains("Unique analyses: 0"))
}

// ---------------------------------------------------------------------------
// Full round-trip
// ---------------------------------------------------------------------------

test "ips_roundtrip_parse_and_filter" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  // Filter by Pfam, then check unique signatures.
  let pfam = @src.interproscan_filter_by_analysis(records, "Pfam")
  let pfam_sigs = @src.interproscan_get_unique_signatures(pfam)
  assert_eq(pfam_sigs.length(), 3)
  assert_true(pfam_sigs.contains("PF00001"))
  assert_true(pfam_sigs.contains("PF00002"))
  assert_true(pfam_sigs.contains("PF00003"))
}

test "ips_roundtrip_group_and_filter" {
  let data = @src.interproscan_sample_data()
  let records = @src.parse_interproscan(data)
  let groups = @src.interproscan_group_by_protein(records)
  // Filter the first group's records by Pfam.
  let (pid, recs) = groups[0]
  let pfam = @src.interproscan_filter_by_analysis(recs, "Pfam")
  assert_eq(pfam.length(), 1)
  assert_eq(pfam[0].signature_acc(), "PF00001")
}
