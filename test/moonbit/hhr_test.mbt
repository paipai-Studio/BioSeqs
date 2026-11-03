///|
/// Tests for the Biopython Bio.Align.hhr-compatible parser.

///|
fn hhr_minimal_text() -> String {
  "Query         mini_query\n" +
  "Match_columns 4\n" +
  "No_of_seqs    1 out of 1\n" +
  "Neff          1.0\n" +
  "Searched_HMMs 1\n" +
  " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM\n" +
  " 1 target_only 80.0 0.01 0.001 20.0 0.0 3 1-3 2-5 (5)\n" +
  "No 1\n" +
  ">target_only\n" +
  "Probab=80.0 E-value=0.01 Score=20.0 Aligned_cols=3 Identities=66.7% Similarity=0.5 Sum_probs=2.4\n" +
  "Q mini_query 1 AC-D 3 (4)\n" +
  "T target_only 2 ACGD 5 (5)\n" +
  "Done!\n"
}

///|
fn hhr_raises(text : String) -> Bool {
  try {
    ignore(@src.hhr_parse(text))
    false
  } catch {
    HhrError(_) => true
  }
}

///|
test "hhr parses file metadata" {
  let record = @src.hhr_parse(@src.hhr_sample_text())
  assert_eq(record.metadata.query_name, "demo_query")
  assert_eq(record.metadata.match_columns, 12)
  assert_eq(record.metadata.sequence_count, 24)
  assert_eq(record.metadata.sequence_count_total, 40)
  assert_true((record.metadata.neff - 3.5).abs() < 1.0e-12)
  match record.metadata.template_neff {
    Some(value) => assert_true((value - 2.1).abs() < 1.0e-12)
    None => abort("expected Template_Neff")
  }
  assert_eq(record.metadata.searched_hmms, 2)
  assert_eq(record.metadata.run_date, "Tue Aug 4 12:00:00 2026")
  assert_eq(record.metadata.command_line, "hhsearch -i demo.a3m -d pdb70")
}

///|
test "hhr parses ranked hit summaries" {
  let record = @src.hhr_parse(@src.hhr_sample_text())
  assert_eq(record.hit_summaries.length(), 2)
  let first = record.hit_summaries[0]
  assert_eq(first.rank, 1)
  assert_eq(first.target_name, "target_A")
  assert_eq(first.description, "Alpha beta enzyme")
  assert_true((first.probability - 99.8).abs() < 1.0e-12)
  assert_true((first.evalue - 1.2e-20).abs() < 1.0e-30)
  assert_true((first.pvalue - 4.0e-25).abs() < 1.0e-34)
  assert_eq(first.aligned_columns, 10)
  assert_eq(first.query_start, 0)
  assert_eq(first.query_end, 12)
  assert_eq(first.target_start, 4)
  assert_eq(first.target_end, 14)
  assert_eq(first.target_length, 40)
}

///|
test "hhr joins multi-block profile alignments" {
  let record = @src.hhr_parse(@src.hhr_sample_text())
  let alignment = record.alignments[0]
  assert_eq(alignment.query_sequence, "ACDEFGHIKLMN")
  assert_eq(alignment.target_sequence, "AC-EFGH-KLMN")
  assert_eq(alignment.query_consensus, "AcDEfGHiKLMn")
  assert_eq(alignment.target_consensus, "AC-EfGH-KLMN")
  assert_eq(alignment.alignment_length(), 12)
  assert_eq(alignment.aligned_columns, 10)
  assert_eq(alignment.query_start, 0)
  assert_eq(alignment.query_end, 12)
  assert_eq(alignment.target_start, 4)
  assert_eq(alignment.target_end, 14)
}

///|
test "hhr preserves per-column annotations" {
  let alignment = @src.hhr_parse(@src.hhr_sample_text()).alignments[0]
  assert_eq(alignment.query_secondary_structure, "CCHHHHHHEECC")
  assert_eq(alignment.target_secondary_structure, "CCHHHHH-EECC")
  assert_eq(alignment.target_dssp, "CCEEEEE-EECC")
  assert_eq(alignment.column_score, "||+|.||.||||")
  assert_eq(alignment.confidence, "998879976899")
}

///|
test "hhr exposes alignment statistics and sequence helpers" {
  let alignment = @src.hhr_parse(@src.hhr_sample_text()).alignments[0]
  assert_true((alignment.probability - 99.8).abs() < 1.0e-12)
  assert_true((alignment.evalue - 1.2e-20).abs() < 1.0e-30)
  assert_true((alignment.score - 85.4).abs() < 1.0e-12)
  assert_true((alignment.identities - 72.7).abs() < 1.0e-12)
  assert_true((alignment.similarity - 1.45).abs() < 1.0e-12)
  assert_true((alignment.sum_probabilities - 10.8).abs() < 1.0e-12)
  assert_eq(alignment.ungapped_query(), "ACDEFGHIKLMN")
  assert_eq(alignment.ungapped_target(), "ACEFGHKLMN")
  assert_eq(alignment.identity_count(), 10)
  assert_true((alignment.query_coverage() - 1.0).abs() < 1.0e-12)
  assert_true((alignment.target_coverage() - 0.25).abs() < 1.0e-12)
}

///|
test "hhr maps query coordinates through target gaps" {
  let alignment = @src.hhr_parse(@src.hhr_sample_text()).alignments[0]
  assert_eq(alignment.query_to_target(0), Some(4))
  assert_true(alignment.query_to_target(2) is None)
  assert_eq(alignment.query_to_target(3), Some(6))
  assert_true(alignment.query_to_target(7) is None)
  assert_eq(alignment.query_to_target(11), Some(13))
  assert_true(alignment.query_to_target(12) is None)
}

///|
test "hhr returns per-column coordinate pairs" {
  let pairs = @src.hhr_parse(@src.hhr_sample_text()).alignments[0].aligned_pairs()
  assert_eq(pairs.length(), 12)
  match pairs[0] {
    (Some(query), Some(target)) => {
      assert_eq(query, 0)
      assert_eq(target, 4)
    }
    _ => abort("expected paired first column")
  }
  match pairs[2] {
    (Some(query), None) => assert_eq(query, 2)
    _ => abort("expected target gap")
  }
  match pairs[11] {
    (Some(query), Some(target)) => {
      assert_eq(query, 11)
      assert_eq(target, 13)
    }
    _ => abort("expected paired final column")
  }
}

///|
test "hhr supports record queries and filters" {
  let record = @src.hhr_parse(@src.hhr_sample_text())
  assert_eq(record.num_hits(), 2)
  assert_true(record.get(-1) is None)
  assert_true(record.get(2) is None)
  match record.get(1) {
    Some(alignment) => assert_eq(alignment.target_name, "target_B")
    None => abort("expected second alignment")
  }
  match record.find_target("target_A") {
    Some(alignment) => assert_eq(alignment.rank, 1)
    None => abort("expected target_A")
  }
  assert_true(record.find_target("missing") is None)
  assert_eq(record.alignments_for_target("target_A").length(), 1)
  assert_eq(record.filter_by_probability(90.0).length(), 1)
  assert_eq(record.filter_by_evalue(0.01).length(), 2)
  match record.best_alignment() {
    Some(alignment) => assert_eq(alignment.rank, 1)
    None => abort("expected best alignment")
  }
}

///|
test "hhr keeps repeated target alignments" {
  let repeated = @src.hhr_sample_text().replace_all(
    old="target_B",
    new="target_A",
  )
  let record = @src.hhr_parse(repeated)
  assert_eq(record.alignments_for_target("target_A").length(), 2)
  match record.find_target("target_A") {
    Some(alignment) => assert_eq(alignment.rank, 1)
    None => abort("expected repeated target")
  }
}

///|
test "hhr parses official layout without blank separators" {
  let text = @src.hhr_sample_text()
  assert_true(!text.contains("pdb70\n\n"))
  assert_true(!text.contains("(20)\n\nNo 1"))
  assert_eq(@src.hhr_parse(text).num_hits(), 2)
}

///|
test "hhr accepts CRLF input" {
  let text = @src.hhr_sample_text().replace_all(old="\n", new="\r\n")
  let record = @src.hhr_parse(text)
  assert_eq(record.metadata.query_name, "demo_query")
  assert_eq(record.num_hits(), 2)
}

///|
test "hhr parses optional metadata and annotation-free alignments" {
  let record = @src.hhr_parse(hhr_minimal_text())
  assert_true(record.metadata.template_neff is None)
  assert_eq(record.metadata.run_date, "")
  let alignment = record.alignments[0]
  assert_eq(alignment.target_description, "")
  assert_eq(alignment.query_sequence, "AC-D")
  assert_eq(alignment.target_sequence, "ACGD")
  assert_eq(alignment.query_consensus, "")
  assert_eq(alignment.target_consensus, "")
  assert_eq(alignment.confidence, "")
  assert_eq(alignment.aligned_columns, 3)
  assert_true((alignment.query_coverage() - 0.75).abs() < 1.0e-12)
  assert_true((alignment.target_coverage() - 0.8).abs() < 1.0e-12)
}

///|
test "hhr parses zero-hit documents" {
  let text = "Query         empty_query\n" +
    "Match_columns 0\n" +
    "No_of_seqs    0 out of 0\n" +
    "Neff          0.0\n" +
    "Searched_HMMs 0\n" +
    " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM\n" +
    "Done!\n"
  let record = @src.hhr_parse(text)
  assert_eq(record.num_hits(), 0)
  assert_true(record.best_alignment() is None)
}

///|
test "hhr serialization round trip preserves the record" {
  let original = @src.hhr_parse(@src.hhr_sample_text())
  let serialized = original.to_string()
  assert_true(serialized.contains("Query         demo_query"))
  assert_true(serialized.contains("Done!"))
  let reparsed = @src.hhr_parse(serialized)
  assert_true(reparsed == original)
}

///|
test "hhr accepts EOF after a complete alignment" {
  let text = hhr_minimal_text().replace_all(old="Done!\n", new="")
  assert_eq(@src.hhr_parse(text).num_hits(), 1)
}

///|
test "hhr rejects empty input" {
  assert_true(hhr_raises(""))
}

///|
test "hhr requires core metadata" {
  let text = hhr_minimal_text().replace_all(old="Match_columns 4\n", new="")
  assert_true(hhr_raises(text))
}

///|
test "hhr rejects unknown header fields" {
  let text = hhr_minimal_text().replace_all(
    old="Match_columns 4\n",
    new="Unexpected    value\nMatch_columns 4\n",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr rejects malformed sequence counts" {
  let text = hhr_minimal_text().replace_all(old="1 out of 1", new="1 of 1")
  assert_true(hhr_raises(text))
}

///|
test "hhr validates the hit table header" {
  let text = hhr_minimal_text().replace_all(old="P-value", new="Pvalue")
  assert_true(hhr_raises(text))
}

///|
test "hhr requires contiguous summary ranks" {
  let text = @src.hhr_sample_text().replace_all(
    old=" 2 target_B",
    new=" 3 target_B",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr requires contiguous detail ranks" {
  let text = @src.hhr_sample_text().replace_all(old="No 2\n", new="No 3\n")
  assert_true(hhr_raises(text))
}

///|
test "hhr validates summary and detail counts" {
  let row = " 2 target_B Membrane protein          72.5 0.004 2.0E-6 31.2 0.0 8 3-10 2-9 (20)\n"
  let text = @src.hhr_sample_text().replace_all(old=row, new="")
  assert_true(hhr_raises(text))
}

///|
test "hhr validates summary and detail targets" {
  let text = @src.hhr_sample_text().replace_all(
    old=">target_A Alpha beta enzyme",
    new=">other Alpha beta enzyme",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr validates aligned column counts" {
  let text = @src.hhr_sample_text().replace_all(
    old="Aligned_cols=10",
    new="Aligned_cols=9",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr validates summary coordinates" {
  let text = @src.hhr_sample_text().replace_all(
    old="10 1-12 5-14 (40)",
    new="10 1-12 5-13 (40)",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr validates sequence coordinates" {
  let text = @src.hhr_sample_text().replace_all(
    old="T target_A       5 AC-EFG   9 (40)",
    new="T target_A       5 AC-EFG  10 (40)",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr validates multi-block coordinate continuity" {
  let text = @src.hhr_sample_text().replace_all(
    old="Q demo_query     7 HIKLMN",
    new="Q demo_query     8 HIKLMN",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr validates annotation widths" {
  let text = @src.hhr_sample_text().replace_all(
    old="Q ss_pred          CCHHHH",
    new="Q ss_pred          CCHHH",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr requires alignment statistics" {
  let statistics = "Probab=80.0 E-value=0.01 Score=20.0 Aligned_cols=3 Identities=66.7% Similarity=0.5 Sum_probs=2.4\n"
  assert_true(
    hhr_raises(hhr_minimal_text().replace_all(old=statistics, new="")),
  )
}

///|
test "hhr requires both alignment rows" {
  let text = hhr_minimal_text().replace_all(
    old="T target_only 2 ACGD 5 (5)\n",
    new="",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr rejects duplicate target headers" {
  let text = hhr_minimal_text().replace_all(
    old=">target_only\n",
    new=">target_only\n>target_only\n",
  )
  assert_true(hhr_raises(text))
}

///|
test "hhr rejects extra content after Done" {
  let text = hhr_minimal_text().replace_all(
    old="Done!\n",
    new="Done!\ntrailing data\n",
  )
  assert_true(hhr_raises(text))
}
