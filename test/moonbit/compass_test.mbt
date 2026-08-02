///|
/// Tests for Bio.Compass profile-profile alignment output parser.

// ---------------------------------------------------------------------------
// CompassRecord construction and accessors
// ---------------------------------------------------------------------------

test "compass_record_creation" {
  let r = @src.CompassRecord::new(
    "query1.msa",
    "template1.msa",
    5,
    4,
    120,
    110,
    185.5,
    3.5e-10,
    28.5,
    1,
    120,
    1,
    110,
    "MALKSLVRLFG",
    "MGVKSAVKT",
    ":  :: :::",
  )
  assert_eq(r.query_name(), "query1.msa")
  assert_eq(r.template_name(), "template1.msa")
  assert_eq(r.query_n_seqs(), 5)
  assert_eq(r.template_n_seqs(), 4)
  assert_eq(r.query_n_cols(), 120)
  assert_eq(r.template_n_cols(), 110)
  assert_eq(r.sw_score(), 185.5)
  assert_eq(r.e_value(), 3.5e-10)
  assert_eq(r.percentage_identity(), 28.5)
  assert_eq(r.query_start(), 1)
  assert_eq(r.query_end(), 120)
  assert_eq(r.template_start(), 1)
  assert_eq(r.template_end(), 110)
  assert_eq(r.aligned_query(), "MALKSLVRLFG")
  assert_eq(r.aligned_template(), "MGVKSAVKT")
  assert_eq(r.consensus_line(), ":  :: :::")
}

test "compass_record_accessors_individual" {
  let r = @src.CompassRecord::new(
    "q", "t", 1, 2, 10, 20, 50.0, 0.001, 35.0,
    5, 50, 10, 60, "ACGT", "ACGT", ":::",
  )
  assert_eq(r.query_name(), "q")
  assert_eq(r.template_name(), "t")
  assert_eq(r.query_n_seqs(), 1)
  assert_eq(r.template_n_seqs(), 2)
  assert_eq(r.query_n_cols(), 10)
  assert_eq(r.template_n_cols(), 20)
  assert_eq(r.sw_score(), 50.0)
  assert_eq(r.e_value(), 0.001)
  assert_eq(r.percentage_identity(), 35.0)
  assert_eq(r.query_start(), 5)
  assert_eq(r.query_end(), 50)
  assert_eq(r.template_start(), 10)
  assert_eq(r.template_end(), 60)
  assert_eq(r.aligned_query(), "ACGT")
  assert_eq(r.aligned_template(), "ACGT")
  assert_eq(r.consensus_line(), ":::")
}

test "compass_record_to_string" {
  let r = @src.CompassRecord::new(
    "query1.msa", "template1.msa",
    5, 4, 120, 110,
    185.5, 3.5e-10, 28.5,
    1, 120, 1, 110,
    "MALKSLVRLFG", "MGVKSAVKT", ":  :: :::",
  )
  let s = r.to_string()
  assert_true(s.contains("query1.msa"))
  assert_true(s.contains("template1.msa"))
  assert_true(s.contains("185.5"))
}

// ---------------------------------------------------------------------------
// Parsing: sample data
// ---------------------------------------------------------------------------

test "compass_parse_sample_data_count" {
  let text = @src.compass_sample_data()
  let records = @src.parse_compass(text)
  assert_eq(records.length(), 2)
}

test "compass_parse_first_record_query_name" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].query_name(), "query1.msa")
}

test "compass_parse_first_record_template_name" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].template_name(), "template1.msa")
}

test "compass_parse_first_record_n_seqs" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].query_n_seqs(), 5)
  assert_eq(records[0].template_n_seqs(), 4)
}

test "compass_parse_first_record_n_cols" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].query_n_cols(), 120)
  assert_eq(records[0].template_n_cols(), 110)
}

test "compass_parse_first_record_sw_score" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].sw_score(), 185.5)
}

test "compass_parse_first_record_e_value" {
  let records = @src.parse_compass(@src.compass_sample_data())
  // Use tolerance comparison due to floating-point precision
  assert_true((records[0].e_value() - 3.5e-10).abs() < 1.0e-20)
}

test "compass_parse_first_record_identity" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].percentage_identity(), 28.5)
}

test "compass_parse_first_record_alignment" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].aligned_query(), "MALKSLVRLFG")
  assert_eq(records[0].aligned_template(), "MGVKSAVKT")
}

test "compass_parse_first_record_positions" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].query_start(), 1)
  assert_eq(records[0].query_end(), 120)
  assert_eq(records[0].template_start(), 1)
  assert_eq(records[0].template_end(), 110)
}

test "compass_parse_first_record_consensus" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[0].consensus_line(), ":  :: :::")
}

test "compass_parse_second_record" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(records[1].query_name(), "query2.msa")
  assert_eq(records[1].template_name(), "template2.msa")
  assert_eq(records[1].query_n_seqs(), 3)
  assert_eq(records[1].template_n_seqs(), 6)
  assert_eq(records[1].query_n_cols(), 80)
  assert_eq(records[1].template_n_cols(), 90)
  assert_eq(records[1].sw_score(), 95.2)
  assert_eq(records[1].percentage_identity(), 15.3)
  assert_eq(records[1].aligned_query(), "ACDEFGHIKL")
  assert_eq(records[1].aligned_template(), "ACDEFGHIKL")
}

// ---------------------------------------------------------------------------
// Version extraction
// ---------------------------------------------------------------------------

test "compass_version_extraction" {
  let text = @src.compass_sample_data()
  let version = @src.compass_version(text)
  assert_eq(version, "2.4.2")
}

test "compass_version_empty_input" {
  assert_eq(@src.compass_version(""), "")
}

test "compass_version_no_version_line" {
  let text = "Some other content\nwithout version\n"
  assert_eq(@src.compass_version(text), "")
}

// ---------------------------------------------------------------------------
// Filtering
// ---------------------------------------------------------------------------

test "compass_filter_by_evalue" {
  let records = @src.parse_compass(@src.compass_sample_data())
  // E-values: 3.5e-10 and 1.2e-05
  let filtered = @src.compass_filter_by_evalue(records, 1.0e-06)
  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].query_name(), "query1.msa")
}

test "compass_filter_by_evalue_all_pass" {
  let records = @src.parse_compass(@src.compass_sample_data())
  let filtered = @src.compass_filter_by_evalue(records, 1.0)
  assert_eq(filtered.length(), 2)
}

test "compass_filter_by_evalue_none_pass" {
  let records = @src.parse_compass(@src.compass_sample_data())
  let filtered = @src.compass_filter_by_evalue(records, 1.0e-20)
  assert_eq(filtered.length(), 0)
}

test "compass_filter_by_identity" {
  let records = @src.parse_compass(@src.compass_sample_data())
  // Identities: 28.5 and 15.3
  let filtered = @src.compass_filter_by_identity(records, 20.0)
  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].query_name(), "query1.msa")
}

test "compass_filter_by_identity_all_pass" {
  let records = @src.parse_compass(@src.compass_sample_data())
  let filtered = @src.compass_filter_by_identity(records, 10.0)
  assert_eq(filtered.length(), 2)
}

// ---------------------------------------------------------------------------
// Alignment length
// ---------------------------------------------------------------------------

test "compass_alignment_length" {
  let records = @src.parse_compass(@src.compass_sample_data())
  assert_eq(@src.compass_alignment_length(records[0]), 11)
  assert_eq(@src.compass_alignment_length(records[1]), 10)
}

test "compass_alignment_length_empty" {
  let r = @src.CompassRecord::new(
    "", "", 0, 0, 0, 0, 0.0, 0.0, 0.0,
    0, 0, 0, 0, "", "", "",
  )
  assert_eq(@src.compass_alignment_length(r), 0)
}

// ---------------------------------------------------------------------------
// Summary
// ---------------------------------------------------------------------------

test "compass_summary_content" {
  let records = @src.parse_compass(@src.compass_sample_data())
  let summary = @src.compass_summary(records)
  assert_true(summary.contains("COMPASS Summary"))
  assert_true(summary.contains("Total records: 2"))
  assert_true(summary.contains("query1.msa"))
  assert_true(summary.contains("query2.msa"))
  assert_true(summary.contains("185.5"))
  assert_true(summary.contains("28.5"))
}

test "compass_summary_empty" {
  let records : Array[@src.CompassRecord] = []
  let summary = @src.compass_summary(records)
  assert_true(summary.contains("Total records: 0"))
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "compass_parse_empty_input" {
  let records = @src.parse_compass("")
  assert_eq(records.length(), 0)
}

test "compass_parse_whitespace_only" {
  let records = @src.parse_compass("   \n  \n  \n")
  assert_eq(records.length(), 0)
}

test "compass_parse_single_record" {
  let text = "COMPASS version 2.4.2\n" +
    "Query alignment: single.msa\n" +
    "Template alignment: tmpl.msa\n" +
    " Number of sequences in query alignment: 2\n" +
    " Number of sequences in template alignment: 3\n" +
    " Number of columns in query alignment: 50\n" +
    " Number of columns in template alignment: 60\n" +
    " Smith-Waterman score: 75.0\n" +
    " E-value: 1e-08\n" +
    " Percentage identity: 22.0\n" +
    " Query:      1 ACDEFG 50\n" +
    "                  : :  :\n" +
    " Template:   1 ACDEFG 60\n"
  let records = @src.parse_compass(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].query_name(), "single.msa")
  assert_eq(records[0].template_name(), "tmpl.msa")
  assert_eq(records[0].query_n_seqs(), 2)
  assert_eq(records[0].template_n_seqs(), 3)
  assert_eq(records[0].sw_score(), 75.0)
  assert_eq(records[0].aligned_query(), "ACDEFG")
  assert_eq(records[0].aligned_template(), "ACDEFG")
  assert_eq(records[0].consensus_line(), ": :  :")
}

test "compass_parse_missing_fields" {
  // Record with only version and query name (missing other fields)
  let text = "COMPASS version 1.0.0\n" +
    "Query alignment: partial.msa\n"
  let records = @src.parse_compass(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].query_name(), "partial.msa")
  assert_eq(records[0].template_name(), "")
  assert_eq(records[0].query_n_seqs(), 0)
  assert_eq(records[0].template_n_seqs(), 0)
  assert_eq(records[0].query_n_cols(), 0)
  assert_eq(records[0].template_n_cols(), 0)
  assert_eq(records[0].sw_score(), 0.0)
  assert_eq(records[0].e_value(), 0.0)
  assert_eq(records[0].percentage_identity(), 0.0)
  assert_eq(records[0].aligned_query(), "")
  assert_eq(records[0].aligned_template(), "")
  assert_eq(records[0].consensus_line(), "")
}

test "compass_parse_no_consensus_line" {
  // Alignment block without a consensus line between Query and Template
  let text = "COMPASS version 2.4.2\n" +
    "Query alignment: q.msa\n" +
    "Template alignment: t.msa\n" +
    " Number of sequences in query alignment: 1\n" +
    " Number of sequences in template alignment: 1\n" +
    " Number of columns in query alignment: 10\n" +
    " Number of columns in template alignment: 10\n" +
    " Smith-Waterman score: 50.0\n" +
    " E-value: 0.001\n" +
    " Percentage identity: 50.0\n" +
    " Query:      1 MALK 10\n" +
    " Template:   1 MGVK 10\n"
  let records = @src.parse_compass(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].aligned_query(), "MALK")
  assert_eq(records[0].aligned_template(), "MGVK")
  assert_eq(records[0].consensus_line(), "")
}

test "compass_parse_three_records" {
  let mut text = ""
  text = text + "COMPASS version 2.4.2\n"
  text = text + "Query alignment: q1.msa\n"
  text = text + "Template alignment: t1.msa\n"
  text = text + " E-value: 1e-05\n"
  text = text + "\n"
  text = text + "COMPASS version 2.4.2\n"
  text = text + "Query alignment: q2.msa\n"
  text = text + "Template alignment: t2.msa\n"
  text = text + " E-value: 1e-10\n"
  text = text + "\n"
  text = text + "COMPASS version 2.4.2\n"
  text = text + "Query alignment: q3.msa\n"
  text = text + "Template alignment: t3.msa\n"
  text = text + " E-value: 1e-02\n"
  let records = @src.parse_compass(text)
  assert_eq(records.length(), 3)
  assert_eq(records[0].query_name(), "q1.msa")
  assert_eq(records[1].query_name(), "q2.msa")
  assert_eq(records[2].query_name(), "q3.msa")
}

test "compass_parse_version_only" {
  let text = "COMPASS version 3.0.0\n"
  let records = @src.parse_compass(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].query_name(), "")
}

test "compass_sample_data_has_version" {
  let text = @src.compass_sample_data()
  assert_true(text.contains("COMPASS version 2.4.2"))
}

test "compass_sample_data_has_evalue" {
  let text = @src.compass_sample_data()
  assert_true(text.contains("3.5e-10"))
  assert_true(text.contains("1.2e-05"))
}

test "compass_filter_empty_records" {
  let records : Array[@src.CompassRecord] = []
  let by_evalue = @src.compass_filter_by_evalue(records, 1.0)
  let by_identity = @src.compass_filter_by_identity(records, 1.0)
  assert_eq(by_evalue.length(), 0)
  assert_eq(by_identity.length(), 0)
}
