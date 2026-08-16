///|
/// Tests for the Biopython Bio.SearchIO.InfernalIO-compatible parsers.

///|
fn infernal_format1_text() -> String {
  "# target name accession query name accession mdl mdl from mdl to seq from seq to strand trunc pass gc bias score E-value inc description of target\n" +
  "# ----------- --------- ---------- --------- --- -------- -------- -------- -------- ------ ----- ---- -- ---- ----- ------- --- ---------------------\n" +
  "seq1 SACC RFQ1 QACC cm 1 70 10 79 + no 1 0.50 0.0 20.0 1.0e-05 ! format one target\n"
}

///|
fn infernal_format2_text() -> String {
  "# idx target name accession query name accession clan name mdl mdl from mdl to seq from seq to strand trunc pass gc bias score E-value inc olp anyidx afrct1 afrct2 winidx wfrct1 wfrct2 mdl len seq len description of target\n" +
  "# --- ----------- --------- ---------- --------- --------- --- -------- -------- -------- -------- ------ ----- ---- -- ---- ----- ------- --- --- ------ ------ ------ ------ ------ ------ ------- ------- ---------------------\n" +
  "1 chr2 ACC2 RFQ2 RFQ2.1 CL001 cm 2 60 900 841 - 5' 3 0.44 0.2 55.0 3.0e-18 ! ^ 7 0.25 0.50 2 0.75 0.80 65 2000 format two target\n"
}

///|
fn infernal_noali_text() -> String {
  "# cmsearch :: search a covariance model against a sequence database\n" +
  "# INFERNAL 1.1.5 (Sep 2023)\n" +
  "# - - -\n" +
  "# target sequence database: noali.fa\n" +
  "# show alignments in output: no\n" +
  "# - - -\n" +
  "Query: RFNOALI [M=50]\n" +
  "Accession: RFNOALI.1\n" +
  "Description: score-only query\n" +
  "Hit scores:\n" +
  " ----   --------- ------ ----- ---- ---------------- -------- -------- --- --- ----- -- --------------\n" +
  " 1 ! 1.0e-08 33.0 0.1 seqA 100 81 - cm no 0.45 reverse score hit\n" +
  " 2 ? 2.0e-04 20.0 0.0 seqA 200 219 + hmm 3' 0.51 second domain\n" +
  " 3 ! 5.0e-03 18.0 0.2 seqB 10 29 + cm no 0.39 other target\n" +
  "Internal CM pipeline statistics summary:\n" +
  "//\n" +
  "[ok]\n"
}

///|
fn infernal_no_hit_text() -> String {
  "# cmscan :: search a sequence against a covariance model database\n" +
  "# INFERNAL 1.1.5 (Sep 2023)\n" +
  "# - - -\n" +
  "# target CM database: Rfam.cm\n" +
  "# show alignments in output: yes\n" +
  "# - - -\n" +
  "Query: empty_query [L=25]\n" +
  "Hit scores:\n" +
  "   [No hits detected that satisfy reporting thresholds]\n" +
  "Internal CM pipeline statistics summary:\n" +
  "//\n" +
  "[ok]\n"
}

///|
fn infernal_reverse_text() -> String {
  @src.infernal_text_sample()
  .replace_all(old="100 118 + .. 0.96 no 0.50", new="200 182 - .. 0.96 no 0.50")
  .replace_all(
    old="targetA 100 ACGUACGUGGGGGAACCGG 118",
    new="targetA 200 ACGUACGUGGGGGAACCGG 182",
  )
}

///|
fn infernal_duplicate_alignment_text() -> String {
  let second = ">> targetA  synthetic RNA target\n" +
    " ----   --------- ------ ----- ---- --- -------- -------- --- -------- -------- --- --- ---- ----- --\n" +
    " 2 ? 4.0e-06 29.0 0.0 hmm 2 18 .. 300 316 + .. 0.88 3' 0.42\n" +
    "CS ((..............))\n" +
    "RFTEST 2 CGUACGUAACCGGAAAAC 18\n" +
    "targetA 300 CGUACGUAACCGGAAAAC 316\n" +
    "PP 88888888888888888\n"
  @src.infernal_text_sample().replace_all(
    old="Internal CM pipeline statistics summary:\n",
    new=second + "Internal CM pipeline statistics summary:\n",
  )
}

///|
fn infernal_tab_raises(content : String, format : Int) -> Bool {
  try {
    ignore(@src.infernal_parse_tabular(content, format~))
    false
  } catch {
    InfernalError(_) => true
  }
}

///|
fn infernal_text_raises(content : String) -> Bool {
  try {
    ignore(@src.infernal_parse_text(content))
    false
  } catch {
    InfernalError(_) => true
  }
}

///|
test "InfernalIO detects tabular formats 1, 2, and 3" {
  assert_true(@src.infernal_tabular_format(infernal_format1_text()) is Format1)
  assert_true(@src.infernal_tabular_format(infernal_format2_text()) is Format2)
  assert_true(
    @src.infernal_tabular_format(@src.infernal_tabular_sample()) is Format3,
  )
}

///|
test "InfernalIO supports an explicit tabular format without a header" {
  let row = "seq1 SACC RFQ1 QACC cm 1 70 10 79 + no 1 0.50 0.0 20.0 1.0e-05 ! explicit\n"
  let results = @src.infernal_parse_tabular(row, format=1)
  assert_eq(results.length(), 1)
  assert_eq(results[0].id, "RFQ1")
}

///|
test "InfernalIO parses format 1 core fields" {
  let query = @src.infernal_parse_tabular(infernal_format1_text())[0]
  assert_eq(query.id, "RFQ1")
  assert_eq(query.accession, "QACC")
  assert_eq(query.sequence_length, 0)
  assert_eq(query.hits[0].id, "seq1")
  assert_eq(query.hits[0].accession, "SACC")
  assert_eq(query.hits[0].description, "format one target")
  assert_eq(query.hits[0].hsps[0].model, "cm")
}

///|
test "InfernalIO parses format 2 clan lengths and overlap fields" {
  let query = @src.infernal_parse_tabular(infernal_format2_text())[0]
  let hit = query.hits[0]
  let hsp = hit.hsps[0]
  assert_eq(query.clan, "CL001")
  assert_eq(query.sequence_length, 2000)
  assert_eq(hit.model_length, 65)
  assert_eq(hsp.overlap, "^")
  assert_eq(hsp.any_index, Some(7))
  assert_eq(hsp.winner_index, Some(2))
  assert_eq(hsp.pipeline_pass, 3)
}

///|
test "InfernalIO parses format 2 overlap fractions" {
  let hsp = @src.infernal_parse_tabular(infernal_format2_text())[0].hits[0].hsps[0]
  assert_eq(hsp.any_overlap_fraction, Some(0.25))
  assert_eq(hsp.reciprocal_overlap_fraction, Some(0.50))
  assert_eq(hsp.winner_overlap_fraction, Some(0.75))
  assert_eq(hsp.winner_reciprocal_fraction, Some(0.80))
}

///|
test "InfernalIO accepts missing format 2 overlap values" {
  let text = infernal_format2_text().replace_all(
    old="^ 7 0.25 0.50 2 0.75 0.80",
    new="= - - - \" \" \"",
  )
  let hsp = @src.infernal_parse_tabular(text)[0].hits[0].hsps[0]
  assert_eq(hsp.overlap, "=")
  assert_true(hsp.any_index is None)
  assert_true(hsp.any_overlap_fraction is None)
  assert_true(hsp.winner_index is None)
  assert_true(hsp.winner_overlap_fraction is None)
}

///|
test "InfernalIO parses format 3 lengths and descriptions" {
  let query = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0]
  assert_eq(query.sequence_length, 1000)
  assert_eq(query.hits[0].model_length, 71)
  assert_eq(query.hits[0].description, "bacterial 5S RNA locus")
}

///|
test "InfernalIO normalizes forward coordinates to zero-based half-open" {
  let fragment = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_start, 0)
  assert_eq(fragment.query_end, 71)
  assert_eq(fragment.hit_start, 100)
  assert_eq(fragment.hit_end, 171)
  assert_true(fragment.hit_strand == @src.strand_plus())
}

///|
test "InfernalIO normalizes reverse coordinates and retains strand" {
  let fragment = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0].hits[0].hsps[1].fragments[0]
  assert_eq(fragment.query_start, 1)
  assert_eq(fragment.query_end, 70)
  assert_eq(fragment.hit_start, 411)
  assert_eq(fragment.hit_end, 480)
  assert_true(fragment.hit_strand == @src.strand_minus())
}

///|
test "InfernalIO aggregates repeated target rows into HSPs" {
  let query = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0]
  assert_eq(query.hits.length(), 1)
  assert_eq(query.hits[0].hsps.length(), 2)
  assert_eq(query.count_hsps(), 2)
}

///|
test "InfernalIO preserves query and hit first-seen order" {
  let base = @src.infernal_tabular_sample()
  let extra = "chrC RFSEQ3 RF00001 RF00001 cm 1 20 700 719 + no 1 0.40 0.0 10.0 0.01 ? 71 1000 later hit\n"
  let text = base + extra
  let results = @src.infernal_parse_tabular(text)
  assert_eq(results.length(), 2)
  assert_eq(results[0].id, "RF00001")
  assert_eq(results[1].id, "RF00002")
  assert_eq(results[0].hits[0].id, "chrA")
  assert_eq(results[0].hits[1].id, "chrC")
}

///|
test "InfernalIO parses inclusion and scientific E-values" {
  let hsps = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0].hits[0].hsps
  assert_true(hsps[0].is_included)
  assert_false(hsps[1].is_included)
  assert_true((hsps[0].evalue - 2.0e-14).abs() < 1.0e-24)
}

///|
test "InfernalIO rejects unsupported explicit tabular formats" {
  assert_true(infernal_tab_raises(infernal_format1_text(), 4))
}

///|
test "InfernalIO requires a header for automatic tabular detection" {
  let row = "seq1 SACC RFQ1 QACC cm 1 70 10 79 + no 1 0.50 0.0 20.0 1.0e-05 ! explicit\n"
  assert_true(infernal_tab_raises(row, 0))
}

///|
test "InfernalIO rejects short tabular rows" {
  let text = "# target name accession query name accession mdl description of target\nshort row\n"
  assert_true(infernal_tab_raises(text, 1))
}

///|
test "InfernalIO rejects invalid integer fields" {
  let text = infernal_format1_text().replace_all(old="cm 1 70", new="cm x 70")
  assert_true(infernal_tab_raises(text, 0))
}

///|
test "InfernalIO rejects invalid floating-point fields" {
  let text = infernal_format1_text().replace_all(
    old="20.0 1.0e-05",
    new="bad 1.0e-05",
  )
  assert_true(infernal_tab_raises(text, 0))
}

///|
test "InfernalIO rejects invalid strand fields" {
  let text = infernal_format1_text().replace_all(
    old="10 79 + no",
    new="10 79 ? no",
  )
  assert_true(infernal_tab_raises(text, 0))
}

///|
test "InfernalIO rejects non-positive coordinates" {
  let text = infernal_format1_text().replace_all(old="cm 1 70", new="cm 0 70")
  assert_true(infernal_tab_raises(text, 0))
}

///|
test "InfernalIO parses plain-text metadata and query fields" {
  let query = @src.infernal_parse_text(@src.infernal_text_sample())[0]
  assert_eq(query.id, "RFTEST")
  assert_eq(query.accession, "RFTEST.1")
  assert_eq(query.description, "synthetic structured RNA")
  assert_eq(query.sequence_length, 19)
  assert_eq(query.program, "cmsearch")
  assert_eq(query.version, "1.1.5")
  assert_eq(query.target, "transcripts.fa")
}

///|
test "InfernalIO parses plain-text CM HSP statistics" {
  let hsp = @src.infernal_parse_text(@src.infernal_text_sample())[0].hits[0].hsps[0]
  assert_eq(hsp.model, "cm")
  assert_eq(hsp.truncated, "no")
  assert_eq(hsp.query_end_type, "[]")
  assert_eq(hsp.hit_end_type, "..")
  assert_true((hsp.bitscore - 45.0).abs() < 1.0e-12)
  assert_true((hsp.average_accuracy - 0.96).abs() < 1.0e-12)
  assert_true((hsp.gc - 0.50).abs() < 1.0e-12)
}

///|
test "InfernalIO parses HMM-only pipeline output" {
  let text = @src.infernal_text_sample()
    .replace_all(old="0.1 cm 1 19", new="0.1 hmm 1 19")
    .replace_all(
      old="Internal CM pipeline statistics summary:",
      new="Internal HMM-only pipeline statistics summary:",
    )
  let hsp = @src.infernal_parse_text(text)[0].hits[0].hsps[0]
  assert_eq(hsp.model, "hmm")
  assert_eq(hsp.fragments.length(), 2)
}

///|
test "InfernalIO preserves alignment sequences and annotations" {
  let fragments = @src.infernal_parse_text(@src.infernal_text_sample())[0].hits[0].hsps[0].fragments
  assert_eq(fragments[0].query_sequence, "ACGUACGU")
  assert_eq(fragments[0].hit_sequence, "ACGUACGU")
  assert_eq(fragments[0].consensus_structure, "(((....)")
  assert_eq(fragments[0].posterior_probability, "99999999")
  assert_eq(fragments[1].query_sequence, "AACCGG")
  assert_eq(fragments[1].hit_sequence, "AACCGG")
}

///|
test "InfernalIO splits local-end markers into multiple fragments" {
  let hsp = @src.infernal_parse_text(@src.infernal_text_sample())[0].hits[0].hsps[0]
  assert_eq(hsp.fragments.length(), 2)
  assert_eq(hsp.fragments[0].model_omission_before, 0)
  assert_eq(hsp.fragments[1].model_omission_before, 5)
  assert_eq(hsp.fragments[1].sequence_omission_before, 5)
  assert_eq(hsp.alignment_span(), 14)
}

///|
test "InfernalIO advances coordinates across local ends" {
  let fragments = @src.infernal_parse_text(@src.infernal_text_sample())[0].hits[0].hsps[0].fragments
  assert_eq(fragments[0].query_start, 0)
  assert_eq(fragments[0].query_end, 8)
  assert_eq(fragments[0].hit_start, 99)
  assert_eq(fragments[0].hit_end, 107)
  assert_eq(fragments[1].query_start, 13)
  assert_eq(fragments[1].query_end, 19)
  assert_eq(fragments[1].hit_start, 112)
  assert_eq(fragments[1].hit_end, 118)
}

///|
test "InfernalIO advances reverse-strand local-end coordinates downward" {
  let hsp = @src.infernal_parse_text(infernal_reverse_text())[0].hits[0].hsps[0]
  assert_eq(hsp.fragments.length(), 2)
  assert_eq(hsp.fragments[0].hit_start, 192)
  assert_eq(hsp.fragments[0].hit_end, 200)
  assert_eq(hsp.fragments[1].hit_start, 181)
  assert_eq(hsp.fragments[1].hit_end, 187)
  assert_true(hsp.fragments[1].hit_strand == @src.strand_minus())
  assert_eq(hsp.hit_start(), 181)
  assert_eq(hsp.hit_end(), 200)
}

///|
test "InfernalIO parses noali score tables" {
  let query = @src.infernal_parse_text(infernal_noali_text())[0]
  assert_eq(query.id, "RFNOALI")
  assert_eq(query.target, "noali.fa")
  assert_eq(query.hits.length(), 2)
  assert_eq(query.count_hsps(), 3)
  let fragment = query.hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_start, -1)
  assert_eq(fragment.hit_start, 80)
  assert_eq(fragment.hit_end, 100)
  assert_true(fragment.hit_strand == @src.strand_minus())
}

///|
test "InfernalIO merges duplicate noali targets and retains HSP order" {
  let hit = @src.infernal_parse_text(infernal_noali_text())[0].hits[0]
  assert_eq(hit.id, "seqA")
  assert_eq(hit.description, "reverse score hit")
  assert_eq(hit.hsps.length(), 2)
  assert_true((hit.hsps[0].evalue - 1.0e-08).abs() < 1.0e-18)
  assert_eq(hit.hsps[1].model, "hmm")
}

///|
test "InfernalIO merges repeated alignment target sections" {
  let query = @src.infernal_parse_text(infernal_duplicate_alignment_text())[0]
  assert_eq(query.hits.length(), 1)
  assert_eq(query.hits[0].hsps.length(), 2)
  assert_eq(query.hits[0].hsps[1].model, "hmm")
}

///|
test "InfernalIO returns empty hits for no-hit queries" {
  let query = @src.infernal_parse_text(infernal_no_hit_text())[0]
  assert_eq(query.id, "empty_query")
  assert_eq(query.sequence_length, 25)
  assert_eq(query.hits.length(), 0)
}

///|
test "InfernalIO parses multiple plain-text queries" {
  let second = @src.infernal_text_sample()
    .replace_all(old="RFTEST", new="RFSECOND")
    .replace_all(old="targetA", new="targetB")
  let text = @src.infernal_text_sample().replace_all(old="[ok]\n", new="") +
    second
  let results = @src.infernal_parse_text(text)
  assert_eq(results.length(), 2)
  assert_eq(results[0].id, "RFTEST")
  assert_eq(results[1].id, "RFSECOND")
}

///|
test "InfernalIO chooses best HSP by E-value then score" {
  let hit = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0].hits[0]
  match hit.best_hsp() {
    Some(hsp) => {
      assert_true((hsp.evalue - 2.0e-14).abs() < 1.0e-24)
      assert_true((hsp.bitscore - 48.5).abs() < 1.0e-12)
    }
    None => abort("expected a best Infernal HSP")
  }
}

///|
test "InfernalIO filters by E-value and inclusion state immutably" {
  let query = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0]
  let filtered = query.filter(maximum_evalue=1.0e-10, included_only=true)
  assert_eq(filtered.hits.length(), 1)
  assert_eq(filtered.count_hsps(), 1)
  assert_eq(query.count_hsps(), 2)
  assert_true(filtered.hits[0].hsps[0].is_included)
}

///|
test "InfernalIO rejects invalid filter thresholds" {
  let query = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0]
  let raised = try {
    ignore(query.filter(maximum_evalue=-1.0))
    false
  } catch {
    InfernalError(_) => true
  }
  assert_true(raised)
}

///|
test "InfernalIO converts to the generic SearchIO hierarchy" {
  let source = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0]
  let query = source.to_searchio()
  assert_eq(query.id, "RF00001")
  assert_eq(query.seq_len, 1000)
  assert_eq(query.hits.length(), 1)
  assert_eq(query.hits[0].id, "chrA")
  assert_eq(query.hits[0].accession, "RFSEQ1")
  assert_eq(query.hits[0].hsps.length(), 2)
  assert_eq(query.hits[0].hsps[1].fragments[0].hit_start, 411)
  assert_true(
    query.hits[0].hsps[1].fragments[0].hit_strand == @src.strand_minus(),
  )
}

///|
test "InfernalIO summarizes query hierarchy counts" {
  let query = @src.infernal_parse_tabular(@src.infernal_tabular_sample())[0]
  assert_eq(
    query.summary(),
    "InfernalQueryResult(query=RF00001, hits=1, hsps=2, program=)",
  )
}

///|
test "InfernalIO rejects malformed text alignment score rows" {
  let malformed = @src.infernal_text_sample().replace_all(
    old=" 1 ! 2.0e-10 45.0 0.1 cm 1 19 [] 100 118 + .. 0.96 no 0.50\n",
    new=" 1 ! 2.0e-10 45.0 0.1 cm 1 19 [] 100 118 + .. 0.96 no\n",
  )
  assert_true(infernal_text_raises(malformed))
}
