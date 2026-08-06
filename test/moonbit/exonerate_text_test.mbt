// Black-box tests for Bio.SearchIO.ExonerateIO.exonerate_text support.

///|
fn exonerate_text_test_parse(text : String) -> @src.ExonerateTextDocument {
  @src.exonerate_text_parse(text) catch {
    ExonerateTextError(message) =>
      abort("valid Exonerate text failed: " + message)
  }
}

///|
fn exonerate_text_test_rejects(text : String) -> Bool {
  try {
    ignore(@src.exonerate_text_parse(text))
    false
  } catch {
    ExonerateTextError(_) => true
  }
}

///|
fn exonerate_text_test_wrap(
  query : String,
  query_description : String,
  target : String,
  target_description : String,
  model : String,
  score : Int,
  query_start : Int,
  query_end : Int,
  target_start : Int,
  target_end : Int,
  body : String,
) -> String {
  "Command line: [exonerate -m " +
  model +
  " query.fa target.fa]\n" +
  "Hostname: [test-host]\n\n" +
  exonerate_text_test_alignment(
    query, query_description, target, target_description, model, score, query_start,
    query_end, target_start, target_end, body,
  ) +
  "-- completed exonerate analysis\n"
}

///|
fn exonerate_text_test_alignment(
  query : String,
  query_description : String,
  target : String,
  target_description : String,
  model : String,
  score : Int,
  query_start : Int,
  query_end : Int,
  target_start : Int,
  target_end : Int,
  body : String,
) -> String {
  "C4 Alignment:\n" +
  "------------\n" +
  "         Query: " +
  query +
  (if query_description.length() > 0 { " " + query_description } else { "" }) +
  "\n" +
  "        Target: " +
  target +
  (if target_description.length() > 0 { " " + target_description } else { "" }) +
  "\n" +
  "         Model: " +
  model +
  "\n" +
  "     Raw score: " +
  score.to_string() +
  "\n" +
  "   Query range: " +
  query_start.to_string() +
  " -> " +
  query_end.to_string() +
  "\n" +
  "  Target range: " +
  target_start.to_string() +
  " -> " +
  target_end.to_string() +
  "\n\n" +
  body +
  "\n"
}

///|
fn exonerate_text_test_dna_report() -> String {
  exonerate_text_test_wrap(
    "query1",
    "demo query",
    "target1",
    "demo target",
    "affine:local:dna2dna",
    31,
    0,
    6,
    100,
    107,
    "   1 : ACGT-AC :    6\n" + "       |||| ||\n" + " 101 : ACGTTAC :  107\n",
  )
}

///|
fn exonerate_text_test_reverse_report() -> String {
  exonerate_text_test_wrap(
    "query1",
    "",
    "target1",
    "reverse target:[revcomp]",
    "ungapped:dna2dna",
    18,
    0,
    6,
    200,
    193,
    "   1 : ACGTAC- :    6\n" + "       |||||| \n" + " 200 : ACGTACT :  194\n",
  )
}

///|
fn exonerate_text_test_center(value : String, width : Int) -> String {
  let left = (width - value.length()) / 2
  " ".repeat(left) + value + " ".repeat(width - value.length() - left)
}

///|
fn exonerate_text_test_intron_report() -> String {
  let marker = "  >>>> Target Intron 1 >>>>  "
  let splice = "gt" + ".".repeat(marker.length() - 4) + "ag"
  let query_row = "ACG" + marker + "TAC"
  let hit_row = "ACG" + splice + "TAC"
  let similarity = "|||" +
    exonerate_text_test_center("7 bp", marker.length()) +
    "|||"
  exonerate_text_test_wrap(
    "transcript",
    "",
    "chromosome",
    "",
    "est2genome",
    42,
    0,
    6,
    100,
    113,
    "   1 : " +
    query_row +
    " :    6\n" +
    "       " +
    similarity +
    "\n" +
    " 101 : " +
    hit_row +
    " :  113\n",
  )
}

///|
fn exonerate_text_test_ner_report() -> String {
  let query_marker = "--<   4   >--"
  let hit_marker = "--<   6   >--"
  exonerate_text_test_wrap(
    "query",
    "",
    "target",
    "",
    "NER:affine:local:dna2dna",
    25,
    0,
    10,
    100,
    112,
    "   1 : AAA" +
    query_marker +
    "CCC :   10\n" +
    "       |||--< NER 1 >--|||\n" +
    " 101 : AAA" +
    hit_marker +
    "CCC :  112\n",
  )
}

///|
fn exonerate_text_test_protein_report() -> String {
  exonerate_text_test_wrap(
    "protein_query",
    "",
    "dna_target",
    "",
    "protein2dna:local",
    55,
    0,
    3,
    100,
    109,
    "   1 : MetGly<->Lys :    3\n" +
    "       ||||||   |||\n" +
    "       MetGly<->Lys\n" +
    " 101 : ATGGGC---AAA :  109\n",
  )
}

///|
fn exonerate_text_test_dna_to_protein_report() -> String {
  exonerate_text_test_wrap(
    "dna_query",
    "",
    "protein_target",
    "",
    "ungapped:dna2protein",
    44,
    0,
    9,
    10,
    13,
    "   1 : ATGGGCAAA :    9\n" +
    "       MetGlyLys\n" +
    "       |||||||||\n" +
    "  11 : MetGlyLys :   13\n",
  )
}

///|
fn exonerate_text_test_frameshift_report() -> String {
  exonerate_text_test_wrap(
    "protein_query",
    "",
    "dna_target",
    "",
    "protein2dna:local",
    27,
    0,
    2,
    100,
    108,
    "   1 : Asp--Ile :    2\n" +
    "       |||##|||\n" +
    "       Asp##Ile\n" +
    " 101 : GATCCATT :  108\n",
  )
}

///|
fn exonerate_text_test_wrapped_report() -> String {
  exonerate_text_test_wrap(
    "wrapped_query",
    "",
    "wrapped_target",
    "",
    "affine:local:dna2dna",
    36,
    0,
    8,
    100,
    108,
    "   1 : ACGT :    4\n" +
    "       ||||\n" +
    " 101 : ACGT :  104\n\n" +
    "   5 : TGCA :    8\n" +
    "       ||||\n" +
    " 105 : TGCA :  108\n",
  )
}

///|
fn exonerate_text_test_joint_intron_report() -> String {
  let marker = "  >>>> Joint Intron 1 >>>>  "
  let similarity = "|||" +
    exonerate_text_test_center("5 bp // 7 bp", marker.length()) +
    "|||"
  exonerate_text_test_wrap(
    "joint_query",
    "",
    "joint_target",
    "",
    "est2genome",
    48,
    0,
    11,
    100,
    113,
    "   1 : AAA" +
    marker +
    "CCC :   11\n" +
    "       " +
    similarity +
    "\n" +
    " 101 : AAA" +
    marker +
    "CCC :  113\n",
  )
}

///|
fn exonerate_text_test_reverse_intron_report() -> String {
  let marker = "  >>>> Target Intron 1 >>>>  "
  let splice = "gt" + ".".repeat(marker.length() - 4) + "ag"
  let similarity = "|||" +
    exonerate_text_test_center("7 bp", marker.length()) +
    "|||"
  exonerate_text_test_wrap(
    "reverse_query",
    "query strand:[revcomp]",
    "reverse_target",
    "target strand:[revcomp]",
    "est2genome",
    41,
    10,
    4,
    200,
    187,
    "  10 : ACG" +
    marker +
    "TAC :    5\n" +
    "       " +
    similarity +
    "\n" +
    " 200 : ACG" +
    splice +
    "TAC :  188\n",
  )
}

///|
fn exonerate_text_test_split_codon_report() -> String {
  let marker = "  >>>> Target Intron 1 >>>>  "
  let splice = "gt" + ".".repeat(marker.length() - 4) + "ag"
  let similarity = "|||{||}" +
    exonerate_text_test_center("7 bp", marker.length()) +
    "{|}|||"
  exonerate_text_test_wrap(
    "protein_query",
    "",
    "genome_target",
    "",
    "protein2genome:local",
    73,
    0,
    3,
    100,
    116,
    "   1 : Gly{Th}" +
    marker +
    "{r}Ala :    3\n" +
    "       " +
    similarity +
    "\n" +
    "       Gly{Th}" +
    " ".repeat(marker.length()) +
    "{r}Ala\n" +
    " 101 : GGT{AC}" +
    splice +
    "{G}GCT :  116\n",
  )
}

///|
fn exonerate_text_test_coding_report() -> String {
  exonerate_text_test_wrap(
    "coding_query",
    "",
    "coding_target",
    "",
    "coding2coding",
    52,
    0,
    6,
    100,
    106,
    "   1 : ATGGAA :    6\n" +
    "       MetGlu\n" +
    "       ||||||\n" +
    "       MetGlu\n" +
    " 101 : ATGGAA :  106\n",
  )
}

///|
fn exonerate_text_test_coding_frameshift_report() -> String {
  exonerate_text_test_wrap(
    "coding_query",
    "",
    "coding_target",
    "",
    "coding2coding",
    39,
    0,
    6,
    100,
    106,
    "   1 : ATGAC-A :    6\n" +
    "       Met#Thr\n" +
    "       |||#|||\n" +
    "       Met-Thr\n" +
    " 101 : ATG-CTA :  106\n",
  )
}

///|
fn exonerate_text_test_special_protein_report() -> String {
  exonerate_text_test_wrap(
    "special_protein",
    "",
    "special_dna",
    "",
    "protein2dna:local",
    29,
    0,
    3,
    100,
    109,
    "   1 : ***UnkSec :    3\n" +
    "       |||||||||\n" +
    "       ***UnkSec\n" +
    " 101 : TAANNNUGA :  109\n",
  )
}

///|
fn exonerate_text_test_multi_report() -> String {
  "Command line: [exonerate -m affine:local:dna2dna queries.fa targets.fa]\n" +
  "Hostname: [aggregate-host]\n\n" +
  exonerate_text_test_alignment(
    "query1", "first query", "target1", "first target", "affine:local:dna2dna", 10,
    0, 3, 100, 103, "   1 : AAA :    3\n       |||\n 101 : AAA :  103\n",
  ) +
  exonerate_text_test_alignment(
    "query1", "first query", "target1", "first target", "affine:local:dna2dna", 20,
    3, 6, 200, 203, "   4 : CCC :    6\n       |||\n 201 : CCC :  203\n",
  ) +
  exonerate_text_test_alignment(
    "query1", "first query", "target2", "second target", "affine:local:dna2dna",
    15, 0, 3, 300, 303, "   1 : GGG :    3\n       |||\n 301 : GGG :  303\n",
  ) +
  exonerate_text_test_alignment(
    "query2", "second query", "target3", "third target", "affine:local:dna2dna",
    30, 0, 3, 400, 403, "   1 : TTT :    3\n       |||\n 401 : TTT :  403\n",
  ) +
  "-- completed exonerate analysis\n"
}

///|
test "Bio.SearchIO.ExonerateIO text parses metadata" {
  let document = exonerate_text_test_parse(exonerate_text_test_dna_report())
  assert_eq(document.metadata.program, "exonerate")
  assert_eq(document.metadata.hostname, "test-host")
  assert_true(document.metadata.command_line.contains("affine:local:dna2dna"))
}

///|
test "Bio.SearchIO.ExonerateIO text aggregates document shape" {
  let document = exonerate_text_test_parse(exonerate_text_test_dna_report())
  assert_eq(document.num_queries(), 1)
  assert_eq(document.num_hits(), 1)
  assert_eq(document.num_hsps(), 1)
}

///|
test "Bio.SearchIO.ExonerateIO text preserves identifiers and descriptions" {
  let query = exonerate_text_test_parse(exonerate_text_test_dna_report()).queries[0]
  assert_eq(query.id, "query1")
  assert_eq(query.description, "demo query")
  assert_eq(query.hits[0].id, "target1")
  assert_eq(query.hits[0].description, "demo target")
}

///|
test "Bio.SearchIO.ExonerateIO text parses DNA alignment rows" {
  let fragment = exonerate_text_test_parse(exonerate_text_test_dna_report()).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_sequence, "ACGT-AC")
  assert_eq(fragment.similarity, "|||| ||")
  assert_eq(fragment.hit_sequence, "ACGTTAC")
}

///|
test "Bio.SearchIO.ExonerateIO text computes DNA coordinates" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_dna_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.query_start, 0)
  assert_eq(hsp.query_end, 6)
  assert_eq(hsp.hit_start, 100)
  assert_eq(hsp.hit_end, 107)
  assert_eq(hsp.fragments[0].query_start, 0)
  assert_eq(hsp.fragments[0].query_end, 6)
}

///|
test "Bio.SearchIO.ExonerateIO text reports strands" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_dna_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.query_strand, @src.ExonerateTextForward)
  assert_eq(hsp.hit_strand, @src.ExonerateTextForward)
}

///|
test "Bio.SearchIO.ExonerateIO text recognizes reverse complement suffix" {
  let query = exonerate_text_test_parse(exonerate_text_test_reverse_report()).queries[0]
  let hsp = query.hits[0].hsps[0]
  assert_eq(query.hits[0].description, "reverse target")
  assert_eq(hsp.hit_strand, @src.ExonerateTextReverse)
  assert_eq(hsp.hit_start, 193)
  assert_eq(hsp.hit_end, 200)
}

///|
test "Bio.SearchIO.ExonerateIO text maps forward positions" {
  let fragment = exonerate_text_test_parse(exonerate_text_test_dna_report()).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_position(0), Some(0))
  assert_eq(fragment.query_position(4), None)
  assert_eq(fragment.query_position(5), Some(4))
  assert_eq(fragment.hit_position(6), Some(106))
}

///|
test "Bio.SearchIO.ExonerateIO text maps reverse positions" {
  let fragment = exonerate_text_test_parse(exonerate_text_test_reverse_report()).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.hit_position(0), Some(199))
  assert_eq(fragment.hit_position(5), Some(194))
}

///|
test "Bio.SearchIO.ExonerateIO text splits target introns" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_intron_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.num_fragments(), 2)
  assert_eq(hsp.fragments[0].query_sequence, "ACG")
  assert_eq(hsp.fragments[1].query_sequence, "TAC")
}

///|
test "Bio.SearchIO.ExonerateIO text computes intron ranges" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_intron_report()).queries[0].hits[0].hsps[0]
  let query_ranges = hsp.query_ranges()
  let hit_ranges = hsp.hit_ranges()
  assert_eq(query_ranges[0], @src.ExonerateTextRange::create(0, 3))
  assert_eq(query_ranges[1], @src.ExonerateTextRange::create(3, 6))
  assert_eq(hit_ranges[0], @src.ExonerateTextRange::create(100, 103))
  assert_eq(hit_ranges[1], @src.ExonerateTextRange::create(110, 113))
}

///|
test "Bio.SearchIO.ExonerateIO text exposes intron inter-ranges" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_intron_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.query_inter_ranges()[0], @src.ExonerateTextRange::create(3, 3))
  assert_eq(
    hsp.hit_inter_ranges()[0],
    @src.ExonerateTextRange::create(103, 110),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text splits NER blocks" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_ner_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.num_fragments(), 2)
  assert_eq(hsp.fragments[0].query_sequence, "AAA")
  assert_eq(hsp.fragments[1].hit_sequence, "CCC")
}

///|
test "Bio.SearchIO.ExonerateIO text applies independent NER lengths" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_ner_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.query_inter_ranges()[0], @src.ExonerateTextRange::create(3, 7))
  assert_eq(
    hsp.hit_inter_ranges()[0],
    @src.ExonerateTextRange::create(103, 109),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text converts protein2dna triplets" {
  let fragment = exonerate_text_test_parse(exonerate_text_test_protein_report()).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_sequence, "MG-K")
  assert_eq(fragment.hit_sequence, "MG-K")
  assert_eq(fragment.query_strand, @src.ExonerateTextProtein)
  assert_eq(fragment.hit_strand, @src.ExonerateTextForward)
}

///|
test "Bio.SearchIO.ExonerateIO text retains translated DNA annotation" {
  let fragment = exonerate_text_test_parse(exonerate_text_test_protein_report()).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.hit_annotation, Some("ATGGGC---AAA"))
  assert_eq(fragment.query_step, 1)
  assert_eq(fragment.hit_step, 3)
}

///|
test "Bio.SearchIO.ExonerateIO text converts dna2protein triplets" {
  let fragment = exonerate_text_test_parse(
      exonerate_text_test_dna_to_protein_report(),
    ).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_sequence, "MGK")
  assert_eq(fragment.hit_sequence, "MGK")
  assert_eq(fragment.query_annotation, Some("ATGGGCAAA"))
  assert_eq(fragment.hit_strand, @src.ExonerateTextProtein)
}

///|
test "Bio.SearchIO.ExonerateIO text splits frameshifts" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_frameshift_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.num_fragments(), 2)
  assert_eq(hsp.fragments[0].query_sequence, "D")
  assert_eq(hsp.fragments[1].query_sequence, "I")
  assert_eq(hsp.fragments[0].hit_sequence, "D")
  assert_eq(hsp.fragments[1].hit_sequence, "I")
}

///|
test "Bio.SearchIO.ExonerateIO text records frameshift coordinate gap" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_frameshift_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.query_inter_ranges()[0], @src.ExonerateTextRange::create(1, 1))
  assert_eq(
    hsp.hit_inter_ranges()[0],
    @src.ExonerateTextRange::create(103, 105),
  )
  assert_eq(hsp.fragments[0].hit_frame, 2)
  assert_eq(hsp.fragments[1].hit_frame, 1)
}

///|
test "Bio.SearchIO.ExonerateIO text computes summary counts" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_dna_report()).queries[0].hits[0].hsps[0]
  let counts = hsp.counts()
  assert_eq(counts.fragments, 1)
  assert_eq(counts.alignment_columns, 7)
  assert_eq(counts.identities, 6)
  assert_eq(counts.gap_columns, 1)
  assert_eq(counts.query_gap_opens, 1)
}

///|
test "Bio.SearchIO.ExonerateIO text stitches wrapped physical blocks" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_wrapped_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.num_fragments(), 1)
  assert_eq(hsp.fragments[0].query_sequence, "ACGTTGCA")
  assert_eq(hsp.fragments[0].hit_sequence, "ACGTTGCA")
  assert_eq(hsp.counts().alignment_columns, 8)
}

///|
test "Bio.SearchIO.ExonerateIO text parses joint intron lengths" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_joint_intron_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.num_fragments(), 2)
  assert_eq(hsp.query_inter_ranges()[0], @src.ExonerateTextRange::create(3, 8))
  assert_eq(
    hsp.hit_inter_ranges()[0],
    @src.ExonerateTextRange::create(103, 110),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text computes reverse intron coordinates" {
  let hsp = exonerate_text_test_parse(
      exonerate_text_test_reverse_intron_report(),
    ).queries[0].hits[0].hsps[0]
  assert_eq(hsp.query_strand, @src.ExonerateTextReverse)
  assert_eq(hsp.hit_strand, @src.ExonerateTextReverse)
  assert_eq(hsp.query_ranges(), [
    @src.ExonerateTextRange::create(7, 10),
    @src.ExonerateTextRange::create(4, 7),
  ])
  assert_eq(
    hsp.hit_inter_ranges()[0],
    @src.ExonerateTextRange::create(190, 197),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text handles split protein codons" {
  let hsp = exonerate_text_test_parse(exonerate_text_test_split_codon_report()).queries[0].hits[0].hsps[0]
  assert_eq(hsp.num_fragments(), 2)
  assert_eq(hsp.fragments[0].query_sequence, "GX")
  assert_eq(hsp.fragments[1].query_sequence, "XA")
  assert_eq(hsp.fragments[0].phase, 0)
  assert_eq(hsp.fragments[1].phase, 1)
  assert_eq(hsp.query_ranges(), [
    @src.ExonerateTextRange::create(0, 1),
    @src.ExonerateTextRange::create(2, 3),
  ])
  assert_eq(hsp.hit_ranges(), [
    @src.ExonerateTextRange::create(100, 105),
    @src.ExonerateTextRange::create(112, 116),
  ])
  assert_eq(hsp.hit_split_codons, [
    @src.ExonerateTextRange::create(103, 105),
    @src.ExonerateTextRange::create(112, 113),
  ])
}

///|
test "Bio.SearchIO.ExonerateIO text flips five-row coding alignments" {
  let fragment = exonerate_text_test_parse(exonerate_text_test_coding_report()).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_sequence, "ATGGAA")
  assert_eq(fragment.hit_sequence, "ATGGAA")
  assert_eq(fragment.query_annotation, Some("MetGlu"))
  assert_eq(fragment.hit_annotation, Some("MetGlu"))
}

///|
test "Bio.SearchIO.ExonerateIO text assigns five-row frameshifts" {
  let hsp = exonerate_text_test_parse(
      exonerate_text_test_coding_frameshift_report(),
    ).queries[0].hits[0].hsps[0]
  assert_eq(hsp.num_fragments(), 2)
  assert_eq(hsp.fragments[0].query_sequence, "ATG")
  assert_eq(hsp.fragments[1].query_sequence, "C-A")
  assert_eq(hsp.query_inter_ranges()[0], @src.ExonerateTextRange::create(3, 4))
  assert_eq(
    hsp.hit_inter_ranges()[0],
    @src.ExonerateTextRange::create(103, 103),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text maps special amino acids" {
  let fragment = exonerate_text_test_parse(
      exonerate_text_test_special_protein_report(),
    ).queries[0].hits[0].hsps[0].fragments[0]
  assert_eq(fragment.query_sequence, "*XU")
  assert_eq(fragment.hit_sequence, "*XU")
}

///|
test "Bio.SearchIO.ExonerateIO text aggregates query hit and HSP levels" {
  let document = exonerate_text_test_parse(exonerate_text_test_multi_report())
  assert_eq(document.num_queries(), 2)
  assert_eq(document.num_hits(), 3)
  assert_eq(document.num_hsps(), 4)
  assert_eq(document.queries[0].hits.length(), 2)
  assert_eq(document.queries[0].hits[0].hsps.length(), 2)
  assert_eq(document.best_hsp().unwrap().score, 30)
}

///|
test "Bio.SearchIO.ExonerateIO text returns defensive hit copies" {
  let document = exonerate_text_test_parse(exonerate_text_test_multi_report())
  let hit = document.queries[0].find_hit("target1").unwrap()
  hit.hsps.clear()
  assert_eq(document.queries[0].hits[0].hsps.length(), 2)
}

///|
test "Bio.SearchIO.ExonerateIO text bundled example parses offline" {
  let document = exonerate_text_test_parse(@src.exonerate_text_example())
  assert_eq(document.num_queries(), 1)
  assert_eq(document.num_hsps(), 1)
  assert_eq(document.queries[0].hits[0].hsps[0].num_fragments(), 2)
}

///|
test "Bio.SearchIO.ExonerateIO text supports no-result reports" {
  let document = exonerate_text_test_parse(
    "Command line: [exonerate query.fa target.fa]\n" +
    "Hostname: [none]\n\n" +
    "-- completed exonerate analysis\n",
  )
  assert_eq(document.num_queries(), 0)
  assert_eq(document.num_hits(), 0)
  assert_eq(document.best_hsp(), None)
}

///|
test "Bio.SearchIO.ExonerateIO text supports CRLF" {
  let text = exonerate_text_test_dna_report().replace_all(old="\n", new="\r\n")
  assert_eq(exonerate_text_test_parse(text).num_hsps(), 1)
}

///|
test "Bio.SearchIO.ExonerateIO text returns defensive query copies" {
  let document = exonerate_text_test_parse(exonerate_text_test_dna_report())
  let copy = document.find_query("query1").unwrap()
  copy.hits.clear()
  assert_eq(document.queries[0].hits.length(), 1)
}

///|
test "Bio.SearchIO.ExonerateIO text returns defensive best-HSP copies" {
  let document = exonerate_text_test_parse(exonerate_text_test_dna_report())
  let copy = document.best_hsp().unwrap()
  copy.fragments.clear()
  assert_eq(document.queries[0].hits[0].hsps[0].fragments.length(), 1)
}

///|
test "Bio.SearchIO.ExonerateIO text summarizes report shape" {
  let summary = exonerate_text_test_parse(exonerate_text_test_dna_report()).summary()
  assert_true(summary.contains("queries=1"))
  assert_true(summary.contains("hits=1"))
  assert_true(summary.contains("hsps=1"))
}

///|
test "Bio.SearchIO.ExonerateIO text rejects missing metadata" {
  assert_true(exonerate_text_test_rejects("-- completed exonerate analysis\n"))
}

///|
test "Bio.SearchIO.ExonerateIO text rejects unfinished reports" {
  assert_true(
    exonerate_text_test_rejects(
      exonerate_text_test_dna_report().replace_all(
        old="-- completed exonerate analysis",
        new="",
      ),
    ),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text rejects incomplete C4 headers" {
  assert_true(
    exonerate_text_test_rejects(
      exonerate_text_test_dna_report().replace_all(
        old="         Model: affine:local:dna2dna\n",
        new="",
      ),
    ),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text rejects malformed ranges" {
  assert_true(
    exonerate_text_test_rejects(
      exonerate_text_test_dna_report().replace_all(old="0 -> 6", new="0 to 6"),
    ),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text rejects inconsistent range spans" {
  assert_true(
    exonerate_text_test_rejects(
      exonerate_text_test_dna_report().replace_all(old="0 -> 6", new="0 -> 7"),
    ),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text rejects reverse orientation mismatch" {
  assert_true(
    exonerate_text_test_rejects(
      exonerate_text_test_reverse_report().replace_all(
        old="200 -> 193",
        new="193 -> 200",
      ),
    ),
  )
}

///|
test "Bio.SearchIO.ExonerateIO text rejects truncated physical rows" {
  assert_true(
    exonerate_text_test_rejects(
      exonerate_text_test_dna_report().replace_all(
        old="       |||| ||\n",
        new="       |||\n",
      ),
    ),
  )
}
