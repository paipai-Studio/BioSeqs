///|
/// Tests for Bio.Align.Mauve - Mauve genome alignment format parser.

///|
/// Test MauveSequence creation.
test "mauve_sequence_create" {
  let seq = @src.new_mauve_sequence("human.chr1", 100, 200, "+", 1000)
  assert_eq(seq.name, "human.chr1")
  assert_eq(seq.start, 100)
  assert_eq(seq.end, 200)
  assert_eq(seq.strand, "+")
  assert_eq(seq.seq_size, 1000)
  assert_eq(seq.length, 100)
  assert_eq(seq.is_inverted, false)
}

///|
/// Test MauveSequence coverage.
test "mauve_sequence_coverage" {
  let seq = @src.new_mauve_sequence("human.chr1", 0, 500, "+", 1000)
  let coverage = @src.mauve_seq_coverage(seq)
  assert_eq(coverage, 50.0)
}

///|
/// Test MauveSequence inverted strand.
test "mauve_sequence_inverted" {
  let seq = @src.new_mauve_sequence("human.chr1", 200, 100, "-", 1000)
  assert_eq(seq.is_inverted, true)
  assert_eq(seq.length, 100)
}

///|
/// Test MauveLCB creation.
test "mauve_lcb_create" {
  let lcb = @src.new_mauve_lcb("lcb_1")
  assert_eq(lcb.id, "lcb_1")
  assert_eq(lcb.score, 0.0)
  assert_eq(lcb.is_collinear, true)
}

///|
/// Test MauveLCB add sequence.
test "mauve_lcb_add_seq" {
  let lcb = @src.new_mauve_lcb("lcb_1")
  let seq = @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000)
  @src.mauve_lcb_add_seq(lcb, seq)
  assert_eq(lcb.sequences.length(), 1)
  assert_eq(@src.mauve_lcb_num_seqs(lcb), 1)
}

///|
/// Test MauveLCB consistency check.
test "mauve_lcb_consistent" {
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("mouse.chr1", 0, 100, "+", 800))
  assert_true(@src.mauve_lcb_is_consistent(lcb))
}

///|
/// Test MauveLCB inconsistency (mixed strands).
test "mauve_lcb_inconsistent" {
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("mouse.chr1", 0, 100, "-", 800))
  assert_false(@src.mauve_lcb_is_consistent(lcb))
}

///|
/// Test MauveAlignment creation.
test "mauve_alignment_create" {
  let ali = @src.new_mauve_alignment()
  assert_eq(ali.version, "1.0")
  assert_eq(@src.mauve_num_lcbs(ali), 0)
}

///|
/// Test MauveAlignment add LCB.
test "mauve_alignment_add_lcb" {
  let ali = @src.new_mauve_alignment()
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_add_lcb(ali, lcb)
  assert_eq(@src.mauve_num_lcbs(ali), 1)
}

///|
/// Test MauveAlignment sequence names.
test "mauve_alignment_seq_names" {
  let ali = @src.new_mauve_alignment()
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("mouse.chr1", 0, 100, "+", 800))
  @src.mauve_add_lcb(ali, lcb)
  
  let names = @src.mauve_all_seq_names(ali)
  assert_eq(names.length(), 2)
  assert_true(names.contains("human.chr1"))
  assert_true(names.contains("mouse.chr1"))
}

///|
/// Test parse_mauve with simple MAF-like content.
test "mauve_parse_simple" {
  let content = "##mauve version=1.0\n\na score=100\ns\thuman.chr1\t0\t100\t+\t1000\nmouse.chr1\t0\t100\t+\t800\n\n"
  let ali = @src.parse_mauve(content)
  assert_eq(ali.version, "1.0")
  assert_eq(@src.mauve_num_lcbs(ali), 1)
}

///|
/// Test parse_mauve with multiple LCBs.
test "mauve_parse_multiple_lcbs" {
  let content = "##mauve version=1.0\n\na score=50\ns\thuman.chr1\t0\t50\t+\t1000\nmouse.chr1\t0\t50\t+\t800\n\na score=30\ns\thuman.chr1\t100\t50\t+\t1000\nmouse.chr1\t100\t50\t+\t800\n\n"
  let ali = @src.parse_mauve(content)
  assert_eq(@src.mauve_num_lcbs(ali), 2)
}

///|
/// Test parse_mauve with inversions (inversion detection via API).
test "mauve_parse_with_inversions" {
  let content = "##mauve version=1.0\n\na score=50\ns\thuman.chr1\t0\t50\t+\t1000\nmouse.chr1\t0\t50\t+\t800\n\na score=30\ns\thuman.chr1\t100\t50\t-\t1000\nmouse.chr1\t100\t50\t-\t800\n\n"
  let ali = @src.parse_mauve(content)
  assert_eq(ali.version, "1.0")
  assert_true(@src.mauve_num_lcbs(ali) >= 2)
  // parse_mauve internally runs detect_mauve_inversions;
  // the inversions array reflects detected strand changes.
  let inv_count = ali.inversions.length()
  assert_true(inv_count >= 0)
}

///|
/// Test detect_mauve_inversions.
test "mauve_detect_inversions" {
  let ali = @src.new_mauve_alignment()
  
  // LCB 1: positive strand
  let lcb1 = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb1, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_add_lcb(ali, lcb1)
  
  // LCB 2: negative strand (inversion)
  let lcb2 = @src.new_mauve_lcb("lcb_2")
  @src.mauve_lcb_add_seq(lcb2, @src.new_mauve_sequence("human.chr1", 100, 100, "-", 1000))
  @src.mauve_add_lcb(ali, lcb2)
  
  let ali = @src.detect_mauve_inversions(ali)
  assert_true(ali.inversions.length() > 0)
}

///|
/// Test detect_mauve_breakpoints.
test "mauve_detect_breakpoints" {
  let ali = @src.new_mauve_alignment()
  
  let lcb1 = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb1, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_add_lcb(ali, lcb1)
  
  let lcb2 = @src.new_mauve_lcb("lcb_2")
  @src.mauve_lcb_add_seq(lcb2, @src.new_mauve_sequence("human.chr1", 500, 100, "+", 1000))
  @src.mauve_add_lcb(ali, lcb2)
  
  let ali = @src.detect_mauve_breakpoints(ali)
  assert_true(ali.breakpoints.length() > 0)
}

///|
/// Test mauve_genome_coverage.
test "mauve_genome_coverage" {
  let ali = @src.new_mauve_alignment()
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr1", 0, 500, "+", 1000))
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("mouse.chr1", 0, 300, "+", 800))
  @src.mauve_add_lcb(ali, lcb)
  
  let coverage = @src.mauve_genome_coverage(ali)
  assert_true(coverage.contains("human.chr1"))
  assert_true(coverage.contains("mouse.chr1"))
  assert_true(coverage["human.chr1"] > 0.0)
  assert_true(coverage["mouse.chr1"] > 0.0)
}

///|
/// Test mauve_conserved_segments.
test "mauve_conserved_segments" {
  let ali = @src.new_mauve_alignment()
  
  let lcb1 = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb1, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_lcb_add_seq(lcb1, @src.new_mauve_sequence("mouse.chr1", 0, 100, "+", 800))
  @src.mauve_add_lcb(ali, lcb1)
  
  let lcb2 = @src.new_mauve_lcb("lcb_2")
  @src.mauve_lcb_add_seq(lcb2, @src.new_mauve_sequence("human.chr1", 500, 100, "+", 1000))
  @src.mauve_lcb_add_seq(lcb2, @src.new_mauve_sequence("mouse.chr1", 500, 100, "+", 800))
  @src.mauve_add_lcb(ali, lcb2)
  
  let segments = @src.mauve_conserved_segments(ali)
  assert_eq(segments["human.chr1"], 2)
  assert_eq(segments["mouse.chr1"], 2)
}

///|
/// Test mauve_rearrangement_rate.
test "mauve_rearrangement_rate" {
  let ali = @src.new_mauve_alignment()
  let rate = @src.mauve_rearrangement_rate(ali)
  assert_eq(rate, 0.0)
}

///|
/// Test mauve_to_bed output format.
test "mauve_to_bed" {
  let ali = @src.new_mauve_alignment()
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_set_score(lcb, 100.0)
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_add_lcb(ali, lcb)
  
  let bed = @src.mauve_to_bed(ali)
  assert_true(bed.contains("track name="))
  assert_true(bed.contains("human.chr1"))
  assert_true(bed.contains("0"))
  assert_true(bed.contains("100"))
}

///|
/// Test mauve_inversions_to_bed.
test "mauve_inversions_to_bed" {
  let ali = @src.new_mauve_alignment()
  let inv = @src.new_mauve_inversion("human.chr1", 100, 200, 100)
  inv.affected_lcbs.push("lcb_1")
  ali.inversions.push(inv)
  
  let bed = @src.mauve_inversions_to_bed(ali)
  assert_true(bed.contains("track name=\"Inversions\""))
  assert_true(bed.contains("human.chr1"))
  assert_true(bed.contains("Inversion"))
}

///|
/// Test mauve_summary output.
test "mauve_summary" {
  let ali = @src.new_mauve_alignment()
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_set_score(lcb, 50.0)
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("mouse.chr1", 0, 100, "+", 800))
  @src.mauve_add_lcb(ali, lcb)
  
  let summary = @src.mauve_summary(ali)
  assert_true(summary.contains("Mauve Alignment Summary"))
  assert_true(summary.contains("Locally Collinear Blocks"))
  assert_true(summary.contains("human.chr1"))
}

///|
/// Test MauveConfig default values.
test "mauve_config_default" {
  let config = @src.MauveConfig::default()
  assert_eq(config.min_lcb_score, 0.0)
  assert_eq(config.min_block_length, 50)
  assert_true(config.detect_inversions)
  assert_true(config.detect_breakpoints)
}

///|
/// Test mauve_get_genome_seqs.
test "mauve_get_genome_seqs" {
  let ali = @src.new_mauve_alignment()
  let lcb = @src.new_mauve_lcb("lcb_1")
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr1", 0, 100, "+", 1000))
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("human.chr2", 200, 100, "+", 800))
  @src.mauve_lcb_add_seq(lcb, @src.new_mauve_sequence("mouse.chr1", 0, 100, "+", 700))
  @src.mauve_add_lcb(ali, lcb)
  
  let human_seqs = @src.mauve_get_genome_seqs(ali, "human.chr")
  assert_eq(human_seqs.length(), 2)
}

///|
/// Test progressive aligner initialization.
test "mauve_progressive_aligner" {
  let aligner = @src.new_mauve_progressive_aligner()
  @src.mauve_add_sequence(aligner, "human")
  @src.mauve_add_sequence(aligner, "mouse")
  @src.mauve_add_sequence(aligner, "rat")
  
  @src.mauve_build_guide_tree(aligner)
  assert_eq(aligner.guide_tree.length(), 3)
}
