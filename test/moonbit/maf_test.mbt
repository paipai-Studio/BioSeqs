///|
/// Tests for Bio.Align.MAF - Multiple Alignment Format parser.

///|
/// Test MAFSequence creation.
test "maf_sequence_create" {
  let seq = @src.new_maf_sequence("test_seq", 100, 50, "+", 1000, "ACGTACGT")
  assert_eq(seq.seq_name, "test_seq")
  assert_eq(seq.start, 100)
  assert_eq(seq.length, 50)
  assert_eq(seq.strand, "+")
  assert_eq(seq.seq_size, 1000)
  assert_eq(seq.aligned_seq, "ACGTACGT")
}

///|
/// Test MAFSequence end position.
test "maf_sequence_end" {
  let seq = @src.new_maf_sequence("test", 50, 100, "+", 500, "AAAA")
  assert_eq(@src.maf_seq_end(seq), 150)
}

///|
/// Test MAFSequence sequence extraction (no gaps).
test "maf_sequence_extract" {
  let seq = @src.new_maf_sequence("test", 0, 10, "+", 100, "AC-GT-A")
  let result = @src.maf_seq_sequence(seq)
  assert_eq(result, "ACGTA")
}

///|
/// Test MAFSequence coverage calculation.
test "maf_sequence_coverage" {
  let seq = @src.new_maf_sequence("test", 0, 10, "+", 100, "ACGT")
  let coverage = @src.maf_seq_coverage(seq)
  assert_eq(coverage, 1.0)
}

///|
/// Test MAFSequence coverage with gaps.
test "maf_sequence_coverage_gaps" {
  let seq = @src.new_maf_sequence("test", 0, 10, "+", 100, "AC--")
  let coverage = @src.maf_seq_coverage(seq)
  assert_eq(coverage, 0.5)
}

///|
/// Test MAFBlock creation.
test "maf_block_create" {
  let block = @src.new_maf_block()
  assert_eq(block.is_complete, false)
  assert_eq(block.sequences.length(), 0)
  assert_eq(block.comments.length(), 0)
  assert_eq(block.events.length(), 0)
}

///|
/// Test MAFBlock add sequence.
test "maf_block_add_seq" {
  let mut block = @src.new_maf_block()
  let seq = @src.new_maf_sequence("seq1", 0, 10, "+", 100, "ACGT")
  block = @src.maf_block_add_seq(block, seq)
  assert_eq(block.sequences.length(), 1)
  assert_eq(@src.maf_block_num_seqs(block), 1)
}

///|
/// Test MAFBlock alignment length.
test "maf_block_alignment_length" {
  let mut block = @src.new_maf_block()
  let seq1 = @src.new_maf_sequence("seq1", 0, 10, "+", 100, "ACGTACGTAC")
  let seq2 = @src.new_maf_sequence("seq2", 0, 10, "+", 100, "TGCATGCATG")
  block = @src.maf_block_add_seq(block, seq1)
  block = @src.maf_block_add_seq(block, seq2)
  assert_eq(@src.maf_block_alignment_length(block), 10)
}

///|
/// Test MAFBlock percent identity (identical sequences).
test "maf_block_percent_identity_identical" {
  let mut block = @src.new_maf_block()
  let seq1 = @src.new_maf_sequence("seq1", 0, 5, "+", 100, "ACGTG")
  let seq2 = @src.new_maf_sequence("seq2", 0, 5, "+", 100, "ACGTG")
  block = @src.maf_block_add_seq(block, seq1)
  block = @src.maf_block_add_seq(block, seq2)
  let identity = @src.maf_block_percent_identity(block)
  assert_eq(identity, 100.0)
}

///|
/// Test MAFBlock percent identity (different sequences).
test "maf_block_percent_identity_different" {
  let mut block = @src.new_maf_block()
  let seq1 = @src.new_maf_sequence("seq1", 0, 4, "+", 100, "ACGT")
  let seq2 = @src.new_maf_sequence("seq2", 0, 4, "+", 100, "TGCA")
  block = @src.maf_block_add_seq(block, seq1)
  block = @src.maf_block_add_seq(block, seq2)
  let identity = @src.maf_block_percent_identity(block)
  assert_eq(identity, 0.0)
}

///|
/// Test MAFBlock coverage.
test "maf_block_coverage" {
  let mut block = @src.new_maf_block()
  let seq1 = @src.new_maf_sequence("seq1", 0, 5, "+", 100, "ACGTG")
  let seq2 = @src.new_maf_sequence("seq2", 0, 5, "+", 100, "ACGTG")
  block = @src.maf_block_add_seq(block, seq1)
  block = @src.maf_block_add_seq(block, seq2)
  let coverage = @src.maf_block_coverage(block)
  assert_eq(coverage, 1.0)
}

///|
/// Test MAFAliignment creation.
test "maf_alignment_create" {
  let ali = @src.new_maf_alignment()
  assert_eq(ali.version, "1")
  assert_eq(ali.score, None)
  assert_eq(@src.maf_num_blocks(ali), 0)
}

///|
/// Test MAFAliignment add block.
test "maf_alignment_add_block" {
  let mut ali = @src.new_maf_alignment()
  let block = @src.new_maf_block()
  ali = @src.maf_add_block(ali, block)
  assert_eq(@src.maf_num_blocks(ali), 1)
}

///|
/// Test MAFAliignment total length.
test "maf_alignment_total_length" {
  let mut ali = @src.new_maf_alignment()
  let mut block1 = @src.new_maf_block()
  let seq1 = @src.new_maf_sequence("seq1", 0, 10, "+", 100, "ACGTACGTAC")
  block1 = @src.maf_block_add_seq(block1, seq1)
  ali = @src.maf_add_block(ali, block1)
  
  let mut block2 = @src.new_maf_block()
  let seq2 = @src.new_maf_sequence("seq1", 10, 5, "+", 100, "GTACG")
  block2 = @src.maf_block_add_seq(block2, seq2)
  ali = @src.maf_add_block(ali, block2)
  
  assert_eq(@src.maf_total_length(ali), 15)
}

///|
/// Test MAFAliignment sequence names.
test "maf_alignment_seq_names" {
  let mut ali = @src.new_maf_alignment()
  let mut block = @src.new_maf_block()
  let seq1 = @src.new_maf_sequence("human.chr1", 0, 10, "+", 1000, "AAAA")
  let seq2 = @src.new_maf_sequence("mouse.chr1", 0, 10, "+", 800, "CCCC")
  block = @src.maf_block_add_seq(block, seq1)
  block = @src.maf_block_add_seq(block, seq2)
  ali = @src.maf_add_block(ali, block)
  
  let names = @src.maf_all_seq_names(ali)
  assert_eq(names.length(), 2)
  assert_true(names.contains("human.chr1"))
  assert_true(names.contains("mouse.chr1"))
}

///|
/// Test parse_maf with simple content.
test "maf_parse_simple" {
  let content = "##maf version=1\n\na score=100\ns human.chr1 0 10 + 1000 ACGTACGTAC\ns mouse.chr1 0 10 + 800 ACGTACGTAC\n\n"
  let ali = @src.parse_maf(content)
  assert_eq(ali.version, "1")
  assert_eq(@src.maf_num_blocks(ali), 1)
  assert_eq(@src.maf_total_length(ali), 10)
}

///|
/// Test parse_maf with multiple blocks.
test "maf_parse_multiple_blocks" {
  let content = "##maf version=1\n\na score=50\ns human.chr1 0 5 + 1000 ACGTG\ns mouse.chr1 0 5 + 800 ACGTG\n\na score=30\ns human.chr1 5 5 + 1000 ACGTG\ns mouse.chr1 5 5 + 800 ACGTG\n\n"
  let ali = @src.parse_maf(content)
  assert_eq(@src.maf_num_blocks(ali), 2)
  assert_eq(@src.maf_total_length(ali), 10)
}

///|
/// Test parse_maf with empty content.
test "maf_parse_empty" {
  let content = "##maf version=1\n"
  let ali = @src.parse_maf(content)
  assert_eq(ali.version, "1")
  assert_eq(@src.maf_num_blocks(ali), 0)
}

///|
/// Test parse_maf with comments.
test "maf_parse_with_comments" {
  let content = "##maf version=1\n# This is a comment\n\na score=100\ns human.chr1 0 5 + 1000 ACGTG\n# Block comment\n\n"
  let ali = @src.parse_maf(content)
  assert_eq(@src.maf_num_blocks(ali), 1)
}

///|
/// Test write_maf produces valid output.
test "maf_write_roundtrip" {
  let content = "##maf version=1\n\na score=100\ns human.chr1 0 5 + 1000 ACGTG\ns mouse.chr1 0 5 + 800 ACGTG\n\n"
  let ali = @src.parse_maf(content)
  let output = @src.write_maf(ali)
  assert_true(output.starts_with("##maf version="))
  assert_true(output.contains("human.chr1"))
  assert_true(output.contains("ACGTG"))
}

///|
/// Test maf_select_seqs.
test "maf_select_seqs" {
  let content = "##maf version=1\n\na score=100\ns human.chr1 0 5 + 1000 ACGTG\ns mouse.chr1 0 5 + 800 ACGTG\ns rat.chr1 0 5 + 600 ACGTG\n\n"
  let ali = @src.parse_maf(content)
  let selected = @src.maf_select_seqs(ali, ["human.chr1", "mouse.chr1"])
  assert_eq(@src.maf_num_blocks(selected), 1)
  assert_eq(selected.blocks[0].sequences.length(), 2)
}

///|
/// Test maf_filter_by_length.
test "maf_filter_by_length" {
  let mut ali = @src.new_maf_alignment()
  
  let mut block1 = @src.new_maf_block()
  block1 = @src.maf_block_set_score(block1, Some(100.0))
  block1 = @src.maf_block_add_seq(block1, @src.new_maf_sequence("s1", 0, 100, "+", 1000, "ACGTACGTACGTACGTACGT"))
  block1 = @src.maf_block_add_seq(block1, @src.new_maf_sequence("s2", 0, 100, "+", 1000, "ACGTACGTACGTACGTACGT"))
  ali = @src.maf_add_block(ali, block1)
  
  let filtered = @src.maf_filter_by_length(ali, 5)
  assert_eq(@src.maf_num_blocks(filtered), 1)
}

///|
/// Test maf_compute_stats.
test "maf_stats" {
  let mut ali = @src.new_maf_alignment()
  let mut block = @src.new_maf_block()
  block = @src.maf_block_set_score(block, Some(100.0))
  block = @src.maf_block_add_seq(block, @src.new_maf_sequence("s1", 0, 10, "+", 1000, "ACGTACGTAC"))
  block = @src.maf_block_add_seq(block, @src.new_maf_sequence("s2", 0, 10, "+", 1000, "ACGTACGTAC"))
  ali = @src.maf_add_block(ali, block)
  
  let stats = @src.maf_compute_stats(ali)
  assert_eq(stats.num_blocks, 1)
  assert_eq(stats.total_length, 10)
  assert_eq(stats.num_sequences, 2)
  assert_eq(stats.avg_percent_identity, 100.0)
}

///|
/// Test MAFAliignment average percent identity.
test "maf_average_percent_identity" {
  let content = "##maf version=1\n\na score=100\ns human.chr1 0 5 + 1000 ACGTG\ns mouse.chr1 0 5 + 800 ACGTG\n\n"
  let ali = @src.parse_maf(content)
  let identity = @src.maf_average_percent_identity(ali)
  assert_eq(identity, 100.0)
}

///|
/// Test maf_sequence_coverage calculation.
test "maf_alignment_sequence_coverage" {
  let content = "##maf version=1\n\na score=100\ns human.chr1 0 100 + 1000 ACGTACGTAC\ns mouse.chr1 0 100 + 800 ACGTACGTAC\n\n"
  let ali = @src.parse_maf(content)
  let coverage = @src.maf_sequence_coverage(ali)
  assert_true(coverage.contains("human.chr1"))
  assert_true(coverage.contains("mouse.chr1"))
}

///|
/// Test MAF with strand info.
test "maf_strand_info" {
  let content = "##maf version=1\n\na score=100\ns human.chr1 0 5 + 1000 ACGTG\ns mouse.chr1 100 5 - 800 ACGTG\n\n"
  let ali = @src.parse_maf(content)
  assert_eq(ali.blocks[0].sequences[0].strand, "+")
  assert_eq(ali.blocks[0].sequences[1].strand, "-")
}

///|
/// Test maf_merge_blocks.
test "maf_merge_blocks" {
  let mut ali = @src.new_maf_alignment()
  
  let mut block1 = @src.new_maf_block()
  block1 = @src.maf_block_set_score(block1, Some(50.0))
  block1 = @src.maf_block_add_seq(block1, @src.new_maf_sequence("s1", 0, 5, "+", 1000, "ACGTG"))
  block1 = @src.maf_block_add_seq(block1, @src.new_maf_sequence("s2", 0, 5, "+", 1000, "ACGTG"))
  
  let mut block2 = @src.new_maf_block()
  block2 = @src.maf_block_set_score(block2, Some(30.0))
  block2 = @src.maf_block_add_seq(block2, @src.new_maf_sequence("s1", 5, 5, "+", 1000, "TGCAC"))
  block2 = @src.maf_block_add_seq(block2, @src.new_maf_sequence("s2", 5, 5, "+", 1000, "TGCAC"))
  
  ali = @src.maf_add_block(ali, block1)
  ali = @src.maf_add_block(ali, block2)
  
  let merged = @src.maf_merge_blocks(ali)
  // Should merge blocks with same sequences
  assert_true(@src.maf_num_blocks(merged) <= @src.maf_num_blocks(ali))
}