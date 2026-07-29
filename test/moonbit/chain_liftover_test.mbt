///|
/// Test chain content: simple single-block chain.
let chain_simple = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"

///|
/// Test chain content: multi-block chain with gaps.
let chain_multi_blocks = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  500 10 20\n  300 5 15\n  200 0 0"

///|
/// Test chain content: reverse query strand.
let chain_reverse_strand = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 - 5000 1 5\n  5000 0 0"

///|
/// Test chain content: two-chain file.
let chain_two_records = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0\nchain 3000 chr1 249250621 + 6000 10000 chr3 8000 + 500 1500 6\n  3000 0 0"

///|
/// Test chain content: three-chain complex file.
let chain_complex = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  400 30 40\n  300 0 0\nchain 3000 chr5 181538259 + 100 200 chr7 159345973 + 300 500 6\n  100 0 0\nchain 2000 chr12 133275309 + 500 800 chr14 107043718 + 700 1000 7\n  200 5 10\n  100 0 0"

///|
/// Test chain content: with blank lines between records.
let chain_with_blank_lines = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n\n  5000 0 0\n\nchain 3000 chr5 181538259 + 100 200 chr7 159345973 + 300 500 6\n\n  100 0 0"

///|
/// Test chain content: tab-separated format.
let chain_tab_separated = "chain\t5000\tchr1\t249250621\t+\t1\t5000\tchr2\t10000\t+\t1\t100\t5\n\t5000\t0\t0"

///|
/// Test chain content: single nucleotide block.
let chain_single_nuc = "chain 100 chr1 1000 + 100 200 chr2 500 + 50 150 1\n  1 0 0"

///|
/// Test chain content: chain with target minus strand.
let chain_target_minus = "chain 5000 chr1 249250621 - 5000 1 chr2 10000 + 1 5000 5\n  5000 0 0"

///|
/// Test chain content: for block boundary tests (1100 total, two blocks of 500+500 with 10 gap).
let chain_block_boundary = "chain 5000 chr1 249250621 + 1 1100 chr2 10000 + 1 1100 5\n  500 10 20\n  500 0 0"

// ============================================================
// cl_parse_chain tests
// ============================================================

///|
test "cl_parse_chain_simple" {
  let cf = @src.cl_parse_chain(chain_simple)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].score, 5000)
  assert_eq(cf.alignments[0].t_name, "chr1")
  assert_eq(cf.alignments[0].t_start, 1)
  assert_eq(cf.alignments[0].t_end, 5000)
  assert_eq(cf.alignments[0].q_name, "chr2")
  assert_eq(cf.alignments[0].q_start, 1)
  assert_eq(cf.alignments[0].q_end, 5000)
  assert_eq(cf.alignments[0].id, 5)
  assert_eq(cf.alignments[0].t_plus, true)
  assert_eq(cf.alignments[0].blocks.length(), 1)
  assert_eq(cf.alignments[0].blocks[0].size, 5000)
  assert_eq(cf.alignments[0].blocks[0].dt, 0)
  assert_eq(cf.alignments[0].blocks[0].dq, 0)
}

///|
test "cl_parse_chain_multiple_blocks" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].blocks.length(), 3)
  assert_eq(cf.alignments[0].blocks[0].size, 500)
  assert_eq(cf.alignments[0].blocks[0].dt, 10)
  assert_eq(cf.alignments[0].blocks[0].dq, 20)
  assert_eq(cf.alignments[0].blocks[1].size, 300)
  assert_eq(cf.alignments[0].blocks[1].dt, 5)
  assert_eq(cf.alignments[0].blocks[1].dq, 15)
  assert_eq(cf.alignments[0].blocks[2].size, 200)
  assert_eq(cf.alignments[0].blocks[2].dt, 0)
  assert_eq(cf.alignments[0].blocks[2].dq, 0)
}

///|
test "cl_parse_chain_empty" {
  let cf = @src.cl_parse_chain("")
  assert_eq(cf.alignments.length(), 0)
}

///|
test "cl_parse_chain_reverse_query_strand" {
  let cf = @src.cl_parse_chain(chain_reverse_strand)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].q_strand, "-")
  assert_eq(cf.alignments[0].q_start, 5000)
  assert_eq(cf.alignments[0].q_end, 1)
}

///|
test "cl_parse_chain_reverse_target_strand" {
  let cf = @src.cl_parse_chain(chain_target_minus)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].t_plus, false)
  assert_eq(cf.alignments[0].t_start, 5000)
  assert_eq(cf.alignments[0].t_end, 1)
}

///|
test "cl_parse_chain_multiple_records" {
  let cf = @src.cl_parse_chain(chain_two_records)
  assert_eq(cf.alignments.length(), 2)
  assert_eq(cf.alignments[0].id, 5)
  assert_eq(cf.alignments[0].q_name, "chr2")
  assert_eq(cf.alignments[1].id, 6)
  assert_eq(cf.alignments[1].q_name, "chr3")
  assert_eq(cf.alignments[0].blocks.length(), 1)
  assert_eq(cf.alignments[1].blocks.length(), 1)
}

///|
test "cl_parse_chain_complex" {
  let cf = @src.cl_parse_chain(chain_complex)
  assert_eq(cf.alignments.length(), 3)
  assert_eq(cf.alignments[0].blocks.length(), 3)
  assert_eq(cf.alignments[1].blocks.length(), 1)
  assert_eq(cf.alignments[2].blocks.length(), 2)
  assert_eq(cf.alignments[0].id, 5)
  assert_eq(cf.alignments[1].id, 6)
  assert_eq(cf.alignments[2].id, 7)
}

///|
test "cl_parse_chain_with_blank_lines" {
  let cf = @src.cl_parse_chain(chain_with_blank_lines)
  assert_eq(cf.alignments.length(), 2)
  assert_eq(cf.alignments[0].blocks.length(), 1)
  assert_eq(cf.alignments[1].blocks.length(), 1)
}

///|
test "cl_parse_chain_tab_separated" {
  let cf = @src.cl_parse_chain(chain_tab_separated)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].blocks.length(), 1)
  assert_eq(cf.alignments[0].blocks[0].size, 5000)
}

///|
test "cl_parse_chain_single_nucleotide_block" {
  let cf = @src.cl_parse_chain(chain_single_nuc)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].blocks.length(), 1)
  assert_eq(cf.alignments[0].blocks[0].size, 1)
  assert_eq(cf.alignments[0].id, 1)
  let result = @src.cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(50))
}

///|
test "cl_parse_chain_id_preserved" {
  let cf = @src.cl_parse_chain("chain 100 chr1 1000 + 100 200 chr2 500 + 50 150 99\n  100 0 0")
  assert_eq(cf.alignments[0].id, 99)
}

// ============================================================
// cl_liftover_pos tests
// ============================================================

///|
test "cl_liftover_pos_forward_basic" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(100))
}

///|
test "cl_liftover_pos_forward_at_start" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_pos(cf, "chr1", 1)
  assert_eq(result, Some(1))
}

///|
test "cl_liftover_pos_forward_at_end" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_pos(cf, "chr1", 4999)
  assert_eq(result, Some(4999))
}

///|
test "cl_liftover_pos_forward_multi_block_first" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(100))
}

///|
test "cl_liftover_pos_forward_multi_block_second" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_liftover_pos(cf, "chr1", 600)
  assert_eq(result, Some(610))
}

///|
test "cl_liftover_pos_forward_multi_block_third" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_liftover_pos(cf, "chr1", 900)
  assert_eq(result, Some(920))
}

///|
test "cl_liftover_pos_reverse_query_strand" {
  let cf = @src.cl_parse_chain(chain_reverse_strand)
  let result = @src.cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(4901))
}

///|
test "cl_liftover_pos_reverse_query_at_start" {
  let cf = @src.cl_parse_chain(chain_reverse_strand)
  let result = @src.cl_liftover_pos(cf, "chr1", 1)
  assert_eq(result, Some(5000))
}

///|
test "cl_liftover_pos_reverse_query_at_end" {
  let cf = @src.cl_parse_chain(chain_reverse_strand)
  let result = @src.cl_liftover_pos(cf, "chr1", 4999)
  assert_eq(result, Some(2))
}

///|
test "cl_liftover_pos_in_gap_returns_none" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_liftover_pos(cf, "chr1", 501)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_out_of_range_before" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_pos(cf, "chr1", 0)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_out_of_range_after" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_pos(cf, "chr1", 6000)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_wrong_chromosome" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_pos(cf, "chrX", 100)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_empty_chains" {
  let cf = @src.cl_parse_chain("")
  let result = @src.cl_liftover_pos(cf, "chr1", 100)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_multi_chain" {
  let cf = @src.cl_parse_chain(chain_two_records)
  let result1 = @src.cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result1, Some(100))
  let result2 = @src.cl_liftover_pos(cf, "chr1", 6500)
  assert_eq(result2, Some(1000))
}

///|
test "cl_liftover_pos_block_boundary" {
  let cf = @src.cl_parse_chain(chain_block_boundary)
  let r1 = @src.cl_liftover_pos(cf, "chr1", 500)
  assert_eq(r1, Some(500))
  let r2 = @src.cl_liftover_pos(cf, "chr1", 501)
  assert_true(r2 is None)
  let r3 = @src.cl_liftover_pos(cf, "chr1", 511)
  assert_eq(r3, Some(521))
}

// ============================================================
// cl_liftover_range tests
// ============================================================

///|
test "cl_liftover_range_basic" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_range(cf, "chr1", 100, 200)
  assert_eq(result, Some((100, 200)))
}

///|
test "cl_liftover_range_across_blocks" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_liftover_range(cf, "chr1", 400, 700)
  assert_eq(result, Some((400, 710)))
}

///|
test "cl_liftover_range_reverse_strand" {
  let cf = @src.cl_parse_chain(chain_reverse_strand)
  let result = @src.cl_liftover_range(cf, "chr1", 100, 200)
  assert_eq(result, Some((4801, 4901)))
}

///|
test "cl_liftover_range_out_of_range" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_range(cf, "chr1", 4900, 6000)
  assert_true(result is None)
}

///|
test "cl_liftover_range_invalid" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_range(cf, "chr1", 200, 100)
  assert_true(result is None)
}

///|
test "cl_liftover_range_spans_multiple_chains" {
  let cf = @src.cl_parse_chain(chain_two_records)
  let result = @src.cl_liftover_range(cf, "chr1", 4900, 7000)
  assert_true(result is None)
}

///|
test "cl_liftover_range_sub_block" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_liftover_range(cf, "chr1", 50, 150)
  assert_eq(result, Some((50, 150)))
}

// ============================================================
// cl_get_chain_summary tests
// ============================================================

///|
test "cl_get_chain_summary_basic" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let summary = @src.cl_get_chain_summary(cf.alignments[0])
  assert_eq(summary["total_block_size"], 1000)
  assert_eq(summary["num_blocks"], 3)
  assert_eq(summary["total_gap_target"], 15)
  assert_eq(summary["total_gap_query"], 35)
  assert_eq(summary["aligned_target"], 4999)
  assert_eq(summary["aligned_query"], 4999)
}

///|
test "cl_get_chain_summary_empty_blocks" {
  let cf = @src.cl_parse_chain("chain 5000 chr1 249250621 + 1 500 chr2 10000 + 1 500 43")
  let summary = @src.cl_get_chain_summary(cf.alignments[0])
  assert_eq(summary["total_block_size"], 0)
  assert_eq(summary["num_blocks"], 0)
  assert_eq(summary["aligned_target"], 499)
  assert_eq(summary["aligned_query"], 499)
}

///|
test "cl_get_chain_summary_single_block" {
  let cf = @src.cl_parse_chain(chain_simple)
  let summary = @src.cl_get_chain_summary(cf.alignments[0])
  assert_eq(summary["total_block_size"], 5000)
  assert_eq(summary["num_blocks"], 1)
  assert_eq(summary["total_gap_target"], 0)
  assert_eq(summary["total_gap_query"], 0)
}

///|
test "cl_get_chain_summary_gap_statistics" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let summary = @src.cl_get_chain_summary(cf.alignments[0])
  assert_eq(summary["total_gap_target"], 15)
  assert_eq(summary["total_gap_query"], 35)
  assert_eq(summary["num_blocks"], 3)
}

// ============================================================
// cl_chain_to_string tests
// ============================================================

///|
test "cl_chain_to_string_roundtrip" {
  let cf = @src.cl_parse_chain(chain_simple)
  let serialized = @src.cl_chain_to_string(cf.alignments[0])
  assert_true(serialized.length() > 0)
  assert_true(serialized.contains("chain 5000"))
  assert_true(serialized.contains("chr1"))
  assert_true(serialized.contains("chr2"))
  assert_true(serialized.contains("5000 0 0"))
}

///|
test "cl_chain_to_string_format" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let s = @src.cl_chain_to_string(cf.alignments[0])
  assert_true(s.contains("chain 5000"))
  assert_true(s.contains("chr1"))
  assert_true(s.contains("500 10 20"))
  assert_true(s.contains("300 5 15"))
  assert_true(s.contains("200 0 0"))
}

///|
test "cl_chain_to_string_preserves_id" {
  let cf = @src.cl_parse_chain("chain 7777 chr1 1000 + 100 200 chr2 500 + 50 150 99\n  100 0 0")
  let s = @src.cl_chain_to_string(cf.alignments[0])
  assert_true(s.contains("chain 7777"))
  assert_true(s.contains("99"))
}

// ============================================================
// cl_find_chain_for_pos tests
// ============================================================

///|
test "cl_find_chain_for_pos_found" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_find_chain_for_pos(cf, "chr1", 100)
  assert_true(result is Some(_))
  match result {
    Some(ch) => assert_eq(ch.id, 5)
    None => ()
  }
}

///|
test "cl_find_chain_for_pos_not_found_wrong_chrom" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_find_chain_for_pos(cf, "chrX", 100)
  assert_true(result is None)
}

///|
test "cl_find_chain_for_pos_not_found_out_of_range" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_find_chain_for_pos(cf, "chr1", 6000)
  assert_true(result is None)
}

///|
test "cl_find_chain_for_pos_empty_chains" {
  let cf = @src.cl_parse_chain("")
  let result = @src.cl_find_chain_for_pos(cf, "chr1", 100)
  assert_true(result is None)
}

///|
test "cl_find_chain_for_pos_multiple_chains" {
  let cf = @src.cl_parse_chain(chain_two_records)
  let r1 = @src.cl_find_chain_for_pos(cf, "chr1", 100)
  match r1 {
    Some(ch) => assert_eq(ch.id, 5)
    None => assert_true(false)
  }
  let r2 = @src.cl_find_chain_for_pos(cf, "chr1", 6500)
  match r2 {
    Some(ch) => assert_eq(ch.id, 6)
    None => assert_true(false)
  }
}

///|
test "cl_find_chain_for_pos_first_block" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_find_chain_for_pos(cf, "chr1", 100)
  assert_true(result is Some(_))
}

///|
test "cl_find_chain_for_pos_second_block" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_find_chain_for_pos(cf, "chr1", 600)
  assert_true(result is Some(_))
}

// ============================================================
// Edge case tests
// ============================================================

///|
test "cl_parse_chain_whitespace_handling" {
  let content = "  chain  5000  chr1  249250621  +  1  5000  chr2  10000  +  1  100  5\n    5000  0  0  "
  let cf = @src.cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].blocks.length(), 1)
}

///|
test "cl_parse_chain_header_only_no_blocks" {
  let content = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5"
  let cf = @src.cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].blocks.length(), 0)
}

///|
test "cl_liftover_pos_at_exact_chain_end" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_pos(cf, "chr1", 5000)
  assert_true(result is None)
}

///|
test "cl_liftover_range_single_point" {
  let cf = @src.cl_parse_chain(chain_simple)
  let result = @src.cl_liftover_range(cf, "chr1", 100, 101)
  assert_eq(result, Some((100, 101)))
}

///|
test "cl_parse_chain_multiple_records_with_boundary" {
  let cf = @src.cl_parse_chain(chain_two_records)
  let result = @src.cl_liftover_pos(cf, "chr1", 5000)
  assert_true(result is None)
  let result2 = @src.cl_liftover_pos(cf, "chr1", 6000)
  assert_eq(result2, Some(500))
}

///|
test "cl_chain_to_string_reverse_strand" {
  let cf = @src.cl_parse_chain(chain_reverse_strand)
  let s = @src.cl_chain_to_string(cf.alignments[0])
  assert_true(s.contains("-"))
  assert_true(s.contains("5000"))
}

///|
test "cl_parse_chain_large_id" {
  let cf = @src.cl_parse_chain("chain 100 chr1 1000 + 100 200 chr2 500 + 50 150 99999\n  100 0 0")
  assert_eq(cf.alignments[0].id, 99999)
}

///|
test "cl_liftover_pos_complex_multi_block" {
  let content = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  1000 50 100\n  2000 30 40\n  2000 0 0"
  let cf = @src.cl_parse_chain(content)
  let r1 = @src.cl_liftover_pos(cf, "chr1", 500)
  assert_eq(r1, Some(500))
  let r2 = @src.cl_liftover_pos(cf, "chr1", 1500)
  assert_eq(r2, Some(1550))
  let r3 = @src.cl_liftover_pos(cf, "chr1", 3500)
  assert_eq(r3, Some(3560))
  let r4 = @src.cl_liftover_pos(cf, "chr1", 4500)
  assert_eq(r4, Some(4560))
}

///|
test "cl_liftover_range_three_blocks" {
  let cf = @src.cl_parse_chain(chain_multi_blocks)
  let result = @src.cl_liftover_range(cf, "chr1", 100, 900)
  assert_eq(result, Some((100, 920)))
}

///|
test "cl_find_chain_for_pos_at_boundary" {
  let cf = @src.cl_parse_chain(chain_two_records)
  let r1 = @src.cl_find_chain_for_pos(cf, "chr1", 4999)
  assert_true(r1 is Some(_))
  let r2 = @src.cl_find_chain_for_pos(cf, "chr1", 5000)
  assert_true(r2 is None)
  let r3 = @src.cl_find_chain_for_pos(cf, "chr1", 6000)
  assert_true(r3 is Some(_))
}