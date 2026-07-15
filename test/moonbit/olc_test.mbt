///|
/// Tests for Overlap-Layout-Consensus algorithm.

///|
test "olc_compute_overlap" {
  let overlap = @src.compute_overlap("ABCDEF", "DEFXYZ", 2)
  assert_eq(overlap, 3)
  
  let overlap2 = @src.compute_overlap("ABC", "DEF", 1)
  assert_eq(overlap2, 0)
  
  let overlap3 = @src.compute_overlap("AAAAA", "AAABBB", 2)
  assert_eq(overlap3, 3)
}

///|
test "olc_compute_all_overlaps" {
  let reads = ["ABCDEF", "DEFXYZ", "XYZ123"]
  let overlaps = @src.compute_all_overlaps(reads, 2)
  assert_eq(overlaps.length(), 2)
}

///|
test "olc_find_layout_path" {
  let reads = ["ABCDEF", "DEFXYZ", "XYZ123"]
  let path = @src.find_layout_path(reads, 2)
  assert_eq(path.length(), 3)
}

///|
test "olc_generate_consensus" {
  let reads = ["ABCDEF", "DEFXYZ", "XYZ123"]
  let path = [0, 1, 2]
  let consensus = @src.generate_consensus(reads, path)
  assert_eq(consensus, "ABCDEFXYZ123")
}

///|
test "olc_assemble_simple" {
  let reads = ["ABCDEF", "DEFXYZ", "XYZ123"]
  let assembly = @src.olc_assemble(reads, 2)
  assert_eq(assembly, "ABCDEFXYZ123")
}

///|
test "olc_assemble_banana" {
  let reads = ["ban", "ana", "nan", "ana"]
  let assembly = @src.olc_assemble(reads, 2)
  assert_true(assembly.length() >= 5)
}

///|
test "olc_assemble_dna" {
  let reads = ["ACGTACGT", "CGTACGTA", "GTACGTAC"]
  let assembly = @src.olc_assemble(reads, 5)
  assert_true(assembly.length() > 0)
}

///|
test "olc_overlap_graph_to_dot" {
  let reads = ["ABCDEF", "DEFXYZ"]
  let dot = @src.overlap_graph_to_dot(reads, 2)
  assert_true(dot.length() > 0)
}

///|
test "olc_single_read" {
  let reads = ["ABCDEF"]
  let assembly = @src.olc_assemble(reads, 2)
  assert_eq(assembly, "ABCDEF")
}

///|
test "olc_empty_reads" {
  let reads : Array[String] = Array::new()
  let assembly = @src.olc_assemble(reads, 2)
  assert_eq(assembly, "")
}

///|
test "olc_correct_errors" {
  let reads = ["ABCDEF", "DEFXYZ"]
  let corrected = @src.correct_errors(reads, 2)
  assert_eq(corrected.length(), 2)
}

///|
test "olc_overlap_struct" {
  let overlap = @src.Overlap::new(0, 1, 3)
  assert_eq(overlap.from(), 0)
  assert_eq(overlap.to(), 1)
  assert_eq(overlap.overlap_len(), 3)
}
