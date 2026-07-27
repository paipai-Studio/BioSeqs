///|
/// Tests for Bio.PDB.FragmentMapper module.

///|
test "fm_classify_ss_helix" {
  let result = @src.fm_classify_secondary_structure("HGI")
  assert_eq(result.length(), 3)
  assert_eq(result[0], @src.ft_helix())
  assert_eq(result[1], @src.ft_helix())
  assert_eq(result[2], @src.ft_helix())
}

///|
test "fm_classify_ss_sheet" {
  let result = @src.fm_classify_secondary_structure("EE")
  assert_eq(result.length(), 2)
  assert_eq(result[0], @src.ft_sheet())
  assert_eq(result[1], @src.ft_sheet())
}

///|
test "fm_classify_ss_loop_coil" {
  let result = @src.fm_classify_secondary_structure("LC ")
  assert_eq(result.length(), 3)
  assert_eq(result[0], @src.ft_loop())
  assert_eq(result[1], @src.ft_coil())
  assert_eq(result[2], @src.ft_coil())
}

///|
test "fm_classify_ss_bridge" {
  let result = @src.fm_classify_secondary_structure("b")
  assert_eq(result.length(), 1)
  assert_eq(result[0], @src.ft_bridge())
}

///|
test "fm_classify_ss_turn" {
  let result = @src.fm_classify_secondary_structure("TS")
  assert_eq(result.length(), 2)
  assert_eq(result[0], @src.ft_turn())
  assert_eq(result[1], @src.ft_turn())
}

///|
test "fm_classify_ss_mixed" {
  let result = @src.fm_classify_secondary_structure("HLT")
  assert_eq(result.length(), 3)
  assert_eq(result[0], @src.ft_helix())
  assert_eq(result[1], @src.ft_loop())
  assert_eq(result[2], @src.ft_turn())
}

///|
test "fm_classify_ss_unknown_char" {
  let result = @src.fm_classify_secondary_structure("X")
  assert_eq(result.length(), 1)
  assert_eq(result[0], @src.ft_coil())
}

///|
test "fm_classify_ss_empty" {
  let result = @src.fm_classify_secondary_structure("")
  assert_eq(result.length(), 0)
}

///|
test "fm_assign_fragments_basic" {
  let result = @src.fm_assign_fragments("ALPHA", "HHHHH")
  assert_eq(result.sequence_length, 5)
  assert_eq(result.num_fragments(), 1)
  assert_eq(result.fragments[0].fragment_type, @src.ft_helix())
  assert_eq(result.fragments[0].start_residue, 1)
  assert_eq(result.fragments[0].end_residue, 5)
}

///|
test "fm_assign_fragments_multiple" {
  let result = @src.fm_assign_fragments("ALPHABET", "HHHHLLLL")
  assert_eq(result.sequence_length, 8)
  assert_eq(result.num_fragments(), 2)
  assert_eq(result.fragments[0].fragment_type, @src.ft_helix())
  assert_eq(result.fragments[0].start_residue, 1)
  assert_eq(result.fragments[0].end_residue, 4)
  assert_eq(result.fragments[1].fragment_type, @src.ft_loop())
  assert_eq(result.fragments[1].start_residue, 5)
  assert_eq(result.fragments[1].end_residue, 8)
}

///|
test "fm_assign_fragments_sequence" {
  let result = @src.fm_assign_fragments("ALPHA", "HHHHH")
  let seq = @src.fm_fragment_sequence(result.fragments[0])
  assert_eq(seq, "ALPHA")
}

///|
test "fm_assign_fragments_id" {
  let result = @src.fm_assign_fragments("A", "H")
  assert_eq(result.fragments[0].id, "fragment_0")
}

///|
test "fm_get_fragment_found" {
  let result = @src.fm_assign_fragments("ALPHA", "HHHHH")
  let frag = @src.fm_get_fragment(result, "fragment_0")
  assert_true(frag is Some(_))
  assert_eq(frag.unwrap().start_residue, 1)
}

///|
test "fm_get_fragment_not_found" {
  let result = @src.fm_assign_fragments("ALPHA", "HHHHH")
  let frag = @src.fm_get_fragment(result, "nonexistent")
  assert_true(frag is None)
}

///|
test "fm_fragment_type" {
  let frag = @src.Fragment::new("test", 1, 10, fragment_type=@src.ft_helix())
  assert_eq(@src.fm_fragment_type(frag), @src.ft_helix())
}

///|
test "fm_fragment_boundaries" {
  let frag = @src.Fragment::new("test", 5, 20)
  let (start, end) = @src.fm_fragment_boundaries(frag)
  assert_eq(start, 5)
  assert_eq(end, 20)
}

///|
test "fm_fragment_length" {
  let frag = @src.Fragment::new("test", 3, 10)
  assert_eq(frag.length(), 8)
}

///|
test "fm_coverage_basic" {
  let result = @src.fm_assign_fragments("AC", "HL")
  let cov = @src.fm_coverage(result)
  assert_eq(cov.get(0).unwrap_or(0.0), 1.0)
  assert_eq(cov.get(1).unwrap_or(0.0), 1.0)
}

///|
test "fm_merge_fragments_adjacent_same_type" {
  let frag1 = @src.Fragment::new("f1", 1, 5,
    fragment_type=@src.ft_helix(),
    residues=[@src.FragmentResidue::new("A", 1), @src.FragmentResidue::new("L", 2)])
  let frag2 = @src.Fragment::new("f2", 6, 10,
    fragment_type=@src.ft_helix(),
    residues=[@src.FragmentResidue::new("P", 6), @src.FragmentResidue::new("H", 7)])
  let merged = @src.fm_merge_fragments([frag1, frag2])
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].start_residue, 1)
  assert_eq(merged[0].end_residue, 10)
  assert_eq(merged[0].residues.length(), 4)
}

///|
test "fm_merge_fragments_different_types" {
  let frag1 = @src.Fragment::new("f1", 1, 5, fragment_type=@src.ft_helix())
  let frag2 = @src.Fragment::new("f2", 6, 10, fragment_type=@src.ft_loop())
  let merged = @src.fm_merge_fragments([frag1, frag2])
  assert_eq(merged.length(), 2)
}

///|
test "fm_merge_fragments_not_adjacent" {
  let frag1 = @src.Fragment::new("f1", 1, 5, fragment_type=@src.ft_helix())
  let frag2 = @src.Fragment::new("f2", 10, 15, fragment_type=@src.ft_helix())
  let merged = @src.fm_merge_fragments([frag1, frag2])
  assert_eq(merged.length(), 2)
}

///|
test "fm_merge_fragments_empty" {
  let merged = @src.fm_merge_fragments([])
  assert_eq(merged.length(), 0)
}

///|
test "fm_filter_by_length_basic" {
  let frag1 = @src.Fragment::new("f1", 1, 5)
  let frag2 = @src.Fragment::new("f2", 6, 8)
  let filtered = @src.fm_filter_by_length([frag1, frag2], 4)
  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].id, "f1")
}

///|
test "fm_filter_by_length_all_pass" {
  let frag1 = @src.Fragment::new("f1", 1, 10)
  let frag2 = @src.Fragment::new("f2", 11, 20)
  let filtered = @src.fm_filter_by_length([frag1, frag2], 5)
  assert_eq(filtered.length(), 2)
}

///|
test "fm_filter_by_length_none_pass" {
  let frag1 = @src.Fragment::new("f1", 1, 2)
  let frag2 = @src.Fragment::new("f2", 3, 4)
  let filtered = @src.fm_filter_by_length([frag1, frag2], 10)
  assert_eq(filtered.length(), 0)
}

///|
test "fm_filter_by_length_empty" {
  let filtered = @src.fm_filter_by_length([], 5)
  assert_eq(filtered.length(), 0)
}

///|
test "fm_parse_fragments_basic" {
  let content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      3  C   ALA A   1      3.000   0.000   0.000  1.00  0.00           C\n" +
    "END\n"
  let result = @src.fm_parse_fragments(content)
  assert_eq(result.sequence_length, 1)
  assert_eq(result.num_fragments(), 1)
  assert_eq(result.protein_sequence, "A")
}

///|
test "fm_parse_fragments_multiple_residues" {
  let content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      3  N   SER A   2      3.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      4  CA  SER A   2      4.500   0.000   0.000  1.00  0.00           C\n" +
    "END\n"
  let result = @src.fm_parse_fragments(content)
  assert_eq(result.sequence_length, 2)
  assert_eq(result.protein_sequence, "AS")
  assert_eq(result.num_fragments(), 1)
}

///|
test "fm_parse_fragments_multiple_chains" {
  let content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      3  N   SER B   1      5.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      4  CA  SER B   1      6.500   0.000   0.000  1.00  0.00           C\n" +
    "END\n"
  let result = @src.fm_parse_fragments(content)
  assert_eq(result.num_fragments(), 2)
}

///|
test "fm_parse_fragments_atom_counts" {
  let content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      3  C   ALA A   1      3.000   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      4  O   ALA A   1      4.500   0.000   0.000  1.00  0.00           O\n" +
    "END\n"
  let result = @src.fm_parse_fragments(content)
  let frag = result.fragments[0]
  assert_eq(frag.residues[0].atom_count, 4)
  assert_true(frag.residues[0].is_c_alpha)
  assert_true(frag.residues[0].is_sidechain)
}

///|
test "fm_parse_fragments_empty" {
  let result = @src.fm_parse_fragments("")
  assert_eq(result.sequence_length, 0)
  assert_eq(result.num_fragments(), 0)
}

///|
test "fm_fragment_type_to_string" {
  assert_eq(@src.ft_helix().to_string(), "Helix")
  assert_eq(@src.ft_sheet().to_string(), "Sheet")
  assert_eq(@src.ft_loop().to_string(), "Loop")
  assert_eq(@src.ft_coil().to_string(), "Coil")
  assert_eq(@src.ft_bridge().to_string(), "Bridge")
  assert_eq(@src.ft_turn().to_string(), "Turn")
}

///|
test "fm_fragment_type_to_char" {
  assert_eq(@src.ft_helix().to_char(), "H")
  assert_eq(@src.ft_sheet().to_char(), "E")
  assert_eq(@src.ft_loop().to_char(), "L")
  assert_eq(@src.ft_coil().to_char(), "C")
  assert_eq(@src.ft_bridge().to_char(), "b")
  assert_eq(@src.ft_turn().to_char(), "T")
}

///|
test "fm_assign_fragments_mixed_types" {
  let result = @src.fm_assign_fragments("ABCDEF", "HELTCB")
  assert_eq(result.num_fragments(), 6)
  assert_eq(result.fragments[0].fragment_type, @src.ft_helix())
  assert_eq(result.fragments[1].fragment_type, @src.ft_sheet())
  assert_eq(result.fragments[2].fragment_type, @src.ft_loop())
  assert_eq(result.fragments[3].fragment_type, @src.ft_turn())
  assert_eq(result.fragments[4].fragment_type, @src.ft_coil())
  assert_eq(result.fragments[5].fragment_type, @src.ft_bridge())
}

///|
test "fm_residue_creation" {
  let res = @src.FragmentResidue::new("ALA", 1, atom_count=4, is_c_alpha=true)
  assert_eq(res.res_name, "ALA")
  assert_eq(res.res_num, 1)
  assert_eq(res.atom_count, 4)
  assert_true(res.is_c_alpha)
  assert_false(res.is_sidechain)
}

///|
test "fm_fragment_creation_defaults" {
  let frag = @src.Fragment::new("test", 1, 10)
  assert_eq(frag.id, "test")
  assert_eq(frag.start_residue, 1)
  assert_eq(frag.end_residue, 10)
  assert_eq(frag.chain_id, "A")
  assert_eq(frag.fragment_type, @src.ft_coil())
  assert_eq(frag.residues.length(), 0)
}

///|
test "fm_result_num_fragments" {
  let result = @src.fm_assign_fragments("ALPHA", "HHHLL")
  assert_eq(result.num_fragments(), 2)
}

///|
test "fm_parse_fragments_hetatm" {
  let content = "HETATM    1  O   HOH A   1      0.000   0.000   0.000  1.00  0.00           O\n" +
    "END\n"
  let result = @src.fm_parse_fragments(content)
  assert_eq(result.sequence_length, 1)
}

///|
test "fm_assign_fragments_coverage" {
  let result = @src.fm_assign_fragments("ACGT", "HLTC")
  let cov = @src.fm_coverage(result)
  assert_eq(cov.get(0).unwrap_or(0.0), 1.0)
  assert_eq(cov.get(1).unwrap_or(0.0), 1.0)
  assert_eq(cov.get(2).unwrap_or(0.0), 1.0)
  assert_eq(cov.get(3).unwrap_or(0.0), 1.0)
}

///|
test "fm_merge_fragments_three_adjacent" {
  let frag1 = @src.Fragment::new("f1", 1, 3, fragment_type=@src.ft_helix())
  let frag2 = @src.Fragment::new("f2", 4, 6, fragment_type=@src.ft_helix())
  let frag3 = @src.Fragment::new("f3", 7, 9, fragment_type=@src.ft_helix())
  let merged = @src.fm_merge_fragments([frag1, frag2, frag3])
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].start_residue, 1)
  assert_eq(merged[0].end_residue, 9)
}

///|
test "fm_parse_fragments_sequence_mapping" {
  let content = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      2  CA  ALA A   1      1.500   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      3  N   GLY A   2      3.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      4  CA  GLY A   2      4.500   0.000   0.000  1.00  0.00           C\n" +
    "ATOM      5  N   PHE A   3      6.000   0.000   0.000  1.00  0.00           N\n" +
    "ATOM      6  CA  PHE A   3      7.500   0.000   0.000  1.00  0.00           C\n" +
    "END\n"
  let result = @src.fm_parse_fragments(content)
  assert_eq(result.protein_sequence, "AGF")
}