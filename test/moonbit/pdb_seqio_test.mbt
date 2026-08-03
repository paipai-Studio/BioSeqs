///|
/// Tests for the PDB sequence IO module.

///|
test "pdb_three_to_one_aa_standard" {
  assert_eq(@src.three_to_one_aa("ALA"), "A")
  assert_eq(@src.three_to_one_aa("GLY"), "G")
  assert_eq(@src.three_to_one_aa("VAL"), "V")
  assert_eq(@src.three_to_one_aa("LEU"), "L")
  assert_eq(@src.three_to_one_aa("MET"), "M")
  assert_eq(@src.three_to_one_aa("PHE"), "F")
  assert_eq(@src.three_to_one_aa("TRP"), "W")
  assert_eq(@src.three_to_one_aa("SER"), "S")
  assert_eq(@src.three_to_one_aa("TYR"), "Y")
  assert_eq(@src.three_to_one_aa("CYS"), "C")
}

///|
test "pdb_three_to_one_aa_special" {
  assert_eq(@src.three_to_one_aa("SEC"), "U")
  assert_eq(@src.three_to_one_aa("PYL"), "O")
}

///|
test "pdb_three_to_one_aa_unknown" {
  assert_eq(@src.three_to_one_aa("XXX"), "X")
}

///|
test "pdb_three_to_one_aa_empty" {
  assert_eq(@src.three_to_one_aa(""), "X")
}

///|
test "pdb_three_to_one_aa_lowercase" {
  assert_eq(@src.three_to_one_aa("ala"), "A")
  assert_eq(@src.three_to_one_aa("gly"), "G")
}

///|
test "pdb_parse_seqres_line_standard" {
  let line = "SEQRES   1 A    4  ALA GLY VAL LEU"
  let (chain_id, residues) = @src.parse_seqres_line(line)
  assert_eq(chain_id, "A")
  assert_eq(residues.length(), 4)
  assert_eq(residues[0], "ALA")
  assert_eq(residues[1], "GLY")
  assert_eq(residues[2], "VAL")
  assert_eq(residues[3], "LEU")
}

///|
test "pdb_parse_seqres_line_single_residue" {
  let line = "SEQRES   1 A    1  MET"
  let (chain_id, residues) = @src.parse_seqres_line(line)
  assert_eq(chain_id, "A")
  assert_eq(residues.length(), 1)
  assert_eq(residues[0], "MET")
}

///|
test "pdb_parse_seqres_line_chain_b" {
  let line = "SEQRES   1 B    3  ALA GLY VAL"
  let (chain_id, residues) = @src.parse_seqres_line(line)
  assert_eq(chain_id, "B")
  assert_eq(residues.length(), 3)
}

///|
test "pdb_seqres_parser_basic" {
  let pdb_text = "SEQRES   1 A    4  ALA GLY VAL LEU\n//\n"
  let records = @src.pdb_seqres_parser(pdb_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AGVL")
  assert_eq(records[0].name, "A")
}

///|
test "pdb_seqres_parser_with_pdb_id" {
  let header = "HEADER    PROTEIN" + " ".repeat(33) + "01-JAN-24   1XXX"
  let pdb_text = header + "\nSEQRES   1 A    4  ALA GLY VAL LEU\n//\n"
  let records = @src.pdb_seqres_parser(pdb_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AGVL")
  assert_eq(records[0].id, "1XXX:A")
}

///|
test "pdb_seqres_parser_multiple_chains" {
  let pdb_text = "SEQRES   1 A    2  ALA GLY\nSEQRES   1 B    2  VAL LEU\n//\n"
  let records = @src.pdb_seqres_parser(pdb_text)
  assert_eq(records.length(), 2)
  assert_eq(records[0].seq.to_string(), "AG")
  assert_eq(records[1].seq.to_string(), "VL")
}

///|
test "pdb_seqres_parser_multi_line_chain" {
  let pdb_text = "SEQRES   1 A    5  ALA GLY VAL LEU MET\nSEQRES   2 A    5  PHE TRP SER TYR CYS\n//\n"
  let records = @src.pdb_seqres_parser(pdb_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AGVLMFWSYC")
}

///|
test "pdb_seqres_parser_empty" {
  let records = @src.pdb_seqres_parser("")
  assert_eq(records.length(), 0)
}

///|
test "pdb_extract_sequence_from_pdb_lines_basic" {
  let lines = [
    "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N  ",
    "ATOM      2  CA  ALA A   1      11.000  21.000  31.000  1.00 20.00           C  ",
    "ATOM      3  N   GLY A   2      12.000  22.000  32.000  1.00 20.00           N  ",
    "ATOM      4  CA  GLY A   2      13.000  23.000  33.000  1.00 20.00           C  ",
    "ATOM      5  N   VAL A   3      14.000  24.000  34.000  1.00 20.00           N  ",
    "ATOM      6  CA  VAL A   3      15.000  25.000  35.000  1.00 20.00           C  ",
    "ATOM      7  N   LEU A   4      16.000  26.000  36.000  1.00 20.00           N  ",
    "ATOM      8  CA  LEU A   4      17.000  27.000  37.000  1.00 20.00           C  ",
  ]
  let seq = @src.extract_sequence_from_pdb_lines(lines, "A")
  assert_eq(seq, "AGVL")
}

///|
test "pdb_extract_sequence_from_pdb_lines_wrong_chain" {
  let lines = [
    "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N  ",
    "ATOM      2  N   GLY B   1      12.000  22.000  32.000  1.00 20.00           N  ",
  ]
  let seq_a = @src.extract_sequence_from_pdb_lines(lines, "A")
  assert_eq(seq_a, "A")
  let seq_b = @src.extract_sequence_from_pdb_lines(lines, "B")
  assert_eq(seq_b, "G")
}

///|
test "pdb_extract_sequence_from_pdb_lines_empty" {
  let seq = @src.extract_sequence_from_pdb_lines([], "A")
  assert_eq(seq, "")
}

///|
test "pdb_extract_sequence_from_pdb_lines_skips_water" {
  let lines = [
    "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N  ",
    "ATOM      2  O   HOH A   2      10.000  20.000  30.000  1.00 20.00           O  ",
    "ATOM      3  N   GLY A   3      12.000  22.000  32.000  1.00 20.00           N  ",
  ]
  let seq = @src.extract_sequence_from_pdb_lines(lines, "A")
  assert_eq(seq, "AG")
}

///|
test "pdb_atom_parser_basic" {
  let pdb_text = "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N  \nATOM      2  CA  ALA A   1      11.000  21.000  31.000  1.00 20.00           C  \nATOM      3  N   GLY A   2      12.000  22.000  32.000  1.00 20.00           N  \nATOM      4  CA  GLY A   2      13.000  23.000  33.000  1.00 20.00           C  \nATOM      5  N   VAL A   3      14.000  24.000  34.000  1.00 20.00           N  \nATOM      6  CA  VAL A   3      15.000  25.000  35.000  1.00 20.00           C  \nATOM      7  N   LEU A   4      16.000  26.000  36.000  1.00 20.00           N  \nATOM      8  CA  LEU A   4      17.000  27.000  37.000  1.00 20.00           C  \nTER       9      LEU A   4\nEND\n"
  let records = @src.pdb_atom_parser(pdb_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AGVL")
  assert_eq(records[0].name, "A")
}

///|
test "pdb_atom_parser_with_pdb_id" {
  let header = "HEADER    PROTEIN" + " ".repeat(33) + "01-JAN-24   1XXX"
  let pdb_text = header + "\nATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N  \nATOM      2  N   GLY A   2      12.000  22.000  32.000  1.00 20.00           N  \nTER       3      GLY A   2\nEND\n"
  let records = @src.pdb_atom_parser(pdb_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AG")
  assert_eq(records[0].id, "1XXX:A")
}

///|
test "pdb_atom_parser_empty" {
  let records = @src.pdb_atom_parser("")
  assert_eq(records.length(), 0)
}

///|
test "pdb_to_seqrecords_atom_source" {
  let pdb_text = "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N  \nATOM      2  N   GLY A   2      12.000  22.000  32.000  1.00 20.00           N  \nTER       3      GLY A   2\nEND\n"
  let records = @src.pdb_to_seqrecords(pdb_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AG")
}

///|
test "pdb_to_seqrecords_seqres_source" {
  let pdb_text = "SEQRES   1 A    4  ALA GLY VAL LEU\n//\n"
  let records = @src.pdb_to_seqrecords(pdb_text, source="seqres")
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "AGVL")
}

///|
test "pdb_to_seqrecords_default_source" {
  let pdb_text = "ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00 20.00           N  \nTER       2      ALA A   1\nEND\n"
  let records = @src.pdb_to_seqrecords(pdb_text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "A")
}

///|
test "pdb_write_pdb_seqrecords_round_trip" {
  let pdb_text = "SEQRES   1 A    4  ALA GLY VAL LEU\n//\n"
  let records = @src.pdb_seqres_parser(pdb_text)
  let output = @src.write_pdb_seqrecords(records)
  assert_true(output.contains("SEQRES"))
  assert_true(output.contains("ALA"))
  assert_true(output.contains("GLY"))
  assert_true(output.contains("VAL"))
  assert_true(output.contains("LEU"))
}

///|
test "pdb_write_pdb_seqrecords_basic" {
  let record = @src.SeqRecord::new(
    @src.Seq::new("AGVL"),
    id="A",
    name="A",
  )
  let output = @src.write_pdb_seqrecords([record])
  assert_true(output.contains("SEQRES"))
  assert_true(output.contains("ALA"))
  assert_true(output.contains("GLY"))
  assert_true(output.contains("VAL"))
  assert_true(output.contains("LEU"))
}

///|
test "pdb_write_pdb_seqrecords_multiple_records" {
  let r1 = @src.SeqRecord::new(
    @src.Seq::new("AG"),
    id="A",
    name="A",
  )
  let r2 = @src.SeqRecord::new(
    @src.Seq::new("VL"),
    id="B",
    name="B",
  )
  let output = @src.write_pdb_seqrecords([r1, r2])
  assert_true(output.contains("SEQRES"))
  assert_true(output.contains("ALA"))
  assert_true(output.contains("GLY"))
  assert_true(output.contains("VAL"))
  assert_true(output.contains("LEU"))
}

///|
test "pdb_write_pdb_seqrecords_empty" {
  let output = @src.write_pdb_seqrecords([])
  assert_eq(output, "")
}

///|
test "pdb_write_pdb_seqrecords_unknown_residue" {
  let record = @src.SeqRecord::new(
    @src.Seq::new("AXG"),
    id="A",
    name="A",
  )
  let output = @src.write_pdb_seqrecords([record])
  assert_true(output.contains("SEQRES"))
  assert_true(output.contains("ALA"))
  assert_true(output.contains("UNK"))
  assert_true(output.contains("GLY"))
}
