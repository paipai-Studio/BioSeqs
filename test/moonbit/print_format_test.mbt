///|
/// Tests for the Bio.Restriction.PrintFormat port. All three formats (list,
/// number, map) are checked byte-for-byte against Biopython 1.88 output.

///|
fn pf_sp(n : Int) -> String {
  let b = StringBuilder::new()
  for i = 0; i < n; i = i + 1 {
    b.write_char(' ')
  }
  b.to_string()
}

///|
fn pf_dct() -> Map[String, Array[Int]] {
  let m : Map[String, Array[Int]] = Map([], capacity=4)
  m.set("EcoRI", [50])
  m.set("BamHI", [])
  m.set("ApaI", [12])
  m
}

///|
fn pf_seq() -> String {
  "GGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTC"
}

///|
test "print_format_list" {
  let out = @src.print_format_list(
    dct=pf_dct(),
    title="My analysis:\n",
    s1="No site:\n",
  )
  assert_eq(
    out,
    "My analysis:\nApaI       :  12.\nEcoRI      :  50.\nNo site:\nBamHI     \n",
  )
}

///|
test "print_format_number" {
  let out = @src.print_format_number(
    dct=pf_dct(),
    title="My analysis:\n",
    s1="No site:\n",
  )
  assert_eq(
    out,
    "My analysis:\n\n\nenzymes which cut 1 times :\n\nApaI       :  12.\nEcoRI      :  50.\nNo site:\nBamHI     \n",
  )
}

///|
test "print_format_map" {
  let out = @src.print_format_map(
    dct=pf_dct(),
    sequence=pf_seq(),
    title="",
    s1="",
  )
  let expected =
    "           12 ApaI\n" +
    "           |" + pf_sp(48) + "\n" +
    "           |" + pf_sp(37) + "50 EcoRI\n" +
    "           |" + pf_sp(37) + "|" + pf_sp(10) + "\n" +
    "GGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTC\n" +
    "||||||||||||||||||||||||||||||||||||||||||||||||||||||\n" +
    "CCATGGCCCGGGGGGGAGCTCCAGCTGCCATAGCTATTCGAACTATAGCTTAAG\n" +
    "1" + pf_sp(51) + "54\n" +
    "\n" +
    "\n" +
    "   Enzymes which do not cut the sequence.\n" +
    "\n" +
    "BamHI" + pf_sp(5) + "\n"
  assert_eq(out, expected)
}

///|
test "print_format_empty_results" {
  let m : Map[String, Array[Int]] = Map([], capacity=0)
  assert_eq(@src.print_format_list(dct=m, title="T", s1="S"), "TS")
}