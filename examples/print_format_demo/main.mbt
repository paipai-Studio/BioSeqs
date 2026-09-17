///|
/// Bio.Restriction.PrintFormat Demo.
///
/// Demonstrates the PrintFormat port: rendering a restriction search result
/// (a map of enzyme name -> 1-based cutting positions, as returned by
/// `RestrictionBatch::search`) in the three supported layouts.

///|
fn main raise {
  let seq = "GGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTC"
  // This map has the same shape as RestrictionBatch::search(seq).
  let dct : Map[String, Array[Int]] = Map([], capacity=3)
  dct.set("EcoRI", [50])
  dct.set("BamHI", [])
  dct.set("ApaI", [12])

  println("=== list (alphabetical) ===")
  println(@src.print_format_list(dct=dct, title="My analysis:\n", s1="No site:\n"))

  println("=== number (grouped by number of cuts) ===")
  println(@src.print_format_number(dct=dct, title="My analysis:\n", s1="No site:\n"))

  println("=== map (graphical) ===")
  println(@src.print_format_map(dct=dct, sequence=seq, title="", s1=""))
}