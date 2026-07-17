///|
/// Overlap-Layout-Consensus Demo - Sequence Assembly
///
/// This demo shows how to use the OLC algorithm for sequence assembly,
/// including overlap detection, layout path finding, and consensus generation.
fn main {
  println("=== Overlap-Layout-Consensus Demo ===")
  println("")

  println("1. Basic overlap detection:")
  let overlap = @src.compute_overlap("ABCDEF", "DEFXYZ", 2)
  println("   Overlap between \"ABCDEF\" and \"DEFXYZ\": \{overlap}")
  println("")

  println("2. Simple sequence assembly:")
  let simple_reads = ["ABCDEF", "DEFXYZ", "XYZ123"]
  println("   Reads:")
  for i = 0; i < simple_reads.length(); i = i + 1 {
    println("     \{i}: \{simple_reads[i]}")
  }
  let simple_assembly = @src.olc_assemble(simple_reads, 2)
  println("   Assembled: \{simple_assembly}")
  println("")

  println("3. DNA sequence assembly:")
  let dna_reads = ["ACGTACGTCCG", "CGTCCGATGCA", "ATGCATGCTGA", "TGCTGACGTTT"]
  println("   Reads:")
  for i = 0; i < dna_reads.length(); i = i + 1 {
    println("     \{i}: \{dna_reads[i]}")
  }
  let dna_assembly = @src.olc_assemble(dna_reads, 5)
  println("   Assembled: \{dna_assembly}")
  println("")

  println("4. Banana example:")
  let banana_reads = ["ban", "ana", "nan", "ana"]
  println("   Reads:")
  for i = 0; i < banana_reads.length(); i = i + 1 {
    println("     \{i}: \{banana_reads[i]}")
  }
  let banana_assembly = @src.olc_assemble(banana_reads, 2)
  println("   Assembled: \{banana_assembly}")
  println("")

  println("5. Find all overlaps:")
  let overlaps = @src.compute_all_overlaps(dna_reads, 5)
  println("   Number of overlaps found: \{overlaps.length()}")
  for i = 0; i < overlaps.length(); i = i + 1 {
    let o = overlaps[i]
    println(
      "     Read \{o.from()} -> Read \{o.to()}: overlap length = \{o.overlap_len()}",
    )
  }
  println("")

  println("6. Overlap graph visualization (Graphviz DOT):")
  let dot = @src.overlap_graph_to_dot(simple_reads, 2)
  println(dot)

  println("7. Layout path finding:")
  let path = @src.find_layout_path(simple_reads, 2)
  println("   Layout path:")
  let path_str = StringBuilder::new()
  for i = 0; i < path.length(); i = i + 1 {
    path_str.write_string(path[i].to_string())
    if i < path.length() - 1 {
      path_str.write_string(" -> ")
    }
  }
  println("   \{path_str.to_string()}")
  println("")

  println("8. Consensus generation from path:")
  let consensus = @src.generate_consensus(simple_reads, path)
  println("   Consensus: \{consensus}")
  println("")

  println("=== Demo Complete ===")
}
