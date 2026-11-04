///|
/// De Bruijn Graph Demo - Sequence Assembly
///
/// This demo shows how to use the De Bruijn graph for sequence assembly.
fn main {
  println("=== De Bruijn Graph Demo ===")
  println("")

  println("1. Basic graph construction:")
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABCDE")
  println(graph.to_string())
  println("")

  println("2. Single sequence assembly:")
  let seq = "ABCDEFGHIJ"
  println("   Original: " + seq)
  match @src.assemble_single(seq, 4) {
    Some(assembled) => println("   Assembled: " + assembled)
    None => println("   Assembly failed")
  }
  println("")

  println("3. Multiple overlapping reads assembly:")
  let reads = ["ABCDE", "CDEFG", "EFGHI"]
  println("   Reads:")
  for i = 0; i < reads.length(); i = i + 1 {
    println("     " + reads[i])
  }
  match @src.assemble_sequences(reads, 3) {
    Some(assembled) => println("   Assembled: " + assembled)
    None => println("   Assembly failed")
  }
  println("")

  println("4. Graph with branches:")
  let graph2 = @src.DeBruijnGraph::new(3)
  graph2.add_sequence("ABCDE")
  graph2.add_sequence("ABXDE")
  println(graph2.to_string())
  println("")

  println("5. Circular sequence (Eulerian circuit):")
  let circular = "ABCA"
  println("   Original: " + circular)
  let graph3 = @src.DeBruijnGraph::new(3)
  graph3.add_sequence(circular)
  println(graph3.to_string())
  match graph3.find_eulerian_path() {
    Some(path) => println("   Eulerian path: " + path)
    None => println("   No Eulerian path found")
  }
}
