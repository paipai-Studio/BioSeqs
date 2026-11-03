///|
/// Nexus format parsing example.
///
/// Run: moon run examples/nexus_demo/main.mbt
fn main {
  let nexus_content = @bio.sample_nexus_content()

  println("=== Nexus Parser Demo ===")
  println("")

  let nf = @bio.parse_nexus(nexus_content)
  println("Blocks found: \{nf.blocks.length()}")
  println("Taxa found: \{nf.taxa.length()}")

  let matrix = nf.get_data_matrix()
  println("Data matrix rows: \{matrix.length()}")

  let trees = nf.get_trees()
  println("Trees found: \{trees.length()}")

  let dist_matrix = nf.get_distance_matrix()
  println("Distance matrix entries: \{dist_matrix.length()}")

  println("")
  println("=== Done ===")
}
