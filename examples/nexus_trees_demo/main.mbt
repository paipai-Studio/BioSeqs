///|
/// nexus_trees demo - NEXUS TREES block parsing and tree collection.
/// Inspired by Biopython Bio.Nexus.Trees and Bio.Nexus.Nodes.

fn main {
  println("=== nexus_trees Demo: NEXUS Tree Collection ===")
  println("")

  // Step 1: Parse NEXUS TREES block
  println("--- Step 1: Parse NEXUS TREES Block ---")
  let content = "#NEXUS\nBEGIN TREES;\n  TREE tree1 = ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\n  TREE tree2 = ((E:0.15,F:0.25):0.35,G:0.45);\nEND;\n"
  let collection = @src.nexus_trees_parse(content)
  println("Trees parsed: " + collection.trees.length().to_string())
  let mut i = 0
  while i < collection.trees.length() {
    let tree = collection.trees[i]
    println("  " + tree.name + ": " + tree.node_count().to_string() + " nodes, " + tree.n_taxa.to_string() + " taxa, rooted=" + tree.rooted.to_string())
    i = i + 1
  }
  println("")

  // Step 2: Round-trip Newick
  println("--- Step 2: Round-trip Newick ---")
  let tree = collection.trees[0]
  println("Original: " + tree.to_newick())
  println("")

  // Step 3: Tree statistics
  println("--- Step 3: Tree Statistics ---")
  let leaves = tree.get_leaves()
  println("Leaves: " + leaves.length().to_string())
  let mut j = 0
  while j < leaves.length() {
    println("  " + leaves[j].taxon + " (node " + leaves[j].node_id.to_string() + ")")
    j = j + 1
  }
  println("Total branch length: " + tree.total_length().to_string())
  println("")

  // Step 4: Patristic distances
  println("--- Step 4: Patristic Distances ---")
  let d_ab = tree.patristic_distance("A", "B")
  let d_ac = tree.patristic_distance("A", "C")
  let d_cd = tree.patristic_distance("C", "D")
  println("d(A,B) = " + d_ab.to_string())
  println("d(A,C) = " + d_ac.to_string())
  println("d(C,D) = " + d_cd.to_string())
  println("")

  // Step 5: Rooted flag
  println("--- Step 5: Rooted Flag ---")
  let rooted_collection = @src.nexus_trees_from_newicks([("t1", "[&R]((A:0.1,B:0.2):0.3,C:0.4);")])
  println("Rooted tree: " + rooted_collection.trees[0].rooted.to_string())
  println("")

  // Step 6: TRANSLATE table
  println("--- Step 6: TRANSLATE Table ---")
  let translate_content = "#NEXUS\nBEGIN TREES;\nTRANSLATE\n  1 Homo_sapiens,\n  2 Mus_musculus,\n  3 Rattus_norvegicus\n  ;\n  TREE tree1 = ((1:0.1,2:0.2):0.3,3:0.4);\nEND;\n"
  let translate_collection = @src.nexus_trees_parse(translate_content)
  println("Translated taxa:")
  for k in translate_collection.translate.keys() {
    match translate_collection.translate.get(k) {
      Some(v) => println("  " + k + " -> " + v)
      None => ()
    }
  }
  println("Resolve '1' -> '" + translate_collection.resolve_taxon("1") + "'")
  println("")

  println("=== nexus_trees Demo Complete ===")
}