///|
/// treeio demo - Phylogenetic tree I/O parsers and utilities.
/// Inspired by Bioconductor treeio package.

fn main {
  println("=== treeio Demo: Phylogenetic Tree I/O ===")
  println("")

  // Step 1: Newick parsing
  println("--- Step 1: Newick Format Parsing ---")
  let newick = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
  println("Input: " + newick)
  let tree = @src.treeio_parse_newick(newick)
  let (n_total, n_tips, n_internal) = @src.treeio_tree_stats(tree)
  println("Nodes: " + n_total.to_string() + " (tips=" + n_tips.to_string() + ", internal=" + n_internal.to_string() + ")")
  println("Rooted: " + @src.treeio_is_rooted(tree).to_string())
  println("Binary: " + @src.treeio_is_binary(tree).to_string())
  println("Total branch length: " + @src.treeio_total_length(tree).to_string())
  println("Max depth: " + @src.treeio_max_depth(tree).to_string())
  println("Tip labels: " + tree.tip_labels.to_string())
  println("")

  // Step 2: Round-trip write
  println("--- Step 2: Newick Round-trip ---")
  let rewritten = @src.treeio_write_newick(tree)
  println("Rewritten: " + rewritten)
  println("")

  // Step 3: Beast annotated tree
  println("--- Step 3: Beast Annotated Tree ---")
  let beast_newick = "((A:1.0[&height=1.0],B:1.0[&height=1.0]):2.0[&height=3.0],C:3.0[&height=3.0]);"
  println("Input: " + beast_newick)
  let beast_tree = @src.treeio_parse_beast(beast_newick)
  let (bt_total, bt_tips, _) = @src.treeio_tree_stats(beast_tree)
  println("Parsed Beast tree: " + bt_total.to_string() + " nodes, " + bt_tips.to_string() + " tips")
  println("")

  // Step 4: jplace format
  println("--- Step 4: jplace Format ---")

  let jplace = "{\"tree\":\"((A:0.1,B:0.2):0.3,C:0.4);\",\"placements\":[{\"p\":[[1,0.5,0.01,0.02,0.03]],\"n\":[\"seq1\"]}],\"fields\":[\"edge_num\",\"likelihood\",\"like_weight_ratio\",\"pendant_length\",\"distal_length\"]}"
  let jp = @src.treeio_parse_jplace(jplace)
  println("Parsed jplace: " + jp.placements.length().to_string() + " placements, " + jp.n_places.to_string() + " places")
  println("")

  // Step 5: r8s format
  println("--- Step 5: r8s Format ---")
  let r8s = "tree1 = ((A:0.1,B:0.2):0.3,C:0.4);"
  let r8s_trees = @src.treeio_parse_r8s(r8s)
  println("Parsed r8s: " + r8s_trees.length().to_string() + " trees")
  println("")

  // Step 6: Phylip format
  println("--- Step 6: Phylip Format ---")
  let phylip = "4\nA    0.1\nB    0.2\nC    0.3\nD    0.4"
  let phylip_tree = @src.treeio_parse_phylip(phylip)
  let (pt_total, pt_tips, _) = @src.treeio_tree_stats(phylip_tree)
  println("Parsed phylip: " + pt_total.to_string() + " nodes, " + pt_tips.to_string() + " tips")
  println("")

  // Step 7: Example tree
  println("--- Step 7: Example Tree ---")
  let ex = @src.treeio_example_tree()
  let (ex_total, ex_tips, ex_internal) = @src.treeio_tree_stats(ex)
  println("Example tree: " + ex_total.to_string() + " nodes (tips=" + ex_tips.to_string() + ", internal=" + ex_internal.to_string() + ")")
  println("Newick: " + @src.treeio_write_newick(ex))
  println("")

  println("=== treeio Demo Complete ===")
}