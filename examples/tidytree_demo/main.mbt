///|
/// tidytree demo - Tidy interface for phylogenetic tree manipulation.
/// Inspired by Bioconductor tidytree package.

fn main {
  println("=== tidytree Demo: Tidy Tree Interface ===")
  println("")

  // Step 1: Build tree and convert to tibble
  println("--- Step 1: as_tibble ---")
  let tree = @src.treeio_example_tree()
  let tibble = @src.tidy_as_tibble(tree)
  println("Tibble rows: " + tibble.rows.length().to_string())
  println("Nodes: " + tibble.n_nodes.to_string() + ", Tips: " + tibble.n_tips.to_string())
  println("")

  // Step 2: Tree topology queries
  println("--- Step 2: Topology Queries ---")
  let tips = @src.tidy_tips(tree)
  let internal = @src.tidy_internal(tree)
  println("Tips: " + tips.to_string())
  println("Internal nodes: " + internal.to_string())
  println("Is leaf (node 1): " + @src.tidy_is_leaf(tree, 1).to_string())
  println("Node depth (root): " + @src.tidy_node_depth(tree, tree.root).to_string())
  println("")

  // Step 3: Ancestor / offspring
  println("--- Step 3: Ancestor / Offspring ---")
  let first_tip = tips[0]
  let ancestors = @src.tidy_ancestors(tree, first_tip)
  println("Ancestors of " + @src.tidy_get_label(tree, first_tip) + ": " + ancestors.to_string())
  let offspring = @src.tidy_offspring(tree, tree.root)
  println("Offspring of root: " + offspring.length().to_string() + " nodes")
  let children = @src.tidy_child(tree, tree.root)
  println("Children of root: " + children.to_string())
  let parent = @src.tidy_parent(tree, first_tip)
  println("Parent of " + @src.tidy_get_label(tree, first_tip) + ": " + parent.to_string())
  println("")

  // Step 4: MRCA
  println("--- Step 4: Most Recent Common Ancestor ---")
  if tips.length() >= 2 {
    let mrca = @src.tidy_mrca(tree, tips[0], tips[1])
    println("MRCA(" + @src.tidy_get_label(tree, tips[0]) + ", " + @src.tidy_get_label(tree, tips[1]) + ") = node " + mrca.to_string())
  }
  println("")

  // Step 5: Branch lengths
  println("--- Step 5: Branch Lengths ---")
  let bl = @src.tidy_branch_length(tree, first_tip)
  println("Branch length of " + @src.tidy_get_label(tree, first_tip) + ": " + bl.to_string())
  let dr = @src.tidy_dist_to_root(tree, first_tip)
  println("Distance to root: " + dr.to_string())
  println("")

  // Step 6: Patristic distance
  println("--- Step 6: Patristic Distance ---")
  if tips.length() >= 2 {
    let pd = @src.tidy_patristic_distance(tree, tips[0], tips[1])
    println("Patristic distance(" + @src.tidy_get_label(tree, tips[0]) + ", " + @src.tidy_get_label(tree, tips[1]) + ") = " + pd.to_string())
  }
  println("")

  // Step 7: Traversal orders
  println("--- Step 7: Traversal Orders ---")
  let post_order = @src.tidy_post_order(tree)
  let pre_order = @src.tidy_pre_order(tree)
  println("Post-order: " + post_order.to_string())
  println("Pre-order: " + pre_order.to_string())
  println("")

  // Step 8: Filter / select
  println("--- Step 8: Filter and Select ---")
  let tips_only = @src.tidy_filter_tips(tibble)
  let internal_only = @src.tidy_filter_internal(tibble)
  println("Tips-only tibble: " + tips_only.rows.length().to_string() + " rows")
  println("Internal-only tibble: " + internal_only.rows.length().to_string() + " rows")
  let selected = @src.tidy_select(tibble, [tree.root])
  println("Selected root: " + selected.rows.length().to_string() + " row")
  println("")

  // Step 9: Print tibble
  println("--- Step 9: Print Tibble ---")
  println(@src.tidy_print(tibble))
  println("")

  println("=== tidytree Demo Complete ===")
}