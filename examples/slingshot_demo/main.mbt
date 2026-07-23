///|
/// slingshot Single-cell Trajectory Inference Demo.
fn main {
  println("=== slingshot Single-cell Trajectory Inference Demo ===\n")

  // Simulate single-cell data with 3 clusters
  println("1. Creating simulated single-cell data...")
  let cell_coords : Array[Array[Double]] = [
    [0.0, 0.0], [0.2, 0.1], [0.4, 0.2],   // Cluster 1 (root)
    [1.0, 0.5], [1.5, 0.8], [2.0, 1.0],   // Cluster 2 (intermediate)
    [3.0, 1.2], [3.5, 1.3], [4.0, 1.5],   // Cluster 3 (terminal)
    [2.5, 2.0], [3.0, 2.5], [3.5, 3.0]    // Cluster 4 (alternative branch)
  ]
  let cluster_labels = ["c1", "c1", "c1", "c2", "c2", "c2", "c3", "c3", "c3", "c4", "c4", "c4"]
  println("   Cells: \{cell_coords.length()}, Clusters: 4")

  // Run slingshot
  println("\n2. Running slingshot trajectory inference...")
  let params = @src.SlingshotParams::new()
  let result = @src.slingshot_run(cell_coords, cluster_labels, params)
  println("   Nodes: \{result.nodes.length()}")
  println("   Edges: \{result.edges.length()}")
  println("   Curves: \{result.curves.length()}")

  // Compute pseudotime
  println("\n3. Cell pseudotime values:")
  let mut i = 0
  while i < result.cell_pseudotime.length() && i < 6 {
    let pt = result.cell_pseudotime[i].to_string()
    println("   Cell \{i}: pseudotime=\{pt}")
    i = i + 1
  }

  // Sorted pseudotime
  println("\n4. Cells sorted by pseudotime:")
  let sorted = result.sling_sorted_pseudotime()
  let mut j = 0
  while j < sorted.length() && j < 5 {
    let (pt, cell_idx) = sorted[j]
    let pt_str = pt.to_string()
    println("   Rank \{j}: Cell \{cell_idx}, pseudotime=\{pt_str}, cluster=\{cluster_labels[cell_idx]}")
    j = j + 1
  }

  // MST edges
  println("\n5. Minimum spanning tree edges:")
  let mut e = 0
  while e < result.edges.length() {
    let w = result.edges[e].weight.to_string()
    println("   \{result.edges[e].from_id} -> \{result.edges[e].to_id} (weight=\{w})")
    e = e + 1
  }

  println("\n=== slingshot Demo Complete ===")
}
