///|
fn main {
  println("=== Graph Demo ===")

  let ppi_graph = @src.create_example_graph()

  println("\n1. PPI Graph Summary:")
  println(@src.graph_summarize(ppi_graph))

  println("\n2. Get neighbors of TP53:")
  let tp53_neighbors = @src.graph_get_neighbors(ppi_graph, "TP53")
  let mut i = 0
  while i < tp53_neighbors.length() {
    println(tp53_neighbors[i])
    if i < tp53_neighbors.length() - 1 {
      println(", ")
    }
    i = i + 1
  }
  println("")

  println("\n3. Get degree of BRCA1:")
  println("  Degree: " + @src.graph_get_degree(ppi_graph, "BRCA1").to_string())

  println("\n4. Check edge between TP53 and BRCA1:")
  println(
    "  Has edge: " + @src.graph_has_edge(ppi_graph, "TP53", "BRCA1").to_string(),
  )

  println("\n5. Shortest path from TP53 to CHEK2:")
  let path = @src.graph_shortest_path(ppi_graph, "TP53", "CHEK2")
  if path.is_some() {
    let path_arr = path.unwrap()
    i = 0
    while i < path_arr.length() {
      println(path_arr[i])
      if i < path_arr.length() - 1 {
        println(" -> ")
      }
      i = i + 1
    }
    println("")
  }

  println("\n6. Connected components:")
  let components = @src.graph_connected_components(ppi_graph)
  println("  Number of components: " + components.length().to_string())

  println("\n=== Directed Pathway Graph ===")
  let pathway_graph = @src.create_pathway_graph()

  println("\n7. Pathway Graph Summary:")
  println(@src.graph_summarize(pathway_graph))

  println("\n8. Shortest path from G1 to M phase:")
  let cell_cycle_path = @src.graph_shortest_path(pathway_graph, "G1", "M")
  if cell_cycle_path.is_some() {
    let cc_path_arr = cell_cycle_path.unwrap()
    i = 0
    while i < cc_path_arr.length() {
      println(cc_path_arr[i])
      if i < cc_path_arr.length() - 1 {
        println(" -> ")
      }
      i = i + 1
    }
    println("")
  }

  println("\n9. Generate DOT format:")
  let dot = @src.graph_to_dot(ppi_graph)
  println(dot)

  println("\nGraph demo completed!")
}
