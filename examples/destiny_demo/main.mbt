///|
/// destiny Demo - Single-cell Diffusion Map Dimensionality Reduction.
/// Inspired by Bioconductor destiny package.

///|
fn main {
  println("=== destiny Diffusion Map Demo ===")
  println("Inspired by Bioconductor destiny package")
  println("")
  
  println("1. Creating single-cell data...")
  let cells = @src.create_example_sc_data(20, 5)
  println("   Cells: \{cells.length().to_string()}")
  println("   Genes per cell: \{cells[0].expression.length().to_string()}")
  println("")
  
  println("2. Computing distance matrix...")
  let dist_matrix = @src.compute_distance_matrix(cells, "euclidean")
  println("   Distance matrix size: \{dist_matrix.n_cells.to_string()}x\{dist_matrix.n_cells.to_string()}")
  println("")
  
  println("3. Computing diffusion map (2 components)...")
  let result = @src.compute_diffusion_map(cells, 2, 1.0)
  println("   Cells embedded: \{result.cell_ids.length().to_string()}")
  println("   Components: \{result.eigenvalues.length().to_string()}")
  println("")
  
  println("4. Eigenvalues:")
  let mut i = 0
  while i < result.eigenvalues.length() {
    println("   Component \{i.to_string()}: eigenvalue = \{result.eigenvalues[i].to_string()}, explained variance = \{(result.explained_variance[i] * 100.0).to_string()}%")
    i = i + 1
  }
  println("")
  
  println("5. Embedding coordinates (first 5 cells):")
  let coords = @src.destiny_plot_coordinates(result, 0, 1)
  let mut c = 0
  while c < 5 && c < coords.length() {
    let (cell_id, x, y) = coords[c]
    println("   \{cell_id}: DC1 = \{x.to_string()}, DC2 = \{y.to_string()}")
    c = c + 1
  }
  println("")
  
  println("6. Automatic embedding with optimal sigma...")
  let auto_result = @src.destiny_create_embedding(cells, 2)
  println("   Cells embedded: \{auto_result.cell_ids.length().to_string()}")
  println("")
  
  println("7. Diffusion map summary:")
  let summary = @src.destiny_summary(result)
  println("   \{summary}")
  
  println("=== destiny demo completed ===")
}
