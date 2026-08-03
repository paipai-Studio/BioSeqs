///|
/// pheatmap_demo - Example demonstrating enhanced heatmap with clustering
/// Inspired by Bioconductor pheatmap package.

///|
fn main {
  println("=== pheatmap Demo ===")
  println("Enhanced heatmap visualization with clustering")
  println("")

  // Create a sample gene expression matrix
  let mat = [
    [5.2, 3.1, 8.4, 2.5, 7.8],
    [1.2, 6.7, 4.3, 9.1, 3.6],
    [7.9, 2.4, 5.6, 1.8, 8.2],
    [3.5, 8.9, 1.7, 6.2, 4.1],
    [9.3, 4.5, 7.1, 3.8, 5.9],
  ]

  let row_names = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]
  let col_names = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]

  // Create heatmap input
  let input = @src.pheatmap_input(mat, row_names, col_names)
  println("Created heatmap input with \{input.mat.length()} rows x \{input.mat[0].length()} cols")

  // Set title
  let with_title = @src.pheatmap_set_title(input, "Gene Expression Heatmap")
  println("Set title: \{with_title.main_title}")

  // Set font size
  let nice_font = @src.pheatmap_set_fontsize(with_title, 12.0)
  println("Set font size to \{nice_font.font_size}")

  // Compute distance matrix
  let dist = @src.pheatmap_distance_matrix(mat, @src.distance_method_euclidean())
  println("Distance matrix computed (\{dist.length()}x\{dist[0].length()})")

  // Compute column distance
  let col_dist = @src.pheatmap_column_distance_matrix(mat, @src.distance_method_euclidean())
  println("Column distance matrix computed (\{col_dist.length()}x\{col_dist[0].length()})")

  // Render the heatmap
  let result = @src.pheatmap_render(nice_font)
  println("")
  println("Heatmap rendered!")
  println("Row order: \{result.row_order}")
  println("Col order: \{result.col_order}")
  println("")
  println("Heatmap text output:")
  println(result.text)

  println("")
  println("=== pheatmap Demo Complete ===")
}
