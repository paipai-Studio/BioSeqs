///|
/// factoextra_demo - Example demonstrating PCA/factor analysis
/// Inspired by Bioconductor factoextra package.

///|
fn main {
  println("=== factoextra Demo ===")
  println("PCA/factor analysis utilities")
  println("")

  // Create a sample gene expression data matrix
  // 4 samples x 6 genes
  let expression_data = [
    [2.5, 3.1, 4.2, 1.8, 5.3, 2.9],
    [3.8, 2.4, 5.1, 2.7, 4.6, 3.2],
    [1.9, 4.3, 3.5, 2.1, 6.0, 2.1],
    [4.1, 3.7, 4.8, 3.2, 5.8, 3.6],
  ]

  // Perform PCA
  println("Performing PCA on expression data...")
  let pca_result = @src.facto_pca(expression_data)
  println("")

  // Get eigenvalues
  let eigenvalues = @src.facto_get_eigenvalue(pca_result)
  println("Eigenvalues:")
  for ev in eigenvalues {
    println("  Dim \{ev.dimension}: eigenvalue=\{ev.eigenvalue}, variance=\{ev.variance_explained}%, cumulative=\{ev.cumulative_variance}%")
  }

  // Get PCA summary
  let summary = @src.facto_summary(pca_result)
  println("")
  println("PCA Summary:")
  println("  Individuals: \{summary.n_individuals}")
  println("  Variables: \{summary.n_variables}")
  println("  Dimensions: \{summary.n_dimensions}")
  println("  Total inertia: \{summary.total_inertia}")
  println("  First dimension variance: \{summary.first_dim_variance}%")

  // Get individual coordinates
  let ind = @src.facto_get_pca_ind(pca_result)
  println("")
  println("Individual coordinates (Dim1, Dim2):")
  let mut i = 0
  while i < ind.coord.length() {
    println("  Sample \{i + 1}: (\{ind.coord[i][0]}, \{ind.coord[i][1]})")
    i = i + 1
  }

  // Get variable coordinates
  let pca_var = @src.facto_get_pca_var(pca_result)
  println("")
  println("Variable loadings (Dim1, Dim2):")
  let mut j = 0
  while j < pca_var.coord.length() {
    println("  Gene \{j + 1}: (\{pca_var.coord[j][0]}, \{pca_var.coord[j][1]})")
    j = j + 1
  }

  // Get contributions to dimension 1
  println("")
  println("Top variable contributions to Dim 1:")
  let contribs = @src.facto_var_contrib_dim(pca_result, 1)
  let mut k = 0
  while k < 3 {
    println("  \{contribs[k].0}: \{contribs[k].1}%")
    k = k + 1
  }

  // Number of dimensions for 80% variance
  let nb_dim = @src.facto_nb_dim(pca_result, 80.0)
  println("")
  println("Dimensions needed for 80% variance: \{nb_dim}")

  // Print full summary
  println("")
  println(@src.facto_print_summary(pca_result))

  println("=== factoextra Demo Complete ===")
}
