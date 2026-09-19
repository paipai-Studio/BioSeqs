///|
/// scikit-bio ordination + mantel demo.
///
/// Shows Gower centering of a distance matrix, Correspondence Analysis of a
/// contingency table, and the Mantel test between two distance matrices.

///|
fn main {
  println("=== Canonical Correspondence Analysis (CCA, scaling=1) ===")
  let y = [[1.0, 0.0, 2.0], [0.0, 3.0, 1.0], [2.0, 1.0, 0.0], [3.0, 1.0, 2.0], [1.0, 2.0, 3.0], [2.0, 0.0, 2.0]]
  let x = [[0.5, 1.0], [1.0, 0.5], [1.5, 2.0], [2.5, 1.5], [1.0, 3.0], [2.0, 2.5]]
  let cca = @src.skbio_cca(y, x, 1)
  println("eigvals   : \{cca.eigvals}")
  println("prop. expl: \{cca.proportion_explained}")

  println("=== Redundancy Analysis (RDA, scaling=1) ===")
  let y = [[1.0, 0.0, 2.0], [0.0, 3.0, 1.0], [2.0, 1.0, 0.0], [3.0, 1.0, 2.0], [1.0, 2.0, 3.0], [2.0, 0.0, 2.0]]
  let x = [[0.5, 1.0], [1.0, 0.5], [1.5, 2.0], [2.5, 1.5], [1.0, 3.0], [2.0, 2.5]]
  let rda = @src.skbio_rda(y, x, 1)
  println("eigvals   : \{rda.eigvals}")
  println("prop. expl: \{rda.proportion_explained}")
  println("sample 0  : \{rda.sample_scores[0]}")


  println("=== center_distance_matrix ===")
  let dm = [[0.0, 1.0, 2.0], [1.0, 0.0, 1.0], [2.0, 1.0, 0.0]]
  let c = @src.skbio_center_distance_matrix(dm)
  for row in c {
    println(row)
  }

  println("\n=== Correspondence Analysis (scaling=1) ===")
  let x = [
    [1.0, 0.0, 2.0, 1.0],
    [0.0, 3.0, 1.0, 0.0],
    [2.0, 1.0, 0.0, 1.0],
  ]
  let ca = @src.skbio_ca(x, 1)
  println("eigvals      : \{ca.eigvals}")
  println("prop. expl.  : \{ca.proportion_explained}")
  println("sample row 0 : \{ca.sample_scores[0]}")

  println("\n=== Mantel test ===")
  let y = [[0.0, 2.0, 7.0], [2.0, 0.0, 5.0], [7.0, 5.0, 0.0]]
  let m = @src.skbio_mantel(dm, y, 999, 42)
  println("corr coeff   : \{m.corr_coeff}")
  println("p-value      : \{m.p_value}")
  println("n            : \{m.n}")
}