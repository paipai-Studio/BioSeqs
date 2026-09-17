///|
/// pcaMethods PPCA + BPCA demo.
///
/// Runs probabilistic PCA (Tipping-Bishop) and Bayesian PCA (Oba 2003) on a
/// small centered matrix with two missing (NaN) entries, then prints the
/// estimated noise variance / precision, loadings, and the imputed values.

///|
fn nan() -> Double {
  0.0 / 0.0
}

///|
fn main raise {
  let data : Array[Array[Double]] = [
    [0.508, 0.242, -0.142, 0.320, 0.012],
    [-0.234, 0.469, nan(), -0.501, -0.751],
    [0.759, 0.024, -0.698, 0.706, 0.520],
    [-1.055, -0.725, 0.143, -0.418, 0.411],
    [0.458, 0.750, 0.325, -0.050, -0.646],
    [nan(), -0.939, -0.357, 0.048, 0.940],
    [1.120, 0.371, -0.668, 0.703, 0.176],
    [-0.860, 0.007, 0.805, -0.779, -0.513],
  ]
  println("=== PPCA (probabilistic PCA, 2 components) ===")
  let p = @src.ppca(matrix=data, n_pcs=2, max_iterations=2000, threshold=1.0e-5)
  println("noise variance (ss) = \{p.sigma2}")
  println("loadings row 0    = \{p.loadings[0][0]}, \{p.loadings[0][1]}")
  println("cumulative R2     = \{p.r2cum[0]}, \{p.r2cum[1]}")

  println("\n=== BPCA (Bayesian PCA, 2 components) ===")
  let b = @src.bpca(matrix=data, n_pcs=2, max_steps=2000, threshold=1.0e-4)
  println("precision (tau)   = \{b.tau}")
  println("loadings row 0    = \{b.loadings[0][0]}, \{b.loadings[0][1]}")
  println("imputed data[1][2] = \{b.imputed[1][2]}")
  println("imputed data[5][0] = \{b.imputed[5][0]}")
}