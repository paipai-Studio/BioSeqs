///|
/// BiocSingular Demo: Singular Value Decomposition for dimensionality reduction.
/// Demonstrates exact SVD, IRLBA, and randomized SVD methods.
///
/// Run: moon run examples/bioc_singular_demo/main.mbt
fn main {
  println("=== BiocSingular Demo ===")
  println("")

  // === 1. Create test matrix ===
  println("--- 1. Create Test Matrix ---")
  let (matrix, nrow, ncol) = @bio.create_test_matrix()
  println("  Matrix size: \{nrow} x \{ncol}")
  for i = 0; i < nrow; i = i + 1 {
    let mut row_str = ""
    for j = 0; j < ncol; j = j + 1 {
      let val = matrix[i * ncol + j]
      row_str = "\{row_str} \{val}"
    }
    println("  Row \{i}:[\{row_str} ]")
  }
  println("")

  // === 2. SVDOptions with different methods ===
  println("--- 2. SVDOptions Construction ---")

  let opts_default = @bio.svd_options_new()
  println("  Default options: method=ExactSVD, rank=\{opts_default.rank}, tol=\{opts_default.tol}, max_iter=\{opts_default.max_iter}")

  let opts_exact = @bio.svd_options_with_method(@bio.exact_svd_method())
  println("  ExactSVD options: method=ExactSVD")

  let opts_irlba = @bio.svd_options_with_method(@bio.irlba_method())
  println("  IRLBA options: method=IRLBA")

  let opts_rand = @bio.svd_options_with_method(@bio.randomized_method())
  println("  Randomized options: method=Randomized")

  let opts_full = @bio.svd_options_full(3, 1.0e-8, 2000, 15, @bio.exact_svd_method())
  println("  Full options: rank=\{opts_full.rank}, tol=\{opts_full.tol}, max_iter=\{opts_full.max_iter}, n_oversamples=\{opts_full.n_oversamples}")
  println("")

  // === 3. Run exact SVD ===
  println("--- 3. Exact SVD ---")
  let result_exact = @bio.run_svd(matrix, nrow, ncol, opts_exact)
  let d_exact = @bio.sv_result_d(result_exact)
  let u_exact = @bio.sv_result_u(result_exact)
  let v_exact = @bio.sv_result_v(result_exact)
  println("  Singular values (d):")
  for i = 0; i < d_exact.length(); i = i + 1 {
    println("    d[\{i}] = \{d_exact[i]}")
  }

  // Check sorted descending
  let mut sorted_ok = true
  for i = 0; i < d_exact.length() - 1; i = i + 1 {
    if d_exact[i] < d_exact[i + 1] {
      sorted_ok = false
    }
  }
  println("  Singular values sorted descending: \{sorted_ok}")
  println("  Rank: \{@bio.sv_result_rank(result_exact)}")
  println("")

  // === 4. Run IRLBA (falls back to exact SVD when rank >= min(nrow,ncol)) ===
  println("--- 4. IRLBA (Implicitly Restarted Lanczos Bidiagonalization) ---")
  let mn = if nrow < ncol { nrow } else { ncol }
  println("  min(nrow, ncol) = \{mn}; using rank=\{mn} to ensure valid decomposition")
  let result_irlba = @bio.run_svd(matrix, nrow, ncol, @bio.svd_options_full(mn, 1.0e-7, 1000, 10, @bio.irlba_method()))
  let d_irlba = @bio.sv_result_d(result_irlba)
  println("  Singular values:")
  for i = 0; i < d_irlba.length(); i = i + 1 {
    println("    d[\{i}] = \{d_irlba[i]}")
  }
  println("  Rank: \{@bio.sv_result_rank(result_irlba)}")
  println("  Note: When rank >= min(nrow,ncol), IRLBA falls back to exact SVD")
  println("  opts_irlba default rank: \{opts_irlba.rank}, method: IRLBA")
  println("")

  // === 5. Run Randomized SVD ===
  println("--- 5. Randomized SVD ---")
  let result_rand = @bio.run_svd(matrix, nrow, ncol, @bio.svd_options_full(mn, 1.0e-7, 1000, 10, @bio.randomized_method()))
  let d_rand = @bio.sv_result_d(result_rand)
  println("  Singular values:")
  for i = 0; i < d_rand.length(); i = i + 1 {
    println("    d[\{i}] = \{d_rand[i]}")
  }
  println("  Rank: \{@bio.sv_result_rank(result_rand)}")
  println("  Note: When rank >= min(nrow,ncol), Randomized SVD falls back to exact SVD")
  println("  opts_rand default rank: \{opts_rand.rank}, method: Randomized")
  println("")

  // === 6. Compare results across methods ===
  println("--- 6. Compare SVD Methods ---")
  println("  Comparing singular values from ExactSVD, IRLBA, and Randomized:")
  println("  Method           d[0]         d[1]         d[2]         d[3]")
  println("  ----             ------       ------       ------       ------")
  println("  ExactSVD         \{d_exact[0]}  \{d_exact[1]}  \{d_exact[2]}  \{d_exact[3]}")
  println("  IRLBA            \{d_irlba[0]}  \{d_irlba[1]}  \{d_irlba[2]}  \{d_irlba[3]}")
  println("  Randomized       \{d_rand[0]}  \{d_rand[1]}  \{d_rand[2]}  \{d_rand[3]}")

  // Verify all methods produce the same results
  let mut all_match = true
  for i = 0; i < d_exact.length(); i = i + 1 {
    if (d_exact[i] - d_irlba[i]).abs() > 0.01 { all_match = false }
    if (d_exact[i] - d_rand[i]).abs() > 0.01 { all_match = false }
  }
  println("")
  println("  All methods produce consistent results: \{all_match}")
  println("")

  // === 7. Demonstrate reconstruction: A = U * Sigma * V^T ===
  println("--- 7. Reconstruction: A = U * Sigma * V^T ---")
  let rank = @bio.sv_result_rank(result_exact)

  // Compute reconstruction directly: result[i][j] = sum_k d[k] * U[i,k] * V[j,k]
  // U is stored as column vectors: u_exact[k] is column k (length nrow)
  // V is stored as column vectors: v_exact[k] is column k (length ncol)
  let reconstruction : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let row = Array::make(ncol, 0.0)
    for j = 0; j < ncol; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < rank; k = k + 1 {
        sum = sum + d_exact[k] * u_exact[k][i] * v_exact[k][j]
      }
      row[j] = sum
    }
    reconstruction.push(row)
  }

  println("  Original matrix vs Reconstructed matrix:")
  let mut total_error = 0.0
  for i = 0; i < nrow; i = i + 1 {
    let mut orig_row = ""
    let mut recon_row = ""
    for j = 0; j < ncol; j = j + 1 {
      let orig_val = matrix[i * ncol + j]
      let recon_val = reconstruction[i][j]
      orig_row = "\{orig_row} \{orig_val}"
      recon_row = "\{recon_row} \{recon_val}"
      let diff = (orig_val - recon_val).abs()
      total_error = total_error + diff
    }
    println("  Row \{i}:")
    println("    Original:   [\{orig_row} ]")
    println("    Reconstruct: [\{recon_row} ]")
  }
  println("  Total reconstruction error (L1): \{total_error}")
  println("  Reconstruction accurate: \{total_error < 0.01}")
  println("")

  // === 8. Show orthonormality of U and V ===
  println("--- 8. Orthonormality of U and V ---")

  // U is stored as column vectors: u_exact[col][row]
  // U^T * U element (i,j) = sum_k U[k,i] * U[k,j] = sum_k u_exact[i][k] * u_exact[j][k]
  println("  Checking U^T * U = I:")
  let mut u_orth_error = 0.0
  for i = 0; i < rank; i = i + 1 {
    for j = 0; j < rank; j = j + 1 {
      let mut dot = 0.0
      for k = 0; k < nrow; k = k + 1 {
        dot = dot + u_exact[i][k] * u_exact[j][k]
      }
      let expected = if i == j { 1.0 } else { 0.0 }
      u_orth_error = u_orth_error + (dot - expected).abs()
    }
  }
  println("    U orthonormality error (L1): \{u_orth_error}")
  println("    U is orthonormal: \{u_orth_error < 0.01}")

  // V is stored as column vectors: v_exact[col][row]
  // V^T * V element (i,j) = sum_k V[k,i] * V[k,j] = sum_k v_exact[i][k] * v_exact[j][k]
  println("  Checking V^T * V = I:")
  let mut v_orth_error = 0.0
  for i = 0; i < rank; i = i + 1 {
    for j = 0; j < rank; j = j + 1 {
      let mut dot = 0.0
      for k = 0; k < ncol; k = k + 1 {
        dot = dot + v_exact[i][k] * v_exact[j][k]
      }
      let expected = if i == j { 1.0 } else { 0.0 }
      v_orth_error = v_orth_error + (dot - expected).abs()
    }
  }
  println("    V orthonormality error (L1): \{v_orth_error}")
  println("    V is orthonormal: \{v_orth_error < 0.01}")
  println("")

  // === 9. Explained variance contribution ===
  println("--- 9. Explained Variance Contribution ---")
  let mut total_var = 0.0
  for i = 0; i < d_exact.length(); i = i + 1 {
    total_var = total_var + d_exact[i] * d_exact[i]
  }
  println("  Total variance (sum of d^2): \{total_var}")
  println("  Individual contribution:")
  let mut cum_var = 0.0
  for i = 0; i < d_exact.length(); i = i + 1 {
    let d_sq = d_exact[i] * d_exact[i]
    let contrib = d_sq / total_var * 100.0
    cum_var = cum_var + contrib
    println("    Component \{i}: d^2 = \{d_sq}, variance = \{contrib}%, cumulative = \{cum_var}%")
  }
  println("")

  // === 10. Truncated SVD convenience function ===
  println("--- 10. Truncated SVD (top 3 components) ---")
  let result_trunc = @bio.run_svd_truncated(matrix, nrow, ncol, 3)
  let d_trunc = @bio.sv_result_d(result_trunc)
  println("  Truncated singular values:")
  for i = 0; i < d_trunc.length(); i = i + 1 {
    println("    d[\{i}] = \{d_trunc[i]}")
  }
  println("  Truncated rank: \{@bio.sv_result_rank(result_trunc)}")
  println("")

  // === 11. Accessor functions ===
  println("--- 11. SVResult Accessor Functions ---")
  println("  sv_result_nrow: \{@bio.sv_result_nrow(result_exact)}")
  println("  sv_result_ncol: \{@bio.sv_result_ncol(result_exact)}")
  println("  sv_result_rank: \{@bio.sv_result_rank(result_exact)}")
  println("  sv_result_d length: \{@bio.sv_result_d(result_exact).length()}")
  println("  sv_result_u dimensions: \{@bio.sv_result_u(result_exact).length()} x \{@bio.sv_result_u(result_exact)[0].length()}")
  println("  sv_result_v dimensions: \{@bio.sv_result_v(result_exact).length()} x \{@bio.sv_result_v(result_exact)[0].length()}")
  println("")

  println("=== Demo Complete ===")
}