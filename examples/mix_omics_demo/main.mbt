///|
/// mixOmics Demo: Multi-omics integration methods.
/// Demonstrates PLS, sparse PLS (sPLS), and DIABLO algorithms
/// for integrating multiple omics data modalities.
///
/// Run: moon run examples/mix_omics_demo/main.mbt
fn main {
  println("=== mixOmics Demo ===")
  println("")

  // === 1. Create simulated multi-omics data ===
  println("--- 1. Create Simulated Multi-Omics Data ---")

  // Gene expression data: 20 samples x 10 genes
  let x_data : Array[Array[Double]] = [
    [2.1, 3.2, 1.5, 4.0, 2.8, 3.5, 1.9, 4.2, 3.1, 2.7],
    [1.8, 2.9, 1.2, 3.8, 2.5, 3.2, 1.7, 4.0, 2.9, 2.4],
    [2.5, 3.6, 1.8, 4.3, 3.0, 3.8, 2.2, 4.5, 3.4, 2.9],
    [1.2, 2.4, 0.8, 3.2, 2.0, 2.7, 1.3, 3.5, 2.3, 1.8],
    [3.0, 4.1, 2.2, 4.8, 3.5, 4.3, 2.6, 5.0, 3.9, 3.3],
    [0.8, 1.9, 0.5, 2.8, 1.5, 2.2, 0.9, 3.1, 1.8, 1.4],
    [2.7, 3.8, 2.0, 4.5, 3.2, 4.0, 2.4, 4.7, 3.6, 3.1],
    [1.5, 2.6, 1.0, 3.5, 2.2, 2.9, 1.5, 3.7, 2.6, 2.0],
    [3.2, 4.3, 2.4, 5.0, 3.7, 4.5, 2.8, 5.2, 4.1, 3.5],
    [1.0, 2.1, 0.7, 3.0, 1.7, 2.4, 1.1, 3.3, 2.0, 1.6],
    [2.3, 3.4, 1.6, 4.1, 2.9, 3.6, 2.0, 4.3, 3.2, 2.8],
    [1.6, 2.7, 1.1, 3.6, 2.3, 3.0, 1.6, 3.8, 2.7, 2.1],
    [2.8, 3.9, 2.1, 4.6, 3.3, 4.1, 2.5, 4.8, 3.7, 3.2],
    [1.3, 2.5, 0.9, 3.3, 2.1, 2.8, 1.4, 3.6, 2.4, 1.9],
    [3.1, 4.2, 2.3, 4.9, 3.6, 4.4, 2.7, 5.1, 4.0, 3.4],
    [0.9, 2.0, 0.6, 2.9, 1.6, 2.3, 1.0, 3.2, 1.9, 1.5],
    [2.6, 3.7, 1.9, 4.4, 3.1, 3.9, 2.3, 4.6, 3.5, 3.0],
    [1.4, 2.5, 1.0, 3.4, 2.2, 2.9, 1.5, 3.6, 2.5, 2.0],
    [2.9, 4.0, 2.2, 4.7, 3.4, 4.2, 2.6, 4.9, 3.8, 3.3],
    [1.1, 2.3, 0.8, 3.1, 1.8, 2.5, 1.2, 3.4, 2.2, 1.7],
  ]

  // Clinical traits data: 20 samples x 3 traits
  let y_data : Array[Array[Double]] = [
    [5.2, 3.1, 10.5],
    [4.8, 2.9, 9.8],
    [5.5, 3.4, 11.2],
    [3.9, 2.3, 8.5],
    [6.0, 3.8, 12.0],
    [3.2, 1.9, 7.2],
    [5.7, 3.5, 11.5],
    [4.3, 2.6, 9.0],
    [6.2, 4.0, 12.5],
    [3.5, 2.1, 7.8],
    [5.3, 3.2, 10.8],
    [4.5, 2.7, 9.2],
    [5.9, 3.7, 11.8],
    [3.8, 2.2, 8.2],
    [6.1, 3.9, 12.2],
    [3.3, 2.0, 7.5],
    [5.6, 3.4, 11.3],
    [4.1, 2.5, 8.8],
    [6.0, 3.8, 12.0],
    [3.6, 2.2, 8.0],
  ]

  let n_samples = x_data.length()
  let n_genes = 10
  let n_traits = 3
  println("  Gene expression (X): \{n_samples} samples x \{n_genes} genes")
  println("  Clinical traits (Y): \{n_samples} samples x \{n_traits} traits")
  println("  X[0] (first 5 genes): [\{x_data[0][0]}, \{x_data[0][1]}, \{x_data[0][2]}, \{x_data[0][3]}, \{x_data[0][4]}]")
  println("  Y[0]: [\{y_data[0][0]}, \{y_data[0][1]}, \{y_data[0][2]}]")
  println("")

  // === 2. Run PLS regression ===
  println("--- 2. PLS Regression (2 components) ---")

  let pls_opts = @bio.pls_options_full(2, "regression", true, 500, 1.0e-6)
  println("  Options: ncomp=\{pls_opts.ncomp}, mode=\{pls_opts.mode}, scale=\{pls_opts.scale}")

  let pls_result = @bio.run_pls(x_data, y_data, pls_opts)

  // Display scores
  println("")
  println("  X Scores (first 5 samples):")
  for i = 0; i < 5; i = i + 1 {
    println("    Sample \{i}: [\{pls_result.x_scores[i][0]}, \{pls_result.x_scores[i][1]}]")
  }

  println("")
  println("  Y Scores (first 5 samples):")
  for i = 0; i < 5; i = i + 1 {
    println("    Sample \{i}: [\{pls_result.y_scores[i][0]}, \{pls_result.y_scores[i][1]}]")
  }

  // Display loadings
  println("")
  println("  X Loadings (gene weights):")
  println("    Gene     Comp1        Comp2")
  println("    ----     ------       ------")
  let gene_names = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ"]
  for j = 0; j < n_genes; j = j + 1 {
    println("    \{gene_names[j]}  \{pls_result.x_loadings[j][0]}  \{pls_result.x_loadings[j][1]}")
  }

  println("")
  println("  Y Loadings (trait weights):")
  println("    Trait    Comp1        Comp2")
  println("    ----     ------       ------")
  let trait_names = ["Trait1", "Trait2", "Trait3"]
  for k = 0; k < n_traits; k = k + 1 {
    println("    \{trait_names[k]}  \{pls_result.y_loadings[k][0]}  \{pls_result.y_loadings[k][1]}")
  }

  // Explained variance
  println("")
  println("  Explained Variance:")
  for comp = 0; comp < pls_result.ncomp; comp = comp + 1 {
    let (x_var, y_var) = pls_result.explained_variance[comp]
    println("    Component \{comp + 1}: X = \{x_var} (\{x_var * 100.0}%), Y = \{y_var} (\{y_var * 100.0}%)")
  }
  println("")

  // === 3. Run sparse PLS ===
  println("--- 3. Sparse PLS (sPLS) with Variable Selection ---")

  let spls_opts = @bio.spls_options_full(2, [5, 5], [2, 2], true, 500, 1.0e-6)
  println("  Options: ncomp=\{spls_opts.ncomp}, keep_x=\{spls_opts.keep_x}, keep_y=\{spls_opts.keep_y}")

  let spls_result = @bio.run_spls(x_data, y_data, spls_opts)

  // Show selected variables
  println("")
  println("  Selected X variables (genes):")
  for comp = 0; comp < spls_result.ncomp; comp = comp + 1 {
    let selected = Array::make(0, 0)
    for j = 0; j < n_genes; j = j + 1 {
      if spls_result.x_selected[j][comp] {
        selected.push(j)
      }
    }
    let mut names = ""
    for i = 0; i < selected.length(); i = i + 1 {
      let idx = selected[i]
      let comma = if i < selected.length() - 1 { ", " } else { "" }
      names = "\{names}\{gene_names[idx]}\{comma}"
    }
    println("    Component \{comp + 1}: [\{names}]")
  }

  println("")
  println("  Selected Y variables (traits):")
  for comp = 0; comp < spls_result.ncomp; comp = comp + 1 {
    let selected = Array::make(0, 0)
    for k = 0; k < n_traits; k = k + 1 {
      if spls_result.y_selected[k][comp] {
        selected.push(k)
      }
    }
    let mut names = ""
    for i = 0; i < selected.length(); i = i + 1 {
      let idx = selected[i]
      let comma = if i < selected.length() - 1 { ", " } else { "" }
      names = "\{names}\{trait_names[idx]}\{comma}"
    }
    println("    Component \{comp + 1}: [\{names}]")
  }

  // Display sPLS loadings
  println("")
  println("  sPLS X Loadings (gene weights):")
  println("    Gene     Comp1        Comp2")
  println("    ----     ------       ------")
  for j = 0; j < n_genes; j = j + 1 {
    println("    \{gene_names[j]}  \{spls_result.x_loadings[j][0]}  \{spls_result.x_loadings[j][1]}")
  }

  // Explained variance
  println("")
  println("  sPLS Explained Variance:")
  for comp = 0; comp < spls_result.ncomp; comp = comp + 1 {
    let (x_var, y_var) = spls_result.explained_variance[comp]
    println("    Component \{comp + 1}: X = \{x_var} (\{x_var * 100.0}%), Y = \{y_var} (\{y_var * 100.0}%)")
  }
  println("")

  // === 4. Run DIABLO ===
  println("--- 4. DIABLO Multi-Block Integration ---")

  // Create 2 blocks: genomics (5 features) + transcriptomics (5 features)
  let block1 : Array[Array[Double]] = [
    [1.2, 2.3, 0.8, 3.1, 2.5],
    [1.0, 2.0, 0.6, 2.9, 2.2],
    [1.5, 2.6, 1.0, 3.3, 2.8],
    [0.8, 1.7, 0.4, 2.4, 1.9],
    [1.8, 2.9, 1.2, 3.6, 3.1],
    [0.5, 1.4, 0.3, 2.0, 1.5],
    [1.6, 2.7, 1.1, 3.4, 2.9],
    [0.9, 1.9, 0.5, 2.7, 2.1],
    [2.0, 3.1, 1.4, 3.8, 3.3],
    [0.6, 1.5, 0.3, 2.2, 1.7],
    [1.4, 2.5, 0.9, 3.2, 2.7],
    [1.1, 2.1, 0.7, 2.8, 2.3],
    [1.7, 2.8, 1.2, 3.5, 3.0],
    [0.7, 1.6, 0.4, 2.5, 1.8],
    [1.9, 3.0, 1.3, 3.7, 3.2],
    [0.5, 1.4, 0.3, 2.1, 1.6],
    [1.5, 2.6, 1.0, 3.3, 2.8],
    [0.8, 1.8, 0.5, 2.6, 2.0],
    [1.8, 2.9, 1.2, 3.6, 3.1],
    [0.6, 1.5, 0.3, 2.3, 1.7],
  ]

  let block2 : Array[Array[Double]] = [
    [3.1, 5.2, 2.8, 4.5, 3.9],
    [2.9, 4.8, 2.5, 4.2, 3.6],
    [3.4, 5.5, 3.1, 4.8, 4.2],
    [2.3, 3.9, 1.8, 3.5, 3.0],
    [3.8, 6.0, 3.5, 5.2, 4.5],
    [1.9, 3.2, 1.5, 2.8, 2.4],
    [3.5, 5.7, 3.2, 4.9, 4.3],
    [2.6, 4.3, 2.1, 3.7, 3.2],
    [4.0, 6.2, 3.7, 5.4, 4.7],
    [2.1, 3.5, 1.7, 3.0, 2.6],
    [3.2, 5.3, 2.9, 4.6, 4.0],
    [2.7, 4.5, 2.3, 3.9, 3.3],
    [3.7, 5.9, 3.4, 5.1, 4.4],
    [2.2, 3.7, 1.8, 3.2, 2.7],
    [3.9, 6.1, 3.6, 5.3, 4.6],
    [2.0, 3.3, 1.6, 2.9, 2.5],
    [3.4, 5.6, 3.1, 4.8, 4.2],
    [2.5, 4.1, 2.0, 3.6, 3.1],
    [3.8, 6.0, 3.5, 5.2, 4.5],
    [2.2, 3.6, 1.8, 3.1, 2.6],
  ]

  let blocks = [block1, block2]
  let block_names = ["Genomics", "Transcriptomics"]

  // DIABLO options with variable selection
  let diablo_opts = @bio.diablo_options_full(
    2,
    [[3, 3], [3, 3]],
    [[0.0, 1.0], [1.0, 0.0]],
    true,
    500,
    1.0e-6,
  )
  println("  Options: ncomp=\{diablo_opts.ncomp}, n_blocks=2, scale=\{diablo_opts.scale}")
  println("  Keep variables per component:")
  for comp = 0; comp < diablo_opts.ncomp; comp = comp + 1 {
    let kv = diablo_opts.keep_variables[comp]
    println("    Component \{comp + 1}: \{kv}")
  }

  let diablo_result = @bio.run_diablo(blocks, block_names, diablo_opts)

  println("")
  println("  Block scores (first 5 samples):")
  for b = 0; b < diablo_result.n_blocks; b = b + 1 {
    println("    \{block_names[b]}:")
    for i = 0; i < 5; i = i + 1 {
      println("      Sample \{i}: [\{diablo_result.block_scores[b][i][0]}, \{diablo_result.block_scores[b][i][1]}]")
    }
  }

  println("")
  println("  Block loadings:")
  for b = 0; b < diablo_result.n_blocks; b = b + 1 {
    println("    \{block_names[b]}:")
    let pb = diablo_result.block_sizes[b]
    let feature_names = if b == 0 {
      ["SNP_A", "SNP_B", "SNP_C", "SNP_D", "SNP_E"]
    } else {
      ["Exp_A", "Exp_B", "Exp_C", "Exp_D", "Exp_E"]
    }
    println("      Feature   Comp1        Comp2")
    println("      ------    ------       ------")
    for j = 0; j < pb; j = j + 1 {
      println("      \{feature_names[j]}  \{diablo_result.block_loadings[b][j][0]}  \{diablo_result.block_loadings[b][j][1]}")
    }
  }

  println("")
  println("  Variable selection per block:")
  for b = 0; b < diablo_result.n_blocks; b = b + 1 {
    println("    \{block_names[b]}:")
    let pb = diablo_result.block_sizes[b]
    for comp = 0; comp < diablo_result.ncomp; comp = comp + 1 {
      let selected = Array::make(0, 0)
      for j = 0; j < pb; j = j + 1 {
        if diablo_result.block_selected[b][j][comp] {
          selected.push(j)
        }
      }
      let mut names = ""
      for i = 0; i < selected.length(); i = i + 1 {
        let idx = selected[i]
        let comma = if i < selected.length() - 1 { ", " } else { "" }
        let fname = if b == 0 {
          ["SNP_A", "SNP_B", "SNP_C", "SNP_D", "SNP_E"][idx]
        } else {
          ["Exp_A", "Exp_B", "Exp_C", "Exp_D", "Exp_E"][idx]
        }
        names = "\{names}\{fname}\{comma}"
      }
      println("      Component \{comp + 1}: [\{names}]")
    }
  }
  println("")

  // === 5. Compare PLS vs sPLS ===
  println("--- 5. Compare PLS vs sPLS ---")

  println("  X Loadings comparison (first 5 genes, Component 1):")
  println("    Gene     PLS          sPLS")
  println("    ----     ---          ----")
  for j = 0; j < 5; j = j + 1 {
    println("    \{gene_names[j]}  \{pls_result.x_loadings[j][0]}  \{spls_result.x_loadings[j][0]}")
  }

  println("")
  println("  Variance explained comparison:")
  println("    Component  PLS_X       sPLS_X     PLS_Y       sPLS_Y")
  println("    ---------  -------     -------    -------     -------")
  for comp = 0; comp < pls_result.ncomp; comp = comp + 1 {
    let (pls_x, pls_y) = pls_result.explained_variance[comp]
    let (spls_x, spls_y) = spls_result.explained_variance[comp]
    println("    Comp \{comp + 1}     \{pls_x}  \{spls_x}  \{pls_y}  \{spls_y}")
  }

  println("")
  println("  Key differences:")
  println("    - PLS uses all variables, sPLS selects top-k per component")
  println("    - sPLS sets unselected variable loadings to exactly 0")
  println("    - sPLS provides variable selection via x_selected/y_selected flags")
  println("    - PLS may explain slightly more variance but includes noise variables")

  println("")
  println("=== Demo Complete ===")
}