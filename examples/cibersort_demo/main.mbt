///|
/// CIBERSORT module demo - immune cell deconvolution from bulk RNA-seq.
///
/// Run: moon run examples/cibersort_demo
fn main {
  println("=== CIBERSORT Module Demo ===")
  println("Cell type deconvolution via Non-negative Least Squares (NNLS)")
  println("Reference: Newman AM et al. (2015) Nature Methods 12:453-457")
  println("")

  // ------------------------------------------------------------------
  // 1. Load the built-in LM22-style signature matrix
  // ------------------------------------------------------------------
  println("=== 1. Load signature matrix ===")
  let sig = @bio.cib_default_signature()
  println("Genes: \{sig.gene_names().length()}")
  println("Cell types: \{sig.cell_types().length()}")
  println("Cell types: \{sig.cell_types()}")
  println("")

  // ------------------------------------------------------------------
  // 2. Build a synthetic mixture matrix
  // ------------------------------------------------------------------
  println("=== 2. Build synthetic mixture samples ===")
  let mix = @bio.cib_sample_mixture()
  println("Samples: \{mix.sample_names().length()}")
  println("Sample names: \{mix.sample_names()}")
  println("Mixture genes: \{mix.gene_names().length()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Run CIBERSORT deconvolution
  // ------------------------------------------------------------------
  println("=== 3. Run CIBERSORT deconvolution ===")
  let result = @bio.cib_run(sig, mix)
  println("Deconvolution completed for \{result.sample_names().length()} samples")
  println("")

  // ------------------------------------------------------------------
  // 4. Display per-sample cell type fractions
  // ------------------------------------------------------------------
  println("=== 4. Per-sample cell type fractions ===")
  for r in result.results() {
    println("--- Sample: \{r.sample_id()} ---")
    let mut i = 0
    while i < r.cell_types().length() {
      let frac = r.fractions()[i]
      if frac > 0.01 {
        let pct = frac * 100.0
        let bar_len = (pct / 5.0).to_int()
        let mut bar = ""
        let mut bi = 0
        while bi < bar_len {
          bar = bar + "#"
          bi = bi + 1
        }
        println(
          "  \{r.cell_types()[i]}: \{pct.to_string()}% \{bar}",
        )
      }
      i = i + 1
    }
    println("  Fit quality: Pearson r=\{r.pearson_r()}, RMSE=\{r.rmse()}")
    println("")
  }

  // ------------------------------------------------------------------
  // 5. Use accessor functions
  // ------------------------------------------------------------------
  println("=== 5. Query specific cell type fractions ===")
  let s1 = result.results()[0]
  let t_cd8_frac = s1.get_fraction("T cells CD8")
  let b_frac = s1.get_fraction("B cells")
  let nk_frac = s1.get_fraction("NK cells")
  println("Sample1 composition:")
  println("  T cells CD8: \{(t_cd8_frac * 100.0).to_string()}%")
  println("  B cells:     \{(b_frac * 100.0).to_string()}%")
  println("  NK cells:    \{(nk_frac * 100.0).to_string()}%")
  println("")

  // ------------------------------------------------------------------
  // 6. Format the full result as a table
  // ------------------------------------------------------------------
  println("=== 6. Full result table ===")
  println(result.to_string())

  // ------------------------------------------------------------------
  // 7. Custom signature example
  // ------------------------------------------------------------------
  println("=== 7. Custom signature example ===")
  let custom_sig = @bio.CibSignatureMatrix::new(
    ["G1", "G2", "G3"],
    ["TypeA", "TypeB"],
    [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]],
  )
  let custom_mix = @bio.CibMixtureMatrix::new(
    ["G1", "G2", "G3"],
    ["Mix1"],
    [[0.6], [0.5], [0.4]],
  )
  let custom_result = @bio.cib_run(custom_sig, custom_mix)
  println("Custom deconvolution for \{custom_result.sample_names()[0]}:")
  let cr = custom_result.results()[0]
  let mut ci = 0
  while ci < cr.cell_types().length() {
    println("  \{cr.cell_types()[ci]}: \{cr.fractions()[ci]}")
    ci = ci + 1
  }
  println("")

  println("=== Demo Complete ===")
}
