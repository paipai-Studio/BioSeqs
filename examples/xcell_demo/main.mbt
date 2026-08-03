///|
/// xCell Demo.
///
/// This example demonstrates immune/stromal cell type scoring using
/// single-sample GSEA (ssGSEA) as in Bioconductor's xCell package, including:
///   1. Understanding the default 64 immune+stromal signatures
///   2. Single-sample GSEA scoring on custom gene sets
///   3. Full xCell pipeline on expression matrix
///   4. Querying scores by cell type, sample, or category
///   5. Immune / Stromal / microenvironment scores summary
///
/// Inspired by Aran, Hu, Butte (2017); xCell Bioconductor package.

fn top_n(arr : Array[Double], n : Int) -> Array[(Int, Double)] {
  let pairs : Array[(Int, Double)] = []
  for i = 0; i < arr.length(); i = i + 1 {
    pairs.push((i, arr[i]))
  }
  // Bubble-sort descending by score
  for i = 0; i < pairs.length(); i = i + 1 {
    for j = i + 1; j < pairs.length(); j = j + 1 {
      if pairs[j].1 > pairs[i].1 {
        let tmp = pairs[i]
        pairs[i] = pairs[j]
        pairs[j] = tmp
      }
    }
  }
  let result : Array[(Int, Double)] = []
  let k = if n < pairs.length() { n } else { pairs.length() }
  for i = 0; i < k; i = i + 1 {
    result.push(pairs[i])
  }
  result
}

fn format_score(s : Double) -> String {
  s.to_string().substring(start=0, end=if s < 0.0 { 7 } else { 6 })
}

fn pad_right(s : String, n : Int) -> String {
  if s.length() >= n { s } else { s + " ".repeat(n - s.length()) }
}
fn pad_left(s : String, n : Int) -> String {
  if s.length() >= n { s } else { " ".repeat(n - s.length()) + s }
}
fn string_index_of(s : String, sub : String) -> Int {
  let n = s.length()
  let kn = sub.length()
  if kn == 0 || kn > n { return -1 }
  let mut i = 0
  while i <= n - kn {
    let mut matched = true
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched { return i }
    i = i + 1
  }
  -1
}
fn array_index_of(arr : Array[String], target : String) -> Int {
  for i = 0; i < arr.length(); i = i + 1 {
    if arr[i] == target { return i }
  }
  -1
}
fn fmt_fixed(d : Double, decimals : Int) -> String {
  let s = d.to_string()
  let dot_idx = string_index_of(s, ".")
  if dot_idx >= 0 {
    let end = if dot_idx + 1 + decimals < s.length() { dot_idx + 1 + decimals } else { s.length() }
    s.substring(start=0, end=end)
  } else {
    let mut r = s + "."
    while r.length() < s.length() + 1 + decimals { r = r + "0" }
    r
  }
}

///|
fn main {
  println("=== xCell Demo ===")
  println("Immune and Stromal Cell Type Enrichment Scoring (ssGSEA)")
  println("Inspired by Aran et al. xCell Bioconductor package")
  println("")

  println("=== 1. Explore Default xCell Signatures ===")
  let sigs = @bio.xcell_default_signatures()
  println("  Total built-in signatures: \{sigs.length()}")

  let categories : Map[String, Int] = Map([], capacity=0)
  for idx = 0; idx < sigs.length(); idx = idx + 1 {
    let s = sigs[idx]
    let cat = s.category
    let c = categories.get_or_default(cat, 0)
    categories.set(cat, c + 1)
  }
  let cat_keys = categories.keys()
  println("\n  Signatures by category:")
  cat_keys.each(fn(cat) {
    println("    \{cat}: \{categories.get_or_default(cat, 0)} signatures")
  })

  println("\n  First 6 signatures sample:")
  for i = 0; i < 6; i = i + 1 {
    let s = sigs[i]
    println(
      "    [\{i}] \{s.cell_type} (\{s.category}): " +
      "\{s.genes.length()} marker genes, e.g. [\{s.genes[0]}, ...]",
    )
  }

  println("\n=== 2. Single-Sample ssGSEA Scoring (Custom Gene Sets) ===")
  // Mock gene expression profile
  let expression_profile = [
    120.0, 80.0, 220.0, 150.0, 95.0, 75.0, 310.0, 280.0, // GeneA-H
    140.0, 50.0, 340.0, 180.0, 45.0, 220.0, 290.0, 380.0, // GeneI-P
    260.0, 110.0, 420.0, 350.0, // GeneQ-T
  ]
  let genes_profile = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH",
    "GeneI", "GeneJ", "GeneK", "GeneL", "GeneM", "GeneN", "GeneO", "GeneP",
    "GeneQ", "GeneR", "GeneS", "GeneT",
  ]
  // Mock signature: up-regulated genes are GeneC, GeneG, GeneK, GeneM, GeneQ
  let set_genes = ["GeneC", "GeneG", "GeneK", "GeneM", "GeneQ", "GeneZ_notexpressed"]

  let single_score = @bio.xcell_ssgsea_single(
    expression_profile, genes_profile, set_genes, alpha=0.25,
  )
  println("  Expression: \{genes_profile.length()} genes")
  println("  Signature set: \{set_genes.length()} genes (\{set_genes.join(", ")})")
  println("  ssGSEA enrichment score (alpha=0.25): \{format_score(single_score)}")

  // With higher alpha
  let single_score2 = @bio.xcell_ssgsea_single(
    expression_profile, genes_profile, set_genes, alpha=0.75,
  )
  println("  ssGSEA enrichment score (alpha=0.75): \{format_score(single_score2)}")

  println("\n=== 3. Build Custom xCell Signatures ===")
  let custom_sigs = [
    @bio.XcellSignature::new(
      "Tcells_Custom", "Lymphoid",
      ["CD3D","CD3E","CD3G","CD4","CD8A","CD8B","LCK","ZAP70"],
    ),
    @bio.XcellSignature::new(
      "Bcells_Custom", "Lymphoid",
      ["CD19","MS4A1","CD79A","CD79B","BLNK","PAX5","IGHM","IGHD"],
    ),
    @bio.XcellSignature::new(
      "Myeloid_Custom", "Myeloid",
      ["LYZ","CD14","CSF1R","FCGR3A","CD68","S100A8","S100A9","MPO"],
    ),
    @bio.XcellSignature::new(
      "Fibroblasts_Custom", "Stromal",
      ["COL1A1","COL1A2","COL3A1","FN1","VIM","SPARC","DCN","LUM"],
    ),
    @bio.XcellSignature::new(
      "Endothelial_Custom", "Stromal",
      ["PECAM1","CDH5","VWF","KDR","FLT1","ENG","TEK","ESAM"],
    ),
  ]
  println("  Custom signatures built: \{custom_sigs.length()} types")
  for idx = 0; idx < custom_sigs.length(); idx = idx + 1 {
    let s = custom_sigs[idx]
    println("    - \{s.cell_type} [\{s.category}]: \{s.genes.length()} markers")
  }

  println("\n=== 4. Synthetic Bulk RNA-seq Expression Matrix ===")
  // 5 samples x 30 genes matrix, simulating different cell-type mixtures
  let gene_names_all = [
    "CD3D","CD3E","CD3G","CD4","CD8A","CD8B","LCK","ZAP70",        // T cell
    "CD19","MS4A1","CD79A","CD79B","BLNK","PAX5","IGHM","IGHD",      // B cell
    "LYZ","CD14","CSF1R","FCGR3A","CD68","S100A8","S100A9","MPO",    // Myeloid
    "COL1A1","COL1A2","VIM","FN1","PECAM1","CDH5",                   // Stromal/Endo
  ]
  let sample_names_x = [
    "Tumor_ImmuneHigh", "Tumor_ImmuneLow", "Tumor_Desmoplastic",
    "Blood_PBMC", "Normal_LymphNode",
  ]
  // Each row: one gene. Columns: samples.
  let exprs_data : Array[Array[Double]] = [
    // Sample 0: Immune high (high lymphoid+myeloid)
    // Sample 1: Immune low
    // Sample 2: Desmoplastic (high stromal)
    // Sample 3: PBMC (high all immune, low stromal)
    // Sample 4: Lymph node (very high lymphoid)
    // T cell genes:
    [280.0, 20.0, 35.0, 190.0, 350.0], // CD3D
    [260.0, 18.0, 30.0, 180.0, 340.0], // CD3E
    [180.0, 15.0, 25.0, 150.0, 290.0], // CD3G
    [140.0, 10.0, 20.0, 100.0, 200.0], // CD4
    [120.0, 8.0, 15.0, 90.0, 180.0],   // CD8A
    [90.0, 6.0, 12.0, 70.0, 140.0],    // CD8B
    [220.0, 16.0, 30.0, 170.0, 300.0], // LCK
    [200.0, 14.0, 25.0, 160.0, 280.0], // ZAP70
    // B cell genes:
    [120.0, 12.0, 25.0, 150.0, 220.0], // CD19
    [140.0, 10.0, 20.0, 160.0, 250.0], // MS4A1=CD20
    [90.0, 8.0, 15.0, 110.0, 180.0],   // CD79A
    [85.0, 7.0, 14.0, 105.0, 170.0],   // CD79B
    [75.0, 6.0, 12.0, 95.0, 150.0],    // BLNK
    [100.0, 9.0, 18.0, 120.0, 190.0],  // PAX5
    [200.0, 15.0, 30.0, 230.0, 280.0], // IGHM
    [110.0, 9.0, 16.0, 140.0, 200.0],  // IGHD
    // Myeloid genes:
    [320.0, 30.0, 80.0, 250.0, 160.0], // LYZ
    [240.0, 25.0, 60.0, 200.0, 120.0], // CD14
    [160.0, 15.0, 50.0, 140.0, 70.0],  // CSF1R
    [180.0, 18.0, 55.0, 160.0, 80.0],  // FCGR3A=CD16
    [220.0, 22.0, 70.0, 190.0, 95.0],  // CD68
    [280.0, 28.0, 90.0, 220.0, 100.0], // S100A8
    [260.0, 26.0, 85.0, 200.0, 95.0],  // S100A9
    [150.0, 20.0, 60.0, 180.0, 80.0],  // MPO
    // Stromal/Endothelial:
    [40.0, 30.0, 600.0, 10.0, 60.0],   // COL1A1
    [35.0, 28.0, 540.0, 9.0, 55.0],    // COL1A2
    [80.0, 40.0, 480.0, 25.0, 70.0],   // VIM
    [50.0, 32.0, 440.0, 12.0, 45.0],   // FN1
    [15.0, 10.0, 200.0, 40.0, 90.0],   // PECAM1=CD31
    [18.0, 12.0, 180.0, 35.0, 85.0],   // CDH5=VE-cad
  ]

  println("  Expression matrix:")
  println("    Genes: \{exprs_data.length()}")
  println("    Samples: \{sample_names_x.length()}")
  let mut header_line = ""
  for idx = 0; idx < sample_names_x.length(); idx = idx + 1 {
    let s = sample_names_x[idx]
    header_line = header_line + "    | " + pad_right(s, 22)
  }
  println(header_line)
  for i = 0; i < gene_names_all.length(); i = i + 1 {
    let mut row_line = "    " + pad_right(gene_names_all[i], 10) + " |"
    for j = 0; j < sample_names_x.length(); j = j + 1 {
      let v = exprs_data[i][j]
      row_line = row_line + " " + pad_left(v.to_string(), 6) + " |"
    }
    println(row_line)
  }

  println("\n=== 5. Run xCell Scoring (Custom Signatures) ===")
  let params = @bio.XcellParams::new()
  println("  ssGSEA alpha: \{params.alpha}")
  println("  (Using custom 5-signature set above)")
  let custom_result = @bio.xcell_score(
    exprs_data, gene_names_all, sample_names_x, custom_sigs, params,
  )
  println("  Result matrix: cell types=\{custom_result.cell_types.length()}, " +
    "samples=\{sample_names_x.length()}")

  // Pretty-print score matrix
  println("\n  --- xCell scores (custom signatures) ---")
  let mut line = "    " + pad_right("Cell type", 28)
  for idx = 0; idx < sample_names_x.length(); idx = idx + 1 {
    let s = sample_names_x[idx]
    line = line + " | " + pad_left(s, 14)
  }
  println(line)
  println("    " + "-".repeat(28 + 16 * sample_names_x.length()))
  for ct = 0; ct < custom_result.cell_types.length(); ct = ct + 1 {
    let name = custom_result.cell_types[ct]
    line = "    " + pad_right(name, 28)
    for j = 0; j < sample_names_x.length(); j = j + 1 {
      line = line + " | " + pad_left(format_score(custom_result.scores[ct][j]), 14)
    }
    println(line)
  }

  println("\n=== 6. Query Scores By Cell Type / Sample ===")
  let tcell_scores = custom_result.get_cell_type_scores("Tcells_Custom")
  println("  Tcells_Custom scores across samples:")
  for j = 0; j < sample_names_x.length(); j = j + 1 {
    println("    - \{sample_names_x[j]}: \{format_score(tcell_scores[j])}")
  }

  let sample1_name = sample_names_x[2]
  let sample_scores = custom_result.get_sample_scores(sample1_name)
  println("\n  Sample '\{sample1_name}' scores across cell types:")
  for ct = 0; ct < custom_result.cell_types.length(); ct = ct + 1 {
    println(
      "    - \{custom_result.cell_types[ct]}: \{format_score(sample_scores[ct])}",
    )
  }

  println("\n=== 7. Top Enriched Cell Types per Sample ===")
  for j = 0; j < sample_names_x.length(); j = j + 1 {
    let sname = sample_names_x[j]
    let top = custom_result.get_top_cell_types(sname, 3)
    println("  Sample '\{sname}' top 3 cell types:")
    for k = 0; k < top.length(); k = k + 1 {
      let (ct_name, ct_score) = top[k]
      println("    \{k + 1}. \{ct_name}: \{format_score(ct_score)}")
    }
  }

  println("\n=== 8. Scores Grouped By Category (Lymphoid/Myeloid/Stromal) ===")
  let cat_scores = custom_result.scores_by_category()
  let cat_names_ordered : Array[String] = ["Lymphoid", "Myeloid", "Stromal"]
  line = "    " + pad_right("Category", 16)
  for idx = 0; idx < sample_names_x.length(); idx = idx + 1 {
    let s = sample_names_x[idx]
    line = line + " | " + pad_left(s, 14)
  }
  println(line)
  println("    " + "-".repeat(16 + 16 * sample_names_x.length()))
  for idx = 0; idx < cat_names_ordered.length(); idx = idx + 1 {
    let cat = cat_names_ordered[idx]
    match cat_scores.get(cat) {
      Some(arr) => {
        line = "    " + pad_right(cat, 16)
        for j = 0; j < arr.length(); j = j + 1 {
          line = line + " | " + pad_left(format_score(arr[j]), 14)
        }
        println(line)
      }
      None => println("    " + cat + ": (no data)")
    }
  }

  println("\n=== 9. Run Default xCell (64 Immune+Stromal) Pipeline ===")
  println("  Using built-in 64 signatures on the same expression data...")
  let default_result = @bio.xcell_run_default(
    exprs_data, gene_names_all, sample_names_x,
  )
  println("  Completed! Cell type signatures scored: \{default_result.cell_types.length()}")

  // Top 5 per sample
  println("\n  Top 5 cell types (from default 64) per sample:")
  for j = 0; j < sample_names_x.length(); j = j + 1 {
    let sname = sample_names_x[j]
    let top5 = default_result.get_top_cell_types(sname, 5)
    println("  --- \{sname} ---")
    for k = 0; k < top5.length(); k = k + 1 {
      let (ct_name, ct_score) = top5[k]
      let found_idx = array_index_of(default_result.cell_types, ct_name)
      let cat_idx = if found_idx >= 0 {
        default_result.categories[found_idx]
      } else {
        "?"
      }
      println(
        "    \{k + 1}. " + pad_right(ct_name, 26) + " [" + pad_right(cat_idx, 10) + "]" +
        " score=\{format_score(ct_score)}",
      )
    }
  }

  println("\n=== 10. Immune / Stromal / Microenvironment Summary (default pipeline) ===")
  println("  " + pad_right("Sample", 26) +
    " | Immune | Stromal | Microenv | ImmuneHigh?")
  println("  " + "-".repeat(26 + 54))
  for j = 0; j < sample_names_x.length(); j = j + 1 {
    let sname = sample_names_x[j]
    let imm = default_result.immune_scores[j]
    let str = default_result.stromal_scores[j]
    let menv = default_result.microenvironment_scores[j]
    let high = if imm > 0.4 { "YES" } else { "no" }
    println(
      "  " + pad_right(sname, 26) + " | " + pad_left(format_score(imm), 6) + " |" +
      " " + pad_left(format_score(str), 7) + " | " + pad_left(format_score(menv), 8) + " |" +
      "  " + high,
    )
  }

  println("\n=== xCell Demo Complete ===")
}
