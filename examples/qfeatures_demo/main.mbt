///|
/// QFeatures module demo - quantitative proteomics data management.
///
/// Run: moon run examples/qfeatures_demo
fn main {
  println("=== QFeatures Module Demo ===")
  println("Quantitative proteomics: PSM -> peptide -> protein aggregation")
  println("Reference: Gatto L et al. (2020) QFeatures, bioRxiv 2020.10.27.358273")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample QFeatures data
  // ------------------------------------------------------------------
  println("=== 1. Load sample QFeatures data ===")
  let qf = @bio.qf_sample_data()
  println(@bio.qf_summary(qf))
  println("")

  // ------------------------------------------------------------------
  // 2. Print summary (n_assays, assay_names, dimensions)
  // ------------------------------------------------------------------
  println("=== 2. Assay overview ===")
  println("n_assays: \{qf.n_assays()}")
  let names = qf.assay_names()
  let mut i = 0
  while i < names.length() {
    match qf.get_assay(names[i]) {
      Some(a) =>
        println(
          "  \{names[i]}: \{a.n_rows()} features x \{a.n_cols()} samples",
        )
      None => println("  \{names[i]}: not found")
    }
    i = i + 1
  }
  println("Total links: \{qf.links().n_links()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show PSM assay details (first few rows)
  // ------------------------------------------------------------------
  println("=== 3. PSM assay details ===")
  let psm = qf.get_assay("psm").unwrap()
  println(@bio.qf_assay_summary(psm))
  let col_names = psm.col_names()
  let header = "  feature\t" + col_names.join("\t")
  println(header)
  let n_show = if psm.n_rows() < 4 { psm.n_rows() } else { 4 }
  let mut r = 0
  while r < n_show {
    let mut line = "  " + psm.row_names()[r]
    let mut c = 0
    while c < psm.n_cols() {
      line = line + "\t" + psm.get(r, c).to_string()
      c = c + 1
    }
    println(line)
    r = r + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Show links between PSM and peptide levels
  // ------------------------------------------------------------------
  println("=== 4. Links between PSM and peptide levels ===")
  let mut p = 0
  while p < psm.n_rows() {
    let links = qf.links().get_links_from("psm", p)
    let mut l = 0
    while l < links.length() {
      let link = links[l]
      let peptide = qf.get_assay("peptide").unwrap()
      let pep_name = peptide.row_names()[link.to_row()]
      println(
        "  PSM \{psm.row_names()[p]} (row \{p}) -> \{link.to_assay()} row \{link.to_row()} (\{pep_name})",
      )
      l = l + 1
    }
    p = p + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Apply log2 transform to PSM assay, show transformed values
  // ------------------------------------------------------------------
  println("=== 5. log2 transform of PSM assay ===")
  let psm_log = @bio.qf_log_transform(psm, 2.0, 1.0)
  println("Before (row 0): " + row_to_string(psm, 0))
  println("After  (row 0): " + row_to_string(psm_log, 0))
  println("")

  // ------------------------------------------------------------------
  // 6. Filter low abundance features
  // ------------------------------------------------------------------
  println("=== 6. Filter low abundance PSM features ===")
  let qf_filtered = @bio.qf_filter_low_abundance(qf, "psm", 15.0)
  let psm_filtered = qf_filtered.get_assay("psm").unwrap()
  println(
    "PSM features: \{psm.n_rows()} -> \{psm_filtered.n_rows()} (threshold=15.0)",
  )
  println("Remaining PSM rows:")
  let mut f = 0
  while f < psm_filtered.n_rows() {
    println("  \{psm_filtered.row_names()[f]}")
    f = f + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Normalize using center/scale
  // ------------------------------------------------------------------
  println("=== 7. Center/scale normalization (peptide assay) ===")
  let peptide = qf.get_assay("peptide").unwrap()
  let peptide_norm = @bio.qf_normalize_center_scale(peptide, true, true)
  println("Peptide before normalization (row 0): " + row_to_string(peptide, 0))
  println("Peptide after  normalization (row 0): " + row_to_string(peptide_norm, 0))
  println("")

  // ------------------------------------------------------------------
  // 8. Show final protein-level assay values
  // ------------------------------------------------------------------
  println("=== 8. Protein-level assay ===")
  let protein = qf.get_assay("protein").unwrap()
  println(@bio.qf_assay_summary(protein))
  let header2 = "  protein\t" + protein.col_names().join("\t")
  println(header2)
  let mut pr = 0
  while pr < protein.n_rows() {
    let mut line = "  " + protein.row_names()[pr]
    let mut pc = 0
    while pc < protein.n_cols() {
      line = line + "\t" + protein.get(pr, pc).to_string()
      pc = pc + 1
    }
    println(line)
    pr = pr + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 9. Convert protein assay to long format and show first few entries
  // ------------------------------------------------------------------
  println("=== 9. Protein assay in long format ===")
  let long = @bio.qf_to_long_format(qf, "protein")
  println("Total long-format entries: \{long.length()}")
  let n_long = if long.length() < 6 { long.length() } else { 6 }
  println("  protein\tsample\tvalue")
  let mut g = 0
  while g < n_long {
    let (feat, sample, value) = long[g]
    println("  \{feat}\t\{sample}\t\{value}")
    g = g + 1
  }

  println("=== Demo complete ===")
}

/// Format a single assay row as a tab-separated string of values.
fn row_to_string(assay : @bio.QfAssay, row : Int) -> String {
  let mut s = ""
  let mut c = 0
  while c < assay.n_cols() {
    s = s + assay.get(row, c).to_string()
    if c < assay.n_cols() - 1 {
      s = s + ", "
    }
    c = c + 1
  }
  s
}
