///|
/// infercnv module example - Copy Number Variation inference from scRNA-seq.
///
/// Run: moon run examples/infercnv_demo
fn main {
  println("=== infercnv Module Demo ===")
  println("Copy Number Variation Inference for Single-Cell RNA-seq")
  println("Reference: Tirosh I et al. (2016) Science, Broad Institute inferCNV")
  println("")

  // ------------------------------------------------------------------
  // 1. Generate synthetic single-cell data
  // ------------------------------------------------------------------
  println("=== 1. Sample synthetic tumour / normal dataset ===")
  let n_tumour = 30
  let n_normal = 20
  let n_chr = 5
  let n_genes_per_chr = 40
  let (input, tumour_cats, normal_cats) = @bio.infercnv_sample_data(
    n_tumour=n_tumour, n_normal=n_normal, n_chr=n_chr,
    n_genes_per_chr=n_genes_per_chr, seed=2025,
  )
  println("Cells: \{input.n_cells}, Genes: \{input.n_genes}, Chromosomes: \{n_chr}")
  println("Tumour categories: \{tumour_cats}")
  println("Normal categories: \{normal_cats}")
  println("")

  // ------------------------------------------------------------------
  // 2. Method 1: global mean as reference (no prior knowledge)
  // ------------------------------------------------------------------
  println("=== 2. Method 1: Global-mean reference ===")
  let params_default = @bio.make_cnv_params(
    21, 1.5, 0.05, @bio.ref_method_global_mean(), [],
  )
  let res_global = @bio.run_infercnv(
    input, params=params_default, excluded_chromosomes=[],
  )
  println("\{res_global.summary()}")
  let mut iter = 0
  while iter < res_global.per_cluster_cnv_score.length() {
    let (cat, score) = res_global.per_cluster_cnv_score[iter]
    println("  Cluster '\{cat}' mean CNV score = \{score}")
    iter = iter + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Method 2: known normal cell types as reference (recommended)
  // ------------------------------------------------------------------
  println("=== 3. Method 2: Reference-categories (use normal cells) ===")
  let params_cat = @bio.make_cnv_params(
    21, 1.5, 0.05, @bio.ref_method_reference_categories(), normal_cats,
  )
  let res_cat = @bio.run_infercnv(
    input, params=params_cat, excluded_chromosomes=["chrX", "chrY", "chrM"],
  )
  println("\{res_cat.summary()}")
  let t_score = res_cat.cluster_score("Tumour")
  let n_score = res_cat.cluster_score("Normal")
  println("  Tumour score: \{t_score}")
  println("  Normal score: \{n_score}")
  println("  Fold difference: \{t_score / n_score}x")
  println("  Chromosomes kept after exclusion: \{res_cat.chromosome_list().length()}")
  println("")

  // ------------------------------------------------------------------
  // 4. Per-chromosome / boundary overview
  // ------------------------------------------------------------------
  println("=== 4. Chromosome boundaries (ordered genes) ===")
  let mut k = 0
  while k < res_cat.chromosome_list().length() {
    let (chr, s, e) = res_cat.chromosome_list()[k]
    println("  \{chr}: genes index \{s} .. \{e} (span = \{e - s + 1})")
    k = k + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Tumour / normal cell prediction
  // ------------------------------------------------------------------
  println("=== 5. Predict malignant cells from CNV score ===")
  let flagged = res_cat.predict_tumour_cells("Normal", threshold_factor=1.3)
  println("  Cells flagged as putatively tumour-like: \{flagged.length()}/\{res_cat.n_cells}")
  // Ground truth: first n_tumour indices are tumour
  let mut true_pos = 0
  let mut j = 0
  while j < flagged.length() {
    if flagged[j] < n_tumour { true_pos = true_pos + 1 }
    j = j + 1
  }
  let false_pos = flagged.length() - true_pos
  println("  Tumour cells correctly flagged: \{true_pos}/\{n_tumour} (TP)")
  println("  Normal cells mis-flagged: \{false_pos}/\{n_normal} (FP)")
  println("")

  // ------------------------------------------------------------------
  // 6. Spot-check CNV matrix around simulated deletion (chr1 first half)
  // ------------------------------------------------------------------
  println("=== 6. Spot-check simulated CNV signal ===")
  // chr1 is chromosome_list()[0]; second half of chr2 is copy gain
  let boundaries = res_cat.chromosome_list()
  if boundaries.length() >= 2 {
    let (chr1, s1, e1) = boundaries[0]
    let (chr2, s2, e2) = boundaries[1]
    let mid1 = (s1 + e1) / 2          // midpoint of chr1 (deletion side)
    let mid2 = (s2 + e2) / 2 + 5      // latter part of chr2 (amplification)
    // Average a few tumour cells vs normal cells
    let mut tumour_del = 0.0
    let mut tumour_amp = 0.0
    let mut normal_del = 0.0
    let mut normal_amp = 0.0
    let n_sample = 5
    let mut s = 0
    while s < n_sample {
      tumour_del = tumour_del + res_cat.at(s, mid1)
      tumour_amp = tumour_amp + res_cat.at(s, mid2)
      normal_del = normal_del + res_cat.at(n_tumour + s, mid1)
      normal_amp = normal_amp + res_cat.at(n_tumour + s, mid2)
      s = s + 1
    }
    let ns = n_sample.to_double()
    println("  Tumour avg chr1-first-half (deletion site) CNV = \{tumour_del / ns}")
    println("  Normal avg chr1-first-half CNV = \{normal_del / ns}")
    println("  Tumour avg chr2-second-half (amplification site) CNV = \{tumour_amp / ns}")
    println("  Normal avg chr2-second-half CNV = \{normal_amp / ns}")
    let _ = chr1
    let _ = chr2
  }
  println("")

  println("=== Demo complete ===")
}
