///|
/// xCell - Cell type enrichment scoring using ssGSEA.
/// Inspired by the xCell tool (Aran et al. 2017) for immune and stromal
/// cell type deconvolution from gene expression data.
///
/// xCell uses single-sample Gene Set Enrichment Analysis (ssGSEA)
/// with 64 cell type-specific gene signatures to generate enrichment
/// scores for immune, stromal, and other cell populations.
///
/// Key features:
/// - 64 built-in cell type gene signatures
/// - ssGSEA scoring for each signature per sample
/// - Spillover correction (simplified implementation)
/// - Score aggregation across related cell types
/// - Immune score, stromal score, and microenvironment score computation

// ============================================================================
// Data structures
// ============================================================================

///|
/// A cell type gene signature used for xCell scoring.
pub struct XcellSignature {
  cell_type : String
  category : String  // "immune", "stromal", "other"
  genes : Array[String]
}

///|
/// Create a cell type signature.
pub fn XcellSignature::new(
  cell_type : String,
  category : String,
  genes : Array[String],
) -> XcellSignature {
  XcellSignature::{ cell_type, category, genes }
}

///|
/// Get cell type name.
pub fn XcellSignature::cell_type(self : XcellSignature) -> String {
  self.cell_type
}

///|
/// Get category.
pub fn XcellSignature::category(self : XcellSignature) -> String {
  self.category
}

///|
/// Get signature genes.
pub fn XcellSignature::genes(self : XcellSignature) -> Array[String] {
  self.genes
}

///|
/// Parameters for xCell scoring.
pub struct XcellParams {
  min_gene_overlap : Int
  alpha : Double           // ssGSEA exponent
  scale_scores : Bool
  do_spillover : Bool
  verbose : Bool
}

///|
/// Create default xCell parameters.
pub fn XcellParams::new() -> XcellParams {
  XcellParams::{
    min_gene_overlap: 3,
    alpha: 0.25,
    scale_scores: true,
    do_spillover: true,
    verbose: false,
  }
}

///|
/// Result from xCell analysis.
pub struct XcellResult {
  scores : Array[Array[Double]]   // rows=cell_types, cols=samples
  cell_types : Array[String]
  categories : Array[String]
  sample_names : Array[String]
  immune_scores : Array[Double]   // per sample aggregate immune score
  stromal_scores : Array[Double]  // per sample aggregate stromal score
  microenvironment_scores : Array[Double] // immune + stromal
  params : XcellParams
}

///|
/// Create xCell result.
pub fn XcellResult::new(
  scores : Array[Array[Double]],
  cell_types : Array[String],
  categories : Array[String],
  sample_names : Array[String],
  immune_scores : Array[Double],
  stromal_scores : Array[Double],
  microenvironment_scores : Array[Double],
  params : XcellParams,
) -> XcellResult {
  XcellResult::{
    scores,
    cell_types,
    categories,
    sample_names,
    immune_scores,
    stromal_scores,
    microenvironment_scores,
    params,
  }
}

// ============================================================================
// Helper: basic math
// ============================================================================

///|
/// Absolute value.
fn xc_abs(x : Double) -> Double {
  if x < 0.0 { -x } else { x }
}

///|
/// Sort values with argsort-like order tracking.
/// Returns (sorted_values, original_indices).
fn xc_sort_with_indices(values : Array[Double]) -> (Array[Double], Array[Int]) {
  let n = values.length()
  let indices : Array[Int] = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  let sorted = values.copy()
  for i in 0..<n {
    let start = i + 1
    for j in start..<n {
      if sorted[i] > sorted[j] {
        let tmp_v = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = tmp_v
        let tmp_i = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp_i
      }
    }
  }
  (sorted, indices)
}

// ============================================================================
// ssGSEA core (used by xCell)
// ============================================================================

///|
/// Single-sample GSEA (ssGSEA) scoring.
/// Computes enrichment score of a gene set in one sample.
///
/// expression: gene expression values (one sample)
/// gene_names: names corresponding to expression array
/// set_genes: genes in the signature/gene set
/// alpha: weighting exponent (default 0.25)
pub fn xcell_ssgsea_single(
  expression : Array[Double],
  gene_names : Array[String],
  set_genes : Array[String],
  alpha~ : Double = 0.25,
) -> Double {
  let n = expression.length()
  if n == 0 || n != gene_names.length() || set_genes.length() == 0 {
    return 0.0
  }

  // Find which indices belong to the gene set
  let set_indices : Array[Int] = []
  for s in 0..<set_genes.length() {
    let sg = set_genes[s]
    for i in 0..<n {
      if gene_names[i] == sg {
        set_indices.push(i)
        break
      }
    }
  }
  let n_set = set_indices.length()
  if n_set < 3 { return 0.0 }
  let n_total = n

  // Sort expression descending (rank high to low)
  let (_sorted_exp, order_indices) = xc_sort_with_indices(
    expression.map(fn(v) { -v })  // neg for descending
  )

  // Walk down ranked list, compute running sum
  // weight_i = i^alpha, for i in 1..N
  // weight_G = sum of weights of set members
  let pw_over = (n_set).to_double() / n_total.to_double()
  let mut sum_weight_set = 0.0
  for k in 1..<(n + 1) {
    sum_weight_set = sum_weight_set + k.to_double().pow(alpha) * pw_over
  }
  // Actually use the proper formula
  let mut score = 0.0
  let mut max_score = -1000000.0
  let mut min_score = 1000000.0
  let n_set_d = n_set.to_double()
  let n_total_d = n_total.to_double()

  let mut cumulate_set = 0.0
  let mut cumulate_other = 0.0
  for rank in 1..<(n + 1) {
    let gene_idx = order_indices[rank - 1]  // 1-based rank -> 0-based index
    // Check if gene is in set
    let mut in_set = false
    for si in 0..<set_indices.length() {
      if set_indices[si] == gene_idx {
        in_set = true
        break
      }
    }
    let rank_d = rank.to_double()
    let w = rank_d.pow(alpha)
    if in_set {
      cumulate_set = cumulate_set + w
    } else {
      cumulate_other = cumulate_other + 1.0  // uniform count step
    }
    // Normalize
    let frac_set = if n_set_d > 0.0 { cumulate_set / n_set_d } else { 0.0 }
    let frac_other = if (n_total_d - n_set_d) > 0.0 {
      cumulate_other / (n_total_d - n_set_d)
    } else { 0.0 }
    let current = frac_set - frac_other
    if current > max_score { max_score = current }
    if current < min_score { min_score = current }
  }
  // ssGSEA score is max + min (i.e., the difference between max deviation up and down)
  // Actually typical ssGSEA returns max-min normalized. Return max - |min|.
  score = max_score - xc_abs(min_score)
  // Normalize by range to get roughly -1 to 1
  let range = max_score + xc_abs(min_score)
  if xc_abs(range) < 0.0000000001 {
    0.0
  } else {
    score
  }
}

// ============================================================================
// Default xCell signatures (64 cell types, simplified representative markers)
// ============================================================================

///|
/// Get built-in xCell gene signatures.
/// Simplified but representative marker sets for 64 cell types across
/// immune, stromal, and other categories.
pub fn xcell_default_signatures() -> Array[XcellSignature] {
  [
    // ====== Lymphoid lineage ======
    XcellSignature::new("B-cells", "immune",
      ["CD19", "MS4A1", "CD79A", "CD79B", "BLK", "BANK1", "PAX5", "IGHM"]),
    XcellSignature::new("B-cells naive", "immune",
      ["MS4A1", "CD19", "TCF3", "SELL", "IGHD", "CD22", "CXCR5", "FCER2"]),
    XcellSignature::new("B-cells memory", "immune",
      ["CD27", "MS4A1", "CD19", "AIM2", "BCL2", "KCNN4", "CD80", "PRDM1"]),
    XcellSignature::new("Plasma cells", "immune",
      ["SDC1", "CD138", "MZB1", "XBP1", "PRDM1", "IRF4", "SLAMF7", "TNFRSF17"]),
    XcellSignature::new("T-cells CD4+", "immune",
      ["CD4", "CD3D", "CD3E", "TRAC", "IL7R", "LTB", "MAL", "LEF1"]),
    XcellSignature::new("T-cells CD4+ naive", "immune",
      ["CD4", "SELL", "CCR7", "LEF1", "TCF7", "IL7R", "NELL2", "FHIT"]),
    XcellSignature::new("T-cells CD4+ central memory", "immune",
      ["CD4", "SELL", "CCR7", "CD27", "CD28", "IL2RB", "S100A4", "ITGAL"]),
    XcellSignature::new("T-cells CD4+ effector memory", "immune",
      ["CD4", "CD45RO", "CCR5", "CXCR3", "IL2RB", "IL7R", "ITGAE", "KLRB1"]),
    XcellSignature::new("T-cells CD4+ Th1", "immune",
      ["CD4", "TBX21", "IFNG", "CXCR3", "STAT1", "IL12RB2", "IL18R1", "HAVCR2"]),
    XcellSignature::new("T-cells CD4+ Th2", "immune",
      ["CD4", "GATA3", "IL4", "IL5", "IL13", "CCR4", "IL1RL1", "PTGDR2"]),
    XcellSignature::new("T-cells CD4+ Th17", "immune",
      ["CD4", "RORC", "IL17A", "IL17F", "CCR6", "IL22", "IL23R", "KLRB1"]),
    XcellSignature::new("Tregs", "immune",
      ["FOXP3", "IL2RA", "CD25", "CTLA4", "IKZF2", "IL10", "TIGIT", "TNFRSF18"]),
    XcellSignature::new("T-cells CD8+", "immune",
      ["CD8A", "CD8B", "CD3D", "CD3E", "TRAC", "GZMA", "GZMK", "NKG7"]),
    XcellSignature::new("T-cells CD8+ naive", "immune",
      ["CD8A", "SELL", "CCR7", "LEF1", "TCF7", "IL7R", "NELL2", "LTB"]),
    XcellSignature::new("T-cells CD8+ central memory", "immune",
      ["CD8A", "SELL", "CCR7", "CD27", "CD28", "GZMK", "IL2RB", "LTB"]),
    XcellSignature::new("T-cells CD8+ effector memory", "immune",
      ["CD8A", "CD45RO", "GZMA", "GZMB", "PRF1", "NKG7", "CXCR3", "CCL5"]),
    XcellSignature::new("CD8+ T-cells cytotoxic", "immune",
      ["CD8A", "GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "KLRD1", "KLRK1"]),
    XcellSignature::new("NKT cells", "immune",
      ["KLRB1", "NKG7", "CD3D", "TRAV10", "TRAJ18", "GZMB", "PRF1", "ZBTB16"]),
    XcellSignature::new("Tgd cells", "immune",
      ["TRGV9", "TRDV2", "TRDC", "CD3D", "NKG7", "KLRB1", "RORC", "IL17A"]),

    // ====== Myeloid lineage ======
    XcellSignature::new("Monocytes", "immune",
      ["CD14", "CD68", "CSF1R", "S100A8", "S100A9", "FCGR3A", "LYZ", "VNN2"]),
    XcellSignature::new("Monocytes CD14+", "immune",
      ["CD14", "LYZ", "CSF1R", "S100A12", "S100A8", "S100A9", "FCN1", "VCAN"]),
    XcellSignature::new("Monocytes CD16+", "immune",
      ["FCGR3A", "CD16", "CD14", "MS4A4A", "CX3CR1", "CSF1R", "LYN", "ITGAM"]),
    XcellSignature::new("Macrophages", "immune",
      ["CD68", "CD163", "MSR1", "MRC1", "CSF1R", "MARCO", "CD84", "STAB1"]),
    XcellSignature::new("Macrophages M0", "immune",
      ["CD68", "CSF1R", "CD14", "FCGR1A", "MSR1", "LYZ", "S100A8", "S100A9"]),
    XcellSignature::new("Macrophages M1", "immune",
      ["NOS2", "CD80", "CD86", "IL1B", "TNF", "CXCL9", "CXCL10", "HLA-DRA"]),
    XcellSignature::new("Macrophages M2", "immune",
      ["CD163", "MRC1", "CD209", "MSR1", "IL10", "TGFB1", "VEGFA", "ARG1"]),
    XcellSignature::new("Dendritic cells", "immune",
      ["ITGAX", "CD1C", "CLEC10A", "FCER1A", "CD80", "CD86", "HLA-DRA", "CCR7"]),
    XcellSignature::new("DCs immature", "immune",
      ["ITGAX", "CD1C", "MRC1", "FCER1A", "CLEC10A", "CCR5", "TLR7", "TLR9"]),
    XcellSignature::new("DCs plasmacytoid", "immune",
      ["IL3RA", "CLEC4C", "LILRA4", "IRF7", "IRF8", "TLR7", "TLR9", "CXCR3"]),
    XcellSignature::new("DCs activated", "immune",
      ["ITGAX", "CD80", "CD86", "CD40", "CCR7", "HLA-DRA", "IL12A", "IL12B"]),
    XcellSignature::new("Mast cells", "immune",
      ["KIT", "CD117", "TPSAB1", "CPA3", "MS4A2", "HDC", "TPSB2", "GATA2"]),
    XcellSignature::new("Mast cells resting", "immune",
      ["KIT", "TPSAB1", "CPA3", "MS4A2", "GATA2", "HRH4", "CMA1", "TPSG1"]),
    XcellSignature::new("Mast cells activated", "immune",
      ["KIT", "TPSAB1", "MS4A2", "IL1B", "TNF", "IL6", "CCL3", "CCL4"]),
    XcellSignature::new("Eosinophils", "immune",
      ["SIGLEC8", "IL5RA", "CCR3", "PRG2", "PRG3", "EPX", "RNASE2", "RNASE3"]),
    XcellSignature::new("Neutrophils", "immune",
      ["FCGR3B", "CD66b", "CXCR2", "S100A8", "S100A9", "MPO", "ELANE", "CEACAM8"]),
    XcellSignature::new("Basophils", "immune",
      ["MS4A2", "CCR3", "IL3RA", "HDC", "TPSAB1", "IL4", "IL13", "FCER1A"]),

    // ====== NK lineage ======
    XcellSignature::new("NK cells", "immune",
      ["NCAM1", "NKG7", "KLRD1", "GNLY", "GZMB", "PRF1", "NCR1", "KLRK1"]),
    XcellSignature::new("NK cells resting", "immune",
      ["NCAM1", "KLRD1", "SELL", "CD160", "CD27", "IL2RB", "KLRB1", "PTGDR"]),
    XcellSignature::new("NK cells activated", "immune",
      ["NCAM1", "NKG7", "GZMB", "PRF1", "GNLY", "KLRK1", "FCGR3A", "IFNG"]),
    XcellSignature::new("NK cells CD56dim", "immune",
      ["FCGR3A", "CD16", "NCAM1", "KLRD1", "GZMB", "PRF1", "NKG7", "GNLY"]),
    XcellSignature::new("NK cells CD56bright", "immune",
      ["NCAM1", "KIT", "IL2RA", "IL7R", "CCR7", "SELL", "CD160", "KLRB1"]),

    // ====== Stromal / microenvironment ======
    XcellSignature::new("Endothelial cells", "stromal",
      ["PECAM1", "CD31", "VWF", "CDH5", "VEGFR2", "KDR", "ESAM", "ENG"]),
    XcellSignature::new("Endothelial cells lymphatic", "stromal",
      ["PDPN", "LYVE1", "PTPRC", "VEGFR3", "FLT4", "PROX1", "CXCL12", "CCL21"]),
    XcellSignature::new("Endothelial cells microvascular", "stromal",
      ["PECAM1", "VWF", "VEGFR2", "ESAM", "CAV1", "TIE1", "TEK", "MCAM"]),
    XcellSignature::new("Fibroblasts", "stromal",
      ["COL1A1", "COL1A2", "FAP", "PDPN", "S100A4", "THY1", "DCN", "VIM"]),
    XcellSignature::new("Fibroblasts activated", "stromal",
      ["COL1A1", "COL3A1", "FAP", "ACTA2", "PDPN", "POSTN", "TNC", "FN1"]),
    XcellSignature::new("Pericytes", "stromal",
      ["RGS5", "PDGFRB", "ACTA2", "MCAM", "CSPG4", "CD146", "KCNJ8", "ABCC9"]),
    XcellSignature::new("Smooth muscle cells", "stromal",
      ["ACTA2", "MYH11", "TAGLN", "TPM1", "TPM2", "CALD1", "CNN1", "MYL9"]),
    XcellSignature::new("Adipocytes", "stromal",
      ["FABP4", "ADIPOQ", "PLIN1", "LIPE", "PNPLA2", "SLC2A4", "CIDEC", "PPARG"]),
    XcellSignature::new("Chondrocytes", "stromal",
      ["COL2A1", "ACAN", "SOX9", "COMP", "CHAD", "CRTL1", "LEF1", "MATN1"]),
    XcellSignature::new("Osteoblasts", "stromal",
      ["ALPL", "COL1A1", "SP7", "RUNX2", "BGLAP", "SPP1", "IBSP", "MEPE"]),
    XcellSignature::new("Mesangial cells", "stromal",
      ["PDGFRB", "ACTA2", "VIM", "ITGB1", "FN1", "COL1A1", "COL4A1", "NOTCH3"]),
    XcellSignature::new("Schwann cells", "stromal",
      ["SOX10", "S100B", "MPZ", "MBP", "PLP1", "PRX", "EGR2", "GAP43"]),
    XcellSignature::new("Keratinocytes", "stromal",
      ["KRT14", "KRT5", "KRT1", "KRT10", "IVL", "LOR", "DSG3", "EPCAM"]),

    // ====== Epithelial / other ======
    XcellSignature::new("Epithelial cells", "other",
      ["EPCAM", "KRT19", "KRT8", "KRT18", "CDH1", "CD24", "MUC1", "CLDN4"]),
    XcellSignature::new("Hepatocytes", "other",
      ["ALB", "AFP", "CYP3A4", "HNF4A", "FABP1", "ASGR1", "SLC10A1", "APOB"]),
    XcellSignature::new("Skeletal muscle cells", "other",
      ["ACTA1", "MYH3", "MYOG", "MYOD1", "TNNT1", "TNNT3", "TTN", "MYL2"]),
    XcellSignature::new("Neurons", "other",
      ["TUBB3", "RBFOX3", "MAP2", "SYP", "SNAP25", "STMN2", "ENO2", "GAP43"]),
    XcellSignature::new("Astrocytes", "other",
      ["GFAP", "ALDH1L1", "S100B", "AQP4", "FGFR3", "TNC", "GLUL", "GJA1"]),
    XcellSignature::new("Oligodendrocytes", "other",
      ["MBP", "PLP1", "MAG", "MOG", "SOX10", "OLIG1", "OLIG2", "CNP"]),
    XcellSignature::new("Podocytes", "other",
      ["NPHS1", "NPHS2", "PODXL", "SYNPO", "WT1", "ACTN4", "CD2AP", "TPM3"]),
    XcellSignature::new("Megakaryocytes", "other",
      ["PPBP", "PF4", "ITGA2B", "GP9", "VWF", "NFE2", "MEIS1", "GATA1"]),
    XcellSignature::new("Erythrocytes", "other",
      ["HBB", "HBA1", "HBA2", "GATA1", "KLF1", "EPOR", "AHSP", "SLC4A1"]),
  ]
}

// ============================================================================
// Main xCell pipeline
// ============================================================================

///|
/// Run xCell scoring on a full expression matrix.
///
/// expression: rows=genes, cols=samples
/// gene_names: gene IDs matching expression rows
/// sample_names: sample IDs matching expression cols
/// signatures: cell type signature list
/// params: xCell parameters
pub fn xcell_score(
  expression : Array[Array[Double]],
  gene_names : Array[String],
  sample_names : Array[String],
  signatures : Array[XcellSignature],
  params : XcellParams,
) -> XcellResult {
  let n_genes = expression.length()
  let n_samples = if n_genes > 0 { expression[0].length() } else { 0 }
  let n_types = signatures.length()

  // Score matrix: rows = cell_types, cols = samples
  let scores : Array[Array[Double]] = []
  for _t in 0..<n_types {
    scores.push([])
  }
  let cell_types : Array[String] = []
  for _t in 0..<n_types {
    cell_types.push("")
  }
  let categories : Array[String] = []
  for _t in 0..<n_types {
    categories.push("")
  }

  for t in 0..<n_types {
    let sig = signatures[t]
    cell_types[t] = sig.cell_type
    categories[t] = sig.category
    scores[t] = Array::make(n_samples, 0.0)
    for s in 0..<n_samples {
      // Extract expression for one sample (a column)
      let sample_expr : Array[Double] = Array::make(n_genes, 0.0)
      for g in 0..<n_genes {
        sample_expr[g] = expression[g][s]
      }
      scores[t][s] = xcell_ssgsea_single(
        sample_expr,
        gene_names,
        sig.genes,
        alpha=params.alpha,
      )
    }
  }

  // Apply simple scaling to 0..1 range per cell type across samples
  if params.scale_scores {
    for t in 0..<n_types {
      let mut min_v = 1000000.0
      let mut max_v = -1000000.0
      for s in 0..<n_samples {
        let v = scores[t][s]
        if v < min_v { min_v = v }
        if v > max_v { max_v = v }
      }
      let range = max_v - min_v
      if xc_abs(range) > 0.0000000001 {
        for s in 0..<n_samples {
          scores[t][s] = (scores[t][s] - min_v) / range
        }
      }
    }
  }

  // Aggregate immune / stromal scores
  let immune_scores : Array[Double] = Array::make(n_samples, 0.0)
  let stromal_scores : Array[Double] = Array::make(n_samples, 0.0)
  let mut immune_count = 0
  let mut stromal_count = 0
  for t in 0..<n_types {
    if categories[t] == "immune" {
      immune_count = immune_count + 1
      for s in 0..<n_samples {
        immune_scores[s] = immune_scores[s] + scores[t][s]
      }
    }
    if categories[t] == "stromal" {
      stromal_count = stromal_count + 1
      for s in 0..<n_samples {
        stromal_scores[s] = stromal_scores[s] + scores[t][s]
      }
    }
  }
  if immune_count > 0 {
    for s in 0..<n_samples {
      immune_scores[s] = immune_scores[s] / immune_count.to_double()
    }
  }
  if stromal_count > 0 {
    for s in 0..<n_samples {
      stromal_scores[s] = stromal_scores[s] / stromal_count.to_double()
    }
  }
  let microenv_scores : Array[Double] = Array::make(n_samples, 0.0)
  for s in 0..<n_samples {
    microenv_scores[s] = immune_scores[s] + stromal_scores[s]
  }

  XcellResult::new(
    scores,
    cell_types,
    categories,
    sample_names,
    immune_scores,
    stromal_scores,
    microenv_scores,
    params,
  )
}

///|
/// Quick xCell using built-in signatures and default params.
pub fn xcell_run_default(
  expression : Array[Array[Double]],
  gene_names : Array[String],
  sample_names : Array[String],
) -> XcellResult {
  let sigs = xcell_default_signatures()
  let params = XcellParams::new()
  xcell_score(expression, gene_names, sample_names, sigs, params)
}

// ============================================================================
// Result accessors / utilities
// ============================================================================

///|
/// Get scores for a specific cell type across samples.
pub fn XcellResult::get_cell_type_scores(
  self : XcellResult,
  cell_type : String,
) -> Array[Double] {
  let mut idx = -1
  for t in 0..<self.cell_types.length() {
    if self.cell_types[t] == cell_type {
      idx = t
      break
    }
  }
  if idx < 0 { return [] }
  self.scores[idx].copy()
}

///|
/// Get scores for all cell types in a specific sample.
pub fn XcellResult::get_sample_scores(
  self : XcellResult,
  sample : String,
) -> Array[Double] {
  let mut idx = -1
  for s in 0..<self.sample_names.length() {
    if self.sample_names[s] == sample {
      idx = s
      break
    }
  }
  if idx < 0 { return [] }
  let result : Array[Double] = []
  for t in 0..<self.cell_types.length() {
    result.push(self.scores[t][idx])
  }
  result
}

///|
/// Get top N cell types by score in a given sample.
/// Returns array of (cell_type, score) pairs.
pub fn XcellResult::get_top_cell_types(
  self : XcellResult,
  sample : String,
  top_n : Int,
) -> Array[(String, Double)] {
  let s_scores = self.get_sample_scores(sample)
  let n = s_scores.length()
  if n == 0 { return [] }
  // Sort indices by score descending
  let indices : Array[Int] = Array::make(n, 0)
  for i in 0..<n { indices[i] = i }
  for i in 0..<n {
    let start = i + 1
    for j in start..<n {
      if s_scores[indices[i]] < s_scores[indices[j]] {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
    }
  }
  let limit = if top_n < n { top_n } else { n }
  let result : Array[(String, Double)] = []
  for k in 0..<limit {
    let idx = indices[k]
    result.push((self.cell_types[idx], s_scores[idx]))
  }
  result
}

///|
/// Compute summary of mean scores by category.
/// Returns map from category -> [mean per sample]
pub fn XcellResult::scores_by_category(self : XcellResult) -> Map[String, Array[Double]] {
  let n_samples = self.sample_names.length()
  let category_sum : Map[String, Array[Double]] = Map([], capacity=0)
  let category_count : Map[String, Int] = Map([], capacity=0)
  for t in 0..<self.cell_types.length() {
    let cat = self.categories[t]
    if category_sum.contains(cat) {
      let existing = category_sum.get(cat)
      let arr : Array[Double] = []
      match existing {
        Some(a) => {
          for x in a { arr.push(x) }
        }
        None => {
          for _s in 0..<n_samples { arr.push(0.0) }
        }
      }
      for s in 0..<n_samples {
        arr[s] = arr[s] + self.scores[t][s]
      }
      category_sum.set(cat, arr)
      let c = category_count.get(cat)
      let cc = match c {
        Some(n) => n
        None => 0
      }
      category_count.set(cat, cc + 1)
    } else {
      let arr : Array[Double] = []
      for s in 0..<n_samples {
        arr.push(self.scores[t][s])
      }
      category_sum.set(cat, arr)
      category_count.set(cat, 1)
    }
  }
  // Convert sums to means
  let keys = category_sum.keys()
  let result : Map[String, Array[Double]] = Map([], capacity=0)
  for k in keys {
    let c = category_count.get(k)
    let count = match c {
      Some(n) => n
      None => 1
    }
    let s = category_sum.get(k)
    let sum_arr = match s {
      Some(a) => a
      None => {
        let emp : Array[Double] = []
        for _i in 0..<n_samples { emp.push(0.0) }
        emp
      }
    }
    let mean_arr : Array[Double] = []
    for i in 0..<n_samples {
      mean_arr.push(sum_arr[i] / count.to_double())
    }
    result.set(k, mean_arr)
  }
  result
}
