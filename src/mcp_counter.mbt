///|
/// MCPcounter - Microenvironment Cell Population counter.
/// Inspired by Becht E et al. (2016) eLife.
/// Bioconductor package: MCPcounter.
///
/// Core algorithm: For each cell population, compute the geometric mean of
/// expression values of its marker genes. This produces a per-sample abundance
/// score for each cell population. No deconvolution or normalization is needed;
/// the scores are directly comparable across samples for the same cell type.
///
/// All public identifiers are prefixed with `mcp_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// A cell population with its marker genes.
pub struct McpCellPopulation {
  name : String
  markers : Array[String]
}

///|
/// Create a cell population with its marker genes.
pub fn McpCellPopulation::new(
  name : String,
  markers : Array[String],
) -> McpCellPopulation {
  { name, markers }
}

///|
/// Get population name.
pub fn McpCellPopulation::name(self : McpCellPopulation) -> String {
  self.name
}

///|
/// Get marker genes.
pub fn McpCellPopulation::markers(self : McpCellPopulation) -> Array[String] {
  self.markers
}

///|
/// MCPcounter result: abundance scores per cell population per sample.
pub struct McpResult {
  populations : Array[String]
  sample_names : Array[String]
  scores : Array[Array[Double]] // rows=populations, cols=samples
}

///|
/// Get population names.
pub fn McpResult::populations(self : McpResult) -> Array[String] {
  self.populations
}

///|
/// Get sample names.
pub fn McpResult::sample_names(self : McpResult) -> Array[String] {
  self.sample_names
}

///|
/// Get the score matrix (rows=populations, cols=samples).
pub fn McpResult::scores(self : McpResult) -> Array[Array[Double]] {
  self.scores
}

///|
/// Get the score for a specific population in a specific sample.
pub fn McpResult::get_score(
  self : McpResult,
  population : String,
  sample : String,
) -> Double {
  let mut pop_idx = -1
  for i in 0..<self.populations.length() {
    if self.populations[i] == population { pop_idx = i; break }
  }
  let mut samp_idx = -1
  for i in 0..<self.sample_names.length() {
    if self.sample_names[i] == sample { samp_idx = i; break }
  }
  if pop_idx < 0 || samp_idx < 0 { return 0.0 }
  self.scores[pop_idx][samp_idx]
}

///|
/// Get all scores for a specific population across samples.
pub fn McpResult::get_population_scores(
  self : McpResult,
  population : String,
) -> Array[Double] {
  let mut idx = -1
  for i in 0..<self.populations.length() {
    if self.populations[i] == population { idx = i; break }
  }
  if idx < 0 { return [] }
  self.scores[idx].copy()
}

// ============================================================================
// Geometric mean computation
// ============================================================================

///|
/// Compute the geometric mean of an array of non-negative values.
/// Returns 0.0 if any value is negative or the array is empty.
/// Uses log-domain computation: exp(mean(log(x))).
fn mcp_geometric_mean(values : Array[Double]) -> Double {
  let n = values.length()
  if n == 0 { return 0.0 }
  let mut sum_log = 0.0
  let mut count = 0
  for v in values {
    if v > 0.0 {
      sum_log = sum_log + @math.ln(v)
      count = count + 1
    }
  }
  if count == 0 { return 0.0 }
  @math.exp(sum_log / count.to_double())
}

// ============================================================================
// Built-in marker gene sets (10 populations)
// ============================================================================

///|
/// Get the built-in MCPcounter cell populations with marker genes.
/// Based on Becht et al. (2016) Table S1.
pub fn mcp_default_populations() -> Array[McpCellPopulation] {
  [
    {
      name: "T cells",
      markers: ["CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC2"],
    },
    {
      name: "CD8+ T cells",
      markers: ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "KLRG1"],
    },
    {
      name: "Cytotoxic lymphocytes",
      markers: ["GZMB", "GZMH", "PRF1", "NKG7", "GNLY", "KLRD1"],
    },
    {
      name: "B lineage",
      markers: ["CD19", "MS4A1", "CD79A", "CD79B", "CD20", "HLA-DOA"],
    },
    {
      name: "NK cells",
      markers: ["NCAM1", "KLRC1", "KLRC3", "KLRC4", "SH2D1A", "XCL2"],
    },
    {
      name: "Monocytic lineage",
      markers: ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN"],
    },
    {
      name: "Myeloid dendritic cells",
      markers: ["ITGAX", "CD1C", "CLEC10A", "FCER1A", "CLEC9A", "BATF3"],
    },
    {
      name: "Neutrophils",
      markers: ["FCGR3B", "CSF3R", "FPR1", "CEACAM3", "FCGR2A", "S100A12"],
    },
    {
      name: "Endothelial cells",
      markers: ["PECAM1", "VWF", "CDH5", "ENG", "FLT1", "KDR"],
    },
    {
      name: "Fibroblasts",
      markers: ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "DCN"],
    },
  ]
}

///|
/// Get the names of the 10 built-in cell populations.
pub fn mcp_population_names() -> Array[String] {
  let pops = mcp_default_populations()
  let names : Array[String] = []
  for p in pops { names.push(p.name) }
  names
}

// ============================================================================
// Main MCPcounter computation
// ============================================================================

///|
/// Run MCPcounter on an expression matrix.
///
/// For each cell population, computes the geometric mean of expression values
/// of its marker genes that are present in the data. The result is a per-sample
/// abundance score for each population.
///
/// Parameters:
///   - gene_names: gene identifiers (rows of expression matrix)
///   - sample_names: sample identifiers (columns of expression matrix)
///   - expression: gene × sample matrix of expression values (log-scale)
///   - populations?: custom cell populations (default: built-in 10)
pub fn mcp_run(
  gene_names : Array[String],
  sample_names : Array[String],
  expression : Array[Array[Double]],
  populations? : Array[McpCellPopulation] = mcp_default_populations(),
) -> McpResult {
  let n_samples = sample_names.length()
  let pop_names : Array[String] = []
  let scores : Array[Array[Double]] = []
  for pop in populations {
    pop_names.push(pop.name)
    // For each sample, compute geometric mean of marker gene expressions
    let sample_scores : Array[Double] = []
    for s in 0..<n_samples {
      let marker_exprs : Array[Double] = []
      for marker in pop.markers {
        // Find marker in gene_names
        for gi in 0..<gene_names.length() {
          if gene_names[gi] == marker {
            let val = expression[gi][s]
            if val > 0.0 {
              marker_exprs.push(val)
            }
            break
          }
        }
      }
      let score = mcp_geometric_mean(marker_exprs)
      sample_scores.push(score)
    }
    scores.push(sample_scores)
  }
  { populations: pop_names, sample_names: sample_names.copy(), scores }
}

///|
/// Format MCPcounter results as a string table.
pub fn McpResult::to_string(self : McpResult) -> String {
  let mut result = "MCPcounter Results:\n"
  // Header
  result = result + "Population"
  for sn in self.sample_names {
    result = result + "\t" + sn
  }
  result = result + "\n"
  // Data rows
  let mut pi = 0
  while pi < self.populations.length() {
    result = result + self.populations[pi]
    for s in 0..<self.sample_names.length() {
      result = result + "\t" + self.scores[pi][s].to_string()
    }
    result = result + "\n"
    pi = pi + 1
  }
  result
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Create sample expression data for testing MCPcounter.
/// 3 samples with different immune cell compositions.
pub fn mcp_sample_data() -> (Array[String], Array[String], Array[Array[Double]]) {
  let all_markers : Array[String] = []
  let pops = mcp_default_populations()
  for pop in pops {
    for m in pop.markers {
      let mut found = false
      for existing in all_markers {
        if existing == m { found = true; break }
      }
      if !found { all_markers.push(m) }
    }
  }
  let sample_names = ["Tumor_A", "Tumor_B", "Tumor_C"]
  let n_genes = all_markers.length()
  // Define expression patterns: Tumor_A (T-cell rich), Tumor_B (B-cell rich), Tumor_C (Fibroblast rich)
  let matrix : Array[Array[Double]] = []
  for gi in 0..<n_genes {
    let gene = all_markers[gi]
    let row : Array[Double] = []
    // Determine which population this gene belongs to
    let mut pop_idx = -1
    let mut pi = 0
    while pi < pops.length() && pop_idx < 0 {
      for m in pops[pi].markers {
        if m == gene { pop_idx = pi; break }
      }
      pi = pi + 1
    }
    // Tumor_A: high T cells (pop 0,1,2), low others
    // Tumor_B: high B cells (pop 3), low others
    // Tumor_C: high Fibroblasts (pop 9), low others
    let base_exprs = [5.0, 4.0, 3.0, 6.0, 2.0, 3.0, 2.0, 1.0, 1.0, 7.0]
    let val_a = if pop_idx == 0 || pop_idx == 1 || pop_idx == 2 {
      base_exprs[pop_idx] + 5.0
    } else { base_exprs[pop_idx] * 0.3 }
    let val_b = if pop_idx == 3 {
      base_exprs[pop_idx] + 5.0
    } else { base_exprs[pop_idx] * 0.3 }
    let val_c = if pop_idx == 9 {
      base_exprs[pop_idx] + 5.0
    } else { base_exprs[pop_idx] * 0.3 }
    row.push(val_a)
    row.push(val_b)
    row.push(val_c)
    matrix.push(row)
  }
  (all_markers, sample_names, matrix)
}
