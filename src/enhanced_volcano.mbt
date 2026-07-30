///|
/// EnhancedVolcano module - Volcano plot visualization for differential expression results.
/// Inspired by Bioconductor EnhancedVolcano package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html
///   - Blighe K, Sferra G, Rabbie S. (2024) EnhancedVolcano: Publication-ready volcano plots
///     with enhanced colouring and labeling. R package version 1.22.0.

///|
/// Classification categories for volcano plot points.
pub enum VolcanoClassification {
  Up
  Down
  NonSig
}

///|
/// Create an Up classification.
pub fn volcano_classification_up() -> VolcanoClassification {
  VolcanoClassification::Up
}

///|
/// Create a Down classification.
pub fn volcano_classification_down() -> VolcanoClassification {
  VolcanoClassification::Down
}

///|
/// Create a NonSig classification.
pub fn volcano_classification_nonsig() -> VolcanoClassification {
  VolcanoClassification::NonSig
}

///|
/// A single gene result for volcano plot.
pub struct VolcanoGene {
  gene_id : String
  label : String
  log2fc : Double
  p_val : Double
  classification : VolcanoClassification
  significant : Bool
}

///|
/// Volcano plot result container.
pub struct VolcanoResult {
  genes : Array[VolcanoGene]
  p_cutoff : Double
  fc_cutoff : Double
  title : String
  x_label : String
  y_label : String
  n_up : Int
  n_down : Int
  n_nonsig : Int
}

///|
/// Classification methods.
pub fn VolcanoClassification::to_string(self : VolcanoClassification) -> String {
  match self {
    VolcanoClassification::Up => "Up"
    VolcanoClassification::Down => "Down"
    VolcanoClassification::NonSig => "NS"
  }
}

///|
/// Create a VolcanoResult from arrays of gene data.
pub fn volcano_plot(
  gene_ids : Array[String],
  log2fc_values : Array[Double],
  p_values : Array[Double],
  p_cutoff? : Double = 0.05,
  fc_cutoff? : Double = 1.0,
  title? : String = "",
  x_label? : String = "log2 Fold Change",
  y_label? : String = "-log10(p-value)",
) -> VolcanoResult {
  let n = gene_ids.length()
  let genes : Array[VolcanoGene] = Array::new()
  let mut n_up = 0
  let mut n_down = 0
  let mut n_nonsig = 0

  let mut i = 0
  while i < n {
    let gene_id = gene_ids[i]
    let lfc = log2fc_values[i]
    let pval = p_values[i]
    let abs_lfc = if lfc < 0.0 { -lfc } else { lfc }
    let significant = pval < p_cutoff && abs_lfc >= fc_cutoff

    let classification = if pval < p_cutoff && lfc >= fc_cutoff {
      n_up = n_up + 1
      VolcanoClassification::Up
    } else if pval < p_cutoff && lfc <= -fc_cutoff {
      n_down = n_down + 1
      VolcanoClassification::Down
    } else {
      n_nonsig = n_nonsig + 1
      VolcanoClassification::NonSig
    }

    genes.push(VolcanoGene::{
      gene_id,
      label: gene_id,
      log2fc: lfc,
      p_val: pval,
      classification,
      significant,
    })
    i = i + 1
  }

  VolcanoResult::{
    genes,
    p_cutoff,
    fc_cutoff,
    title,
    x_label,
    y_label,
    n_up,
    n_down,
    n_nonsig,
  }
}

///|
/// Create a volcano result with default labels.
pub fn volcano_plot_default(
  gene_ids : Array[String],
  log2fc_values : Array[Double],
  p_values : Array[Double],
) -> VolcanoResult {
  volcano_plot(
    gene_ids,
    log2fc_values,
    p_values,
    p_cutoff=0.05,
    fc_cutoff=1.0,
    title="Volcano Plot",
    x_label="log2(fold change)",
    y_label="-log10(p-value)",
  )
}

///|
/// Get the number of genes.
pub fn VolcanoResult::get_n_genes(self : VolcanoResult) -> Int {
  self.genes.length()
}

///|
/// Get the number of up-regulated genes.
pub fn VolcanoResult::get_n_up(self : VolcanoResult) -> Int {
  self.n_up
}

///|
/// Get the number of down-regulated genes.
pub fn VolcanoResult::get_n_down(self : VolcanoResult) -> Int {
  self.n_down
}

///|
/// Get the number of non-significant genes.
pub fn VolcanoResult::get_n_nonsig(self : VolcanoResult) -> Int {
  self.n_nonsig
}

///|
/// Get the p-value cutoff.
pub fn VolcanoResult::get_p_cutoff(self : VolcanoResult) -> Double {
  self.p_cutoff
}

///|
/// Get the fold-change cutoff.
pub fn VolcanoResult::get_fc_cutoff(self : VolcanoResult) -> Double {
  self.fc_cutoff
}

///|
/// Get all genes.
pub fn VolcanoResult::get_genes(self : VolcanoResult) -> Array[VolcanoGene] {
  self.genes
}

///|
/// Get significant genes only.
pub fn VolcanoResult::get_significant_genes(self : VolcanoResult) -> Array[VolcanoGene] {
  let sig : Array[VolcanoGene] = Array::new()
  let mut i = 0
  while i < self.genes.length() {
    if self.genes[i].significant {
      sig.push(self.genes[i])
    }
    i = i + 1
  }
  sig
}

///|
/// Get top significant genes by p-value.
pub fn VolcanoResult::get_top_genes(self : VolcanoResult, n : Int) -> Array[VolcanoGene] {
  let sig = self.get_significant_genes()
  let sorted = sig.copy()
  let len = sorted.length()
  let mut i = 0
  while i < len - 1 {
    let mut j = 0
    while j < len - 1 - i {
      if sorted[j].p_val > sorted[j + 1].p_val {
        let tmp = sorted[j]
        sorted.set(j, sorted[j + 1])
        sorted.set(j + 1, tmp)
      }
      j = j + 1
    }
    i = i + 1
  }

  let result : Array[VolcanoGene] = Array::new()
  let take = if n < len { n } else { len }
  let mut k = 0
  while k < take {
    result.push(sorted[k])
    k = k + 1
  }
  result
}

///|
/// Get up-regulated genes.
pub fn VolcanoResult::get_up_genes(self : VolcanoResult) -> Array[VolcanoGene] {
  let up : Array[VolcanoGene] = Array::new()
  let mut i = 0
  while i < self.genes.length() {
    if self.genes[i].classification is VolcanoClassification::Up {
      up.push(self.genes[i])
    }
    i = i + 1
  }
  up
}

///|
/// Get down-regulated genes.
pub fn VolcanoResult::get_down_genes(self : VolcanoResult) -> Array[VolcanoGene] {
  let down : Array[VolcanoGene] = Array::new()
  let mut i = 0
  while i < self.genes.length() {
    if self.genes[i].classification is VolcanoClassification::Down {
      down.push(self.genes[i])
    }
    i = i + 1
  }
  down
}

///|
/// Calculate -log10(p-value), handling p=0 as a very small number.
pub fn neg_log10_p(p_val : Double) -> Double {
  let p = if p_val <= 0.0 { 1.0e-10 } else { p_val }
  -@math.log10(p)
}

///|
/// Format the volcano result as an ASCII visualization.
pub fn VolcanoResult::to_ascii(self : VolcanoResult, width? : Int = 60, height? : Int = 20) -> String {
  let w = if width < 20 { 60 } else { width }
  let h = if height < 10 { 20 } else { height }

  // Find ranges
  let mut min_lfc = 0.0
  let mut max_lfc = 0.0
  let mut max_neglog = 0.0

  let mut i = 0
  while i < self.genes.length() {
    let g = self.genes[i]
    if g.log2fc < min_lfc { min_lfc = g.log2fc }
    if g.log2fc > max_lfc { max_lfc = g.log2fc }
    let nlp = neg_log10_p(g.p_val)
    if nlp > max_neglog { max_neglog = nlp }
    i = i + 1
  }

  // Add padding
  let lfc_range = max_lfc - min_lfc
  if lfc_range < 1.0e-9 {
    min_lfc = min_lfc - 1.0
    max_lfc = max_lfc + 1.0
  }
  let nlp_range = if max_neglog < 0.5 { 1.0 } else { max_neglog }

  // Build plot
  let plot : Array[Array[String]] = Array::new()
  let mut r = 0
  while r < h {
    let row : Array[String] = Array::new()
    let mut c = 0
    while c < w {
      row.push(".")
      c = c + 1
    }
    plot.push(row)
    r = r + 1
  }

  // Mark significant cutoffs
  let x_zero = ((0.0 - min_lfc) / (max_lfc - min_lfc) * (w.to_double() - 1.0)).to_int()
  let y_cutoff_val = neg_log10_p(self.p_cutoff)
  let y_cutoff_row = h - 1 - ((y_cutoff_val / nlp_range) * (h.to_double() - 1.0)).to_int()

  // Draw vertical dashed line at x=0
  let mut r2 = 0
  while r2 < h {
    if x_zero >= 0 && x_zero < w {
      plot[r2].set(x_zero, "|")
    }
    r2 = r2 + 1
  }

  // Draw horizontal dashed line at p_cutoff
  if y_cutoff_row >= 0 && y_cutoff_row < h {
    let mut c = 0
    while c < w {
      if plot[y_cutoff_row][c] == "." {
        plot[y_cutoff_row].set(c, "-")
      }
      c = c + 1
    }
  }

  // Plot genes
  let mut i2 = 0
  while i2 < self.genes.length() {
    let g = self.genes[i2]
    let x_pos = ((g.log2fc - min_lfc) / (max_lfc - min_lfc) * (w.to_double() - 1.0)).to_int()
    let nlp = neg_log10_p(g.p_val)
    let y_pos = h - 1 - ((nlp / nlp_range) * (h.to_double() - 1.0)).to_int()

    if x_pos >= 0 && x_pos < w && y_pos >= 0 && y_pos < h {
      let ch = if g.classification is VolcanoClassification::Up {
        "+"
      } else if g.classification is VolcanoClassification::Down {
        "-"
      } else {
        "o"
      }
      plot[y_pos].set(x_pos, ch)
    }
    i2 = i2 + 1
  }

  // Build output string
  let mut result = ""
  result = result + self.title + "\n"
  result = result + self.x_label + "\n"

  let mut r3 = 0
  while r3 < h {
    result = result + plot[r3].join("") + "\n"
    r3 = r3 + 1
  }

  result = result + self.y_label + "\n"
  result = result + "Cutoffs: p < " + self.p_cutoff.to_string() + ", |log2FC| >= " + self.fc_cutoff.to_string() + "\n"
  result = result + "Up: " + self.n_up.to_string() + ", Down: " + self.n_down.to_string() + ", NS: " + self.n_nonsig.to_string() + "\n"
  result
}

///|
/// Generate a summary string of the volcano result.
pub fn VolcanoResult::summary(self : VolcanoResult) -> String {
  let mut s = "=== Volcano Plot Summary ===\n".to_string()
  s = s + "Title: " + self.title + "\n"
  s = s + "Total genes: " + self.genes.length().to_string() + "\n"
  s = s + "Significant (p < " + self.p_cutoff.to_string() + "): " + (self.n_up + self.n_down).to_string() + "\n"
  s = s + "  Up-regulated: " + self.n_up.to_string() + "\n"
  s = s + "  Down-regulated: " + self.n_down.to_string() + "\n"
  s = s + "  Non-significant: " + self.n_nonsig.to_string() + "\n"
  s = s + "Cutoffs: |log2FC| >= " + self.fc_cutoff.to_string() + "\n"
  s
}

///|
/// Create sample data for demo purposes.
pub fn volcano_sample(n_genes : Int) -> VolcanoResult {
  let n = if n_genes < 10 { 10 } else { n_genes }
  let ids : Array[String] = Array::new()
  let lfcs : Array[Double] = Array::new()
  let pvals : Array[Double] = Array::new()

  // Generate deterministic sample data
  let names = ["TP53", "BRCA1", "BRCA2", "MYC", "KRAS", "EGFR", "PTEN", "APC", "RB1", "VHL",
    "STAT3", "AKT1", "MTOR", "JAK2", "NOTCH1", "WNT3", "CTNNB1", "SMAD4", "NFKB1", "BCL2"]

  let mut i = 0
  while i < n {
    let idx = i % names.length()
    ids.push(names[idx] + "_" + i.to_string())

    // Create varied log2FC and p-values
    let lfc_val = if i % 4 == 0 { 2.5 + (i.to_double() * 0.1) }
      else if i % 4 == 1 { -2.8 - (i.to_double() * 0.1) }
      else if i % 4 == 2 { 0.3 + (i.to_double() * 0.05) }
      else { -0.2 - (i.to_double() * 0.03) }

    let p_val = if i % 3 == 0 { 0.001 / (1.0 + i.to_double()) }
      else if i % 3 == 1 { 0.01 + (i.to_double() * 0.005) }
      else { 0.05 + (i.to_double() * 0.02) }

    lfcs.push(lfc_val)
    pvals.push(p_val)
    i = i + 1
  }

  volcano_plot_default(ids, lfcs, pvals)
}