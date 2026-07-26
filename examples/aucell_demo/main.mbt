///|
/// AUCell example: gene set scoring for single-cell data using AUC.
///
/// AUCell uses Area Under Curve (AUC) to score gene sets per cell.
/// It ranks genes by expression and computes the enrichment of gene sets.
///
/// Run: moon run examples/aucell_demo/main.mbt
fn main {
  // === Create synthetic single-cell expression data ===
  let gene_names = ["TP53", "BRCA1", "MYC", "CCND1", "BCL2", "BAX", "CASP3", "CDK1", "E2F1", "RB1"]
  let cell_names = ["CellTypeA_1", "CellTypeA_2", "CellTypeA_3", "CellTypeB_1", "CellTypeB_2", "CellTypeB_3"]

  // Expression matrix (cells × genes)
  let expression : Array[Array[Double]] = [
    // TP53  BRCA1  MYC   CCND1  BCL2  BAX   CASP3  CDK1  E2F1  RB1
    [10.0, 8.0, 2.0, 1.0, 9.0, 7.0, 8.0, 1.0, 2.0, 3.0],   // CellTypeA_1
    [9.0, 7.0, 3.0, 2.0, 8.0, 6.0, 7.0, 2.0, 3.0, 4.0],    // CellTypeA_2
    [11.0, 9.0, 1.0, 0.5, 10.0, 8.0, 9.0, 0.5, 1.0, 2.0],  // CellTypeA_3
    [2.0, 3.0, 10.0, 9.0, 1.0, 2.0, 1.0, 10.0, 9.0, 8.0],  // CellTypeB_1
    [1.0, 2.0, 11.0, 10.0, 0.5, 1.0, 0.5, 11.0, 10.0, 9.0], // CellTypeB_2
    [3.0, 4.0, 8.0, 7.0, 2.0, 3.0, 2.0, 8.0, 7.0, 6.0],    // CellTypeB_3
  ]

  // === Define gene sets (signatures) ===
  let apoptosis_genes = ["TP53", "BAX", "CASP3", "BCL2"]
  let apoptosis_set = @bio.GeneSet::new("Apoptosis", apoptosis_genes)

  let cell_cycle_genes = ["MYC", "CCND1", "CDK1", "E2F1", "RB1"]
  let cell_cycle_set = @bio.GeneSet::new("Cell_Cycle", cell_cycle_genes)

  let dna_repair_genes = ["BRCA1", "TP53"]
  let dna_repair_set = @bio.GeneSet::new("DNA_Repair", dna_repair_genes)

  let gene_sets = [apoptosis_set, cell_cycle_set, dna_repair_set]

  // === Create AUCell data and run analysis ===
  println("=== AUCell Gene Set Scoring ===")
  println("  Cells: \{cell_names.length()}")
  println("  Genes: \{gene_names.length()}")
  println("  Gene sets: \{gene_sets.length()}")

  let data = @bio.AUCellData::new(expression, cell_names, gene_names, gene_sets)
  let result = @bio.run_aucell(data)

  // === Display raw scores ===
  println("\n=== Raw AUCell Scores ===")
  let mut header = "Cell\t"
  for gs_name in result.gene_set_names {
    header = header + gs_name + "\t"
  }
  println(header)

  let mut i = 0
  while i < result.cell_names.length() {
    let mut line = result.cell_names[i] + "\t"
    let mut j = 0
    while j < result.score_matrix[i].length() {
      let val = result.score_matrix[i][j]
      line = line + val.to_string() + "\t"
      j = j + 1
    }
    println(line)
    i = i + 1
  }

  // === Normalize scores ===
  println("\n=== Normalized AUCell Scores (min-max per gene set) ===")
  let normalized = @bio.normalize_aucell_scores(result)
  let mut header2 = "Cell\t"
  for gs_name in normalized.gene_set_names {
    header2 = header2 + gs_name + "\t"
  }
  println(header2)

  let mut c = 0
  while c < normalized.cell_names.length() {
    let mut line = normalized.cell_names[c] + "\t"
    let mut s = 0
    while s < normalized.score_matrix[c].length() {
      let val = normalized.score_matrix[c][s]
      line = line + val.to_string() + "\t"
      s = s + 1
    }
    println(line)
    c = c + 1
  }

  // === Get cell scores ===
  println("\n=== Gene Set Scores for CellTypeA_1 ===")
  let cell_scores = @bio.get_cell_scores(normalized, "CellTypeA_1")
  let mut cs = 0
  while cs < normalized.gene_set_names.length() {
    println("  \{normalized.gene_set_names[cs]}: \{cell_scores[cs].to_string()}")
    cs = cs + 1
  }

  // === Get gene set scores across cells ===
  println("\n=== Apoptosis Scores Across Cells ===")
  let apoptosis_scores = @bio.get_set_scores(normalized, "Apoptosis")
  let mut a = 0
  while a < normalized.cell_names.length() {
    println("  \{normalized.cell_names[a]}: \{apoptosis_scores[a].to_string()}")
    a = a + 1
  }

  println("\nDone!")
}
