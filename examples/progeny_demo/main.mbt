///|
/// PROGENy example: infer pathway activity from gene expression data.
///
/// PROGENy uses L2-regularized linear regression (Ridge regression)
/// to estimate pathway activity from gene expression data.
///
/// Run: moon run examples/progeny_demo/main.mbt
fn main {
  // === Create synthetic gene expression data ===
  let gene_names = ["TP53", "BRCA1", "BRCA2", "MDM2", "CDKN1A", "E2F1", "RB1", "MYC", "CCND1", "BCL2"]
  let sample_names = ["Sample_A", "Sample_B", "Sample_C", "Sample_D"]

  // Expression matrix (samples × genes)
  let expression : Array[Array[Double]] = [
    // TP53  BRCA1  BRCA2  MDM2  CDKN1A  E2F1  RB1  MYC  CCND1  BCL2
    [10.0, 8.0, 7.0, 2.0, 5.0, 6.0, 4.0, 9.0, 3.0, 1.0],   // Sample_A
    [5.0, 3.0, 4.0, 9.0, 2.0, 8.0, 7.0, 2.0, 6.0, 10.0],  // Sample_B
    [8.0, 6.0, 5.0, 3.0, 7.0, 4.0, 3.0, 10.0, 2.0, 2.0],  // Sample_C
    [3.0, 2.0, 1.0, 7.0, 4.0, 9.0, 6.0, 1.0, 5.0, 8.0],   // Sample_D
  ]

  // === Define pathway gene sets with weights ===
  let p53_genes = ["TP53", "MDM2", "CDKN1A"]
  let p53_weights = [1.0, -1.0, 1.0]
  let p53_pathway = @bio.PathwayGeneSet::new("p53", p53_genes, p53_weights)

  let cell_cycle_genes = ["E2F1", "RB1", "CCND1", "MYC"]
  let cell_cycle_weights = [1.0, -1.0, 1.0, 0.5]
  let cell_cycle_pathway = @bio.PathwayGeneSet::new("Cell_Cycle", cell_cycle_genes, cell_cycle_weights)

  let repair_genes = ["BRCA1", "BRCA2", "TP53"]
  let repair_weights = [1.0, 1.0, 0.5]
  let repair_pathway = @bio.PathwayGeneSet::new("DNA_Repair", repair_genes, repair_weights)

  let apoptosis_genes = ["BCL2", "TP53", "MYC"]
  let apoptosis_weights = [-1.0, 0.5, 0.5]
  let apoptosis_pathway = @bio.PathwayGeneSet::new("Apoptosis", apoptosis_genes, apoptosis_weights)

  let pathways = [p53_pathway, cell_cycle_pathway, repair_pathway, apoptosis_pathway]

  // === Create PROGENy data and run analysis ===
  println("=== PROGENy Pathway Activity Inference ===")
  println("  Genes: \{gene_names.length()}")
  println("  Samples: \{sample_names.length()}")
  println("  Pathways: \{pathways.length()}")

  let lambda = 0.1
  let data = @bio.ProgenyData::new(expression, sample_names, gene_names, pathways, lambda)
  let result = @bio.run_progeny(data)

  // === Display results ===
  println("\n=== Pathway Activity Matrix ===")
  let mut header = "Sample\t"
  for pathway_name in result.pathway_names {
    header = header + pathway_name + "\t"
  }
  println(header)

  let mut i = 0
  while i < result.sample_names.length() {
    let mut line = result.sample_names[i] + "\t"
    let mut j = 0
    while j < result.activity_matrix[i].length() {
      let val = result.activity_matrix[i][j]
      line = line + val.to_string() + "\t"
      j = j + 1
    }
    println(line)
    i = i + 1
  }

  // === Get specific pathway activity ===
  println("\n=== p53 Pathway Activity per Sample ===")
  let p53_activity = @bio.get_pathway_activity(result, "p53")
  let mut s = 0
  while s < p53_activity.length() {
    println("  \{result.sample_names[s]}: \{p53_activity[s].to_string()}")
    s = s + 1
  }

  // === Get sample profile ===
  println("\n=== Sample_A Pathway Profile ===")
  let profile = @bio.get_sample_profile(result, "Sample_A")
  let mut p = 0
  while p < result.pathway_names.length() {
    println("  \{result.pathway_names[p]}: \{profile[p].to_string()}")
    p = p + 1
  }

  // === Compare regularization effects ===
  println("\n=== Effect of Regularization (lambda) ===")
  let lambdas = [0.01, 0.1, 1.0, 10.0]
  for lam in lambdas {
    let d = @bio.ProgenyData::new(expression, sample_names, gene_names, pathways, lam)
    let r = @bio.run_progeny(d)
    let first_activity = r.activity_matrix[0][0]
    println("  lambda=\{lam.to_string()}: Sample_A/p53 = \{first_activity.to_string()}")
  }

  println("\nDone!")
}
