///|
fn demo_data() -> (
  Array[Array[Double]],
  Array[String],
  Array[String],
  Array[String],
  Array[String],
  Array[String],
  Array[String],
) {
  let genes = [
    "response_up", "stable", "T_cell_marker", "donor_signal", "response_down",
  ]
  let counts : Array[Array[Double]] = [[], [], [], [], []]
  let cells : Array[String] = []
  let samples : Array[String] = []
  let clusters : Array[String] = []
  let treatments : Array[String] = []
  let donors : Array[String] = []
  for donor in 0..<4 {
    for treatment in 0..<2 {
      let sample = "D" +
        (donor + 1).to_string() +
        (if treatment == 0 { "_C" } else { "_T" })
      for cluster in 0..<2 {
        let cluster_name = if cluster == 0 { "T_cell" } else { "Monocyte" }
        for replicate in 0..<3 {
          cells.push(
            sample + "_" + cluster_name + "_" + (replicate + 1).to_string(),
          )
          samples.push(sample)
          clusters.push(cluster_name)
          treatments.push(if treatment == 0 { "Control" } else { "Treated" })
          donors.push("D" + (donor + 1).to_string())
          counts[0].push(
            (10 + treatment * 10 + donor + (if replicate == 2 { 1 } else { 0 })).to_double(),
          )
          counts[1].push((8 + replicate % 2).to_double())
          counts[2].push(
            (if cluster == 0 { 18 + treatment } else { 4 + treatment }).to_double(),
          )
          counts[3].push((6 + donor * 2 + replicate).to_double())
          counts[4].push(
            (18 - treatment * 8 + cluster + (if replicate == 1 { 1 } else { 0 })).to_double(),
          )
        }
      }
    }
  }
  (counts, genes, cells, samples, clusters, treatments, donors)
}

///|
fn main {
  println("=== Bioconductor dreamlet Demo ===")
  let (counts, genes, cells, samples, clusters, treatments, donors) = demo_data()

  println("\n1. Build a SingleCellExperiment with cohort metadata")
  let sce = @src.SingleCellExperiment::new(counts, genes, cells)
  sce.col_data["sample_id"] = samples
  sce.col_data["cluster_id"] = clusters
  sce.col_data["Treatment"] = treatments
  sce.col_data["Donor"] = donors
  println(
    "   " +
    @src.sce_get_n_genes(sce).to_string() +
    " genes x " +
    @src.sce_get_n_cells(sce).to_string() +
    " cells",
  )

  println("\n2. Aggregate raw counts by sample and cell type")
  let pseudobulk = @src.dreamlet_aggregate_sce(sce, "sample_id", "cluster_id", metadata_fields=[
    "Treatment", "Donor",
  ]) catch {
    _ => abort("failed to aggregate the SingleCellExperiment")
  }
  println(pseudobulk.summary())

  println("3. Apply filtering, TMM normalization, and voom weights")
  let model = @src.dreamlet_model([
    @src.dreamlet_categorical_effect("Treatment", reference="Control"),
    @src.dreamlet_random_effect("Donor"),
  ]) catch {
    _ => abort("failed to build the repeated-measures model")
  }
  let config = @src.DreamletProcessConfig::create(
    min_cells=3,
    min_count=2.0,
    min_samples=4,
    min_prop=0.5,
    min_total_count=10.0,
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("failed to build dreamlet processing controls")
  }
  let processed = @src.dreamlet_process_assays(pseudobulk, model, config~) catch {
    _ => abort("failed to process pseudobulk assays")
  }
  println(processed.summary())
  for assay in processed.assays {
    println(
      "   " +
      assay.cluster_id +
      ": " +
      assay.gene_names.length().to_string() +
      " genes, TMM reference=" +
      assay.tmm_reference_sample,
    )
  }

  println("4. Fit the donor random-intercept treatment contrast")
  let result = @src.dreamlet(
    processed,
    "Treatment:Treated",
    contrast_name="Treated-Control",
    max_iterations=80,
    tolerance=0.01,
  ) catch {
    _ => abort("failed to run dreamlet differential expression")
  }
  println(result.summary())

  println("5. Inspect the study-wide top table")
  let top = result.top_table() catch {
    _ => abort("failed to create the dreamlet top table")
  }
  let limit = if top.length() < 6 { top.length() } else { 6 }
  for index in 0..<limit {
    let gene = top[index]
    println(
      "   " +
      gene.cluster_id +
      "/" +
      gene.gene_name +
      ": logFC=" +
      gene.estimate.to_string() +
      ", p=" +
      gene.p_value.to_string() +
      ", within-FDR=" +
      gene.within_cluster_fdr.to_string() +
      ", study-FDR=" +
      gene.study_wide_fdr.to_string(),
    )
  }

  println("\n=== Demo Complete ===")
}
