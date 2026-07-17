///|
/// MultiAssayExperiment Demo - Multi-omics data coordination.
fn main {
  println("=== MultiAssayExperiment Demo ===\n")

  println("1. Creating multi-omics experiment from multiple data sources...")
  let rna_counts = [
    [10.0, 20.0, 30.0, 40.0],
    [50.0, 60.0, 70.0, 80.0],
    [5.0, 10.0, 15.0, 20.0],
  ]
  let rna_genes = ["TP53", "BRCA1", "EGFR"]
  let rna_samples = ["rna_s1", "rna_s2", "rna_s3", "rna_s4"]
  let chip_scores = [[1.5, 2.5, 3.5], [4.0, 5.0, 6.0]]
  let chip_peaks = ["peak1", "peak2"]
  let chip_samples = ["chip_s1", "chip_s2", "chip_s3"]
  let sample_ids = ["patient1", "patient2", "patient3", "patient4"]
  let mae = @src.create_multi_omics(
    rna_counts, rna_genes, rna_samples, chip_scores, chip_peaks, chip_samples, sample_ids,
  )
  println("   \{@src.MultiAssayExperiment::summary(mae)}")

  println("\n2. Listing experiments...")
  let names = @src.MultiAssayExperiment::experiment_names(mae)
  println("   Experiments: \{names.join(", ")}")
  let dims = @src.MultiAssayExperiment::dims(mae)
  let mut i = 0
  while i < dims.length() {
    let d = dims[i]
    println("   \{d.0}: \{d.1} rows x \{d.2} cols")
    i = i + 1
  }

  println("\n3. Getting RNA-seq experiment details...")
  match @src.MultiAssayExperiment::get_experiment(mae, "RNA-seq") {
    Some(exp) => {
      println("   Assay type: \{exp.assay_type}")
      println("   Rows: \{@src.Experiment::nrow(exp)}")
      println("   Cols: \{@src.Experiment::ncol(exp)}")
      match @src.MultiAssayExperiment::rownames(mae, "RNA-seq") {
        Some(r) => println("   Genes: \{r.join(", ")}")
        None => println("   No row names")
      }
      match @src.MultiAssayExperiment::colnames(mae, "RNA-seq") {
        Some(c) => println("   Samples: \{c.join(", ")}")
        None => println("   No col names")
      }
    }
    None => println("   RNA-seq not found")
  }

  println("\n4. Getting assay data...")
  match @src.MultiAssayExperiment::assay(mae, "RNA-seq") {
    Some(data) => {
      println("   RNA-seq data (\{data.length()} genes):")
      let mut r = 0
      while r < data.length() {
        let row_str = data[r].map(fn(x) { x.to_string() }).join(", ")
        match @src.MultiAssayExperiment::rownames(mae, "RNA-seq") {
          Some(names) => println("     \{names[r]}: [\{row_str}]")
          None => println("     [\{row_str}]")
        }
        r = r + 1
      }
    }
    None => println("   RNA-seq assay not found")
  }

  println("\n5. Finding common samples across experiments...")
  let common = @src.MultiAssayExperiment::common_samples(mae)
  println("   Common samples: \{common.length()}")

  println("\n6. Subsetting by experiment...")
  let subset_exp = @src.MultiAssayExperiment::subset_by_experiment(mae, [
    "RNA-seq",
  ])
  println(
    "   After subsetting: \{@src.MultiAssayExperiment::num_experiments(subset_exp)} experiments",
  )

  println("\n7. Adding metadata...")
  let mae2 = @src.MultiAssayExperiment::add_metadata(mae, "study", "TCGA-BRCA")
  let mae3 = @src.MultiAssayExperiment::add_metadata(
    mae2, "platform", "Illumina NovaSeq",
  )
  match @src.MultiAssayExperiment::get_metadata(mae3, "study") {
    Some(v) => println("   Study: \{v}")
    None => println("   Study not found")
  }
  match @src.MultiAssayExperiment::get_metadata(mae3, "platform") {
    Some(v) => println("   Platform: \{v}")
    None => println("   Platform not found")
  }

  println("\n8. Creating custom MultiAssayExperiment from scratch...")
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[100.0, 200.0], [300.0, 400.0]],
    row_names=["GAPDH", "TP53"],
    col_names=["sample1", "sample2"],
  )
  let snv_exp = @src.Experiment::new(
    name="SNV",
    assay_type="variants",
    data=[[1.0, 0.0], [0.0, 1.0]],
    row_names=["rs123", "rs456"],
    col_names=["sample1", "sample2"],
  )
  let sm = @src.SampleMap::new(
    assays=["RNA-seq", "RNA-seq", "SNV", "SNV"],
    primary_ids=["p1", "p2", "p1", "p2"],
    colnames=["sample1", "sample2", "sample1", "sample2"],
  )
  let col_data : Array[(String, Map[String, String])] = [
    ("p1", Map([("tissue", "tumor"), ("stage", "III")])),
    ("p2", Map([("tissue", "normal"), ("stage", "0")])),
  ]
  let custom_mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp, snv_exp],
    sample_map=sm,
    col_data~,
  )
  println("   \{@src.MultiAssayExperiment::summary(custom_mae)}")

  println("\n9. Looking up sample column data...")
  match @src.MultiAssayExperiment::sample_col_data(custom_mae, "p1") {
    Some(d) => {
      match d.get("tissue") {
        Some(v) => println("   p1 tissue: \{v}")
        None => println("   tissue not found")
      }
      match d.get("stage") {
        Some(v) => println("   p1 stage: \{v}")
        None => println("   stage not found")
      }
    }
    None => println("   p1 not found")
  }

  println("\n=== MultiAssayExperiment Demo Complete ===")
}
