///|
fn make_taxa_tree() -> @src.Tree {
  let taxon_a = @src.Clade::new(name=Some("TaxonA"))
  let taxon_b = @src.Clade::new(name=Some("TaxonB"))
  let taxon_c = @src.Clade::new(name=Some("TaxonC"))
  let taxon_d = @src.Clade::new(name=Some("TaxonD"))
  let firmicutes = @src.Clade::new(name=Some("Firmicutes"), clades=[
    taxon_a, taxon_b,
  ])
  let bacteroidota = @src.Clade::new(name=Some("Bacteroidota"), clades=[
    taxon_c, taxon_d,
  ])
  let bacteria = @src.Clade::new(name=Some("Bacteria"), clades=[
    firmicutes, bacteroidota,
  ])
  @src.Tree::new(bacteria, rooted=true, name=Some("taxonomy"))
}

///|
fn make_sample_tree() -> @src.Tree {
  let control_1 = @src.Clade::new(name=Some("Control1"))
  let control_2 = @src.Clade::new(name=Some("Control2"))
  let treated_1 = @src.Clade::new(name=Some("Treated1"))
  let treated_2 = @src.Clade::new(name=Some("Treated2"))
  let controls = @src.Clade::new(name=Some("Controls"), clades=[
    control_1, control_2,
  ])
  let treated = @src.Clade::new(name=Some("Treated"), clades=[
    treated_1, treated_2,
  ])
  let samples = @src.Clade::new(name=Some("Samples"), clades=[controls, treated])
  @src.Tree::new(samples, rooted=true, name=Some("sample_groups"))
}

///|
fn print_assay(data : Array[Array[Double]]) -> Unit {
  for row in data {
    println("    [" + row.map(fn(value) { value.to_string() }).join(", ") + "]")
  }
}

///|
fn main {
  println("=== TreeSummarizedExperiment Demo ===")

  let experiment = @src.summarized_experiment(
    Map([
      (
        "counts",
        [
          [10.0, 12.0, 20.0, 24.0],
          [5.0, 7.0, 15.0, 17.0],
          [30.0, 32.0, 18.0, 20.0],
          [8.0, 10.0, 14.0, 16.0],
        ],
      ),
    ]),
    [("TaxonA", 0, 0), ("TaxonB", 0, 0), ("TaxonC", 0, 0), ("TaxonD", 0, 0)],
    [
      Map([("sample", "Control1")]),
      Map([("sample", "Control2")]),
      Map([("sample", "Treated1")]),
      Map([("sample", "Treated2")]),
    ],
    Map([("study", "microbiome-treatment")]),
  )

  let tse = @src.TreeSummarizedExperiment::new(
    experiment~,
    row_tree=Some(make_taxa_tree()),
    row_node_labels=["TaxonA", "TaxonB", "TaxonC", "TaxonD"],
    col_tree=Some(make_sample_tree()),
    col_node_labels=["Control1", "Control2", "Treated1", "Treated2"],
    reference_sequences=["ACGT", "AAGT", "CCGT", "TTGT"],
  )

  println("\n1. Container and links")
  println("  Valid: \{tse.is_valid()}")
  println("  Dimensions: \{tse.nrow()} rows x \{tse.ncol()} columns")
  for link in tse.row_links() {
    println(
      "  \{link.node_alias()} -> \{link.node_label()}, leaf=\{link.is_leaf()}",
    )
  }

  println("\n2. Tree queries")
  let descendants = @src.tse_find_descendants(
    make_taxa_tree(),
    "Firmicutes",
    leaves_only=true,
  )
  println("  Firmicutes descendants: " + descendants.join(", "))
  println(
    "  TaxonC ancestors: " +
    @src.tse_find_ancestors(make_taxa_tree(), "TaxonC").join(" -> "),
  )

  println("\n3. Subset by an internal tree node")
  match tse.subset_by_row_nodes(["Firmicutes"]) {
    Some(subset) =>
      match subset.assay("counts") {
        Some(data) => print_assay(data)
        None => ()
      }
    None => println("  Row tree is unavailable")
  }

  println("\n4. Aggregate taxa to phylum level (sum)")
  match tse.aggregate_rows(["Firmicutes", "Bacteroidota"]) {
    Some(aggregated) =>
      match aggregated.assay("counts") {
        Some(data) => print_assay(data)
        None => ()
      }
    None => println("  No matching row-tree nodes")
  }

  println("\n5. Aggregate samples to treatment groups (mean)")
  match
    tse.aggregate_cols(
      ["Controls", "Treated"],
      aggregation=@src.tse_aggregation_mean(),
    ) {
    Some(aggregated) =>
      match aggregated.assay("counts") {
        Some(data) => print_assay(data)
        None => ()
      }
    None => println("  No matching column-tree nodes")
  }

  println("\n=== Demo Complete ===")
}
