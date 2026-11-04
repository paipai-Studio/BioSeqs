///|
fn print_ranges(names : Array[String], ranges : @src.GRanges) -> Unit {
  for index = 0; index < @src.granges_length(ranges); index = index + 1 {
    println(
      "  \{names[index]}: \{ranges.seqnames[index]}:\{ranges.starts[index]}-\{ranges.ends[index]}",
    )
  }
}

///|
fn main {
  println("=== RangedSummarizedExperiment Demo ===")

  let experiment = @src.RangedSummarizedExperiment::new(
    assays=Map([
      (
        "counts",
        [
          [30.0, 31.0, 32.0],
          [10.0, 11.0, 12.0],
          [20.0, 21.0, 22.0],
          [40.0, 41.0, 42.0],
        ],
      ),
      (
        "logcounts",
        [
          [3.40, 3.43, 3.47],
          [2.40, 2.48, 2.56],
          [3.04, 3.09, 3.14],
          [3.71, 3.74, 3.76],
        ],
      ),
    ]),
    row_ranges=@src.granges(
      ["chr2", "chr1", "chr1", "chr3"],
      [(300, 349), (100, 149), (200, 249), (50, 99)],
      [
        @src.strand_plus(),
        @src.strand_plus(),
        @src.strand_minus(),
        @src.strand_star(),
      ],
    ),
    col_data=[
      Map([("sample", "S1"), ("condition", "control")]),
      Map([("sample", "S2"), ("condition", "treated")]),
      Map([("sample", "S3"), ("condition", "treated")]),
    ],
    row_data=[
      Map([("biotype", "protein_coding")]),
      Map([("biotype", "lncRNA")]),
      Map([("biotype", "protein_coding")]),
      Map([("biotype", "protein_coding")]),
    ],
    row_names=["geneC", "geneA", "geneB", "geneD"],
    metadata=Map([("study", "range-aware-expression")]),
  ) catch {
    _ => abort("failed to construct RangedSummarizedExperiment")
  }

  println("\n1. Container and row ranges")
  println("  " + experiment.summary())
  print_ranges(experiment.row_names(), experiment.row_ranges())

  let regions = @src.granges(["chr1", "chr1"], [(120, 170), (220, 230)], [
    @src.strand_plus(),
    @src.strand_minus(),
  ])
  println("\n2. Strand-aware overlaps")
  for hit in experiment.find_overlaps(regions) {
    println("  \{experiment.row_names()[hit.0]} overlaps query \{hit.1}")
  }
  let overlapping = experiment.subset_by_overlaps(regions) catch {
    _ => abort("failed to subset overlapping rows")
  }
  println("  Retained rows: " + overlapping.row_names().join(", "))

  println("\n3. Nearest ranges and interval distances")
  let nearest = experiment.nearest(regions)
  let distances = experiment.distance_to_nearest(regions)
  for index = 0; index < experiment.nrow(); index = index + 1 {
    println(
      "  \{experiment.row_names()[index]}: subject=\{nearest[index]}, distance=\{distances[index]}",
    )
  }

  println("\n4. Strand-aware promoter ranges")
  let promoter_experiment = experiment.promoters(upstream=100, downstream=20) catch {
    _ => abort("failed to create promoter ranges")
  }
  print_ranges(
    promoter_experiment.row_names(),
    promoter_experiment.row_ranges(),
  )

  println("\n5. Sort ranges and assays together")
  let sorted = experiment.sort() catch {
    _ => abort("failed to sort experiment")
  }
  println("  Row order: " + sorted.row_names().join(", "))
  match sorted.assay("counts") {
    Some(counts) =>
      println(
        "  First-sample counts: " +
        counts.map(fn(row) { row[0].to_string() }).join(", "),
      )
    None => println("  counts assay is unavailable")
  }

  println("\n=== Demo Complete ===")
}
