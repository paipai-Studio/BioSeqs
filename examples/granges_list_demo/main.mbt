///|
fn print_element(name : String, ranges : @src.GRanges) -> Unit {
  println("  \{name} (\{@src.granges_length(ranges)} ranges)")
  for index = 0; index < @src.granges_length(ranges); index = index + 1 {
    println(
      "    \{ranges.seqnames[index]}:\{ranges.starts[index]}-\{ranges.ends[index]}",
    )
  }
}

///|
fn main {
  println("=== Bioconductor GRangesList Demo ===")

  let transcripts = @src.granges_list(
    [
      @src.granges(["chr1", "chr1"], [(100, 109), (200, 209)], [
        @src.strand_plus(),
        @src.strand_plus(),
      ]),
      @src.granges(["chr1", "chr1"], [(150, 159), (300, 309)], [
        @src.strand_minus(),
        @src.strand_minus(),
      ]),
      @src.granges_single("chr2", 50, 59, @src.strand_star()),
    ],
    names=["txA", "txB", "txC"],
    element_metadata=[
      Map([("gene", "geneA")]),
      Map([("gene", "geneB")]),
      Map([("gene", "geneC")]),
    ],
  ) catch {
    _ => abort("failed to construct transcript ranges")
  }

  println("\n1. Compound genomic features")
  println("  " + transcripts.summary())
  let names = transcripts.names()
  for index = 0; index < transcripts.length(); index = index + 1 {
    match transcripts.get(index) {
      Some(ranges) => print_element(names[index], ranges)
      None => abort("missing transcript ranges")
    }
  }

  println("\n2. Unlist and relist")
  let flat = transcripts.unlist()
  println("  Flattened range count: \{@src.granges_length(flat)}")
  let rebuilt = transcripts.relist(flat) catch {
    _ => abort("failed to relist flattened ranges")
  }
  println(
    "  Rebuilt element sizes: " +
    rebuilt.element_lengths().map(fn(size) { size.to_string() }).join(", "),
  )

  println("\n3. Feature-level overlaps use member exons")
  let queries = @src.granges(["chr1", "chr1"], [(120, 130), (205, 206)], [
    @src.strand_plus(),
    @src.strand_plus(),
  ])
  println("  Query 0 lies inside txA bounds but only in its intron.")
  for hit in transcripts.find_overlaps(queries) {
    println("  \{names[hit.0]} overlaps query \{hit.1}")
  }

  println("\n4. Grouped RangedSummarizedExperiment")
  let experiment = @src.RangedSummarizedExperiment::new_with_range_groups(
    assays=Map([("counts", [[12.0, 18.0], [25.0, 30.0], [8.0, 11.0]])]),
    row_range_groups=transcripts,
    col_data=[
      Map([("sample", "S1"), ("condition", "control")]),
      Map([("sample", "S2"), ("condition", "treated")]),
    ],
    row_names=names,
    metadata=Map([("organism", "human")]),
  ) catch {
    _ => abort("failed to construct grouped experiment")
  }
  println("  " + experiment.summary())
  let overlapping = experiment.subset_by_overlaps(queries) catch {
    _ => abort("failed to subset grouped experiment")
  }
  println("  Overlapping assay rows: " + overlapping.row_names().join(", "))

  let shifted = experiment.shift(10) catch {
    _ => abort("failed to shift grouped experiment")
  }
  match shifted.row_range_groups() {
    Some(groups) =>
      match groups.get(0) {
        Some(ranges) =>
          println(
            "  Shifted txA exon starts: " +
            ranges.starts.map(fn(start) { start.to_string() }).join(", "),
          )
        None => abort("missing shifted transcript")
      }
    None => abort("grouping was not preserved")
  }

  println("\n5. Coverage uses all member ranges")
  let coverage = experiment.coverage(Map([("chr1", 320), ("chr2", 80)]))
  println(
    "  chr1 coverage at exon/intron positions 100/125: \{coverage["chr1"][99]}/\{coverage["chr1"][124]}",
  )

  println("\n=== Demo Complete ===")
}
