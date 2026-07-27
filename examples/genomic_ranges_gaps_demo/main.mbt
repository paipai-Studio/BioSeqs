///|
/// GRanges advanced operations demo.
///
/// Demonstrates GRangesLocal construction, gaps, reduce, disjoin,
/// set operations (setdiff, intersect, union), coverage computation,
/// and promoter region calculation.
///
/// Run: moon run examples/genomic_ranges_gaps_demo/main.mbt
fn main {
  println("=== GRanges Advanced Operations Demo ===")
  println("")

  println("1. Creating GRangesLocal with sample ranges")
  let seqnames = ["chr1", "chr1", "chr1", "chr2", "chr2"]
  let starts = [1, 50, 150, 1, 100]
  let ends = [10, 80, 200, 50, 150]
  let strands = ["+", "+", "-", "+", "-"]
  let gr = @bio.grg_granges(seqnames, starts, ends, strands=strands)
  println("   Created \{@bio.grg_length(gr)} ranges")
  let mut i = 0
  while i < gr.seqnames.length() {
    println(
      "     \{gr.seqnames[i]}:\{gr.starts[i]}-\{gr.ends[i]} (strand=\{gr.strands[i]})",
    )
    i = i + 1
  }

  println("")
  println("2. Finding gaps between ranges")
  let seqlengths : Map[String, Int] = Map([], capacity=2)
  seqlengths["chr1"] = 300
  seqlengths["chr2"] = 200
  let gaps = @bio.grg_gaps(gr, seqlengths)
  println("   Gaps found: \{gaps.seqnames.length()}")
  i = 0
  while i < gaps.seqnames.length() {
    println(
      "     \{gaps.seqnames[i]}:\{gaps.starts[i]}-\{gaps.ends[i]} (length=\{gaps.gap_lengths[i]})",
    )
    i = i + 1
  }

  println("")
  println("3. Reduce - merge overlapping/adjacent ranges")
  // Create ranges that overlap
  let gr_overlap = @bio.grg_granges(
    ["chr1", "chr1", "chr1"],
    [1, 5, 20],
    [10, 15, 30],
  )
  println("   Input ranges:")
  i = 0
  while i < gr_overlap.seqnames.length() {
    println(
      "     \{gr_overlap.seqnames[i]}:\{gr_overlap.starts[i]}-\{gr_overlap.ends[i]}",
    )
    i = i + 1
  }
  let reduced = @bio.grg_reduce(gr_overlap)
  println("   After reduce:")
  i = 0
  while i < reduced.seqnames.length() {
    println(
      "     \{reduced.seqnames[i]}:\{reduced.starts[i]}-\{reduced.ends[i]} (merged \{reduced.merged_count[i]} ranges)",
    )
    i = i + 1
  }

  println("")
  println("4. Disjoin - split into non-overlapping segments")
  let disjoined = @bio.grg_disjoin(gr_overlap)
  println("   Disjoined into segments: \{disjoined.seqnames.length()}")
  i = 0
  while i < disjoined.seqnames.length() {
    println(
      "     \{disjoined.seqnames[i]}:\{disjoined.starts[i]}-\{disjoined.ends[i]} (segment id=\{disjoined.segment_ids[i]})",
    )
    i = i + 1
  }

  println("")
  println("5. Set operations")
  let gr_a = @bio.grg_granges(
    ["chr1", "chr1"],
    [1, 50],
    [30, 80],
    strands=["+", "+"],
  )
  let gr_b = @bio.grg_granges(
    ["chr1"],
    [10],
    [20],
    strands=["-"],
  )

  println("   Set difference (A \\ B):")
  let setdiff_result = @bio.grg_setdiff(gr_a, gr_b)
  i = 0
  while i < setdiff_result.seqnames.length() {
    println(
      "     \{setdiff_result.seqnames[i]}:\{setdiff_result.starts[i]}-\{setdiff_result.ends[i]}",
    )
    i = i + 1
  }

  println("   Intersection (A ∩ B):")
  let intersect_result = @bio.grg_intersect(gr_a, gr_b)
  i = 0
  while i < intersect_result.seqnames.length() {
    println(
      "     \{intersect_result.seqnames[i]}:\{intersect_result.starts[i]}-\{intersect_result.ends[i]}",
    )
    i = i + 1
  }

  println("   Union (A ∪ B):")
  let union_result = @bio.grg_union(gr_a, gr_b)
  i = 0
  while i < union_result.seqnames.length() {
    println(
      "     \{union_result.seqnames[i]}:\{union_result.starts[i]}-\{union_result.ends[i]}",
    )
    i = i + 1
  }

  println("")
  println("6. Computing coverage")
  let gr_cov = @bio.grg_granges(
    ["chr1", "chr1", "chr1"],
    [1, 3, 8],
    [5, 8, 12],
  )
  let cov_sl : Map[String, Int] = Map([], capacity=1)
  cov_sl["chr1"] = 15
  let coverage = @bio.grg_compute_coverage(gr_cov, cov_sl)
  let cov_arr = coverage["chr1"]
  println("   Coverage array for chr1 (length 15):")
  println("     \{cov_arr}")

  println("")
  println("7. Computing promoter regions")
  let gr_genes = @bio.grg_granges(
    ["chr1", "chr1", "chr2"],
    [100, 500, 200],
    [200, 600, 300],
    strands=["+", "-", "*"],
  )
  println("   Gene ranges:")
  i = 0
  while i < gr_genes.seqnames.length() {
    println(
      "     \{gr_genes.seqnames[i]}:\{gr_genes.starts[i]}-\{gr_genes.ends[i]} (strand=\{gr_genes.strands[i]})",
    )
    i = i + 1
  }

  let promoters = @bio.grg_promoters(gr_genes, upstream=50, downstream=20)
  println("   Promoter regions (upstream=50, downstream=20):")
  i = 0
  while i < promoters.seqnames.length() {
    println(
      "     \{promoters.seqnames[i]}:\{promoters.starts[i]}-\{promoters.ends[i]} (strand=\{promoters.strands[i]})",
    )
    i = i + 1
  }

  println("")
  println("8. Trimming ranges to chromosome bounds")
  let gr_trim = @bio.grg_granges(
    ["chr1", "chr1", "chr2"],
    [-5, 8, 1],
    [10, 250, 300],
  )
  let trim_sl : Map[String, Int] = Map([], capacity=2)
  trim_sl["chr1"] = 100
  trim_sl["chr2"] = 200
  let trimmed = @bio.grg_trim(gr_trim, trim_sl)
  println("   After trimming to chromosome bounds:")
  i = 0
  while i < trimmed.seqnames.length() {
    println(
      "     \{trimmed.seqnames[i]}:\{trimmed.starts[i]}-\{trimmed.ends[i]}",
    )
    i = i + 1
  }

  println("")
  println("=== Demo Complete ===")
}