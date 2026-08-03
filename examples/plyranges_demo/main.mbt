///|
/// plyranges_demo - Example demonstrating tidy genomic data manipulation
/// Inspired by Bioconductor plyranges package.

///|
fn main {
  println("=== plyranges Demo ===")
  println("Tidy genomic data manipulation inspired by Bioconductor plyranges")
  println("")

  // Create GRanges object
  let gr = @src.tgr_create(
    ["chr1", "chr1", "chr2", "chr2", "chr3"],
    [100, 500, 200, 800, 300],
    [200, 600, 300, 900, 400],
    ["+", "-", "+", "-", "+"],
  )
  println("Created GRanges with \{gr.seqnames.length()} ranges")

  // Filter by seqname
  let chr1_ranges = @src.tgr_filter_seqname(gr, "chr1")
  println("Filtered chr1 ranges: \{chr1_ranges.seqnames.length()} ranges")

  // Filter by strand
  let plus_strand = @src.tgr_filter_strand(gr, "+")
  println("Plus strand ranges: \{plus_strand.seqnames.length()} ranges")

  // Filter by start position
  let after_300 = @src.tgr_filter_start(gr, 300)
  println("Ranges starting after 300: \{after_300.seqnames.length()}")

  // Add metadata
  let with_gc = @src.tgr_mutate(gr, "gc_content", [0.45, 0.62, 0.38, 0.71, 0.52])
  println("Added gc_content metadata")

  let with_gene = @src.tgr_mutate_str(with_gc, "gene_name", ["BRCA1", "TP53", "EGFR", "MYC", "KRAS"])
  println("Added gene_name metadata")

  // Arrange by start position
  let arranged = @src.tgr_arrange(with_gene)
  println("Arranged by start position")
  println("First start: \{arranged.starts[0]}, Last start: \{arranged.starts[arranged.starts.length() - 1]}")

  // Select columns
  let selected = @src.tgr_select(arranged, ["gc_content", "gene_name"])
  println("Selected metadata columns")

  // Rename
  let renamed = @src.tgr_rename(selected, "gc_content", "GC")
  println("Renamed gc_content to GC")

  // Summary
  let summary = @src.tgr_summary(gr)
  println("")
  println("Summary:")
  println("  Ranges: \{summary.n_ranges}")
  println("  Seqnames: \{summary.n_seqnames}")
  println("  Total bases: \{summary.total_bases}")
  println("  Mean width: \{summary.mean_width}")

  println("")
  println("=== plyranges Demo Complete ===")
}
