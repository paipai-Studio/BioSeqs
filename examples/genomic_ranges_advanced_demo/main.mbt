///|
/// GenomicRanges Advanced Demo - Tiling, Windows, and Set Operations
///
/// This demo shows how to:
/// 1. Tile genomic ranges into non-overlapping sub-ranges
/// 2. Generate sliding windows over ranges
/// 3. Tile an entire genome by chromosome
/// 4. Compute coverage in sliding windows
/// 5. Perform set operations (gaps, subtract)
/// 6. Create promoter regions

fn main {
  println("=== GenomicRanges Advanced: Tiling and Windowing Demo ===\n")

  // Step 1: Create GRanges
  println("Step 1: Creating GRanges...")
  let gr = @bio.granges(
    ["chr1", "chr1", "chr2"],
    [(100, 300), (500, 900), (200, 600)],
    [@bio.strand_plus(), @bio.strand_minus(), @bio.strand_star()],
  )
  println("  Original ranges: \{gr.seqnames.length()}")
  for i = 0; i < gr.seqnames.length(); i = i + 1 {
    println("  \{gr.seqnames[i]}: \{gr.starts[i]}-\{gr.ends[i]} (\{@bio.strand_to_char(gr.strands[i])})")
  }

  // Step 2: Tile by fixed width
  println("\nStep 2: Tiling by width=100...")
  let tiled = @bio.granges_tile(gr, width=100)
  println("  Tiled ranges: \{tiled.seqnames.length()}")
  for i = 0; i < tiled.seqnames.length(); i = i + 1 {
    println("  \{tiled.seqnames[i]}: \{tiled.starts[i]}-\{tiled.ends[i]}")
  }

  // Step 3: Tile by number
  println("\nStep 3: Tiling by n=3 (3 tiles per range)...")
  let tiled_n = @bio.granges_tile(gr, n=3)
  println("  Tiled ranges: \{tiled_n.seqnames.length()}")
  for i = 0; i < tiled_n.seqnames.length(); i = i + 1 {
    println("  \{tiled_n.seqnames[i]}: \{tiled_n.starts[i]}-\{tiled_n.ends[i]}")
  }

  // Step 4: Sliding windows
  println("\nStep 4: Sliding windows (width=100, step=50)...")
  let windows = @bio.granges_sliding_windows(gr, 100, 50)
  println("  Window count: \{windows.seqnames.length()}")
  for i = 0; i < windows.seqnames.length(); i = i + 1 {
    println("  \{windows.seqnames[i]}: \{windows.starts[i]}-\{windows.ends[i]}")
  }

  // Step 5: Tile entire genome
  println("\nStep 5: Tiling entire genome (width=500)...")
  let seqlengths = [("chr1", 2000), ("chr2", 1500)]
  let genome_tiles = @bio.tile_genome(seqlengths, 500)
  println("  Genome tiles: \{genome_tiles.seqnames.length()}")
  for i = 0; i < genome_tiles.seqnames.length(); i = i + 1 {
    println("  \{genome_tiles.seqnames[i]}: \{genome_tiles.starts[i]}-\{genome_tiles.ends[i]}")
  }

  // Step 6: Tile by ranges (exactly n tiles per range)
  println("\nStep 6: Tile by ranges (n=4 per range)...")
  let tile_by_r = @bio.tile_by_ranges(gr, 4)
  println("  Tiled ranges: \{tile_by_r.seqnames.length()}")
  for i = 0; i < tile_by_r.seqnames.length(); i = i + 1 {
    println("  \{tile_by_r.seqnames[i]}: \{tile_by_r.starts[i]}-\{tile_by_r.ends[i]} (width=\{tile_by_r.widths[i]})")
  }

  // Step 7: Coverage by window
  println("\nStep 7: Coverage by window (chr1, width=100, step=50)...")
  let coverage = @bio.coverage_by_window(gr, "chr1", 1000, 100, 50)
  println("  Windows: \{coverage.length()}")
  for entry in coverage {
    let (start, end, count) = entry
    println("  chr1: \{start}-\{end} -> coverage=\{count}")
  }

  // Step 8: Bin genome
  println("\nStep 8: Bin genome (width=200)...")
  let bins = @bio.bin_genome(gr, seqlengths, 200)
  println("  Bins: \{bins.length()}")
  for entry in bins {
    let (chrom, start, end, count, density) = entry
    if count > 0 {
      println("  \{chrom}: \{start}-\{end} -> count=\{count}, density=\{density.to_string()}")
    }
  }

  // Step 9: Promoters
  println("\nStep 9: Creating promoter regions (upstream=200, downstream=50)...")
  let promoters = @bio.granges_promoters(gr, 200, 50)
  println("  Promoter ranges: \{promoters.seqnames.length()}")
  for i = 0; i < promoters.seqnames.length(); i = i + 1 {
    println("  \{promoters.seqnames[i]}: \{promoters.starts[i]}-\{promoters.ends[i]} (\{@bio.strand_to_char(promoters.strands[i])})")
  }

  // Step 10: Gaps
  println("\nStep 10: Finding gaps in chr1 (1 to 1000)...")
  let gr_chr1 = @bio.granges(
    ["chr1", "chr1"],
    [(100, 200), (400, 600)],
    [@bio.strand_plus(), @bio.strand_plus()],
  )
  let gaps = @bio.granges_gaps(gr_chr1, start_limit=1, end_limit=1000)
  println("  Gap ranges: \{gaps.seqnames.length()}")
  for i = 0; i < gaps.seqnames.length(); i = i + 1 {
    println("  \{gaps.seqnames[i]}: \{gaps.starts[i]}-\{gaps.ends[i]}")
  }

  // Step 11: Subtract
  println("\nStep 11: Subtracting ranges...")
  let gr_main = @bio.granges(["chr1"], [(100, 500)], [@bio.strand_plus()])
  let gr_sub = @bio.granges(["chr1"], [(200, 300)], [@bio.strand_minus()])
  let subtracted = @bio.granges_subtract(gr_main, gr_sub)
  println("  After subtraction: \{subtracted.seqnames.length()} ranges")
  for i = 0; i < subtracted.seqnames.length(); i = i + 1 {
    println("  \{subtracted.seqnames[i]}: \{subtracted.starts[i]}-\{subtracted.ends[i]}")
  }

  // Step 12: Resize advanced
  println("\nStep 12: Resizing ranges (width=100, fix=center)...")
  let resized = @bio.granges_resize_advanced(gr, 100, "center")
  println("  Resized ranges: \{resized.seqnames.length()}")
  for i = 0; i < resized.seqnames.length(); i = i + 1 {
    println("  \{resized.seqnames[i]}: \{resized.starts[i]}-\{resized.ends[i]} (width=\{resized.widths[i]})")
  }

  println("\n=== GenomicRanges Advanced Demo Complete ===")
}
