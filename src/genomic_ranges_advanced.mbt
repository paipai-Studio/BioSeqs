///|
/// GenomicRanges Advanced - Tiling and Windowing Operations.
/// Inspired by Bioconductor GenomicRanges package.
///
/// Provides:
/// - tile: Break each range into non-overlapping sub-ranges
/// - sliding_windows: Generate overlapping windows of specified width and step
/// - tile_genome: Tile entire genome by chromosome lengths
/// - tile_by_ranges: Tile each range individually
/// - coverage_by_window: Compute coverage in sliding windows

///|
/// Tile a GRanges object into non-overlapping sub-ranges.
/// Each input range is divided into tiles of approximately equal width.
///
/// Parameters:
/// - gr: Input GRanges object
/// - n: Number of tiles per range (if specified, width is computed)
/// - width: Width of each tile (if specified, n is computed)
///
/// If both n and width are specified, width takes precedence.
/// Returns a new GRanges with the tiled sub-ranges.
pub fn granges_tile(
  gr : GRanges,
  n? : Int = 0,
  width? : Int = 0,
) -> GRanges {
  let n_ranges = gr.seqnames.length()
  let result_seqnames : Array[String] = Array::new()
  let result_ranges : Array[(Int, Int)] = Array::new()
  let result_strands : Array[Strand] = Array::new()

  for i = 0; i < n_ranges; i = i + 1 {
    let seqname = gr.seqnames[i]
    let start = gr.starts[i]
    let end = gr.ends[i]
    let strand = gr.strands[i]
    let range_width = end - start + 1

    if width > 0 {
      // Tile by fixed width
      let mut pos = start
      while pos <= end {
        let tile_end = if pos + width - 1 < end {
          pos + width - 1
        } else {
          end
        }
        result_seqnames.push(seqname)
        result_ranges.push((pos, tile_end))
        result_strands.push(strand)
        pos = pos + width
      }
    } else if n > 0 {
      // Tile by number of tiles
      let tile_width = (range_width + n - 1) / n // ceiling division
      let mut pos = start
      let mut tile_count = 0
      while pos <= end && tile_count < n {
        let tile_end = if pos + tile_width - 1 < end {
          pos + tile_width - 1
        } else {
          end
        }
        result_seqnames.push(seqname)
        result_ranges.push((pos, tile_end))
        result_strands.push(strand)
        pos = pos + tile_width
        tile_count = tile_count + 1
      }
    } else {
      // Default: single tile (the range itself)
      result_seqnames.push(seqname)
      result_ranges.push((start, end))
      result_strands.push(strand)
    }
  }

  granges(result_seqnames, result_ranges, result_strands)
}

///|
/// Generate sliding windows over each range in a GRanges object.
/// Windows overlap by (width - step) positions.
///
/// Parameters:
/// - gr: Input GRanges object
/// - width: Width of each window
/// - step: Step size between consecutive window starts
///
/// Returns a new GRanges with the sliding windows.
pub fn granges_sliding_windows(
  gr : GRanges,
  width : Int,
  step : Int,
) -> GRanges {
  let n_ranges = gr.seqnames.length()
  let result_seqnames : Array[String] = Array::new()
  let result_ranges : Array[(Int, Int)] = Array::new()
  let result_strands : Array[Strand] = Array::new()

  for i = 0; i < n_ranges; i = i + 1 {
    let seqname = gr.seqnames[i]
    let start = gr.starts[i]
    let end = gr.ends[i]
    let strand = gr.strands[i]

    let mut pos = start
    while pos <= end {
      let win_end = if pos + width - 1 < end {
        pos + width - 1
      } else {
        end
      }
      // Only add windows that are at least half the width
      let actual_width = win_end - pos + 1
      if actual_width >= width / 2 || pos == start {
        result_seqnames.push(seqname)
        result_ranges.push((pos, win_end))
        result_strands.push(strand)
      }
      pos = pos + step
    }
  }

  granges(result_seqnames, result_ranges, result_strands)
}

///|
/// Tile an entire genome by chromosome.
/// Creates non-overlapping tiles spanning each chromosome from 1 to the
/// specified length.
///
/// Parameters:
/// - seqlengths: Array of (chromosome_name, length) tuples
/// - width: Width of each tile
///
/// Returns a GRanges with tiles covering the entire genome.
pub fn tile_genome(
  seqlengths : Array[(String, Int)],
  width : Int,
) -> GRanges {
  let result_seqnames : Array[String] = Array::new()
  let result_ranges : Array[(Int, Int)] = Array::new()
  let result_strands : Array[Strand] = Array::new()

  for pair in seqlengths {
    let (chrom, length) = pair
    let mut pos = 1
    while pos <= length {
      let tile_end = if pos + width - 1 < length {
        pos + width - 1
      } else {
        length
      }
      result_seqnames.push(chrom)
      result_ranges.push((pos, tile_end))
      result_strands.push(Strand::Star)
      pos = pos + width
    }
  }

  granges(result_seqnames, result_ranges, result_strands)
}

///|
/// Tile each range into a specified number of equal-width sub-ranges.
/// Unlike granges_tile, this guarantees exactly n tiles per range,
/// with the last tile potentially being shorter.
///
/// Parameters:
/// - gr: Input GRanges object
/// - n: Number of tiles per range
///
/// Returns a new GRanges with exactly n tiles per input range.
pub fn tile_by_ranges(gr : GRanges, n : Int) -> GRanges {
  let n_ranges = gr.seqnames.length()
  let result_seqnames : Array[String] = Array::new()
  let result_ranges : Array[(Int, Int)] = Array::new()
  let result_strands : Array[Strand] = Array::new()

  for i = 0; i < n_ranges; i = i + 1 {
    let start = gr.starts[i]
    let end = gr.ends[i]
    let range_width = end - start + 1
    let tile_width = range_width / n
    let remainder = range_width % n

    let mut pos = start
    for k = 0; k < n; k = k + 1 {
      // Distribute remainder among first few tiles
      let extra = if k < remainder { 1 } else { 0 }
      let tw = tile_width + extra
      let tile_end = pos + tw - 1
      result_seqnames.push(gr.seqnames[i])
      result_ranges.push((pos, tile_end))
      result_strands.push(gr.strands[i])
      pos = pos + tw
    }
  }

  granges(result_seqnames, result_ranges, result_strands)
}

///|
/// Compute coverage in sliding windows.
/// For each window, compute the number of ranges that overlap it.
///
/// Parameters:
/// - gr: Input GRanges with ranges to compute coverage from
/// - seqname: Chromosome name to compute coverage on
/// - seqlength: Length of the chromosome
/// - width: Window width
/// - step: Step size between windows
///
/// Returns array of (window_start, window_end, coverage_count).
pub fn coverage_by_window(
  gr : GRanges,
  seqname : String,
  seqlength : Int,
  width : Int,
  step : Int,
) -> Array[(Int, Int, Int)] {
  let result : Array[(Int, Int, Int)] = Array::new()

  // Filter ranges for this chromosome
  let chrom_ranges : Array[(Int, Int)] = Array::new()
  for i = 0; i < gr.seqnames.length(); i = i + 1 {
    if gr.seqnames[i] == seqname {
      chrom_ranges.push((gr.starts[i], gr.ends[i]))
    }
  }

  let mut pos = 1
  while pos <= seqlength {
    let win_end = if pos + width - 1 < seqlength {
      pos + width - 1
    } else {
      seqlength
    }

    // Count overlapping ranges
    let mut count = 0
    for range in chrom_ranges {
      let (r_start, r_end) = range
      // Check overlap: r_start <= win_end && r_end >= pos
      if r_start <= win_end && r_end >= pos {
        count = count + 1
      }
    }

    result.push((pos, win_end, count))
    pos = pos + step
  }

  result
}

///|
/// Bin a GRanges object into fixed-width bins across the genome.
/// Similar to tile_genome but also computes per-bin statistics.
///
/// Parameters:
/// - gr: Input GRanges with data ranges
/// - seqlengths: Array of (chromosome_name, length) tuples
/// - bin_width: Width of each bin
///
/// Returns array of (chrom, bin_start, bin_end, count, density).
pub fn bin_genome(
  gr : GRanges,
  seqlengths : Array[(String, Int)],
  bin_width : Int,
) -> Array[(String, Int, Int, Int, Double)] {
  let result : Array[(String, Int, Int, Int, Double)] = Array::new()

  // Group ranges by chromosome
  let chrom_map : Map[String, Array[(Int, Int)]] = Map([], capacity=seqlengths.length())
  for i = 0; i < gr.seqnames.length(); i = i + 1 {
    let chrom = gr.seqnames[i]
    let ranges = chrom_map.get(chrom)
    match ranges {
      Some(arr) => arr.push((gr.starts[i], gr.ends[i]))
      None => {
        let arr : Array[(Int, Int)] = [(gr.starts[i], gr.ends[i])]
        chrom_map.set(chrom, arr)
      }
    }
  }

  let total_ranges = gr.seqnames.length()

  for pair in seqlengths {
    let (chrom, length) = pair
    let chrom_ranges = chrom_map.get(chrom).unwrap_or([])

    let mut pos = 1
    while pos <= length {
      let bin_end = if pos + bin_width - 1 < length {
        pos + bin_width - 1
      } else {
        length
      }

      let mut count = 0
      for range in chrom_ranges {
        let (r_start, r_end) = range
        if r_start <= bin_end && r_end >= pos {
          count = count + 1
        }
      }

      let density = if total_ranges > 0 {
        count.to_double() / total_ranges.to_double()
      } else {
        0.0
      }

      result.push((chrom, pos, bin_end, count, density))
      pos = pos + bin_width
    }
  }

  result
}

///|
/// Resize ranges to a specified width, fixing the start, end, or center.
/// Extends the resize() function with more options.
///
/// Parameters:
/// - gr: Input GRanges
/// - width: New width for all ranges
/// - fix: "start", "end", or "center"
pub fn granges_resize_advanced(
  gr : GRanges,
  width : Int,
  fix : String,
) -> GRanges {
  let n = gr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)

  for i = 0; i < n; i = i + 1 {
    let start = gr.starts[i]
    let end = gr.ends[i]
    match fix {
      "start" => {
        new_starts[i] = start
        new_ends[i] = start + width - 1
      }
      "end" => {
        new_starts[i] = end - width + 1
        new_ends[i] = end
      }
      "center" => {
        let center = (start + end) / 2
        new_starts[i] = center - width / 2
        new_ends[i] = new_starts[i] + width - 1
      }
      _ => {
        new_starts[i] = start
        new_ends[i] = start + width - 1
      }
    }
  }

  GRanges::{
    seqnames: gr.seqnames.copy(),
    starts: new_starts,
    ends: new_ends,
    strands: gr.strands.copy(),
    widths: Array::make(n, width),
  }
}

///|
/// Create promoters around each range.
/// Extends each range upstream and downstream by specified amounts,
/// considering strand direction.
///
/// Parameters:
/// - gr: Input GRanges (typically TSS or gene ranges)
/// - upstream: Number of bases to extend upstream
/// - downstream: Number of bases to extend downstream
pub fn granges_promoters(
  gr : GRanges,
  upstream : Int,
  downstream : Int,
) -> GRanges {
  let n = gr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)

  for i = 0; i < n; i = i + 1 {
    let start = gr.starts[i]
    let end = gr.ends[i]
    match gr.strands[i] {
      Strand::Minus => {
        // For minus strand, promoter is downstream of the end
        new_starts[i] = end - downstream + 1
        new_ends[i] = end + upstream
      }
      _ => {
        // For plus strand and star, promoter is upstream of the start
        new_starts[i] = start - upstream
        new_ends[i] = start + downstream - 1
      }
    }
  }

  let new_widths : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    new_widths[i] = new_ends[i] - new_starts[i] + 1
  }

  GRanges::{
    seqnames: gr.seqnames.copy(),
    starts: new_starts,
    ends: new_ends,
    strands: gr.strands.copy(),
    widths: new_widths,
  }
}

///|
/// Compute gaps between ranges within each chromosome.
/// Returns the gaps (regions not covered by any range) for each chromosome
/// present in the input.
///
/// Parameters:
/// - gr: Input GRanges
/// - start_limit: Minimum position (default: 1)
/// - end_limit: Maximum position (0 means no limit)
pub fn granges_gaps(
  gr : GRanges,
  start_limit? : Int = 1,
  end_limit? : Int = 0,
) -> GRanges {
  // Group ranges by chromosome
  let chrom_map : Map[String, Array[(Int, Int)]] = Map([], capacity=4)
  for i = 0; i < gr.seqnames.length(); i = i + 1 {
    let chrom = gr.seqnames[i]
    let ranges = chrom_map.get(chrom)
    match ranges {
      Some(arr) => arr.push((gr.starts[i], gr.ends[i]))
      None => {
        let arr : Array[(Int, Int)] = [(gr.starts[i], gr.ends[i])]
        chrom_map.set(chrom, arr)
      }
    }
  }

  let result_seqnames : Array[String] = Array::new()
  let result_ranges : Array[(Int, Int)] = Array::new()
  let result_strands : Array[Strand] = Array::new()

  // For each chromosome, sort ranges and find gaps
  let chrom_keys = chrom_map.keys().collect()
  let mut ck = 0
  while ck < chrom_keys.length() {
    let chrom = chrom_keys[ck]
    let ranges = chrom_map.get(chrom).unwrap_or([])

    // Sort by start position
    let sorted = sort_ranges(ranges)

    let mut prev_end = start_limit - 1
    let mut ri = 0
    while ri < sorted.length() {
      let (r_start, r_end) = sorted[ri]
      if r_start > prev_end + 1 {
        result_seqnames.push(chrom)
        result_ranges.push((prev_end + 1, r_start - 1))
        result_strands.push(Strand::Star)
      }
      if r_end > prev_end {
        prev_end = r_end
      }
      ri = ri + 1
    }

    // Gap after last range
    if end_limit == 0 || prev_end < end_limit {
      let gap_end = if end_limit > 0 { end_limit } else { prev_end + 1 }
      if gap_end > prev_end {
        result_seqnames.push(chrom)
        result_ranges.push((prev_end + 1, gap_end))
        result_strands.push(Strand::Star)
      }
    }
    ck = ck + 1
  }

  granges(result_seqnames, result_ranges, result_strands)
}

///|
/// Sort an array of (start, end) ranges by start position.
fn sort_ranges(ranges : Array[(Int, Int)]) -> Array[(Int, Int)] {
  let result = ranges.copy()
  let n = result.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - 1 - i {
      if result[j].0 > result[j + 1].0 {
        let tmp = result[j]
        result[j] = result[j + 1]
        result[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Subtract ranges: find regions in gr1 not covered by gr2.
pub fn granges_subtract(
  gr1 : GRanges,
  gr2 : GRanges,
) -> GRanges {
  let result_seqnames : Array[String] = Array::new()
  let result_ranges : Array[(Int, Int)] = Array::new()
  let result_strands : Array[Strand] = Array::new()

  for i = 0; i < gr1.seqnames.length(); i = i + 1 {
    let chrom = gr1.seqnames[i]
    let mut segments : Array[(Int, Int)] = [(gr1.starts[i], gr1.ends[i])]

    // Subtract each range from gr2 that overlaps
    for j = 0; j < gr2.seqnames.length(); j = j + 1 {
      if gr2.seqnames[j] == chrom {
        let new_segments : Array[(Int, Int)] = Array::new()
        for seg in segments {
          let (s, e) = seg
          let sub_s = gr2.starts[j]
          let sub_e = gr2.ends[j]
          if sub_e < s || sub_s > e {
            // No overlap
            new_segments.push((s, e))
          } else {
            // Overlap: may split into up to 2 segments
            if sub_s > s {
              new_segments.push((s, sub_s - 1))
            }
            if sub_e < e {
              new_segments.push((sub_e + 1, e))
            }
          }
        }
        segments = new_segments
      }
    }

    for seg in segments {
      result_seqnames.push(chrom)
      result_ranges.push(seg)
      result_strands.push(gr1.strands[i])
    }
  }

  granges(result_seqnames, result_ranges, result_strands)
}
