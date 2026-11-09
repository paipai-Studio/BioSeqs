// Bioframe third-round: Hi-C interval transforms (items 1-5)
// arrows, stripes, dots, insulation_score, boundary_calling

///|
pub struct Bf3Interval {
  chrom : String
  start : Int32
  end_ : Int32
  name : String
  score : Double
  strand : String
}

///|
pub fn bf3_interval(chrom : String, start : Int32, end_ : Int32) -> Bf3Interval {
  { chrom, start, end_, name: "", score: 0.0, strand: ".", }
}

// --- Arrows: Hi-C pixel → arrow interval representation ---

///|
pub struct Bf3Arrow {
  chrom : String
  start : Int32 // start of arrow (bin1)
  end_ : Int32 // end of arrow (bin2)
  direction : String // "up" or "down"
  score : Double
}

///|
pub fn bf3_arrows(
  pixels : Array[(String, Int32, Int32, Double)],
  bin_size : Int32,
) -> Array[Bf3Arrow] {
  // Convert Hi-C pixels (chrom, bin1_pos, bin2_pos, count) to arrow intervals
  // Each pixel in upper triangle becomes an arrow from bin1 to bin2
  let arrows : Array[Bf3Arrow] = []
  for i in 0..<pixels.length() {
    let (chrom, bin1, bin2, count) = pixels[i]
    if bin1 <= bin2 {
      arrows.push({
        chrom,
        start: bin1,
        end_: bin2 + bin_size,
        direction: "up",
        score: count,
      })
    } else {
      arrows.push({
        chrom,
        start: bin2,
        end_: bin1 + bin_size,
        direction: "down",
        score: count,
      })
    }
  }
  arrows
}

// --- Stripes: Hi-C pixel → row/column stripe intervals ---

///|
pub struct Bf3Stripe {
  chrom : String
  start : Int32
  end_ : Int32
  stripe_type : String // "row" or "column"
  anchor : Int32 // anchor position
  score : Double
}

///|
pub fn bf3_stripes(
  pixels : Array[(String, Int32, Int32, Double)],
  bin_size : Int32,
  stripe_type : String,
) -> Array[Bf3Stripe] {
  // Extract row or column stripes from Hi-C pixels
  let stripes : Array[Bf3Stripe] = []
  for i in 0..<pixels.length() {
    let (chrom, bin1, bin2, count) = pixels[i]
    if stripe_type == "row" {
      // Row stripe: fixed bin1, varying bin2
      stripes.push({
        chrom,
        start: bin2,
        end_: bin2 + bin_size,
        stripe_type: "row",
        anchor: bin1,
        score: count,
      })
    } else {
      // Column stripe: fixed bin2, varying bin1
      stripes.push({
        chrom,
        start: bin1,
        end_: bin1 + bin_size,
        stripe_type: "column",
        anchor: bin2,
        score: count,
      })
    }
  }
  stripes
}

// --- Dots: significant contact points as point intervals ---

///|
pub struct Bf3Dot {
  chrom : String
  pos : Int32 // position of the dot
  score : Double // significance score
  is_significant : Bool
}

///|
pub fn bf3_dots(
  pixels : Array[(String, Int32, Int32, Double)],
  threshold : Double,
) -> Array[Bf3Dot] {
  // Extract significant contact dots (pixels above threshold)
  let dots : Array[Bf3Dot] = []
  for i in 0..<pixels.length() {
    let (chrom, bin1, bin2, count) = pixels[i]
    if count >= threshold {
      let midpoint = (bin1 + bin2) / 2
      dots.push({ chrom, pos: midpoint, score: count, is_significant: true, })
    }
  }
  dots
}

// --- Insulation score ---

///|
pub fn bf3_insulation_score(
  pixels : Array[(Int32, Int32, Double)],
  chrom_size : Int32,
  bin_size : Int32,
  window_bins : Int32,
) -> Array[(Int32, Double)] {
  // Calculate insulation score along the diagonal
  // pixels: (bin1_idx, bin2_idx, count) — upper triangle
  // For each position, sum contacts in a window around the diagonal
  let num_bins = chrom_size / bin_size
  let scores : Array[(Int32, Double)] = []
  // Build contact matrix (sparse)
  let mut max_idx = 0
  for i in 0..<pixels.length() {
    let (b1, b2, _) = pixels[i]
    if b2 > max_idx {
      max_idx = b2
    }
    if b1 > max_idx {
      max_idx = b1
    }
  }
  let _max_bin = if max_idx > num_bins { max_idx } else { num_bins }
  for bin_idx in 0..<num_bins {
    let w_start = bin_idx - window_bins
    let w_end = bin_idx + window_bins
    let lo = if w_start < 0 { 0 } else { w_start }
    let hi = if w_end > num_bins { num_bins } else { w_end }
    let mut sum_contacts = 0.0
    let mut count_pixels = 0
    for i in 0..<pixels.length() {
      let (b1, b2, count) = pixels[i]
      // Check if pixel is in the window around the diagonal at bin_idx
      if b1 >= lo && b1 < hi && b2 >= lo && b2 < hi {
        // Pixel is within the window
        if (b1 < bin_idx && b2 >= bin_idx) || (b2 < bin_idx && b1 >= bin_idx) {
          // Pixel crosses the boundary at bin_idx — this is what insulation measures
          sum_contacts = sum_contacts + count
          count_pixels = count_pixels + 1
        }
      }
    }
    // Insulation score = mean contacts crossing the boundary
    let score = if count_pixels > 0 {
      sum_contacts / count_pixels.to_double()
    } else {
      0.0
    }
    scores.push((bin_idx * bin_size, score))
  }
  scores
}

// --- Boundary calling from insulation score ---

///|
pub struct Bf3Boundary {
  pos : Int32
  strength : Double // how deep the local minimum is
  is_boundary : Bool
}

///|
pub fn bf3_boundary_calling(
  insulation : Array[(Int32, Double)],
  min_strength : Double,
  window : Int32,
) -> Array[Bf3Boundary] {
  // Find local minima in insulation score (boundaries = deep minima)
  let boundaries : Array[Bf3Boundary] = []
  for i in 0..<insulation.length() {
    let (pos, score) = insulation[i]
    let lo = if i - window < 0 { 0 } else { i - window }
    let hi = if i + window >= insulation.length() {
      insulation.length() - 1
    } else {
      i + window
    }
    let mut is_min = true
    let mut max_neighbor = score
    for j in lo..<=hi {
      if j != i {
        let (_, neighbor_score) = insulation[j]
        if neighbor_score < score {
          is_min = false
        }
        if neighbor_score > max_neighbor {
          max_neighbor = neighbor_score
        }
      }
    }
    if is_min {
      let strength = max_neighbor - score
      if strength >= min_strength {
        boundaries.push({ pos, strength, is_boundary: true, })
      }
    }
  }
  boundaries
}

///|
pub fn bf3_boundaries_to_intervals(
  boundaries : Array[Bf3Boundary],
  chrom : String,
  bin_size : Int32,
) -> Array[Bf3Interval] {
  let intervals : Array[Bf3Interval] = []
  for i in 0..<boundaries.length() {
    let b = boundaries[i]
    intervals.push({
      chrom,
      start: b.pos,
      end_: b.pos + bin_size,
      name: "boundary",
      score: b.strength,
      strand: ".",
    })
  }
  intervals
}
