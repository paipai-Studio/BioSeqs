/// pysam_bam_gc_content.mbt — GC content, base composition, read length, MAPQ distribution
/// Items 14-17

///|
/// GC content result for a read.
pub struct PysamGCResult {
  qname : String
  gc_count : Int32
  total : Int32
  gc_fraction : Double
}

///|
/// Calculate GC content for a single sequence.
pub fn pysam_gc_content(seq : String) -> PysamGCResult {
  let mut gc = 0
  let mut total = 0
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i].to_int()
    if c == 71 || c == 67 || c == 103 || c == 99 { // G, C, g, c
      gc = gc + 1
      total = total + 1
    } else if c == 65 || c == 84 || c == 97 || c == 116 { // A, T, a, t
      total = total + 1
    }
  }
  let frac = if total > 0 { gc.to_double() / total.to_double() } else { 0.0 }
  { qname: "", gc_count: gc, total, gc_fraction: frac, }
}

///|
/// GC content distribution histogram across all reads.
pub fn pysam_gc_distribution(
  sam_text : String,
  bins~ : Int32 = 20,
) -> Array[(Int32, Int32)] {
  let lines = sam_text.split("\n").to_array()
  let hist : Array[Int32] = Array::make(bins + 1, 0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 10 {
      continue
    }
    let seq = fields[9].to_owned()
    let gc = pysam_gc_content(seq)
    let bin = (gc.gc_fraction * bins.to_double()).to_int()
    let bin_idx = if bin > bins { bins } else { bin }
    hist[bin_idx] = hist[bin_idx] + 1
  }
  let result : Array[(Int32, Int32)] = Array(capacity=0)
  for b = 0; b <= bins; b = b + 1 {
    result.push((b, hist[b]))
  }
  result
}

///|
/// Base composition matrix: position x base (A/C/G/T/N) counts.
pub fn pysam_base_composition_matrix(
  sam_text : String,
  max_pos~ : Int32 = 200,
) -> Array[Array[Int32]] {
  // 5 rows: A, C, G, T, N; max_pos columns
  let matrix : Array[Array[Int32]] = Array::make(5, Array::make(max_pos, 0))
  let lines = sam_text.split("\n").to_array()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 10 {
      continue
    }
    let seq = fields[9].to_owned()
    for p = 0; p < seq.length() && p < max_pos; p = p + 1 {
      let base_idx = match9_base(seq[p].to_int())
      matrix[base_idx][p] = matrix[base_idx][p] + 1
    }
  }
  matrix
}

///|
/// Map base to index 0-4.
fn match9_base(c : Int32) -> Int32 {
  match c {
    65 | 97 => 0 // A, a
    67 | 99 => 1 // C, c
    71 | 103 => 2 // G, g
    84 | 116 => 3 // T, t
    _ => 4 // N and others
  }
}

///|
/// Read length distribution.
pub struct PysamReadLenDist {
  histogram : Array[(Int32, Int32)] // (length, count)
  mean : Double
  n50 : Int32
  total_reads : Int32
}

///|
/// Calculate read length distribution.
pub fn pysam_read_length_dist(sam_text : String) -> PysamReadLenDist {
  let lines = sam_text.split("\n").to_array()
  let lengths : Array[Int32] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 10 {
      continue
    }
    let seq = fields[9].to_owned()
    if seq != "*" {
      lengths.push(seq.length())
    }
  }
  if lengths.length() == 0 {
    return { histogram: [], mean: 0.0, n50: 0, total_reads: 0, }
  }
  lengths.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  // Build histogram
  let hist : Array[(Int32, Int32)] = Array(capacity=0)
  let mut i = 0
  while i < lengths.length() {
    let len = lengths[i]
    let mut count = 0
    while i < lengths.length() && lengths[i] == len {
      count = count + 1
      i = i + 1
    }
    hist.push((len, count))
  }
  // Mean
  let mut total = 0
  for j = 0; j < lengths.length(); j = j + 1 {
    total = total + lengths[j]
  }
  let mean = total.to_double() / lengths.length().to_double()
  // N50: sort descending, find length at 50% of total bases
  let sorted_desc = lengths.copy()
  sorted_desc.sort_by(fn(a, b) {
    if a > b {
      -1
    } else if a < b {
      1
    } else {
      0
    }
  })
  let half = total / 2
  let mut cumsum = 0
  let mut n50 = 0
  for j = 0; j < sorted_desc.length(); j = j + 1 {
    cumsum = cumsum + sorted_desc[j]
    if cumsum >= half {
      n50 = sorted_desc[j]
      break
    }
  }
  { histogram: hist, mean, n50, total_reads: lengths.length(), }
}

///|
/// MAPQ distribution.
pub struct PysamMapqDist {
  histogram : Array[(Int32, Int32)] // (mapq, count)
  mean : Double
  median : Int32
  total : Int32
}

///|
/// Calculate MAPQ distribution.
pub fn pysam_mapq_dist(sam_text : String) -> PysamMapqDist {
  let lines = sam_text.split("\n").to_array()
  let mapqs : Array[Int32] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 5 {
      continue
    }
    let mapq = parse_int(fields[4].to_owned())
    mapqs.push(mapq)
  }
  if mapqs.length() == 0 {
    return { histogram: [], mean: 0.0, median: 0, total: 0, }
  }
  // Histogram (0-60)
  let hist : Array[Int32] = Array::make(61, 0)
  for i = 0; i < mapqs.length(); i = i + 1 {
    let q = if mapqs[i] > 60 { 60 } else { mapqs[i] }
    hist[q] = hist[q] + 1
  }
  let hist_result : Array[(Int32, Int32)] = Array(capacity=0)
  for q = 0; q <= 60; q = q + 1 {
    if hist[q] > 0 {
      hist_result.push((q, hist[q]))
    }
  }
  // Mean
  let mut total = 0
  for i = 0; i < mapqs.length(); i = i + 1 {
    total = total + mapqs[i]
  }
  let mean = total.to_double() / mapqs.length().to_double()
  // Median
  let sorted = mapqs.copy()
  sorted.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let median = sorted[sorted.length() / 2]
  { histogram: hist_result, mean, median, total: mapqs.length(), }
}
