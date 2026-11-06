/// biopython_chaos_game.mbt — Chaos Game Representation for DNA sequences
/// Inspired by Biopython Bio.Graphics and chaos game representation

///|
/// Chaos Game Representation (CGR): map a DNA sequence to 2D coordinates.
/// Each nucleotide maps to a corner of a square:
/// A -> (0,0), C -> (0,1), G -> (1,1), T -> (1,0)
/// Each position: new_point = (old_point + corner) / 2
pub fn biopython_chaos_game_representation(
  sequence : String,
  n_points~ : Int32 = 0,
) -> Array[(Double, Double)] {
  let n = if n_points > 0 { n_points } else { sequence.length() }
  let points : Array[(Double, Double)] = Array::make(n, (0.0, 0.0))
  let mut x = 0.5
  let mut y = 0.5
  for i = 0; i < n; i = i + 1 {
    let ch = sequence[i].to_int()
    // Get corner coordinates
    let (cx, cy) = match ch {
      65 | 97 => (0.0, 0.0) // A
      67 | 99 => (0.0, 1.0) // C
      71 | 103 => (1.0, 1.0) // G
      84 | 116 => (1.0, 0.0) // T
      _ => (0.5, 0.5) // Unknown
    }
    // CGR formula: new = (old + corner) / 2
    x = (x + cx) / 2.0
    y = (y + cy) / 2.0
    points[i] = (x, y)
  }
  points
}

///|
/// CGR frequency matrix: divide the unit square into a k x k grid
/// and count how many CGR points fall in each cell.
pub fn biopython_cgr_frequency_matrix(
  sequence : String,
  k : Int32,
) -> Array[Array[Int32]] {
  let points = biopython_chaos_game_representation(sequence)
  let matrix : Array[Array[Int32]] = Array::make(k, Array::make(k, 0))
  for i = 0; i < points.length(); i = i + 1 {
    let (px, py) = points[i]
    let row = (py * k.to_double()).floor().to_int()
    let col = (px * k.to_double()).floor().to_int()
    let r = if row >= k { k - 1 } else if row < 0 { 0 } else { row }
    let c = if col >= k { k - 1 } else if col < 0 { 0 } else { col }
    matrix[r][c] = matrix[r][c] + 1
  }
  matrix
}

///|
/// USM (Universal Sequence Map): CGR with higher-order mapping.
/// Maps k-mers instead of individual nucleotides.
pub fn biopython_universal_sequence_map(
  sequence : String,
  kmer_size~ : Int32 = 2,
) -> Array[(Double, Double)] {
  let n = sequence.length()
  if n < kmer_size {
    return []
  }
  let n_kmers = n - kmer_size + 1
  let points : Array[(Double, Double)] = Array::make(n_kmers, (0.0, 0.0))
  let mut x = 0.5
  let mut y = 0.5
  for i = 0; i < n_kmers; i = i + 1 {
    // Process k-mer
    for j = 0; j < kmer_size; j = j + 1 {
      let ch = sequence[i + j].to_int()
      let (cx, cy) = match ch {
        65 | 97 => (0.0, 0.0)
        67 | 99 => (0.0, 1.0)
        71 | 103 => (1.0, 1.0)
        84 | 116 => (1.0, 0.0)
        _ => (0.5, 0.5)
      }
      x = (x + cx) / 2.0
      y = (y + cy) / 2.0
    }
    points[i] = (x, y)
  }
  points
}
