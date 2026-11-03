/// Sequence gap operations.
///
/// scikit-bio equivalent: `Sequence.degap()`, `has_gaps()`,
/// `gap_vector()`, `gap_map()`
///
/// Explicit gap manipulation for aligned sequences.

///|
/// Remove all gap characters from a sequence.
/// Default gap character is '-'.
pub fn skbio_degap(seq : String, gap~ : Char = '-') -> String {
  let sb = StringBuilder::new()
  for c in seq {
    if c != gap {
      sb.write_char(c)
    }
  }
  sb.to_string()
}

///|
/// Check if a sequence contains any gap characters.
pub fn skbio_has_gaps(seq : String, gap~ : Char = '-') -> Bool {
  for c in seq {
    if c == gap {
      return true
    }
  }
  false
}

///|
/// Return a boolean vector indicating gap positions.
/// True = gap, False = non-gap
pub fn skbio_gap_vector(seq : String, gap~ : Char = '-') -> Array[Bool] {
  let result : Array[Bool] = Array(capacity=0)
  for c in seq {
    result.push(c == gap)
  }
  result
}

///|
/// Count the number of gap characters in a sequence.
pub fn skbio_count_gaps(seq : String, gap~ : Char = '-') -> Int32 {
  let mut count : Int32 = 0
  for c in seq {
    if c == gap {
      count = count + 1
    }
  }
  count
}

///|
/// Compute the fraction of gap characters in a sequence.
pub fn skbio_gap_fraction(seq : String, gap~ : Char = '-') -> Double {
  if seq.length() == 0 {
    return 0.0
  }
  skbio_count_gaps(seq, gap~).to_double() / seq.length().to_double()
}

///|
/// Map ungapped positions to gapped positions.
/// Returns an array where result[ungapped_index] = gapped_index.
/// Useful for transferring coordinates from ungapped to aligned sequence.
pub fn skbio_gap_map(seq : String, gap~ : Char = '-') -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  let mut gapped_pos : Int32 = 0
  for c in seq {
    if c != gap {
      result.push(gapped_pos)
    }
    gapped_pos = gapped_pos + 1
  }
  result
}

///|
/// Map gapped positions to ungapped positions.
/// Returns an array where result[gapped_index] = ungapped_index, or -1 if gap.
pub fn skbio_gap_map_inverse(seq : String, gap~ : Char = '-') -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  let mut ungapped_pos : Int32 = 0
  for c in seq {
    if c == gap {
      result.push(-1)
    } else {
      result.push(ungapped_pos)
      ungapped_pos = ungapped_pos + 1
    }
  }
  result
}

///|
/// Find all gap runs (contiguous gap segments).
/// Returns array of (start, end) pairs (0-based, end exclusive).
pub fn skbio_find_gap_runs(
  seq : String,
  gap~ : Char = '-',
) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = Array(capacity=0)
  let mut in_gap = false
  let mut gap_start : Int32 = 0
  let mut pos : Int32 = 0
  for c in seq {
    if c == gap {
      if !in_gap {
        gap_start = pos
        in_gap = true
      }
    } else if in_gap {
      result.push((gap_start, pos))
      in_gap = false
    }
    pos = pos + 1
  }
  if in_gap {
    result.push((gap_start, pos))
  }
  result
}

///|
/// Insert gaps into a sequence at specified positions.
/// positions: sorted array of 0-based positions where gaps should be inserted
pub fn skbio_insert_gaps(
  seq : String,
  positions : Array[Int32],
  gap~ : Char = '-',
) -> String {
  let sb = StringBuilder::new()
  let pos_set : Set[Int32] = Set::from_array(positions)
  let mut i : Int32 = 0
  for c in seq {
    while pos_set.contains(i) {
      sb.write_char(gap)
      i = i + 1
    }
    sb.write_char(c)
    i = i + 1
  }
  while pos_set.contains(i) {
    sb.write_char(gap)
    i = i + 1
  }
  sb.to_string()
}
