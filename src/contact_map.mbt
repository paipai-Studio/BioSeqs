///|
/// Distance matrix and contact map analysis.
///
/// Biotite provides `distance()` and `contact_map()` functions for
/// pairwise atom/residue distance analysis.

///|
/// Compute a full distance matrix between two sets of coordinates.
///
/// Returns an n×m matrix where element [i][j] = |coords1[i] - coords2[j]|.
pub fn coordinate_distance_matrix(
  coords1 : Array[Vector3],
  coords2 : Array[Vector3],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for c1 in coords1 {
    let row : Array[Double] = []
    for c2 in coords2 {
      row.push(c1.distance(c2))
    }
    result.push(row)
  }
  result
}

///|
/// Compute a self-distance matrix (coords vs coords).
///
/// Only returns the upper triangle (i < j) to avoid redundancy.
/// The result is a flat array of distances.
pub fn self_distance_matrix(coords : Array[Vector3]) -> Array[Double] {
  let n = coords.length()
  let result : Array[Double] = []
  for i in 0..<n {
    for j in (i + 1)..<n {
      result.push(coords[i].distance(coords[j]))
    }
  }
  result
}

///|
/// Generate a binary contact map.
///
/// Two atoms are "in contact" if their distance is <= threshold.
pub fn contact_map(
  coords : Array[Vector3],
  threshold~ : Double = 8.0,
) -> Array[Array[Bool]] {
  let n = coords.length()
  let result : Array[Array[Bool]] = []
  for i in 0..<n {
    let row = Array::make(n, false)
    for j in 0..<n {
      if i != j {
        row[j] = coords[i].distance(coords[j]) <= threshold
      }
    }
    result.push(row)
  }
  result
}

///|
/// Generate a residue-level contact map from C-alpha coordinates.
///
/// Two residues are in contact if their C-alpha distance is <= threshold.
/// Default threshold is 8.0 Å (commonly used for residue contacts).
pub fn residue_contact_map(
  ca_coords : Array[Vector3],
  threshold~ : Double = 8.0,
) -> Array[Array[Bool]] {
  contact_map(ca_coords, threshold~)
}

///|
/// Count total contacts in a contact map.
pub fn count_contacts(map : Array[Array[Bool]]) -> Int32 {
  let mut count = 0
  let n = map.length()
  for i in 0..<n {
    for j in (i + 1)..<n {
      if i < map.length() && j < map[i].length() && map[i][j] {
        count = count + 1
      }
    }
  }
  count
}

///|
/// Find contacts between two different sets of atoms.
///
/// Returns pairs (i, j) where coords1[i] is within threshold of coords2[j].
pub fn inter_contact_pairs(
  coords1 : Array[Vector3],
  coords2 : Array[Vector3],
  threshold~ : Double = 5.0,
) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = []
  for i in 0..<coords1.length() {
    for j in 0..<coords2.length() {
      if coords1[i].distance(coords2[j]) <= threshold {
        result.push((i, j))
      }
    }
  }
  result
}

///|
/// Calculate contact order (a measure of folding rate).
///
/// Contact order = average sequence separation of contacting residues
/// divided by chain length.
pub fn contact_order(
  ca_coords : Array[Vector3],
  threshold~ : Double = 8.0,
) -> Double {
  let n = ca_coords.length()
  if n < 2 {
    return 0.0
  }
  let mut sum_sep = 0.0
  let mut count = 0
  for i in 0..<n {
    for j in (i + 1)..<n {
      if ca_coords[i].distance(ca_coords[j]) <= threshold {
        sum_sep = sum_sep + (j - i).to_double()
        count = count + 1
      }
    }
  }
  if count == 0 {
    return 0.0
  }
  // Contact order = average separation / chain length
  sum_sep / count.to_double() / n.to_double()
}

///|
/// Long-range order parameter.
///
/// Fraction of contacts between residues separated by > 12 in sequence.
pub fn long_range_order(
  ca_coords : Array[Vector3],
  threshold~ : Double = 8.0,
  min_sep~ : Int32 = 12,
) -> Double {
  let n = ca_coords.length()
  if n < 2 {
    return 0.0
  }
  let mut long_contacts = 0
  let mut total_contacts = 0
  for i in 0..<n {
    for j in (i + 1)..<n {
      if ca_coords[i].distance(ca_coords[j]) <= threshold {
        total_contacts = total_contacts + 1
        if j - i > min_sep {
          long_contacts = long_contacts + 1
        }
      }
    }
  }
  if total_contacts == 0 {
    return 0.0
  }
  long_contacts.to_double() / total_contacts.to_double()
}
