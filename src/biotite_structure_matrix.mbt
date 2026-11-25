///|
/// biotite_structure_matrix.mbt — distance / angle / dihedral matrices
/// and sequence identity, following biotite.structure.matrix.

///|
/// Full pairwise coordinate distance matrix.
pub fn biotite_get_distance_matrix(arr : AtomArray) -> Array[Array[Double]] {
  let n = arr.len()
  let result : Array[Array[Double]] = Array(capacity=n)
  for i in 0..<n {
    let row : Array[Double] = Array(capacity=n)
    for j in 0..<n {
      row.push(arr.get_coord(i).distance(arr.get_coord(j)))
    }
    result.push(row)
  }
  result
}

///|
/// Angle (degrees) for each (i, j, k) triple, with j as the vertex.
pub fn biotite_get_angle_matrix(
  arr : AtomArray,
  triples : Array[(Int32, Int32, Int32)],
) -> Array[Double] {
  let result : Array[Double] = Array(capacity=triples.length())
  for (i, j, k) in triples {
    if i >= 0 &&
      i < arr.len() &&
      j >= 0 &&
      j < arr.len() &&
      k >= 0 &&
      k < arr.len() {
      let a = arr.get_coord(i).sub(arr.get_coord(j))
      let b = arr.get_coord(k).sub(arr.get_coord(j))
      result.push(a.angle_degrees(b))
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
/// Dihedral angle (degrees) for each (i, j, k, l) quartet.
pub fn biotite_get_dihedral_matrix(
  arr : AtomArray,
  quartets : Array[Array[Int32]],
) -> Array[Double] {
  let result : Array[Double] = Array(capacity=quartets.length())
  for q in quartets {
    if q.length() < 4 {
      result.push(0.0)
      continue
    }
    let i = q[0]
    let j = q[1]
    let k = q[2]
    let l = q[3]
    if i >= 0 &&
      i < arr.len() &&
      j >= 0 &&
      j < arr.len() &&
      k >= 0 &&
      k < arr.len() &&
      l >= 0 &&
      l < arr.len() {
      result.push(
        biotite_dihedral_degrees(
          arr.get_coord(i),
          arr.get_coord(j),
          arr.get_coord(k),
          arr.get_coord(l),
        ),
      )
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
/// Fraction of identical one-letter residues between two structures.
/// Chain separators emitted by AtomArray::to_sequence are ignored.
pub fn biotite_get_sequence_identity(a : AtomArray, b : AtomArray) -> Double {
  let s1 = biotite_strip_chain_separators(a.to_sequence())
  let s2 = biotite_strip_chain_separators(b.to_sequence())
  let n = if s1.length() < s2.length() { s1.length() } else { s2.length() }
  if n == 0 {
    return 0.0
  }
  let mut same = 0
  for i in 0..<n {
    if s1[i] == s2[i] {
      same = same + 1
    }
  }
  same.to_double() / n.to_double()
}

///|
fn biotite_strip_chain_separators(s : String) -> String {
  let sb = StringBuilder()
  for ch in s {
    if ch != '/' {
      let _ = sb.write_char(ch)
    }
  }
  sb.to_string()
}

///|
/// Dihedral angle in degrees, using the atan2 convention.
pub fn biotite_dihedral_degrees(
  p1 : Vector3,
  p2 : Vector3,
  p3 : Vector3,
  p4 : Vector3,
) -> Double {
  let b1 = p2.sub(p1)
  let b2 = p3.sub(p2)
  let b3 = p4.sub(p3)
  let n1 = b1.cross(b2)
  let n2 = b2.cross(b3)
  let b2n = b2.normalize()
  let x = n1.dot(n2)
  let y = n1.cross(n2).dot(b2n)
  @math.atan2(y, x) * 180.0 / 3.141592653589793
}
