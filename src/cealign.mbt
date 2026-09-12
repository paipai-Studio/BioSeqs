///|
/// Protein structure alignment by Combinatorial Extension (CE).
///
/// This module follows the core semantics of Biopython's
/// `Bio.PDB.cealign`: structures are represented by ordered CA or C4' guide
/// atoms, aligned fragment pairs are extended into monotonic paths, and the
/// longest paths are ranked by their post-superposition RMSD.

///|
/// Error raised for invalid CE parameters or structures.
pub suberror CeAlignError {
  CeAlignError(String)
}

///|
/// A guide atom used by CEAligner.
///
/// Proteins contribute CA atoms and nucleic acids contribute C4' atoms when
/// no CA atom is present in the residue.
pub struct CeGuideAtom {
  model_id : Int
  chain_id : Char
  residue_number : Int
  atom_name : String
  coordinate : Vector3
}

///|
/// Create a CE guide atom.
pub fn CeGuideAtom::new(
  model_id~ : Int,
  chain_id~ : Char,
  residue_number~ : Int,
  atom_name~ : String,
  coordinate~ : Vector3,
) -> CeGuideAtom {
  CeGuideAtom::{ model_id, chain_id, residue_number, atom_name, coordinate, }
}

///|
/// One CE alignment path expanded to individual guide atom indices.
pub struct CeAlignmentPath {
  reference_indices : Array[Int]
  mobile_indices : Array[Int]
  fragment_count : Int
  similarity : Double
  z_score : Double
}

///|
/// Number of aligned guide atom pairs.
pub fn CeAlignmentPath::aligned_length(self : CeAlignmentPath) -> Int {
  self.reference_indices.length()
}

///|
/// Return aligned index pairs.
pub fn CeAlignmentPath::aligned_pairs(
  self : CeAlignmentPath,
) -> Array[(Int, Int)] {
  let pairs : Array[(Int, Int)] = []
  let mut index = 0
  while index < self.reference_indices.length() {
    pairs.push((self.reference_indices[index], self.mobile_indices[index]))
    index = index + 1
  }
  pairs
}

///|
/// Return a defensive copy of the reference indices.
pub fn CeAlignmentPath::reference_path(self : CeAlignmentPath) -> Array[Int] {
  self.reference_indices.copy()
}

///|
/// Return a defensive copy of the mobile indices.
pub fn CeAlignmentPath::mobile_path(self : CeAlignmentPath) -> Array[Int] {
  self.mobile_indices.copy()
}

///|
/// Result of aligning a mobile structure onto a reference structure.
pub struct CeAlignmentResult {
  path : CeAlignmentPath
  rmsd : Double
  rotation : Array[Array[Double]]
  translation : Array[Double]
  reference_guide_count : Int
  mobile_guide_count : Int
  transformed : Bool
  optimization_applied : Bool
  structure : Structure
}

///|
/// Fraction of reference guide atoms included in the alignment.
pub fn CeAlignmentResult::reference_coverage(
  self : CeAlignmentResult,
) -> Double {
  if self.reference_guide_count == 0 {
    0.0
  } else {
    self.path.aligned_length().to_double() /
    self.reference_guide_count.to_double()
  }
}

///|
/// Fraction of mobile guide atoms included in the alignment.
pub fn CeAlignmentResult::mobile_coverage(self : CeAlignmentResult) -> Double {
  if self.mobile_guide_count == 0 {
    0.0
  } else {
    self.path.aligned_length().to_double() / self.mobile_guide_count.to_double()
  }
}

///|
/// Return a defensive copy of the rotation matrix.
pub fn CeAlignmentResult::rotation_matrix(
  self : CeAlignmentResult,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for row in self.rotation {
    result.push(row.copy())
  }
  result
}

///|
/// Return a defensive copy of the translation vector.
pub fn CeAlignmentResult::translation_vector(
  self : CeAlignmentResult,
) -> Array[Double] {
  self.translation.copy()
}

///|
/// Human-readable alignment summary.
pub fn CeAlignmentResult::summary(self : CeAlignmentResult) -> String {
  "CEAlignment(aligned=" +
  self.path.aligned_length().to_string() +
  ", fragments=" +
  self.path.fragment_count.to_string() +
  ", rmsd=" +
  self.rmsd.to_string() +
  ", z_score=" +
  self.path.z_score.to_string() +
  ", reference_coverage=" +
  self.reference_coverage().to_string() +
  ", mobile_coverage=" +
  self.mobile_coverage().to_string() +
  ")"
}

///|
/// Immutable CE aligner configuration and reference guide atom set.
pub struct CEAligner {
  window_size : Int
  max_gap : Int
  reference_atoms : Array[CeGuideAtom]
}

///|
/// Create a CEAligner.
pub fn CEAligner::new(
  window_size? : Int = 8,
  max_gap? : Int = 30,
) -> CEAligner raise CeAlignError {
  if window_size <= 0 {
    raise CeAlignError("window_size must be greater than zero")
  }
  if max_gap < 0 {
    raise CeAlignError("max_gap must be non-negative")
  }
  CEAligner::{ window_size, max_gap, reference_atoms: [], }
}

///|
/// Number of guide atoms in the current reference.
pub fn CEAligner::reference_length(self : CEAligner) -> Int {
  self.reference_atoms.length()
}

///|
/// Return a copy of the current reference coordinates.
pub fn CEAligner::reference_coordinates(self : CEAligner) -> Array[Vector3] {
  let coordinates : Array[Vector3] = []
  for atom in self.reference_atoms {
    coordinates.push(atom.coordinate)
  }
  coordinates
}

///|
/// Extract ordered CA/C4' guide atoms from a structure.
///
/// Ordering by model, chain, and residue number makes CE path construction
/// deterministic even when input chains are not stored in lexical order.
pub fn cealign_get_guide_atoms(
  structure : Structure,
) -> Array[CeGuideAtom] raise CeAlignError {
  let atoms : Array[CeGuideAtom] = []
  for model in structure.models {
    for chain in model.chains {
      for residue in chain.residues {
        let guide = match residue.find_atom("CA") {
          Some(atom) => Some(atom)
          None => residue.find_atom("C4'")
        }
        match guide {
          Some(atom) =>
            atoms.push(
              CeGuideAtom::new(
                model_id=model.id,
                chain_id=chain.id,
                residue_number=residue.resseq,
                atom_name=atom.name,
                coordinate=atom.coord,
              ),
            )
          None => ()
        }
      }
    }
  }
  atoms.sort_by(fn(left : CeGuideAtom, right : CeGuideAtom) -> Int {
    if left.model_id < right.model_id {
      -1
    } else if left.model_id > right.model_id {
      1
    } else if left.chain_id.to_int() < right.chain_id.to_int() {
      -1
    } else if left.chain_id.to_int() > right.chain_id.to_int() {
      1
    } else if left.residue_number < right.residue_number {
      -1
    } else if left.residue_number > right.residue_number {
      1
    } else {
      0
    }
  })
  if atoms.length() == 0 {
    raise CeAlignError(
      "Structure " + structure.id + " does not have CA or C4' guide atoms",
    )
  }
  atoms
}

///|
/// Set the immutable reference structure.
pub fn CEAligner::set_reference(
  self : CEAligner,
  structure : Structure,
) -> CEAligner raise CeAlignError {
  let atoms = cealign_get_guide_atoms(structure)
  if atoms.length() < self.window_size * 2 {
    raise CeAlignError(
      "Too few guide atoms in reference structure: " +
      atoms.length().to_string() +
      "; need at least " +
      (self.window_size * 2).to_string(),
    )
  }
  CEAligner::{ ..self, reference_atoms: atoms, }
}

///|
/// Calculate a symmetric intramolecular distance matrix.
pub fn cealign_distance_matrix(
  coordinates : Array[Vector3],
) -> Array[Array[Double]] {
  let length = coordinates.length()
  let matrix : Array[Array[Double]] = []
  let mut row = 0
  while row < length {
    matrix.push(Array::make(length, 0.0))
    row = row + 1
  }
  row = 0
  while row < length {
    let mut column = row
    while column < length {
      let distance = coordinates[row].distance(coordinates[column])
      matrix[row][column] = distance
      matrix[column][row] = distance
      column = column + 1
    }
    row = row + 1
  }
  matrix
}

///|
/// CE distance measure II for one aligned fragment pair.
///
/// Values are non-positive; zero represents identical internal geometry.
pub fn cealign_fragment_similarity(
  reference_distances : Array[Array[Double]],
  mobile_distances : Array[Array[Double]],
  reference_start : Int,
  mobile_start : Int,
  window_size : Int,
) -> Double {
  if window_size <= 1 {
    return 0.0
  }
  if window_size == 2 {
    return -(reference_distances[reference_start][reference_start + 1] -
    mobile_distances[mobile_start][mobile_start + 1]).abs()
  }
  let mut difference = 0.0
  let mut first = 0
  while first < window_size - 2 {
    let mut second = first + 2
    while second < window_size {
      difference = difference +
        (reference_distances[reference_start + first][reference_start + second] -
        mobile_distances[mobile_start + first][mobile_start + second]).abs()
      second = second + 1
    }
    first = first + 1
  }
  let term_count = (window_size - 1) * (window_size - 2) / 2
  -difference / term_count.to_double()
}

///|
/// Build the CE aligned-fragment similarity matrix.
pub fn cealign_similarity_matrix(
  reference : Array[Vector3],
  mobile : Array[Vector3],
  window_size : Int,
) -> Array[Array[Double]] raise CeAlignError {
  if window_size <= 0 {
    raise CeAlignError("window_size must be greater than zero")
  }
  if reference.length() < window_size || mobile.length() < window_size {
    raise CeAlignError("coordinate sets must contain at least one full window")
  }
  let reference_distances = cealign_distance_matrix(reference)
  let mobile_distances = cealign_distance_matrix(mobile)
  let row_count = reference.length() - window_size + 1
  let column_count = mobile.length() - window_size + 1
  let similarity : Array[Array[Double]] = []
  let mut row = 0
  while row < row_count {
    let values : Array[Double] = []
    let mut column = 0
    while column < column_count {
      values.push(
        cealign_fragment_similarity(
          reference_distances, mobile_distances, row, column, window_size,
        ),
      )
      column = column + 1
    }
    similarity.push(values)
    row = row + 1
  }
  similarity
}

///|
struct CePathCandidate {
  reference_starts : Array[Int]
  mobile_starts : Array[Int]
  fragment_count : Int
  similarity : Double
  z_score : Double
}

///|
fn cealign_extension_similarity(
  reference_distances : Array[Array[Double]],
  mobile_distances : Array[Array[Double]],
  reference_start_i : Int,
  mobile_start_i : Int,
  reference_start_j : Int,
  mobile_start_j : Int,
  window_size : Int,
) -> Double {
  let last = window_size - 1
  let mut difference = (reference_distances[reference_start_i][reference_start_j] -
    mobile_distances[mobile_start_i][mobile_start_j]).abs() +
    (reference_distances[reference_start_i + last][reference_start_j + last] -
    mobile_distances[mobile_start_i + last][mobile_start_j + last]).abs()
  let mut offset = 1
  while offset < last {
    difference = difference +
      (reference_distances[reference_start_i + offset][reference_start_j +
      last -
      offset] -
      mobile_distances[mobile_start_i + offset][mobile_start_j + last - offset]).abs()
    offset = offset + 1
  }
  -difference / window_size.to_double()
}

///|
fn cealign_normal_survival(z : Double) -> Double {
  if z <= 0.0 {
    return 0.5
  }
  let t = 1.0 / (1.0 + 0.2316419 * z)
  let density = 0.3989422804014327 * @math.exp(-z * z / 2.0)
  density *
  t *
  (
    0.319381530 +
    t *
    (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
  )
}

///|
fn cealign_inverse_survival(probability : Double) -> Double {
  if probability >= 0.5 {
    return 0.0
  }
  let p = if probability < 1.0e-300 { 1.0e-300 } else { probability }
  let eta = (-2.0 * @math.ln(p)).sqrt()
  eta -
  (2.515517 + 0.802853 * eta + 0.010328 * eta * eta) /
  (1.0 + 1.432788 * eta + 0.189269 * eta * eta + 0.001308 * eta * eta * eta)
}

///|
fn cealign_similarity_z_score(
  fragment_count : Int,
  similarity : Double,
) -> Double {
  if fragment_count < 1 {
    return 0.0
  }
  let averages = [
    2.54, 2.51, 2.72, 3.01, 3.31, 3.61, 3.90, 4.19, 4.47, 4.74, 4.99, 5.22, 5.46,
    5.70, 5.94, 6.13, 6.36, 6.52, 6.68, 6.91,
  ]
  let deviations = [
    1.33, 0.88, 0.73, 0.71, 0.74, 0.80, 0.86, 0.92, 0.98, 1.04, 1.08, 1.10, 1.15,
    1.19, 1.23, 1.25, 1.32, 1.34, 1.36, 1.45,
  ]
  let average = if fragment_count <= 20 {
    averages[fragment_count - 1]
  } else {
    0.209874 * fragment_count.to_double() + 2.944714
  }
  let deviation = if fragment_count <= 20 {
    deviations[fragment_count - 1]
  } else {
    0.039487 * fragment_count.to_double() + 0.675735
  }
  let positive_similarity = -similarity
  if positive_similarity > average {
    0.0
  } else {
    (average - positive_similarity) / deviation
  }
}

///|
fn cealign_gap_z_score(fragment_count : Int, gap_count : Int) -> Double {
  if fragment_count <= 1 {
    return 0.0
  }
  let averages = [
    0.00, 11.50, 23.32, 35.95, 49.02, 62.44, 76.28, 90.26, 104.86, 119.97, 134.86,
    150.54, 164.86, 179.57, 194.39, 209.38, 224.74, 238.96, 253.72, 270.79,
  ]
  let deviations = [
    0.00, 9.88, 14.34, 17.99, 21.10, 23.89, 26.55, 29.00, 31.11, 33.10, 35.02, 36.03,
    37.19, 38.82, 41.04, 43.35, 45.45, 48.41, 50.87, 52.27,
  ]
  let average = if fragment_count <= 20 {
    averages[fragment_count - 1]
  } else {
    14.949173 * fragment_count.to_double() - 14.581193
  }
  let deviation = if fragment_count <= 20 {
    deviations[fragment_count - 1]
  } else {
    2.045067 * fragment_count.to_double() + 13.191095
  }
  if gap_count.to_double() > average || deviation == 0.0 {
    0.0
  } else {
    (average - gap_count.to_double()) / deviation
  }
}

///|
fn cealign_z_score(
  window_size : Int,
  fragment_count : Int,
  similarity : Double,
  gap_count : Int,
) -> Double {
  if window_size != 8 {
    return 0.0
  }
  let similarity_z = cealign_similarity_z_score(fragment_count, similarity)
  let gap_z = cealign_gap_z_score(fragment_count, gap_count)
  let combined_probability = cealign_normal_survival(similarity_z) *
    cealign_normal_survival(gap_z)
  cealign_inverse_survival(combined_probability)
}

///|
fn cealign_store_candidate(
  candidates : Array[CePathCandidate],
  candidate : CePathCandidate,
) -> Unit {
  candidates.push(candidate)
  candidates.sort_by(fn(
    left : CePathCandidate,
    right : CePathCandidate,
  ) -> Int {
    if left.fragment_count > right.fragment_count {
      -1
    } else if left.fragment_count < right.fragment_count {
      1
    } else if left.similarity > right.similarity {
      -1
    } else if left.similarity < right.similarity {
      1
    } else {
      0
    }
  })
  if candidates.length() > 20 {
    ignore(candidates.pop())
  }
}

///|
fn cealign_expand_candidate(
  candidate : CePathCandidate,
  window_size : Int,
) -> CeAlignmentPath {
  let reference_indices : Array[Int] = []
  let mobile_indices : Array[Int] = []
  let mut fragment = 0
  while fragment < candidate.fragment_count {
    let mut offset = 0
    while offset < window_size {
      reference_indices.push(candidate.reference_starts[fragment] + offset)
      mobile_indices.push(candidate.mobile_starts[fragment] + offset)
      offset = offset + 1
    }
    fragment = fragment + 1
  }
  CeAlignmentPath::{
    reference_indices,
    mobile_indices,
    fragment_count: candidate.fragment_count,
    similarity: candidate.similarity,
    z_score: candidate.z_score,
  }
}

///|
/// Find up to 20 CE paths, ordered by fragment count and similarity.
pub fn cealign_find_paths(
  reference : Array[Vector3],
  mobile : Array[Vector3],
  window_size : Int,
  max_gap : Int,
) -> Array[CeAlignmentPath] raise CeAlignError {
  if window_size <= 0 {
    raise CeAlignError("window_size must be greater than zero")
  }
  if max_gap < 0 {
    raise CeAlignError("max_gap must be non-negative")
  }
  if reference.length() < window_size || mobile.length() < window_size {
    raise CeAlignError("coordinate sets must contain at least one full window")
  }
  let reference_distances = cealign_distance_matrix(reference)
  let mobile_distances = cealign_distance_matrix(mobile)
  let similarity = cealign_similarity_matrix(reference, mobile, window_size)
  let candidates : Array[CePathCandidate] = []
  let row_count = similarity.length()
  let column_count = similarity[0].length()
  let mut reference_start = 0
  while reference_start < row_count {
    let mut mobile_start = 0
    while mobile_start < column_count {
      if similarity[reference_start][mobile_start] > -3.0 {
        let reference_starts = [reference_start]
        let mobile_starts = [mobile_start]
        let mut fragment_count = 1
        let mut path_similarity = similarity[reference_start][mobile_start]
        let mut extending = true
        while extending {
          let previous_reference = reference_starts[fragment_count - 1]
          let previous_mobile = mobile_starts[fragment_count - 1]
          let mut best_reference = -1
          let mut best_mobile = -1
          let mut best_extension_similarity = -1000000.0
          let mut gap_slot = 0
          while gap_slot < max_gap * 2 + 1 {
            let mut next_reference = previous_reference + window_size
            let mut next_mobile = previous_mobile + window_size
            if (gap_slot + 1) % 2 == 0 {
              next_reference = next_reference + (gap_slot + 1) / 2
            } else {
              next_mobile = next_mobile + (gap_slot + 1) / 2
            }
            if next_reference < row_count &&
              next_mobile < column_count &&
              similarity[next_reference][next_mobile] > -3.0 {
              let mut extension_similarity = 0.0
              let mut existing = 0
              while existing < fragment_count {
                extension_similarity = extension_similarity +
                  cealign_extension_similarity(
                    reference_distances,
                    mobile_distances,
                    reference_starts[existing],
                    mobile_starts[existing],
                    next_reference,
                    next_mobile,
                    window_size,
                  )
                existing = existing + 1
              }
              extension_similarity = extension_similarity /
                fragment_count.to_double()
              if extension_similarity > -4.0 &&
                extension_similarity > best_extension_similarity {
                best_extension_similarity = extension_similarity
                best_reference = next_reference
                best_mobile = next_mobile
              }
            }
            gap_slot = gap_slot + 1
          }
          if best_reference < 0 {
            extending = false
          } else {
            let count = fragment_count.to_double()
            let current_terms = count + count * (count - 1.0) / 2.0
            let new_terms = count + 1.0 + count * (count + 1.0) / 2.0
            let new_similarity = (
                current_terms * path_similarity +
                count * best_extension_similarity +
                similarity[best_reference][best_mobile]
              ) /
              new_terms
            if new_similarity > -4.0 {
              reference_starts.push(best_reference)
              mobile_starts.push(best_mobile)
              fragment_count = fragment_count + 1
              path_similarity = new_similarity
            } else {
              extending = false
            }
          }
        }
        let mut gap_count = 0
        let mut fragment = 1
        while fragment < fragment_count {
          gap_count = gap_count +
            reference_starts[fragment] -
            reference_starts[fragment - 1] -
            1
          gap_count = gap_count +
            mobile_starts[fragment] -
            mobile_starts[fragment - 1] -
            1
          fragment = fragment + 1
        }
        cealign_store_candidate(candidates, CePathCandidate::{
          reference_starts,
          mobile_starts,
          fragment_count,
          similarity: path_similarity,
          z_score: cealign_z_score(
            window_size, fragment_count, path_similarity, gap_count,
          ),
        })
      }
      mobile_start = mobile_start + 1
    }
    reference_start = reference_start + 1
  }
  let paths : Array[CeAlignmentPath] = []
  for candidate in candidates {
    paths.push(cealign_expand_candidate(candidate, window_size))
  }
  paths
}

///|
fn cealign_fit_path(
  path : CeAlignmentPath,
  reference : Array[Vector3],
  mobile : Array[Vector3],
) -> QCPSuperimposer {
  let fixed : Array[QCPAtomCoordinate] = []
  let moving : Array[QCPAtomCoordinate] = []
  let mut index = 0
  while index < path.reference_indices.length() {
    let reference_coordinate = reference[path.reference_indices[index]]
    let mobile_coordinate = mobile[path.mobile_indices[index]]
    fixed.push(
      QCPAtomCoordinate::new(
        reference_coordinate.x,
        reference_coordinate.y,
        reference_coordinate.z,
      ),
    )
    moving.push(
      QCPAtomCoordinate::new(
        mobile_coordinate.x,
        mobile_coordinate.y,
        mobile_coordinate.z,
      ),
    )
    index = index + 1
  }
  QCPSuperimposer::new().set(fixed, moving)
}

///|
fn cealign_path_with_indices(
  path : CeAlignmentPath,
  reference_indices : Array[Int],
  mobile_indices : Array[Int],
) -> CeAlignmentPath {
  CeAlignmentPath::{
    ..path,
    reference_indices: reference_indices.copy(),
    mobile_indices: mobile_indices.copy(),
  }
}

///|
fn cealign_optimize_path(
  path : CeAlignmentPath,
  reference : Array[Vector3],
  mobile : Array[Vector3],
  window_size : Int,
  initial_fit : QCPSuperimposer,
) -> (CeAlignmentPath, QCPSuperimposer) {
  let reference_indices = path.reference_indices.copy()
  let mobile_indices = path.mobile_indices.copy()
  let mut best_fit = initial_fit
  let half_window = window_size / 2
  let mut side = 0
  while side < 2 {
    let indices = if side == 0 { reference_indices } else { mobile_indices }
    let coordinate_count = if side == 0 {
      reference.length()
    } else {
      mobile.length()
    }
    let mut index = 1
    while index + 1 < indices.length() {
      let center = indices[index]
      let left = indices[index - 1]
      let right = indices[index + 1]
      let mut minimum_shift = -half_window
      if left - center + 1 > minimum_shift {
        minimum_shift = left - center + 1
      }
      if -center > minimum_shift {
        minimum_shift = -center
      }
      let mut maximum_shift = half_window
      if right - center - 1 < maximum_shift {
        maximum_shift = right - center - 1
      }
      if coordinate_count - 1 - center < maximum_shift {
        maximum_shift = coordinate_count - 1 - center
      }
      let mut best_shift = 0
      let mut shift = minimum_shift
      while shift <= maximum_shift {
        if shift != 0 {
          indices[index] = center + shift
          let candidate_path = cealign_path_with_indices(
            path, reference_indices, mobile_indices,
          )
          let candidate_fit = cealign_fit_path(
            candidate_path, reference, mobile,
          )
          if candidate_fit.rmsd + 1.0e-12 < best_fit.rmsd {
            best_shift = shift
            best_fit = candidate_fit
          }
        }
        indices[index] = center
        shift = shift + 1
      }
      indices[index] = center + best_shift
      index = index + 1
    }
    side = side + 1
  }
  (cealign_path_with_indices(path, reference_indices, mobile_indices), best_fit)
}

///|
fn cealign_validate_transform(
  rotation : Array[Array[Double]],
  translation : Array[Double],
) -> Unit raise CeAlignError {
  if rotation.length() != 3 ||
    rotation[0].length() != 3 ||
    rotation[1].length() != 3 ||
    rotation[2].length() != 3 {
    raise CeAlignError("rotation must be a 3 x 3 matrix")
  }
  if translation.length() != 3 {
    raise CeAlignError("translation must contain three values")
  }
}

///|
fn cealign_transform_vector(
  coordinate : Vector3,
  rotation : Array[Array[Double]],
  translation : Array[Double],
) -> Vector3 {
  Vector3::new(
    rotation[0][0] * coordinate.x +
    rotation[0][1] * coordinate.y +
    rotation[0][2] * coordinate.z +
    translation[0],
    rotation[1][0] * coordinate.x +
    rotation[1][1] * coordinate.y +
    rotation[1][2] * coordinate.z +
    translation[1],
    rotation[2][0] * coordinate.x +
    rotation[2][1] * coordinate.y +
    rotation[2][2] * coordinate.z +
    translation[2],
  )
}

///|
/// Apply a rigid transform to every atom and return a new structure.
pub fn cealign_transform_structure(
  structure : Structure,
  rotation : Array[Array[Double]],
  translation : Array[Double],
) -> Structure raise CeAlignError {
  cealign_validate_transform(rotation, translation)
  let models : Array[Model] = []
  for model in structure.models {
    let chains : Array[Chain] = []
    for chain in model.chains {
      let residues : Array[Residue] = []
      for residue in chain.residues {
        let atoms : Array[Atom] = []
        for atom in residue.atoms {
          atoms.push(Atom::{
            ..atom,
            coord: cealign_transform_vector(atom.coord, rotation, translation),
          })
        }
        residues.push(Residue::{ ..residue, atoms, })
      }
      chains.push(Chain::{ ..chain, residues, })
    }
    models.push(Model::{ ..model, chains, })
  }
  Structure::{ ..structure, models, }
}

///|
/// Align a mobile structure onto the configured reference.
///
/// When `transform` is true, `result.structure` contains a transformed copy.
/// The input structure is never modified.
pub fn CEAligner::align(
  self : CEAligner,
  mobile_structure : Structure,
  transform? : Bool = true,
  final_optimization? : Bool = true,
) -> CeAlignmentResult raise CeAlignError {
  if self.reference_atoms.length() == 0 {
    raise CeAlignError("reference structure has not been set")
  }
  let mobile_atoms = cealign_get_guide_atoms(mobile_structure)
  if mobile_atoms.length() < self.window_size * 2 {
    raise CeAlignError(
      "Too few guide atoms in mobile structure: " +
      mobile_atoms.length().to_string() +
      "; need at least " +
      (self.window_size * 2).to_string(),
    )
  }
  let reference_coordinates : Array[Vector3] = []
  for atom in self.reference_atoms {
    reference_coordinates.push(atom.coordinate)
  }
  let mobile_coordinates : Array[Vector3] = []
  for atom in mobile_atoms {
    mobile_coordinates.push(atom.coordinate)
  }
  let paths = cealign_find_paths(
    reference_coordinates,
    mobile_coordinates,
    self.window_size,
    self.max_gap,
  )
  if paths.length() == 0 {
    raise CeAlignError("failed to find a suitable CE alignment path")
  }
  let longest_fragment_count = paths[0].fragment_count
  let mut best_path = paths[0]
  let mut best_fit = cealign_fit_path(
    best_path, reference_coordinates, mobile_coordinates,
  )
  let mut candidate_index = 1
  while candidate_index < paths.length() &&
        paths[candidate_index].fragment_count == longest_fragment_count {
    let candidate_fit = cealign_fit_path(
      paths[candidate_index],
      reference_coordinates,
      mobile_coordinates,
    )
    if candidate_fit.rmsd < best_fit.rmsd {
      best_path = paths[candidate_index]
      best_fit = candidate_fit
    }
    candidate_index = candidate_index + 1
  }
  let optimization_applied = final_optimization && best_path.z_score >= 3.5
  if optimization_applied {
    let optimized = cealign_optimize_path(
      best_path,
      reference_coordinates,
      mobile_coordinates,
      self.window_size,
      best_fit,
    )
    best_path = optimized.0
    best_fit = optimized.1
  }
  let rotation = best_fit.get_rot()
  let translation = best_fit.get_tran()
  let result_structure = if transform {
    cealign_transform_structure(mobile_structure, rotation, translation)
  } else {
    mobile_structure
  }
  CeAlignmentResult::{
    path: best_path,
    rmsd: best_fit.rmsd,
    rotation,
    translation,
    reference_guide_count: self.reference_atoms.length(),
    mobile_guide_count: mobile_atoms.length(),
    transformed: transform,
    optimization_applied,
    structure: result_structure,
  }
}

///|
fn cealign_sample_structure(
  id : String,
  angle : Double,
  translation : Vector3,
) -> Structure {
  let residues : Array[Residue] = []
  let cosine = @math.cos(angle)
  let sine = @math.sin(angle)
  let mut index = 0
  while index < 32 {
    let parameter = index.to_double() * 0.72
    let reference = Vector3::new(
      @math.cos(parameter) * 4.0,
      @math.sin(parameter) * 4.0,
      index.to_double() * 1.35,
    )
    let coordinate = Vector3::new(
      cosine * reference.x - sine * reference.y + translation.x,
      sine * reference.x + cosine * reference.y + translation.y,
      reference.z + translation.z,
    )
    let ca = Atom::new(
      name="CA",
      coord=coordinate,
      resname="ALA",
      chainid='A',
      resseq=index + 1,
      element="C",
    )
    let side = Atom::new(
      name="CB",
      coord=Vector3::new(
        coordinate.x + cosine,
        coordinate.y + sine,
        coordinate.z + 0.3,
      ),
      resname="ALA",
      chainid='A',
      resseq=index + 1,
      element="C",
    )
    residues.push(
      Residue::new(resname="ALA", chainid='A', resseq=index + 1, atoms=[
        ca, side,
      ]),
    )
    index = index + 1
  }
  Structure::new(id~, models=[
    Model::new(id=0, chains=[Chain::new(id='A', residues~)]),
  ])
}

///|
/// Deterministic sample structures for tests and examples.
pub fn cealign_sample_structures() -> (Structure, Structure) {
  (
    cealign_sample_structure("reference", 0.0, Vector3::new(0.0, 0.0, 0.0)),
    cealign_sample_structure("mobile", 0.65, Vector3::new(8.0, -3.0, 4.0)),
  )
}
