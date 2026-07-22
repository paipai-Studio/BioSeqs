///|
/// StructureAlignment module - Multiple protein structure alignment.
/// Inspired by Biopython Bio.PDB.StructureAlignment.
///
/// Provides functions for aligning multiple protein structures using
/// dynamic programming and structural similarity measures.

///|
/// Represents 3D coordinates for StructureAlignment.
pub struct SAPoint3D {
  x : Double
  y : Double
  z : Double
}

///|
/// Create a new SAPoint3D.
pub fn SAPoint3D::new(x : Double, y : Double, z : Double) -> SAPoint3D {
  SAPoint3D::{ x, y, z }
}

///|
/// Calculate distance between two points.
pub fn SAPoint3D::distance(self : SAPoint3D, other : SAPoint3D) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// A residue for StructureAlignment.
pub struct SAResidue {
  resname : String
  resseq : Int
  ca_coord : SAPoint3D
}

///|
/// Create a new SAResidue.
pub fn SAResidue::new(resname : String, resseq : Int, ca_coord : SAPoint3D) -> SAResidue {
  SAResidue::{ resname, resseq, ca_coord }
}

///|
/// A protein structure.
pub struct SAStructure {
  id : String
  residues : Array[SAResidue]
}

///|
/// Create a new SAStructure.
pub fn SAStructure::new(id : String, residues : Array[SAResidue]) -> SAStructure {
  SAStructure::{ id, residues }
}

///|
/// Structure alignment result.
pub struct SAAlignmentResult {
  structure_ids : Array[String]
  aligned_residues : Array[Array[SAResidue?]]
  rmsd : Double
  tm_score : Double
  aligned_length : Int
}

///|
/// Create a new SAAlignmentResult.
pub fn SAAlignmentResult::new(
  structure_ids : Array[String],
  aligned_residues : Array[Array[SAResidue?]],
  rmsd : Double,
  tm_score : Double,
  aligned_length : Int
) -> SAAlignmentResult {
  SAAlignmentResult::{ structure_ids, aligned_residues, rmsd, tm_score, aligned_length }
}

///|
/// Calculate the RMSD between two aligned residue sets.
pub fn sa_calculate_rmsd(residues1 : Array[SAPoint3D], residues2 : Array[SAPoint3D]) -> Double {
  if residues1.length() == 0 || residues2.length() == 0 {
    return 0.0
  }
  
  let len = if residues1.length() < residues2.length() { residues1.length() } else { residues2.length() }
  let mut sum_sq = 0.0
  
  let mut i = 0
  while i < len {
    let dx = residues1[i].x - residues2[i].x
    let dy = residues1[i].y - residues2[i].y
    let dz = residues1[i].z - residues2[i].z
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    i = i + 1
  }
  
  (sum_sq / len.to_double()).sqrt()
}

///|
/// Calculate the TM-score (Template Modeling score).
pub fn sa_calculate_tm_score(residues1 : Array[SAPoint3D], residues2 : Array[SAPoint3D], l_target : Int) -> Double {
  if residues1.length() == 0 || residues2.length() == 0 {
    return 0.0
  }
  
  let len = if residues1.length() < residues2.length() { residues1.length() } else { residues2.length() }
  let d0 = if l_target >= 15 {
    1.24 * @math.pow(l_target.to_double() - 15.0, 1.0 / 3.0) - 1.8
  } else {
    1.24 * @math.pow(15.0 - l_target.to_double(), 1.0 / 3.0) * 0.5
  }
  
  let mut sum = 0.0
  
  let mut i = 0
  while i < len {
    let dist = residues1[i].distance(residues2[i])
    sum = sum + 1.0 / (1.0 + (dist / d0) * (dist / d0))
    i = i + 1
  }
  
  sum / len.to_double()
}

///|
/// Calculate the distance matrix between two structures.
pub fn sa_calculate_distance_matrix(struct1 : SAStructure, struct2 : SAStructure) -> Array[Array[Double]] {
  let len1 = struct1.residues.length()
  let len2 = struct2.residues.length()
  
  let matrix : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < len1 {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < len2 {
      let dist = struct1.residues[i].ca_coord.distance(struct2.residues[j].ca_coord)
      row.push(dist)
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  
  matrix
}

///|
/// Perform dynamic programming alignment between two structures.
pub fn sa_dp_alignment(struct1 : SAStructure, struct2 : SAStructure, gap_open? : Double = 5.0, gap_extend? : Double = 1.0) -> (Array[(Int, Int)], Double) {
  let dist_matrix = sa_calculate_distance_matrix(struct1, struct2)
  let len1 = struct1.residues.length()
  let len2 = struct2.residues.length()
  
  let score : Array[Array[Double]] = Array::new()
  let trace : Array[Array[Int]] = Array::new()
  
  let mut i = 0
  while i <= len1 {
    let score_row : Array[Double] = Array::new()
    let trace_row : Array[Int] = Array::new()
    let mut j = 0
    while j <= len2 {
      score_row.push(0.0)
      trace_row.push(0)
      j = j + 1
    }
    score.push(score_row)
    trace.push(trace_row)
    i = i + 1
  }
  
  i = 1
  while i <= len1 {
    score[i].set(0, score[i-1][0] - gap_open)
    trace[i].set(0, 1)
    i = i + 1
  }
  
  let mut j = 1
  while j <= len2 {
    score[0].set(j, score[0][j-1] - gap_open)
    trace[0].set(j, 2)
    j = j + 1
  }
  
  i = 1
  while i <= len1 {
    j = 1
    while j <= len2 {
      let gap_open1 = if trace[i-1][j] == 1 { gap_extend } else { gap_open }
      let gap_open2 = if trace[i][j-1] == 2 { gap_extend } else { gap_open }
      
      let match_score = score[i-1][j-1] - dist_matrix[i-1][j-1]
      let delete_score = score[i-1][j] - gap_open1
      let insert_score = score[i][j-1] - gap_open2
      
      if match_score >= delete_score && match_score >= insert_score {
        score[i].set(j, match_score)
        trace[i].set(j, 0)
      } else if delete_score >= insert_score {
        score[i].set(j, delete_score)
        trace[i].set(j, 1)
      } else {
        score[i].set(j, insert_score)
        trace[i].set(j, 2)
      }
      j = j + 1
    }
    i = i + 1
  }
  
  let alignment : Array[(Int, Int)] = Array::new()
  i = len1
  j = len2
  
  while i > 0 && j > 0 {
    let t = trace[i][j]
    if t == 0 {
      alignment.push((i-1, j-1))
      i = i - 1
      j = j - 1
    } else if t == 1 {
      i = i - 1
    } else {
      j = j - 1
    }
  }
  
  let reversed : Array[(Int, Int)] = Array::new()
  let align_len = alignment.length()
  let mut k = 0
  while k < align_len {
    reversed.push(alignment[align_len - 1 - k])
    k = k + 1
  }
  
  (reversed, score[len1][len2])
}

///|
/// Align two protein structures.
pub fn sa_align_two_structures(struct1 : SAStructure, struct2 : SAStructure) -> SAAlignmentResult {
  let (alignment, _) = sa_dp_alignment(struct1, struct2)
  
  let coords1 : Array[SAPoint3D] = Array::new()
  let coords2 : Array[SAPoint3D] = Array::new()
  let aligned_residues : Array[Array[SAResidue?]] = Array::new()
  
  let align_len = alignment.length()
  let mut idx = 0
  while idx < align_len {
    let pair = alignment[idx]
    let i = pair.0
    let j = pair.1
    
    coords1.push(struct1.residues[i].ca_coord)
    coords2.push(struct2.residues[j].ca_coord)
    
    let row : Array[SAResidue?] = Array::new()
    row.push(Some(struct1.residues[i]))
    row.push(Some(struct2.residues[j]))
    aligned_residues.push(row)
    
    idx = idx + 1
  }
  
  let rmsd = sa_calculate_rmsd(coords1, coords2)
  let tm_score = sa_calculate_tm_score(coords1, coords2, struct1.residues.length())
  
  SAAlignmentResult::new(
    [struct1.id, struct2.id],
    aligned_residues,
    rmsd,
    tm_score,
    alignment.length()
  )
}

///|
/// Perform multiple structure alignment using progressive approach.
pub fn sa_align_multiple_structures(structures : Array[SAStructure]) -> SAAlignmentResult {
  if structures.length() == 0 {
    return SAAlignmentResult::new([], [], 0.0, 0.0, 0)
  }
  
  if structures.length() == 1 {
    let aligned_residues : Array[Array[SAResidue?]] = Array::new()
    let mut i = 0
    while i < structures[0].residues.length() {
      let row : Array[SAResidue?] = Array::new()
      row.push(Some(structures[0].residues[i]))
      aligned_residues.push(row)
      i = i + 1
    }
    return SAAlignmentResult::new([structures[0].id], aligned_residues, 0.0, 1.0, structures[0].residues.length())
  }
  
  let mut current_alignment = sa_align_two_structures(structures[0], structures[1])
  
  let mut k = 2
  while k < structures.length() {
    current_alignment = sa_align_structure_to_profile(current_alignment, structures[k])
    k = k + 1
  }
  
  current_alignment
}

///|
/// Align a structure to an existing alignment profile.
pub fn sa_align_structure_to_profile(alignment : SAAlignmentResult, structure : SAStructure) -> SAAlignmentResult {
  let num_structures = alignment.structure_ids.length()
  let aligned_length = alignment.aligned_length
  
  let consensus_coords : Array[SAPoint3D] = Array::new()
  
  let mut i = 0
  while i < aligned_length {
    let mut sum_x = 0.0
    let mut sum_y = 0.0
    let mut sum_z = 0.0
    let mut count = 0
    
    let mut j = 0
    while j < num_structures {
      match alignment.aligned_residues[i][j] {
        Some(res) => {
          sum_x = sum_x + res.ca_coord.x
          sum_y = sum_y + res.ca_coord.y
          sum_z = sum_z + res.ca_coord.z
          count = count + 1
        }
        None => ignore({ () })
      }
      j = j + 1
    }
    
    if count > 0 {
      consensus_coords.push(SAPoint3D::new(sum_x / count.to_double(), sum_y / count.to_double(), sum_z / count.to_double()))
    }
    i = i + 1
  }
  
  let consensus_struct = SAStructure::new("consensus", sa_residues_from_coords(consensus_coords))
  
  let (new_alignment, _) = sa_dp_alignment(consensus_struct, structure)
  
  let new_aligned_residues : Array[Array[SAResidue?]] = Array::new()
  
  let new_align_len = new_alignment.length()
  let mut idx = 0
  while idx < new_align_len {
    let pair = new_alignment[idx]
    let i = pair.0
    let j = pair.1
    
    let row : Array[SAResidue?] = Array::new()
    
    let mut k = 0
    while k < num_structures {
      let val = if i < alignment.aligned_length { alignment.aligned_residues[i][k] } else { None }
      row.push(val)
      k = k + 1
    }
    
    row.push(Some(structure.residues[j]))
    new_aligned_residues.push(row)
    
    idx = idx + 1
  }
  
  let new_ids : Array[String] = Array::new()
  let mut m = 0
  while m < alignment.structure_ids.length() {
    new_ids.push(alignment.structure_ids[m])
    m = m + 1
  }
  new_ids.push(structure.id)
  
  let coords1 : Array[SAPoint3D] = Array::new()
  let coords2 : Array[SAPoint3D] = Array::new()
  
  idx = 0
  while idx < new_align_len {
    let pair = new_alignment[idx]
    let i = pair.0
    let j = pair.1
    
    if i < consensus_coords.length() {
      coords1.push(consensus_coords[i])
    }
    coords2.push(structure.residues[j].ca_coord)
    
    idx = idx + 1
  }
  
  let rmsd = sa_calculate_rmsd(coords1, coords2)
  let tm_score = sa_calculate_tm_score(coords1, coords2, consensus_struct.residues.length())
  
  SAAlignmentResult::new(new_ids, new_aligned_residues, rmsd, tm_score, new_alignment.length())
}

///|
/// Create residues from coordinates.
pub fn sa_residues_from_coords(coords : Array[SAPoint3D]) -> Array[SAResidue] {
  let residues : Array[SAResidue] = Array::new()
  
  let mut i = 0
  while i < coords.length() {
    residues.push(SAResidue::new("XXX", i + 1, coords[i]))
    i = i + 1
  }
  
  residues
}

///|
/// Calculate pairwise RMSD matrix for multiple structures.
pub fn sa_calculate_pairwise_rmsd(structures : Array[SAStructure]) -> Array[Array[Double]] {
  let n = structures.length()
  
  let matrix : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      if i == j {
        row.push(0.0)
      } else {
        let alignment = sa_align_two_structures(structures[i], structures[j])
        row.push(alignment.rmsd)
      }
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  
  matrix
}

///|
/// Bio API functions.

pub fn bio_structure_align_two(struct1 : SAStructure, struct2 : SAStructure) -> SAAlignmentResult {
  sa_align_two_structures(struct1, struct2)
}

pub fn bio_structure_align_multiple(structures : Array[SAStructure]) -> SAAlignmentResult {
  sa_align_multiple_structures(structures)
}

pub fn bio_structure_calculate_rmsd(struct1 : SAStructure, struct2 : SAStructure) -> Double {
  let alignment = sa_align_two_structures(struct1, struct2)
  alignment.rmsd
}

pub fn bio_structure_calculate_tm_score(struct1 : SAStructure, struct2 : SAStructure) -> Double {
  let alignment = sa_align_two_structures(struct1, struct2)
  alignment.tm_score
}

pub fn bio_structure_pairwise_rmsd(structures : Array[SAStructure]) -> Array[Array[Double]] {
  sa_calculate_pairwise_rmsd(structures)
}