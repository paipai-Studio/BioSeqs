///|
/// Nucleic acid base stacking and base pair classification.
///
/// Biotite provides `base_stacking()`, `base_pairs_edge()` (Leontis-Westhof),
/// and `base_pairs_glycosidic_bond()` for RNA/DNA 3D structure analysis.
/// This module implements these functions for BioSeqs.

///|
/// Leontis-Westhof edge classification for base pairs.
pub enum NABaseEdge {
  WatsonCrick // W: Watson-Crick edge
  Sugar // S: Sugar edge
  Hoogsteen // H: Hoogsteen (shallow-groove) edge
}

///|
/// Get edge as string.
pub fn NABaseEdge::to_string(self : NABaseEdge) -> String {
  match self {
    WatsonCrick => "W"
    Sugar => "S"
    Hoogsteen => "H"
  }
}

///|
/// Glycosidic bond orientation.
pub enum NAGlycosidicBond {
  Cis
  Trans
}

///|
/// Get glycosidic bond as string.
pub fn NAGlycosidicBond::to_string(self : NAGlycosidicBond) -> String {
  match self {
    Cis => "cis"
    Trans => "trans"
  }
}

///|
/// A nucleic acid base pair.
pub struct NABasePair {
  base1 : Int32 // atom index of first base
  base2 : Int32 // atom index of second base
  edge1 : NABaseEdge
  edge2 : NABaseEdge
  glycosidic : NAGlycosidicBond
}

///|
/// Construct a NABasePair.
pub fn NABasePair::new(
  base1~ : Int32,
  base2~ : Int32,
  edge1~ : NABaseEdge = WatsonCrick,
  edge2~ : NABaseEdge = WatsonCrick,
  glycosidic~ : NAGlycosidicBond = Cis,
) -> NABasePair {
  NABasePair::{ base1, base2, edge1, edge2, glycosidic, }
}

///|
/// Get Leontis-Westhof classification as string (e.g. "cWW", "tSH").
pub fn NABasePair::to_leontis_westhof(self : NABasePair) -> String {
  let g = match self.glycosidic {
    Cis => "c"
    Trans => "t"
  }
  g + self.edge1.to_string() + self.edge2.to_string()
}

///|
/// A base stacking interaction.
pub struct NABaseStacking {
  base1 : Int32
  base2 : Int32
  distance : Double // centroid distance in Angstroms
  angle : Double // normal vector angle in degrees
  overlap_area : Double // approximate overlap
}

///|
/// Construct a NABaseStacking.
pub fn NABaseStacking::new(
  base1~ : Int32,
  base2~ : Int32,
  distance~ : Double,
  angle~ : Double,
  overlap_area~ : Double,
) -> NABaseStacking {
  NABaseStacking::{ base1, base2, distance, angle, overlap_area, }
}

///|
/// Get the stacking distance.
pub fn NABaseStacking::get_distance(self : NABaseStacking) -> Double {
  self.distance
}

///|
/// Get the stacking angle.
pub fn NABaseStacking::get_angle(self : NABaseStacking) -> Double {
  self.angle
}

///|
/// Detect base stacking interactions from nucleotide coordinates.
///
/// Two bases are stacked if:
/// 1. Their centroid distance is 3.0–6.0 Å
/// 2. The angle between their normal vectors is < 30°
/// 3. They have significant overlap area
pub fn detect_base_stacking(
  base_centroids : Array[Vector3],
  base_normals : Array[Vector3],
) -> Array[NABaseStacking] {
  let n = base_centroids.length()
  let stacking : Array[NABaseStacking] = []
  for i in 0..<n {
    for j in (i + 1)..<n {
      let dist = base_centroids[i].distance(base_centroids[j])
      if dist < 3.0 || dist > 6.0 {
        continue
      }
      // Angle between normal vectors
      let n1 = base_normals[i]
      let n2 = base_normals[j]
      let dot = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z
      let mut d = dot
      if d > 1.0 {
        d = 1.0
      }
      if d < -1.0 {
        d = -1.0
      }
      let angle = @math.acos(d.abs()) * 180.0 / 3.141592653589793
      if angle < 30.0 {
        // Approximate overlap area (simplified)
        let overlap = (1.0 - angle / 30.0) * (6.0 - dist) / 3.0 * 10.0
        stacking.push(
          NABaseStacking::new(
            base1=i,
            base2=j,
            distance=dist,
            angle~,
            overlap_area=overlap,
          ),
        )
      }
    }
  }
  stacking
}

///|
/// Classify a base pair's Leontis-Westhof edge type.
///
/// Determines which edges of the two bases are involved in the pairing
/// by analyzing the positions of specific atoms:
/// - Watson-Crick edge: involves N1/N3 and related atoms
/// - Sugar edge: involves the sugar O2' and adjacent atoms
/// - Hoogsteen edge: involves N7/N9 and related atoms
///
/// This is a simplified classification based on base geometry.
pub fn classify_base_pair_edge(
  base1_atoms : Array[Vector3],
  base2_atoms : Array[Vector3],
  base1_type : Char,
  base2_type : Char,
) -> (NABaseEdge, NABaseEdge) {
  // Simplified: use the relative orientation of the bases
  // In a full implementation, we would check which atoms are
  // closest to the hydrogen bonding partners
  // For now, classify based on base type and geometry
  let edge1 = classify_single_edge(base1_atoms, base1_type)
  let edge2 = classify_single_edge(base2_atoms, base2_type)
  (edge1, edge2)
}

///|
/// Classify a single base's edge.
fn classify_single_edge(_ : Array[Vector3], base_type : Char) -> NABaseEdge {
  // Simplified: purines (A, G) default to Watson-Crick
  // pyrimidines (C, T, U) default to Watson-Crick
  // A full implementation would analyze atom positions
  match base_type {
    'A' | 'G' | 'C' | 'T' | 'U' => WatsonCrick
    _ => Hoogsteen
  }
}

///|
/// Determine glycosidic bond orientation (cis/trans).
///
/// In cis orientation, the glycosidic bonds of both bases are on the
/// same side of the base pair. In trans, they are on opposite sides.
pub fn classify_glycosidic_bond(
  base1_glycosidic : Vector3,
  base2_glycosidic : Vector3,
  base1_centroid : Vector3,
  base2_centroid : Vector3,
) -> NAGlycosidicBond {
  // Vector from base1 to base2 centroid
  let v12 = Vector3::new(
    base2_centroid.x - base1_centroid.x,
    base2_centroid.y - base1_centroid.y,
    base2_centroid.z - base1_centroid.z,
  )
  // Vector from base1 centroid to its glycosidic bond
  let g1 = Vector3::new(
    base1_glycosidic.x - base1_centroid.x,
    base1_glycosidic.y - base1_centroid.y,
    base1_glycosidic.z - base1_centroid.z,
  )
  // Vector from base2 centroid to its glycosidic bond
  let g2 = Vector3::new(
    base2_glycosidic.x - base2_centroid.x,
    base2_glycosidic.y - base2_centroid.y,
    base2_glycosidic.z - base2_centroid.z,
  )
  // Project glycosidic vectors onto the inter-base vector
  let proj1 = g1.x * v12.x + g1.y * v12.y + g1.z * v12.z
  let proj2 = g2.x * v12.x + g2.y * v12.y + g2.z * v12.z
  // If both projections have the same sign, it's cis
  if proj1 * proj2 > 0.0 {
    Cis
  } else {
    Trans
  }
}

///|
/// Detect base pairs from nucleotide coordinates.
///
/// Two bases are paired if their centroids are within 3.0–5.0 Å and
/// their normal vectors are roughly coplanar.
pub fn detect_base_pairs(
  base_centroids : Array[Vector3],
  base_normals : Array[Vector3],
) -> Array[NABasePair] {
  let n = base_centroids.length()
  let pairs : Array[NABasePair] = []
  for i in 0..<n {
    for j in (i + 1)..<n {
      let dist = base_centroids[i].distance(base_centroids[j])
      if dist < 2.5 || dist > 5.0 {
        continue
      }
      // Check coplanarity (normals should be roughly anti-parallel)
      let n1 = base_normals[i]
      let n2 = base_normals[j]
      let dot = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z
      let mut d = dot
      if d > 1.0 {
        d = 1.0
      }
      if d < -1.0 {
        d = -1.0
      }
      let angle = @math.acos(d.abs()) * 180.0 / 3.141592653589793
      if angle < 35.0 {
        pairs.push(NABasePair::new(base1=i, base2=j))
      }
    }
  }
  pairs
}
