///|
/// Disulfide bond and metal coordination site detection.
///
/// Biotite provides `disulfide_bonds()` and metal coordination analysis.
/// This module implements these for BioSeqs.

///|
/// A disulfide bond between two cysteine residues.
pub struct SSBond {
  residue1 : Int32 // cysteine residue index
  residue2 : Int32 // cysteine residue index
  sg_distance : Double // S-S distance in Angstroms
  chi1_1 : Double // chi1 angle of first cysteine
  chi1_2 : Double // chi1 angle of second cysteine
}

///|
/// Get the S-S distance.
pub fn SSBond::get_distance(self : SSBond) -> Double {
  self.sg_distance
}

///|
/// Detect disulfide bonds in a structure.
///
/// Two cysteine residues form a disulfide bond if their SG atoms are
/// within 2.5 Å of each other (typical S-S bond length is ~2.05 Å).
pub fn detect_disulfide_bonds(
  cys_residues : Array[(Int32, Vector3)],
  max_distance~ : Double = 2.5,
) -> Array[SSBond] {
  let bonds : Array[SSBond] = []
  let n = cys_residues.length()
  for i in 0..<n {
    for j in (i + 1)..<n {
      let (res1, sg1) = cys_residues[i]
      let (res2, sg2) = cys_residues[j]
      let dist = sg1.distance(sg2)
      if dist <= max_distance {
        bonds.push(SSBond::{
          residue1: res1,
          residue2: res2,
          sg_distance: dist,
          chi1_1: 0.0, // would require full backbone calculation
          chi1_2: 0.0,
        })
      }
    }
  }
  bonds
}

///|
/// Check if a disulfide bond is in the typical range.
///
/// Typical S-S bond length: 2.02–2.08 Å
pub fn is_typical_disulfide(bond : SSBond) -> Bool {
  bond.sg_distance >= 1.9 && bond.sg_distance <= 2.2
}

///|
/// Find all cysteine residues with SG atom coordinates.
pub fn find_cysteine_sg(residues : Array[Residue]) -> Array[(Int32, Vector3)] {
  let result : Array[(Int32, Vector3)] = []
  for i in 0..<residues.length() {
    let r = residues[i]
    if r.resname == "CYS" {
      for a in r.get_atoms() {
        if a.name == "SG" {
          result.push((i, a.get_coord()))
          break
        }
      }
    }
  }
  result
}

///|
/// Metal coordination geometry classification.
pub enum CoordinationGeometry {
  Octahedral // 6 ligands
  Tetrahedral // 4 ligands
  SquarePlanar // 4 ligands, planar
  TrigonalBipyramidal // 5 ligands
  SquarePyramidal // 5 ligands
  Linear // 2 ligands
  Unknown
}

///|
/// Get geometry as string.
pub fn CoordinationGeometry::to_string(self : CoordinationGeometry) -> String {
  match self {
    Octahedral => "octahedral"
    Tetrahedral => "tetrahedral"
    SquarePlanar => "square planar"
    TrigonalBipyramidal => "trigonal bipyramidal"
    SquarePyramidal => "square pyramidal"
    Linear => "linear"
    Unknown => "unknown"
  }
}

///|
/// A metal coordination site.
pub struct MetalCoordSite {
  metal_index : Int32 // atom index of metal
  metal_type : String // element symbol
  ligand_indices : Array[Int32] // coordinating atom indices
  ligand_distances : Array[Double] // M-L distances
  geometry : CoordinationGeometry
}

///|
/// Get the number of ligands.
pub fn MetalCoordSite::coordination_number(self : MetalCoordSite) -> Int32 {
  self.ligand_indices.length()
}

///|
/// Detect metal coordination sites.
///
/// Finds metal ions and their coordinating atoms (typically N, O, S
/// within 3.0 Å of the metal).
pub fn detect_metal_sites(
  coords : Array[Vector3],
  elements : Array[String],
  _atom_names : Array[String],
  max_distance~ : Double = 3.0,
) -> Array[MetalCoordSite] {
  let metals = ["FE", "ZN", "CU", "MG", "CA", "MN", "NI", "CO", "MO", "W"]
  let sites : Array[MetalCoordSite] = []
  for i in 0..<coords.length() {
    let mut is_metal = false
    for m in metals {
      if elements[i] == m {
        is_metal = true
        break
      }
    }
    if !is_metal {
      continue
    }
    let ligand_indices : Array[Int32] = []
    let ligand_distances : Array[Double] = []
    for j in 0..<coords.length() {
      if i == j {
        continue
      }
      // Only consider coordinating atoms (N, O, S)
      let is_donor = elements[j] == "N" ||
        elements[j] == "O" ||
        elements[j] == "S"
      if !is_donor {
        continue
      }
      let dist = coords[i].distance(coords[j])
      if dist <= max_distance {
        ligand_indices.push(j)
        ligand_distances.push(dist)
      }
    }
    if ligand_indices.length() >= 2 {
      let geom = classify_coordination_geometry(ligand_indices.length())
      sites.push(MetalCoordSite::{
        metal_index: i,
        metal_type: elements[i],
        ligand_indices,
        ligand_distances,
        geometry: geom,
      })
    }
  }
  sites
}

///|
/// Classify coordination geometry based on ligand count.
fn classify_coordination_geometry(n_ligands : Int32) -> CoordinationGeometry {
  match n_ligands {
    2 => Linear
    4 => Tetrahedral // could also be square planar
    5 => TrigonalBipyramidal // could also be square pyramidal
    6 => Octahedral
    _ => Unknown
  }
}

///|
/// Calculate average metal-ligand distance.
pub fn average_metal_ligand_distance(site : MetalCoordSite) -> Double {
  if site.ligand_distances.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for d in site.ligand_distances {
    sum = sum + d
  }
  sum / site.ligand_distances.length().to_double()
}

///|
fn _coord_geom_unused_constructors() -> Unit {
  ignore(CoordinationGeometry::SquarePlanar)
  ignore(CoordinationGeometry::SquarePyramidal)
}
