///|
/// ProDy-style PDBEnsemble with atom mapping.
///
/// Handles multi-conformation ensembles where conformations may have
/// different atom compositions (missing atoms, different ordering).
/// Atom mapping matches atoms by (chainid, resseq, atom_name).
///
/// Features:
///   - Build from reference + multiple mobile structures
///   - Atom mapping with missing-atom handling
///   - Iterative superimposition to reference
///   - Pairwise RMSD matrix, per-atom RMSF
///   - Mean structure, coordinate deviation
///   - PCA on ensemble fluctuations

///|
/// Atom mapping between a reference and a mobile structure.
/// ref_indices[i] and mob_indices[i] are paired atoms.
pub struct AtomMap {
  ref_indices : Array[Int]
  mob_indices : Array[Int]
  n_mapped : Int
}

///|
/// Build an atom map by matching (chainid, resseq, atom_name).
/// Returns pairs of (ref_idx, mob_idx) for matched atoms.
pub fn AtomMap::new(
  ref_atoms : Array[Atom],
  mob_atoms : Array[Atom],
) -> AtomMap {
  let ref_indices : Array[Int] = Array(capacity=0)
  let mob_indices : Array[Int] = Array(capacity=0)
  // For each reference atom, find matching mobile atom
  for i = 0; i < ref_atoms.length(); i = i + 1 {
    let ra = ref_atoms[i]
    for j = 0; j < mob_atoms.length(); j = j + 1 {
      let ma = mob_atoms[j]
      if ra.chainid == ma.chainid &&
         ra.resseq == ma.resseq &&
         ra.name == ma.name {
        ref_indices.push(i)
        mob_indices.push(j)
        break
      }
    }
  }
  AtomMap::{ ref_indices, mob_indices, n_mapped: ref_indices.length() }
}

///|
/// Number of mapped atom pairs.
pub fn AtomMap::num_mapped(self : AtomMap) -> Int {
  self.n_mapped
}

///|
/// Get the reference indices.
pub fn AtomMap::ref_indices(self : AtomMap) -> Array[Int] {
  self.ref_indices
}

///|
/// Get the mobile indices.
pub fn AtomMap::mob_indices(self : AtomMap) -> Array[Int] {
  self.mob_indices
}

///|
/// Coverage: fraction of reference atoms mapped.
pub fn AtomMap::coverage(self : AtomMap, n_ref : Int) -> Double {
  if n_ref == 0 {
    return 0.0
  }
  self.n_mapped.to_double() / n_ref.to_double()
}

// ===== PDBEnsemble =====

///|
/// A multi-conformation ensemble with atom mapping.
/// Stores a reference structure and multiple coordinate sets,
/// each mapped to the reference atom ordering.
pub struct PDBEnsemble {
  n_atoms : Int
  // Reference atom metadata for mapping
  atom_names : Array[String]
  atom_chains : Array[Char]
  atom_resseqs : Array[Int]
  // Reference coordinates
  ref_coords : Array[Vector3]
  // All conformation coordinates (mapped to reference ordering)
  // coordsets[conf][i] = coordinate of reference atom i in conformation conf
  // Missing atoms stored as (NaN, NaN, NaN)
  coordsets : Array[Array[Vector3]]
  // Superimposed coordinates
  aligned : Array[Array[Vector3]]
  is_aligned : Bool
  // Weights for superimposition (default: all 1.0)
  weights : Array[Double]
}

///|
/// Build a PDBEnsemble from a reference chain.
/// Uses all atoms (or CA atoms if specified).
pub fn PDBEnsemble::new(ref_chain : Chain, calpha_only : Bool) -> PDBEnsemble? {
  let ref_atoms = if calpha_only {
    get_ca_atoms(ref_chain)
  } else {
    ref_chain.get_atoms()
  }
  let n = ref_atoms.length()
  if n < 3 {
    return None
  }
  let names : Array[String] = Array(capacity=n)
  let chains : Array[Char] = Array(capacity=n)
  let resseqs : Array[Int] = Array(capacity=n)
  let coords : Array[Vector3] = Array(capacity=n)
  let wts : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    names.push(ref_atoms[i].name)
    chains.push(ref_atoms[i].chainid)
    resseqs.push(ref_atoms[i].resseq)
    coords.push(ref_atoms[i].coord)
    wts.push(1.0)
  }
  Some(PDBEnsemble::{
    n_atoms: n,
    atom_names: names,
    atom_chains: chains,
    atom_resseqs: resseqs,
    ref_coords: coords,
    coordsets: [],
    aligned: [],
    is_aligned: false,
    weights: wts,
  })
}

///|
/// Extract CA atoms from a chain.
fn get_ca_atoms(chain : Chain) -> Array[Atom] {
  let result : Array[Atom] = Array(capacity=0)
  for res in chain.get_residues() {
    match res.get_ca() {
      Some(atom) => result.push(atom)
      None => ()
    }
  }
  result
}

///|
/// Add a conformation from a chain.
/// Atoms are mapped to the reference by (chainid, resseq, name).
/// Missing atoms are filled with NaN.
pub fn PDBEnsemble::add_chain(self : PDBEnsemble, chain : Chain) -> PDBEnsemble {
  let mob_atoms = chain.get_atoms()
  let coords : Array[Vector3] = Array::make(self.n_atoms, Vector3::new(999999.999, 999999.999, 999999.999))
  // Build a quick lookup: (chainid, resseq, name) → coord
  for j = 0; j < mob_atoms.length(); j = j + 1 {
    let ma = mob_atoms[j]
    for i = 0; i < self.n_atoms; i = i + 1 {
      if self.atom_chains[i] == ma.chainid &&
         self.atom_resseqs[i] == ma.resseq &&
         self.atom_names[i] == ma.name {
        coords[i] = ma.coord
        break
      }
    }
  }
  let new_coordsets = Array(capacity=0)
  for c in self.coordsets {
    new_coordsets.push(c)
  }
  new_coordsets.push(coords)
  PDBEnsemble::{ ..self, coordsets: new_coordsets, is_aligned: false }
}

///|
/// Add a conformation from raw coordinates (already in reference ordering).
pub fn PDBEnsemble::add_coordset(self : PDBEnsemble, coords : Array[Vector3]) -> PDBEnsemble {
  let n = if coords.length() < self.n_atoms {
    coords.length()
  } else {
    self.n_atoms
  }
  let full : Array[Vector3] = Array::make(self.n_atoms, Vector3::new(999999.999, 999999.999, 999999.999))
  for i = 0; i < n; i = i + 1 {
    full[i] = coords[i]
  }
  let new_coordsets = Array(capacity=0)
  for c in self.coordsets {
    new_coordsets.push(c)
  }
  new_coordsets.push(full)
  PDBEnsemble::{ ..self, coordsets: new_coordsets, is_aligned: false }
}

///|
/// Add a conformation using an explicit AtomMap.
pub fn PDBEnsemble::add_mapped(
  self : PDBEnsemble,
  mob_atoms : Array[Atom],
  map : AtomMap,
) -> PDBEnsemble {
  let coords : Array[Vector3] = Array::make(self.n_atoms, Vector3::new(999999.999, 999999.999, 999999.999))
  for k = 0; k < map.n_mapped; k = k + 1 {
    let ri = map.ref_indices[k]
    let mi = map.mob_indices[k]
    if ri < self.n_atoms && mi < mob_atoms.length() {
      coords[ri] = mob_atoms[mi].coord
    }
  }
  let new_coordsets = Array(capacity=0)
  for c in self.coordsets {
    new_coordsets.push(c)
  }
  new_coordsets.push(coords)
  PDBEnsemble::{ ..self, coordsets: new_coordsets, is_aligned: false }
}

// ===== Superimposition =====

///|
/// Superimpose all conformations onto the reference.
/// Uses weighted Kabsch algorithm.
pub fn PDBEnsemble::superimpose(self : PDBEnsemble) -> PDBEnsemble {
  if self.coordsets.length() == 0 {
    return PDBEnsemble::{ ..self, is_aligned: true }
  }
  let aligned : Array[Array[Vector3]] = Array(capacity=0)
  for conf = 0; conf < self.coordsets.length(); conf = conf + 1 {
    let (mapped_ref, mapped_mob) = get_mapped_coords(self, conf)
    if mapped_ref.length() < 3 {
      // Not enough atoms, just copy
      aligned.push(copy_coords(self.coordsets[conf]))
    } else {
      let (rot, trans) = vector_superposition_transform(mapped_mob, mapped_ref)
      let transformed : Array[Vector3] = Array::make(self.n_atoms, Vector3::new(999999.999, 999999.999, 999999.999))
      for i = 0; i < self.n_atoms; i = i + 1 {
        let c = self.coordsets[conf][i]
        if is_valid_coord(c) {
          transformed[i] = rot.transform(c).add(trans)
        }
      }
      aligned.push(transformed)
    }
  }
  PDBEnsemble::{ ..self, aligned, is_aligned: true }
}

///|
/// Get mapped (non-NaN) coordinates for a conformation.
fn get_mapped_coords(ens : PDBEnsemble, conf : Int) -> (Array[Vector3], Array[Vector3]) {
  let ref : Array[Vector3] = Array(capacity=0)
  let mob : Array[Vector3] = Array(capacity=0)
  for i = 0; i < ens.n_atoms; i = i + 1 {
    let c = ens.coordsets[conf][i]
    if is_valid_coord(c) {
      ref.push(ens.ref_coords[i])
      mob.push(c)
    }
  }
  (ref, mob)
}

///|
/// Check if a coordinate is valid (not NaN/placeholder).
fn is_valid_coord(v : Vector3) -> Bool {
  v.x < 100000.0 && v.y < 100000.0 && v.z < 100000.0
}

///|
fn copy_coords(coords : Array[Vector3]) -> Array[Vector3] {
  let result : Array[Vector3] = Array(capacity=0)
  for c in coords {
    result.push(c)
  }
  result
}

// ===== Analysis =====

///|
/// Get the aligned coordinates for a conformation.
/// If not yet superimposed, returns original coordinates.
pub fn PDBEnsemble::get_coords(self : PDBEnsemble, conf : Int) -> Array[Vector3] {
  if self.is_aligned && conf < self.aligned.length() {
    self.aligned[conf]
  } else if conf < self.coordsets.length() {
    self.coordsets[conf]
  } else {
    Array::make(self.n_atoms, Vector3::zero())
  }
}

///|
/// Get reference coordinates.
pub fn PDBEnsemble::get_ref_coords(self : PDBEnsemble) -> Array[Vector3] {
  self.ref_coords
}

///|
/// Number of atoms in the reference.
pub fn PDBEnsemble::num_atoms(self : PDBEnsemble) -> Int {
  self.n_atoms
}

///|
/// Number of conformations.
pub fn PDBEnsemble::num_conformations(self : PDBEnsemble) -> Int {
  self.coordsets.length()
}

///|
/// Check if aligned.
pub fn PDBEnsemble::is_aligned(self : PDBEnsemble) -> Bool {
  self.is_aligned
}

///|
/// RMSD of each conformation from the reference.
/// If aligned, uses aligned coords; otherwise uses raw coords.
pub fn PDBEnsemble::get_rmsd(self : PDBEnsemble) -> Array[Double] {
  let result : Array[Double] = Array(capacity=0)
  for conf = 0; conf < self.coordsets.length(); conf = conf + 1 {
    let coords = self.get_coords(conf)
    let mut sum_sq = 0.0
    let mut n_valid = 0
    for i = 0; i < self.n_atoms; i = i + 1 {
      if is_valid_coord(coords[i]) {
        let d = self.ref_coords[i].sub(coords[i])
        sum_sq = sum_sq + d.dot(d)
        n_valid = n_valid + 1
      }
    }
    if n_valid > 0 {
      result.push((sum_sq / n_valid.to_double()).sqrt())
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
/// Pairwise RMSD matrix between all conformations.
/// Returns n_conf × n_conf matrix.
pub fn PDBEnsemble::get_rmsd_matrix(self : PDBEnsemble) -> Array[Array[Double]] {
  let n_conf = self.coordsets.length()
  let matrix : Array[Array[Double]] = Array::makei(n_conf, _ => Array::make(n_conf, 0.0))
  for i = 0; i < n_conf; i = i + 1 {
    let ci = self.get_coords(i)
    for j = i + 1; j < n_conf; j = j + 1 {
      let cj = self.get_coords(j)
      let mut sum_sq = 0.0
      let mut n_valid = 0
      for k = 0; k < self.n_atoms; k = k + 1 {
        if is_valid_coord(ci[k]) && is_valid_coord(cj[k]) {
          let d = ci[k].sub(cj[k])
          sum_sq = sum_sq + d.dot(d)
          n_valid = n_valid + 1
        }
      }
      let rmsd = if n_valid > 0 {
        (sum_sq / n_valid.to_double()).sqrt()
      } else {
        0.0
      }
      matrix[i][j] = rmsd
      matrix[j][i] = rmsd
    }
  }
  matrix
}

///|
/// Per-atom RMSF (root mean square fluctuation).
/// Requires the ensemble to be aligned.
pub fn PDBEnsemble::get_rmsf(self : PDBEnsemble) -> Array[Double] {
  let n_conf = self.coordsets.length()
  let rmsf : Array[Double] = Array::make(self.n_atoms, 0.0)
  if n_conf == 0 {
    return rmsf
  }
  // Compute mean position
  let mean : Array[Vector3] = Array::make(self.n_atoms, Vector3::zero())
  let counts : Array[Int] = Array::make(self.n_atoms, 0)
  for conf = 0; conf < n_conf; conf = conf + 1 {
    let coords = self.get_coords(conf)
    for i = 0; i < self.n_atoms; i = i + 1 {
      if is_valid_coord(coords[i]) {
        mean[i] = mean[i].add(coords[i])
        counts[i] = counts[i] + 1
      }
    }
  }
  for i = 0; i < self.n_atoms; i = i + 1 {
    if counts[i] > 0 {
      mean[i] = mean[i].scale(1.0 / counts[i].to_double())
    }
  }
  // Compute RMSF
  for conf = 0; conf < n_conf; conf = conf + 1 {
    let coords = self.get_coords(conf)
    for i = 0; i < self.n_atoms; i = i + 1 {
      if is_valid_coord(coords[i]) && counts[i] > 0 {
        let d = coords[i].sub(mean[i])
        rmsf[i] = rmsf[i] + d.dot(d)
      }
    }
  }
  for i = 0; i < self.n_atoms; i = i + 1 {
    if counts[i] > 0 {
      rmsf[i] = (rmsf[i] / counts[i].to_double()).sqrt()
    }
  }
  rmsf
}

///|
/// Mean structure (average coordinates over all conformations).
pub fn PDBEnsemble::get_mean_structure(self : PDBEnsemble) -> Array[Vector3] {
  let n_conf = self.coordsets.length()
  let mean : Array[Vector3] = Array::make(self.n_atoms, Vector3::zero())
  if n_conf == 0 {
    return self.ref_coords
  }
  let counts : Array[Int] = Array::make(self.n_atoms, 0)
  for conf = 0; conf < n_conf; conf = conf + 1 {
    let coords = self.get_coords(conf)
    for i = 0; i < self.n_atoms; i = i + 1 {
      if is_valid_coord(coords[i]) {
        mean[i] = mean[i].add(coords[i])
        counts[i] = counts[i] + 1
      }
    }
  }
  for i = 0; i < self.n_atoms; i = i + 1 {
    if counts[i] > 0 {
      mean[i] = mean[i].scale(1.0 / counts[i].to_double())
    } else {
      mean[i] = self.ref_coords[i]
    }
  }
  mean
}

///|
/// Per-atom coordinate deviation (standard deviation).
pub fn PDBEnsemble::get_coord_deviations(self : PDBEnsemble) -> Array[Double] {
  let n_conf = self.coordsets.length()
  let dev : Array[Double] = Array::make(self.n_atoms, 0.0)
  if n_conf == 0 {
    return dev
  }
  let mean = self.get_mean_structure()
  let counts : Array[Int] = Array::make(self.n_atoms, 0)
  for conf = 0; conf < n_conf; conf = conf + 1 {
    let coords = self.get_coords(conf)
    for i = 0; i < self.n_atoms; i = i + 1 {
      if is_valid_coord(coords[i]) {
        let d = coords[i].sub(mean[i])
        dev[i] = dev[i] + d.dot(d)
        counts[i] = counts[i] + 1
      }
    }
  }
  for i = 0; i < self.n_atoms; i = i + 1 {
    if counts[i] > 1 {
      dev[i] = (dev[i] / (counts[i] - 1).to_double()).sqrt()
    }
  }
  dev
}

///|
/// Set weights for superimposition.
pub fn PDBEnsemble::set_weights(self : PDBEnsemble, weights : Array[Double]) -> PDBEnsemble {
  let n = if weights.length() < self.n_atoms {
    weights.length()
  } else {
    self.n_atoms
  }
  let wts : Array[Double] = Array::make(self.n_atoms, 1.0)
  for i = 0; i < n; i = i + 1 {
    wts[i] = weights[i]
  }
  PDBEnsemble::{ ..self, weights: wts }
}

///|
/// Get the weights.
pub fn PDBEnsemble::get_weights(self : PDBEnsemble) -> Array[Double] {
  self.weights
}

///|
/// Remove a conformation by index.
pub fn PDBEnsemble::remove_conformation(self : PDBEnsemble, index : Int) -> PDBEnsemble {
  let new_coordsets : Array[Array[Vector3]] = Array(capacity=0)
  for i = 0; i < self.coordsets.length(); i = i + 1 {
    if i != index {
      new_coordsets.push(self.coordsets[i])
    }
  }
  PDBEnsemble::{ ..self, coordsets: new_coordsets, is_aligned: false }
}

///|
/// Get the atom names.
pub fn PDBEnsemble::get_atom_names(self : PDBEnsemble) -> Array[String] {
  self.atom_names
}

///|
/// Get the atom chains.
pub fn PDBEnsemble::get_atom_chains(self : PDBEnsemble) -> Array[Char] {
  self.atom_chains
}

///|
/// Get the atom residue numbers.
pub fn PDBEnsemble::get_atom_resseqs(self : PDBEnsemble) -> Array[Int] {
  self.atom_resseqs
}

///|
/// Coverage of each conformation (fraction of atoms present).
pub fn PDBEnsemble::get_coverage(self : PDBEnsemble) -> Array[Double] {
  let result : Array[Double] = Array(capacity=0)
  for conf = 0; conf < self.coordsets.length(); conf = conf + 1 {
    let mut n_valid = 0
    for i = 0; i < self.n_atoms; i = i + 1 {
      if is_valid_coord(self.coordsets[conf][i]) {
        n_valid = n_valid + 1
      }
    }
    result.push(n_valid.to_double() / self.n_atoms.to_double())
  }
  result
}

///|
/// Build a PDBEnsemble from multiple chains.
/// First chain is the reference; remaining chains are added as conformations.
pub fn PDBEnsemble::from_chains(chains : Array[Chain], calpha_only : Bool) -> PDBEnsemble? {
  if chains.length() == 0 {
    return None
  }
  match PDBEnsemble::new(chains[0], calpha_only) {
    Some(ens) => {
      let mut result = ens
      for i = 1; i < chains.length(); i = i + 1 {
        result = result.add_chain(chains[i])
      }
      Some(result)
    }
    None => None
  }
}

///|
/// Iterative superimposition: align, compute mean, re-align to mean.
pub fn PDBEnsemble::iterative_superimpose(self : PDBEnsemble, n_iter : Int) -> PDBEnsemble {
  let mut result = self.superimpose()
  for iter = 0; iter < n_iter; iter = iter + 1 {
    let mean = result.get_mean_structure()
    // Re-align each conformation to the mean
    let new_aligned : Array[Array[Vector3]] = Array(capacity=0)
    for conf = 0; conf < result.coordsets.length(); conf = conf + 1 {
      let (mapped_ref, mapped_mob) = get_mapped_coords_to_mean(result, conf, mean)
      if mapped_ref.length() < 3 {
        new_aligned.push(copy_coords(result.coordsets[conf]))
      } else {
        let (rot, trans) = vector_superposition_transform(mapped_mob, mapped_ref)
        let transformed : Array[Vector3] = Array::make(result.n_atoms, Vector3::new(999999.999, 999999.999, 999999.999))
        for i = 0; i < result.n_atoms; i = i + 1 {
          let c = result.coordsets[conf][i]
          if is_valid_coord(c) {
            transformed[i] = rot.transform(c).add(trans)
          }
        }
        new_aligned.push(transformed)
      }
    }
    result = PDBEnsemble::{ ..result, aligned: new_aligned, is_aligned: true }
  }
  result
}

///|
fn get_mapped_coords_to_mean(
  ens : PDBEnsemble,
  conf : Int,
  mean : Array[Vector3],
) -> (Array[Vector3], Array[Vector3]) {
  let ref : Array[Vector3] = Array(capacity=0)
  let mob : Array[Vector3] = Array(capacity=0)
  for i = 0; i < ens.n_atoms; i = i + 1 {
    let c = ens.coordsets[conf][i]
    if is_valid_coord(c) {
      ref.push(mean[i])
      mob.push(c)
    }
  }
  (ref, mob)
}
