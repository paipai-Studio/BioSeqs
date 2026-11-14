///|
/// AtomArray — unified NumPy-style structure representation.
///
/// Biotite's core abstraction: all atoms stored in parallel arrays
/// (coordinates, names, residue IDs, etc.) for efficient vectorized
/// operations. This module provides a MoonBit equivalent.

///|
/// Unified atom array — stores atom data in parallel arrays for
/// efficient filtering, indexing, and batch operations.
pub struct AtomArray {
  coords : Array[Vector3]
  names : Array[String]
  res_names : Array[String]
  chain_ids : Array[Char]
  res_ids : Array[Int32]
  elements : Array[String]
  b_factors : Array[Double]
  occupancies : Array[Double]
  altlocs : Array[Char]
}

///|
/// Create an empty AtomArray.
pub fn AtomArray::new() -> AtomArray {
  AtomArray::{
    coords: [],
    names: [],
    res_names: [],
    chain_ids: [],
    res_ids: [],
    elements: [],
    b_factors: [],
    occupancies: [],
    altlocs: [],
  }
}

///|
/// Create an AtomArray from an array of Atoms.
pub fn AtomArray::from_atoms(atoms : Array[Atom]) -> AtomArray {
  let coords : Array[Vector3] = []
  let names : Array[String] = []
  let res_names : Array[String] = []
  let chain_ids : Array[Char] = []
  let res_ids : Array[Int32] = []
  let elements : Array[String] = []
  let b_factors : Array[Double] = []
  let occupancies : Array[Double] = []
  let altlocs : Array[Char] = []
  for a in atoms {
    coords.push(a.get_coord())
    names.push(a.name)
    res_names.push(a.resname)
    chain_ids.push(a.chainid)
    res_ids.push(a.resseq)
    elements.push(a.element)
    b_factors.push(a.bfactor)
    occupancies.push(a.occupancy)
    altlocs.push(a.altloc)
  }
  AtomArray::{
    coords,
    names,
    res_names,
    chain_ids,
    res_ids,
    elements,
    b_factors,
    occupancies,
    altlocs,
  }
}

///|
/// Number of atoms.
pub fn AtomArray::len(self : AtomArray) -> Int32 {
  self.coords.length()
}

///|
/// Get coordinates of atom at index.
pub fn AtomArray::get_coord(self : AtomArray, i : Int32) -> Vector3 {
  self.coords[i]
}

///|
/// Get name of atom at index.
pub fn AtomArray::get_name(self : AtomArray, i : Int32) -> String {
  self.names[i]
}

///|
/// Get residue name of atom at index.
pub fn AtomArray::get_res_name(self : AtomArray, i : Int32) -> String {
  self.res_names[i]
}

///|
/// Get chain ID of atom at index.
pub fn AtomArray::get_chain_id(self : AtomArray, i : Int32) -> Char {
  self.chain_ids[i]
}

///|
/// Filter atoms by a boolean mask.
///
/// Returns a new AtomArray containing only atoms where mask[i] is true.
pub fn AtomArray::filter(self : AtomArray, mask : Array[Bool]) -> AtomArray {
  let coords : Array[Vector3] = []
  let names : Array[String] = []
  let res_names : Array[String] = []
  let chain_ids : Array[Char] = []
  let res_ids : Array[Int32] = []
  let elements : Array[String] = []
  let b_factors : Array[Double] = []
  let occupancies : Array[Double] = []
  let altlocs : Array[Char] = []
  for i in 0..<self.len() {
    if i < mask.length() && mask[i] {
      coords.push(self.coords[i])
      names.push(self.names[i])
      res_names.push(self.res_names[i])
      chain_ids.push(self.chain_ids[i])
      res_ids.push(self.res_ids[i])
      elements.push(self.elements[i])
      b_factors.push(self.b_factors[i])
      occupancies.push(self.occupancies[i])
      altlocs.push(self.altlocs[i])
    }
  }
  AtomArray::{
    coords,
    names,
    res_names,
    chain_ids,
    res_ids,
    elements,
    b_factors,
    occupancies,
    altlocs,
  }
}

///|
/// Filter to only C-alpha atoms.
pub fn AtomArray::filter_ca(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(self.names[i] == "CA")
  }
  self.filter(mask)
}

///|
/// Filter to only backbone atoms (N, CA, C, O).
pub fn AtomArray::filter_backbone(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    let n = self.names[i]
    mask.push(n == "N" || n == "CA" || n == "C" || n == "O")
  }
  self.filter(mask)
}

///|
/// Filter to a specific chain.
pub fn AtomArray::filter_chain(self : AtomArray, chain_id : Char) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(self.chain_ids[i] == chain_id)
  }
  self.filter(mask)
}

///|
/// Filter to heavy atoms (non-hydrogen).
pub fn AtomArray::filter_heavy(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(self.elements[i] != "H" && self.elements[i] != "D")
  }
  self.filter(mask)
}

///|
/// Concatenate two AtomArrays.
pub fn AtomArray::concatenate(self : AtomArray, other : AtomArray) -> AtomArray {
  let coords = aa_concat_v3(self.coords, other.coords)
  let names = aa_concat_str(self.names, other.names)
  let res_names = aa_concat_str(self.res_names, other.res_names)
  let chain_ids = aa_concat_char(self.chain_ids, other.chain_ids)
  let res_ids = aa_concat_i32(self.res_ids, other.res_ids)
  let elements = aa_concat_str(self.elements, other.elements)
  let b_factors = aa_concat_dbl(self.b_factors, other.b_factors)
  let occupancies = aa_concat_dbl(self.occupancies, other.occupancies)
  let altlocs = aa_concat_char(self.altlocs, other.altlocs)
  AtomArray::{
    coords,
    names,
    res_names,
    chain_ids,
    res_ids,
    elements,
    b_factors,
    occupancies,
    altlocs,
  }
}

///|
/// Convert AtomArray back to Array[Atom].
pub fn AtomArray::to_atoms(self : AtomArray) -> Array[Atom] {
  let atoms : Array[Atom] = []
  for i in 0..<self.len() {
    atoms.push(
      Atom::new(
        name=self.names[i],
        coord=self.coords[i],
        resname=self.res_names[i],
        chainid=self.chain_ids[i],
        resseq=self.res_ids[i],
        element=self.elements[i],
        bfactor=self.b_factors[i],
        occupancy=self.occupancies[i],
      ),
    )
  }
  atoms
}

///|
fn aa_concat_v3(a : Array[Vector3], b : Array[Vector3]) -> Array[Vector3] {
  let r = Array(capacity=a.length() + b.length())
  for v in a {
    r.push(v)
  }
  for v in b {
    r.push(v)
  }
  r
}

///|
fn aa_concat_str(a : Array[String], b : Array[String]) -> Array[String] {
  let r = Array(capacity=a.length() + b.length())
  for v in a {
    r.push(v)
  }
  for v in b {
    r.push(v)
  }
  r
}

///|
fn aa_concat_char(a : Array[Char], b : Array[Char]) -> Array[Char] {
  let r = Array(capacity=a.length() + b.length())
  for v in a {
    r.push(v)
  }
  for v in b {
    r.push(v)
  }
  r
}

///|
fn aa_concat_i32(a : Array[Int32], b : Array[Int32]) -> Array[Int32] {
  let r = Array(capacity=a.length() + b.length())
  for v in a {
    r.push(v)
  }
  for v in b {
    r.push(v)
  }
  r
}

///|
fn aa_concat_dbl(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let r = Array(capacity=a.length() + b.length())
  for v in a {
    r.push(v)
  }
  for v in b {
    r.push(v)
  }
  r
}

// ============================================================
// P0 — Residue / Chain level tool functions (Biotite pattern)
// ============================================================

///|
/// Get starting atom index of each residue.
/// A residue boundary is detected when chain_id or res_id changes.
pub fn AtomArray::get_residue_starts(self : AtomArray) -> Array[Int32] {
  let n = self.len()
  let starts : Array[Int32] = Array(capacity=n)
  if n == 0 {
    return starts
  }
  starts.push(0)
  for i in 1..<n {
    if self.chain_ids[i] != self.chain_ids[i - 1] ||
      self.res_ids[i] != self.res_ids[i - 1] {
      starts.push(i)
    }
  }
  starts
}

///|
/// Get number of residues in the atom array.
pub fn AtomArray::get_residue_count(self : AtomArray) -> Int32 {
  if self.len() == 0 {
    return 0
  }
  self.get_residue_starts().length()
}

///|
/// Get the number of atoms in each residue.
pub fn AtomArray::get_residue_lengths(self : AtomArray) -> Array[Int32] {
  let starts = self.get_residue_starts()
  let n = self.len()
  let lengths : Array[Int32] = Array(capacity=starts.length())
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    lengths.push(e - s)
  }
  lengths
}

///|
/// Get end atom index (exclusive) of each residue.
pub fn AtomArray::get_residue_ends(self : AtomArray) -> Array[Int32] {
  let starts = self.get_residue_starts()
  let n = self.len()
  let ends : Array[Int32] = Array(capacity=starts.length())
  for i in 0..<starts.length() {
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    ends.push(e)
  }
  ends
}

///|
/// Spread per-residue Double data to per-atom level.
/// residue_data must have exactly one element per residue.
pub fn spread_residue_double(
  arr : AtomArray,
  residue_data : Array[Double],
) -> Array[Double] {
  let starts = arr.get_residue_starts()
  let n = arr.len()
  let result : Array[Double] = Array(capacity=n)
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    let val = residue_data[i]
    for _ in s..<e {
      result.push(val)
    }
  }
  result
}

///|
/// Spread per-residue Bool data to per-atom level.
pub fn spread_residue_bool(
  arr : AtomArray,
  residue_data : Array[Bool],
) -> Array[Bool] {
  let starts = arr.get_residue_starts()
  let n = arr.len()
  let result : Array[Bool] = Array(capacity=n)
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    let val = residue_data[i]
    for _ in s..<e {
      result.push(val)
    }
  }
  result
}

///|
/// Get a boolean mask for each residue.
/// Returns Array[Array[Bool]], one mask per residue.
pub fn AtomArray::get_residue_masks(self : AtomArray) -> Array[Array[Bool]] {
  let starts = self.get_residue_starts()
  let n = self.len()
  let masks : Array[Array[Bool]] = Array(capacity=starts.length())
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    let mask : Array[Bool] = Array(capacity=n)
    for _ in 0..<n {
      mask.push(false)
    }
    for j in s..<e {
      mask[j] = true
    }
    masks.push(mask)
  }
  masks
}

// --- Chain-level ---

///|
/// Get starting atom index of each chain.
/// A chain boundary is detected when chain_id changes.
pub fn AtomArray::get_chain_starts(self : AtomArray) -> Array[Int32] {
  let n = self.len()
  let starts : Array[Int32] = Array(capacity=n)
  if n == 0 {
    return starts
  }
  starts.push(0)
  for i in 1..<n {
    if self.chain_ids[i] != self.chain_ids[i - 1] {
      starts.push(i)
    }
  }
  starts
}

///|
/// Get number of chains in the atom array.
pub fn AtomArray::get_chain_count(self : AtomArray) -> Int32 {
  if self.len() == 0 {
    return 0
  }
  self.get_chain_starts().length()
}

///|
/// Get the number of atoms in each chain.
pub fn AtomArray::get_chain_lengths(self : AtomArray) -> Array[Int32] {
  let starts = self.get_chain_starts()
  let n = self.len()
  let lengths : Array[Int32] = Array(capacity=starts.length())
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    lengths.push(e - s)
  }
  lengths
}

///|
/// Get end atom index (exclusive) of each chain.
pub fn AtomArray::get_chain_ends(self : AtomArray) -> Array[Int32] {
  let starts = self.get_chain_starts()
  let n = self.len()
  let ends : Array[Int32] = Array(capacity=starts.length())
  for i in 0..<starts.length() {
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    ends.push(e)
  }
  ends
}

///|
/// Get a boolean mask for each chain.
pub fn AtomArray::get_chain_masks(self : AtomArray) -> Array[Array[Bool]] {
  let starts = self.get_chain_starts()
  let n = self.len()
  let masks : Array[Array[Bool]] = Array(capacity=starts.length())
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    let mask : Array[Bool] = Array(capacity=n)
    for _ in 0..<n {
      mask.push(false)
    }
    for j in s..<e {
      mask[j] = true
    }
    masks.push(mask)
  }
  masks
}

// ============================================================
// P0 — Filter family (Biotite pattern)
// ============================================================

///|
/// 20 canonical amino acid three-letter codes.
fn is_canonical_aa(res_name : String) -> Bool {
  match res_name {
    "ALA"
    | "ARG"
    | "ASN"
    | "ASP"
    | "CYS"
    | "GLN"
    | "GLU"
    | "GLY"
    | "HIS"
    | "ILE"
    | "LEU"
    | "LYS"
    | "MET"
    | "PHE"
    | "PRO"
    | "SER"
    | "THR"
    | "TRP"
    | "TYR"
    | "VAL" => true
    _ => false
  }
}

///|
/// All amino acid three-letter codes (canonical + selenocysteine + pyrrolysine + oxidized MSE).
fn is_any_amino_acid(res_name : String) -> Bool {
  match res_name {
    "ALA"
    | "ARG"
    | "ASN"
    | "ASP"
    | "CYS"
    | "GLN"
    | "GLU"
    | "GLY"
    | "HIS"
    | "ILE"
    | "LEU"
    | "LYS"
    | "MET"
    | "PHE"
    | "PRO"
    | "SER"
    | "THR"
    | "TRP"
    | "TYR"
    | "VAL"
    | "MSE"
    | "SEC"
    | "PYL" => true
    _ => false
  }
}

///|
/// Standard nucleotide residue names (DNA + RNA).
fn is_canonical_nucleotide(res_name : String) -> Bool {
  match res_name {
    "DA" | "DC" | "DG" | "DT" | "DU" => true
    _ => false
  }
}

///|
/// All nucleotide residue names (canonical + one-letter + modified).
fn is_any_nucleotide(res_name : String) -> Bool {
  match res_name {
    "DA"
    | "DC"
    | "DG"
    | "DT"
    | "DU"
    | "DI"
    | "A"
    | "C"
    | "G"
    | "T"
    | "U"
    | "I" => true
    _ => false
  }
}

///|
/// Water / solvent residue names.
fn is_solvent(res_name : String) -> Bool {
  match res_name {
    "WAT" | "HOH" | "H2O" | "TIP3" | "TIP4" | "SOL" | "DOD" | "HO4" | "W01" =>
      true
    _ => false
  }
}

///|
/// Monoatomic ion residue names.
fn is_monoatomic_ion(res_name : String) -> Bool {
  match res_name {
    "NA"
    | "CL"
    | "K"
    | "ZN"
    | "MG"
    | "CA"
    | "FE"
    | "MN"
    | "CU"
    | "CO"
    | "NI"
    | "CS"
    | "RB"
    | "SR"
    | "BA"
    | "AL"
    | "LI" => true
    _ => false
  }
}

///|
/// Carbohydrate residue names.
fn is_carbohydrate(res_name : String) -> Bool {
  match res_name {
    "GLC"
    | "FRU"
    | "GAL"
    | "MAL"
    | "SUC"
    | "TRE"
    | "RIB"
    | "XYL"
    | "ARA"
    | "NAG"
    | "BMA"
    | "BGC"
    | "GCU"
    | "GCR"
    | "MAN"
    | "G6P"
    | "F6P" => true
    _ => false
  }
}

///|
/// Peptide backbone atom names (protein).
fn is_peptide_backbone(atom_name : String) -> Bool {
  match atom_name {
    "N" | "CA" | "C" | "O" => true
    _ => false
  }
}

///|
/// Nucleotide phosphate/sugar backbone atom names.
fn is_phosphate_backbone(atom_name : String) -> Bool {
  match atom_name {
    "P" | "OP1" | "OP2" | "O5'" | "O3'" | "O4'" => true
    _ => false
  }
}

///|
/// Filter to only canonical 20 amino acids (excludes hetero / modified residues).
pub fn AtomArray::filter_canonical_amino_acids(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_canonical_aa(self.res_names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to all amino acids including modified (MSE, SEC, PYL).
pub fn AtomArray::filter_amino_acids(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_any_amino_acid(self.res_names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to only canonical nucleotides (DA/DC/DG/DT/DU).
pub fn AtomArray::filter_canonical_nucleotides(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_canonical_nucleotide(self.res_names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to all nucleotides.
pub fn AtomArray::filter_nucleotides(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_any_nucleotide(self.res_names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to peptide backbone atoms (N, CA, C, O).
pub fn AtomArray::filter_peptide_backbone(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_peptide_backbone(self.names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to nucleotide phosphate/sugar backbone atoms.
pub fn AtomArray::filter_phosphate_backbone(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_phosphate_backbone(self.names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to solvent / water molecules.
pub fn AtomArray::filter_solvent(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_solvent(self.res_names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to monoatomic ions.
pub fn AtomArray::filter_monoatomic_ions(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_monoatomic_ion(self.res_names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to carbohydrate residues.
pub fn AtomArray::filter_carbohydrates(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    mask.push(is_carbohydrate(self.res_names[i]))
  }
  self.filter(mask)
}

///|
/// Filter to first altloc per residue.
/// Uses altloc char: space (no altloc) or 'A' is kept, others removed.
pub fn AtomArray::filter_first_altloc(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    let a = self.altlocs[i]
    mask.push(a == ' ' || a == 'A')
  }
  self.filter(mask)
}

///|
/// Filter to highest-occupancy altloc per residue.
/// For each residue group, finds the altloc with max occupancy and keeps only those atoms.
pub fn AtomArray::filter_highest_occupancy_altloc(
  self : AtomArray,
) -> AtomArray {
  let starts = self.get_residue_starts()
  let n = self.len()
  let mask : Array[Bool] = Array(capacity=n)
  for _ in 0..<n {
    mask.push(true)
  }
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    let mut has_altloc = false
    let mut max_occ = -1.0
    let mut best_altloc = ' '
    for j in s..<e {
      let a = self.altlocs[j]
      if a != ' ' {
        has_altloc = true
        if self.occupancies[j] > max_occ {
          max_occ = self.occupancies[j]
          best_altloc = a
        }
      }
    }
    if has_altloc {
      for j in s..<e {
        let a = self.altlocs[j]
        if a != ' ' && a != best_altloc {
          mask[j] = false
        }
      }
    }
  }
  self.filter(mask)
}

///|
/// Filter atoms present in both AtomArrays.
/// An atom matches if chain_id, res_id, res_name and name all equal.
pub fn AtomArray::filter_intersection(
  self : AtomArray,
  other : AtomArray,
) -> AtomArray {
  let n_other = other.len()
  let keys_other : Array[String] = Array(capacity=n_other)
  for idx in 0..<n_other {
    let key = "${other.chain_ids[idx]}-${other.res_ids[idx]}-${other.res_names[idx]}-${other.names[idx]}"
    keys_other.push(key)
  }
  let mask : Array[Bool] = []
  for idx2 in 0..<self.len() {
    let key = "${self.chain_ids[idx2]}-${self.res_ids[idx2]}-${self.res_names[idx2]}-${self.names[idx2]}"
    let mut found = false
    for k in keys_other {
      if k == key {
        found = true
        break
      }
    }
    mask.push(found)
  }
  self.filter(mask)
}

// ============================================================
// P0 — AtomArrayStack (multi-model structure)
// ============================================================

///|
/// A stack of AtomArrays sharing the same annotation arrays
/// but having different coordinate sets.
/// Used for NMR multi-model structures, MD trajectories, etc.
pub struct AtomArrayStack {
  coords_stack : Array[Array[Vector3]]
  names : Array[String]
  res_names : Array[String]
  chain_ids : Array[Char]
  res_ids : Array[Int32]
  elements : Array[String]
  b_factors : Array[Double]
  occupancies : Array[Double]
  altlocs : Array[Char]
}

///|
/// Create an AtomArrayStack from an AtomArray (single model).
pub fn AtomArrayStack::from_atom_array(arr : AtomArray) -> AtomArrayStack {
  let coords : Array[Vector3] = Array(capacity=arr.coords.length())
  for c in arr.coords {
    coords.push(c)
  }
  AtomArrayStack::{
    coords_stack: [coords],
    names: arr.names,
    res_names: arr.res_names,
    chain_ids: arr.chain_ids,
    res_ids: arr.res_ids,
    elements: arr.elements,
    b_factors: arr.b_factors,
    occupancies: arr.occupancies,
    altlocs: arr.altlocs,
  }
}

///|
/// Create an AtomArrayStack from annotation arrays and multiple coordinate sets.
pub fn AtomArrayStack::new(
  coords_stack : Array[Array[Vector3]],
  names : Array[String],
  res_names : Array[String],
  chain_ids : Array[Char],
  res_ids : Array[Int32],
  elements : Array[String],
  b_factors : Array[Double],
  occupancies : Array[Double],
  altlocs : Array[Char],
) -> AtomArrayStack {
  AtomArrayStack::{
    coords_stack,
    names,
    res_names,
    chain_ids,
    res_ids,
    elements,
    b_factors,
    occupancies,
    altlocs,
  }
}

///|
/// Number of models in the stack.
pub fn AtomArrayStack::n_models(self : AtomArrayStack) -> Int32 {
  self.coords_stack.length()
}

///|
/// Number of atoms per model.
pub fn AtomArrayStack::len(self : AtomArrayStack) -> Int32 {
  self.names.length()
}

///|
/// Extract one model as an AtomArray.
pub fn AtomArrayStack::get_model(
  self : AtomArrayStack,
  idx : Int32,
) -> AtomArray {
  let coords = self.coords_stack[idx]
  AtomArray::{
    coords,
    names: self.names,
    res_names: self.res_names,
    chain_ids: self.chain_ids,
    res_ids: self.res_ids,
    elements: self.elements,
    b_factors: self.b_factors,
    occupancies: self.occupancies,
    altlocs: self.altlocs,
  }
}

///|
/// Append a new model (coordinate set) to the stack.
pub fn AtomArrayStack::append_model(
  self : AtomArrayStack,
  coords : Array[Vector3],
) -> AtomArrayStack {
  let new_stack : Array[Array[Vector3]] = []
  for c in self.coords_stack {
    new_stack.push(c)
  }
  new_stack.push(coords)
  AtomArrayStack::{
    coords_stack: new_stack,
    names: self.names,
    res_names: self.res_names,
    chain_ids: self.chain_ids,
    res_ids: self.res_ids,
    elements: self.elements,
    b_factors: self.b_factors,
    occupancies: self.occupancies,
    altlocs: self.altlocs,
  }
}

///|
/// Concatenate multiple AtomArrayStacks along the model axis.
/// All stacks must have the same atom annotations.
pub fn AtomArrayStack::concatenate_models(
  stacks : Array[AtomArrayStack],
) -> AtomArrayStack {
  let first = stacks[0]
  let new_stack : Array[Array[Vector3]] = []
  for s in stacks {
    for c in s.coords_stack {
      new_stack.push(c)
    }
  }
  AtomArrayStack::{
    coords_stack: new_stack,
    names: first.names,
    res_names: first.res_names,
    chain_ids: first.chain_ids,
    res_ids: first.res_ids,
    elements: first.elements,
    b_factors: first.b_factors,
    occupancies: first.occupancies,
    altlocs: first.altlocs,
  }
}

// ============================================================
// P1 — Structure analysis metrics (Biotite pattern)
// ============================================================

// --- P1-1: Mass center + radius of gyration ---

///|
/// Compute the total mass of all atoms using atomic element masses.
pub fn AtomArray::total_mass(self : AtomArray) -> Double {
  let mp = calculate_mass_properties(self.coords, self.elements)
  mp.get_total_mass()
}

///|
/// Compute the mass center (center of mass) of all atoms.
/// Uses atomic element masses from mass_charge.
pub fn AtomArray::mass_center(self : AtomArray) -> Vector3 {
  let mp = calculate_mass_properties(self.coords, self.elements)
  mp.get_center_of_mass()
}

///|
/// Compute the radius of gyration of all atoms.
/// Rg = sqrt(Σ mi * |ri - rc|² / Σ mi).
pub fn AtomArray::gyration_radius(self : AtomArray) -> Double {
  let mp = calculate_mass_properties(self.coords, self.elements)
  mp.get_radius_of_gyration()
}

// --- P1-2: Secondary structure assignment (PSEA adapter) ---

///|
/// Extract alpha-carbon atoms and run PSEA 3-state secondary structure
/// prediction. Residues without a CA atom are skipped.
pub fn AtomArray::annotate_sse(self : AtomArray) -> PseaResult {
  // Group CA atoms by chain + residue order
  let chain_order : Array[Char] = []
  let res_seen : Array[(Char, Int32)] = []
  let cas : Array[(Int32, Int32)] = [] // (res_idx, atom_idx) pairs
  let starts = self.get_residue_starts()
  for i in 0..<starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { self.len() }
    let chain = self.chain_ids[s]
    let rid = self.res_ids[s]
    // Look for CA in this residue
    let mut ca_idx = -1
    for j in s..<e {
      if self.names[j] == "CA" {
        ca_idx = j
        break
      }
    }
    if ca_idx >= 0 {
      cas.push((i, ca_idx))
      // Track chain order
      if !(chain_order.length() > 0 &&
        chain_order[chain_order.length() - 1] == chain) {
        chain_order.push(chain)
      }
      res_seen.push((chain, rid))
    }
  }
  // Build PseaAtom list, grouping by chain
  let psea_atoms : Array[PseaAtom] = Array(capacity=cas.length())
  let mut prev_chain : Char = if cas.length() > 0 {
    self.chain_ids[cas[0].1]
  } else {
    ' '
  }
  let mut seq_num : Int32 = 0
  let mut k = 0
  while k < cas.length() {
    let (res_i, atom_i) = cas[k]
    let chain = self.chain_ids[atom_i]
    if chain != prev_chain {
      seq_num = 0 // restart chain-local numbering
      prev_chain = chain
    }
    let coord = self.coords[atom_i]
    let rn = self.res_names[atom_i]
    let rid = self.res_ids[atom_i]
    // Prefer actual residue seq number when available
    let resnum = if rid > 0 { rid } else { seq_num + 1 }
    psea_atoms.push(PseaAtom::new(rn, resnum, "CA", coord.x, coord.y, coord.z))
    seq_num = seq_num + 1
    k = k + 1
  }
  psea_run(psea_atoms)
}

// --- P1-4: Orient by principal component analysis ---

///|
/// Reorient the AtomArray so that its principal axes align with
/// the Cartesian x/y/z axes (PCA alignment).
///
/// Uses unweighted covariance (all atoms contribute equally).
/// The center of mass is moved to the origin and the coordinate
/// system is rotated so that:
///   x-axis ← largest variance direction
///   y-axis ← second largest variance direction
///   z-axis ← smallest variance direction
///
/// Returns a new AtomArray with updated coordinates.
pub fn AtomArray::orient_principal_components(self : AtomArray) -> AtomArray {
  let n = self.len()
  if n < 3 {
    return self
  }
  // 1. Compute centroid
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for i in 0..<n {
    cx = cx + self.coords[i].x
    cy = cy + self.coords[i].y
    cz = cz + self.coords[i].z
  }
  let nf = n.to_double()
  cx = cx / nf
  cy = cy / nf
  cz = cz / nf
  // 2. Build 3x3 covariance matrix (unweighted)
  // C = (1/n) * Σ (r - c)(r - c)^T
  let mut cxx = 0.0
  let mut cyy = 0.0
  let mut czz = 0.0
  let mut cxy = 0.0
  let mut cxz = 0.0
  let mut cyz = 0.0
  for i in 0..<n {
    let dx = self.coords[i].x - cx
    let dy = self.coords[i].y - cy
    let dz = self.coords[i].z - cz
    cxx = cxx + dx * dx
    cyy = cyy + dy * dy
    czz = czz + dz * dz
    cxy = cxy + dx * dy
    cxz = cxz + dx * dz
    cyz = cyz + dy * dz
  }
  cxx = cxx / nf
  cyy = cyy / nf
  czz = czz / nf
  cxy = cxy / nf
  cxz = cxz / nf
  cyz = cyz / nf
  let matrix : Array[Array[Double]] = [
    [cxx, cxy, cxz],
    [cxy, cyy, cyz],
    [cxz, cyz, czz],
  ]
  // 3. Jacobi eigendecomposition (3x3)
  let (eigenvalues, eigenvectors) = linalg_jacobi(matrix, 3, 100)
  // 4. Sort descending
  let (_, v_sorted) = linalg_sort_eigen(eigenvalues, eigenvectors, 3)
  // eigenvectors columns are eigenvectors. Build rotation:
  // V^T rows = original V columns
  // v_sorted[row][col] — we want row [0][*], [1][*], [2][*] to be the
  // eigenvectors (V^T). linalg_sort_eigen reorders columns. v_sorted[i][j]
  // is element row i, column j of sorted eigenvector matrix.
  // So eigenvector for j-th largest eigenvalue = [v_sorted[0][j], v_sorted[1][j], v_sorted[2][j]]
  // For rotation V^T, we use these as rows.
  let e0x = v_sorted[0][0]
  let e0y = v_sorted[1][0]
  let e0z = v_sorted[2][0]
  let e1x = v_sorted[0][1]
  let e1y = v_sorted[1][1]
  let e1z = v_sorted[2][1]
  let e2x = v_sorted[0][2]
  let e2y = v_sorted[1][2]
  let e2z = v_sorted[2][2]
  // 5. Rotate + translate to origin
  let new_coords : Array[Vector3] = Array(capacity=n)
  for i in 0..<n {
    let dx = self.coords[i].x - cx
    let dy = self.coords[i].y - cy
    let dz = self.coords[i].z - cz
    let nx = e0x * dx + e0y * dy + e0z * dz
    let ny = e1x * dx + e1y * dy + e1z * dz
    let nz = e2x * dx + e2y * dy + e2z * dz
    new_coords.push(Vector3::new(nx, ny, nz))
  }
  AtomArray::{
    coords: new_coords,
    names: self.names,
    res_names: self.res_names,
    chain_ids: self.chain_ids,
    res_ids: self.res_ids,
    elements: self.elements,
    b_factors: self.b_factors,
    occupancies: self.occupancies,
    altlocs: self.altlocs,
  }
}

// ============================================================
// P2 — Bond connectivity (Biotite pattern)
// ============================================================

///|
/// Infer bond connectivity from coordinates + elements using covalent
/// radii distance cutoffs. Wraps `bondlist_from_coordinates`.
pub fn AtomArray::compute_bonds(self : AtomArray) -> BondList {
  bondlist_from_coordinates(self.coords, self.elements)
}

///|
/// Find rotatable single bonds (non-terminal atoms, non-rings).
/// Infers connectivity first if not precomputed.
pub fn AtomArray::find_rotatable_bonds(self : AtomArray) -> Array[Bond] {
  let bl = bondlist_from_coordinates(self.coords, self.elements)
  bl.find_rotatable_bonds()
}

///|
/// Get all neighbors of an atom via inferred bond graph.
pub fn AtomArray::get_neighbors(
  self : AtomArray,
  atom_idx : Int32,
) -> Array[Int32] {
  let bl = bondlist_from_coordinates(self.coords, self.elements)
  bl.get_neighbors(atom_idx)
}

// ============================================================
// P2 — Box / unit cell helpers (Biotite pattern)
// ============================================================

///|
/// Compute the axis-aligned bounding box of the AtomArray.
/// Returns a `BiotiteBox` from orthogonal edge vectors.
pub fn AtomArray::get_bounding_box(self : AtomArray) -> BiotiteBox {
  let n = self.len()
  if n == 0 {
    return BiotiteBox::from_vectors([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [
      0.0, 0.0, 0.0,
    ])
  }
  let mut minx = self.coords[0].x
  let mut miny = self.coords[0].y
  let mut minz = self.coords[0].z
  let mut maxx = minx
  let mut maxy = miny
  let mut maxz = minz
  let mut i = 1
  while i < n {
    let c = self.coords[i]
    if c.x < minx {
      minx = c.x
    } else if c.x > maxx {
      maxx = c.x
    }
    if c.y < miny {
      miny = c.y
    } else if c.y > maxy {
      maxy = c.y
    }
    if c.z < minz {
      minz = c.z
    } else if c.z > maxz {
      maxz = c.z
    }
    i = i + 1
  }
  let ax = maxx - minx
  let by = maxy - miny
  let cz = maxz - minz
  BiotiteBox::from_vectors([ax, 0.0, 0.0], [0.0, by, 0.0], [0.0, 0.0, cz])
}

///|
/// Translate the AtomArray so the geometric center moves to the origin.
pub fn AtomArray::center(self : AtomArray) -> AtomArray {
  let n = self.len()
  if n == 0 {
    return self
  }
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  let mut i = 0
  while i < n {
    cx = cx + self.coords[i].x
    cy = cy + self.coords[i].y
    cz = cz + self.coords[i].z
    i = i + 1
  }
  let nf = n.to_double()
  cx = cx / nf
  cy = cy / nf
  cz = cz / nf
  let new_coords : Array[Vector3] = Array(capacity=n)
  let mut j = 0
  while j < n {
    new_coords.push(
      Vector3::new(
        self.coords[j].x - cx,
        self.coords[j].y - cy,
        self.coords[j].z - cz,
      ),
    )
    j = j + 1
  }
  AtomArray::{
    coords: new_coords,
    names: self.names,
    res_names: self.res_names,
    chain_ids: self.chain_ids,
    res_ids: self.res_ids,
    elements: self.elements,
    b_factors: self.b_factors,
    occupancies: self.occupancies,
    altlocs: self.altlocs,
  }
}
