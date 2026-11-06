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
