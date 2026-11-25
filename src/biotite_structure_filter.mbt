/// biotite_structure_filter.mbt — Filter structures by atom properties
/// Inspired by biotite.structure.filter_amino_acids, filter_backbone, etc.

///|
/// Atom metadata for filtering.
pub(all) struct BiotiteAtomInfo {
  atom_name : String
  res_name : String
  chain_id : String
  res_id : Int32
  hetero : Bool
}

///|
/// Default constructor for BiotiteAtomInfo.
pub fn BiotiteAtomInfo::default() -> BiotiteAtomInfo {
  { atom_name: "", res_name: "", chain_id: "", res_id: 0, hetero: false, }
}

///|
/// Standard amino acid three-letter codes.
fn biotite_is_amino_acid(res_name : String) -> Bool {
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
/// Standard nucleotide residue names.
fn biotite_is_nucleotide(res_name : String) -> Bool {
  match res_name {
    "DA"
    | "DC"
    | "DG"
    | "DT"
    | "DI"
    | "DU"
    | "A"
    | "C"
    | "G"
    | "T"
    | "I"
    | "U" => true
    _ => false
  }
}

///|
/// Backbone atom names.
fn biotite_is_backbone_atom(atom_name : String) -> Bool {
  match atom_name {
    "N" | "CA" | "C" | "O" | "CB" | "H" | "HA" | "OXT" => true
    _ => false
  }
}

///|
/// Filter to keep only amino acid residues.
pub fn biotite_filter_amino_acids(
  atoms : Array[BiotiteAtomInfo],
) -> Array[BiotiteAtomInfo] {
  let result : Array[BiotiteAtomInfo] = Array(capacity=0)
  for i = 0; i < atoms.length(); i = i + 1 {
    if biotite_is_amino_acids(atoms[i].res_name) && not(atoms[i].hetero) {
      result.push(atoms[i])
    }
  }
  result
}

///|
/// Filter to keep only backbone atoms of amino acid residues.
pub fn biotite_filter_backbone(
  atoms : Array[BiotiteAtomInfo],
) -> Array[BiotiteAtomInfo] {
  let result : Array[BiotiteAtomInfo] = Array(capacity=0)
  for i = 0; i < atoms.length(); i = i + 1 {
    if biotite_is_amino_acid(atoms[i].res_name) &&
      biotite_is_backbone_atom(atoms[i].atom_name) {
      result.push(atoms[i])
    }
  }
  result
}

///|
/// Filter to keep only nucleotide residues.
pub fn biotite_filter_nucleotides(
  atoms : Array[BiotiteAtomInfo],
) -> Array[BiotiteAtomInfo] {
  let result : Array[BiotiteAtomInfo] = Array(capacity=0)
  for i = 0; i < atoms.length(); i = i + 1 {
    if biotite_is_nucleotide(atoms[i].res_name) && not(atoms[i].hetero) {
      result.push(atoms[i])
    }
  }
  result
}

///|
/// Filter to keep only hetero atoms (non-standard residues).
pub fn biotite_filter_hetero(
  atoms : Array[BiotiteAtomInfo],
) -> Array[BiotiteAtomInfo] {
  let result : Array[BiotiteAtomInfo] = Array(capacity=0)
  for i = 0; i < atoms.length(); i = i + 1 {
    if atoms[i].hetero {
      result.push(atoms[i])
    }
  }
  result
}

///|
/// Filter by chain ID.
pub fn biotite_filter_chain(
  atoms : Array[BiotiteAtomInfo],
  chain_id : String,
) -> Array[BiotiteAtomInfo] {
  let result : Array[BiotiteAtomInfo] = Array(capacity=0)
  for i = 0; i < atoms.length(); i = i + 1 {
    if atoms[i].chain_id == chain_id {
      result.push(atoms[i])
    }
  }
  result
}

///|
/// Filter by residue name.
pub fn biotite_filter_residue(
  atoms : Array[BiotiteAtomInfo],
  res_name : String,
) -> Array[BiotiteAtomInfo] {
  let result : Array[BiotiteAtomInfo] = Array(capacity=0)
  for i = 0; i < atoms.length(); i = i + 1 {
    if atoms[i].res_name == res_name {
      result.push(atoms[i])
    }
  }
  result
}

///|
/// Get the starting indices of each chain.
pub fn biotite_get_chain_starts(atoms : Array[BiotiteAtomInfo]) -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  if atoms.length() == 0 {
    return result
  }
  result.push(0)
  for i = 1; i < atoms.length(); i = i + 1 {
    if atoms[i].chain_id != atoms[i - 1].chain_id {
      result.push(i)
    }
  }
  result
}

///|
/// Get the starting indices of each residue.
pub fn biotite_get_residue_starts(
  atoms : Array[BiotiteAtomInfo],
) -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  if atoms.length() == 0 {
    return result
  }
  result.push(0)
  for i = 1; i < atoms.length(); i = i + 1 {
    if atoms[i].res_id != atoms[i - 1].res_id ||
      atoms[i].chain_id != atoms[i - 1].chain_id {
      result.push(i)
    }
  }
  result
}

///|
/// Get the starting indices of each segment (contiguous run of same chain+res_name).
pub fn biotite_get_segment_starts(
  atoms : Array[BiotiteAtomInfo],
) -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  if atoms.length() == 0 {
    return result
  }
  result.push(0)
  for i = 1; i < atoms.length(); i = i + 1 {
    if atoms[i].chain_id != atoms[i - 1].chain_id ||
      atoms[i].res_name != atoms[i - 1].res_name {
      result.push(i)
    }
  }
  result
}

///|
/// Helper: negate boolean (not() is deprecated)
fn not(b : Bool) -> Bool {
  !b
}

///|
/// Helper: check amino acid with boolean return
fn biotite_is_amino_acids(res_name : String) -> Bool {
  biotite_is_amino_acid(res_name)
}

// ============================================================
// AtomArray adapters — reuse the repository's single structure model
// instead of maintaining a second one.
// ============================================================

///|
/// Project an AtomArray into the lightweight BiotiteAtomInfo view.
pub fn biotite_atom_info_from_array(arr : AtomArray) -> Array[BiotiteAtomInfo] {
  let result : Array[BiotiteAtomInfo] = Array(capacity=arr.len())
  for i in 0..<arr.len() {
    result.push({
      atom_name: arr.get_name(i),
      res_name: arr.get_res_name(i),
      chain_id: arr.get_chain_id(i).to_string(),
      res_id: arr.res_ids[i],
      hetero: false,
    })
  }
  result
}

///|
pub fn biotite_filter_amino_acids_array(arr : AtomArray) -> AtomArray {
  arr.filter_amino_acids()
}

///|
pub fn biotite_filter_backbone_array(arr : AtomArray) -> AtomArray {
  arr.filter_backbone()
}

///|
pub fn biotite_filter_nucleotides_array(arr : AtomArray) -> AtomArray {
  arr.filter_nucleotides()
}

///|
pub fn biotite_filter_chain_array(
  arr : AtomArray,
  chain_id : Char,
) -> AtomArray {
  arr.filter_chain(chain_id)
}

///|
pub fn biotite_filter_residue_array(
  arr : AtomArray,
  res_name : String,
) -> AtomArray {
  let mask : Array[Bool] = Array(capacity=arr.len())
  for i in 0..<arr.len() {
    mask.push(arr.get_res_name(i) == res_name)
  }
  arr.filter(mask)
}

///|
pub fn biotite_get_chain_starts_array(arr : AtomArray) -> Array[Int32] {
  arr.get_chain_starts()
}

///|
pub fn biotite_get_residue_starts_array(arr : AtomArray) -> Array[Int32] {
  arr.get_residue_starts()
}

///|
/// Segment starts: contiguous runs of the same chain + residue name.
pub fn biotite_get_segment_starts_array(arr : AtomArray) -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  if arr.len() == 0 {
    return result
  }
  result.push(0)
  for i in 1..<arr.len() {
    if arr.chain_ids[i] != arr.chain_ids[i - 1] ||
      arr.res_names[i] != arr.res_names[i - 1] {
      result.push(i)
    }
  }
  result
}
