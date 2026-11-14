///|
/// Structure check + repair for AtomArray.
///
/// Inspired by Biotite's `structure.repair` submodule.
/// All functions operate on the AtomArray directly (in-place returns
/// a new AtomArray, since MoonBit structs are values).

// ============================================================
// P2-1: Infer element from atom name
// ============================================================

///|
/// Infer an element symbol from a PDB atom name.
///
/// PDB convention: atom names encode element by position.
/// - Column 1 (index 0): element symbol right-aligned if len(name) <= 3
/// - Otherwise the first non-digit character after stripping digits
///
/// Returns "X" if unknown.
pub fn infer_element(atom_name : String) -> String {
  // Try exact metal/halogen/phosphorus 2-letter names first
  // Strip trailing digits
  let mut stripped = atom_name
  while stripped.length() > 0 {
    let last = stripped[stripped.length() - 1]
    if last >= '0' && last <= '9' {
      stripped = stripped[:stripped.length() - 1].to_string()
    } else {
      break
    }
  }
  // Match common 2-letter element names
  match stripped {
    "FE"
    | "ZN"
    | "MG"
    | "NA"
    | "CL"
    | "BR"
    | "CU"
    | "CO"
    | "NI"
    | "MN"
    | "CA"
    | "CS"
    | "RB"
    | "SR"
    | "BA"
    | "AL"
    | "LI"
    | "SI"
    | "SE"
    | "AS"
    | "GE"
    | "GA"
    | "AU"
    | "AG"
    | "HG"
    | "CD"
    | "PD"
    | "PT"
    | "OS"
    | "IR"
    | "RH"
    | "RU"
    | "TC"
    | "MO"
    | "TI"
    | "CR"
    | "PB"
    | "BI"
    | "SB"
    | "H"
    | "C"
    | "N"
    | "O"
    | "S"
    | "P"
    | "F"
    | "V" => stripped
    _ => {
      // First letter uppercase is usually the element
      if stripped.length() == 0 {
        return "X"
      }
      let first = stripped[0]
      if first >= '1' && first <= '9' {
        // PDB atom name like "1HG" → element starts at index 1
        if stripped.length() >= 2 {
          let f2 = stripped[1]
          return f2.to_string()
        }
        return "X"
      }
      first.to_string()
    }
  }
}

///|
/// Fill in missing element fields (empty or "X") using atom name inference.
/// Returns a new AtomArray with element array updated.
pub fn AtomArray::infer_elements(self : AtomArray) -> AtomArray {
  let new_elements : Array[String] = Array(capacity=self.len())
  for i in 0..<self.len() {
    let e = self.elements[i]
    if e == "" || e == "X" {
      new_elements.push(infer_element(self.names[i]))
    } else {
      new_elements.push(e)
    }
  }
  AtomArray::{
    coords: self.coords,
    names: self.names,
    res_names: self.res_names,
    chain_ids: self.chain_ids,
    res_ids: self.res_ids,
    elements: new_elements,
    b_factors: self.b_factors,
    occupancies: self.occupancies,
    altlocs: self.altlocs,
  }
}

// ============================================================
// P2-2: Residue ID continuity
// ============================================================

///|
/// Check whether residue IDs form a continuous run within each chain.
///
/// Returns an array of (chain_id, gap_size) pairs describing where
/// res_id jumps by more than 1.
pub fn AtomArray::check_res_id_continuity(
  self : AtomArray,
) -> Array[(Char, Int32, Int32, Int32)] {
  // (chain_id, from_res_id, to_res_id, gap)
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  let gaps : Array[(Char, Int32, Int32, Int32)] = []
  let mut i = 1
  while i < n_res {
    let prev_s = starts[i - 1]
    let cur_s = starts[i]
    let prev_chain = self.chain_ids[prev_s]
    let cur_chain = self.chain_ids[cur_s]
    // Only check within same chain
    if prev_chain == cur_chain {
      let prev_rid = self.res_ids[prev_s]
      let cur_rid = self.res_ids[cur_s]
      let diff = cur_rid - prev_rid
      if diff > 1 {
        gaps.push((prev_chain, prev_rid, cur_rid, diff - 1))
      }
    }
    i = i + 1
  }
  gaps
}

///|
/// Reassign residue IDs to be continuous within each chain.
/// The first residue in each chain keeps its original ID; subsequent
/// residues are assigned sequentially (+1 each).
pub fn AtomArray::create_continuous_res_ids(self : AtomArray) -> AtomArray {
  let starts = self.get_residue_starts()
  let n = self.len()
  let new_res_ids : Array[Int32] = Array(capacity=n)
  let new_resnames : Array[String] = Array(capacity=n)
  let new_chainids : Array[Char] = Array(capacity=n)
  // For each residue, build new res_id
  let new_ids : Array[Int32] = []
  let mut prev_chain : Char = if starts.length() > 0 {
    self.chain_ids[starts[0]]
  } else {
    ' '
  }
  let mut prev_id = if starts.length() > 0 {
    self.res_ids[starts[0]] - 1
  } else {
    0
  }
  let mut k = 0
  while k < starts.length() {
    let s = starts[k]
    let chain = self.chain_ids[s]
    if chain != prev_chain {
      prev_id = self.res_ids[s] - 1
      prev_chain = chain
    }
    prev_id = prev_id + 1
    new_ids.push(prev_id)
    k = k + 1
  }
  // Spread to per-atom arrays
  let mut i = 0
  while i < starts.length() {
    let s = starts[i]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { n }
    let rid = new_ids[i]
    let rn = self.res_names[s]
    let ch = self.chain_ids[s]
    let mut j = s
    while j < e {
      new_res_ids.push(rid)
      new_resnames.push(rn)
      new_chainids.push(ch)
      j = j + 1
    }
    i = i + 1
  }
  AtomArray::{
    coords: self.coords,
    names: self.names,
    res_names: new_resnames,
    chain_ids: new_chainids,
    res_ids: new_res_ids,
    elements: self.elements,
    b_factors: self.b_factors,
    occupancies: self.occupancies,
    altlocs: self.altlocs,
  }
}

// ============================================================
// P2-3: Backbone continuity check
// ============================================================

///|
/// Standard peptide backbone atom names.
fn standard_backbone_atoms() -> Array[String] {
  ["N", "CA", "C", "O"]
}

///|
/// Standard nucleotide phosphate/sugar backbone atom names.
fn standard_phosphate_backbone_atoms() -> Array[String] {
  ["P", "OP1", "OP2", "O5'", "O3'", "O4'"]
}

///|
/// Check whether each residue has a complete backbone.
///
/// Returns an array of (residue_index, is_complete, missing_atoms)
/// for each residue.
pub fn AtomArray::check_backbone_continuity(
  self : AtomArray,
) -> Array[(Int32, Bool, Array[String])] {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  let result : Array[(Int32, Bool, Array[String])] = []
  let mut i = 0
  while i < n_res {
    let s = starts[i]
    let e = if i + 1 < n_res { starts[i + 1] } else { self.len() }
    let rn = self.res_names[s]
    // Pick which backbone set to check
    let expected : Array[String] = if is_any_amino_acid(rn) {
      standard_backbone_atoms()
    } else {
      standard_phosphate_backbone_atoms()
    }
    // Collect present atom names
    let present : Array[String] = []
    let mut j = s
    while j < e {
      present.push(self.names[j])
      j = j + 1
    }
    // Check expected
    let missing : Array[String] = []
    for exp in expected {
      let mut found = false
      for p in present {
        if p == exp {
          found = true
          break
        }
      }
      if !found {
        missing.push(exp)
      }
    }
    result.push((i, missing.length() == 0, missing))
    i = i + 1
  }
  result
}

///|
/// Count residues with missing backbone atoms.
pub fn AtomArray::count_incomplete_backbone(self : AtomArray) -> Int32 {
  let checks = self.check_backbone_continuity()
  let mut n = 0
  for c in checks {
    let (_, ok, _) = c
    if !ok {
      n = n + 1
    }
  }
  n
}

// ============================================================
// P2-4: Filter to protein / nucleic acid only (exclude solvents/ions)
// ============================================================

///|
/// Remove solvent molecules and monoatomic ions from the AtomArray.
pub fn AtomArray::filter_biomolecule(self : AtomArray) -> AtomArray {
  let mask : Array[Bool] = []
  for i in 0..<self.len() {
    let rn = self.res_names[i]
    mask.push(!is_solvent(rn) && !is_monoatomic_ion(rn))
  }
  self.filter(mask)
}
