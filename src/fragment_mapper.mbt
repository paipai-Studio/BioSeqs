///|
/// Bio.PDB.FragmentMapper module, inspired by Biopython Bio.PDB.FragmentMapper.
/// Provides protein structure fragment mapping and classification.

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

///|
/// Fragment classification types based on DSSP assignments.
pub enum FragmentType {
  Helix
  Sheet
  Loop
  Coil
  Bridge
  Turn
} derive(Eq, Debug)

///| Construct a Helix FragmentType.
pub fn ft_helix() -> FragmentType { Helix }
///| Construct a Sheet FragmentType.
pub fn ft_sheet() -> FragmentType { Sheet }
///| Construct a Loop FragmentType.
pub fn ft_loop() -> FragmentType { Loop }
///| Construct a Coil FragmentType.
pub fn ft_coil() -> FragmentType { Coil }
///| Construct a Bridge FragmentType.
pub fn ft_bridge() -> FragmentType { Bridge }
///| Construct a Turn FragmentType.
pub fn ft_turn() -> FragmentType { Turn }

///|
/// A single residue within a fragment.
pub struct FragmentResidue {
  pub res_name : String
  pub res_num : Int
  pub atom_count : Int
  pub is_c_alpha : Bool
  pub is_sidechain : Bool
} derive(Eq, Debug)

///|
/// A protein structure fragment.
pub struct Fragment {
  pub id : String
  pub start_residue : Int
  pub end_residue : Int
  pub chain_id : String
  pub fragment_type : FragmentType
  pub residues : Array[FragmentResidue]
} derive(Eq, Debug)

///|
/// Result of fragment mapping operation.
pub struct FragmentMapperResult {
  pub fragments : Array[Fragment]
  pub protein_sequence : String
  pub sequence_length : Int
  pub coverage : Map[Int, Double]
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

///|
/// Map 3-letter amino acid code to 1-letter code.
fn _three_to_one(res_name : String) -> String {
  match res_name {
    "ALA" => "A"
    "ARG" => "R"
    "ASN" => "N"
    "ASP" => "D"
    "CYS" => "C"
    "GLN" => "Q"
    "GLU" => "E"
    "GLY" => "G"
    "HIS" => "H"
    "ILE" => "I"
    "LEU" => "L"
    "LYS" => "K"
    "MET" => "M"
    "PHE" => "F"
    "PRO" => "P"
    "SER" => "S"
    "THR" => "T"
    "TRP" => "W"
    "TYR" => "Y"
    "VAL" => "V"
    _ => "X"
  }
}

///|
/// Trim whitespace from a string (spaces only).
fn trim_spaces(s : String) -> String {
  let mut start = 0
  while start < s.length() && s.unsafe_get(start).to_int() == 32 {
    start = start + 1
  }
  let mut end = s.length()
  while end > start && s.unsafe_get(end - 1).to_int() == 32 {
    end = end - 1
  }
  if start >= end {
    ""
  } else {
    s[start:end].to_owned()
  }
}

///|
/// Classify a single DSSP character to FragmentType.
fn classify_char(c : Int) -> FragmentType {
  if c == 72 || c == 71 || c == 73 {
    return Helix
  }
  if c == 69 {
    return Sheet
  }
  if c == 76 {
    return Loop
  }
  if c == 66 || c == 98 {
    return Bridge
  }
  if c == 84 || c == 83 {
    return Turn
  }
  Coil
}

// ---------------------------------------------------------------------------
// FragmentType helpers
// ---------------------------------------------------------------------------

///|
/// Convert a FragmentType to its string representation.
pub fn FragmentType::to_string(self : FragmentType) -> String {
  match self {
    Helix => "Helix"
    Sheet => "Sheet"
    Loop => "Loop"
    Coil => "Coil"
    Bridge => "Bridge"
    Turn => "Turn"
  }
}

///|
/// Convert a FragmentType to a single-letter DSSP-like code.
pub fn FragmentType::to_char(self : FragmentType) -> String {
  match self {
    Helix => "H"
    Sheet => "E"
    Loop => "L"
    Coil => "C"
    Bridge => "b"
    Turn => "T"
  }
}

// ---------------------------------------------------------------------------
// FragmentResidue helpers
// ---------------------------------------------------------------------------

///|
/// Create a new FragmentResidue.
pub fn FragmentResidue::new(
  res_name : String,
  res_num : Int,
  atom_count? : Int = 0,
  is_c_alpha? : Bool = false,
  is_sidechain? : Bool = false,
) -> FragmentResidue {
  FragmentResidue::{ res_name, res_num, atom_count, is_c_alpha, is_sidechain }
}

// ---------------------------------------------------------------------------
// Fragment helpers
// ---------------------------------------------------------------------------

///|
/// Create a new Fragment.
pub fn Fragment::new(
  id : String,
  start_residue : Int,
  end_residue : Int,
  chain_id? : String = "A",
  fragment_type? : FragmentType = Coil,
  residues? : Array[FragmentResidue] = [],
) -> Fragment {
  Fragment::{
    id,
    start_residue,
    end_residue,
    chain_id,
    fragment_type,
    residues,
  }
}

///|
/// Get the length (number of residues) of a fragment.
pub fn Fragment::length(self : Fragment) -> Int {
  self.end_residue - self.start_residue + 1
}

// ---------------------------------------------------------------------------
// FragmentMapperResult helpers
// ---------------------------------------------------------------------------

///|
/// Get the number of fragments in the result.
pub fn FragmentMapperResult::num_fragments(self : FragmentMapperResult) -> Int {
  self.fragments.length()
}

// ---------------------------------------------------------------------------
// Core API functions (fm_ prefix)
// ---------------------------------------------------------------------------

///|
/// Classify a DSSP secondary structure string into an array of FragmentType.
///
/// DSSP assignment mapping:
///   - 'H', 'G', 'I' → Helix
///   - 'E', 'B' → Sheet
///   - 'L', 'C', ' ' → Loop
///   - 'b'           → Bridge
///   - 'T', 'S'      → Turn
///   - other         → Coil
pub fn fm_classify_secondary_structure(
  ss_string : String,
) -> Array[FragmentType] {
  let result : Array[FragmentType] = Array::new()
  let mut i = 0
  while i < ss_string.length() {
    let c = ss_string.unsafe_get(i).to_int()
    result.push(classify_char(c))
    i = i + 1
  }
  result
}

///|
/// Assign fragments from an amino acid sequence and a DSSP secondary
/// structure string. Consecutive positions of the same fragment type
/// are merged into a single fragment.
pub fn fm_assign_fragments(
  sequence : String,
  secondary_structure : String,
) -> FragmentMapperResult {
  let ss_types = fm_classify_secondary_structure(secondary_structure)
  let seq_len = sequence.length()
  let frag_len = ss_types.length()
  let total = if seq_len > frag_len { seq_len } else { frag_len }

  let fragments : Array[Fragment] = Array::new()
  let coverage : Map[Int, Double] = Map([], capacity=total)

  let mut i = 0
  let mut fragment_idx = 0
  while i < frag_len {
    let current_type = ss_types[i]
    let start = i
    let mut end = i
    while end + 1 < frag_len && ss_types[end + 1] == current_type {
      end = end + 1
    }

    let residues : Array[FragmentResidue] = Array::new()
    let mut j = start
    while j <= end {
      let res_name = if j < seq_len {
        char_at(sequence, j)
      } else {
        "X"
      }
      residues.push(FragmentResidue::{
        res_name,
        res_num: j + 1,
        atom_count: 0,
        is_c_alpha: false,
        is_sidechain: false,
      })
      j = j + 1
    }

    let id = "fragment_" + fragment_idx.to_string()
    fragments.push(Fragment::{
      id,
      start_residue: start + 1,
      end_residue: end + 1,
      chain_id: "A",
      fragment_type: current_type,
      residues,
    })

    let mut k = start
    while k <= end {
      coverage[k] = 1.0
      k = k + 1
    }

    i = end + 1
    fragment_idx = fragment_idx + 1
  }

  let protein_seq = if seq_len > 0 { sequence } else { "" }

  FragmentMapperResult::{
    fragments,
    protein_sequence: protein_seq,
    sequence_length: seq_len,
    coverage,
  }
}

///|
/// Parse PDB fragment assignment from content string containing ATOM/HETATM
/// records. Groups residues by chain and creates fragments covering each chain.
pub fn fm_parse_fragments(content : String) -> FragmentMapperResult {
  let lines = split_by_char(content, 10)

  let residue_map : Map[String, FragmentResidue] = Map([], capacity=256)
  let chain_map : Map[String, Array[FragmentResidue]] = Map([], capacity=16)
  let mut sequence = ""

  for line in lines {
    if line.length() < 54 {
      continue
    }
    let record_type = substring(line, 0, 6)
    if record_type == "ATOM  " || record_type == "HETATM" {
      let atom_name = trim_spaces(substring(line, 12, 4))
      let res_name = trim_spaces(substring(line, 17, 3))
      let chain_id = substring(line, 21, 1)
      let res_num_str = trim_spaces(substring(line, 22, 4))

      if res_num_str.length() == 0 {
        continue
      }
      let res_num = parse_int(res_num_str)
      let key = chain_id + res_num.to_string()

      let existing = residue_map.get(key)
      match existing {
        Some(prev) => {
          let updated = FragmentResidue::{
            res_name: prev.res_name,
            res_num: prev.res_num,
            atom_count: prev.atom_count + 1,
            is_c_alpha: prev.is_c_alpha || atom_name == "CA",
            is_sidechain: prev.is_sidechain ||
              !(atom_name == "N" || atom_name == "CA" || atom_name == "C"),
          }
          residue_map[key] = updated
          let chain_res = chain_map.get(chain_id)
          match chain_res {
            Some(arr) => {
              let new_arr = arr.copy()
              let mut idx = 0
              while idx < new_arr.length() {
                if new_arr[idx].res_num == res_num {
                  new_arr[idx] = updated
                  break
                }
                idx = idx + 1
              }
              chain_map[chain_id] = new_arr
            }
            None => ()
          }
        }
        None => {
          let is_ca = atom_name == "CA"
          let is_sc = !(atom_name == "N" || atom_name == "CA" || atom_name == "C")
          let new_res = FragmentResidue::{
            res_name,
            res_num,
            atom_count: 1,
            is_c_alpha: is_ca,
            is_sidechain: is_sc,
          }
          residue_map[key] = new_res

          let chain_res = chain_map.get(chain_id)
          match chain_res {
            Some(arr) => {
              let new_arr = arr.copy()
              new_arr.push(new_res)
              chain_map[chain_id] = new_arr
            }
            None => {
              chain_map[chain_id] = [new_res]
            }
          }

          sequence = sequence + _three_to_one(res_name)
        }
      }
    }
  }

  let fragments : Array[Fragment] = Array::new()
  let coverage : Map[Int, Double] = Map([], capacity=sequence.length())

  let mut frag_idx = 0
  for chain_id, chain_residues in chain_map {
    if chain_residues.length() == 0 {
      continue
    }

    let sorted_residues = _sort_by_resnum(chain_residues)
    let start_res = sorted_residues[0].res_num
    let end_res = sorted_residues[sorted_residues.length() - 1].res_num

    let id = "fragment_" + frag_idx.to_string()
    fragments.push(Fragment::{
      id,
      start_residue: start_res,
      end_residue: end_res,
      chain_id,
      fragment_type: Coil,
      residues: sorted_residues,
    })

    let mut pos = 0
    for _ in sorted_residues {
      coverage[pos] = 1.0
      pos = pos + 1
    }

    frag_idx = frag_idx + 1
  }

  FragmentMapperResult::{
    fragments,
    protein_sequence: sequence,
    sequence_length: sequence.length(),
    coverage,
  }
}

///|
/// Sort an array of FragmentResidue by res_num (ascending).
fn _sort_by_resnum(residues : Array[FragmentResidue]) -> Array[FragmentResidue] {
  let n = residues.length()
  let result = residues.copy()
  let mut i = 1
  while i < n {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j].res_num > key.res_num {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
    i = i + 1
  }
  result
}

///|
/// Get a fragment by its ID string. Returns None if not found.
pub fn fm_get_fragment(
  result : FragmentMapperResult,
  id : String,
) -> Fragment? {
  for frag in result.fragments {
    if frag.id == id {
      return Some(frag)
    }
  }
  None
}

///|
/// Classify the fragment type based on its residues and fragment metadata.
/// Returns the fragment's fragment_type field.
pub fn fm_fragment_type(fragment : Fragment) -> FragmentType {
  fragment.fragment_type
}

///|
/// Get the amino acid sequence (1-letter codes) of a fragment.
pub fn fm_fragment_sequence(fragment : Fragment) -> String {
  let mut result = ""
  for res in fragment.residues {
    result = result + res.res_name
  }
  result
}

///|
/// Get per-residue coverage as a map from residue index to coverage value.
pub fn fm_coverage(result : FragmentMapperResult) -> Map[Int, Double] {
  result.coverage
}

///|
/// Get the start and end positions (1-based) of a fragment as a tuple.
pub fn fm_fragment_boundaries(fragment : Fragment) -> (Int, Int) {
  (fragment.start_residue, fragment.end_residue)
}

///|
/// Merge adjacent fragments of the same type on the same chain into
/// larger fragments. Non-adjacent or differently-typed fragments are
/// left unchanged.
pub fn fm_merge_fragments(fragments : Array[Fragment]) -> Array[Fragment] {
  if fragments.length() == 0 {
    return []
  }

  let sorted = _sort_fragments_by_position(fragments)
  let result : Array[Fragment] = Array::new()
  let mut current = sorted[0]

  let mut i = 1
  while i < sorted.length() {
    let next = sorted[i]
    let is_adjacent = next.start_residue == current.end_residue + 1
    let same_type = next.fragment_type == current.fragment_type
    let same_chain = next.chain_id == current.chain_id

    if is_adjacent && same_type && same_chain {
      let merged_residues = _concat_arrays(current.residues, next.residues)
      current = Fragment::{
        id: current.id,
        start_residue: current.start_residue,
        end_residue: next.end_residue,
        chain_id: current.chain_id,
        fragment_type: current.fragment_type,
        residues: merged_residues,
      }
    } else {
      result.push(current)
      current = next
    }
    i = i + 1
  }
  result.push(current)
  result
}

///|
/// Sort fragments by start_residue (ascending), then by chain_id.
fn _sort_fragments_by_position(fragments : Array[Fragment]) -> Array[Fragment] {
  let n = fragments.length()
  let result = fragments.copy()
  let mut i = 1
  while i < n {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 {
      let should_swap =
        result[j].start_residue > key.start_residue ||
        (result[j].start_residue == key.start_residue &&
          result[j].chain_id > key.chain_id)
      if should_swap {
        result[j + 1] = result[j]
        j = j - 1
      } else {
        break
      }
    }
    result[j + 1] = key
    i = i + 1
  }
  result
}

///|
/// Concatenate two arrays of FragmentResidue.
fn _concat_arrays(
  a : Array[FragmentResidue],
  b : Array[FragmentResidue],
) -> Array[FragmentResidue] {
  let result : Array[FragmentResidue] = Array::new()
  for item in a {
    result.push(item)
  }
  for item in b {
    result.push(item)
  }
  result
}

///|
/// Filter fragments by minimum length (number of residues).
/// Fragments shorter than min_length are excluded.
pub fn fm_filter_by_length(
  fragments : Array[Fragment],
  min_length : Int,
) -> Array[Fragment] {
  let result : Array[Fragment] = Array::new()
  for frag in fragments {
    if frag.length() >= min_length {
      result.push(frag)
    }
  }
  result
}