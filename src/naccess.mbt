// Bio.PDB.NACCESS - Parser for NACCESS output files, mirroring BioPython's
// Bio.PDB.NACCESS.
//
// NACCESS (Atomic Solvent Accessible Area Calculations) is a standalone
// program that calculates the solvent accessible surface area of atoms,
// residues, and chains in a PDB structure. Biopython's Bio.PDB.NACCESS
// module runs the external NACCESS program and parses its output files.
//
// NACCESS produces three output files:
//   .asa  - Atomic solvent accessible surface areas
//   .rsa  - Residue solvent accessible surface areas (relative + absolute)
//   .log  - Log file
//
// The .rsa file format (per-residue):
//   RES  NUM  CHN  ABS  REL  ...
//   Where:
//     RES  - Residue name (3-letter code)
//     NUM  - Residue number
//     CHN  - Chain identifier
//     ABS  - Absolute accessible area (Å²)
//     REL  - Relative accessible area (% of standard)
//
// The .asa file format (per-atom):
//   ATM  NUM  RES  CHN  ...  ABS  ...

// ============================================================================
// Data structures
// ============================================================================

///|
/// Atomic solvent accessible area record.
pub struct NaccessAtom {
  pub atom_name : String // e.g. "CA", "N", "O"
  pub atom_num : Int
  pub res_name : String // 3-letter residue code
  pub res_num : Int
  pub chain_id : String
  pub abs_asa : Double // absolute accessible surface area (Å²)
  pub rel_asa : Double // relative accessible surface area (%)
  pub element : String
}

///|
/// Construct a NaccessAtom.
pub fn NaccessAtom::new(
  atom_name~ : String,
  atom_num~ : Int,
  res_name~ : String,
  res_num~ : Int,
  chain_id~ : String,
  abs_asa~ : Double,
  rel_asa~ : Double = 0.0,
  element~ : String = "",
) -> NaccessAtom {
  NaccessAtom::{
    atom_name,
    atom_num,
    res_name,
    res_num,
    chain_id,
    abs_asa,
    rel_asa,
    element,
  }
}

///|
/// Residue solvent accessible area record.
pub struct NaccessResidue {
  pub res_name : String // 3-letter code
  pub res_num : Int
  pub chain_id : String
  pub abs_all : Double // absolute ASA for all atoms
  pub abs_side_chain : Double // absolute ASA for side chain
  pub abs_main_chain : Double // absolute ASA for main chain
  pub rel_all : Double // relative ASA for all atoms (%)
  pub rel_side_chain : Double // relative ASA for side chain (%)
  pub rel_main_chain : Double // relative ASA for main chain (%)
  pub non_polar_abs : Double // non-polar absolute ASA
  pub non_polar_rel : Double // non-polar relative ASA
  pub polar_abs : Double // polar absolute ASA
  pub polar_rel : Double // polar relative ASA
}

///|
/// Construct a NaccessResidue.
pub fn NaccessResidue::new(
  res_name~ : String,
  res_num~ : Int,
  chain_id~ : String,
  abs_all~ : Double = 0.0,
  abs_side_chain~ : Double = 0.0,
  abs_main_chain~ : Double = 0.0,
  rel_all~ : Double = 0.0,
  rel_side_chain~ : Double = 0.0,
  rel_main_chain~ : Double = 0.0,
  non_polar_abs~ : Double = 0.0,
  non_polar_rel~ : Double = 0.0,
  polar_abs~ : Double = 0.0,
  polar_rel~ : Double = 0.0,
) -> NaccessResidue {
  NaccessResidue::{
    res_name,
    res_num,
    chain_id,
    abs_all,
    abs_side_chain,
    abs_main_chain,
    rel_all,
    rel_side_chain,
    rel_main_chain,
    non_polar_abs,
    non_polar_rel,
    polar_abs,
    polar_rel,
  }
}

///|
/// A parsed NACCESS result.
pub struct NaccessResult {
  pub structure_id : String
  residues : Array[NaccessResidue]
  atoms : Array[NaccessAtom]
  pub mut total_abs_asa : Double
  pub total_rel_asa : Double
  pub chain_totals : Map[String, Double]
}

///|
/// Construct an empty NaccessResult.
pub fn NaccessResult::new() -> NaccessResult {
  NaccessResult::{
    structure_id: "",
    residues: [],
    atoms: [],
    total_abs_asa: 0.0,
    total_rel_asa: 0.0,
    chain_totals: Map([], capacity=0),
  }
}

///|
/// Add a residue record.
pub fn NaccessResult::add_residue(self : NaccessResult, r : NaccessResidue) -> Unit {
  self.residues.push(r)
}

///|
/// Add an atom record.
pub fn NaccessResult::add_atom(self : NaccessResult, a : NaccessAtom) -> Unit {
  self.atoms.push(a)
}

///|
/// Get all residue records.
pub fn NaccessResult::get_residues(self : NaccessResult) -> Array[NaccessResidue] {
  self.residues
}

///|
/// Get all atom records.
pub fn NaccessResult::get_atoms(self : NaccessResult) -> Array[NaccessAtom] {
  self.atoms
}

///|
/// Get the number of residue records.
pub fn NaccessResult::get_num_residues(self : NaccessResult) -> Int {
  self.residues.length()
}

///|
/// Get the number of atom records.
pub fn NaccessResult::get_num_atoms(self : NaccessResult) -> Int {
  self.atoms.length()
}

// ============================================================================
// RSA file parsing (residue-level)
// ============================================================================

///|
/// Parse an NACCESS .rsa file content.
///
/// RSA file format (space-separated columns):
///   RES  NUM  CHN  ABS_ALL  REL_ALL  ABS_SIDE  REL_SIDE  ABS_MAIN  REL_MAIN
///   NON_POL_ABS  NON_POL_REL  POL_ABS  POL_REL
pub fn naccess_parse_rsa(content : String) -> Array[NaccessResidue] {
  let residues = Array::new()
  let lines = content.split("\n").to_array()
  for line in lines {
    let trimmed = line.to_string()
    if trimmed.length() < 10 {
      continue
    }
    // Lines starting with "RES" are residue records
    if trimmed.starts_with("RES ") {
      let parts = trimmed.split(" ").to_array()
      // Filter out empty strings from multiple spaces
      let filtered = Array::new()
      for p in parts {
        let s = p.to_string()
        if s.length() > 0 {
          filtered.push(s)
        }
      }
      if filtered.length() >= 5 {
        // Skip header lines: the residue-number field must be numeric.
        match parse_double(filtered[2]) {
          None => continue
          Some(_) => ()
        }
        let res_name = filtered[1]
        let res_num = parse_int(filtered[2])
        let chain_id = filtered[3]
        let abs_all = parse_double(filtered[4]).unwrap_or(0.0)
        let rel_all = if filtered.length() > 5 { parse_double(filtered[5]).unwrap_or(0.0) } else { 0.0 }
        let abs_side = if filtered.length() > 6 { parse_double(filtered[6]).unwrap_or(0.0) } else { 0.0 }
        let rel_side = if filtered.length() > 7 { parse_double(filtered[7]).unwrap_or(0.0) } else { 0.0 }
        let abs_main = if filtered.length() > 8 { parse_double(filtered[8]).unwrap_or(0.0) } else { 0.0 }
        let rel_main = if filtered.length() > 9 { parse_double(filtered[9]).unwrap_or(0.0) } else { 0.0 }
        let np_abs = if filtered.length() > 10 { parse_double(filtered[10]).unwrap_or(0.0) } else { 0.0 }
        let np_rel = if filtered.length() > 11 { parse_double(filtered[11]).unwrap_or(0.0) } else { 0.0 }
        let p_abs = if filtered.length() > 12 { parse_double(filtered[12]).unwrap_or(0.0) } else { 0.0 }
        let p_rel = if filtered.length() > 13 { parse_double(filtered[13]).unwrap_or(0.0) } else { 0.0 }
        residues.push(
          NaccessResidue::new(
            res_name=res_name,
            res_num=res_num,
            chain_id=chain_id,
            abs_all=abs_all,
            rel_all=rel_all,
            abs_side_chain=abs_side,
            rel_side_chain=rel_side,
            abs_main_chain=abs_main,
            rel_main_chain=rel_main,
            non_polar_abs=np_abs,
            non_polar_rel=np_rel,
            polar_abs=p_abs,
            polar_rel=p_rel,
          ),
        )
      }
    }
  }
  residues
}

// ============================================================================
// ASA file parsing (atom-level)
// ============================================================================

///|
/// Parse an NACCESS .asa file content.
///
/// ASA file format (space-separated columns):
///   ATM  NUM  RES  CHN  ATOM_NAME  ...  ABS_ASA  REL_ASA
pub fn naccess_parse_asa(content : String) -> Array[NaccessAtom] {
  let atoms = Array::new()
  let lines = content.split("\n").to_array()
  for line in lines {
    let trimmed = line.to_string()
    if trimmed.length() < 10 {
      continue
    }
    // Lines starting with "ATM" are atom records
    if trimmed.starts_with("ATM ") || trimmed.starts_with("HETATM") {
      let parts = trimmed.split(" ").to_array()
      let filtered = Array::new()
      for p in parts {
        let s = p.to_string()
        if s.length() > 0 {
          filtered.push(s)
        }
      }
      if filtered.length() >= 5 {
        // Skip header lines: the atom-number field must be numeric.
        match parse_double(filtered[1]) {
          None => continue
          Some(_) => ()
        }
        let is_het = filtered[0] == "HETATM"
        let prefix_len = if is_het { 0 } else { 0 }
        let atom_name = if filtered.length() > 4 { filtered[4] } else { "" }
        let atom_num = parse_int(filtered[1 + prefix_len])
        let res_name = filtered[2 + prefix_len]
        let res_num = if is_het {
          parse_int(filtered[3 + prefix_len])
        } else {
          parse_int(filtered[3])
        }
        let chain_id = filtered[3 + prefix_len]
        let abs_asa = if filtered.length() > 5 + prefix_len {
          parse_double(filtered[filtered.length() - 2]).unwrap_or(0.0)
        } else {
          0.0
        }
        let rel_asa = if filtered.length() > 6 + prefix_len {
          parse_double(filtered[filtered.length() - 1]).unwrap_or(0.0)
        } else {
          0.0
        }
        atoms.push(
          NaccessAtom::new(
            atom_name=atom_name,
            atom_num=atom_num,
            res_name=res_name,
            res_num=res_num,
            chain_id=chain_id,
            abs_asa=abs_asa,
            rel_asa=rel_asa,
          ),
        )
      }
    }
  }
  atoms
}

// ============================================================================
// Combined parsing
// ============================================================================

///|
/// Parse NACCESS output from both RSA and ASA content.
pub fn naccess_parse(rsa_content : String, asa_content : String) -> NaccessResult {
  let result = NaccessResult::new()
  let residues = naccess_parse_rsa(rsa_content)
  let atoms = naccess_parse_asa(asa_content)
  for r in residues {
    result.add_residue(r)
  }
  for a in atoms {
    result.add_atom(a)
  }
  // Compute totals
  let mut total_abs = 0.0
  for r in residues {
    total_abs = total_abs + r.abs_all
    match result.chain_totals.get(r.chain_id) {
      Some(v) => result.chain_totals[r.chain_id] = v + r.abs_all
      None => result.chain_totals[r.chain_id] = r.abs_all
    }
  }
  result.total_abs_asa = total_abs
  result
}

// ============================================================================
// Serialization
// ============================================================================

///|
/// Serialize a NaccessResidue to RSA format line.
pub fn naccess_residue_to_rsa(r : NaccessResidue) -> String {
  "RES  " +
  r.res_name +
  "  " +
  r.res_num.to_string() +
  "  " +
  r.chain_id +
  "  " +
  r.abs_all.to_string() +
  "  " +
  r.rel_all.to_string() +
  "  " +
  r.abs_side_chain.to_string() +
  "  " +
  r.rel_side_chain.to_string() +
  "  " +
  r.abs_main_chain.to_string() +
  "  " +
  r.rel_main_chain.to_string() +
  "  " +
  r.non_polar_abs.to_string() +
  "  " +
  r.non_polar_rel.to_string() +
  "  " +
  r.polar_abs.to_string() +
  "  " +
  r.polar_rel.to_string()
}

///|
/// Serialize a NaccessAtom to ASA format line.
pub fn naccess_atom_to_asa(a : NaccessAtom) -> String {
  "ATM  " +
  a.atom_num.to_string() +
  "  " +
  a.res_name +
  "  " +
  a.chain_id +
  "  " +
  a.atom_name +
  "  " +
  a.abs_asa.to_string() +
  "  " +
  a.rel_asa.to_string()
}

///|
/// Serialize a NaccessResult to RSA file content.
pub fn naccess_write_rsa(result : NaccessResult) -> String {
  let sb = StringBuilder::new()
  sb.write_string("NACCESS  Version 2.1\n")
  sb.write_string("RES     NUM  CHN  ABS_ALL  REL_ALL  ABS_SIDE  REL_SIDE  ABS_MAIN  REL_MAIN  NP_ABS  NP_REL  P_ABS  P_REL\n")
  for r in result.residues {
    sb.write_string(naccess_residue_to_rsa(r) + "\n")
  }
  sb.to_string()
}

///|
/// Serialize a NaccessResult to ASA file content.
pub fn naccess_write_asa(result : NaccessResult) -> String {
  let sb = StringBuilder::new()
  sb.write_string("NACCESS  Version 2.1\n")
  sb.write_string("ATM     NUM  RES  CHN  ATOM_NAME  ABS_ASA  REL_ASA\n")
  for a in result.atoms {
    sb.write_string(naccess_atom_to_asa(a) + "\n")
  }
  sb.to_string()
}

// ============================================================================
// Analysis functions
// ============================================================================

///|
/// Get residues with relative ASA above a threshold (buried vs exposed).
pub fn naccess_exposed_residues(
  result : NaccessResult,
  threshold : Double,
) -> Array[NaccessResidue] {
  let exposed = Array::new()
  for r in result.residues {
    if r.rel_all > threshold {
      exposed.push(r)
    }
  }
  exposed
}

///|
/// Get buried residues (relative ASA below threshold).
pub fn naccess_buried_residues(
  result : NaccessResult,
  threshold : Double,
) -> Array[NaccessResidue] {
  let buried = Array::new()
  for r in result.residues {
    if r.rel_all < threshold {
      buried.push(r)
    }
  }
  buried
}

///|
/// Get the average relative ASA for all residues.
pub fn naccess_mean_rel_asa(result : NaccessResult) -> Double {
  if result.residues.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for r in result.residues {
    sum = sum + r.rel_all
  }
  sum / result.residues.length().to_double()
}

///|
/// Get residues for a specific chain.
pub fn naccess_get_chain_residues(
  result : NaccessResult,
  chain_id : String,
) -> Array[NaccessResidue] {
  let chain_res = Array::new()
  for r in result.residues {
    if r.chain_id == chain_id {
      chain_res.push(r)
    }
  }
  chain_res
}

///|
/// Count residues by type (exposed vs buried).
pub fn naccess_count_exposed(
  result : NaccessResult,
  threshold : Double,
) -> Int {
  let mut count = 0
  for r in result.residues {
    if r.rel_all > threshold {
      count = count + 1
    }
  }
  count
}

///|
/// Get total ASA for a specific chain.
pub fn naccess_chain_total_asa(result : NaccessResult, chain_id : String) -> Double {
  let mut total = 0.0
  for r in result.residues {
    if r.chain_id == chain_id {
      total = total + r.abs_all
    }
  }
  total
}

///|
/// Classify residue as exposed (True) or buried (False).
pub fn naccess_is_exposed(r : NaccessResidue, threshold : Double) -> Bool {
  r.rel_all > threshold
}

// ============================================================================
// Standard ASA reference values (used for relative ASA calculation)
// ============================================================================

///|
/// Standard ASA values for amino acid residues (maximum accessible area).
/// Used by NACCESS to compute relative ASA = (observed / standard) * 100.
pub fn naccess_standard_asa(res_name : String) -> Double {
  match res_name {
    "ALA" => 113.0
    "ARG" => 241.0
    "ASN" => 158.0
    "ASP" => 151.0
    "CYS" => 140.0
    "GLN" => 189.0
    "GLU" => 183.0
    "GLY" => 85.0
    "HIS" => 194.0
    "ILE" => 182.0
    "LEU" => 180.0
    "LYS" => 211.0
    "MET" => 204.0
    "PHE" => 218.0
    "PRO" => 143.0
    "SER" => 122.0
    "THR" => 146.0
    "TRP" => 259.0
    "TYR" => 229.0
    "VAL" => 160.0
    _ => 0.0
  }
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Create a sample NACCESS RSA content for testing.
pub fn naccess_sample_rsa() -> String {
  "NACCESS  Version 2.1\nRES     NUM  CHN  ABS_ALL  REL_ALL  ABS_SIDE  REL_SIDE  ABS_MAIN  REL_MAIN  NP_ABS  NP_REL  P_ABS  P_REL\nRES MET 1 A 45.3 22.1 30.1 18.5 15.2 35.2 25.1 20.1 20.2 25.1\nRES ALA 2 A 20.1 17.8 5.1 8.2 15.0 34.9 10.0 8.9 10.1 14.2\nRES CYS 3 A 5.2 3.7 1.1 1.0 4.1 9.5 3.0 2.1 2.2 3.1\nRES HOH 1 B 25.0 0.0 25.0 0.0 0.0 0.0 25.0 0.0 0.0 0.0\n"
}

///|
/// Create a sample NACCESS ASA content for testing.
pub fn naccess_sample_asa() -> String {
  "NACCESS  Version 2.1\nATM     NUM  RES  CHN  ATOM_NAME  ABS_ASA  REL_ASA\nATM     1   MET  A    N        15.2    35.2\nATM     2   MET  A    CA       10.1    22.1\nATM     3   ALA  A    N         8.0    18.6\nATM     4   ALA  A    CA        5.1     8.2\nATM     5   CYS  A    N         2.1     4.9\nATM     6   HOH  B    O        25.0     0.0\n"
}

///|
/// Create a sample NaccessResult for testing.
pub fn naccess_sample() -> NaccessResult {
  naccess_parse(naccess_sample_rsa(), naccess_sample_asa())
}

///|
/// Generate a summary string.
pub fn naccess_summary(result : NaccessResult) -> String {
  let sb = StringBuilder::new()
  sb.write_string("NACCESS Result Summary:\n")
  sb.write_string("  Residues: " + result.get_num_residues().to_string() + "\n")
  sb.write_string("  Atoms: " + result.get_num_atoms().to_string() + "\n")
  sb.write_string("  Total ASA: " + result.total_abs_asa.to_string() + " Å²\n")
  sb.write_string("  Mean Rel ASA: " + naccess_mean_rel_asa(result).to_string() + "%\n")
  sb.write_string("  Exposed (>20%): " + naccess_count_exposed(result, 20.0).to_string() + "\n")
  sb.write_string("  Buried (≤20%): " + (result.get_num_residues() - naccess_count_exposed(result, 20.0)).to_string() + "\n")
  sb.write_string("  Chains: " + result.chain_totals.size().to_string() + "\n")
  for chain_id in result.chain_totals.keys() {
    match result.chain_totals.get(chain_id) {
      Some(v) =>
        sb.write_string("    Chain " + chain_id + ": " + v.to_string() + " Å²\n")
      None => ()
    }
  }
  sb.to_string()
}
