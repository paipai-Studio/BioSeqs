///|
/// Chemical Component Dictionary (CCD) — PDB chemical component lookup.
///
/// The CCD is a library of standard residue and ligand definitions used
/// by the PDB. Each entry describes the atoms, bonds, bond orders, and
/// stereochemistry of a chemical component.

///|
/// A chemical component entry from the CCD.
pub struct CcdEntry {
  id : String
  name : String
  comp_type : String
  formula : String
  atoms : Array[CcdAtom]
  bonds : Array[CcdBond]
  parent : String
  one_letter_code : String
}

///|
pub(all) struct CcdAtom {
  name : String
  element : String
  charge : Int
  leaving_atom : Bool
  aromatic : Bool
}

///|
pub(all) struct CcdBond {
  atom1 : String
  atom2 : String
  order : String
  aromatic : Bool
}

///|
/// Look up a standard CCD entry by its 3-letter code.
pub fn ccd_lookup(code : String) -> CcdEntry? {
  ccd_standard_table().get(code)
}

///|
/// Get the one-letter code for a 3-letter amino acid code.
pub fn ccd_one_letter(three_letter : String) -> String {
  match ccd_standard_table().get(three_letter) {
    Some(entry) => entry.one_letter_code
    None => "X"
  }
}

///|
/// Check if a 3-letter code is a standard amino acid.
pub fn ccd_is_amino_acid(code : String) -> Bool {
  match ccd_standard_table().get(code) {
    Some(entry) => entry.comp_type == "peptide linking" || entry.comp_type == "L-peptide linking" || entry.comp_type == "D-peptide linking"
    None => false
  }
}

///|
/// Check if a 3-letter code is a nucleotide.
pub fn ccd_is_nucleotide(code : String) -> Bool {
  match ccd_standard_table().get(code) {
    Some(entry) => entry.comp_type == "DNA linking" || entry.comp_type == "RNA linking"
    None => false
  }
}

///|
/// Check if a 3-letter code is a water molecule.
pub fn ccd_is_water(code : String) -> Bool {
  code == "HOH" || code == "WAT" || code == "DOD" || code == "H2O"
}

///|
/// Get the element for an atom in a residue.
pub fn ccd_atom_element(residue : String, atom_name : String) -> String {
  match ccd_standard_table().get(residue) {
    Some(entry) => {
      for i = 0; i < entry.atoms.length(); i = i + 1 {
        if entry.atoms[i].name == atom_name {
          return entry.atoms[i].element
        }
      }
      ccd_guess_element(atom_name)
    }
    None => ccd_guess_element(atom_name)
  }
}

///|
/// Guess the element from an atom name.
/// Checks 2-char metal names first, then falls back to 1-char elements.
fn ccd_guess_element(atom_name : String) -> String {
  if atom_name.length() == 0 {
    return "C"
  }
  // Check 2-char elements first (metals, halogens with 2-letter symbols)
  if atom_name.length() >= 2 {
    let two = atom_name.substring(start=0, end=2)
    let two_char = ["FE", "CL", "BR", "ZN", "MG", "MN", "CU", "NA", "NI", "CO",
      "SE", "MO", "RU", "OS", "IR", "PT", "AU", "HG", "PB", "SR",
      "BA", "AL", "SI", "LI", "BE", "NE", "AR", "KR", "XE", "CS",
      "RA", "FR", "AC", "AG", "CD", "GA", "GE", "AS", "IN", "SN",
      "SB", "TE", "BI", "PO", "AT", "RN", "RE", "RH", "PD", "TA",
      "TI", "SC", "ZR", "NB", "TC", "ND", "TB", "DY", "HO", "ER",
      "TM", "YB", "LU", "HF", "LA", "CE", "PR", "PM", "SM", "EU",
      "GD", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
      "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS",
      "RG", "CN", "FL", "LV", "HE"]
    for i = 0; i < two_char.length(); i = i + 1 {
      if two == two_char[i] {
        return two
      }
    }
  }
  // Fall back to 1-char elements
  let c0 = atom_name.unsafe_get(0).unsafe_to_char()
  if c0 == 'H' { return "H" }
  if c0 == 'C' { return "C" }
  if c0 == 'N' { return "N" }
  if c0 == 'O' { return "O" }
  if c0 == 'F' { return "F" }
  if c0 == 'P' { return "P" }
  if c0 == 'S' { return "S" }
  if c0 == 'K' { return "K" }
  if c0 == 'I' { return "I" }
  if c0 == 'B' { return "B" }
  if c0 == 'U' { return "U" }
  if c0 == 'V' { return "V" }
  if c0 == 'W' { return "W" }
  if c0 == 'Y' { return "Y" }
  "C"
}

///|
fn ccd_standard_table() -> Map[String, CcdEntry] {
  let m : Map[String, CcdEntry] = {}
  ccd_add_aa(m, "ALA", "ALANINE", "A", "C3 H7 N O2")
  ccd_add_aa(m, "ARG", "ARGININE", "R", "C6 H14 N4 O2")
  ccd_add_aa(m, "ASN", "ASPARAGINE", "N", "C4 H8 N2 O3")
  ccd_add_aa(m, "ASP", "ASPARTIC ACID", "D", "C4 H7 N O4")
  ccd_add_aa(m, "CYS", "CYSTEINE", "C", "C3 H7 N O2 S")
  ccd_add_aa(m, "GLN", "GLUTAMINE", "Q", "C5 H10 N2 O3")
  ccd_add_aa(m, "GLU", "GLUTAMIC ACID", "E", "C5 H9 N O4")
  ccd_add_aa(m, "GLY", "GLYCINE", "G", "C2 H5 N O2")
  ccd_add_aa(m, "HIS", "HISTIDINE", "H", "C6 H9 N3 O2")
  ccd_add_aa(m, "ILE", "ISOLEUCINE", "I", "C6 H13 N O2")
  ccd_add_aa(m, "LEU", "LEUCINE", "L", "C6 H13 N O2")
  ccd_add_aa(m, "LYS", "LYSINE", "K", "C6 H14 N2 O2")
  ccd_add_aa(m, "MET", "METHIONINE", "M", "C5 H11 N O2 S")
  ccd_add_aa(m, "PHE", "PHENYLALANINE", "F", "C9 H11 N O2")
  ccd_add_aa(m, "PRO", "PROLINE", "P", "C5 H9 N O2")
  ccd_add_aa(m, "SER", "SERINE", "S", "C3 H7 N O3")
  ccd_add_aa(m, "THR", "THREONINE", "T", "C4 H9 N O3")
  ccd_add_aa(m, "TRP", "TRYPTOPHAN", "W", "C11 H12 N2 O2")
  ccd_add_aa(m, "TYR", "TYROSINE", "Y", "C9 H11 N O3")
  ccd_add_aa(m, "VAL", "VALINE", "V", "C5 H11 N O2")
  ccd_add_special(m, "SEC", "SELENOCYSTEINE", "U", "C3 H7 N O2 Se", "L-peptide linking")
  ccd_add_special(m, "PYL", "PYRROLYSINE", "O", "C12 H21 N3 O3", "L-peptide linking")
  ccd_add_special(m, "ASX", "ASN/ASP", "B", "C4 H8 N2 O3", "peptide linking")
  ccd_add_special(m, "GLX", "GLN/GLU", "Z", "C5 H10 N2 O3", "peptide linking")
  ccd_add_special(m, "XAA", "UNKNOWN", "X", "C0 H0 N O", "peptide linking")
  ccd_add_nuc(m, "DA", "DEOXYADENOSINE", "A")
  ccd_add_nuc(m, "DC", "DEOXYCYTIDINE", "C")
  ccd_add_nuc(m, "DG", "DEOXYGUANOSINE", "G")
  ccd_add_nuc(m, "DT", "THYMIDINE", "T")
  ccd_add_nuc(m, "A", "ADENOSINE", "A")
  ccd_add_nuc(m, "C", "CYTIDINE", "C")
  ccd_add_nuc(m, "G", "GUANOSINE", "G")
  ccd_add_nuc(m, "U", "URIDINE", "U")
  ccd_add_special(m, "HOH", "WATER", "", "H2 O", "non-polymer")
  ccd_add_special(m, "WAT", "WATER", "", "H2 O", "non-polymer")
  ccd_add_special(m, "CA", "CALCIUM ION", "", "Ca", "non-polymer")
  ccd_add_special(m, "MG", "MAGNESIUM ION", "", "Mg", "non-polymer")
  ccd_add_special(m, "ZN", "ZINC ION", "", "Zn", "non-polymer")
  ccd_add_special(m, "NA", "SODIUM ION", "", "Na", "non-polymer")
  ccd_add_special(m, "K", "POTASSIUM ION", "", "K", "non-polymer")
  ccd_add_special(m, "FE", "IRON ION", "", "Fe", "non-polymer")
  ccd_add_special(m, "MN", "MANGANESE ION", "", "Mn", "non-polymer")
  ccd_add_special(m, "CU", "COPPER ION", "", "Cu", "non-polymer")
  ccd_add_special(m, "CL", "CHLORIDE ION", "", "Cl", "non-polymer")
  ccd_add_special(m, "HEM", "HEME", "", "C34 H32 Fe N4 O4", "non-polymer")
  ccd_add_special(m, "ATP", "ATP", "", "C10 H16 N5 O13 P3", "non-polymer")
  ccd_add_special(m, "GTP", "GTP", "", "C10 H16 N5 O14 P3", "non-polymer")
  ccd_add_special(m, "FAD", "FAD", "", "C27 H33 N9 O15 P2", "non-polymer")
  ccd_add_special(m, "NAD", "NAD", "", "C21 H27 N7 O14 P2", "non-polymer")
  ccd_add_special(m, "COA", "COENZYME A", "", "C21 H36 N7 O16 P3 S", "non-polymer")
  ccd_add_special(m, "PLP", "PYRIDOXAL PHOSPHATE", "", "C8 H10 N O6 P", "non-polymer")
  ccd_add_special(m, "RET", "RETINAL", "", "C20 H28 O", "non-polymer")
  m
}

///|
fn ccd_add_aa(
  m : Map[String, CcdEntry],
  id : String,
  name : String,
  one : String,
  formula : String,
) -> Unit {
  m[id] = CcdEntry::{
    id,
    name,
    comp_type: "L-peptide linking",
    formula,
    atoms: ccd_aa_atoms(id),
    bonds: [],
    parent: "",
    one_letter_code: one,
  }
}

///|
fn ccd_add_nuc(
  m : Map[String, CcdEntry],
  id : String,
  name : String,
  one : String,
) -> Unit {
  let is_dna = id.length() == 2 && id.unsafe_get(0).unsafe_to_char() == 'D'
  m[id] = CcdEntry::{
    id,
    name,
    comp_type: if is_dna { "DNA linking" } else { "RNA linking" },
    formula: "",
    atoms: [],
    bonds: [],
    parent: "",
    one_letter_code: one,
  }
}

///|
fn ccd_add_special(
  m : Map[String, CcdEntry],
  id : String,
  name : String,
  one : String,
  formula : String,
  comp_type : String,
) -> Unit {
  m[id] = CcdEntry::{
    id,
    name,
    comp_type,
    formula,
    atoms: [],
    bonds: [],
    parent: "",
    one_letter_code: one,
  }
}

///|
fn ccd_aa_atoms(id : String) -> Array[CcdAtom] {
  let atoms : Array[CcdAtom] = [
    CcdAtom::{ name: "N", element: "N", charge: 0, leaving_atom: true, aromatic: false },
    CcdAtom::{ name: "CA", element: "C", charge: 0, leaving_atom: false, aromatic: false },
    CcdAtom::{ name: "C", element: "C", charge: 0, leaving_atom: false, aromatic: false },
    CcdAtom::{ name: "O", element: "O", charge: 0, leaving_atom: true, aromatic: false },
    CcdAtom::{ name: "OXT", element: "O", charge: 0, leaving_atom: true, aromatic: false },
  ]
  match id {
    "ALA" => atoms.push(CcdAtom::{ name: "CB", element: "C", charge: 0, leaving_atom: false, aromatic: false })
    "CYS" => {
      atoms.push(CcdAtom::{ name: "CB", element: "C", charge: 0, leaving_atom: false, aromatic: false })
      atoms.push(CcdAtom::{ name: "SG", element: "S", charge: 0, leaving_atom: false, aromatic: false })
    }
    "SER" => {
      atoms.push(CcdAtom::{ name: "CB", element: "C", charge: 0, leaving_atom: false, aromatic: false })
      atoms.push(CcdAtom::{ name: "OG", element: "O", charge: 0, leaving_atom: false, aromatic: false })
    }
    "THR" => {
      atoms.push(CcdAtom::{ name: "CB", element: "C", charge: 0, leaving_atom: false, aromatic: false })
      atoms.push(CcdAtom::{ name: "OG1", element: "O", charge: 0, leaving_atom: false, aromatic: false })
      atoms.push(CcdAtom::{ name: "CG2", element: "C", charge: 0, leaving_atom: false, aromatic: false })
    }
    "GLY" => ()
    "PRO" => {
      atoms.push(CcdAtom::{ name: "CB", element: "C", charge: 0, leaving_atom: false, aromatic: false })
      atoms.push(CcdAtom::{ name: "CG", element: "C", charge: 0, leaving_atom: false, aromatic: false })
      atoms.push(CcdAtom::{ name: "CD", element: "C", charge: 0, leaving_atom: false, aromatic: false })
    }
    _ => ()
  }
  atoms
}
