///|
/// bp3_iupac_pdbdata.mbt — Biopython round 3: IUPACData + vdw_radii + PDBData
/// #11 IUPACData — ambiguity codes, monoisotopic/average masses
/// #18 vdw_radii — comprehensive van der Waals radii table
/// #24 PDBData — PDB data tables (residue masses, bond lengths/angles)

///|
/// #11 IUPACData
let bp3_iupac_ambiguous_dna_values : Map[String, String] = Map([
  ("A", "A"),
  ("C", "C"),
  ("G", "G"),
  ("T", "T"),
  ("U", "T"),
  ("R", "AG"),
  ("Y", "CT"),
  ("S", "GC"),
  ("W", "AT"),
  ("K", "GT"),
  ("M", "AC"),
  ("B", "CGT"),
  ("D", "AGT"),
  ("H", "ACT"),
  ("V", "ACG"),
  ("N", "ACGT"),
])

///|
let bp3_iupac_ambiguous_rna_values : Map[String, String] = Map([
  ("A", "A"),
  ("C", "C"),
  ("G", "G"),
  ("U", "U"),
  ("T", "U"),
  ("R", "AG"),
  ("Y", "CU"),
  ("S", "GC"),
  ("W", "AU"),
  ("K", "GU"),
  ("M", "AC"),
  ("B", "CGU"),
  ("D", "AGU"),
  ("H", "ACU"),
  ("V", "ACG"),
  ("N", "ACGU"),
])

///|
let bp3_iupac_extended_dna_values : Map[String, String] = Map([
  ("A", "A"),
  ("C", "C"),
  ("G", "G"),
  ("T", "T"),
  ("B", "CGT"),
  ("D", "AGT"),
  ("H", "ACT"),
  ("V", "ACG"),
  ("N", "ACGT"),
])

///|
let bp3_iupac_ambiguous_dna_letters : String = "ACGTURYSWKMBDHVN"

///|
let bp3_iupac_ambiguous_rna_letters : String = "ACGUURYSWKMBDHVN"

///|
let bp3_iupac_protein_letters : String = "ACDEFGHIKLMNPQRSTVWY"

///|
let bp3_iupac_extended_protein_letters : String = "ACDEFGHIKLMNPQRSTVWYBXZJUO*"

///|
pub fn bp3_iupac_expand(ambiguous : String, is_rna~ : Bool = false) -> String {
  let table = if is_rna {
    bp3_iupac_ambiguous_rna_values
  } else {
    bp3_iupac_ambiguous_dna_values
  }
  let r = table.get(ambiguous.to_upper())
  match r {
    Some(v) => v
    None => ""
  }
}

///|
pub fn bp3_iupac_valid_letter(letter : String, is_rna~ : Bool = false) -> Bool {
  let letters = if is_rna {
    bp3_iupac_ambiguous_rna_letters
  } else {
    bp3_iupac_ambiguous_dna_letters
  }
  for i in 0..<letters.length() {
    if letters[i:i + 1].to_owned() == letter.to_upper() {
      return true
    }
  }
  false
}

///|
/// #11 IUPACData — monoisotopic/average masses delegated to data.mbt
/// protein_residue_mass_monoisotopic() and protein_residue_mass_average().

///|
pub fn bp3_molecular_weight_monoisotopic(sequence : String) -> Double {
  protein_molecular_weight(sequence, is_monoisotopic=true)
}

///|
pub fn bp3_molecular_weight_average(sequence : String) -> Double {
  protein_molecular_weight(sequence, is_monoisotopic=false)
}

///|
/// #18 vdw_radii — comprehensive van der Waals radii table (Bondi 1964 + extensions)
let bp3_vdw_radii : Map[String, Double] = Map([
  ("H", 1.20),
  ("HE", 1.40),
  ("LI", 1.82),
  ("BE", 1.53),
  ("B", 1.92),
  ("C", 1.70),
  ("N", 1.55),
  ("O", 1.52),
  ("F", 1.47),
  ("NE", 1.54),
  ("NA", 2.27),
  ("MG", 1.73),
  ("AL", 1.84),
  ("SI", 2.10),
  ("P", 1.80),
  ("S", 1.80),
  ("CL", 1.75),
  ("AR", 1.88),
  ("K", 2.75),
  ("CA", 2.31),
  ("SC", 2.11),
  ("TI", 2.08),
  ("V", 1.92),
  ("CR", 1.83),
  ("MN", 1.79),
  ("FE", 1.83),
  ("CO", 1.83),
  ("NI", 1.63),
  ("CU", 1.40),
  ("ZN", 1.39),
  ("GA", 1.87),
  ("GE", 2.11),
  ("AS", 1.85),
  ("SE", 1.90),
  ("BR", 1.85),
  ("KR", 2.02),
  ("RB", 2.65),
  ("SR", 2.49),
  ("Y", 2.22),
  ("ZR", 1.72),
  ("NB", 1.52),
  ("MO", 1.72),
  ("TC", 1.72),
  ("RU", 1.72),
  ("RH", 1.72),
  ("PD", 1.63),
  ("AG", 1.72),
  ("CD", 1.58),
  ("IN", 1.93),
  ("SN", 2.17),
  ("SB", 2.06),
  ("TE", 2.06),
  ("I", 1.98),
  ("XE", 2.16),
  ("CS", 3.43),
  ("BA", 2.68),
  ("LA", 2.22),
  ("CE", 2.22),
  ("PR", 2.22),
  ("ND", 2.22),
  ("PM", 2.22),
  ("SM", 2.22),
  ("EU", 2.22),
  ("GD", 2.22),
  ("TB", 2.22),
  ("DY", 2.22),
  ("HO", 2.22),
  ("ER", 2.22),
  ("TM", 2.22),
  ("YB", 2.22),
  ("LU", 2.22),
  ("HF", 1.72),
  ("TA", 1.72),
  ("W", 1.72),
  ("RE", 1.72),
  ("OS", 1.72),
  ("IR", 1.72),
  ("PT", 1.72),
  ("AU", 1.66),
  ("HG", 1.55),
  ("TL", 1.96),
  ("PB", 2.02),
  ("BI", 2.07),
  ("PO", 1.97),
  ("AT", 2.02),
  ("RN", 2.20),
  ("FR", 3.48),
  ("RA", 2.83),
  ("AC", 2.22),
  ("TH", 2.22),
  ("PA", 2.22),
  ("U", 1.86),
  ("NP", 1.86),
  ("PU", 1.86),
])

///|
pub fn bp3_vdw_radius(element : String) -> Double {
  bp3_vdw_radii.get(element.to_upper()).unwrap_or(1.70)
}

///|
pub fn bp3_vdw_radius_list(elements : Array[String]) -> Array[Double] {
  let result = Array(capacity=elements.length())
  for i in 0..<elements.length() {
    result.push(bp3_vdw_radius(elements[i]))
  }
  result
}

///|
/// #24 PDBData — PDB data tables
let bp3_residue_mw : Map[String, Double] = Map([
  ("ALA", 89.09),
  ("ARG", 174.20),
  ("ASN", 132.12),
  ("ASP", 133.10),
  ("CYS", 121.16),
  ("GLN", 146.15),
  ("GLU", 147.13),
  ("GLY", 75.07),
  ("HIS", 155.16),
  ("ILE", 131.17),
  ("LEU", 131.17),
  ("LYS", 146.19),
  ("MET", 149.21),
  ("PHE", 165.19),
  ("PRO", 115.13),
  ("SER", 105.09),
  ("THR", 119.12),
  ("TRP", 204.23),
  ("TYR", 181.19),
  ("VAL", 117.15),
  ("DA", 331.22),
  ("DC", 307.20),
  ("DG", 347.22),
  ("DT", 322.21),
])

///|
let bp3_bond_lengths : Map[String, Double] = Map([
  ("C-N", 1.33),
  ("N-CA", 1.47),
  ("CA-C", 1.53),
  ("C-O", 1.23),
  ("CA-CB", 1.54),
  ("CB-CG", 1.52),
  ("CG-CD", 1.52),
  ("CD-CE", 1.52),
  ("N-H", 1.01),
  ("CA-HA", 1.09),
  ("O-HO", 0.96),
  ("S-S", 2.05),
  ("CA-SG", 1.81),
  ("CB-SG", 1.82),
  ("P-O", 1.60),
  ("P-O5", 1.59),
])

///|
let bp3_bond_angles : Map[String, Double] = Map([
  ("N-CA-C", 111.0),
  ("CA-C-N", 116.0),
  ("C-N-CA", 121.0),
  ("CA-C-O", 121.0),
  ("N-CA-CB", 110.0),
  ("CB-CA-C", 110.0),
  ("C-CA-CB", 110.0),
  ("O-C-N", 123.0),
])

///|
pub fn bp3_residue_mass(res_name : String) -> Double {
  bp3_residue_mw.get(res_name.to_upper()).unwrap_or(110.0)
}

///|
pub fn bp3_bond_length(atom1 : String, atom2 : String) -> Double {
  let key = atom1 + "-" + atom2
  let r = bp3_bond_lengths.get(key)
  match r {
    Some(v) => v
    None => bp3_bond_lengths.get(atom2 + "-" + atom1).unwrap_or(1.53)
  }
}

///|
pub fn bp3_bond_angle(atom1 : String, atom2 : String, atom3 : String) -> Double {
  let key = atom1 + "-" + atom2 + "-" + atom3
  let r = bp3_bond_angles.get(key)
  match r {
    Some(v) => v
    None => 109.5
  }
}

///|
pub fn bp3_protein_letters_3to1() -> Map[String, String] {
  Map([
    ("ALA", "A"),
    ("ARG", "R"),
    ("ASN", "N"),
    ("ASP", "D"),
    ("CYS", "C"),
    ("GLU", "E"),
    ("GLN", "Q"),
    ("GLY", "G"),
    ("HIS", "H"),
    ("ILE", "I"),
    ("LEU", "L"),
    ("LYS", "K"),
    ("MET", "M"),
    ("PHE", "F"),
    ("PRO", "P"),
    ("SER", "S"),
    ("THR", "T"),
    ("TRP", "W"),
    ("TYR", "Y"),
    ("VAL", "V"),
  ])
}

///|
pub fn bp3_protein_letters_1to3() -> Map[String, String] {
  Map([
    ("A", "ALA"),
    ("R", "ARG"),
    ("N", "ASN"),
    ("D", "ASP"),
    ("C", "CYS"),
    ("E", "GLU"),
    ("Q", "GLN"),
    ("G", "GLY"),
    ("H", "HIS"),
    ("I", "ILE"),
    ("L", "LEU"),
    ("K", "LYS"),
    ("M", "MET"),
    ("F", "PHE"),
    ("P", "PRO"),
    ("S", "SER"),
    ("T", "THR"),
    ("W", "TRP"),
    ("Y", "TYR"),
    ("V", "VAL"),
  ])
}

///|
pub fn bp3_iupac_get_extended_dna_values() -> Map[String, String] {
  bp3_iupac_extended_dna_values
}

///|
pub fn bp3_iupac_get_protein_letters() -> String {
  bp3_iupac_protein_letters
}

///|
pub fn bp3_iupac_get_extended_protein_letters() -> String {
  bp3_iupac_extended_protein_letters
}
