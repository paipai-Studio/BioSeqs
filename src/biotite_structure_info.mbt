///|
/// biotite_structure_info.mbt — residue/atom dictionaries used by the
/// Biotite-flavoured helpers.

///|
/// Three-letter to one-letter code for canonical amino acids.
pub fn biotite_one_letter_code(res_name : String) -> String {
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
/// One-letter to three-letter code for canonical amino acids.
pub fn biotite_three_letter_code(one : String) -> String {
  match one {
    "A" => "ALA"
    "R" => "ARG"
    "N" => "ASN"
    "D" => "ASP"
    "C" => "CYS"
    "Q" => "GLN"
    "E" => "GLU"
    "G" => "GLY"
    "H" => "HIS"
    "I" => "ILE"
    "L" => "LEU"
    "K" => "LYS"
    "M" => "MET"
    "F" => "PHE"
    "P" => "PRO"
    "S" => "SER"
    "T" => "THR"
    "W" => "TRP"
    "Y" => "TYR"
    "V" => "VAL"
    _ => "UNK"
  }
}

///|
/// Average residue mass in Da (free amino acid, not peptide-residue mass).
pub fn biotite_residue_mass(res_name : String) -> Double {
  match res_name {
    "ALA" => 89.09
    "ARG" => 174.20
    "ASN" => 132.12
    "ASP" => 133.10
    "CYS" => 121.16
    "GLN" => 146.15
    "GLU" => 147.13
    "GLY" => 75.07
    "HIS" => 155.16
    "ILE" => 131.17
    "LEU" => 131.17
    "LYS" => 146.19
    "MET" => 149.21
    "PHE" => 165.19
    "PRO" => 115.13
    "SER" => 105.09
    "THR" => 119.12
    "TRP" => 204.23
    "TYR" => 181.19
    "VAL" => 117.15
    _ => 0.0
  }
}

///|
/// Standard atomic weight for common elements (0.0 = unknown).
pub fn biotite_atom_mass(element : String) -> Double {
  match to_uppercase(element) {
    "H" => 1.008
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "P" => 30.974
    "S" => 32.06
    "F" => 18.998
    "CL" => 35.45
    "BR" => 79.904
    "I" => 126.904
    "FE" => 55.845
    "ZN" => 65.38
    "MG" => 24.305
    "CA" => 40.078
    "NA" => 22.990
    "K" => 39.098
    _ => 0.0
  }
}

///|
/// Net side-chain charge at pH 7 for the ionisable residues.
pub fn biotite_residue_charge(res_name : String) -> Int32 {
  match res_name {
    "ASP" | "GLU" => -1
    "LYS" | "ARG" => 1
    "HIS" => 0
    _ => 0
  }
}

///|
/// Total mass of an AtomArray using the element dictionary.
pub fn biotite_structure_mass(arr : AtomArray) -> Double {
  let mut total = 0.0
  for i in 0..<arr.len() {
    total = total + biotite_atom_mass(arr.elements[i])
  }
  total
}

///|
/// Van-der-Waals radius (Angstrom) of an element
/// (official `biotite.structure.info.vdw_radius_single`).
///
/// Main group values follow Mantina et al. 2009 (Table 12), the
/// transition metal subset follows RDKit 2024.09.4, exactly like the
/// official `_SINGLE_ATOM_VDW_RADII` table. Element lookup is
/// case-insensitive (the official implementation upper-cases the
/// input); `None` is returned for elements without an entry.
pub fn biotite_vdw_radius_single(element : String) -> Double? {
  match element.to_string().to_upper() {
    "H" => Some(1.10)
    "HE" => Some(1.40)
    "LI" => Some(1.81)
    "BE" => Some(1.53)
    "B" => Some(1.92)
    "C" => Some(1.70)
    "N" => Some(1.55)
    "O" => Some(1.52)
    "F" => Some(1.47)
    "NE" => Some(1.54)
    "NA" => Some(2.27)
    "MG" => Some(1.73)
    "AL" => Some(1.84)
    "SI" => Some(2.10)
    "P" => Some(1.80)
    "S" => Some(1.80)
    "CL" => Some(1.75)
    "AR" => Some(1.88)
    "K" => Some(2.75)
    "CA" => Some(2.31)
    "GA" => Some(1.87)
    "GE" => Some(2.11)
    "AS" => Some(1.85)
    "SE" => Some(1.90)
    "BR" => Some(1.83)
    "KR" => Some(2.02)
    "RB" => Some(3.03)
    "SR" => Some(2.49)
    "IN" => Some(1.93)
    "SN" => Some(2.17)
    "SB" => Some(2.06)
    "TE" => Some(2.06)
    "I" => Some(1.98)
    "XE" => Some(2.16)
    "CS" => Some(3.43)
    "BA" => Some(2.68)
    "TL" => Some(1.96)
    "PB" => Some(2.02)
    "BI" => Some(2.07)
    "PO" => Some(1.97)
    "AT" => Some(2.02)
    "RN" => Some(2.20)
    "FR" => Some(3.48)
    "RA" => Some(2.83)
    "FE" => Some(2.05)
    "CU" => Some(2.00)
    "ZN" => Some(2.10)
    "MN" => Some(2.05)
    "CO" => Some(2.00)
    "NI" => Some(2.00)
    "MO" => Some(2.10)
    "RU" => Some(2.05)
    "W" => Some(2.10)
    "PT" => Some(2.05)
    "AU" => Some(2.10)
    _ => None
  }
}
