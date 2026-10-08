///|
/// Atomic property lookup tables — masses, Van der Waals radii, covalent radii.
/// Ported from biotite.structure.info.masses / radii.
/// Data sources: IUPAC atomic weights (2021), Bondi VdW radii (1964),
/// Cordero covalent radii (2008).

///|
/// Atomic mass (Da) for elements H–U (1–92).
/// Returns 0.0 for unknown elements.
pub fn atom_mass(element : String) -> Double {
  let e = to_uppercase(element)
  match e {
    "H" => 1.008
    "HE" => 4.0026
    "LI" => 6.94
    "BE" => 9.0122
    "B" => 10.81
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "F" => 18.998
    "NE" => 20.180
    "NA" => 22.990
    "MG" => 24.305
    "AL" => 26.982
    "SI" => 28.085
    "P" => 30.974
    "S" => 32.06
    "CL" => 35.45
    "AR" => 39.948
    "K" => 39.098
    "CA" => 40.078
    "SC" => 44.956
    "TI" => 47.867
    "V" => 50.942
    "CR" => 51.996
    "MN" => 54.938
    "FE" => 55.845
    "CO" => 58.933
    "NI" => 58.693
    "CU" => 63.546
    "ZN" => 65.38
    "GA" => 69.723
    "GE" => 72.630
    "AS" => 74.922
    "SE" => 78.971
    "BR" => 79.904
    "KR" => 83.798
    "RB" => 85.468
    "SR" => 87.62
    "Y" => 88.906
    "ZR" => 90.224
    "NB" => 92.906
    "MO" => 95.95
    "TC" => 98.0
    "RU" => 101.07
    "RH" => 102.91
    "PD" => 106.42
    "AG" => 107.87
    "CD" => 112.41
    "IN" => 114.82
    "SN" => 118.71
    "SB" => 121.76
    "TE" => 127.60
    "I" => 126.90
    "XE" => 131.29
    "CS" => 132.91
    "BA" => 137.33
    "LA" => 138.91
    "CE" => 140.12
    "PR" => 140.91
    "ND" => 144.24
    "PM" => 145.0
    "SM" => 150.36
    "EU" => 151.96
    "GD" => 157.25
    "TB" => 158.93
    "DY" => 162.50
    "HO" => 164.93
    "ER" => 167.26
    "TM" => 168.93
    "YB" => 173.05
    "LU" => 174.97
    "HF" => 178.49
    "TA" => 180.95
    "W" => 183.84
    "RE" => 186.21
    "OS" => 190.23
    "IR" => 192.22
    "PT" => 195.08
    "AU" => 196.97
    "HG" => 200.59
    "TL" => 204.38
    "PB" => 207.2
    "BI" => 208.98
    "PO" => 209.0
    "AT" => 210.0
    "RN" => 222.0
    "FR" => 223.0
    "RA" => 226.0
    "AC" => 227.0
    "TH" => 232.04
    "PA" => 231.04
    "U" => 238.03
    _ => 0.0
  }
}

///|
/// Van der Waals radius (Å) — Bondi (1964) extended by Mantina et al. (2009).
/// Named atom_vdw_radius to avoid conflict with existing vdw_radius in pdb_analysis.
/// Returns 0.0 for unknown elements.
pub fn atom_vdw_radius(element : String) -> Double {
  let e = to_uppercase(element)
  match e {
    "H" => 1.20
    "HE" => 1.40
    "LI" => 1.82
    "BE" => 1.53
    "B" => 1.92
    "C" => 1.70
    "N" => 1.55
    "O" => 1.52
    "F" => 1.47
    "NE" => 1.54
    "NA" => 2.27
    "MG" => 1.73
    "AL" => 1.84
    "SI" => 2.10
    "P" => 1.80
    "S" => 1.80
    "CL" => 1.75
    "AR" => 1.88
    "K" => 2.75
    "CA" => 2.31
    "SC" => 2.30
    "TI" => 2.15
    "V" => 2.05
    "CR" => 2.05
    "MN" => 2.05
    "FE" => 2.04
    "CO" => 2.00
    "NI" => 1.97
    "CU" => 1.96
    "ZN" => 2.01
    "GA" => 1.87
    "GE" => 2.11
    "AS" => 1.85
    "SE" => 1.90
    "BR" => 1.85
    "KR" => 2.02
    "RB" => 3.03
    "SR" => 2.49
    "I" => 1.98
    "XE" => 2.16
    "CS" => 3.43
    "BA" => 2.68
    "AU" => 1.66
    "HG" => 1.70
    "PB" => 2.02
    _ => 0.0
  }
}

///|
/// Covalent radius (Å) — Cordero et al. (2008).
/// Returns 0.0 for unknown elements.
pub fn covalent_radius(element : String) -> Double {
  let e = to_uppercase(element)
  match e {
    "H" => 0.31
    "HE" => 0.28
    "LI" => 1.28
    "BE" => 0.96
    "B" => 0.84
    "C" => 0.76
    "N" => 0.71
    "O" => 0.66
    "F" => 0.57
    "NE" => 0.58
    "NA" => 1.66
    "MG" => 1.41
    "AL" => 1.21
    "SI" => 1.11
    "P" => 1.07
    "S" => 1.05
    "CL" => 1.02
    "AR" => 1.06
    "K" => 2.03
    "CA" => 1.76
    "SC" => 1.70
    "TI" => 1.60
    "V" => 1.53
    "CR" => 1.39
    "MN" => 1.39
    "FE" => 1.32
    "CO" => 1.26
    "NI" => 1.24
    "CU" => 1.32
    "ZN" => 1.22
    "GA" => 1.22
    "GE" => 1.20
    "AS" => 1.19
    "SE" => 1.20
    "BR" => 1.20
    "KR" => 1.16
    "RB" => 2.20
    "SR" => 1.95
    "I" => 1.39
    "XE" => 1.40
    "CS" => 2.44
    "BA" => 2.15
    "AU" => 1.36
    "HG" => 1.32
    "PB" => 1.46
    _ => 0.0
  }
}

///|
/// Formal charge for common amino-acid atom types.
/// Returns 0 for neutral/unknown.
pub fn formal_charge(resname : String, atom_name : String) -> Int {
  let r = to_uppercase(resname)
  let a = to_uppercase(atom_name)
  // Negative: carboxylate oxygens
  if (r == "ASP" || r == "GLU") && (a == "OD1" || a == "OD2" || a == "OE1" || a == "OE2") {
    return -1
  }
  // Positive: amine nitrogens
  if a == "NZ" && r == "LYS" {
    return 1
  }
  if a == "NH1" || a == "NH2" && r == "ARG" {
    return 1
  }
  if a == "NE" && r == "ARG" {
    return 1
  }
  if a == "ND1" || a == "NE2" {
    // HIS: depends on protonation state; default neutral
    return 0
  }
  // Terminal groups
  if a == "OXT" {
    return -1
  }
  0
}
