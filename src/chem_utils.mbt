///|
/// Bio.PDB.chem_utils - Chemical utilities for protein structure analysis.
/// Inspired by Biopython's Bio.PDB.chem_utils module.
/// Provides bond length calculations, angle calculations, van der Waals radii,
/// covalent radii, ionic radii, and other chemical property calculations.
///
/// These utilities are fundamental for structure validation, energy calculations,
/// and chemical analysis in structural biology.

///|
/// Van der Waals radii by element (in Angstroms).
/// Values from Mantina et al. (2008) J. Phys. Chem. A.
pub fn chem_vdw_radius(element : String) -> Double {
  let el = element.to_upper()
  if el == "H" { 1.20 }
  else if el == "C" { 1.70 }
  else if el == "N" { 1.55 }
  else if el == "O" { 1.52 }
  else if el == "F" { 1.47 }
  else if el == "S" { 1.80 }
  else if el == "P" { 1.80 }
  else if el == "CL" || el == "Cl" { 1.75 }
  else if el == "BR" || el == "Br" { 1.85 }
  else if el == "I" { 1.98 }
  else if el == "FE" || el == "Fe" { 1.94 }
  else if el == "ZN" || el == "Zn" { 1.39 }
  else if el == "MG" || el == "Mg" { 1.73 }
  else if el == "CA" || el == "Ca" { 2.31 }
  else if el == "MN" || el == "Mn" { 1.97 }
  else if el == "CU" || el == "Cu" { 1.40 }
  else if el == "NA" || el == "Na" { 2.27 }
  else if el == "K" { 2.75 }
  else if el == "BR" { 1.85 }
  else { 1.70 } // Default to carbon
}

///|
/// Covalent radii by element (in Angstroms).
/// Values from Cordero et al. (2008) Dalton Trans.
pub fn chem_covalent_radius(element : String) -> Double {
  let el = element.to_upper()
  if el == "H" { 0.31 }
  else if el == "C" { 0.76 }
  else if el == "N" { 0.71 }
  else if el == "O" { 0.66 }
  else if el == "F" { 0.57 }
  else if el == "S" { 1.05 }
  else if el == "P" { 1.07 }
  else if el == "CL" || el == "Cl" { 1.02 }
  else if el == "BR" || el == "Br" { 1.2 }
  else if el == "I" { 1.39 }
  else if el == "FE" || el == "Fe" { 1.32 }
  else if el == "ZN" || el == "Zn" { 1.22 }
  else if el == "MG" || el == "Mg" { 1.41 }
  else if el == "CA" || el == "Ca" { 1.76 }
  else if el == "MN" || el == "Mn" { 1.39 }
  else if el == "CU" || el == "Cu" { 1.32 }
  else if el == "NA" || el == "Na" { 1.66 }
  else if el == "K" { 2.0 }
  else { 0.77 } // Default to carbon
}

///|
/// Ionic radii by element (in Angstroms) for common oxidation states.
/// Returns a tuple of (radius, oxidation_state).
pub fn chem_ionic_radius(element : String, oxidation_state? : Int = 0) -> Double {
  let el = element.to_upper()
  if el == "H" { 0.0 }
  else if el == "C" { 0.08 }
  else if el == "N" { 0.12 }
  else if el == "O" { 1.40 }
  else if el == "F" { 1.33 }
  else if el == "S" { 1.84 }
  else if el == "P" { 0.17 }
  else if el == "CL" || el == "Cl" { 1.81 }
  else if el == "BR" || el == "Br" { 1.96 }
  else if el == "I" { 2.20 }
  else if el == "FE" || el == "Fe" {
    if oxidation_state == 2 { 0.78 } else { 0.65 }
  }
  else if el == "ZN" || el == "Zn" { 0.74 }
  else if el == "MG" || el == "Mg" { 0.72 }
  else if el == "CA" || el == "Ca" { 1.00 }
  else if el == "MN" || el == "Mn" {
    if oxidation_state == 2 { 0.80 } else { 0.67 }
  }
  else if el == "CU" || el == "Cu" {
    if oxidation_state == 1 { 0.77 } else { 0.73 }
  }
  else if el == "NA" || el == "Na" { 1.02 }
  else if el == "K" { 1.38 }
  else { 0.0 }
}

///|
/// Calculate bond length between two atoms given their coordinates.
/// atom1_coords: [x, y, z] of first atom
/// atom2_coords: [x, y, z] of second atom
pub fn chem_bond_length(atom1_coords : Array[Double], atom2_coords : Array[Double]) -> Double {
  let dx = atom1_coords[0] - atom2_coords[0]
  let dy = atom1_coords[1] - atom2_coords[1]
  let dz = atom1_coords[2] - atom2_coords[2]
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// Calculate bond angle at atom2 between atom1-atom2-atom3 (in degrees).
/// atom1_coords, atom2_coords, atom3_coords: [x, y, z] arrays
pub fn chem_bond_angle(
  atom1_coords : Array[Double],
  atom2_coords : Array[Double],
  atom3_coords : Array[Double],
) -> Double {
  // Vectors from atom2 to atom1 and atom2 to atom3
  let v1x = atom1_coords[0] - atom2_coords[0]
  let v1y = atom1_coords[1] - atom2_coords[1]
  let v1z = atom1_coords[2] - atom2_coords[2]
  let v2x = atom3_coords[0] - atom2_coords[0]
  let v2y = atom3_coords[1] - atom2_coords[1]
  let v2z = atom3_coords[2] - atom2_coords[2]
  let dot = v1x * v2x + v1y * v2y + v1z * v2z
  let mag1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt()
  let mag2 = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt()
  if mag1 == 0.0 || mag2 == 0.0 {
    return 0.0
  }
  let cos_angle = dot / (mag1 * mag2)
  let clamped = if cos_angle > 1.0 { 1.0 } else if cos_angle < -1.0 { -1.0 } else { cos_angle }
  // acos approximation
  acos_approx(clamped) * 180.0 / 3.141592653589793
}

///|
/// Calculate dihedral (torsion) angle for four atoms (in degrees).
/// atom1, atom2, atom3, atom4: [x, y, z] arrays
pub fn chem_dihedral_angle(
  atom1_coords : Array[Double],
  atom2_coords : Array[Double],
  atom3_coords : Array[Double],
  atom4_coords : Array[Double],
) -> Double {
  // Vector b1 = atom2 - atom1
  let b1x = atom2_coords[0] - atom1_coords[0]
  let b1y = atom2_coords[1] - atom1_coords[1]
  let b1z = atom2_coords[2] - atom1_coords[2]
  // Vector b2 = atom3 - atom2
  let b2x = atom3_coords[0] - atom2_coords[0]
  let b2y = atom3_coords[1] - atom2_coords[1]
  let b2z = atom3_coords[2] - atom2_coords[2]
  // Vector b3 = atom4 - atom3
  let b3x = atom4_coords[0] - atom3_coords[0]
  let b3y = atom4_coords[1] - atom3_coords[1]
  let b3z = atom4_coords[2] - atom3_coords[2]
  // Normal n1 = b1 x b2
  let n1x = b1y * b2z - b1z * b2y
  let n1y = b1z * b2x - b1x * b2z
  let n1z = b1x * b2y - b1y * b2x
  // Normal n2 = b2 x b3
  let n2x = b2y * b3z - b2z * b3y
  let n2y = b2z * b3x - b2x * b3z
  let n2z = b2x * b3y - b2y * b3x
  // Normalize b2
  let mag_b2 = (b2x * b2x + b2y * b2y + b2z * b2z).sqrt()
  if mag_b2 == 0.0 {
    return 0.0
  }
  // Unit b2
  let ub2x = b2x / mag_b2
  let ub2y = b2y / mag_b2
  let ub2z = b2z / mag_b2
  // m1 = n1 x ub2
  let m1x = n1y * ub2z - n1z * ub2y
  let m1y = n1z * ub2x - n1x * ub2z
  let m1z = n1x * ub2y - n1y * ub2x
  // x = n1 . n2
  let x = n1x * n2x + n1y * n2y + n1z * n2z
  // y = m1 . n2
  let y = m1x * n2x + m1y * n2y + m1z * n2z
  // atan2(y, x)
  let angle = atan2_approx(y, x)
  angle * 180.0 / 3.141592653589793
}

///|
/// Check if two atoms are likely bonded based on distance and covalent radii.
pub fn chem_is_bonded(
  atom1_element : String,
  atom2_element : String,
  distance : Double,
  tolerance? : Double = 0.4,
) -> Bool {
  let r1 = chem_covalent_radius(atom1_element)
  let r2 = chem_covalent_radius(atom2_element)
  let expected_max = r1 + r2 + tolerance
  let expected_min = r1 + r2 - tolerance
  distance <= expected_max && distance >= expected_min
}

///|
/// Calculate total charge of a molecule given atomic charges.
pub fn chem_total_charge(charges : Array[Double]) -> Double {
  let mut total = 0.0
  for c in charges {
    total = total + c
  }
  total
}

///|
/// Calculate empirical formula from element counts.
pub fn chem_empirical_formula(element_counts : Map[String, Int]) -> String {
  let result = StringBuilder::new()
  // Standard Hill system: C first, H second, then alphabetical
  let elements = element_counts.keys().collect()
  let has_c = elements.contains("C")
  let has_h = elements.contains("H")
  if has_c {
    result.write_string("C")
    let c_count = element_counts["C"]
    if c_count > 1 {
      result.write_string(c_count.to_string())
    }
  }
  if has_h {
    result.write_string("H")
    let h_count = element_counts["H"]
    if h_count > 1 {
      result.write_string(h_count.to_string())
    }
  }
  // Remaining elements alphabetically
  let remaining = Array::new()
  for el in elements {
    if el != "C" && el != "H" {
      remaining.push(el)
    }
  }
  // Sort remaining alphabetically
  let mut sort_again = true
  while sort_again {
    sort_again = false
    let mut i = 0
    while i < remaining.length() - 1 {
      if remaining[i] > remaining[i + 1] {
        let tmp = remaining[i]
        remaining[i] = remaining[i + 1]
        remaining[i + 1] = tmp
        sort_again = true
      }
      i = i + 1
    }
  }
  for el in remaining {
    result.write_string(el)
    let count = element_counts[el]
    if count > 1 {
      result.write_string(count.to_string())
    }
  }
  result.to_string()
}

///|
/// Calculate molecular weight from element counts and atomic weights.
pub fn chem_molecular_weight(element_counts : Map[String, Int]) -> Double {
  let atomic_weights_elements = [
    "H", "C", "N", "O", "F", "P", "S", "CL", "BR", "I",
    "FE", "ZN", "MG", "CA", "MN", "CU", "NA", "K",
  ]
  let atomic_weights_values = [
    1.008, 12.011, 14.007, 15.999, 18.998, 30.974, 32.06, 35.45, 79.904, 126.904,
    55.845, 65.38, 24.305, 40.078, 54.938, 63.546, 22.990, 39.098,
  ]
  let mut weight = 0.0
  let entries = element_counts.to_array()
  let mut i = 0
  while i < entries.length() {
    let entry = entries[i]
    let el = entry.0
    let count = entry.1
    let upper_el = el.to_upper()
    let mut aw = 0.0
    let mut j = 0
    while j < atomic_weights_elements.length() {
      if atomic_weights_elements[j] == upper_el {
        aw = atomic_weights_values[j]
        break
      }
      j = j + 1
    }
    weight = weight + aw * count.to_double()
    i = i + 1
  }
  weight
}

///|
/// Approximate arccosine function.
/// Uses the identity: acos(x) = atan2(sqrt(1-x^2), x)
fn acos_approx(x : Double) -> Double {
  if x >= 1.0 { return 0.0 }
  if x <= -1.0 { return 3.141592653589793 }
  atan2_approx((1.0 - x * x).sqrt(), x)
}

///|
/// Approximate atan2 function using a simple polynomial approximation.
fn atan2_approx(y : Double, x : Double) -> Double {
  if x > 0.0 {
    atan_approx(y / x)
  } else if x < 0.0 {
    if y >= 0.0 {
      atan_approx(y / x) + 3.141592653589793
    } else {
      atan_approx(y / x) - 3.141592653589793
    }
  } else {
    if y > 0.0 { 3.141592653589793 / 2.0 }
    else if y < 0.0 { -3.141592653589793 / 2.0 }
    else { 0.0 }
  }
}

///|
/// Approximate arctangent for a single value.
fn atan_approx(x : Double) -> Double {
  // Range reduction
  if x < -1.0 {
    return -3.141592653589793 / 2.0 + atan_approx(1.0 / x)
  } else if x > 1.0 {
    return 3.141592653589793 / 2.0 - atan_approx(1.0 / x)
  }
  // Polynomial approximation for |x| <= 1
  let x2 = x * x
  let x3 = x2 * x
  let x5 = x3 * x2
  let x7 = x5 * x2
  let x9 = x7 * x2
  x - x3 / 3.0 + x5 / 5.0 - x7 / 7.0 + x9 / 9.0
}

///|
/// Get standard amino acid element composition.
pub fn chem_amino_acid_elements(aa_name : String) -> Map[String, Int] {
  let name = aa_name.to_upper()
  let elements : Map[String, Int] = Map([], capacity=16)
  if name == "ALA" || name == "A" {
    elements.set("C", 3); elements.set("H", 7); elements.set("N", 1); elements.set("O", 2)
  } else if name == "ARG" || name == "R" {
    elements.set("C", 6); elements.set("H", 14); elements.set("N", 4); elements.set("O", 2)
  } else if name == "ASN" || name == "N" {
    elements.set("C", 4); elements.set("H", 8); elements.set("N", 2); elements.set("O", 3)
  } else if name == "ASP" || name == "D" {
    elements.set("C", 4); elements.set("H", 7); elements.set("N", 1); elements.set("O", 4)
  } else if name == "CYS" || name == "C" {
    elements.set("C", 3); elements.set("H", 7); elements.set("N", 1); elements.set("O", 2); elements.set("S", 1)
  } else if name == "GLN" || name == "Q" {
    elements.set("C", 5); elements.set("H", 10); elements.set("N", 2); elements.set("O", 3)
  } else if name == "GLU" || name == "E" {
    elements.set("C", 5); elements.set("H", 9); elements.set("N", 1); elements.set("O", 4)
  } else if name == "GLY" || name == "G" {
    elements.set("C", 2); elements.set("H", 5); elements.set("N", 1); elements.set("O", 2)
  } else if name == "HIS" || name == "H" {
    elements.set("C", 6); elements.set("H", 9); elements.set("N", 3); elements.set("O", 2)
  } else if name == "ILE" || name == "I" {
    elements.set("C", 6); elements.set("H", 13); elements.set("N", 1); elements.set("O", 2)
  } else if name == "LEU" || name == "L" {
    elements.set("C", 6); elements.set("H", 13); elements.set("N", 1); elements.set("O", 2)
  } else if name == "LYS" || name == "K" {
    elements.set("C", 6); elements.set("H", 14); elements.set("N", 2); elements.set("O", 2)
  } else if name == "MET" || name == "M" {
    elements.set("C", 5); elements.set("H", 11); elements.set("N", 1); elements.set("O", 2); elements.set("S", 1)
  } else if name == "PHE" || name == "F" {
    elements.set("C", 9); elements.set("H", 11); elements.set("N", 1); elements.set("O", 2)
  } else if name == "PRO" || name == "P" {
    elements.set("C", 5); elements.set("H", 9); elements.set("N", 1); elements.set("O", 2)
  } else if name == "SER" || name == "S" {
    elements.set("C", 3); elements.set("H", 7); elements.set("N", 1); elements.set("O", 3)
  } else if name == "THR" || name == "T" {
    elements.set("C", 4); elements.set("H", 9); elements.set("N", 1); elements.set("O", 3)
  } else if name == "TRP" || name == "W" {
    elements.set("C", 11); elements.set("H", 12); elements.set("N", 2); elements.set("O", 2)
  } else if name == "TYR" || name == "Y" {
    elements.set("C", 9); elements.set("H", 11); elements.set("N", 1); elements.set("O", 3)
  } else if name == "VAL" || name == "V" {
    elements.set("C", 5); elements.set("H", 11); elements.set("N", 1); elements.set("O", 2)
  }
  elements
}

///|
/// Get the average molecular weight of an amino acid residue.
/// Monoisotopic weight (minus water for peptide bond).
pub fn chem_aa_molecular_weight(aa_name : String) -> Double {
  let elements = chem_amino_acid_elements(aa_name)
  let weight = chem_molecular_weight(elements)
  // Subtract water (18.015) for each peptide bond
  // But here we just return the amino acid MW
  weight
}

///|
/// Get average molecular weight of an amino acid residue (minus water).
pub fn chem_aa_residue_weight(aa_name : String) -> Double {
  let full_weight = chem_aa_molecular_weight(aa_name)
  full_weight - 18.01524 // Subtract water
}

///|
/// Calculate extinction coefficient at 280nm for a protein sequence.
/// Uses Pace et al. (1995) method with Trp, Tyr, and Cys counts.
pub fn chem_extinction_coefficient(sequence : String) -> Double {
  let upper_seq = sequence.to_upper()
  let mut trp_count = 0
  let mut tyr_count = 0
  let mut cys_count = 0
  let mut i = 0
  while i < upper_seq.length() {
    let ch = upper_seq.unsafe_get(i).unsafe_to_char()
    if ch == 'W' { trp_count = trp_count + 1 }
    else if ch == 'Y' { tyr_count = tyr_count + 1 }
    else if ch == 'C' { cys_count = cys_count + 1 }
    i = i + 1
  }
  // Pace et al. values
  let ext_trp = 5500.0
  let ext_tyr = 1490.0
  let ext_cystine = 125.0 // per disulfide bond (divide Cys by 2)
  let n_cystine = cys_count / 2
  ext_trp * trp_count.to_double() + ext_tyr * tyr_count.to_double() + ext_cystine * n_cystine.to_double()
}

///|
/// Calculate isoelectric point (pI) of a protein using pKa values.
/// Implements a simple pI calculation based on charge balance.
pub fn chem_isoelectric_point(
  sequence : String,
  p_nterm? : Double = 7.47,
  p_cterm? : Double = 3.55,
) -> Double {
  let upper_seq = sequence.to_upper()
  let mut pos_charge = 0.0
  let mut neg_charge = 0.0
  // Count charged residues
  let arg_count = count_char(upper_seq, 'R')
  let lys_count = count_char(upper_seq, 'K')
  let his_count = count_char(upper_seq, 'H')
  let asp_count = count_char(upper_seq, 'D')
  let glu_count = count_char(upper_seq, 'E')
  let cys_count = count_char(upper_seq, 'C')
  let tyr_count = count_char(upper_seq, 'Y')
  // pKa values
  let pka_arg = 12.48
  let pka_lys = 10.54
  let pka_his = 6.04
  let pka_asp = 3.57
  let pka_glu = 4.15
  let pka_cys = 8.00
  let pka_tyr = 10.10
  // Use bisection method to find pH where net charge = 0
  let mut low = 0.0
  let mut high = 14.0
  let mut pH = 7.0
  let mut iter = 0
  while iter < 20 {
    pH = (low + high) / 2.0
    // Calculate net charge at this pH
    pos_charge = 0.0
    neg_charge = 0.0
    // Positive charges
    pos_charge = pos_charge + 1.0 / (1.0 + chem_pow10(pH - p_nterm)) // N-term
    pos_charge = pos_charge + arg_count.to_double() / (1.0 + chem_pow10(pH - pka_arg))
    pos_charge = pos_charge + lys_count.to_double() / (1.0 + chem_pow10(pH - pka_lys))
    pos_charge = pos_charge + his_count.to_double() / (1.0 + chem_pow10(pH - pka_his))
    // Negative charges
    neg_charge = neg_charge + 1.0 / (1.0 + chem_pow10(p_cterm - pH)) // C-term
    neg_charge = neg_charge + asp_count.to_double() / (1.0 + chem_pow10(pka_asp - pH))
    neg_charge = neg_charge + glu_count.to_double() / (1.0 + chem_pow10(pka_glu - pH))
    neg_charge = neg_charge + cys_count.to_double() / (1.0 + chem_pow10(pka_cys - pH))
    neg_charge = neg_charge + tyr_count.to_double() / (1.0 + chem_pow10(pka_tyr - pH))
    let net = pos_charge - neg_charge
    if net > 0.0 {
      low = pH
    } else {
      high = pH
    }
    iter = iter + 1
  }
  pH
}

///|
/// Count occurrences of a character in a string.
fn count_char(s : String, ch : Char) -> Int {
  let mut count = 0
  let mut i = 0
  let target = ch.to_int()
  while i < s.length() {
    if s.unsafe_get(i).to_int() == target {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Calculate 10^x using exp(x * ln(10)).
fn chem_pow10(x : Double) -> Double {
  @math.exp(x * @math.ln(10.0))
}

///|
/// Calculate total number of atoms in a protein sequence.
pub fn chem_count_atoms(sequence : String) -> Int {
  let upper_seq = sequence.to_upper()
  let mut total = 0
  let mut i = 0
  while i < upper_seq.length() {
    let aa = upper_seq.unsafe_get(i).unsafe_to_char().to_string()
    let elements = chem_amino_acid_elements(aa)
    let el_entries = elements.to_array()
    let mut ei = 0
    while ei < el_entries.length() {
      total = total + el_entries[ei].1
      ei = ei + 1
    }
    i = i + 1
  }
  // Adjust for peptide bonds (subtract water for each bond)
  let n_bonds = if sequence.length() > 0 { sequence.length() - 1 } else { 0 }
  total - 3 * n_bonds // Remove H2O atoms (2H + 1O) per peptide bond
}

///|
/// Calculate formula weight from a chemical formula string.
/// Parses formulas like "H2O", "C6H12O6", "NaCl", etc.
pub fn chem_formula_weight(formula : String) -> Double {
  let atomic_weights_elements = [
    "H", "C", "N", "O", "F", "P", "S", "CL", "BR", "I",
    "FE", "ZN", "MG", "CA", "MN", "CU", "NA", "K",
  ]
  let atomic_weights_values = [
    1.008, 12.011, 14.007, 15.999, 18.998, 30.974, 32.06, 35.45, 79.904, 126.904,
    55.845, 65.38, 24.305, 40.078, 54.938, 63.546, 22.990, 39.098,
  ]
  let weight_map : Map[String, Double] = Map([], capacity=32)
  let mut i = 0
  while i < atomic_weights_elements.length() {
    weight_map.set(atomic_weights_elements[i], atomic_weights_values[i])
    i = i + 1
  }
  // Parse formula (convert to uppercase for consistent matching)
  let upper = formula.to_upper()
  let mut total_weight = 0.0
  let mut pos = 0
  while pos < upper.length() {
    let ch = upper.unsafe_get(pos).to_int()
    if ch >= 65 && ch <= 90 {
      // Start of element symbol
      let mut el = upper.unsafe_get(pos).unsafe_to_char().to_string()
      pos = pos + 1
      // Try two-letter element first (e.g., "CL", "BR", "FE", "NA", "K" is single letter)
      if pos < upper.length() {
        let next_ch = upper.unsafe_get(pos).to_int()
        if next_ch >= 65 && next_ch <= 90 {
          // Check if two-letter combination is in weight map
          let two_letter = el + upper.unsafe_get(pos).unsafe_to_char().to_string()
          if weight_map.contains(two_letter) {
            el = two_letter
            pos = pos + 1
          }
        }
      }
      // Read count
      let mut count = 0
      while pos < upper.length() {
        let digit = upper.unsafe_get(pos).to_int()
        if digit >= 48 && digit <= 57 {
          count = count * 10 + (digit - 48)
          pos = pos + 1
        } else {
          break
        }
      }
      if count == 0 {
        count = 1
      }
      // Look up atomic weight
      let mut aw = 0.0
      match weight_map.get(el) {
        Some(val) => aw = val
        None => ()
      }
      total_weight = total_weight + aw * count.to_double()
    } else {
      pos = pos + 1
    }
  }
  total_weight
}

///|
/// Calculate typical hydrogen bond length between two atom types.
/// Returns approximate distance in Angstroms.
pub fn chem_hydrogen_bond_length(donor_type : String, acceptor_type : String) -> Double {
  let donor = donor_type.to_upper()
  let acceptor = acceptor_type.to_upper()
  let key = donor + "..." + acceptor
  // Typical hydrogen bond distances (Angstroms)
  if key == "OH...OH" || key == "OH...O" { 2.72 }
  else if key == "NH...OH" || key == "NH...O" { 2.93 }
  else if key == "SH...SH" { 3.45 }
  else if key == "OH...NH" || key == "O...HN" { 2.88 }
  else if key == "NH...N" { 3.10 }
  else if key == "OH...S" { 3.40 }
  else if key == "NH...S" { 3.45 }
  else { 2.90 } // Default
}
