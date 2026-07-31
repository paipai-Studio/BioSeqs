///|
/// Chemical Utilities Demo (Bio.PDB.chem_utils equivalent).
///
/// Demonstrates chemical property calculations for protein structures.
///
/// Run: moon run examples/chem_utils_demo/main.mbt
fn main {
  println("=== Chemical Utilities Demo ===\n")

  // Van der Waals radii
  println("1. Van der Waals Radii (in Angstroms):")
  println("   C: \{@bio.chem_vdw_radius("C")}")
  println("   N: \{@bio.chem_vdw_radius("N")}")
  println("   O: \{@bio.chem_vdw_radius("O")}")
  println("   S: \{@bio.chem_vdw_radius("S")}")
  println("   H: \{@bio.chem_vdw_radius("H")}")

  // Covalent radii
  println("\n2. Covalent Radii (in Angstroms):")
  println("   C: \{@bio.chem_covalent_radius("C")}")
  println("   N: \{@bio.chem_covalent_radius("N")}")
  println("   O: \{@bio.chem_covalent_radius("O")}")

  // Bond length
  println("\n3. Bond Length Calculation:")
  let atom1 = [0.0, 0.0, 0.0]
  let atom2 = [1.54, 0.0, 0.0] // ~C-C bond distance
  let bond_len = @bio.chem_bond_length(atom1, atom2)
  println("   Atom1: (0, 0, 0)")
  println("   Atom2: (1.54, 0, 0)")
  println("   Bond length: \{bond_len} Angstroms")

  // Bond angle
  println("\n4. Bond Angle Calculation:")
  let a = [1.0, 0.0, 0.0]
  let b = [0.0, 0.0, 0.0]
  let c = [0.0, 1.0, 0.0]
  let angle = @bio.chem_bond_angle(a, b, c)
  println("   Angle (right angle test): \{angle} degrees")

  let a2 = [1.0, 0.0, 0.0]
  let b2 = [0.0, 0.0, 0.0]
  let c2 = [-1.0, 0.0, 0.0]
  let angle2 = @bio.chem_bond_angle(a2, b2, c2)
  println("   Angle (collinear test): \{angle2} degrees")

  // Dihedral angle
  println("\n5. Dihedral Angle Calculation:")
  let d1 = [1.0, 0.0, 0.0]
  let d2 = [0.0, 0.0, 0.0]
  let d3 = [0.0, 1.0, 0.0]
  let d4 = [1.0, 1.0, 0.0]
  let dihedral = @bio.chem_dihedral_angle(d1, d2, d3, d4)
  println("   Dihedral angle (0 test): \{dihedral} degrees")

  // Empirical formula
  println("\n6. Empirical Formula:")
  let water_counts : Map[String, Int] = Map([], capacity=8)
  water_counts.set("H", 2)
  water_counts.set("O", 1)
  let formula = @bio.chem_empirical_formula(water_counts)
  println("   Water -> Formula: \{formula}")

  // Formula weight
  println("\n7. Formula Weight Calculation:")
  println("   H2O: \{@bio.chem_formula_weight("H2O")} g/mol")
  println("   C6H12O6 (glucose): \{@bio.chem_formula_weight("C6H12O6")} g/mol")
  println("   NaCl: \{@bio.chem_formula_weight("NaCl")} g/mol")

  // Isoelectric point
  println("\n8. Isoelectric Points (pI):")
  println("   K (Lysine): \{@bio.chem_isoelectric_point("K")}")
  println("   D (Aspartate): \{@bio.chem_isoelectric_point("D")}")
  println("   A (Alanine): \{@bio.chem_isoelectric_point("A")}")

  // Hydrogen bond length
  println("\n9. Hydrogen Bond Lengths:")
  println("   OH...OH: \{@bio.chem_hydrogen_bond_length("OH", "OH")} Angstroms")
  println("   NH...OH: \{@bio.chem_hydrogen_bond_length("NH", "OH")} Angstroms")

  println("\n=== Demo Complete ===")
}
