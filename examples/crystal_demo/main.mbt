///|
/// Bio.Crystal demo: small-molecule crystallographic CIF parsing.
///
/// Demonstrates unit-cell geometry, space-group lookup, coordinate
/// transformations, CIF parsing, and crystal-structure analysis.
fn main {
  println("=== Bio.Crystal Demo ===\n")

  // 1. Unit cell.
  println("1. Unit Cell:")
  let cell = @bio.UnitCell::new(
    a=5.418,
    b=5.418,
    c=7.543,
    alpha=90.0,
    beta=90.0,
    gamma=120.0,
  )
  println("   a=\{cell.a}, b=\{cell.b}, c=\{cell.c}")
  println("   alpha=\{cell.alpha}, beta=\{cell.beta}, gamma=\{cell.gamma}")
  let vol = @bio.unit_cell_volume(cell)
  println("   Volume: \{vol} A^3")

  // 2. Crystal density.
  println("\n2. Crystal Density:")
  let rho = @bio.crystal_density(cell, 2, 78.11)
  println("   Z=2, MW=78.11 g/mol")
  println("   Density: \{rho} g/cm^3")

  // 3. Space group lookup.
  println("\n3. Space Group Lookup:")
  let sg14 = @bio.lookup_space_group(14)
  println("   #14: \{sg14.symbol}")
  let sg19 = @bio.lookup_space_group(19)
  println("   #19: \{sg19.symbol}")
  let sg225 = @bio.lookup_space_group(225)
  println("   #225: \{sg225.symbol}")

  // 4. Coordinate transformations.
  println("\n4. Coordinate Transformations:")
  let(x, y, z) = @bio.fractional_to_cartesian(cell, 0.5, 0.5, 0.5)
  println("   Fractional (0.5, 0.5, 0.5) -> Cartesian (\{x}, \{y}, \{z})")
  let(xf, yf, zf) = @bio.cartesian_to_fractional(cell, x, y, z)
  println("   Cartesian back to Fractional (\{xf}, \{yf}, \{zf})")

  // 5. Parse a CIF file.
  println("\n5. Parse CIF File:")
  let cif_text = @bio.crystal_sample_cif()
  let blocks = @bio.parse_cif(cif_text)
  println("   Parsed \{blocks.length()} block(s)")
  println("   Block name: \{blocks[0].name}")

  // 6. Convert CIF block to CrystalStructure.
  println("\n6. Convert to CrystalStructure:")
  let s = @bio.cif_block_to_crystal(blocks[0])
  println("   Name: \{s.name}")
  println("   Atoms: \{s.atoms.length()}")
  println("   Z = \{s.z}")
  println("   MW = \{s.molecular_weight} g/mol")
  println("   Space group: \{s.space_group.symbol} (#\{s.space_group.number})")
  println("   Volume: \{@bio.unit_cell_volume(s.cell)} A^3")
  let density = @bio.crystal_density(s.cell, s.z, s.molecular_weight)
  println("   Density: \{density} g/cm^3")

  // 7. Atoms.
  println("\n7. Atom Sites:")
  println("   label\telement\tx\ty\tz")
  for i in 0..<s.atoms.length() {
    let a = s.atoms[i]
    println("   \{a.label}\t\{a.element}\t\{a.x}\t\{a.y}\t\{a.z}")
  }

  // 8. Cartesian atoms.
  println("\n8. Atoms in Cartesian Coordinates:")
  let cart_atoms = @bio.atoms_to_cartesian(s)
  for i in 0..<cart_atoms.length() {
    let a = cart_atoms[i]
    println(
      "   \{a.label}: (\{a.x}, \{a.y}, \{a.z}) [is_fractional=\{a.is_fractional}]",
    )
  }

  // 9. Crystal structure from sample.
  println("\n9. Sample Crystal Structure (benzene):")
  let sample = @bio.crystal_sample_structure()
  println("   Name: \{sample.name}")
  println("   Atoms: \{sample.atoms.length()}")
  println("   Bonds: \{sample.bonds.length()}")
  for i in 0..<sample.bonds.length() {
    let b = sample.bonds[i]
    println("   Bond: \{b.atom1} - \{b.atom2}, distance=\{b.distance} A")
  }

  // 10. Center of mass.
  println("\n10. Center of Mass:")
  let(cx, cy, cz) = @bio.center_of_mass(sample)
  println("   Center: (\{cx}, \{cy}, \{cz})")

  // 11. Atom distance.
  println("\n11. Atom Distances:")
  let d12 = @bio.crystal_atom_distance(sample, "C1", "C2")
  println("   C1-C2 distance: \{d12.unwrap_or(-1.0)} A")

  // 12. Format crystal structure.
  println("\n12. Format Crystal Structure:")
  let out = @bio.format_crystal_structure(sample)
  println(out)

  // 13. Crystal summary.
  println("13. Crystal Summary:")
  let summary = @bio.crystal_summary(sample)
  println("   Name: \{summary.name}")
  println("   Atoms: \{summary.n_atoms}")
  println("   Bonds: \{summary.n_bonds}")
  println("   Volume: \{summary.volume} A^3")
  println("   Density: \{summary.density} g/cm^3")
  println("   Space group: \{summary.space_group}")
  println("   Z: \{summary.z}")

  println("\n=== Demo Complete ===")
}
