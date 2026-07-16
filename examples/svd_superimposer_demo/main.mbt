///|
/// SVD-based protein structure superposition example.
/// Inspired by Biopython Bio.PDB.SVDSuperimposer.
///
/// Run: moon run examples/svd_superimposer_demo/main.mbt
fn main {
  println("=== SVD Superimposer Demo ===")
  println("Inspired by Biopython Bio.PDB.SVDSuperimposer")
  println("")

  println("1. Create AtomCoordinate objects")
  let atom1 = @bio.AtomCoordinate::new(1.0, 0.0, 0.0)
  let atom2 = @bio.AtomCoordinate::new(0.0, 1.0, 0.0)
  let atom3 = @bio.AtomCoordinate::new(0.0, 0.0, 1.0)
  println("  Atom1: (\{atom1.x}, \{atom1.y}, \{atom1.z})")
  println("  Atom2: (\{atom2.x}, \{atom2.y}, \{atom2.z})")
  println("  Distance between atoms: \{atom1.distance(atom2)}")
  println("")

  println("2. Superimpose identical structures")
  let atoms1 = [atom1, atom2, atom3]
  let atoms2 = atoms1.copy()
  let (transformed1, rmsd1) = @bio.superimpose(atoms1, atoms2)
  println("  RMSD (identical): \{rmsd1}")
  println("")

  println("3. Superimpose translated structures")
  let translated = [
    @bio.AtomCoordinate::new(3.0, 2.0, 1.0),
    @bio.AtomCoordinate::new(2.0, 3.0, 1.0),
    @bio.AtomCoordinate::new(2.0, 2.0, 2.0),
  ]
  let (transformed2, rmsd2) = @bio.superimpose(atoms1, translated)
  println("  RMSD (translated): \{rmsd2}")
  println("")

  println("4. Superimpose rotated structures")
  let rotated = [
    @bio.AtomCoordinate::new(0.0, 1.0, 0.0),
    @bio.AtomCoordinate::new(-1.0, 0.0, 0.0),
    @bio.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let (transformed3, rmsd3) = @bio.superimpose(atoms1, rotated)
  println("  RMSD (rotated): \{rmsd3}")
  println("")

  println("5. Get rotation and translation matrices")
  let superimposer = @bio.SVDSuperimposer::set(atoms1, translated)
  let rot = superimposer.get_rotation()
  let tran = superimposer.get_translation()
  println("  Rotation matrix:")
  let mut i = 0
  while i < 3 {
    println("    \{rot[i]}")
    i = i + 1
  }
  println("  Translation vector: \{tran}")
  println("")

  println("6. Apply transformation")
  let applied = superimposer.apply(atoms1)
  println("  Applied transformation:")
  i = 0
  while i < applied.length() {
    println(
      "    Atom \{i}: (\{applied[i].x}, \{applied[i].y}, \{applied[i].z})",
    )
    i = i + 1
  }
}
