///|
/// QCP Superimposer demo - Protein structure superposition using Quaternion
/// Characteristic Polynomial algorithm.
/// Inspired by Biopython Bio.PDB.QCPSuperimposer.
///
/// Run: moon run examples/qcp_superimposer_demo/main.mbt
fn main {
  println("=== QCP Superimposer Demo ===")
  println("Inspired by Biopython Bio.PDB.QCPSuperimposer")
  println("")

  println("1. Create QCPAtomCoordinate objects")
  let atom1 = @bio.QCPAtomCoordinate::new(1.0, 0.0, 0.0)
  let atom2 = @bio.QCPAtomCoordinate::new(0.0, 1.0, 0.0)
  let atom3 = @bio.QCPAtomCoordinate::new(0.0, 0.0, 1.0)
  println("  Atom1: (" + atom1.x.to_string() + ", " + atom1.y.to_string() + ", " + atom1.z.to_string() + ")")
  println("  Atom2: (" + atom2.x.to_string() + ", " + atom2.y.to_string() + ", " + atom2.z.to_string() + ")")
  println("  Distance between atoms: " + atom1.distance(atom2).to_string())
  println("")

  println("2. Superimpose identical structures")
  let atoms1 = [atom1, atom2, atom3]
  let atoms2 = atoms1.copy()
  let superimposer1 = @bio.bio_qcp_superimpose(atoms1, atoms2)
  println("  RMSD (identical): " + superimposer1.get_rmsd().to_string())
  println("")

  println("3. Superimpose translated structures")
  let translated = [
    @bio.QCPAtomCoordinate::new(3.0, 2.0, 1.0),
    @bio.QCPAtomCoordinate::new(2.0, 3.0, 1.0),
    @bio.QCPAtomCoordinate::new(2.0, 2.0, 2.0),
  ]
  let superimposer2 = @bio.bio_qcp_superimpose(atoms1, translated)
  println("  RMSD (translated): " + superimposer2.get_rmsd().to_string())
  println("")

  println("4. Apply transformation to moving coordinates")
  let transformed = superimposer2.apply(translated)
  println("  Transformed coordinates:")
  let mut i = 0
  while i < transformed.length() {
    println("    (" + transformed[i].x.to_string() + ", " + transformed[i].y.to_string() + ", " + transformed[i].z.to_string() + ")")
    i = i + 1
  }
  println("")

  println("5. Direct RMSD calculation")
  let rmsd = @bio.bio_qcp_calculate_rmsd(atoms1, translated)
  println("  Calculated RMSD: " + rmsd.to_string())
  println("")

  println("=== QCP Superimposer Demo Complete ===")
}