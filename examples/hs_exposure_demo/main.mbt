///|
/// HSExposure demo.
///
/// Run: moon run examples/hs_exposure_demo/main.mbt
fn main {
  println("=== HSExposure Demo ===\n")

  println("1. Create PDB atoms:")
  let ca = @bio.PDBAtom::new("CA", 0.0, 0.0, 0.0)
  let cb = @bio.PDBAtom::new("CB", 1.5, 0.0, 0.0)
  let n = @bio.PDBAtom::new("N", 0.0, 1.5, 0.0)
  println("   CA: \{ca.name} (\{ca.x}, \{ca.y}, \{ca.z})")
  println("   CB: \{cb.name} (\{cb.x}, \{cb.y}, \{cb.z})")
  println("   N: \{n.name} (\{n.x}, \{n.y}, \{n.z})")

  println("\n2. Prepare all CA atoms:")
  let all_ca : Array[@bio.PDBAtom] = Array::new()
  all_ca.push(ca)
  all_ca.push(@bio.PDBAtom::new("CA2", 5.0, 0.0, 0.0))
  all_ca.push(@bio.PDBAtom::new("CA3", 0.0, 5.0, 0.0))
  all_ca.push(@bio.PDBAtom::new("CA4", 0.0, 0.0, 5.0))
  all_ca.push(@bio.PDBAtom::new("CA5", 3.0, 3.0, 0.0))
  println("   Total CA atoms: \{all_ca.length()}")

  println("\n3. Calculate HSE:")
  let result = @bio.calculate_hse(ca, cb, n, all_ca)
  println("   Residue: \{result.residue_name}")
  println("   HSE up: \{result.hse_up}")
  println("   HSE down: \{result.hse_down}")
  println("   HSE total: \{result.hse_total}")

  println("\n4. Classify exposure:")
  let classification = @bio.classify_exposure(result.hse_total)
  println("   Classification: \{classification}")

  println("\n5. Example atoms:")
  let example_atoms = @bio.create_example_atoms()
  println("   Example atoms count: \{example_atoms.length()}")

  println("\n=== Demo Complete ===")
}