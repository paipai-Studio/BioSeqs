///|
/// Internal Coordinates Demo - Protein structure representation using
/// bond lengths, bond angles, and torsion angles.
/// Inspired by Biopython's Bio.PDB.internal_coords module.

///|
fn main {
  println("=== Internal Coordinates Demo ===")
  println("Inspired by Biopython's Bio.PDB.internal_coords module")
  println("")

  println("1. Dihedral Angle Calculation")
  println("   Computing dihedral angle for four coplanar points...")
  let p1 = [0.0, 0.0, 0.0]
  let p2 = [1.0, 0.0, 0.0]
  let p3 = [1.0, 1.0, 0.0]
  let p4 = [0.0, 1.0, 0.0]
  let angle = @bio.ic_dihedral_angle(p1, p2, p3, p4)
  let angle_deg = @bio.ic_rad_to_deg(angle)
  println("   Dihedral angle (radians): \{angle.to_string()}")
  println("   Dihedral angle (degrees): \{angle_deg.to_string()}")
  println("   PI constant: \{@bio.ic_pi().to_string()}")
  println("")

  println("2. Convert between radians and degrees")
  let rad = @bio.ic_deg_to_rad(180.0)
  println("   180 degrees in radians: \{rad.to_string()}")
  let deg = @bio.ic_rad_to_deg(@bio.ic_pi())
  println("   PI radians in degrees: \{deg.to_string()}")
  println("")

  println("3. Build Extended Polypeptide Chain")
  let sequence = ["ALA", "GLY", "SER"]
  println("   Building extended chain for: ALA-GLY-SER")
  let result = @bio.ic_build_extended_chain(sequence)
  println("   Number of atoms: \{result.x.length().to_string()}")
  println("   Number of residues: \{result.residue_ids.length().to_string()}")
  println("   First atom: \{result.atom_names[0]} at (\{result.x[0].to_string()}, \{result.y[0].to_string()}, \{result.z[0].to_string()})")
  println("")

  println("4. Backbone Torsion Angles (phi/psi)")
  let c_prev = [-1.0, 0.0, 0.0]
  let n_curr = [0.0, 0.0, 0.0]
  let ca_curr = [1.0, 0.0, 0.0]
  let c_curr = [2.0, 0.0, 0.0]
  let n_next = [3.0, 0.5, 0.0]
  let phi = @bio.ic_compute_phi(c_prev, n_curr, ca_curr, c_curr)
  let psi = @bio.ic_compute_psi(n_curr, ca_curr, c_curr, n_next)
  println("   Phi angle (rad): \{phi.to_string()}, (deg): \{@bio.ic_rad_to_deg(phi).to_string()}")
  println("   Psi angle (rad): \{psi.to_string()}, (deg): \{@bio.ic_rad_to_deg(psi).to_string()}")
  println("")

  println("5. Standard Bond Geometry Parameters")
  println("   Peptide bond length (C-N): \{@bio.ic_peptide_bond_length().to_string()} Angstroms")
  println("   CA-C bond length: \{@bio.ic_ca_c_bond_length().to_string()} Angstroms")
  println("   N-CA bond length: \{@bio.ic_n_ca_bond_length().to_string()} Angstroms")
  println("")

  println("6. Ramachandran Plot Data (phi/psi statistics)")
  let phi_vals = [-0.8, -0.5, -0.3, 0.2, 0.5, 0.7]
  let psi_vals = [-0.6, -0.2, 0.1, 0.3, 0.6, 0.8]
  let (mean_phi, mean_psi, var_phi, var_psi) = @bio.ic_ramachandran_data(phi_vals, psi_vals)
  println("   Mean phi: \{mean_phi.to_string()} rad (\{@bio.ic_rad_to_deg(mean_phi).to_string()} deg)")
  println("   Mean psi: \{mean_psi.to_string()} rad (\{@bio.ic_rad_to_deg(mean_psi).to_string()} deg)")
  println("   Variance phi: \{var_phi.to_string()}")
  println("   Variance psi: \{var_psi.to_string()}")
  println("")

  println("7. Rotamer Library for Side Chains")
  let amino_acids = ["LEU", "VAL", "PHE", "GLY", "ALA", "SER"]
  for aa in amino_acids {
    let rotamers = @bio.ic_chi1_rotamers(aa)
    println("   \{aa}: \{rotamers.length().to_string()} rotamers")
    for r in rotamers {
      println("     - \{r.name}: probability = \{r.probability.to_string()}, chi angles = \{r.chi_angles.length().to_string()}")
    }
  }
  println("")

  println("8. InternalCoordChain Operations")
  let chain = @bio.InternalCoordChain::new(chain_id='A')
  let summary = @bio.ic_chain_summary(chain)
  println("   Empty chain summary: \{summary}")
  let residue = @bio.InternalCoordResidue::new(resname="ALA", resseq=1)
  let chain2 = chain.add_residue(residue)
  println("   Chain after adding ALA: \{@bio.ic_chain_summary(chain2)}")
  println("")

  println("=== Demo completed successfully! ===")
  println("")
  println("Key concepts:")
  println("- Internal coordinates (bond lengths, angles, torsions) are")
  println("  an alternative to Cartesian coordinates for protein structure")
  println("- Dihedral angles (phi/psi/omega) define backbone conformation")
  println("- Side-chain rotamers define preferred chi angle combinations")
  println("- Ramachandran plots show allowed phi/psi angle regions")
}
