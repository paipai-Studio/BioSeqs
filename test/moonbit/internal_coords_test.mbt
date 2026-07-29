///|
/// Tests for Bio.PDB.internal_coords module.

///|
test "ic_dihedral_angle_planar" {
  // Four coplanar points in a rectangle give 0 dihedral angle
  let p1 = [0.0, 0.0, 0.0]
  let p2 = [1.0, 0.0, 0.0]
  let p3 = [1.0, 1.0, 0.0]
  let p4 = [0.0, 1.0, 0.0]

  let angle = @src.ic_dihedral_angle(p1, p2, p3, p4)
  assert_true(angle.abs() < 1.0e-10)
}

///|
test "ic_dihedral_angle_pi" {
  // Four coplanar points on opposite sides give PI
  let p1 = [0.0, 1.0, 0.0]
  let p2 = [0.0, 0.0, 0.0]
  let p3 = [1.0, 0.0, 0.0]
  let p4 = [1.0, -1.0, 0.0]

  let angle = @src.ic_dihedral_angle(p1, p2, p3, p4)
  assert_true((angle.abs() - @src.ic_pi()).abs() < 1.0e-10)
}

///|
test "ic_dihedral_angle_zero" {
  // All points on a line should give 0
  let p1 = [0.0, 0.0, 0.0]
  let p2 = [1.0, 0.0, 0.0]
  let p3 = [2.0, 0.0, 0.0]
  let p4 = [3.0, 0.0, 0.0]

  let angle = @src.ic_dihedral_angle(p1, p2, p3, p4)
  assert_true(angle.abs() < 1.0e-10)
}

///|
test "ic_deg_to_rad" {
  let rad = @src.ic_deg_to_rad(180.0)
  assert_true((rad - @src.ic_pi()).abs() < 1.0e-10)
}

///|
test "ic_rad_to_deg" {
  let deg = @src.ic_rad_to_deg(@src.ic_pi())
  assert_true((deg - 180.0).abs() < 1.0e-10)
}

///|
test "ic_vec_norm" {
  let v = [3.0, 4.0, 0.0]
  let norm = @src.ic_vec_norm(v)
  assert_true((norm - 5.0).abs() < 1.0e-10)
}

///|
test "ic_vec_dot" {
  let v1 = [1.0, 0.0, 0.0]
  let v2 = [0.0, 1.0, 0.0]
  let dot = @src.ic_dot(v1, v2)
  assert_true(dot.abs() < 1.0e-10)
}

///|
test "ic_vec_cross" {
  let v1 = [1.0, 0.0, 0.0]
  let v2 = [0.0, 1.0, 0.0]
  let cross = @src.ic_cross(v1, v2)
  assert_true((cross[2] - 1.0).abs() < 1.0e-10)
}

///|
test "ic_distance" {
  let d = @src.ic_distance([0.0, 0.0, 0.0], [3.0, 4.0, 0.0])
  assert_true((d - 5.0).abs() < 1.0e-10)
}

///|
test "ic_torsion_angle_creation" {
  let tau = @src.TorsionAngle::new(
    name="phi",
    value=-0.57,
    atom_names=["C", "N", "CA", "C"],
  )
  assert_eq(tau.name, "phi")
  assert_true((tau.value - (-0.57)).abs() < 1.0e-10)
}

///|
test "ic_bond_length_creation" {
  let bl = @src.BondLength::new(atom1="N", atom2="CA", length=1.47)
  assert_eq(bl.atom1, "N")
  assert_eq(bl.atom2, "CA")
  assert_true((bl.length - 1.47).abs() < 1.0e-10)
}

///|
test "ic_bond_angle_creation" {
  let ba = @src.BondAngle::new(atom1="N", atom2="CA", atom3="C", angle=110.0)
  assert_eq(ba.atom1, "N")
  assert_eq(ba.atom2, "CA")
  assert_eq(ba.atom3, "C")
  assert_true((ba.angle - 110.0).abs() < 1.0e-10)
}

///|
test "ic_residue_creation" {
  let res = @src.InternalCoordResidue::new(resname="ALA", resseq=1)
  assert_eq(res.resname, "ALA")
  assert_eq(res.resseq, 1)
  assert_eq(res.start_atom, "N")
}

///|
test "ic_chain_creation" {
  let chain = @src.InternalCoordChain::new(chain_id='A')
  assert_eq(chain.chain_id, 'A')
  assert_eq(chain.n_residues(), 0)
}

///|
test "ic_add_residue_to_chain" {
  let chain = @src.InternalCoordChain::new(chain_id='A')
  let res = @src.InternalCoordResidue::new(resname="ALA", resseq=1)
  let chain2 = chain.add_residue(res)
  assert_eq(chain2.n_residues(), 1)
}

///|
test "ic_compute_phi" {
  let phi = @src.ic_compute_phi(
    [0.0, 0.0, 0.0],   // prev C
    [1.5, 0.0, 0.0],   // N
    [3.0, 0.0, 0.0],   // CA
    [4.5, 1.0, 0.0],   // C
  )
  assert_true(phi >= -@src.ic_pi() && phi <= @src.ic_pi())
}

///|
test "ic_compute_psi" {
  let psi = @src.ic_compute_psi(
    [0.0, 0.0, 0.0],   // N
    [1.5, 0.0, 0.0],   // CA
    [3.0, 0.0, 0.0],   // C
    [4.5, 1.0, 0.0],   // next N
  )
  assert_true(psi >= -@src.ic_pi() && psi <= @src.ic_pi())
}

///|
test "ic_build_extended_chain" {
  let sequence = ["ALA", "GLY", "SER"]
  let result = @src.ic_build_extended_chain(sequence)
  assert_true(result.atom_names.length() > 0)
  assert_eq(result.atom_names.length(), result.x.length())
  assert_eq(result.atom_names.length(), result.y.length())
  assert_eq(result.atom_names.length(), result.z.length())
}

///|
test "ic_build_empty_chain" {
  let sequence : Array[String] = []
  let result = @src.ic_build_extended_chain(sequence)
  assert_eq(result.atom_names.length(), 0)
}

///|
test "ic_standard_bond_lengths" {
  let peptide = @src.ic_peptide_bond_length()
  assert_true((peptide - 1.33).abs() < 1.0e-10)
  let ca_c = @src.ic_ca_c_bond_length()
  assert_true((ca_c - 1.52).abs() < 1.0e-10)
  let n_ca = @src.ic_n_ca_bond_length()
  assert_true((n_ca - 1.47).abs() < 1.0e-10)
}

///|
test "ic_standard_angles" {
  let omega = @src.ic_omega_angle()
  assert_true((omega - 180.0).abs() < 1.0e-10)
  let tau = @src.ic_tau_angle()
  assert_true((tau - 110.0).abs() < 1.0e-10)
}

///|
test "ic_rotamer_creation" {
  let rotamer = @src.Rotamer::new(
    name="g+",
    chi_angles=[0.5, -0.3],
    probability=0.95,
  )
  assert_eq(rotamer.name, "g+")
  assert_true((rotamer.probability - 0.95).abs() < 1.0e-10)
}

///|
test "ic_rotamer_library_entry" {
  let entry = @src.RotamerLibraryEntry::new(phi=-0.57, psi=-0.45)
  assert_true((entry.phi - (-0.57)).abs() < 1.0e-10)
  assert_true((entry.psi - (-0.45)).abs() < 1.0e-10)
}

///|
test "ic_to_cartesian_result_creation" {
  let result = @src.ICToCartesianResult::new()
  assert_eq(result.atom_names.length(), 0)
  assert_eq(result.x.length(), 0)
  assert_eq(result.y.length(), 0)
  assert_eq(result.z.length(), 0)
}

///|
test "ic_chain_summary" {
  let chain = @src.InternalCoordChain::new(chain_id='A')
  let summary = @src.ic_chain_summary(chain)
  assert_true(summary.length() > 0)
}

///|
test "ic_ramachandran_data" {
  let phi_vals = [-1.0, -0.5, 0.0, 0.5, 1.0]
  let psi_vals = [-1.0, -0.5, 0.0, 0.5, 1.0]
  let (mean_phi, mean_psi, var_phi, var_psi) = @src.ic_ramachandran_data(phi_vals, psi_vals)
  assert_true(mean_phi.abs() < 1.0e-10)
  assert_true(mean_psi.abs() < 1.0e-10)
  assert_true(var_phi > 0.0)
  assert_true(var_psi > 0.0)
}

///|
test "ic_chi1_rotamers_leu" {
  let rotamers = @src.ic_chi1_rotamers("LEU")
  assert_true(rotamers.length() > 0)
}

///|
test "ic_chi1_rotamers_gly" {
  let rotamers = @src.ic_chi1_rotamers("GLY")
  // GLY has one default rotamer
  assert_eq(rotamers.length(), 1)
}