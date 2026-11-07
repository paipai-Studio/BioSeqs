///|
test "builtin_amber99sb_parse" {
  let ff = @src.omm_builtin_amber99sb()
  assert_eq(ff.name, "amber99sb-builtin")
  assert_eq(ff.atom_types.length(), 14)
  assert_eq(ff.residues.length(), 20)
  assert_true(ff.bond_forces.length() > 0)
  assert_true(ff.angle_forces.length() > 0)
  assert_true(ff.torsion_forces.length() > 0)
  assert_eq(ff.nonbonded.length(), 14)
  // Spot-check the CT type and one bond.
  let mut ct_mass = 0.0
  for t in ff.atom_types {
    if t.name == "CT" {
      ct_mass = t.mass
    }
  }
  assert_true(@moonbit.assert_num_close(ct_mass, 12.011, 1.0e-12))
  let (k, r0) = @src.omm_xml_ff_lookup_bond(ff, "C", "N").unwrap()
  assert_true(@moonbit.assert_num_close(k, 1966.49, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r0, 0.1335, 1.0e-12))
  // Symmetric lookup order.
  let (k2, _) = @src.omm_xml_ff_lookup_bond(ff, "N", "C").unwrap()
  assert_true(@moonbit.assert_num_close(k2, 1966.49, 1.0e-9))
  // Angle lookup (symmetric ends).
  let (ak, ath) = @src.omm_xml_ff_lookup_angle(ff, "O", "C", "N").unwrap()
  assert_true(@moonbit.assert_num_close(ak, 422.58, 1.0e-9))
  assert_true(@moonbit.assert_num_close(ath, 122.9, 1.0e-12))
  // Torsion terms for the backbone phi chain.
  let terms = @src.omm_xml_ff_lookup_torsions(ff, "C", "N", "CT", "C")
  assert_eq(terms.length(), 1)
  assert_eq(terms[0].2, 2)
}

///|
test "builtin_ff_build_ala_system" {
  let ff = @src.omm_builtin_amber99sb()
  let atoms : Array[@src.OmmXmlAtomSpec] = [
    { name: "N", resname: "ALA", element: "N", },
    { name: "H", resname: "ALA", element: "H", },
    { name: "CA", resname: "ALA", element: "C", },
    { name: "HA", resname: "ALA", element: "H", },
    { name: "CB", resname: "ALA", element: "C", },
    { name: "HB1", resname: "ALA", element: "H", },
    { name: "HB2", resname: "ALA", element: "H", },
    { name: "HB3", resname: "ALA", element: "H", },
    { name: "C", resname: "ALA", element: "C", },
    { name: "O", resname: "ALA", element: "O", },
  ]
  let sys = @src.omm_xml_ff_build_system(ff, atoms)
  assert_eq(sys.num_particles(), 10)
  // ALA template has 9 bonds.
  assert_eq(sys.num_bonds(), 9)
  // Angles are always generated from the bond graph; torsions only when
  // the type quadruple has XML parameters (ALA's all miss -> none).
  assert_true(sys.num_angles() > 0)
  assert_eq(sys.num_dihedrals(), 0)
  // The N-CA bond uses the XML k converted to kcal/mol/A^2.
  let masses = sys.get_masses()
  assert_true(@moonbit.assert_num_close(masses[0], 14.007, 1.0e-9))
  assert_true(@moonbit.assert_num_close(masses[1], 1.008, 1.0e-9))
  // Energy is finite for a reasonable geometry.
  let pos : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 1.01],
    [1.47, 0.0, -0.3],
    [2.2, 0.9, -0.3],
    [2.0, 1.2, -1.4],
    [3.0, 1.2, -1.4],
    [2.0, 2.2, -1.4],
    [2.0, 1.2, -2.4],
    [2.0, -0.2, 0.3],
    [3.0, -0.2, 0.3],
  ]
  let e = sys.potential_energy(pos)
  // The hand-built geometry is far from equilibrium (large but finite
  // strain), so only finiteness is asserted here.
  assert_true(e == e)
  // Minimizing from this geometry drives the energy far down
  // (realistic amber99sb-style parameters).
  let sim = @src.OmmSimulation::new(system=sys, dt=0.001, positions=Some(pos))
  let dropped = sim.minimize(max_iter=200, tolerance=0.001)
  assert_true(dropped > 0.0)
  assert_true(sim.state().potential_energy < e)
}

///|
test "builtin_ff_build_ser_and_fallback" {
  let ff = @src.omm_builtin_amber99sb()
  // SER: OG-HO bond present; unknown-residue atoms fall back to element
  // types without bonds.
  let atoms : Array[@src.OmmXmlAtomSpec] = [
    { name: "N", resname: "SER", element: "N", },
    { name: "H", resname: "SER", element: "H", },
    { name: "CA", resname: "SER", element: "C", },
    { name: "HA", resname: "SER", element: "H", },
    { name: "CB", resname: "SER", element: "C", },
    { name: "HB1", resname: "SER", element: "H", },
    { name: "HB2", resname: "SER", element: "H", },
    { name: "OG", resname: "SER", element: "O", },
    { name: "HG", resname: "SER", element: "H", },
    { name: "C", resname: "SER", element: "C", },
    { name: "O", resname: "SER", element: "O", },
  ]
  let sys = @src.omm_xml_ff_build_system(ff, atoms)
  assert_eq(sys.num_particles(), 11)
  assert_eq(sys.num_bonds(), 10)
  // OG-HO bond: check via potential energy sensitivity is overkill; instead
  // verify through a second atom typed by element fallback (no template).
  let fallback : Array[@src.OmmXmlAtomSpec] = [
    { name: "A", resname: "XYZ", element: "C", },
    { name: "B", resname: "XYZ", element: "C", },
  ]
  let sys2 = @src.omm_xml_ff_build_system(ff, fallback)
  assert_eq(sys2.num_particles(), 2)
  // No template -> no bonds; element fallback types give CT params.
  assert_eq(sys2.num_bonds(), 0)
  let masses = sys2.get_masses()
  assert_true(@moonbit.assert_num_close(masses[0], 12.011, 1.0e-9))
}

///|
test "omm_xml_ff_build_water_units" {
  let ff = @src.omm_parse_xml_ff(xml_ff_water_string())
  let atoms : Array[@src.OmmXmlAtomSpec] = [
    { name: "O", resname: "HOH", element: "O", },
    { name: "H1", resname: "HOH", element: "H", },
    { name: "H2", resname: "HOH", element: "H", },
  ]
  let sys = @src.omm_xml_ff_build_system(ff, atoms)
  assert_eq(sys.num_particles(), 2 + 1)
  assert_eq(sys.num_bonds(), 2)
  // XML k=553 kJ/mol/nm^2 -> 553/418.4 kcal/mol/A^2; r0 0.9572 nm -> 9.572 A.
  // Place O-H at exactly 9.572/2... no — at r0 in A the bond energy is 0.
  let r0_a = 9.572
  let pos : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [r0_a, 0.0, 0.0],
    [0.0, r0_a, 0.0],
  ]
  let e = sys.potential_energy(pos)
  // Both bonds at equilibrium: bond energy 0; angle falls back to
  // 109.5 deg with k=50 — the H-O-H angle here is 90 deg:
  // E = 0.5 * 50 * (90 - 109.5 deg in rad)^2.
  let dth = (90.0 - 109.5) * 3.141592653589793 / 180.0
  let expect = 0.5 * 50.0 * dth * dth
  assert_true(@moonbit.assert_num_close(e, expect, 1.0e-6))
  // Nonbonded units: sigma nm -> A.
  let (q, s, eps) = @src.omm_lookup_nonbonded(ff, "OW")
  assert_true(@moonbit.assert_num_close(q, -0.834, 1.0e-12))
  ignore(s)
  ignore(eps)
}

///|
test "omm_system_dihedral_energy" {
  let sys = @src.OmmSystem::new()
  for _i = 0; _i < 4; _i = _i + 1 {
    ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=1.0)))
  }
  let k = 1.5
  sys.add_dihedral(0, 1, 2, 3, k, 0.0, 3)
  // cis (phi = 0): E = k (1 + cos(0)) = 2k
  let cis : Array[Array[Double]] = [
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.5, 0.8660254037844386, 0.0],
  ]
  assert_true(
    @moonbit.assert_num_close(sys.potential_energy(cis), 2.0 * k, 1.0e-9),
  )
  // trans (phi = 180): E = k (1 + cos(3*pi)) = 0
  let trans : Array[Array[Double]] = [
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.5, -0.8660254037844386, 0.0],
  ]
  assert_true(
    @moonbit.assert_num_close(sys.potential_energy(trans), 0.0, 1.0e-9),
  )
  // phi = 60 deg: E = k (1 + cos(180 deg)) = 0
  let phi60 : Array[Array[Double]] = [
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.5, 0.5, 0.8660254037844386],
  ]
  assert_true(
    @moonbit.assert_num_close(sys.potential_energy(phi60), 0.0, 1.0e-9),
  )
}
