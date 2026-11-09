///|
test "system_universe_roundtrip" {
  let atoms : Array[@src.UnivAtom] = [
    {
      name: "C1",
      resname: "MOL",
      resid: 1,
      chain: "A",
      segid: "SYS",
      element: "C",
      atom_type: "CT",
      charge: -0.1,
      mass: 12.011,
      is_protein: false,
      is_water: false,
      is_backbone: false,
    },
    {
      name: "H1",
      resname: "MOL",
      resid: 1,
      chain: "A",
      segid: "SYS",
      element: "H",
      atom_type: "HC",
      charge: 0.1,
      mass: 1.008,
      is_protein: false,
      is_water: false,
      is_backbone: false,
    },
  ]
  let frame : Array[Array[Double]] = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
  let u0 = @src.Universe::new(atoms, @src.Trajectory::new([frame]))
  // Universe -> OmmSystem
  let sys = @src.omm_system_from_universe(u0, sigma=3.0, epsilon=0.1)
  assert_eq(sys.num_particles(), 2)
  let charges = sys.get_charges()
  assert_true(@moonbit.assert_num_close(charges[0], -0.1, 1.0e-12))
  let masses = sys.get_masses()
  assert_true(@moonbit.assert_num_close(masses[0], 12.011, 1.0e-12))
  assert_true(@moonbit.assert_num_close(masses[1], 1.008, 1.0e-12))
  sys.set_periodic_box(10.0, 10.0, 10.0)
  assert_true(sys.is_periodic())
  // OmmSystem -> Universe (elements guessed from mass)
  let u1 = @src.omm_system_to_universe(sys, positions=Some(frame))
  assert_eq(u1.n_atoms(), 2)
  let info = u1.atom_info()
  assert_eq(info[0].element, "C")
  assert_eq(info[1].element, "H")
  assert_true(@moonbit.assert_num_close(info[0].charge, -0.1, 1.0e-12))
  assert_eq(u1.frame_index(), 0)
  let back_pos = u1.positions_at(0)
  assert_true(@moonbit.assert_num_close(back_pos[1][0], 1.5, 1.0e-12))
}

///|
test "prmtop_parse_full_and_bridges" {
  let prmtop = "%VERSION  VERSION STAMP = V0001.000  DATE = 11/07/26\n%FLAG POINTERS\n%FORMAT(10I8)\n       2       2       1       0       0       0       0       0       0       1\n%FLAG ATOM_NAME\n%FORMAT(20a4)\nC   H   \n%FLAG CHARGE\n%FORMAT(5E16.8)\n  -1.82000000E+00   1.82000000E+00\n%FLAG MASS\n%FORMAT(5E16.8)\n  1.20110000E+01   1.00800000E+00\n%FLAG ATOM_TYPE_INDEX\n%FORMAT(10I8)\n       1       2\n%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n       6       1\n%FLAG AMBER_ATOM_TYPE\n%FORMAT(20a4)\nCT  HC  \n%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\nMOL \n%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n       1\n%FLAG BOND_FORCE_CONSTANT\n%FORMAT(5E16.8)\n  3.00000000E+02\n%FLAG BOND_EQUIL_VALUE\n%FORMAT(5E16.8)\n  1.50000000E+00\n%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n       0       3       1\n%FLAG LENNARD_JONES_ACOEF\n%FORMAT(5E16.8)\n  1.00000000E+03   8.00000000E+02   6.00000000E+02\n%FLAG LENNARD_JONES_BCOEF\n%FORMAT(5E16.8)\n  5.00000000E+00   4.00000000E+00   3.00000000E+00\n"
  let top = @src.omm_parse_prmtop_full(prmtop)
  assert_eq(top.n_atoms, 2)
  assert_eq(top.n_types, 2)
  assert_eq(top.atom_names.length(), 2)
  assert_eq(top.atom_names[0], "C")
  assert_eq(top.atom_names[1], "H")
  assert_eq(top.atomic_numbers[0], 6)
  assert_eq(top.bonds.length(), 1)
  // bond indices are 0-based: (0, 1, param 0)
  let b = top.bonds[0]
  assert_eq(b.0, 0)
  assert_eq(b.1, 1)
  assert_eq(b.2, 0)
  // charges: internal units / 18.2223 -> e (|q| = 1.82 / 18.2223 = 0.09989)
  assert_true(
    @moonbit.assert_num_close(top.charges[0].abs(), 1.82 / 18.2223, 1.0e-9),
  )
  // LJ diagonal (1,1): sigma = (1000/5)^(1/6), eps = 25/4000
  let lj = @src.omm_prmtop_lj(top, 1, 1)
  assert_true(
    @moonbit.assert_num_close(lj.0, @math.pow(200.0, 1.0 / 6.0), 1.0e-6),
  )
  assert_true(@moonbit.assert_num_close(lj.1, 0.00625, 1.0e-9))
  // Off-diagonal symmetric lookup (2,1) == (1,2)
  let lj2a = @src.omm_prmtop_lj(top, 1, 2)
  let lj2b = @src.omm_prmtop_lj(top, 2, 1)
  assert_true(@moonbit.assert_num_close(lj2a.0, lj2b.0, 1.0e-12))
  // prmtop -> Universe (with inpcrd coordinates)
  let inpcrd = "test system\n        2\n     0.0000000    0.0000000    0.0000000\n     1.5000000    0.0000000    0.0000000\n"
  let u = @src.omm_universe_from_prmtop(prmtop, inpcrd~)
  assert_eq(u.n_atoms(), 2)
  assert_eq(u.n_frames(), 1)
  let info = u.atom_info()
  assert_eq(info[0].resname, "MOL")
  assert_eq(info[0].resid, 1)
  assert_eq(info[0].element, "C")
  assert_eq(info[1].element, "H")
  assert_true(info[0].charge < 0.0)
  let p = u.positions_at(0)
  assert_true(@moonbit.assert_num_close(p[1][0], 1.5, 1.0e-6))
  // prmtop -> OmmSystem
  let sys = @src.omm_system_from_prmtop(prmtop)
  assert_eq(sys.num_particles(), 2)
  assert_eq(sys.num_constraints(), 0)
  let e = sys.potential_energy([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
  // Bond at equilibrium contributes 0; LJ is repulsive here. The pair
  // epsilon is the geometric mean of the per-type epsilons
  // (sqrt(0.00625 * 0.00375) = 0.004841):
  // 4*eps*((s/r)^12-(s/r)^6) = 4*0.004841*(308.294 - 17.558) = 5.6301
  assert_true(@moonbit.assert_num_close(e, 5.6301, 1.0e-3))
  let e2 = sys.potential_energy([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
  // LJ is repulsive at r=2.0 < sigma=2.418 (4*eps*((s/r)^12-(s/r)^6) =
  // 0.1286), plus bond strain 0.5*300*0.25 = 37.5
  assert_true(e2 > 10.0)
  // Box variant
  let sys2 = @src.omm_system_from_prmtop_box(prmtop, box=[20.0, 20.0, 20.0])
  assert_true(sys2.is_periodic())
}

///|
test "element_from_number_and_mass" {
  assert_eq(@src.omm_element_from_number(6), "C")
  assert_eq(@src.omm_element_from_number(8), "O")
  assert_eq(@src.omm_element_from_number(11), "Na")
  assert_eq(@src.omm_element_from_number(16), "S")
  assert_eq(@src.omm_element_from_number(999), "X")
  assert_eq(@src.omm_element_from_mass(1.008), "H")
  assert_eq(@src.omm_element_from_mass(12.011), "C")
  assert_eq(@src.omm_element_from_mass(15.999), "O")
  assert_eq(@src.omm_element_from_mass(22.990), "Na")
  assert_eq(@src.omm_element_from_mass(3.14159), "X")
}

///|
test "prmtop_missing_flags_defaults" {
  // A minimal file with only POINTERS parses to empty arrays.
  let tiny = "%FLAG POINTERS\n%FORMAT(10I8)\n       3       1\n"
  let top = @src.omm_parse_prmtop_full(tiny)
  assert_eq(top.n_atoms, 3)
  assert_eq(top.bonds.length(), 0)
  let u = @src.omm_universe_from_prmtop(tiny)
  assert_eq(u.n_atoms(), 3)
  assert_eq(u.n_frames(), 0)
}
