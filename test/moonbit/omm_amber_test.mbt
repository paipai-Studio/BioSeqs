///|
/// OmmAmberParams cannot be constructed outside the src package, so start from
/// an empty parse result and inject bonded terms via the mutable arrays.
fn omm_amber_make_params(
  bonds : Array[(Int32, Int32, Double, Double)],
  angles : Array[(Int32, Int32, Int32, Double, Double)],
  dihedrals : Array[(Int32, Int32, Int32, Int32, Double, Double, Int32)],
) -> @src.OmmAmberParams {
  let params = @src.omm_parse_prmtop("")
  for b in bonds {
    params.bonds.push(b)
  }
  for a in angles {
    params.angles.push(a)
  }
  for d in dihedrals {
    params.dihedrals.push(d)
  }
  params
}

///|
test "omm amber simple bond ff builder" {
  let bonds = @src.omm_simple_bond_ff(3, [(0, 1), (1, 2)], k=120.0, r0=1.4)
  assert_eq(bonds.length(), 2)
  assert_eq(bonds[0].0, 0)
  assert_eq(bonds[0].1, 1)
  assert_eq(bonds[0].2, 120.0)
  assert_eq(bonds[0].3, 1.4)
  assert_eq(bonds[1].0, 1)
  assert_eq(bonds[1].1, 2)
  assert_eq(bonds[1].2, 120.0)
  assert_eq(bonds[1].3, 1.4)
}

///|
test "omm amber simple bond ff default parameters" {
  let bonds = @src.omm_simple_bond_ff(2, [(0, 1)])
  assert_eq(bonds.length(), 1)
  assert_eq(bonds[0].2, 100.0)
  assert_eq(bonds[0].3, 1.0)
}

///|
test "omm amber energy zero for empty prmtop parse" {
  let params = @src.omm_parse_prmtop("%FLAG POINTERS\n%FORMAT(10I8)\n")
  assert_eq(params.n_atoms(), 0)
  let e = @src.omm_amber_energy([[0.0, 0.0, 0.0]], params)
  assert_eq(e, 0.0)
}

///|
test "omm amber bond energy harmonic at stretched bond" {
  // E = 0.5 * k * (r - r0)^2 with k=100, r0=1, r=2 -> 50
  let params = omm_amber_make_params([(0, 1, 100.0, 1.0)], [], [])
  let e = @src.omm_amber_energy([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], params)
  assert_true(@moonbit.assert_num_close(e, 50.0, 1.0e-10))
}

///|
test "omm amber bond energy zero at equilibrium" {
  let params = omm_amber_make_params([(0, 1, 100.0, 1.5)], [], [])
  let e = @src.omm_amber_energy([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], params)
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-12))
}

///|
test "omm amber angle energy right angle vs linear reference" {
  // atoms: j at origin, i and k at right angle -> theta = pi/2, theta0 = 0
  // E = 0.5 * k * (pi/2)^2
  let params = omm_amber_make_params([], [(0, 1, 2, 1.0, 0.0)], [])
  let positions = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let e = @src.omm_amber_energy(positions, params)
  assert_true(
    @moonbit.assert_num_close(
      e,
      0.5 * 3.141592653589793 * 3.141592653589793 / 4.0,
      1.0e-10,
    ),
  )
}

///|
test "omm amber angle energy zero at equilibrium angle" {
  // theta = pi/2 = theta0 -> zero energy
  let theta0 = 3.141592653589793 / 2.0
  let params = omm_amber_make_params([], [(0, 1, 2, 5.0, theta0)], [])
  let positions = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let e = @src.omm_amber_energy(positions, params)
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-12))
}

///|
test "omm amber dihedral energy trans conformation" {
  // trans dihedral phi = pi: E = k_phi * (1 + cos(n*phi - phi0))
  // k_phi=2, n=1, phi0=0 -> 2*(1 + cos(pi)) = 0
  let params = omm_amber_make_params([], [], [(0, 1, 2, 3, 2.0, 0.0, 1)])
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [2.0, 1.0, 0.0],
  ]
  let e = @src.omm_amber_energy(positions, params)
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-10))
}

///|
test "omm amber dihedral energy cis conformation" {
  // cis dihedral phi = 0: k_phi=2, n=1, phi0=0 -> 2*(1+1) = 4
  let params = omm_amber_make_params([], [], [(0, 1, 2, 3, 2.0, 0.0, 1)])
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let e = @src.omm_amber_energy(positions, params)
  assert_true(@moonbit.assert_num_close(e, 4.0, 1.0e-10))
}

///|
test "omm amber dihedral energy periodicity two" {
  // n=2, phi=0, phi0=0: E = k*(1+cos(0)) = 2k
  let params = omm_amber_make_params([], [], [(0, 1, 2, 3, 1.5, 0.0, 2)])
  let positions = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let e = @src.omm_amber_energy(positions, params)
  assert_true(@moonbit.assert_num_close(e, 3.0, 1.0e-10))
}

///|
test "omm amber energy sums bond and angle contributions" {
  // bond: r=2, r0=1, k=100 -> 50 ; angle: pi/2 vs theta0=0, k=2 -> (pi/2)^2
  let params = omm_amber_make_params(
    [(0, 1, 100.0, 1.0)],
    [(0, 1, 2, 2.0, 0.0)],
    [],
  )
  let positions = [
    [2.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [9.0, 9.0, 9.0],
  ]
  let e = @src.omm_amber_energy(positions, params)
  // angle term: 0.5 * k_theta * (pi/2)^2 = pi^2/4 with k_theta=2
  let expected = 50.0 + 3.141592653589793 * 3.141592653589793 / 4.0
  assert_true(@moonbit.assert_num_close(e, expected, 1.0e-10))
}
