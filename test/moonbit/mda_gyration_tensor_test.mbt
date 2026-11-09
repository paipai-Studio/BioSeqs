///|
/// Tests for mda_gyration_tensor.mbt (gyration tensor & shape descriptors).

///|
test "mda_gyration_tensor_empty" {
  let gt = @src.mda_gyration_tensor([])
  assert_true(@moonbit.assert_num_close(gt.rg, 0.0, 1.0e-12))
  assert_eq(gt.eigenvalues.length(), 3)
  assert_true(@moonbit.assert_num_close(gt.eigenvalues[0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(gt.asphericity, 0.0, 1.0e-12))
}

///|
test "mda_gyration_tensor_single_atom" {
  let gt = @src.mda_gyration_tensor([[2.0, -1.0, 0.5]])
  // A single atom has zero deviations: tensor and eigenvalues are 0.
  let t = gt.tensor()
  assert_true(@moonbit.assert_num_close(t[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[2][2], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(gt.rg, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(gt.asphericity, 0.0, 1.0e-12))
}

///|
test "mda_gyration_tensor_isotropic_cube" {
  // Vertices of a cube of side 2 centered at the origin: S = I.
  // p2_val = 0 exactly -> degenerate branch -> eigenvalues [1,1,1].
  let coords : Array[Array[Double]] = Array(capacity=8)
  for x = -1; x <= 1; x = x + 2 {
    for y = -1; y <= 1; y = y + 2 {
      for z = -1; z <= 1; z = z + 2 {
        coords.push([x.to_double(), y.to_double(), z.to_double()])
      }
    }
  }
  let gt = @src.mda_gyration_tensor(coords)
  let t = gt.tensor()
  assert_true(@moonbit.assert_num_close(t[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[1][1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[2][2], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[0][1], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(gt.eigenvalues[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.eigenvalues[1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.eigenvalues[2], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.rg, 3.0.sqrt(), 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.asphericity, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.acylindricity, 0.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(gt.relative_shape_anisotropy, 0.0, 1.0e-9),
  )
}

///|
test "mda_gyration_tensor_planar_anisotropic" {
  // 4 atoms on the line x=0 with <y^2> = 2, <z^2> = 1, <yz> = 0:
  // S = diag(0, 2, 1). One eigenvalue equals the mean (1), so the
  // characteristic-polynomial solver is exact here.
  let s2 = 2.0.sqrt()
  let coords = [
    [0.0, s2, 1.0],
    [0.0, s2, -1.0],
    [0.0, -s2, 1.0],
    [0.0, -s2, -1.0],
  ]
  let gt = @src.mda_gyration_tensor(coords)
  let t = gt.tensor()
  assert_true(@moonbit.assert_num_close(t[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[1][1], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[2][2], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[0][1], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[1][2], 0.0, 1.0e-12))
  // Eigenvalues sorted descending: [2, 1, 0].
  assert_true(@moonbit.assert_num_close(gt.eigenvalues[0], 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.eigenvalues[1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.eigenvalues[2], 0.0, 1.0e-9))
  // rg = sqrt(trace) = sqrt(3).
  assert_true(@moonbit.assert_num_close(gt.rg, 3.0.sqrt(), 1.0e-9))
  // b = 2 - (1+0)/2 = 1.5; c = 1 - 0 = 1; kappa^2 = 1.5*(2.25+0.75)/9 = 0.5.
  assert_true(@moonbit.assert_num_close(gt.asphericity, 1.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.acylindricity, 1.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(gt.relative_shape_anisotropy, 0.5, 1.0e-9),
  )
}

///|
test "mda_gyration_tensor_mass_weighting" {
  // Atoms (0,0,0) and (2,0,0) with masses 1 and 3: COM at x = 1.5.
  // sxx = (1*2.25 + 3*0.25)/4 = 0.75; all other components 0.
  let gt = @src.mda_gyration_tensor([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], masses=[
    1.0, 3.0,
  ])
  let t = gt.tensor()
  assert_true(@moonbit.assert_num_close(t[0][0], 0.75, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[1][1], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[0][1], 0.0, 1.0e-12))
  // The sum of eigenvalues equals the trace.
  let sum_ev = gt.eigenvalues[0] + gt.eigenvalues[1] + gt.eigenvalues[2]
  assert_true(@moonbit.assert_num_close(sum_ev, 0.75, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.rg, 0.75.sqrt(), 1.0e-9))
}

///|
test "mda_gyration_tensor_offdiagonal_symmetric" {
  // Two atoms at +(1,2,3) and -(1,2,3): tensor [[1,2,3],[2,4,6],[3,6,9]].
  let gt = @src.mda_gyration_tensor([[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]])
  let t = gt.tensor()
  assert_true(@moonbit.assert_num_close(t[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[1][1], 4.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[2][2], 9.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[0][1], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[0][2], 3.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[1][2], 6.0, 1.0e-12))
  // Symmetry.
  assert_true(@moonbit.assert_num_close(t[1][0], t[0][1], 1.0e-12))
  assert_true(@moonbit.assert_num_close(t[2][0], t[0][2], 1.0e-12))
  // Eigenvalue sum always equals the trace (= 14).
  let sum_ev = gt.eigenvalues[0] + gt.eigenvalues[1] + gt.eigenvalues[2]
  assert_true(@moonbit.assert_num_close(sum_ev, 14.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(gt.rg, 14.0.sqrt(), 1.0e-9))
}

///|
test "mda_shape_descriptors_match_tensor" {
  let s2 = 2.0.sqrt()
  let coords = [
    [0.0, s2, 1.0],
    [0.0, s2, -1.0],
    [0.0, -s2, 1.0],
    [0.0, -s2, -1.0],
  ]
  let (b, c, kappa2) = @src.mda_shape_descriptors(coords)
  assert_true(@moonbit.assert_num_close(b, 1.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(c, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(kappa2, 0.5, 1.0e-9))
}
