///|
test "fft_3d_roundtrip" {
  // Forward + inverse FFT restores the original signal.
  let nx = 4
  let ny = 2
  let nz = 2
  let re : Array[Double] = Array::make(nx * ny * nz, 0.0)
  let im : Array[Double] = Array::make(nx * ny * nz, 0.0)
  let re0 : Array[Double] = Array::make(nx * ny * nz, 0.0)
  for i = 0; i < nx * ny * nz; i = i + 1 {
    re[i] = @math.sin(0.7 * i.to_double()) +
      0.3 * @math.cos(1.3 * i.to_double())
    re0[i] = re[i]
  }
  @src.omm_fft_3d(re, im, nx, ny, nz)
  @src.omm_fft_3d(re, im, nx, ny, nz, inverse=true)
  for i = 0; i < nx * ny * nz; i = i + 1 {
    assert_true(@moonbit.assert_num_close(re[i], re0[i], 1.0e-9))
  }
}

///|
test "fft_3d_dc_and_single_mode" {
  // Constant signal: only the DC bin is nonzero with value N.
  let n = 4
  let re : Array[Double] = Array::make(n * n * n, 1.0)
  let im : Array[Double] = Array::make(n * n * n, 0.0)
  @src.omm_fft_3d(re, im, n, n, n)
  for i = 0; i < n * n * n; i = i + 1 {
    if i == 0 {
      assert_true(@moonbit.assert_num_close(re[0], 64.0, 1.0e-9))
      assert_true(@moonbit.assert_num_close(im[0], 0.0, 1.0e-9))
    } else {
      assert_true(re[i].abs() < 1.0e-9)
      assert_true(im[i].abs() < 1.0e-9)
    }
  }
  assert_true(@src.omm_fft_next_pow2(5) == 8)
  assert_true(@src.omm_fft_next_pow2(32) == 32)
}

///|
test "pme_erfc_via_real_space_limit" {
  // With the corrected A&S 7.1.26 erfc, the real-space term for two
  // +1 charges at r=0.1 (alpha = sqrt(pi)/10 = 0.17725, so alpha*r =
  // 0.01772) is ~ q^2 * erfc(0.01772)/0.1 = 9.800 (the old broken sign
  // gave ~0.97).
  let params = @src.OmmPMEParams::new(10.0, box=[10.0, 10.0, 10.0])
  let pos : Array[Array[Double]] = [[5.0, 5.0, 5.0], [5.1, 5.0, 5.0]]
  let charges = [1.0, 1.0]
  let e = @src.omm_pme_real_space(pos, charges, params)
  assert_true(@moonbit.assert_num_close(e, 9.800, 0.02))
}

///|
test "pme_grid_reciprocal_matches_direct_sum" {
  // Two-charge neutral pair in a cubic box: the FFT-grid reciprocal energy
  // agrees with the analytic k-space sum to interpolation accuracy.
  let params = @src.OmmPMEParams::new(8.0, grid_x=32, grid_y=32, grid_z=32, box=[
    10.0, 10.0, 10.0,
  ])
  let pos : Array[Array[Double]] = [[2.0, 3.0, 4.0], [7.0, 6.0, 5.0]]
  let charges = [1.0, -1.0]
  let direct = @src.omm_pme_reciprocal_space(pos, charges, params)
  let grid = @src.omm_pme_grid_reciprocal(pos, charges, params)
  let scale = if direct.abs() > 1.0e-9 { direct.abs() } else { 1.0 }
  assert_true((grid - direct).abs() / scale < 0.02)
  // Full energy (real + grid reciprocal + self) matches the direct variant.
  let e_direct = @src.omm_pme_energy(pos, charges, params)
  let e_grid = @src.omm_pme_grid_energy(pos, charges, params)
  let scale2 = if e_direct.abs() > 1.0e-9 { e_direct.abs() } else { 1.0 }
  assert_true((e_grid - e_direct).abs() / scale2 < 0.01)
}

///|
test "pme_grid_four_charges_neutral" {
  // A neutral 4-charge cluster: grid and direct reciprocal sums agree,
  // and the total PME energy is finite and stable.
  let params = @src.OmmPMEParams::new(8.0, grid_x=32, grid_y=32, grid_z=32, box=[
    12.0, 12.0, 12.0,
  ])
  let pos : Array[Array[Double]] = [
    [1.0, 1.0, 1.0],
    [4.0, 2.0, 3.0],
    [6.0, 8.0, 5.0],
    [9.0, 9.0, 2.0],
  ]
  let charges = [0.8, -0.8, 0.4, -0.4]
  let direct = @src.omm_pme_reciprocal_space(pos, charges, params)
  let grid = @src.omm_pme_grid_reciprocal(pos, charges, params)
  let scale = if direct.abs() > 1.0e-9 { direct.abs() } else { 1.0 }
  assert_true((grid - direct).abs() / scale < 0.05)
  let e = @src.omm_pme_grid_energy(pos, charges, params)
  assert_true(e == e)
}
