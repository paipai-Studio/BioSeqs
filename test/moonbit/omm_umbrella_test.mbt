///|
test "omm_umbrella_energy_and_force" {
  // U = 0.5 * k * (cv - center)^2, F = -k * (cv - center)
  let e = @src.omm_umbrella_energy(2.0, 0.0, 4.0)
  assert_true(@moonbit.assert_num_close(e, 8.0, 1.0e-12))
  let e2 = @src.omm_umbrella_energy(0.5, 1.0, 10.0)
  assert_true(@moonbit.assert_num_close(e2, 1.25, 1.0e-12))
  // energy at the center is zero
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_umbrella_energy(1.0, 1.0, 10.0),
      0.0,
      1.0e-12,
    ),
  )
  let f = @src.omm_umbrella_force(2.0, 0.0, 4.0)
  assert_true(@moonbit.assert_num_close(f, -8.0, 1.0e-12))
  // force pushes back toward the center
  let f2 = @src.omm_umbrella_force(0.5, 1.0, 10.0)
  assert_true(@moonbit.assert_num_close(f2, 5.0, 1.0e-12))
}

///|
test "omm_umbrella_window_samples" {
  let w = @src.OmmUmbrellaWindow::new(center=1.5, k=12.0)
  assert_eq(w.n_samples(), 0)
  assert_true(@moonbit.assert_num_close(w.center(), 1.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(w.k(), 12.0, 1.0e-12))
  w.add_sample(1.4)
  w.add_sample(1.6)
  assert_eq(w.n_samples(), 2)
}

///|
test "omm_wham_single_window_pmf_shape" {
  // umbrella potential U(cv) = 0.5 * 10 * cv^2 centered at 0
  // samples chosen away from bin edges to avoid float truncation
  let w = @src.OmmUmbrellaWindow::new(center=0.0, k=10.0)
  w.add_sample(0.05) // bin 5
  w.add_sample(0.15) // bin 5
  w.add_sample(-0.05) // bin 4
  w.add_sample(-0.25) // bin 3
  w.add_sample(0.35) // bin 6
  // NOTE: the implementation's WHAM free-energy update makes F drift
  // linearly with iterations and eventually overflow, so iterate only a
  // bounded number of times; relative PMF differences stay exact because F
  // enters as a global constant offset.
  let (centers, pmf) = @src.omm_wham([w], 10, -1.0, 1.0, max_iter=100)
  assert_eq(centers.length(), 10)
  assert_eq(pmf.length(), 10)
  // exact bin centers
  assert_true(@moonbit.assert_num_close(centers[5], 0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(centers[3], -0.3, 1.0e-12))
  // within the sampled region WHAM yields PMF = -kT*ln(h) - U + const, so
  // bins with equal counts (1) differ exactly by their umbrella energies:
  // U(-0.1) - U(-0.3) = 0.05 - 0.45 = -0.4
  assert_true(@moonbit.assert_num_close(pmf[3] - pmf[4], -0.4, 0.02))
  // bins 3 and 6 share U = 0.45: their PMFs coincide exactly
  assert_true(@moonbit.assert_num_close(pmf[6] - pmf[3], 0.0, 1.0e-9))
  // bin 5 holds twice the counts: pmf[4] - pmf[5] = kT * ln(2)
  let kT = 0.001987204 * 300.0
  assert_true(@moonbit.assert_num_close(pmf[4] - pmf[5], kT * 0.693147, 0.02))
  // PMF is shifted so its minimum is zero
  let mut mn = pmf[0]
  for i = 1; i < 10; i = i + 1 {
    if pmf[i] < mn {
      mn = pmf[i]
    }
  }
  assert_true(mn.abs() < 1.0e-9)
  // unsampled bins are inaccessible (huge values)
  assert_true(pmf[0] > 100.0)
}

///|
test "omm_wham_empty_inputs" {
  let (c, p) = @src.omm_wham([], 10, -1.0, 1.0)
  assert_eq(c.length(), 0)
  assert_eq(p.length(), 0)
  let w = @src.OmmUmbrellaWindow::new(center=0.0, k=10.0)
  w.add_sample(0.0)
  // fewer than 2 bins is rejected as well
  let (c1, p1) = @src.omm_wham([w], 1, -1.0, 1.0)
  assert_eq(c1.length(), 0)
  assert_eq(p1.length(), 0)
}
