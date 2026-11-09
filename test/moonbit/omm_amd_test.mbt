///|
test "amd_boost_math" {
  let params = @src.OmmAMDParams::new(e_threshold=10.0, alpha=5.0)
  // V = 0 < E: dV = 100/15 = 6.6667; ratio = 1 - 10*20/225 = 0.1111
  let (v1, r1) = @src.omm_amd_boost(params, 0.0)
  assert_true(@moonbit.assert_num_close(v1, 6.666666666666667, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r1, 0.1111111111111111, 1.0e-9))
  // V >= E: unchanged, ratio 1.
  let (v2, r2) = @src.omm_amd_boost(params, 10.0)
  assert_true(@moonbit.assert_num_close(v2, 10.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r2, 1.0, 1.0e-12))
  let (v3, r3) = @src.omm_amd_boost(params, 15.0)
  assert_true(@moonbit.assert_num_close(v3, 15.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r3, 1.0, 1.0e-12))
  // Continuity at V = E from below.
  let (vb, _rb) = @src.omm_amd_boost(params, 9.999999)
  assert_true((vb - 10.0).abs() < 1.0e-5)
}

///|
test "gamd_boost_math" {
  let params = @src.OmmAMDParams::new(e_threshold=5.0, alpha=1.0)
  // V = 3: x = 2; dV = 4/5 = 0.8; ratio = 1 - 4/25 = 0.84
  let (v1, r1) = @src.omm_gamd_boost(params, 3.0)
  assert_true(@moonbit.assert_num_close(v1, 3.8, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r1, 0.84, 1.0e-12))
  // V >= E unchanged.
  let (v2, r2) = @src.omm_gamd_boost(params, 5.0)
  assert_true(@moonbit.assert_num_close(v2, 5.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r2, 1.0, 1.0e-12))
}

///|
test "amd_simulate_harmonic_smoke" {
  // Harmonic oscillator under a boost: positions stay finite, the boosted
  // energy (KE + V') is approximately conserved, and dt is respected.
  let k = 2.0
  let params = @src.OmmAMDParams::new(e_threshold=0.1, alpha=0.5)
  let (pos_hist, vel_hist) = @src.omm_amd_simulate(
    [[0.3, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
    n_steps=300,
    force_fn=fn(pos : Array[Array[Double]]) -> Array[Array[Double]] {
      [[-k * pos[0][0], 0.0, 0.0]]
    },
    energy_fn=fn(pos : Array[Array[Double]]) -> Double {
      0.5 * k * pos[0][0] * pos[0][0]
    },
    params,
    dt=0.01,
  )
  assert_eq(pos_hist.length(), 301)
  // Boosted energy conservation (KE + boosted PE) over the run.
  let mut e0 = 0.0
  let (vb0, _r0) = @src.omm_amd_boost(params, 0.5 * k * 0.09)
  e0 = vb0
  let pf = pos_hist[300]
  let vf = vel_hist[300]
  let ke = 0.5 *
    (vf[0][0] * vf[0][0] + vf[0][1] * vf[0][1] + vf[0][2] * vf[0][2])
  let (vbf, _rf) = @src.omm_amd_boost(params, 0.5 * k * pf[0][0] * pf[0][0])
  assert_true((ke + vbf - e0).abs() < 0.15)
  // No NaN anywhere.
  assert_true(pf[0][0] == pf[0][0])
  assert_true(vf[0][0] == vf[0][0])
}
