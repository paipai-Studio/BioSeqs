///|
test "drude_polarizability_and_spring" {
  // alpha = q^2 * 332.06371 / k.
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_drude_polarizability(1.0, 1000.0),
      0.33206371,
      1.0e-9,
    ),
  )
  assert_true(@src.omm_drude_polarizability(1.0, 0.0) == 0.0)
  // Spring forces: equal and opposite, proportional to separation.
  let pairs : Array[@src.OmmDrudePair] = [
    { core: 0, particle: 1, k_drude: 500.0, },
  ]
  let pos : Array[Array[Double]] = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]
  let f = @src.omm_drude_spring_forces(pos, pairs)
  // The spring pulls each particle toward the other: the core (at the
  // origin) feels +50 along +x, the Drude particle feels -50.
  assert_true(@moonbit.assert_num_close(f[0][0], 50.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f[1][0], -50.0, 1.0e-9))
}

///|
test "drude_oscillator_relaxes_to_boltzmann" {
  // Free core, Drude particle attached by a spring: the equilibrium
  // separation follows sqrt(3 kT / k) (3 relative DOF).
  let k = 500.0
  let pairs : Array[@src.OmmDrudePair] = [{ core: 0, particle: 1, k_drude: k, }]
  let masses = [12.0, 0.4]
  @src.md_seed_rng(2024)
  let (pos_hist, _vel_hist) = @src.omm_drude_langevin_simulate(
    [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    force_fn=fn(_pos : Array[Array[Double]]) -> Array[Array[Double]] {
      @src.md_new_forces(2)
    },
    masses~,
    pairs~,
    n_steps=4000,
    dt=0.0005,
    temperature=300.0,
    friction_core=1.0,
    friction_drude=50.0,
    n_inner=4,
  )
  // Mean separation over the last 1000 steps.
  let mut sum_r = 0.0
  for s = 3000; s < 4000; s = s + 1 {
    let p = pos_hist[s]
    let dx = p[1][0] - p[0][0]
    let dy = p[1][1] - p[0][1]
    let dz = p[1][2] - p[0][2]
    sum_r = sum_r + @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
  }
  let mean_r = sum_r / 1000.0
  let expect = @math.pow(3.0 * 0.001987204 * 300.0 / k, 0.5)
  // Wide statistical band (single short trajectory with a strong OU
  // coupling on a light particle).
  assert_true(mean_r > 0.3 * expect && mean_r < 2.5 * expect)
}

///|
test "drude_nve_reduces_to_verlet" {
  // No Drude pairs + zero friction: plain velocity Verlet on free flight.
  let pairs : Array[@src.OmmDrudePair] = []
  let (pos_hist, vel_hist) = @src.omm_drude_langevin_simulate(
    [[0.0, 0.0, 0.0]],
    [[2.0, 0.0, 0.0]],
    force_fn=fn(_pos : Array[Array[Double]]) -> Array[Array[Double]] {
      @src.md_new_forces(1)
    },
    masses=[1.0],
    pairs~,
    n_steps=10,
    dt=0.1,
    temperature=0.0,
    friction_core=0.0,
    n_inner=1,
  )
  assert_true(@moonbit.assert_num_close(pos_hist[10][0][0], 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vel_hist[10][0][0], 2.0, 1.0e-9))
}
