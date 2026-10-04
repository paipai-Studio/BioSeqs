///|
/// Tests for conformation sampling from normal modes.

test "SamplingConfig default" {
  let cfg = @src.SamplingConfig::default()
  assert_true((cfg.kT - 1.0).abs() < 1.0e-10)
  assert_eq(cfg.num_modes, 0)
  assert_eq(cfg.seed, 42)
  assert_true((cfg.scale - 1.0).abs() < 1.0e-10)
}

test "calc_rmsd identical coords" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 1.0, 1.0),
    @src.Vector3::new(2.0, 0.0, 0.0),
  ]
  let rmsd = @src.calc_rmsd(coords, coords)
  assert_true(rmsd.abs() < 1.0e-10)
}

test "calc_rmsd shifted coords" {
  let coords1 = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let coords2 = [
    @src.Vector3::new(0.0, 1.0, 0.0),
    @src.Vector3::new(1.0, 1.0, 0.0),
  ]
  let rmsd = @src.calc_rmsd(coords1, coords2)
  assert_true((rmsd - 1.0).abs() < 1.0e-10)
}

test "calc_rmsd empty" {
  let coords : Array[@src.Vector3] = []
  let rmsd = @src.calc_rmsd(coords, coords)
  assert_eq(rmsd, 0.0)
}

test "deform_mode zero amplitude ANM" {
  let coords = [
    @src.Vector3::new(1.0, 2.0, 3.0),
    @src.Vector3::new(4.0, 5.0, 6.0),
  ]
  let evec = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
  let result = @src.deform_mode(coords, evec, 0.0)
  assert_true((result[0].x - 1.0).abs() < 1.0e-10)
  assert_true((result[0].y - 2.0).abs() < 1.0e-10)
  assert_true((result[1].x - 4.0).abs() < 1.0e-10)
}

test "deform_mode nonzero amplitude ANM" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 0.0),
  ]
  let evec = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
  let result = @src.deform_mode(coords, evec, 2.0)
  assert_true((result[0].x - 2.0).abs() < 1.0e-10)
  assert_true((result[0].y - 0.0).abs() < 1.0e-10)
  assert_true((result[1].x - 0.0).abs() < 1.0e-10)
  assert_true((result[1].y - 2.0).abs() < 1.0e-10)
}

test "deform_mode GNM isotropic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 1.0, 1.0),
  ]
  let evec = [1.0, -1.0]
  let result = @src.deform_mode(coords, evec, 3.0)
  assert_true((result[0].x - 3.0).abs() < 1.0e-10)
  assert_true((result[0].y - 3.0).abs() < 1.0e-10)
  assert_true((result[0].z - 3.0).abs() < 1.0e-10)
  assert_true((result[1].x - (-2.0)).abs() < 1.0e-10)
  assert_true((result[1].z - (-2.0)).abs() < 1.0e-10)
}

test "deform_mode mismatched evec length" {
  let coords = [
    @src.Vector3::new(1.0, 2.0, 3.0),
  ]
  let evec = [1.0, 2.0]
  let result = @src.deform_mode(coords, evec, 5.0)
  assert_true((result[0].x - 1.0).abs() < 1.0e-10)
  assert_true((result[0].y - 2.0).abs() < 1.0e-10)
}

test "traverse_mode basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
  ]
  let evec = [1.0, 0.0, 0.0]
  let traj = @src.traverse_mode(coords, evec, 5, 1.0)
  assert_eq(traj.length(), 5)
  assert_true((traj[0][0].x - (-1.0)).abs() < 1.0e-10)
  assert_true((traj[2][0].x - 0.0).abs() < 1.0e-10)
  assert_true((traj[4][0].x - 1.0).abs() < 1.0e-10)
}

test "traverse_mode single step" {
  let coords = [
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evec = [1.0, 0.0, 0.0]
  let traj = @src.traverse_mode(coords, evec, 1, 5.0)
  assert_eq(traj.length(), 1)
  assert_true((traj[0][0].x - 1.0).abs() < 1.0e-10)
}

test "traverse_mode zero steps" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
  ]
  let evec = [1.0, 0.0, 0.0]
  let traj = @src.traverse_mode(coords, evec, 0, 1.0)
  assert_eq(traj.length(), 0)
}

test "traverse_mode symmetric" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
  ]
  let evec = [0.0, 1.0, 0.0]
  let traj = @src.traverse_mode(coords, evec, 3, 2.0)
  assert_true((traj[0][0].y - (-2.0)).abs() < 1.0e-10)
  assert_true((traj[1][0].y - 0.0).abs() < 1.0e-10)
  assert_true((traj[2][0].y - 2.0).abs() < 1.0e-10)
}

test "traverse_mode_amplitudes explicit" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
  ]
  let evec = [1.0, 0.0, 0.0]
  let amps = [0.5, 1.5, -0.3]
  let traj = @src.traverse_mode_amplitudes(coords, evec, amps)
  assert_eq(traj.length(), 3)
  assert_true((traj[0][0].x - 0.5).abs() < 1.0e-10)
  assert_true((traj[1][0].x - 1.5).abs() < 1.0e-10)
  assert_true((traj[2][0].x - (-0.3)).abs() < 1.0e-10)
}

test "sample_anm_conformations basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(3.0, 0.0, 0.0),
  ]
  let eigenvalues = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]
  let eigenvectors : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  eigenvectors[0][6] = 0.7071
  eigenvectors[3][6] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, eigenvalues, eigenvectors, 5, cfg)
  assert_eq(confs.length(), 5)
  for c in confs {
    assert_eq(c.get_coords().length(), 2)
    assert_true(c.get_rmsd() >= 0.0)
  }
}

test "sample_anm reproducibility" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(8, 0.0) })
  evecs[0][6] = 0.5
  evecs[3][6] = -0.5
  evecs[1][7] = 0.5
  evecs[4][7] = -0.5
  let cfg = @src.SamplingConfig::default()
  let c1 = @src.sample_anm_conformations(coords, evals, evecs, 3, cfg)
  let c2 = @src.sample_anm_conformations(coords, evals, evecs, 3, cfg)
  for i in 0..<3 {
    assert_true((c1[i].get_rmsd() - c2[i].get_rmsd()).abs() < 1.0e-10)
  }
}

test "sample_anm different seeds" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 0.7071
  evecs[3][6] = -0.7071
  let cfg1 = @src.SamplingConfig::new(seed=1)
  let cfg2 = @src.SamplingConfig::new(seed=2)
  let c1 = @src.sample_anm_conformations(coords, evals, evecs, 1, cfg1)
  let c2 = @src.sample_anm_conformations(coords, evals, evecs, 1, cfg2)
  assert_true((c1[0].get_rmsd() - c2[0].get_rmsd()).abs() > 1.0e-6)
}

test "sample_anm zero scale" {
  let coords = [
    @src.Vector3::new(1.0, 2.0, 3.0),
    @src.Vector3::new(4.0, 5.0, 6.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 1.0
  let cfg = @src.SamplingConfig::new(scale=0.0)
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 3, cfg)
  for c in confs {
    assert_true(c.get_rmsd().abs() < 1.0e-10)
  }
}

test "sample_anm num_modes limit" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(9, 0.0) })
  evecs[0][6] = 0.5
  evecs[1][7] = 0.5
  evecs[2][8] = 0.5
  let cfg = @src.SamplingConfig::new(num_modes=1)
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 2, cfg)
  for c in confs {
    assert_eq(c.get_coefficients().length(), 1)
  }
}

test "sample_gnm basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 1.0, 1.0),
    @src.Vector3::new(2.0, 2.0, 2.0),
  ]
  let evals = [0.0, 1.0, 2.0]
  let evecs : Array[Array[Double]] = Array::makei(3, { _ => Array::make(3, 0.0) })
  evecs[0][1] = 0.5774
  evecs[1][1] = -0.5774
  evecs[2][1] = 0.5774
  evecs[0][2] = 0.8165
  evecs[2][2] = -0.8165
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_gnm_conformations(coords, evals, evecs, 5, cfg)
  assert_eq(confs.length(), 5)
  for c in confs {
    assert_eq(c.get_coords().length(), 3)
    assert_true(c.get_rmsd() >= 0.0)
  }
}

test "sample_gnm reproducibility" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(2, { _ => Array::make(2, 0.0) })
  evecs[0][1] = 0.7071
  evecs[1][1] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let c1 = @src.sample_gnm_conformations(coords, evals, evecs, 3, cfg)
  let c2 = @src.sample_gnm_conformations(coords, evals, evecs, 3, cfg)
  for i in 0..<3 {
    assert_true((c1[i].get_rmsd() - c2[i].get_rmsd()).abs() < 1.0e-10)
  }
}

test "sample_gnm zero scale" {
  let coords = [
    @src.Vector3::new(1.0, 2.0, 3.0),
    @src.Vector3::new(4.0, 5.0, 6.0),
  ]
  let evals = [0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(2, { _ => Array::make(2, 0.0) })
  evecs[0][1] = 0.7071
  evecs[1][1] = -0.7071
  let cfg = @src.SamplingConfig::new(scale=0.0)
  let confs = @src.sample_gnm_conformations(coords, evals, evecs, 3, cfg)
  for c in confs {
    assert_true(c.get_rmsd().abs() < 1.0e-10)
  }
}

test "average_rmsd basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 0.7071
  evecs[3][6] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 10, cfg)
  let avg = @src.average_rmsd(confs)
  assert_true(avg >= 0.0)
}

test "average_rmsd empty" {
  let confs : Array[@src.Conformation] = []
  assert_eq(@src.average_rmsd(confs), 0.0)
}

test "rmsd_stddev basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 0.7071
  evecs[3][6] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 10, cfg)
  let sd = @src.rmsd_stddev(confs)
  assert_true(sd >= 0.0)
}

test "max_rmsd_conformation" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 0.7071
  evecs[3][6] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 10, cfg)
  let max_c = match @src.max_rmsd_conformation(confs) {
    Some(c) => c
    None => @src.Conformation::new([], [], 0.0)
  }
  for c in confs {
    assert_true(c.get_rmsd() <= max_c.get_rmsd() + 1.0e-10)
  }
}

test "min_rmsd_conformation" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 0.7071
  evecs[3][6] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 10, cfg)
  let min_c = match @src.min_rmsd_conformation(confs) {
    Some(c) => c
    None => @src.Conformation::new([], [], 1.0e30)
  }
  for c in confs {
    assert_true(c.get_rmsd() >= min_c.get_rmsd() - 1.0e-10)
  }
}

test "mean_conformation basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 0.7071
  evecs[3][6] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 20, cfg)
  let mean = match @src.mean_conformation(confs) {
    Some(m) => m
    None => []
  }
  assert_eq(mean.length(), 2)
}

test "ensemble_rmsf basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
    @src.Vector3::new(2.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0]
  let evecs : Array[Array[Double]] = Array::makei(9, { _ => Array::make(9, 0.0) })
  evecs[0][6] = 0.5
  evecs[3][6] = -0.5
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 20, cfg)
  let rmsf = @src.ensemble_rmsf(confs)
  assert_eq(rmsf.length(), 3)
  for v in rmsf {
    assert_true(v >= 0.0)
  }
}

test "ensemble_rmsd_matrix basic" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let evals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let evecs : Array[Array[Double]] = Array::makei(6, { _ => Array::make(7, 0.0) })
  evecs[0][6] = 0.7071
  evecs[3][6] = -0.7071
  let cfg = @src.SamplingConfig::default()
  let confs = @src.sample_anm_conformations(coords, evals, evecs, 5, cfg)
  let mat = @src.ensemble_rmsd_matrix(confs)
  assert_eq(mat.length(), 5)
  // Diagonal should be zero
  for i in 0..<5 {
    assert_true(mat[i][i].abs() < 1.0e-10)
  }
  // Symmetric
  for i in 0..<5 {
    for j in 0..<5 {
      assert_true((mat[i][j] - mat[j][i]).abs() < 1.0e-10)
    }
  }
}
