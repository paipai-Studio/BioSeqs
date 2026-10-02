///|
/// Tier 3 ProDy tests: mode comparison, RTB-ANM, essential dynamics,
/// comparative dynamics, signature dynamics.

///|
fn make_t3_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let coords = [
    (1.0, 0.0, 0.0),
    (0.31, 0.95, 1.5),
    (-0.81, 0.59, 3.0),
    (-0.81, -0.59, 4.5),
    (0.31, -0.95, 6.0),
    (1.0, 0.0, 7.5),
    (0.31, 0.95, 9.0),
    (-0.81, 0.59, 10.5),
    (-0.81, -0.59, 12.0),
    (0.31, -0.95, 13.5),
    (1.0, 0.0, 15.0),
    (0.31, 0.95, 16.5),
  ]
  for i = 0; i < 12; i = i + 1 {
    let (x, y, z) = coords[i]
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, y, z), resname="ALA", chainid='A', resseq=i + 1, element="C")
    let n_atom = @src.Atom::new(name="N", coord=@src.Vector3::new(x - 0.5, y, z - 0.3), resname="ALA", chainid='A', resseq=i + 1, element="N")
    let c_atom = @src.Atom::new(name="C", coord=@src.Vector3::new(x + 0.5, y, z + 0.3), resname="ALA", chainid='A', resseq=i + 1, element="C")
    let o_atom = @src.Atom::new(name="O", coord=@src.Vector3::new(x + 0.7, y + 0.2, z + 0.5), resname="ALA", chainid='A', resseq=i + 1, element="O")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=i + 1, atoms=[n_atom, ca, c_atom, o_atom])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

// === Mode Comparison ===

test "mode_overlap_identical" {
  let v = [1.0, 0.0, 0.0, 0.0]
  let ov = @src.mode_overlap(v, v)
  assert_true((ov - 1.0).abs() < 0.001)
}

test "mode_overlap_orthogonal" {
  let v1 = [1.0, 0.0, 0.0, 0.0]
  let v2 = [0.0, 1.0, 0.0, 0.0]
  let ov = @src.mode_overlap(v1, v2)
  assert_true(ov < 0.001)
}

test "mode_overlap_general" {
  let v1 = [0.6, 0.8, 0.0, 0.0]
  let v2 = [0.8, 0.6, 0.0, 0.0]
  let ov = @src.mode_overlap(v1, v2)
  // dot = 0.48 + 0.48 = 0.96
  assert_true((ov - 0.96).abs() < 0.001)
}

test "cumulative_overlap" {
  let mode = [1.0, 0.0, 0.0]
  let set = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let co = @src.cumulative_overlap(mode, set, 2)
  assert_true((co - 1.0).abs() < 0.001)
}

test "scaling_factor" {
  let v1 = [2.0, 4.0]
  let v2 = [1.0, 2.0]
  let s = @src.scaling_factor(v1, v2)
  assert_true((s - 2.0).abs() < 0.001)
}

test "scaling_factor_orthogonal" {
  let v1 = [1.0, 0.0]
  let v2 = [0.0, 1.0]
  let s = @src.scaling_factor(v1, v2)
  assert_true(s.abs() < 0.001)
}

test "overlap_matrix_shape" {
  let m1 = [[1.0, 0.0], [0.0, 1.0]]
  let m2 = [[1.0, 0.0], [0.0, 1.0]]
  let om = @src.overlap_matrix(m1, m2)
  assert_eq(om.length(), 2)
  assert_eq(om[0].length(), 2)
}

test "overlap_matrix_diagonal" {
  let m1 = [[1.0, 0.0], [0.0, 1.0]]
  let m2 = [[1.0, 0.0], [0.0, 1.0]]
  let om = @src.overlap_matrix(m1, m2)
  assert_true((om[0][0] - 1.0).abs() < 0.001)
  assert_true(om[0][1] < 0.001)
}

test "effective_overlap" {
  let mode = [1.0, 0.0, 0.0]
  let set = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let eo = @src.effective_overlap(mode, set, 1)
  assert_true((eo - 1.0).abs() < 0.001)
}

test "subspace_angle_identical" {
  let m1 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let m2 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let angle = @src.subspace_angle(m1, m2)
  assert_true(angle < 1.0)
}

test "subspace_angle_orthogonal" {
  let m1 = [[1.0, 0.0, 0.0]]
  let m2 = [[0.0, 1.0, 0.0]]
  let angle = @src.subspace_angle(m1, m2)
  assert_true((angle - 90.0).abs() < 1.0)
}

// === RTB-ANM ===

test "rtb_anm_construction" {
  let chain = make_t3_chain()
  let residues = chain.get_residues()
  let coords : Array[@src.Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => coords.push(ca.coord)
      None => ()
    }
  }
  match @src.RTBANM::new(coords, block_size=4) {
    Some(rtb) => {
      assert_eq(rtb.num_residues(), 12)
      assert_true(rtb.num_blocks() >= 3)
    }
    None => assert_true(false)
  }
}

test "rtb_anm_blocks" {
  let chain = make_t3_chain()
  let residues = chain.get_residues()
  let coords : Array[@src.Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => coords.push(ca.coord)
      None => ()
    }
  }
  match @src.RTBANM::new(coords, block_size=3) {
    Some(rtb) => {
      assert_eq(rtb.num_blocks(), 4)
      let sizes = rtb.get_block_sizes()
      let total = {
        let mut s = 0
        for sz in sizes {
          s = s + sz
        }
        s
      }
      assert_eq(total, 12)
    }
    None => assert_true(false)
  }
}

test "rtb_anm_eigenvalues" {
  let chain = make_t3_chain()
  let residues = chain.get_residues()
  let coords : Array[@src.Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => coords.push(ca.coord)
      None => ()
    }
  }
  match @src.RTBANM::new(coords, block_size=4) {
    Some(rtb) => {
      let evals = rtb.get_eigenvalues()
      assert_true(evals.length() > 0)
      assert_true(evals[0] >= -0.001)

    }
    None => assert_true(false)
  }
}

test "rtb_anm_predict_bfactors" {
  let chain = make_t3_chain()
  let residues = chain.get_residues()
  let coords : Array[@src.Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => coords.push(ca.coord)
      None => ()
    }
  }
  match @src.RTBANM::new(coords, block_size=4) {
    Some(rtb) => {
      let bf = rtb.predict_bfactors()
      assert_eq(bf.length(), 12)
      let mut all_nonneg = true
      for b in bf {
        if b < -0.001 {
          all_nonneg = false
        }
      }
      assert_true(all_nonneg)
    }
    None => assert_true(false)
  }
}

// === Essential Dynamics ===

test "essential_dynamics_basic" {
  let confs = [
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 2.0, 0.0, 0.0],
    [2.0, 0.0, 0.0, 3.0, 0.0, 0.0],
    [3.0, 0.0, 0.0, 4.0, 0.0, 0.0],
  ]
  match @src.EssentialDynamics::new(confs) {
    Some(eda) => {
      assert_eq(eda.n_conformations, 4)
      assert_eq(eda.n_coords, 6)
    }
    None => assert_true(false)
  }
}

test "essential_dynamics_variance_explained" {
  let confs = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [-2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
  ]
  match @src.EssentialDynamics::new(confs) {
    Some(eda) => {
      let ve = eda.variance_explained(1)
      assert_true(ve > 0.0)
      let ve_all = eda.variance_explained(6)
      assert_true(ve_all > 0.99)
    }
    None => assert_true(false)
  }
}

test "essential_dynamics_project" {
  let confs = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 2.0, 0.0, 0.0],
    [2.0, 0.0, 0.0, 4.0, 0.0, 0.0],
  ]
  match @src.EssentialDynamics::new(confs) {
    Some(eda) => {
      let proj = eda.project(confs[1], 2)
      assert_eq(proj.length(), 2)
    }
    None => assert_true(false)
  }
}

test "essential_dynamics_reconstruct" {
  let confs = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 2.0, 0.0, 0.0],
    [2.0, 0.0, 0.0, 4.0, 0.0, 0.0],
  ]
  match @src.EssentialDynamics::new(confs) {
    Some(eda) => {
      let proj = eda.project(confs[1], 6)
      let recon = eda.reconstruct(proj)
      assert_eq(recon.length(), 6)
    }
    None => assert_true(false)
  }
}

test "essential_dynamics_rmsf" {
  let confs = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 2.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0, -2.0, 0.0, 0.0],
  ]
  match @src.EssentialDynamics::new(confs) {
    Some(eda) => {
      let rmsf = eda.rmsf()
      assert_eq(rmsf.length(), 2)
      assert_true(rmsf[0] >= 0.0)
    }
    None => assert_true(false)
  }
}

test "essential_dynamics_num_essential" {
  let confs = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 2.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, -2.0, 0.0, 0.0],
  ]
  match @src.EssentialDynamics::new(confs) {
    Some(eda) => {
      let n = eda.num_essential_modes(threshold=0.95)
      assert_true(n >= 1)
    }
    None => assert_true(false)
  }
}

// === Comparative Dynamics ===

test "compare_fluctuations_identical" {
  let f = [1.0, 2.0, 3.0, 4.0]
  let (rmsd, corr) = @src.compare_fluctuations(f, f)
  assert_true(rmsd < 0.001)
  assert_true((corr - 1.0).abs() < 0.001)
}

test "compare_fluctuations_anticorrelated" {
  let f1 = [1.0, 2.0, 3.0, 4.0]
  let f2 = [4.0, 3.0, 2.0, 1.0]
  let (_, corr) = @src.compare_fluctuations(f1, f2)
  assert_true(corr < -0.99)
}

test "compare_fluctuations_orthogonal" {
  let f1 = [1.0, -1.0, 1.0, -1.0]
  let f2 = [1.0, 1.0, -1.0, -1.0]
  let (_, corr) = @src.compare_fluctuations(f1, f2)
  assert_true(corr.abs() < 0.01)
}

test "compare_bfactors" {
  let bf1 = [1.0, 2.0, 3.0, 4.0]
  let bf2 = [1.0, 2.0, 3.0, 4.0]
  let (corr, rmsd) = @src.compare_bfactors(bf1, bf2)
  assert_true((corr - 1.0).abs() < 0.001)
  assert_true(rmsd < 0.001)
}

test "average_collectivity" {
  let modes = [
    [0.5, 0.5, 0.5, 0.5],
    [0.7, 0.7, 0.0, 0.0],
  ]
  let ac = @src.average_collectivity(modes, 2)
  assert_true(ac > 0.0 && ac <= 1.01)
}

test "dynamics_similarity_identical" {
  let modes = [[1.0, 0.0, 0.0, 0.0]]
  let fluc = [1.0, 0.0, 0.0, 0.0]
  let sim = @src.dynamics_similarity(modes, modes, fluc, fluc, 1)
  assert_true(sim > 0.5)
}

// === Signature Dynamics ===

test "dynamics_signature_basic" {
  let evals = [0.0, 1.0, 2.0, 3.0]
  let evecs = [
    [0.5, 0.5, 0.5, 0.5],
    [0.5, 0.5, -0.5, -0.5],
    [0.5, -0.5, 0.5, -0.5],
    [0.5, -0.5, -0.5, 0.5],
  ]
  let sig = @src.dynamics_signature(evals, evecs, n_modes=3)
  assert_eq(sig.n_residues, 4)
  assert_eq(sig.msf_profile.length(), 4)
}

test "dynamics_signature_stats" {
  let evals = [0.0, 1.0, 2.0]
  let evecs = [
    [0.5, 0.5, 0.5, 0.5],
    [0.5, 0.5, -0.5, -0.5],
    [0.5, -0.5, 0.5, -0.5],
  ]
  let sig = @src.dynamics_signature(evals, evecs, n_modes=2)
  assert_true(sig.mean_fluctuation >= 0.0)
  assert_true(sig.max_fluctuation >= sig.mean_fluctuation)
  assert_true(sig.min_fluctuation <= sig.mean_fluctuation)
}

test "rank_values" {
  let ranks = @src.rank_values([3.0, 1.0, 4.0, 2.0])
  assert_eq(ranks[0], 3.0) // 3.0 is 3rd
  assert_eq(ranks[1], 1.0) // 1.0 is 1st
  assert_eq(ranks[2], 4.0) // 4.0 is 4th
  assert_eq(ranks[3], 2.0) // 2.0 is 2nd
}

test "rank_values_ties" {
  let ranks = @src.rank_values([1.0, 1.0, 3.0, 3.0])
  // Both 1.0s get average rank (1+2)/2 = 1.5
  // Both 3.0s get average rank (3+4)/2 = 3.5
  assert_true((ranks[0] - 1.5).abs() < 0.001)
  assert_true((ranks[1] - 1.5).abs() < 0.001)
  assert_true((ranks[2] - 3.5).abs() < 0.001)
  assert_true((ranks[3] - 3.5).abs() < 0.001)
}

test "rank_correlation_perfect" {
  let v1 = [1.0, 2.0, 3.0, 4.0]
  let v2 = [2.0, 4.0, 6.0, 8.0]
  let rc = @src.rank_correlation(v1, v2)
  assert_true((rc - 1.0).abs() < 0.001)
}

test "rank_correlation_reverse" {
  let v1 = [1.0, 2.0, 3.0, 4.0]
  let v2 = [4.0, 3.0, 2.0, 1.0]
  let rc = @src.rank_correlation(v1, v2)
  assert_true((rc + 1.0).abs() < 0.001)
}

test "signature_correlation_identical" {
  let evals = [0.0, 1.0, 2.0]
  let evecs = [
    [0.5, 0.5, 0.5, 0.5],
    [0.7, 0.5, 0.3, 0.1],
    [0.1, 0.3, 0.5, 0.7],
  ]
  let sig = @src.dynamics_signature(evals, evecs, n_modes=3)
  let corr = @src.signature_correlation(sig, sig)
  assert_true((corr - 1.0).abs() < 0.001)
}

test "signature_distance_zero" {
  let evals = [0.0, 1.0, 2.0]
  let evecs = [
    [0.5, 0.5, 0.5, 0.5],
    [0.5, 0.5, -0.5, -0.5],
    [0.5, -0.5, 0.5, -0.5],
  ]
  let sig = @src.dynamics_signature(evals, evecs, n_modes=2)
  let dist = @src.signature_distance(sig, sig)
  assert_true(dist < 0.001)
}

test "rank_signature_identical" {
  let evals = [0.0, 1.0, 2.0]
  let evecs = [
    [0.5, 0.5, 0.5, 0.5],
    [0.7, 0.5, 0.3, 0.1],
    [0.1, 0.3, 0.5, 0.7],
  ]
  let sig = @src.dynamics_signature(evals, evecs, n_modes=3)
  let rs = @src.rank_signature(sig, sig)
  assert_true((rs - 1.0).abs() < 0.001)
}
