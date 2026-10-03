///|
/// Tests for betaANM and pfANM — distance-dependent ANM variants.

///|
fn make_test_chain_beta() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let coords = [
    (1.0, 0.0, 0.0),
    (0.31, 0.95, 1.5),
    (-0.81, 0.59, 3.0),
    (-0.81, -0.59, 4.5),
    (0.31, -0.95, 6.0),
  ]
  for i = 0; i < 5; i = i + 1 {
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

///|
fn make_larger_chain_beta() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let coords = [
    (0.0, 0.0, 0.0),
    (3.8, 0.0, 0.0),
    (7.6, 0.0, 0.0),
    (11.4, 0.0, 0.0),
    (15.2, 0.0, 0.0),
    (15.2, 3.8, 0.0),
    (15.2, 7.6, 0.0),
    (11.4, 7.6, 0.0),
    (7.6, 7.6, 0.0),
    (3.8, 7.6, 0.0),
  ]
  for i = 0; i < 10; i = i + 1 {
    let (x, y, z) = coords[i]
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, y, z), resname="GLY", chainid='A', resseq=i + 1, element="C")
    let res = @src.Residue::new(resname="GLY", chainid='A', resseq=i + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

// ===== betaANM tests =====

test "beta_anm_construction" {
  let chain = make_test_chain_beta()
  match @src.BetaANM::new(chain) {
    Some(anm) => {
      assert_eq(@src.BetaANM::num_nodes(anm), 5)
      assert_true(@src.BetaANM::beta_exponent(anm) - 2.0 < 0.001)
    }
    None => assert_true(false)
  }
}

test "beta_anm_custom_beta" {
  let chain = make_test_chain_beta()
  match @src.BetaANM::new(chain, beta=3.0) {
    Some(anm) => {
      assert_true(@src.BetaANM::beta_exponent(anm) - 3.0 < 0.001)
    }
    None => assert_true(false)
  }
}

test "beta_anm_eigenvalues" {
  let chain = make_test_chain_beta()
  match @src.BetaANM::new(chain) {
    Some(anm) => {
      let evals = @src.BetaANM::get_eigenvalues(anm)
      // 3N = 15 eigenvalues
      assert_eq(evals.length(), 15)
      // First 6 should be near-zero (3 trans + 3 rot)
      for i = 0; i < 6; i = i + 1 {
        assert_true(evals[i].abs() < 0.01)
      }
      // Remaining should be positive
      for i = 6; i < 15; i = i + 1 {
        assert_true(evals[i] > -0.001)
      }
    }
    None => assert_true(false)
  }
}

test "beta_anm_bfactors" {
  let chain = make_larger_chain_beta()
  match @src.BetaANM::new(chain) {
    Some(anm) => {
      let bf = @src.BetaANM::predict_bfactors(anm)
      assert_eq(bf.length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(bf[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "beta_anm_fluctuations" {
  let chain = make_larger_chain_beta()
  match @src.BetaANM::new(chain) {
    Some(anm) => {
      let fluc = @src.BetaANM::fluctuations(anm)
      let bf = @src.BetaANM::predict_bfactors(anm)
      // fluctuations == bfactors
      for i = 0; i < 10; i = i + 1 {
        assert_true((fluc[i] - bf[i]).abs() < 1.0e-10)
      }
    }
    None => assert_true(false)
  }
}

test "beta_anm_slow_modes" {
  let chain = make_larger_chain_beta()
  match @src.BetaANM::new(chain) {
    Some(anm) => {
      let modes = @src.BetaANM::slow_modes(anm, 3)
      assert_eq(modes.length(), 3)
      // Each mode: (index, eigenvalue, eigenvector)
      for i = 0; i < 3; i = i + 1 {
        let (_, eval, evec) = modes[i]
        assert_true(eval > -0.001)
        assert_eq(evec.length(), 30)
      }
    }
    None => assert_true(false)
  }
}

test "beta_anm_cross_correlation" {
  let chain = make_larger_chain_beta()
  match @src.BetaANM::new(chain) {
    Some(anm) => {
      let corr = @src.BetaANM::cross_correlation(anm)
      assert_eq(corr.length(), 10)
      assert_eq(corr[0].length(), 10)
      // Diagonal should be positive
      for i = 0; i < 10; i = i + 1 {
        assert_true(corr[i][i] > 0.0)
      }
      // Symmetry: corr[i][j] ≈ corr[j][i]
      for i = 0; i < 10; i = i + 1 {
        for j = i + 1; j < 10; j = j + 1 {
          assert_true((corr[i][j] - corr[j][i]).abs() < 0.01)
        }
      }
    }
    None => assert_true(false)
  }
}

test "beta_anm_hessian_symmetric" {
  let chain = make_test_chain_beta()
  match @src.BetaANM::new(chain) {
    Some(anm) => {
      let h = @src.BetaANM::hessian_matrix(anm)
      let dim = 15
      // Hessian should be symmetric
      for i = 0; i < dim; i = i + 1 {
        for j = i + 1; j < dim; j = j + 1 {
          assert_true((h[i][j] - h[j][i]).abs() < 1.0e-10)
        }
      }
    }
    None => assert_true(false)
  }
}

test "beta_anm_higher_beta_more_local" {
  let chain = make_larger_chain_beta()
  // Higher beta → more local interactions → different B-factor pattern
  match @src.BetaANM::new(chain, beta=1.0) {
    Some(anm1) => {
      match @src.BetaANM::new(chain, beta=4.0) {
        Some(anm4) => {
          let bf1 = @src.BetaANM::predict_bfactors(anm1)
          let bf4 = @src.BetaANM::predict_bfactors(anm4)
          // B-factor patterns should differ
          let mut diff = 0.0
          for i = 0; i < 10; i = i + 1 {
            diff = diff + (bf1[i] - bf4[i]).abs()
          }
          assert_true(diff > 0.001)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "beta_anm_too_few_residues" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA", chainid='A', resseq=1, element="C")
  let res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[ca])
  residues.push(res)
  let chain = @src.Chain::new(id='A', residues~)
  match @src.BetaANM::new(chain) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ===== pfANM tests =====

test "pf_anm_construction" {
  let chain = make_test_chain_beta()
  match @src.PfANM::new(chain) {
    Some(anm) => {
      assert_eq(@src.PfANM::num_nodes(anm), 5)
    }
    None => assert_true(false)
  }
}

test "pf_anm_eigenvalues" {
  let chain = make_test_chain_beta()
  match @src.PfANM::new(chain) {
    Some(anm) => {
      let evals = @src.PfANM::get_eigenvalues(anm)
      assert_eq(evals.length(), 15)
      // First 6 near-zero
      for i = 0; i < 6; i = i + 1 {
        assert_true(evals[i].abs() < 0.01)
      }
      // Remaining positive
      for i = 6; i < 15; i = i + 1 {
        assert_true(evals[i] > -0.001)
      }
    }
    None => assert_true(false)
  }
}

test "pf_anm_bfactors" {
  let chain = make_larger_chain_beta()
  match @src.PfANM::new(chain) {
    Some(anm) => {
      let bf = @src.PfANM::predict_bfactors(anm)
      assert_eq(bf.length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(bf[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "pf_anm_fluctuations" {
  let chain = make_larger_chain_beta()
  match @src.PfANM::new(chain) {
    Some(anm) => {
      let fluc = @src.PfANM::fluctuations(anm)
      let bf = @src.PfANM::predict_bfactors(anm)
      for i = 0; i < 10; i = i + 1 {
        assert_true((fluc[i] - bf[i]).abs() < 1.0e-10)
      }
    }
    None => assert_true(false)
  }
}

test "pf_anm_slow_modes" {
  let chain = make_larger_chain_beta()
  match @src.PfANM::new(chain) {
    Some(anm) => {
      let modes = @src.PfANM::slow_modes(anm, 2)
      assert_eq(modes.length(), 2)
      for i = 0; i < 2; i = i + 1 {
        let (_, eval, evec) = modes[i]
        assert_true(eval > -0.001)
        assert_eq(evec.length(), 30)
      }
    }
    None => assert_true(false)
  }
}

test "pf_anm_cross_correlation" {
  let chain = make_larger_chain_beta()
  match @src.PfANM::new(chain) {
    Some(anm) => {
      let corr = @src.PfANM::cross_correlation(anm)
      assert_eq(corr.length(), 10)
      assert_eq(corr[0].length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(corr[i][i] > 0.0)
      }
      for i = 0; i < 10; i = i + 1 {
        for j = i + 1; j < 10; j = j + 1 {
          assert_true((corr[i][j] - corr[j][i]).abs() < 0.01)
        }
      }
    }
    None => assert_true(false)
  }
}

test "pf_anm_hessian_symmetric" {
  let chain = make_test_chain_beta()
  match @src.PfANM::new(chain) {
    Some(anm) => {
      let h = @src.PfANM::hessian_matrix(anm)
      let dim = 15
      for i = 0; i < dim; i = i + 1 {
        for j = i + 1; j < dim; j = j + 1 {
          assert_true((h[i][j] - h[j][i]).abs() < 1.0e-10)
        }
      }
    }
    None => assert_true(false)
  }
}

test "pf_anm_no_cutoff_needed" {
  // pfANM with very large cutoff should still work
  let chain = make_larger_chain_beta()
  match @src.PfANM::new(chain, cutoff=100.0) {
    Some(anm) => {
      let bf = @src.PfANM::predict_bfactors(anm)
      assert_eq(bf.length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(bf[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "pf_anm_too_few_residues" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA", chainid='A', resseq=1, element="C")
  let res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[ca])
  residues.push(res)
  let chain = @src.Chain::new(id='A', residues~)
  match @src.PfANM::new(chain) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ===== Comparison: betaANM vs ANM vs pfANM =====

test "beta_vs_standard_anm_differs" {
  // With beta=0, betaANM should reduce to standard ANM
  // With beta>0, results should differ
  let chain = make_larger_chain_beta()
  match @src.ANM::new(chain) {
    Some(std_anm) => {
      match @src.BetaANM::new(chain, beta=2.0) {
        Some(beta_anm) => {
          let std_bf = @src.ANM::predict_bfactors(std_anm)
          let beta_bf = @src.BetaANM::predict_bfactors(beta_anm)
          // B-factor patterns should differ
          let mut diff = 0.0
          for i = 0; i < 10; i = i + 1 {
            diff = diff + (std_bf[i] - beta_bf[i]).abs()
          }
          assert_true(diff > 0.001)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "pf_vs_beta_differs" {
  let chain = make_larger_chain_beta()
  match @src.PfANM::new(chain) {
    Some(pf_anm) => {
      match @src.BetaANM::new(chain, beta=2.0) {
        Some(beta_anm) => {
          let pf_bf = @src.PfANM::predict_bfactors(pf_anm)
          let beta_bf = @src.BetaANM::predict_bfactors(beta_anm)
          let mut diff = 0.0
          for i = 0; i < 10; i = i + 1 {
            diff = diff + (pf_bf[i] - beta_bf[i]).abs()
          }
          assert_true(diff > 0.001)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}
