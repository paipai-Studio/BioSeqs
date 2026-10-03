///|
/// Tests for mANM — Membrane ANM.

///|
/// Build a transmembrane-like chain: helix crossing the membrane
/// (z from -12 to +12), with residues 2-7 in the membrane core.
fn make_tm_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let coords = [
    (0.0, 0.0, -15.0),   // 0: outside (lower)
    (0.0, 0.0, -12.0),   // 1: membrane edge
    (0.0, 0.0, -8.0),    // 2: membrane core
    (0.0, 0.0, -4.0),    // 3: membrane core
    (0.0, 0.0, 0.0),     // 4: membrane core
    (0.0, 0.0, 4.0),     // 5: membrane core
    (0.0, 0.0, 8.0),     // 6: membrane core
    (0.0, 0.0, 12.0),    // 7: membrane edge
    (0.0, 0.0, 15.0),    // 8: outside (upper)
    (0.0, 0.0, 18.0),    // 9: outside (upper)
  ]
  for i = 0; i < 10; i = i + 1 {
    let (x, y, z) = coords[i]
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, y, z), resname="ALA", chainid='A', resseq=i + 1, element="C")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=i + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
/// Build a chain with 2D spread (for contact testing).
fn make_spread_chain() -> @src.Chain {
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

// ===== mANM tests =====

test "manm_construction" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]  // residues 2-5 in membrane
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      assert_eq(@src.MANM::num_nodes(anm), 10)
      assert_eq(@src.MANM::num_membrane(anm), 4)
      assert_true(@src.MANM::membrane_restraint(anm) - 1.0 < 0.001)
    }
    None => assert_true(false)
  }
}

test "manm_custom_restraint" {
  let chain = make_spread_chain()
  let mem_res = [3, 4, 5]
  match @src.MANM::new(chain, mem_res, k_membrane=5.0, k_lateral=0.5) {
    Some(anm) => {
      assert_true(@src.MANM::membrane_restraint(anm) - 5.0 < 0.001)
      assert_true(@src.MANM::lateral_restraint(anm) - 0.5 < 0.001)
      assert_eq(@src.MANM::num_membrane(anm), 3)
    }
    None => assert_true(false)
  }
}

test "manm_eigenvalues" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      let evals = @src.MANM::get_eigenvalues(anm)
      assert_eq(evals.length(), 30)
      // First 6 near-zero (3 trans + 3 rot)
      for i = 0; i < 6; i = i + 1 {
        assert_true(evals[i].abs() < 0.5)
      }
      // Remaining non-negative
      for i = 6; i < 30; i = i + 1 {
        assert_true(evals[i] > -0.5)
      }
    }
    None => assert_true(false)
  }
}

test "manm_bfactors" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      let bf = @src.MANM::predict_bfactors(anm)
      assert_eq(bf.length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(bf[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "manm_fluctuations" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      let fluc = @src.MANM::fluctuations(anm)
      let bf = @src.MANM::predict_bfactors(anm)
      for i = 0; i < 10; i = i + 1 {
        assert_true((fluc[i] - bf[i]).abs() < 1.0e-10)
      }
    }
    None => assert_true(false)
  }
}

test "manm_slow_modes" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      let modes = @src.MANM::slow_modes(anm, 3)
      assert_eq(modes.length(), 3)
      for i = 0; i < 3; i = i + 1 {
        let (_, eval, evec) = modes[i]
        assert_true(eval > -0.5)
        assert_eq(evec.length(), 30)
      }
    }
    None => assert_true(false)
  }
}

test "manm_cross_correlation" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      let corr = @src.MANM::cross_correlation(anm)
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

test "manm_hessian_symmetric" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      let h = @src.MANM::hessian_matrix(anm)
      let dim = 30
      for i = 0; i < dim; i = i + 1 {
        for j = i + 1; j < dim; j = j + 1 {
          assert_true((h[i][j] - h[j][i]).abs() < 1.0e-10)
        }
      }
    }
    None => assert_true(false)
  }
}

test "manm_no_membrane_equals_anm" {
  // With no membrane residues, mANM should match standard ANM
  let chain = make_spread_chain()
  let empty_mem : Array[Int] = Array(capacity=0)
  match @src.MANM::new(chain, empty_mem) {
    Some(m_anm) => {
      match @src.ANM::new(chain) {
        Some(std_anm) => {
          let m_bf = @src.MANM::predict_bfactors(m_anm)
          let std_bf = @src.ANM::predict_bfactors(std_anm)
          for i = 0; i < 10; i = i + 1 {
            assert_true((m_bf[i] - std_bf[i]).abs() < 0.001)
          }
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "manm_with_membrane_differs_from_anm" {
  // With membrane restraints, B-factors should differ from standard ANM
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res, k_membrane=2.0) {
    Some(m_anm) => {
      match @src.ANM::new(chain) {
        Some(std_anm) => {
          let m_bf = @src.MANM::predict_bfactors(m_anm)
          let std_bf = @src.ANM::predict_bfactors(std_anm)
          let mut diff = 0.0
          for i = 0; i < 10; i = i + 1 {
            diff = diff + (m_bf[i] - std_bf[i]).abs()
          }
          assert_true(diff > 0.001)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "manm_stronger_restraint_reduces_fluctuation" {
  // Stronger membrane restraint should reduce z-fluctuations of
  // membrane residues relative to non-membrane residues
  let chain = make_spread_chain()
  let mem_res = [4]
  match @src.MANM::new(chain, mem_res, k_membrane=0.1) {
    Some(weak) => {
      match @src.MANM::new(chain, mem_res, k_membrane=10.0) {
        Some(strong) => {
          let bf_weak = @src.MANM::predict_bfactors(weak)
          let bf_strong = @src.MANM::predict_bfactors(strong)
          // Total fluctuation should decrease with stronger restraint
          let mut total_weak = 0.0
          let mut total_strong = 0.0
          for i = 0; i < 10; i = i + 1 {
            total_weak = total_weak + bf_weak[i]
            total_strong = total_strong + bf_strong[i]
          }
          assert_true(total_strong < total_weak)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "manm_lateral_restraint" {
  let chain = make_spread_chain()
  let mem_res = [2, 3, 4, 5]
  match @src.MANM::new(chain, mem_res, k_membrane=1.0, k_lateral=0.5) {
    Some(anm) => {
      let bf = @src.MANM::predict_bfactors(anm)
      assert_eq(bf.length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(bf[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "manm_too_few_residues" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA", chainid='A', resseq=1, element="C")
  let res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[ca])
  residues.push(res)
  let chain = @src.Chain::new(id='A', residues~)
  match @src.MANM::new(chain, [0]) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "manm_residues_by_z" {
  let chain = make_tm_chain()
  let mem_res = @src.membrane_residues_by_z(chain, -12.0, 12.0)
  // Residues 1-7 have z in [-12, 12]
  assert_eq(mem_res.length(), 7)
  // Residue 0 (z=-15) and 8 (z=15), 9 (z=18) should be excluded
  let mem_res2 = @src.membrane_residues_by_z(chain, -10.0, 10.0)
  // Residues 2-6 have z in [-10, 10]
  assert_eq(mem_res2.length(), 5)
}

test "manm_invalid_membrane_indices_ignored" {
  // Out-of-range indices should be silently ignored
  let chain = make_spread_chain()
  let mem_res = [2, 3, 100, -1]  // 100 and -1 are invalid
  match @src.MANM::new(chain, mem_res) {
    Some(anm) => {
      assert_eq(@src.MANM::num_nodes(anm), 10)
      let bf = @src.MANM::predict_bfactors(anm)
      assert_eq(bf.length(), 10)
    }
    None => assert_true(false)
  }
}
