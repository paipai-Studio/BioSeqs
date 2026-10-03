///|
/// Tests for ccANM — Crystal Contact ANM.

///|
fn make_test_chain_cc() -> @src.Chain {
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

///|
fn make_small_chain_cc() -> @src.Chain {
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
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=i + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
fn make_crystal_p21() -> @src.CrystalInfo {
  @src.CrystalInfo::{
    a: 30.0,
    b: 40.0,
    c: 50.0,
    alpha: 90.0,
    beta: 90.0,
    gamma: 90.0,
    sym_ops: @src.spacegroup_p21(),
  }
}

///|
fn make_crystal_p212121() -> @src.CrystalInfo {
  @src.CrystalInfo::{
    a: 30.0,
    b: 40.0,
    c: 50.0,
    alpha: 90.0,
    beta: 90.0,
    gamma: 90.0,
    sym_ops: @src.spacegroup_p212121(),
  }
}

///|
fn make_crystal_p1_large() -> @src.CrystalInfo {
  @src.CrystalInfo::{
    a: 100.0,
    b: 100.0,
    c: 100.0,
    alpha: 90.0,
    beta: 90.0,
    gamma: 90.0,
    sym_ops: @src.spacegroup_p1(),
  }
}

// ===== ccANM tests =====

test "ccanm_construction_p21" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p21()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      assert_eq(@src.CcANM::num_nodes(anm), 10)
    }
    None => assert_true(false)
  }
}

test "ccanm_construction_p212121" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p212121()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      assert_eq(@src.CcANM::num_nodes(anm), 10)
    }
    None => assert_true(false)
  }
}

test "ccanm_eigenvalues" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p21()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      let evals = @src.CcANM::get_eigenvalues(anm)
      assert_eq(evals.length(), 30)
      for i = 0; i < 6; i = i + 1 {
        assert_true(evals[i].abs() < 0.5)
      }
      for i = 6; i < 30; i = i + 1 {
        assert_true(evals[i] > -0.5)
      }
    }
    None => assert_true(false)
  }
}

test "ccanm_bfactors" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p21()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      let bf = @src.CcANM::predict_bfactors(anm)
      assert_eq(bf.length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(bf[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "ccanm_fluctuations" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p212121()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      let fluc = @src.CcANM::fluctuations(anm)
      let bf = @src.CcANM::predict_bfactors(anm)
      for i = 0; i < 10; i = i + 1 {
        assert_true((fluc[i] - bf[i]).abs() < 1.0e-10)
      }
    }
    None => assert_true(false)
  }
}

test "ccanm_slow_modes" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p21()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      let modes = @src.CcANM::slow_modes(anm, 3)
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

test "ccanm_cross_correlation" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p21()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      let corr = @src.CcANM::cross_correlation(anm)
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

test "ccanm_hessian_symmetric" {
  let chain = make_small_chain_cc()
  let crystal = make_crystal_p21()
  match @src.CcANM::new(chain, crystal) {
    Some(anm) => {
      let h = @src.CcANM::hessian_matrix(anm)
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

test "ccanm_p1_equals_standard_anm" {
  let chain = make_test_chain_cc()
  let crystal = make_crystal_p1_large()
  match @src.CcANM::new(chain, crystal) {
    Some(cc_anm) => {
      match @src.ANM::new(chain) {
        Some(std_anm) => {
          let cc_bf = @src.CcANM::predict_bfactors(cc_anm)
          let std_bf = @src.ANM::predict_bfactors(std_anm)
          for i = 0; i < 10; i = i + 1 {
            assert_true((cc_bf[i] - std_bf[i]).abs() < 0.001)
          }
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "ccanm_p21_differs_from_p1" {
  let chain = make_test_chain_cc()
  let crystal_p21 = make_crystal_p21()
  let crystal_p1 = make_crystal_p1_large()
  match @src.CcANM::new(chain, crystal_p21) {
    Some(anm_p21) => {
      match @src.CcANM::new(chain, crystal_p1) {
        Some(anm_p1) => {
          let bf_p21 = @src.CcANM::predict_bfactors(anm_p21)
          let bf_p1 = @src.CcANM::predict_bfactors(anm_p1)
          let mut diff = 0.0
          for i = 0; i < 10; i = i + 1 {
            diff = diff + (bf_p21[i] - bf_p1[i]).abs()
          }
          assert_true(diff > 0.001)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "ccanm_too_few_residues" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA", chainid='A', resseq=1, element="C")
  let res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[ca])
  residues.push(res)
  let chain = @src.Chain::new(id='A', residues~)
  let crystal = make_crystal_p21()
  match @src.CcANM::new(chain, crystal) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "ccanm_spacegroup_helpers" {
  let p1 = @src.spacegroup_p1()
  assert_eq(p1.length(), 1)
  let p21 = @src.spacegroup_p21()
  assert_eq(p21.length(), 2)
  let p212121 = @src.spacegroup_p212121()
  assert_eq(p212121.length(), 4)
  // Identity should always be first
  let id_op = @src.sym_identity()
  assert_true(id_op.rot[0][0] - 1.0 < 0.001)
  assert_true(id_op.trans[0] < 0.001)
}
