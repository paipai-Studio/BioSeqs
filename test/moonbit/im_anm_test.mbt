///|
/// Tests for imANM (mixed ANM) — per-residue gamma values.

///|
fn make_test_chain_im() -> @src.Chain {
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
fn make_larger_chain_im() -> @src.Chain {
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

test "gammas_uniform" {
  let g = @src.gammas_uniform(5, 2.0)
  assert_eq(g.length(), 5)
  for i = 0; i < 5; i = i + 1 {
    assert_true((g[i] - 2.0).abs() < 1.0e-10)
  }
}

test "gammas_by_region basic" {
  let regions = [(0, 3, 2.0), (4, 6, 0.5)]
  let g = @src.gammas_by_region(10, regions, default_gamma=1.0)
  assert_eq(g.length(), 10)
  assert_true((g[0] - 2.0).abs() < 1.0e-10)
  assert_true((g[3] - 2.0).abs() < 1.0e-10)
  assert_true((g[4] - 0.5).abs() < 1.0e-10)
  assert_true((g[6] - 0.5).abs() < 1.0e-10)
  assert_true((g[7] - 1.0).abs() < 1.0e-10)
  assert_true((g[9] - 1.0).abs() < 1.0e-10)
}

test "gammas_by_region default only" {
  let regions : Array[(Int, Int, Double)] = []
  let g = @src.gammas_by_region(5, regions, default_gamma=3.0)
  for i = 0; i < 5; i = i + 1 {
    assert_true((g[i] - 3.0).abs() < 1.0e-10)
  }
}

test "gammas_by_region out of bounds" {
  let regions = [(0, 20, 5.0)]
  let g = @src.gammas_by_region(5, regions, default_gamma=1.0)
  for i = 0; i < 5; i = i + 1 {
    assert_true((g[i] - 5.0).abs() < 1.0e-10)
  }
}

test "gammas_from_bfactors basic" {
  let bf = [10.0, 20.0, 5.0, 10.0]
  let g = @src.gammas_from_bfactors(bf, scale=1.0)
  assert_eq(g.length(), 4)
  assert_true((g[0] - 0.1).abs() < 1.0e-10)
  assert_true((g[1] - 0.05).abs() < 1.0e-10)
  assert_true((g[2] - 0.2).abs() < 1.0e-10)
  assert_true((g[3] - 0.1).abs() < 1.0e-10)
}

test "gammas_from_bfactors with scale" {
  let bf = [10.0, 20.0]
  let g = @src.gammas_from_bfactors(bf, scale=2.0)
  assert_true((g[0] - 0.2).abs() < 1.0e-10)
  assert_true((g[1] - 0.1).abs() < 1.0e-10)
}

test "gammas_from_bfactors zero fallback" {
  let bf = [0.0, 10.0, 0.0]
  let g = @src.gammas_from_bfactors(bf, scale=1.0)
  assert_true((g[1] - 0.1).abs() < 1.0e-10)
  assert_true((g[0] - 0.1).abs() < 1.0e-10)
  assert_true((g[2] - 0.1).abs() < 1.0e-10)
}

test "gammas_by_ss basic" {
  let ss = ["H", "H", "E", "E", "C", "C"]
  let g = @src.gammas_by_secondary_structure(6, ss, helix_gamma=1.0, sheet_gamma=1.0, coil_gamma=0.3)
  assert_eq(g.length(), 6)
  assert_true((g[0] - 1.0).abs() < 1.0e-10)
  assert_true((g[4] - 0.3).abs() < 1.0e-10)
}

test "gammas_by_ss 310 helix" {
  let ss = ["G", "I", "B"]
  let g = @src.gammas_by_secondary_structure(3, ss, helix_gamma=2.0, sheet_gamma=1.5, coil_gamma=0.5)
  assert_true((g[0] - 2.0).abs() < 1.0e-10)
  assert_true((g[1] - 2.0).abs() < 1.0e-10)
  assert_true((g[2] - 1.5).abs() < 1.0e-10)
}

test "gammas_gaussian basic" {
  let g = @src.gammas_gaussian(11, 5, 2.0, base=1.0, amplitude=2.0)
  assert_eq(g.length(), 11)
  assert_true((g[5] - 3.0).abs() < 1.0e-10)
  assert_true(g[0] > 1.0)
  assert_true(g[0] < 2.0)
  assert_true((g[4] - g[6]).abs() < 1.0e-10)
}

test "ImANM construction" {
  let chain = make_test_chain_im()
  let gammas = @src.gammas_uniform(5, 1.0)
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.ArithmeticMean)
  match imanm {
    Some(m) => {
      assert_eq(m.num_nodes(), 5)
      assert_eq(m.get_gammas().length(), 5)
    }
    None => assert_true(false)
  }
}

test "ImANM mismatched gammas returns None" {
  let chain = make_test_chain_im()
  let gammas = [1.0, 2.0, 3.0]
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0)
  match imanm {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

test "ImANM uniform equals standard ANM" {
  let chain = make_test_chain_im()
  let gammas = @src.gammas_uniform(5, 1.0)
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.ArithmeticMean)
  let anm = @src.ANM::new(chain, cutoff=15.0, gamma=1.0)
  match (imanm, anm) {
    (Some(im), Some(a)) => {
      let im_eigs = im.get_eigenvalues()
      let a_eigs = a.get_eigenvalues()
      for i = 0; i < 15; i = i + 1 {
        assert_true((im_eigs[i] - a_eigs[i]).abs() < 1.0e-6)
      }
    }
    _ => assert_true(false)
  }
}

test "ImANM eigenvalues non-negative" {
  let chain = make_larger_chain_im()
  let gammas = @src.gammas_uniform(10, 1.0)
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0)
  match imanm {
    Some(m) => {
      let eigs = m.get_eigenvalues()
      for e in eigs {
        assert_true(e >= -1.0e-6)
      }
    }
    None => assert_true(false)
  }
}

test "ImANM first 6 eigenvalues near zero" {
  let chain = make_larger_chain_im()
  let gammas = @src.gammas_uniform(10, 1.0)
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0)
  match imanm {
    Some(m) => {
      let eigs = m.get_eigenvalues()
      for i = 0; i < 6; i = i + 1 {
        assert_true(eigs[i].abs() < 1.0e-6)
      }
    }
    None => assert_true(false)
  }
}

test "ImANM hetero different from uniform" {
  let chain = make_larger_chain_im()
  let uniform_g = @src.gammas_uniform(10, 1.0)
  let regions = [(4, 6, 0.1)]
  let hetero_g = @src.gammas_by_region(10, regions, default_gamma=1.0)
  let im_uniform = @src.ImANM::new(chain, uniform_g, cutoff=15.0)
  let im_hetero = @src.ImANM::new(chain, hetero_g, cutoff=15.0)
  match (im_uniform, im_hetero) {
    (Some(u), Some(h)) => {
      let u_eigs = u.get_eigenvalues()
      let h_eigs = h.get_eigenvalues()
      let mut diff = 0.0
      for i = 6; i < 30; i = i + 1 {
        diff = diff + (u_eigs[i] - h_eigs[i]).abs()
      }
      assert_true(diff > 1.0e-6)
    }
    _ => assert_true(false)
  }
}

test "ImANM different combine rules" {
  let chain = make_test_chain_im()
  let gammas = [1.0, 2.0, 0.5, 1.5, 1.0]
  let im_arith = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.ArithmeticMean)
  let im_geom = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.GeometricMean)
  let im_prod = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.Product)
  let im_min = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.Minimum)
  let im_harm = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.HarmonicMean)
  match (im_arith, im_geom, im_prod, im_min, im_harm) {
    (Some(a), Some(g), Some(p), Some(mn), Some(hm)) => {
      assert_eq(a.num_nodes(), 5)
      assert_eq(g.num_nodes(), 5)
      assert_eq(p.num_nodes(), 5)
      assert_eq(mn.num_nodes(), 5)
      assert_eq(hm.num_nodes(), 5)
    }
    _ => assert_true(false)
  }
}

test "ImANM predict_bfactors" {
  let chain = make_larger_chain_im()
  let gammas = @src.gammas_uniform(10, 1.0)
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0)
  match imanm {
    Some(m) => {
      let bf = m.predict_bfactors()
      assert_eq(bf.length(), 10)
      for b in bf {
        assert_true(b >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "ImANM slow_modes" {
  let chain = make_larger_chain_im()
  let gammas = @src.gammas_uniform(10, 1.0)
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0)
  match imanm {
    Some(m) => {
      let modes = m.slow_modes(3)
      assert_eq(modes.length(), 3)
      assert_true(modes[0].1 <= modes[1].1)
      assert_true(modes[1].1 <= modes[2].1)
      for (_, _, vec) in modes {
        assert_eq(vec.length(), 30)
      }
    }
    None => assert_true(false)
  }
}

test "ImANM cross_correlation" {
  let chain = make_larger_chain_im()
  let gammas = @src.gammas_uniform(10, 1.0)
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0)
  match imanm {
    Some(m) => {
      let corr = m.cross_correlation()
      assert_eq(corr.length(), 10)
      assert_eq(corr[0].length(), 10)
      for i = 0; i < 10; i = i + 1 {
        assert_true(corr[i][i] >= -1.0e-10)
      }
      for i = 0; i < 10; i = i + 1 {
        for j = 0; j < 10; j = j + 1 {
          assert_true((corr[i][j] - corr[j][i]).abs() < 1.0e-8)
        }
      }
    }
    None => assert_true(false)
  }
}

test "ImANM get_gammas" {
  let chain = make_test_chain_im()
  let gammas = [1.0, 2.0, 0.5, 1.5, 3.0]
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0)
  match imanm {
    Some(m) => {
      let g = m.get_gammas()
      assert_eq(g.length(), 5)
      assert_true((g[0] - 1.0).abs() < 1.0e-10)
      assert_true((g[1] - 2.0).abs() < 1.0e-10)
      assert_true((g[4] - 3.0).abs() < 1.0e-10)
    }
    None => assert_true(false)
  }
}

test "ImANM hessian symmetric" {
  let chain = make_test_chain_im()
  let gammas = [1.0, 2.0, 0.5, 1.5, 1.0]
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.GeometricMean)
  match imanm {
    Some(m) => {
      let h = m.hessian_matrix()
      let dim = 15
      for i = 0; i < dim; i = i + 1 {
        for j = 0; j < dim; j = j + 1 {
          assert_true((h[i][j] - h[j][i]).abs() < 1.0e-10)
        }
      }
    }
    None => assert_true(false)
  }
}

test "ImANM hessian row sums zero" {
  let chain = make_test_chain_im()
  let gammas = [1.0, 2.0, 0.5, 1.5, 1.0]
  let imanm = @src.ImANM::new(chain, gammas, cutoff=15.0, combine=@src.HarmonicMean)
  match imanm {
    Some(m) => {
      let h = m.hessian_matrix()
      let dim = 15
      for i = 0; i < dim; i = i + 1 {
        let mut row_sum = 0.0
        for j = 0; j < dim; j = j + 1 {
          row_sum = row_sum + h[i][j]
        }
        assert_true(row_sum.abs() < 1.0e-8)
      }
    }
    None => assert_true(false)
  }
}
