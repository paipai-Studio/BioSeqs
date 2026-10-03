///|
/// Tests for ENM parameter fitting.

///|
/// Build a test chain: 8 residues in a helix-like arrangement.
fn make_fit_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 8; r = r + 1 {
    let angle = r.to_double() * 1.75
    let x = 1.5 * @math.cos(angle)
    let y = 1.5 * @math.sin(angle)
    let z = r.to_double() * 2.3
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, y, z), resname="ALA", chainid='A', resseq=r + 1, element="C")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
/// Generate synthetic B-factors from ANM with known parameters.
fn make_synthetic_bfactors(chain : @src.Chain, cutoff : Double, gamma : Double) -> Array[Double] {
  match @src.ANM::new(chain, cutoff~, gamma~) {
    Some(anm) => @src.ANM::predict_bfactors(anm)
    None => Array::make(8, 0.0)
  }
}

// ===== fit_gamma =====

test "fit_gamma_basic" {
  let predicted = [1.0, 2.0, 3.0, 4.0, 5.0]
  let observed = [2.0, 4.0, 6.0, 8.0, 10.0]
  let gamma = @src.fit_gamma(predicted, observed)
  // observed = 2.0 * predicted, so gamma = 2.0
  assert_true((gamma - 2.0).abs() < 0.001)
}

test "fit_gamma_identity" {
  let predicted = [1.0, 2.0, 3.0, 4.0, 5.0]
  let observed = [1.0, 2.0, 3.0, 4.0, 5.0]
  let gamma = @src.fit_gamma(predicted, observed)
  assert_true((gamma - 1.0).abs() < 0.001)
}

test "fit_gamma_zero" {
  let predicted = [0.0, 0.0, 0.0]
  let observed = [1.0, 2.0, 3.0]
  let gamma = @src.fit_gamma(predicted, observed)
  assert_true((gamma - 1.0).abs() < 0.001)
}

test "fit_gamma_recovers_known" {
  let chain = make_fit_chain()
  let true_gamma = 3.5
  let true_cutoff = 10.0
  let observed = make_synthetic_bfactors(chain, true_cutoff, true_gamma)
  // Get predicted at gamma=1
  match @src.ANM::new(chain, cutoff=true_cutoff, gamma=1.0) {
    Some(anm) => {
      let pred = @src.ANM::predict_bfactors(anm)
      let fitted_gamma = @src.fit_gamma(pred, observed)
      assert_true((fitted_gamma * true_gamma - 1.0).abs() < 0.01)
    }
    None => assert_true(false)
  }
}

// ===== fit_quality =====

test "fit_quality_perfect" {
  let predicted = [1.0, 2.0, 3.0, 4.0, 5.0]
  let observed = [1.0, 2.0, 3.0, 4.0, 5.0]
  let (corr, rms) = @src.fit_quality(predicted, observed)
  assert_true((corr - 1.0).abs() < 0.001)
  assert_true(rms < 0.001)
}

test "fit_quality_anticorrelated" {
  let predicted = [1.0, 2.0, 3.0, 4.0, 5.0]
  let observed = [5.0, 4.0, 3.0, 2.0, 1.0]
  let (corr, _) = @src.fit_quality(predicted, observed)
  assert_true((corr + 1.0).abs() < 0.001)
}

test "fit_quality_rms" {
  let predicted = [1.0, 2.0, 3.0]
  let observed = [2.0, 3.0, 4.0]
  let (_, rms) = @src.fit_quality(predicted, observed)
  // diff = [1,1,1], rms = 1.0
  assert_true((rms - 1.0).abs() < 0.001)
}

// ===== scale_bfactors =====

test "scale_bfactors" {
  let pred = [1.0, 2.0, 3.0]
  let scaled = @src.scale_bfactors(pred, 2.5)
  assert_true((scaled[0] - 2.5).abs() < 0.001)
  assert_true((scaled[1] - 5.0).abs() < 0.001)
  assert_true((scaled[2] - 7.5).abs() < 0.001)
}

// ===== fit_cutoff =====

test "fit_cutoff_finds_known" {
  let chain = make_fit_chain()
  let true_cutoff = 8.0
  let true_gamma = 2.0
  let observed = make_synthetic_bfactors(chain, true_cutoff, true_gamma)
  let cutoffs = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0]
  match @src.fit_cutoff(chain, observed, cutoffs) {
    Some(result) => {
      // Should find cutoff close to 8.0
      assert_true((@src.ENMFitResult::get_cutoff(result) - true_cutoff).abs() < 2.0)
      assert_true(@src.ENMFitResult::get_correlation(result) > 0.9)
    }
    None => assert_true(false)
  }
}

test "fit_cutoff_returns_best" {
  let chain = make_fit_chain()
  let observed = make_synthetic_bfactors(chain, 10.0, 1.0)
  let cutoffs = [7.0, 10.0, 13.0]
  match @src.fit_cutoff(chain, observed, cutoffs) {
    Some(result) => {
      assert_true(@src.ENMFitResult::get_correlation(result) > 0.0)
    }
    None => assert_true(false)
  }
}

test "fit_cutoff_empty" {
  let chain = make_fit_chain()
  let observed = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  match @src.fit_cutoff(chain, observed, []) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ===== fit_cutoff_gamma =====

test "fit_cutoff_gamma_all" {
  let chain = make_fit_chain()
  let observed = make_synthetic_bfactors(chain, 10.0, 2.0)
  let cutoffs = [7.0, 10.0, 13.0]
  let results = @src.fit_cutoff_gamma(chain, observed, cutoffs)
  assert_true(results.length() >= 1)
  // Best should be at cutoff=10
  match @src.best_fit(results) {
    Some(best) => {
      assert_true(@src.ENMFitResult::get_correlation(best) > 0.5)
    }
    None => assert_true(false)
  }
}

test "best_fit_empty" {
  match @src.best_fit([]) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ===== fit_cutoff_gnm =====

test "fit_cutoff_gnm" {
  let chain = make_fit_chain()
  // Generate B-factors from GNM
  match @src.GNM::new(chain, cutoff=10.0, gamma=1.0) {
    Some(gnm) => {
      let observed = @src.GNM::predict_bfactors(gnm)
      let cutoffs = [7.0, 10.0, 13.0]
      match @src.fit_cutoff_gnm(chain, observed, cutoffs) {
        Some(result) => {
          assert_true(@src.ENMFitResult::get_correlation(result) > 0.5)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

// ===== Mode matching =====

test "match_modes_identity" {
  let modes = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let overlap = @src.match_modes(modes, modes)
  // Identity: diagonal should be 1.0
  assert_true((overlap[0][0] - 1.0).abs() < 0.001)
  assert_true((overlap[1][1] - 1.0).abs() < 0.001)
  assert_true(overlap[0][1] < 0.001)
}

test "match_modes_orthogonal" {
  let modes_a = [[1.0, 0.0, 0.0]]
  let modes_b = [[0.0, 1.0, 0.0]]
  let overlap = @src.match_modes(modes_a, modes_b)
  assert_true(overlap[0][0] < 0.001)
}

test "enm_effective_overlap" {
  let modes = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let overlap = @src.match_modes(modes, modes)
  let eff = @src.enm_effective_overlap(overlap)
  assert_true((eff - 1.0).abs() < 0.001)
}

test "mode_matching_score" {
  let modes = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let overlap = @src.match_modes(modes, modes)
  let score = @src.mode_matching_score(overlap)
  assert_true((score - 1.0).abs() < 0.001)
}

test "enm_cumulative_overlap" {
  let modes_a = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let modes_b = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let cum = @src.enm_cumulative_overlap(modes_a, modes_b)
  assert_eq(cum.length(), 2)
  assert_true((cum[0] - 1.0).abs() < 0.001)
  assert_true((cum[1] - 1.0).abs() < 0.001)
}

test "enm_cumulative_overlap_partial" {
  let modes_a = [
    [0.6, 0.8, 0.0],
  ]
  let modes_b = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  let cum = @src.enm_cumulative_overlap(modes_a, modes_b)
  // Mode a is in span of modes_b, so cumulative overlap = 1
  assert_true((cum[0] - 1.0).abs() < 0.001)
}

// ===== Cross-validation =====

test "loo_cross_validate" {
  let predicted = [1.0, 2.0, 3.0, 4.0, 5.0]
  let observed = [2.0, 4.0, 6.0, 8.0, 10.0]
  let residuals = @src.loo_cross_validate(predicted, observed)
  assert_eq(residuals.length(), 5)
  // Since observed = 2*predicted, LOO should recover gamma=2, residuals ~0
  for i = 0; i < 5; i = i + 1 {
    assert_true(residuals[i].abs() < 0.001)
  }
}

test "loo_cross_validate_noisy" {
  let predicted = [1.0, 2.0, 3.0, 4.0, 5.0]
  let observed = [2.1, 3.9, 6.2, 7.8, 10.1]
  let residuals = @src.loo_cross_validate(predicted, observed)
  assert_eq(residuals.length(), 5)
  // Residuals should be small
  for i = 0; i < 5; i = i + 1 {
    assert_true(residuals[i].abs() < 1.0)
  }
}

test "kfold_cross_validate" {
  let predicted = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let observed = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
  let (corr, rms) = @src.kfold_cross_validate(predicted, observed, 3)
  // Perfect linear relationship → high correlation
  assert_true(corr > 0.99)
  assert_true(rms < 0.1)
}

// ===== fit_beta =====

test "fit_beta_finds_known" {
  let chain = make_fit_chain()
  // Generate B-factors from betaANM with beta=2
  match @src.BetaANM::new(chain, beta=2.0) {
    Some(anm) => {
      let observed = @src.BetaANM::predict_bfactors(anm)
      let betas = [0.0, 1.0, 2.0, 3.0, 4.0]
      match @src.fit_beta(chain, observed, betas) {
        Some((best_beta, corr, _)) => {
          assert_true((best_beta - 2.0).abs() < 0.5)
          assert_true(corr > 0.5)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

// ===== compare_anm_gnm =====

test "compare_anm_gnm" {
  let chain = make_fit_chain()
  let observed = make_synthetic_bfactors(chain, 10.0, 1.0)
  let (anm_corr, gnm_corr) = @src.compare_anm_gnm(chain, observed, 10.0)
  // ANM should correlate well with its own B-factors
  assert_true(anm_corr > 0.5)
}

// ===== correlation_profile =====

test "correlation_profile" {
  let chain = make_fit_chain()
  let observed = make_synthetic_bfactors(chain, 10.0, 1.0)
  let cutoffs = [7.0, 10.0, 13.0]
  let profile = @src.correlation_profile(chain, observed, cutoffs)
  assert_true(profile.length() >= 1)
  for (cutoff, corr) in profile {
    assert_true(corr >= -1.0 && corr <= 1.0)
    assert_true(cutoff > 0.0)
  }
}

// ===== ENMFitResult accessors =====

test "enm_fit_result_accessors" {
  let chain = make_fit_chain()
  let observed = make_synthetic_bfactors(chain, 10.0, 1.0)
  match @src.fit_cutoff(chain, observed, [10.0]) {
    Some(result) => {
      assert_true(@src.ENMFitResult::get_cutoff(result) > 0.0)
      assert_true(@src.ENMFitResult::get_gamma(result) > 0.0)
      assert_true(@src.ENMFitResult::get_correlation(result) >= -1.0)
      assert_true(@src.ENMFitResult::get_rms_error(result) >= 0.0)
      let pred = @src.ENMFitResult::get_predicted_bfactors(result)
      assert_eq(pred.length(), 8)
    }
    None => assert_true(false)
  }
}
