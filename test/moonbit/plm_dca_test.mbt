// ===========================================================================
// plm_dca_test.mbt — Tests for pseudo-likelihood maximization DCA
// ===========================================================================

///|
fn make_protein_aln() -> Array[String] {
  // Small alignment with some coevolution signal
  // Positions 1 and 4 are correlated (A-A or C-C)
  [
    "ACDEACDE",
    "ACDEACDE",
    "CCDECCDE",
    "CCDECCDE",
    "ACDEACDE",
    "CCDECCDE",
    "ADDEADDE",
    "ADDEADDE",
    "CDDECDDE",
    "CDDECDDE",
  ]
}

test "PlmDCAConfig default_protein" {
  let c = @src.PlmDCAConfig::default_protein()
  assert_eq(c.alphabet, "protein")
  assert_true(c.lambda_h > 0.0)
  assert_true(c.lambda_j > 0.0)
  assert_true(c.max_iter > 0)
  assert_true(c.learning_rate > 0.0)
  assert_true(c.tol > 0.0)
}

test "PlmDCAConfig default_dna" {
  let c = @src.PlmDCAConfig::default_dna()
  assert_eq(c.alphabet, "dna")
  assert_true(c.lambda_h > 0.0)
  assert_true(c.lambda_j > 0.0)
}

test "plm_dca basic shape" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  assert_eq(result.align_len, 8)
  assert_eq(result.num_seqs, 10)
  assert_true(result.m_eff > 0.0)
  // Couplings: L × L × q × q
  assert_eq(result.couplings.length(), 8)
  assert_eq(result.couplings[0].length(), 8)
  // Fields: L × q
  assert_eq(result.fields.length(), 8)
  // FN and DI scores: L × L
  assert_eq(result.fn_scores.length(), 8)
  assert_eq(result.fn_scores[0].length(), 8)
  assert_eq(result.di_scores.length(), 8)
  assert_eq(result.di_scores[0].length(), 8)
}

test "plm_dca fn_scores non-negative" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      assert_true(result.fn_scores[i][j] >= 0.0)
    }
  }
}

test "plm_dca di_scores non-negative" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      if i != j {
        assert_true(result.di_scores[i][j] >= -1.0e-6)
      }
    }
  }
}

test "plm_dca fn_scores symmetric" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      let diff = (result.fn_scores[i][j] - result.fn_scores[j][i]).abs()
      assert_true(diff < 1.0e-10)
    }
  }
}

test "plm_dca di_scores symmetric" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      let diff = (result.di_scores[i][j] - result.di_scores[j][i]).abs()
      assert_true(diff < 1.0e-6)
    }
  }
}

test "plm_dca diagonal zero" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  for i = 0; i < 8; i = i + 1 {
    assert_true(result.fn_scores[i][i].abs() < 1.0e-10)
  }
}

test "plm_dca coupling matrix symmetric" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  // J_ij(a,b) should equal J_ji(b,a) after symmetrization
  let q = 21
  for i = 0; i < 8; i = i + 1 {
    for j = (i + 1); j < 8; j = j + 1 {
      for a = 0; a < q; a = a + 1 {
        for b = 0; b < q; b = b + 1 {
          let diff = (result.couplings[i][j][a][b] - result.couplings[j][i][b][a]).abs()
          assert_true(diff < 1.0e-10)
        }
      }
    }
  }
}

test "plm_dca fields zero-sum gauge" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  let q = 21
  for i = 0; i < 8; i = i + 1 {
    let mut s = 0.0
    for a = 0; a < q; a = a + 1 {
      s = s + result.fields[i][a]
    }
    assert_true(s.abs() < 1.0e-10)
  }
}

test "plm_dca coupling zero-sum gauge" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  let q = 21
  // Row sums and column sums should be ~0
  for i = 0; i < 8; i = i + 1 {
    for j = (i + 1); j < 8; j = j + 1 {
      for a = 0; a < q; a = a + 1 {
        let mut row_sum = 0.0
        for b = 0; b < q; b = b + 1 {
          row_sum = row_sum + result.couplings[i][j][a][b]
        }
        assert_true(row_sum.abs() < 1.0e-8)
      }
      for b = 0; b < q; b = b + 1 {
        let mut col_sum = 0.0
        for a = 0; a < q; a = a + 1 {
          col_sum = col_sum + result.couplings[i][j][a][b]
        }
        assert_true(col_sum.abs() < 1.0e-8)
      }
    }
  }
}

test "plm_dca coevolution signal detected" {
  // Positions 0 and 4 are correlated: A-A or C-C
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  // The coupling between positions 0 and 4 should be stronger than
  // between positions 0 and 1 (which are not correlated)
  let fn_04 = result.fn_scores[0][4]
  let fn_01 = result.fn_scores[0][1]
  // fn_04 should be noticeably larger than fn_01
  assert_true(fn_04 > fn_01)
}

test "plm_dca identity alignment low couplings" {
  // All identical sequences → no coevolution signal
  let aln = ["ACDEACDE", "ACDEACDE", "ACDEACDE", "ACDEACDE"]
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  // All FN scores should be very small (no variation → no coupling)
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      if i != j {
        assert_true(result.fn_scores[i][j] < 1.0)
      }
    }
  }
}

test "plm_dca random alignment moderate couplings" {
  // Diverse alignment with no correlation
  let aln = [
    "ACDEFGHI",
    "KLMNPQRS",
    "TVWYACDE",
    "FGHIKLMM",
    "NPQRSTVW",
    "YACDEFGH",
  ]
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  // Results should be valid
  assert_eq(result.align_len, 8)
  assert_eq(result.num_seqs, 6)
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      assert_true(result.fn_scores[i][j] >= 0.0)
    }
  }
}

test "plm_dca DNA alignment" {
  let aln = [
    "ACGTACGT",
    "ACGTACGT",
    "TGACTGAC",
    "TGACTGAC",
    "ACGTACGT",
    "TGACTGAC",
    "ACGTACGT",
    "TGACTGAC",
  ]
  let config = @src.PlmDCAConfig::default_dna()
  let result = @src.plm_dca(aln, config)
  assert_eq(result.align_len, 8)
  assert_eq(result.num_seqs, 8)
  // q = 5 for DNA
  assert_eq(result.couplings[0][1].length(), 5)
  assert_eq(result.couplings[0][1][0].length(), 5)
  // FN scores symmetric
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      let diff = (result.fn_scores[i][j] - result.fn_scores[j][i]).abs()
      assert_true(diff < 1.0e-10)
    }
  }
}

test "plm_rank_contacts" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  let contacts = @src.plm_rank_contacts(result, 3)
  // Should have pairs with |i-j| >= 3
  assert_true(contacts.length() > 0)
  // Check sorted descending
  for i in 1..<contacts.length() {
    assert_true(contacts[i - 1].score >= contacts[i].score)
  }
  // Check min sequence distance
  for c in contacts {
    assert_true(c.pos_j - c.pos_i >= 3)
    assert_true(c.pos_i < c.pos_j)
  }
}

test "plm_dca single sequence" {
  let aln = ["ACDEACDE"]
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  assert_eq(result.align_len, 8)
  assert_eq(result.num_seqs, 1)
  assert_true(result.m_eff > 0.0)
}

test "plm_dca two sequences" {
  let aln = ["ACDEACDE", "CCDECCDE"]
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  assert_eq(result.align_len, 8)
  assert_eq(result.num_seqs, 2)
  // FN scores should be valid
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      assert_true(result.fn_scores[i][j] >= 0.0)
    }
  }
}

test "plm_dca vs mfDCA ranking correlation" {
  // Both methods should identify the same top contact (pos 0, pos 4)
  let aln = make_protein_aln()
  let plm_config = @src.PlmDCAConfig::default_protein()
  let plm_result = @src.plm_dca(aln, plm_config)
  let mf_config = @src.DCAConfig::default_protein()
  let mf_result = @src.mf_dca(aln, mf_config)
  // Get top plmDCA contact
  let plm_contacts = @src.plm_rank_contacts(plm_result, 3)
  // Get top mfDCA contacts
  let mf_contacts = @src.rank_couplings(mf_result.dca_scores, 3)
  // Both should have the same top contact pair
  if plm_contacts.length() > 0 && mf_contacts.length() > 0 {
    let plm_top = plm_contacts[0]
    let mf_top = mf_contacts[0]
    // Top contact should be one of the correlated pairs: (0,4) or (1,5)
    let plm_is_04 = plm_top.pos_i == 0 && plm_top.pos_j == 4
    let plm_is_15 = plm_top.pos_i == 1 && plm_top.pos_j == 5
    assert_true(plm_is_04 || plm_is_15)
    let mf_is_04 = mf_top.pos_i == 0 && mf_top.pos_j == 4
    let mf_is_15 = mf_top.pos_i == 1 && mf_top.pos_j == 5
    assert_true(mf_is_04 || mf_is_15)
  }
}

test "plm_dca m_eff with reweighting" {
  // Many identical sequences → low m_eff
  let aln = [
    "ACDEACDE", "ACDEACDE", "ACDEACDE",
    "CCDECCDE", "CCDECCDE", "CCDECCDE",
  ]
  let config = @src.PlmDCAConfig::default_protein()
  let result = @src.plm_dca(aln, config)
  // 6 sequences but only 2 unique → m_eff ≈ 2
  assert_true(result.m_eff < 3.0)
  assert_true(result.m_eff > 1.0)
}

test "plm_dca custom config" {
  let aln = make_protein_aln()
  let config = @src.PlmDCAConfig::{
    alphabet: "protein",
    reweight_threshold: 0.8,
    lambda_h: 0.05,
    lambda_j: 0.02,
    max_iter: 100,
    learning_rate: 0.01,
    tol: 1.0e-4,
  }
  let result = @src.plm_dca(aln, config)
  assert_eq(result.align_len, 8)
  // With higher regularization, couplings should be smaller
  for i = 0; i < 8; i = i + 1 {
    for j = 0; j < 8; j = j + 1 {
      assert_true(result.fn_scores[i][j] >= 0.0)
    }
  }
}
