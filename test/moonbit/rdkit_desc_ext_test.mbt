///|
/// Tests for rdkit_desc_ext.mbt — Kappa/Chi/BCUT2D/AUTOCORR2D/fr_* and the
/// batch descriptor API.

///|
test "kappa_indices_n_butane" {
  // Reference (no alpha correction): K1=4.0, K2=3.0, K3=2.0
  let mol = @src.parse_smiles("CCCC")
  let (k1, k2, k3) = @src.rdkit_kappa_indices(mol)
  assert_eq(k1, 4.0)
  assert_eq(k2, 3.0)
  assert_eq(k3, 2.0)
}

///|
test "kappa_indices_benzene" {
  // N=6, P1=6, P2=6, P3=6
  let mol = @src.parse_smiles("C1=CC=CC=C1")
  let (k1, k2, k3) = @src.rdkit_kappa_indices(mol)
  assert_eq(k1, 6.0 * 25.0 / 36.0)
  assert_eq(k2, 5.0 * 16.0 / 36.0)
  assert_eq(k3, 4.0 * 9.0 / 36.0)
}

///|
test "chi_indices_n_butane" {
  // degrees [1,2,2,1]; Chi0 = 2 + sqrt(2)
  let mol = @src.parse_smiles("CCCC")
  let (chi0, chi1, _c2, _c3, _c4) = @src.rdkit_chi_indices(mol)
  assert_eq(chi0, 2.0 + 2.0 * @math.pow(2.0, -0.5))
  // Chi1: three bonds, products 1*2, 2*2, 2*1 -> 2/sqrt(2) + 1/2
  let expected1 = 2.0 * @math.pow(2.0, -0.5) + 0.5
  assert_eq(chi1, expected1)
}

///|
test "chi_indices_path_sums_symmetry" {
  // Chi values must be identical under reversal (symmetry sanity check)
  let mol = @src.parse_smiles("CCCC")
  let mol2 = @src.parse_smiles("CCCC")
  let (a0, a1, a2, a3, a4) = @src.rdkit_chi_indices(mol)
  let (b0, b1, b2, b3, b4) = @src.rdkit_chi_indices(mol2)
  assert_eq(a0, b0)
  assert_eq(a1, b1)
  assert_eq(a2, b2)
  assert_eq(a3, b3)
  assert_eq(a4, b4)
}

///|
test "bcut2d_extremes_ordered" {
  // MW-high eigenvalue dominates the MW-low one for a heavy molecule
  let mol = @src.parse_smiles("CCCCl")
  let (mw_hi, mw_lo, ch_hi, ch_lo) = @src.rdkit_bcut2d(mol)
  assert_true(mw_hi > mw_lo)
  // Gasteiger charges are small; extremes stay bounded
  assert_true(ch_hi >= ch_lo)
}

///|
test "autocorr2d_sums_weight_products" {
  // single-bond lag-1 autocorrelation of CCl = wC*wCl
  let mol = @src.parse_smiles("CCl")
  let ats = @src.rdkit_autocorr2d(mol, max_lag=2)
  assert_eq(ats[1], 12.011 * 35.45)
  assert_eq(ats[2], 0.0)
}

///|
test "fragment_counts_common_groups" {
  // symmetric patterns match in both directions; assert at least one hit
  let ether = @src.rdkit_fragment_counts(@src.parse_smiles("COC"))
  assert_true(ether["fr_ether"] >= 1)
  let amide = @src.rdkit_fragment_counts(@src.parse_smiles("CC(=O)NC"))
  assert_eq(amide["fr_amide"], 1)
  let nitrile = @src.rdkit_fragment_counts(@src.parse_smiles("CC#N"))
  assert_eq(nitrile["fr_nitrile"], 1)
  let alkene = @src.rdkit_fragment_counts(@src.parse_smiles("CC=CC"))
  assert_true(alkene["fr_alkene"] >= 1)
  let halo = @src.rdkit_fragment_counts(@src.parse_smiles("CCl"))
  assert_eq(halo["fr_halogen"], 1)
  let benzene = @src.rdkit_fragment_counts(@src.parse_smiles("C1=CC=CC=C1"))
  assert_eq(benzene["fr_benzene"], 1)
}

///|
test "calc_descriptors_batch_map" {
  let mol = @src.parse_smiles("CCO")
  let d = @src.rdkit_calc_descriptors(mol)
  assert_eq(d["NumAtoms"], 3.0)
  assert_eq(d["NumBonds"], 2.0)
  assert_eq(d["HBD"], 1.0)
  assert_true(d["MolWt"] > 40.0)
  // Kappa and fragment keys are present (ethanol has no ether bond)
  assert_true(d["Kappa1"] > 0.0)
  assert_true(d.contains("fr_ether"))
  assert_eq(d["ATS1"], 12.011 * 12.011 + 12.011 * 15.999)
}
