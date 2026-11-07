///|
/// Tests for rdkit_forcefield.mbt — UFF-style energy, optimization and
/// multi-conformer sampling.

///|
test "ff_energy_optimization_decreases" {
  let mol = @src.parse_smiles("CCO")
  let raw = @src.rdkit_generate_conformer("CCO", seed=7)
  let e_before = @src.rdkit_ff_energy(mol, raw.coords)
  let (e_after, coords) = @src.rdkit_ff_optimize(mol, raw.coords, max_iters=150)
  assert_true(e_after < e_before)
  assert_true(coords.length() == 3)
  // determinism: same input gives same output
  let (e_again, _c2) = @src.rdkit_ff_optimize(mol, raw.coords, max_iters=150)
  assert_eq(e_after, e_again)
}

///|
test "ff_bond_lengths_near_ideal_after_optimization" {
  let mol = @src.parse_smiles("CCO")
  let raw = @src.rdkit_generate_conformer("CCO", seed=3)
  let (_e, coords) = @src.rdkit_ff_optimize(mol, raw.coords, max_iters=200)
  // bond 0-1 is a C-C single bond (ideal ~1.52), bond 1-2 C-O (~1.43)
  let d01 = @src.ff_bond_length(coords, 0, 1)
  let d12 = @src.ff_bond_length(coords, 1, 2)
  assert_true(d01 > 1.3 && d01 < 1.75)
  assert_true(d12 > 1.2 && d12 < 1.65)
}

///|
test "conformer_rmsd_zero_for_identical" {
  let conf = @src.rdkit_generate_conformer("CCO", seed=5)
  // identical geometry -> numerically zero after superposition
  assert_true(@src.rdkit_conformer_rmsd(conf, conf) < 0.000001)
}

///|
test "conformer_rmsd_invariant_under_translation" {
  let conf = @src.rdkit_generate_conformer("CCO", seed=5)
  let shifted = @src.RdkitConformer::from_coords(conf.atoms, [
    (conf.coords[0].0 + 5.0, conf.coords[0].1 - 2.0, conf.coords[0].2 + 1.0),
    (conf.coords[1].0 + 5.0, conf.coords[1].1 - 2.0, conf.coords[1].2 + 1.0),
    (conf.coords[2].0 + 5.0, conf.coords[2].1 - 2.0, conf.coords[2].2 + 1.0),
  ])
  assert_true(@src.rdkit_conformer_rmsd(conf, shifted) < 0.000001)
}

///|
test "embed_multiple_conformers_sorted_and_deduped" {
  let confs = @src.rdkit_embed_multiple_conformers(
    "CCO",
    n_confs=4,
    rms_threshold=0.3,
    seed=11,
  )
  assert_true(confs.length() >= 1)
  assert_true(confs.length() <= 4)
  // energies sorted ascending
  for i = 1; i < confs.length(); i = i + 1 {
    assert_true(confs[i].energy >= confs[i - 1].energy)
  }
}
