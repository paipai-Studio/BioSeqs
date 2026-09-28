///|
/// Tests for SMILES parsing + Morgan fingerprints (rdkit.mbt).

fn fp_equal(a : Array[Int], b : Array[Int]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let mut i = 0
  while i < a.length() {
    if a[i] != b[i] {
      return false
    }
    i = i + 1
  }
  true
}

test "parse_smiles_ethanol" {
  let mol = @src.parse_smiles("CCO")
  assert_eq(mol.num_atoms(), 3)
  assert_eq(mol.num_bonds(), 2)
  assert_eq(mol.atom_symbols(), ["C", "C", "O"])
  // 两个单键
  assert_eq(mol.bond_triples()[0], (0, 1, 1))
  assert_eq(mol.bond_triples()[1], (1, 2, 1))
}

test "parse_smiles_double_bond" {
  let mol = @src.parse_smiles("CC=C")
  assert_eq(mol.num_atoms(), 3)
  // 第二个键是双键
  assert_eq(mol.bond_triples()[1].2, 2)
}

test "parse_smiles_ring_and_aromatic" {
  // 环己烷：6 碳 6 键（环闭合）
  let cyc = @src.parse_smiles("C1CCCCC1")
  assert_eq(cyc.num_atoms(), 6)
  assert_eq(cyc.num_bonds(), 6)

  // 苯：6 个芳香碳
  let benz = @src.parse_smiles("c1ccccc1")
  assert_eq(benz.num_atoms(), 6)
  assert_eq(benz.num_bonds(), 6)
}

test "morgan_fingerprint_length_and_stability" {
  let m1 = @src.parse_smiles("CCO")
  let m2 = @src.parse_smiles("CCO")
  let fp1 = @src.morgan_fingerprint(m1, 2)
  let fp2 = @src.morgan_fingerprint(m2, 2)
  // n 原子 × (radius+1) 个哈希
  assert_eq(fp1.length(), 3 * 3)
  // 相同分子 → 相同指纹
  assert_true(fp_equal(fp1, fp2))
}

test "morgan_fingerprint_distinguishes" {
  // 乙烷 CC vs 甲醇片段 CO：radius 1 的环境哈希不同
  let cc = @src.parse_smiles("CC")
  let co = @src.parse_smiles("CO")
  let fp_cc = @src.morgan_fingerprint(cc, 1)
  let fp_co = @src.morgan_fingerprint(co, 1)
  assert_true(not(fp_equal(fp_cc, fp_co)))
}
