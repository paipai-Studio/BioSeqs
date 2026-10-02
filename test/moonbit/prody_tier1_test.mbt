///|
/// Tier 1 ProDy/OpenStructure tests: GNM, ANM, ensemble PCA,
/// clash detection, bond inference.

///|
fn make_test_chain() -> @src.Chain {
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
fn make_clashing_atoms() -> Array[@src.Atom] {
  let atoms : Array[@src.Atom] = Array(capacity=0)
  atoms.push(@src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="CB", coord=@src.Vector3::new(0.5, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="CG", coord=@src.Vector3::new(1.5, 0.0, 0.0), element="C"))
  atoms
}

///|
fn make_bonded_atoms() -> Array[@src.Atom] {
  let atoms : Array[@src.Atom] = Array(capacity=0)
  atoms.push(@src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="C2", coord=@src.Vector3::new(1.54, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="O1", coord=@src.Vector3::new(1.54, 1.43, 0.0), element="O"))
  atoms.push(@src.Atom::new(name="N1", coord=@src.Vector3::new(10.0, 10.0, 10.0), element="N"))
  atoms
}

// GNM tests

test "gnm_construction" {
  let chain = make_test_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => assert_eq(gnm.num_nodes(), 5)
    None => assert_true(false)
  }
}

test "gnm_kirchhoff_shape" {
  let chain = make_test_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let mat = gnm.kirchhoff_matrix()
      assert_eq(mat.length(), 5)
      assert_eq(mat[0].length(), 5)
    }
    None => assert_true(false)
  }
}

test "gnm_kirchhoff_row_sum_zero" {
  let chain = make_test_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let mat = gnm.kirchhoff_matrix()
      for i = 0; i < 5; i = i + 1 {
        let mut sum = 0.0
        for j = 0; j < 5; j = j + 1 { sum = sum + mat[i][j] }
        assert_true(sum.abs() < 1.0e-10)
      }
    }
    None => assert_true(false)
  }
}

test "gnm_eigenvalues" {
  let chain = make_test_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let vals = gnm.get_eigenvalues()
      assert_true(vals[0].abs() < 1.0e-6)
      for i = 0; i < 5; i = i + 1 { assert_true(vals[i] >= -1.0e-6) }
    }
    None => assert_true(false)
  }
}

test "gnm_bfactors" {
  let chain = make_test_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => {
      let bf = gnm.predict_bfactors()
      assert_eq(bf.length(), 5)
      for i = 0; i < 5; i = i + 1 { assert_true(bf[i] >= 0.0) }
    }
    None => assert_true(false)
  }
}

test "gnm_slow_modes" {
  let chain = make_test_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => assert_eq(gnm.slow_modes(2).length(), 2)
    None => assert_true(false)
  }
}

test "gnm_num_contacts" {
  let chain = make_test_chain()
  match @src.GNM::new(chain) {
    Some(gnm) => assert_true(gnm.num_contacts() > 0)
    None => assert_true(false)
  }
}

test "gnm_kirchhoff_function" {
  let coords = [@src.Vector3::new(0.0, 0.0, 0.0), @src.Vector3::new(3.0, 0.0, 0.0), @src.Vector3::new(6.0, 0.0, 0.0)]
  let mat = @src.gnm_kirchhoff(coords)
  assert_eq(mat.length(), 3)
  for i = 0; i < 3; i = i + 1 {
    let mut sum = 0.0
    for j = 0; j < 3; j = j + 1 { sum = sum + mat[i][j] }
    assert_true(sum.abs() < 1.0e-10)
  }
}

// ANM tests

test "anm_construction" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => assert_eq(anm.num_nodes(), 5)
    None => assert_true(false)
  }
}

test "anm_hessian_shape" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let h = anm.hessian_matrix()
      assert_eq(h.length(), 15)
      assert_eq(h[0].length(), 15)
    }
    None => assert_true(false)
  }
}

test "anm_hessian_symmetric" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let h = anm.hessian_matrix()
      for i = 0; i < 15; i = i + 1 {
        for j = 0; j < 15; j = j + 1 { assert_true((h[i][j] - h[j][i]).abs() < 1.0e-10) }
      }
    }
    None => assert_true(false)
  }
}

test "anm_hessian_row_sum_zero" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let h = anm.hessian_matrix()
      for i = 0; i < 15; i = i + 1 {
        let mut sum = 0.0
        for j = 0; j < 15; j = j + 1 { sum = sum + h[i][j] }
        assert_true(sum.abs() < 1.0e-10)
      }
    }
    None => assert_true(false)
  }
}

test "anm_eigenvalues" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let vals = anm.get_eigenvalues()
      assert_eq(vals.length(), 15)
      for i = 0; i < 6; i = i + 1 { assert_true(vals[i].abs() < 1.0) }
    }
    None => assert_true(false)
  }
}

test "anm_bfactors" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let bf = anm.predict_bfactors()
      assert_eq(bf.length(), 5)
      for i = 0; i < 5; i = i + 1 { assert_true(bf[i] >= 0.0) }
    }
    None => assert_true(false)
  }
}

test "anm_slow_modes" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => assert_true(anm.slow_modes(3).length() >= 1)
    None => assert_true(false)
  }
}

test "anm_cross_correlation" {
  let chain = make_test_chain()
  match @src.ANM::new(chain) {
    Some(anm) => {
      let corr = anm.cross_correlation()
      assert_eq(corr.length(), 5)
      // Diagonal values are normalized fluctuations in [0, 1]
      let mut max_diag = 0.0
      for i = 0; i < 5; i = i + 1 {
        assert_true(corr[i][i] >= -0.01)
        assert_true(corr[i][i] <= 1.01)
        if corr[i][i] > max_diag { max_diag = corr[i][i] }
      }
      // Max diagonal should be ~1.0 (normalized)
      assert_true((max_diag - 1.0).abs() < 0.01)
    }
    None => assert_true(false)
  }
}

///|
fn make_shifted_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let coords = [
    (1.5, 0.2, 0.1),
    (0.41, 1.05, 1.6),
    (-0.71, 0.69, 3.1),
    (-0.71, -0.49, 4.6),
    (0.41, -0.85, 6.1),
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

// Ensemble PCA tests

test "ensemble_construction" {
  let chain1 = make_test_chain()
  let chain2 = make_test_chain()
  match @src.Ensemble::new([chain1, chain2]) {
    Some(ens) => {
      assert_eq(ens.num_conformations(), 2)
      assert_eq(ens.num_residues(), 5)
    }
    None => assert_true(false)
  }
}

test "ensemble_superimpose" {
  let chain1 = make_test_chain()
  let chain2 = make_test_chain()
  match @src.Ensemble::new([chain1, chain2]) {
    Some(ens) => {
      let aligned = ens.superimpose()
      assert_eq(aligned.get_aligned(0).length(), 5)
    }
    None => assert_true(false)
  }
}

test "ensemble_pca" {
  let chain1 = make_test_chain()
  let chain2 = make_shifted_chain()
  match @src.Ensemble::new([chain1, chain2]) {
    Some(ens) => {
      match ens.pca() {
        Some(pca) => {
          let varr = pca.get_variance_explained()
          let mut total = 0.0
          for v in varr { total = total + v }
          assert_true((total - 1.0).abs() < 0.01)
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "ensemble_pca_eigenvalues_descending" {
  let chain1 = make_test_chain()
  let chain2 = make_shifted_chain()
  match @src.Ensemble::new([chain1, chain2]) {
    Some(ens) => {
      match ens.pca() {
        Some(pca) => {
          let vals = pca.get_eigenvalues()
          for i = 1; i < vals.length(); i = i + 1 {
            assert_true(vals[i] <= vals[i - 1] + 1.0e-10)
          }
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

test "ensemble_pca_project_mean" {
  let chain1 = make_test_chain()
  let chain2 = make_shifted_chain()
  match @src.Ensemble::new([chain1, chain2]) {
    Some(ens) => {
      match ens.pca() {
        Some(pca) => {
          let mean = pca.mean_structure()
          let projections = pca.project(mean)
          for p in projections { assert_true(p.abs() < 1.0) }
        }
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
}

// Clash detection tests

test "clash_no_clash" {
  let atoms : Array[@src.Atom] = Array(capacity=0)
  atoms.push(@src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="C2", coord=@src.Vector3::new(5.0, 0.0, 0.0), element="C"))
  let result = @src.detect_clashes(atoms)
  assert_eq(result.num_clashes, 0)
}

test "clash_detected" {
  let atoms = make_clashing_atoms()
  let result = @src.detect_clashes(atoms)
  assert_true(result.num_clashes > 0)
}

test "clash_score_positive" {
  let atoms = make_clashing_atoms()
  assert_true(@src.clash_score(atoms) > 0.0)
}

test "clash_count" {
  let atoms = make_clashing_atoms()
  assert_true(@src.count_clashes(atoms) > 0)
}

test "clash_severe" {
  let atoms : Array[@src.Atom] = Array(capacity=0)
  atoms.push(@src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="C2", coord=@src.Vector3::new(0.3, 0.0, 0.0), element="C"))
  assert_true(@src.has_severe_clashes(atoms))
}

test "clash_not_severe" {
  let atoms : Array[@src.Atom] = Array(capacity=0)
  atoms.push(@src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="C2", coord=@src.Vector3::new(5.0, 0.0, 0.0), element="C"))
  assert_false(@src.has_severe_clashes(atoms))
}

// Bond inference tests

test "infer_bonds_basic" {
  let atoms = make_bonded_atoms()
  let bonds = @src.infer_bonds(atoms)
  assert_true(bonds.length() >= 2)
}

test "infer_bonds_count" {
  let atoms = make_bonded_atoms()
  assert_true(@src.count_inferred_bonds(atoms) >= 2)
}

test "infer_bonds_pairs" {
  let atoms = make_bonded_atoms()
  let pairs = @src.inferred_bonds_pairs(atoms)
  assert_true(pairs.length() >= 2)
}

test "infer_bonds_no_false_positives" {
  let atoms : Array[@src.Atom] = Array(capacity=0)
  atoms.push(@src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="N1", coord=@src.Vector3::new(10.0, 0.0, 0.0), element="N"))
  assert_eq(@src.infer_bonds(atoms).length(), 0)
}

test "infer_bonds_cc_distance" {
  let atoms : Array[@src.Atom] = Array(capacity=0)
  atoms.push(@src.Atom::new(name="C1", coord=@src.Vector3::new(0.0, 0.0, 0.0), element="C"))
  atoms.push(@src.Atom::new(name="C2", coord=@src.Vector3::new(1.54, 0.0, 0.0), element="C"))
  let bonds = @src.infer_bonds(atoms)
  assert_eq(bonds.length(), 1)
  assert_true((bonds[0].distance - 1.54).abs() < 0.01)
}

test "infer_bonds_report" {
  let atoms = make_bonded_atoms()
  let bonds = @src.infer_bonds(atoms)
  let report = @src.bond_report(bonds)
  assert_true(report.contains("Inferred"))
}

test "clash_report_string" {
  let atoms = make_clashing_atoms()
  let result = @src.detect_clashes(atoms)
  let report = @src.clash_report(result)
  assert_true(report.contains("Clash"))
}

// Stereochemistry validation

test "stereo_validation_basic" {
  let chain = make_test_chain()
  let result = @src.validate_stereochemistry(chain)
  assert_eq(result.num_chiral_centers, 5)
}
