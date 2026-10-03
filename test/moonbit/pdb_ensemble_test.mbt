///|
/// Tests for PDBEnsemble with atom mapping.

///|
/// Build a reference chain: 5 ALA residues with CA atoms.
fn make_ref_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 5; r = r + 1 {
    let x = r.to_double() * 3.8
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, 0.0, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
/// Build a mobile chain with slight perturbation.
fn make_perturbed_chain(dx : Double, dy : Double) -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 5; r = r + 1 {
    let x = r.to_double() * 3.8 + dx
    let y = dy
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, y, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
/// Build a chain with missing residue 3.
fn make_missing_res_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 5; r = r + 1 {
    if r == 2 {
      continue
    }
    let x = r.to_double() * 3.8 + 0.5
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, 0.1, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
/// Build a chain with reordered residues (5,4,3,2,1).
fn make_reordered_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 4; r >= 0; r = r - 1 {
    let x = r.to_double() * 3.8 + 0.3
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, 0.2, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[ca])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
/// Build a chain with extra atoms (N, CA, C per residue).
fn make_full_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 5; r = r + 1 {
    let x = r.to_double() * 3.8
    let n = @src.Atom::new(name="N", coord=@src.Vector3::new(x - 0.5, 0.0, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="N")
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, 0.0, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    let c = @src.Atom::new(name="C", coord=@src.Vector3::new(x + 0.5, 0.0, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[n, ca, c])
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

// ===== AtomMap tests =====

test "atommap_basic" {
  let ref = make_ref_chain()
  let mob = make_perturbed_chain(0.5, 0.1)
  let map = @src.AtomMap::new(ref.get_atoms(), mob.get_atoms())
  assert_eq(@src.AtomMap::num_mapped(map), 5)
  assert_true(@src.AtomMap::coverage(map, 5) - 1.0 < 0.001)
}

test "atommap_missing" {
  let ref = make_ref_chain()
  let mob = make_missing_res_chain()
  let map = @src.AtomMap::new(ref.get_atoms(), mob.get_atoms())
  assert_eq(@src.AtomMap::num_mapped(map), 4)
  assert_true(@src.AtomMap::coverage(map, 5) - 0.8 < 0.001)
}

test "atommap_reordered" {
  let ref = make_ref_chain()
  let mob = make_reordered_chain()
  let map = @src.AtomMap::new(ref.get_atoms(), mob.get_atoms())
  assert_eq(@src.AtomMap::num_mapped(map), 5)
}

test "atommap_indices" {
  let ref = make_ref_chain()
  let mob = make_reordered_chain()
  let map = @src.AtomMap::new(ref.get_atoms(), mob.get_atoms())
  let ri = @src.AtomMap::ref_indices(map)
  let mi = @src.AtomMap::mob_indices(map)
  assert_eq(ri.length(), 5)
  assert_eq(mi.length(), 5)
  for i = 0; i < 5; i = i + 1 {
    assert_eq(ri[i], i)
  }
}

// ===== PDBEnsemble construction =====

test "pdb_ensemble_construction" {
  let ref = make_ref_chain()
  match @src.PDBEnsemble::new(ref, true) {
    Some(ens) => {
      assert_eq(@src.PDBEnsemble::num_atoms(ens), 5)
      assert_eq(@src.PDBEnsemble::num_conformations(ens), 0)
      assert_true(!@src.PDBEnsemble::is_aligned(ens))
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_add_chain" {
  let ref = make_ref_chain()
  match @src.PDBEnsemble::new(ref, true) {
    Some(ens) => {
      let mob = make_perturbed_chain(0.5, 0.1)
      let ens2 = ens.add_chain(mob)
      assert_eq(@src.PDBEnsemble::num_conformations(ens2), 1)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_from_chains" {
  let ref = make_ref_chain()
  let mob1: @src.Chain = make_perturbed_chain(0.5, 0.1)
  let mob2: @src.Chain = make_perturbed_chain(-0.3, 0.2)
  match @src.PDBEnsemble::from_chains([ref, mob1, mob2], true) {
    Some(ens) => {
      assert_eq(@src.PDBEnsemble::num_atoms(ens), 5)
      assert_eq(@src.PDBEnsemble::num_conformations(ens), 2)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_add_coordset" {
  let ref = make_ref_chain()
  match @src.PDBEnsemble::new(ref, true) {
    Some(ens) => {
      let coords = [
        @src.Vector3::new(0.5, 0.1, 0.0),
        @src.Vector3::new(4.3, 0.1, 0.0),
        @src.Vector3::new(8.1, 0.1, 0.0),
        @src.Vector3::new(11.9, 0.1, 0.0),
        @src.Vector3::new(15.7, 0.1, 0.0),
      ]
      let ens2 = ens.add_coordset(coords)
      assert_eq(@src.PDBEnsemble::num_conformations(ens2), 1)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_add_mapped" {
  let ref = make_ref_chain()
  let mob = make_reordered_chain()
  let map = @src.AtomMap::new(ref.get_atoms(), mob.get_atoms())
  match @src.PDBEnsemble::new(ref, true) {
    Some(ens) => {
      let ens2 = ens.add_mapped(mob.get_atoms(), map)
      assert_eq(@src.PDBEnsemble::num_conformations(ens2), 1)
    }
    None => assert_true(false)
  }
}

// ===== Superimposition =====

test "pdb_ensemble_superimpose" {
  let ref = make_ref_chain()
  let mob = make_perturbed_chain(0.5, 0.1)
  match @src.PDBEnsemble::from_chains([ref, mob], true) {
    Some(ens) => {
      let aligned = ens.superimpose()
      assert_true(@src.PDBEnsemble::is_aligned(aligned))
      let rmsd = @src.PDBEnsemble::get_rmsd(aligned)
      assert_eq(rmsd.length(), 1)
      assert_true(rmsd[0] < 1.0)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_superimpose_identity" {
  let ref = make_ref_chain()
  match @src.PDBEnsemble::from_chains([ref, ref], true) {
    Some(ens) => {
      let aligned = ens.superimpose()
      let rmsd = @src.PDBEnsemble::get_rmsd(aligned)
      assert_true(rmsd[0] < 0.001)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_iterative_superimpose" {
  let ref = make_ref_chain()
  let mob1: @src.Chain = make_perturbed_chain(0.5, 0.1)
  let mob2: @src.Chain = make_perturbed_chain(-0.3, 0.2)
  match @src.PDBEnsemble::from_chains([ref, mob1, mob2], true) {
    Some(ens) => {
      let aligned = ens.iterative_superimpose(3)
      assert_true(@src.PDBEnsemble::is_aligned(aligned))
      let rmsd = @src.PDBEnsemble::get_rmsd(aligned)
      assert_eq(rmsd.length(), 2)
    }
    None => assert_true(false)
  }
}

// ===== RMSD =====

test "pdb_ensemble_rmsd" {
  let ref = make_ref_chain()
  let mob = make_perturbed_chain(0.5, 0.1)
  match @src.PDBEnsemble::from_chains([ref, mob], true) {
    Some(ens) => {
      let rmsd = @src.PDBEnsemble::get_rmsd(ens)
      assert_true(rmsd[0] > 0.01)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_rmsd_matrix" {
  let ref = make_ref_chain()
  let mob1: @src.Chain = make_perturbed_chain(0.5, 0.1)
  let mob2: @src.Chain = make_perturbed_chain(-0.3, 0.2)
  match @src.PDBEnsemble::from_chains([ref, mob1, mob2], true) {
    Some(ens) => {
      let aligned = ens.superimpose()
      let matrix = @src.PDBEnsemble::get_rmsd_matrix(aligned)
      assert_eq(matrix.length(), 2)
      assert_eq(matrix[0].length(), 2)
      assert_true(matrix[0][0] < 0.001)
      assert_true(matrix[1][1] < 0.001)
      assert_true((matrix[0][1] - matrix[1][0]).abs() < 0.001)
    }
    None => assert_true(false)
  }
}

// ===== RMSF =====

test "pdb_ensemble_rmsf" {
  let ref = make_ref_chain()
  let mob1: @src.Chain = make_perturbed_chain(0.5, 0.1)
  let mob2: @src.Chain = make_perturbed_chain(-0.3, 0.2)
  match @src.PDBEnsemble::from_chains([ref, mob1, mob2], true) {
    Some(ens) => {
      let aligned = ens.superimpose()
      let rmsf = @src.PDBEnsemble::get_rmsf(aligned)
      assert_eq(rmsf.length(), 5)
      for i = 0; i < 5; i = i + 1 {
        assert_true(rmsf[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_rmsf_identical" {
  let ref = make_ref_chain()
  match @src.PDBEnsemble::from_chains([ref, ref, ref], true) {
    Some(ens) => {
      let aligned = ens.superimpose()
      let rmsf = @src.PDBEnsemble::get_rmsf(aligned)
      for i = 0; i < 5; i = i + 1 {
        assert_true(rmsf[i] < 0.001)
      }
    }
    None => assert_true(false)
  }
}

// ===== Mean structure =====

test "pdb_ensemble_mean_structure" {
  let ref = make_ref_chain()
  match @src.PDBEnsemble::from_chains([ref, ref], true) {
    Some(ens) => {
      let aligned = ens.superimpose()
      let mean = @src.PDBEnsemble::get_mean_structure(aligned)
      assert_eq(mean.length(), 5)
      let ref_coords = @src.PDBEnsemble::get_ref_coords(ens)
      for i = 0; i < 5; i = i + 1 {
        assert_true((mean[i].x - ref_coords[i].x).abs() < 0.001)
        assert_true((mean[i].y - ref_coords[i].y).abs() < 0.001)
      }
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_coord_deviations" {
  let ref = make_ref_chain()
  let mob = make_perturbed_chain(0.0, 0.5)
  match @src.PDBEnsemble::from_chains([ref, mob], true) {
    Some(ens) => {
      let aligned = ens.superimpose()
      let dev = @src.PDBEnsemble::get_coord_deviations(aligned)
      assert_eq(dev.length(), 5)
      for i = 0; i < 5; i = i + 1 {
        assert_true(dev[i] >= 0.0)
      }
    }
    None => assert_true(false)
  }
}

// ===== Missing atoms =====

test "pdb_ensemble_missing_atoms" {
  let ref = make_ref_chain()
  let mob = make_missing_res_chain()
  match @src.PDBEnsemble::from_chains([ref, mob], true) {
    Some(ens) => {
      let coverage = @src.PDBEnsemble::get_coverage(ens)
      assert_eq(coverage.length(), 1)
      assert_true(coverage[0] - 0.8 < 0.001)
      let aligned = ens.superimpose()
      let rmsd = @src.PDBEnsemble::get_rmsd(aligned)
      assert_eq(rmsd.length(), 1)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_reordered_mapping" {
  let ref = make_ref_chain()
  let mob = make_reordered_chain()
  match @src.PDBEnsemble::from_chains([ref, mob], true) {
    Some(ens) => {
      let coverage = @src.PDBEnsemble::get_coverage(ens)
      assert_true(coverage[0] - 1.0 < 0.001)
    }
    None => assert_true(false)
  }
}

// ===== Full atom ensemble =====

test "pdb_ensemble_full_atoms" {
  let ref = make_full_chain()
  match @src.PDBEnsemble::new(ref, false) {
    Some(ens) => {
      assert_eq(@src.PDBEnsemble::num_atoms(ens), 15)
      let ens2 = ens.add_chain(ref)
      assert_eq(@src.PDBEnsemble::num_conformations(ens2), 1)
    }
    None => assert_true(false)
  }
}

// ===== Utility =====

test "pdb_ensemble_remove_conformation" {
  let ref = make_ref_chain()
  let mob1: @src.Chain = make_perturbed_chain(0.5, 0.1)
  let mob2: @src.Chain = make_perturbed_chain(-0.3, 0.2)
  match @src.PDBEnsemble::from_chains([ref, mob1, mob2], true) {
    Some(ens) => {
      assert_eq(@src.PDBEnsemble::num_conformations(ens), 2)
      let ens2 = ens.remove_conformation(0)
      assert_eq(@src.PDBEnsemble::num_conformations(ens2), 1)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_weights" {
  let ref = make_ref_chain()
  match @src.PDBEnsemble::new(ref, true) {
    Some(ens) => {
      let ens2 = ens.set_weights([1.0, 2.0, 1.0, 2.0, 1.0])
      let wts = @src.PDBEnsemble::get_weights(ens2)
      assert_eq(wts.length(), 5)
      assert_true(wts[1] - 2.0 < 0.001)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_metadata" {
  let ref = make_full_chain()
  match @src.PDBEnsemble::new(ref, false) {
    Some(ens) => {
      let names = @src.PDBEnsemble::get_atom_names(ens)
      let chains = @src.PDBEnsemble::get_atom_chains(ens)
      let resseqs = @src.PDBEnsemble::get_atom_resseqs(ens)
      assert_eq(names.length(), 15)
      assert_eq(chains.length(), 15)
      assert_eq(resseqs.length(), 15)
      assert_eq(names[0], "N")
      assert_eq(names[1], "CA")
      assert_eq(names[2], "C")
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_get_coords" {
  let ref = make_ref_chain()
  let mob = make_perturbed_chain(0.5, 0.1)
  match @src.PDBEnsemble::from_chains([ref, mob], true) {
    Some(ens) => {
      let coords = @src.PDBEnsemble::get_coords(ens, 0)
      assert_eq(coords.length(), 5)
      let ref_coords = @src.PDBEnsemble::get_ref_coords(ens)
      assert_eq(ref_coords.length(), 5)
    }
    None => assert_true(false)
  }
}

test "pdb_ensemble_too_few_atoms" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA", chainid='A', resseq=1, element="C")
  let res = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[ca])
  residues.push(res)
  let chain = @src.Chain::new(id='A', residues~)
  match @src.PDBEnsemble::new(chain, true) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}
