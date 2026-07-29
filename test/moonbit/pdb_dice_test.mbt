///|
/// Tests for pdb_dice module.

// ============================================================================
// Helper: create a test structure
// ============================================================================

fn create_test_structure() -> @src.Structure {
  // Create a simple structure with 2 chains, each with 3 residues
  let atoms_a1 : Array[@src.Atom] = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0), resname="ALA", chainid='A', resseq=1),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 0.0, 0.0), resname="ALA", chainid='A', resseq=1),
    @src.Atom::new(name="C", coord=@src.Vector3::new(2.5, 1.0, 0.0), resname="ALA", chainid='A', resseq=1),
    @src.Atom::new(name="O", coord=@src.Vector3::new(2.5, 2.0, 0.0), resname="ALA", chainid='A', resseq=1),
  ]
  let atoms_a2 : Array[@src.Atom] = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(3.5, 0.0, 0.0), resname="GLY", chainid='A', resseq=2),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(5.0, 0.0, 0.0), resname="GLY", chainid='A', resseq=2),
  ]
  let atoms_a3 : Array[@src.Atom] = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(6.0, 0.0, 0.0), resname="VAL", chainid='A', resseq=3),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(7.5, 0.0, 0.0), resname="VAL", chainid='A', resseq=3),
  ]
  let res_a1 = @src.Residue::new(resname="ALA", chainid='A', resseq=1, atoms=atoms_a1)
  let res_a2 = @src.Residue::new(resname="GLY", chainid='A', resseq=2, atoms=atoms_a2)
  let res_a3 = @src.Residue::new(resname="VAL", chainid='A', resseq=3, atoms=atoms_a3)
  let chain_a = @src.Chain::new(id='A', residues=[res_a1, res_a2, res_a3])

  let atoms_b1 : Array[@src.Atom] = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 5.0, 0.0), resname="SER", chainid='B', resseq=10),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 5.0, 0.0), resname="SER", chainid='B', resseq=10),
  ]
  let atoms_b2 : Array[@src.Atom] = [
    @src.Atom::new(name="N", coord=@src.Vector3::new(3.5, 5.0, 0.0), resname="THR", chainid='B', resseq=11),
    @src.Atom::new(name="CA", coord=@src.Vector3::new(5.0, 5.0, 0.0), resname="THR", chainid='B', resseq=11),
  ]
  let res_b1 = @src.Residue::new(resname="SER", chainid='B', resseq=10, atoms=atoms_b1)
  let res_b2 = @src.Residue::new(resname="THR", chainid='B', resseq=11, atoms=atoms_b2)
  let chain_b = @src.Chain::new(id='B', residues=[res_b1, res_b2])

  let model = @src.Model::new(id=0, chains=[chain_a, chain_b])
  @src.Structure::new(id="TEST", models=[model])
}

// ============================================================================
// ResidueId tests
// ============================================================================

///|
test "residue_id_creation" {
  let rid = @src.ResidueId::new('A', 10)
  assert_eq(rid.chain_id, 'A')
  assert_eq(rid.resseq, 10)
  assert_eq(rid.icode, ' ')
}

///|
test "residue_id_with_icode" {
  let rid = @src.ResidueId::new('A', 10, icode='A')
  assert_eq(rid.chain_id, 'A')
  assert_eq(rid.resseq, 10)
  assert_eq(rid.icode, 'A')
}

// ============================================================================
// DiceResult tests
// ============================================================================

///|
test "dice_result_creation" {
  let result = @src.DiceResult::new("test", 100, 20, 2, 1)
  assert_eq(result.structure_id, "test")
  assert_eq(result.num_atoms, 100)
  assert_eq(result.num_residues, 20)
  assert_eq(result.num_chains, 2)
  assert_eq(result.num_models, 1)
}

// ============================================================================
// Chain extraction tests
// ============================================================================

///|
test "extract_single_chain" {
  let structure = create_test_structure()
  let result = @src.extract_chain(structure, 'A')
  assert_eq(result.models.length(), 1)
  assert_eq(result.models[0].chains.length(), 1)
  assert_eq(result.models[0].chains[0].id, 'A')
}

///|
test "extract_multiple_chains" {
  let structure = create_test_structure()
  let chains = ['A', 'B']
  let result = @src.extract_chains(structure, chains)
  assert_eq(result.models.length(), 1)
  assert_eq(result.models[0].chains.length(), 2)
}

///|
test "extract_nonexistent_chain" {
  let structure = create_test_structure()
  let result = @src.extract_chain(structure, 'Z')
  assert_eq(result.models.length(), 0)
}

// ============================================================================
// Residue range extraction tests
// ============================================================================

///|
test "extract_residue_range" {
  let structure = create_test_structure()
  let result = @src.extract_residue_range(structure, 'A', 1, 2)
  assert_eq(result.models.length(), 1)
  assert_eq(result.models[0].chains.length(), 1)
  assert_eq(result.models[0].chains[0].residues.length(), 2)
  assert_eq(result.models[0].chains[0].residues[0].resseq, 1)
  assert_eq(result.models[0].chains[0].residues[1].resseq, 2)
}

///|
test "extract_residue_range_single" {
  let structure = create_test_structure()
  let result = @src.extract_residue_range(structure, 'A', 2, 2)
  assert_eq(result.models[0].chains[0].residues.length(), 1)
  assert_eq(result.models[0].chains[0].residues[0].resseq, 2)
}

///|
test "extract_residue_range_chain_b" {
  let structure = create_test_structure()
  let result = @src.extract_residue_range(structure, 'B', 10, 11)
  assert_eq(result.models[0].chains.length(), 1)
  assert_eq(result.models[0].chains[0].id, 'B')
  assert_eq(result.models[0].chains[0].residues.length(), 2)
}

// ============================================================================
// Residue extraction by IDs
// ============================================================================

///|
test "extract_residues_by_ids" {
  let structure = create_test_structure()
  let ids = [
    @src.ResidueId::new('A', 1),
    @src.ResidueId::new('B', 10),
  ]
  let result = @src.extract_residues(structure, ids)
  // Should have 2 chains, each with 1 residue
  assert_eq(result.models[0].chains.length(), 2)
}

///|
test "extract_residues_empty_ids" {
  let structure = create_test_structure()
  let empty : Array[@src.ResidueId] = []
  let result = @src.extract_residues(structure, empty)
  assert_eq(result.models.length(), 0)
}

// ============================================================================
// Atom selection tests
// ============================================================================

///|
test "extract_ca_atoms" {
  let structure = create_test_structure()
  let result = @src.extract_ca_atoms(structure)
  // Each residue should have 1 CA atom
  let total_residues = result.get_num_residues()
  let total_atoms = result.get_num_atoms()
  assert_eq(total_atoms, total_residues)
}

///|
test "extract_backbone_atoms" {
  let structure = create_test_structure()
  let result = @src.extract_backbone_atoms(structure)
  assert_true(result.get_num_atoms() > 0)
}

///|
test "extract_atoms_by_name" {
  let structure = create_test_structure()
  let names = ["N", "CA"]
  let result = @src.extract_atoms_by_name(structure, names)
  assert_true(result.get_num_atoms() > 0)
  assert_true(result.get_num_atoms() < structure.get_num_atoms())
}

// ============================================================================
// Model extraction tests
// ============================================================================

///|
test "extract_first_model" {
  let structure = create_test_structure()
  let result = @src.extract_first_model(structure)
  assert_eq(result.models.length(), 1)
}

///|
test "extract_models_by_id" {
  let structure = create_test_structure()
  let ids = [0]
  let result = @src.extract_models(structure, ids)
  assert_eq(result.models.length(), 1)
}

///|
test "extract_models_nonexistent" {
  let structure = create_test_structure()
  let ids = [99]
  let result = @src.extract_models(structure, ids)
  assert_eq(result.models.length(), 0)
}

// ============================================================================
// Het residue tests
// ============================================================================

///|
test "remove_hetero_residues" {
  let structure = create_test_structure()
  let result = @src.remove_hetero_residues(structure)
  // All residues are standard, so should be same count
  assert_eq(result.get_num_residues(), structure.get_num_residues())
}

///|
test "keep_hetero_residues_none" {
  let structure = create_test_structure()
  let result = @src.keep_hetero_residues(structure)
  // No hetero residues in test structure
  assert_eq(result.models.length(), 0)
}

// ============================================================================
// Selection utility tests
// ============================================================================

///|
test "get_unique_chain_ids" {
  let structure = create_test_structure()
  let ids = @src.get_unique_chain_ids(structure)
  assert_eq(ids.length(), 2)
}

///|
test "get_unique_residue_names" {
  let structure = create_test_structure()
  let names = @src.get_unique_residue_names(structure)
  assert_true(names.length() >= 4) // ALA, GLY, VAL, SER, THR
}

///|
test "get_unique_atom_names" {
  let structure = create_test_structure()
  let names = @src.get_unique_atom_names(structure)
  assert_true(names.length() >= 4) // N, CA, C, O
}

///|
test "count_atoms_by_element" {
  let structure = create_test_structure()
  let counts = @src.count_atoms_by_element(structure)
  // Most atoms have no element set, so count "?"
  let unknown = counts.get("?").unwrap_or(0)
  assert_eq(unknown, structure.get_num_atoms())
}

// ============================================================================
// Structure statistics tests
// ============================================================================

///|
test "structure_stats" {
  let structure = create_test_structure()
  let stats = @src.structure_stats(structure)
  assert_eq(stats.structure_id, "TEST")
  assert_eq(stats.num_models, 1)
  assert_eq(stats.num_chains, 2)
  assert_eq(stats.num_residues, 5) // 3 in A, 2 in B
  assert_true(stats.num_atoms > 0)
}

// ============================================================================
// Chain sequence tests
// ============================================================================

///|
test "chain_sequence_chain_a" {
  let structure = create_test_structure()
  let seq = @src.chain_sequence(structure, 'A')
  assert_eq(seq, "AGV")
}

///|
test "chain_sequence_chain_b" {
  let structure = create_test_structure()
  let seq = @src.chain_sequence(structure, 'B')
  assert_eq(seq, "ST")
}

///|
test "chain_sequence_nonexistent" {
  let structure = create_test_structure()
  let seq = @src.chain_sequence(structure, 'Z')
  assert_eq(seq, "")
}

// ============================================================================
// B-factor and occupancy filtering tests
// ============================================================================

///|
test "filter_by_bfactor" {
  let structure = create_test_structure()
  // All atoms have bfactor 0.0, so filter should keep all
  let result = @src.filter_by_bfactor(structure, -1.0, 100.0)
  assert_eq(result.get_num_atoms(), structure.get_num_atoms())
}

///|
test "filter_by_bfactor_none" {
  let structure = create_test_structure()
  // All atoms have bfactor 0.0, so filter with min 10 should remove all
  let result = @src.filter_by_bfactor(structure, 10.0, 100.0)
  assert_eq(result.models.length(), 0)
}

///|
test "filter_by_occupancy" {
  let structure = create_test_structure()
  // All atoms have occupancy 1.0
  let result = @src.filter_by_occupancy(structure, 0.5)
  assert_eq(result.get_num_atoms(), structure.get_num_atoms())
}

///|
test "filter_by_occupancy_high" {
  let structure = create_test_structure()
  let result = @src.filter_by_occupancy(structure, 2.0)
  assert_eq(result.models.length(), 0)
}

// ============================================================================
// Geometric selection tests
// ============================================================================

///|
test "select_atoms_near_point" {
  let structure = create_test_structure()
  let point = @src.Vector3::new(0.0, 0.0, 0.0)
  let atoms = @src.select_atoms_near_point(structure, point, 2.0)
  assert_true(atoms.length() > 0)
}

///|
test "select_atoms_near_point_far" {
  let structure = create_test_structure()
  let point = @src.Vector3::new(100.0, 100.0, 100.0)
  let atoms = @src.select_atoms_near_point(structure, point, 1.0)
  assert_eq(atoms.length(), 0)
}

///|
test "select_residues_near_point" {
  let structure = create_test_structure()
  let point = @src.Vector3::new(0.0, 0.0, 0.0)
  let residues = @src.select_residues_near_point(structure, point, 3.0)
  assert_true(residues.length() > 0)
}

///|
test "select_residues_near_point_far" {
  let structure = create_test_structure()
  let point = @src.Vector3::new(100.0, 100.0, 100.0)
  let residues = @src.select_residues_near_point(structure, point, 1.0)
  assert_eq(residues.length(), 0)
}
