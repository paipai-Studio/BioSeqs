///|
test "structure_alignment_point3d" {
  let p = @src.SAPoint3D::new(1.0, 2.0, 3.0)

  assert_eq(p.x, 1.0)
  assert_eq(p.y, 2.0)
  assert_eq(p.z, 3.0)
}

///|
test "structure_alignment_create_residue" {
  let residue = @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(1.0, 0.0, 0.0))

  assert_eq(residue.resname, "ALA")
  assert_eq(residue.resseq, 1)
}

///|
test "structure_alignment_create_structure" {
  let residues = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let structure = @src.SAStructure::new("PDB1", residues)

  assert_eq(structure.id, "PDB1")
  assert_eq(structure.residues.length(), 2)
}

///|
test "structure_alignment_align_two" {
  let residues1 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
    @src.SAResidue::new("SER", 3, @src.SAPoint3D::new(2.0, 0.0, 0.0)),
  ]
  let residues2 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let struct1 = @src.SAStructure::new("PDB1", residues1)
  let struct2 = @src.SAStructure::new("PDB2", residues2)

  let result = @src.bio_structure_align_two(struct1, struct2)

  assert_eq(result.structure_ids.length(), 2)
  assert_true(result.aligned_length > 0)
}

///|
test "structure_alignment_align_multiple" {
  let residues1 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let residues2 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let residues3 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
  ]
  let structures = [
    @src.SAStructure::new("PDB1", residues1),
    @src.SAStructure::new("PDB2", residues2),
    @src.SAStructure::new("PDB3", residues3),
  ]

  let result = @src.bio_structure_align_multiple(structures)

  assert_eq(result.structure_ids.length(), 3)
}

///|
test "structure_alignment_calculate_rmsd" {
  let residues1 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let residues2 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let struct1 = @src.SAStructure::new("PDB1", residues1)
  let struct2 = @src.SAStructure::new("PDB2", residues2)

  let rmsd = @src.bio_structure_calculate_rmsd(struct1, struct2)

  assert_true(rmsd < 0.1)
}

///|
test "structure_alignment_calculate_tm_score" {
  let residues1 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let residues2 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let struct1 = @src.SAStructure::new("PDB1", residues1)
  let struct2 = @src.SAStructure::new("PDB2", residues2)

  let tm_score = @src.bio_structure_calculate_tm_score(struct1, struct2)

  assert_true(tm_score > 0.5)
}

///|
test "structure_alignment_pairwise_rmsd" {
  let residues1 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let residues2 = [
    @src.SAResidue::new("ALA", 1, @src.SAPoint3D::new(0.0, 0.0, 0.0)),
    @src.SAResidue::new("GLY", 2, @src.SAPoint3D::new(1.0, 0.0, 0.0)),
  ]
  let structures = [
    @src.SAStructure::new("PDB1", residues1),
    @src.SAStructure::new("PDB2", residues2),
  ]

  let matrix = @src.bio_structure_pairwise_rmsd(structures)

  assert_eq(matrix.length(), 2)
  assert_eq(matrix[0][0], 0.0)
}