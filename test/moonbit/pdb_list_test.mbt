///|
/// Tests for PDBList module.

test "PDBList creation" {
  let pdblist = @src.create_example_pdblist()
  assert_eq(pdblist.pdb_dir, "/data/pdb")
}

test "PDBList download_pdb" {
  let pdblist = @src.PDBList::new("/data/pdb")
  let file = pdblist.download_pdb("1XYZ", "pdb")
  assert_true(file.contains("1XYZ"))
}

test "PDBList get_pdb_file" {
  let pdblist = @src.PDBList::new("/data/pdb")
  let file = pdblist.get_pdb_file("1XYZ", "pdb", false)
  assert_true(file.contains("1XYZ"))
}

test "PDBList get_pdb_file obsolete" {
  let pdblist = @src.PDBList::new("/data/pdb")
  let file = pdblist.get_pdb_file("1XYZ", "pdb", true)
  assert_true(file.contains("obsolete"))
}

test "PDBList resolve_obsolete" {
  let pdblist = @src.PDBList::new("/data/pdb")
  let resolved = pdblist.resolve_obsolete("1XYZ")
  assert_eq(resolved, "2XYZ")
}