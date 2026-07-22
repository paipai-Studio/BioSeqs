///|
/// PDBList Demo - Demonstrates PDB structure downloading and management.

fn main {
  let pdblist = @bio.PDBList::new("/data/pdb")
  
  println("=== PDBList Demo ===")
  println("PDB Directory: \{pdblist.pdb_dir}")
  
  let file = pdblist.download_pdb("1XYZ", "pdb")
  println("Download path: \{file}")
  
  let obsolete_file = pdblist.get_pdb_file("1XYZ", "pdb", true)
  println("Obsolete path: \{obsolete_file}")
  
  let resolved = pdblist.resolve_obsolete("1XYZ")
  println("Resolved obsolete: \{resolved}")
}