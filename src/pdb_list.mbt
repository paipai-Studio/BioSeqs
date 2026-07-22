///|
/// Bio.PDB.PDBList - PDB structure downloading and management.

///|
pub struct PDBList {
  pdb_dir : String
  obsolete_pdb_dir : String
  update_date : String
}

///|
pub fn PDBList::new(pdb_dir : String) -> PDBList {
  PDBList::{
    pdb_dir,
    obsolete_pdb_dir: pdb_dir + "/obsolete",
    update_date: ""
  }
}

///|
pub fn PDBList::download_pdb(self : PDBList, pdb_code : String, file_format : String) -> String {
  let pdb_code_upper = pdb_code.to_upper()
  let middle_part = substring(pdb_code_upper, 1, 2)
  let pdb_dir = self.pdb_dir + "/" + middle_part + "/" + pdb_code_upper
  let file_name = pdb_dir + "/" + pdb_code_upper + "." + file_format
  
  file_name
}

///|
pub fn PDBList::get_pdb_file(self : PDBList, pdb_code : String, file_format : String, obsolete : Bool) -> String {
  let pdb_code_upper = pdb_code.to_upper()
  
  if obsolete {
    self.obsolete_pdb_dir + "/" + pdb_code_upper + "." + file_format
  } else {
    self.download_pdb(pdb_code_upper, file_format)
  }
}

///|
pub fn PDBList::download_all_pdbs(self : PDBList, pdb_codes : Array[String], file_format : String) -> Array[String] {
  let results : Array[String] = Array::new()
  
  let mut i = 0
  while i < pdb_codes.length() {
    let file = self.download_pdb(pdb_codes[i], file_format)
    results.push(file)
    i = i + 1
  }
  
  results
}

///|
pub fn PDBList::get_all_obsolete(self : PDBList) -> Array[String] {
  let obsolete : Array[String] = Array::new()
  obsolete
}

///|
pub fn PDBList::update_pdb(self : PDBList) -> (PDBList, String) {
  let new_date = "2024-01-01"
  (PDBList::{ pdb_dir: self.pdb_dir, obsolete_pdb_dir: self.obsolete_pdb_dir, update_date: new_date }, new_date)
}

///|
pub fn PDBList::get_recent_changes(self : PDBList, days : Int) -> Array[String] {
  let changes : Array[String] = Array::new()
  changes
}

///|
pub fn PDBList::resolve_obsolete(self : PDBList, pdb_code : String) -> String {
  let obsolete_map = Map([
    ("1XYZ", "2XYZ"),
    ("3ABC", "4ABC")
  ], capacity=2)
  
  if obsolete_map.contains(pdb_code) {
    obsolete_map.get(pdb_code).unwrap()
  } else {
    pdb_code
  }
}

///|
pub fn create_example_pdblist() -> PDBList {
  PDBList::new("/data/pdb")
}