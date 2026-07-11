///|
/// Bio.PDB - Protein Data Bank data structures.

///|
/// A point in 3D space with x, y, z coordinates.
pub struct Vector3 {
  x : Double
  y : Double
  z : Double
}

///|
/// Construct a Vector3.
pub fn Vector3::new(x : Double, y : Double, z : Double) -> Vector3 {
  Vector3::{ x, y, z }
}

///|
/// Calculate the Euclidean distance between two vectors.
pub fn Vector3::distance(self : Vector3, other : Vector3) -> Double {
  let dx = self.x - other.x
  let dy = self.y - other.y
  let dz = self.z - other.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// An atom in a protein structure.
pub struct Atom {
  name : String
  altloc : Char
  resname : String
  chainid : Char
  resseq : Int
  icode : Char
  coord : Vector3
  occupancy : Double
  bfactor : Double
  element : String
  charge : Int
}

///|
/// Construct an Atom with default values.
pub fn Atom::new(
  name~ : String,
  coord~ : Vector3,
  resname? : String = "UNK",
  chainid? : Char = ' ',
  resseq? : Int = 1,
  occupancy? : Double = 1.0,
  bfactor? : Double = 0.0,
  element? : String = "",
  charge? : Int = 0,
) -> Atom {
  Atom::{
    name,
    altloc: ' ',
    resname,
    chainid,
    resseq,
    icode: ' ',
    coord,
    occupancy,
    bfactor,
    element,
    charge,
  }
}

///|
/// Get the coordinates of the atom.
pub fn Atom::get_coord(self : Atom) -> Vector3 {
  self.coord
}

///|
/// Calculate the distance to another atom.
pub fn Atom::distance(self : Atom, other : Atom) -> Double {
  self.coord.distance(other.coord)
}

///|
/// A residue in a protein chain.
pub struct Residue {
  resname : String
  chainid : Char
  resseq : Int
  icode : Char
  atoms : Array[Atom]
  hetfield : String
}

///|
/// Construct a Residue.
pub fn Residue::new(
  resname~ : String,
  chainid~ : Char,
  resseq~ : Int,
  atoms? : Array[Atom] = [],
  hetfield? : String = " ",
) -> Residue {
  Residue::{ resname, chainid, resseq, icode: ' ', atoms, hetfield }
}

///|
/// Check if this residue is a hetero-residue.
pub fn Residue::is_het(self : Residue) -> Bool {
  self.hetfield != " "
}

///|
/// Get the number of atoms in this residue.
pub fn Residue::get_num_atoms(self : Residue) -> Int {
  self.atoms.length()
}

///|
/// Get all atoms in this residue.
pub fn Residue::get_atoms(self : Residue) -> Array[Atom] {
  self.atoms
}

///|
/// Find an atom by name.
pub fn Residue::find_atom(self : Residue, name : String) -> Atom? {
  for atom in self.atoms {
    if atom.name == name {
      return Some(atom)
    }
  }
  None
}

///|
/// Get the CA (alpha carbon) atom if it exists.
pub fn Residue::get_ca(self : Residue) -> Atom? {
  self.find_atom("CA")
}

///|
/// Get the center of mass of this residue.
pub fn Residue::center_of_mass(self : Residue) -> Vector3? {
  if self.atoms.length() == 0 {
    return None
  }

  let mut x = 0.0
  let mut y = 0.0
  let mut z = 0.0
  for atom in self.atoms {
    x = x + atom.coord.x
    y = y + atom.coord.y
    z = z + atom.coord.z
  }

  let n = self.atoms.length().to_double()
  Some(Vector3::new(x / n, y / n, z / n))
}

///|
/// A chain in a protein model.
pub struct Chain {
  id : Char
  residues : Array[Residue]
}

///|
/// Construct a Chain.
pub fn Chain::new(id~ : Char, residues? : Array[Residue] = []) -> Chain {
  Chain::{ id, residues }
}

///|
/// Get the chain ID.
pub fn Chain::get_id(self : Chain) -> Char {
  self.id
}

///|
/// Get the number of residues in this chain.
pub fn Chain::get_num_residues(self : Chain) -> Int {
  self.residues.length()
}

///|
/// Get all residues in this chain.
pub fn Chain::get_residues(self : Chain) -> Array[Residue] {
  self.residues
}

///|
/// Find a residue by sequence number.
pub fn Chain::find_residue(self : Chain, resseq : Int) -> Residue? {
  for residue in self.residues {
    if residue.resseq == resseq {
      return Some(residue)
    }
  }
  None
}

///|
/// Get all atoms in this chain.
pub fn Chain::get_atoms(self : Chain) -> Array[Atom] {
  let result = Array::new()
  for residue in self.residues {
    for atom in residue.atoms {
      result.push(atom)
    }
  }
  result
}

///|
/// A model in a protein structure.
pub struct Model {
  id : Int
  chains : Array[Chain]
}

///|
/// Construct a Model.
pub fn Model::new(id~ : Int, chains? : Array[Chain] = []) -> Model {
  Model::{ id, chains }
}

///|
/// Get the model ID.
pub fn Model::get_id(self : Model) -> Int {
  self.id
}

///|
/// Get the number of chains in this model.
pub fn Model::get_num_chains(self : Model) -> Int {
  self.chains.length()
}

///|
/// Get all chains in this model.
pub fn Model::get_chains(self : Model) -> Array[Chain] {
  self.chains
}

///|
/// Find a chain by ID.
pub fn Model::find_chain(self : Model, id : Char) -> Chain? {
  for chain in self.chains {
    if chain.id == id {
      return Some(chain)
    }
  }
  None
}

///|
/// Get all residues in this model.
pub fn Model::get_residues(self : Model) -> Array[Residue] {
  let result = Array::new()
  for chain in self.chains {
    for residue in chain.residues {
      result.push(residue)
    }
  }
  result
}

///|
/// Get all atoms in this model.
pub fn Model::get_atoms(self : Model) -> Array[Atom] {
  let result = Array::new()
  for chain in self.chains {
    for residue in chain.residues {
      for atom in residue.atoms {
        result.push(atom)
      }
    }
  }
  result
}

///|
/// A protein structure.
pub struct Structure {
  id : String
  models : Array[Model]
}

///|
/// Construct a Structure.
pub fn Structure::new(id~ : String, models? : Array[Model] = []) -> Structure {
  Structure::{ id, models }
}

///|
/// Get the structure ID.
pub fn Structure::get_id(self : Structure) -> String {
  self.id
}

///|
/// Get the number of models in this structure.
pub fn Structure::get_num_models(self : Structure) -> Int {
  self.models.length()
}

///|
/// Get the number of chains in this structure.
pub fn Structure::get_num_chains(self : Structure) -> Int {
  self.get_chains().length()
}

///|
/// Get the number of atoms in this structure.
pub fn Structure::get_num_atoms(self : Structure) -> Int {
  self.get_atoms().length()
}

///|
/// Get the number of residues in this structure.
pub fn Structure::get_num_residues(self : Structure) -> Int {
  self.get_residues().length()
}

///|
/// Get all models in this structure.
pub fn Structure::get_models(self : Structure) -> Array[Model] {
  self.models
}

///|
/// Find a model by ID.
pub fn Structure::find_model(self : Structure, id : Int) -> Model? {
  for model in self.models {
    if model.id == id {
      return Some(model)
    }
  }
  None
}

///|
/// Get all chains in this structure.
pub fn Structure::get_chains(self : Structure) -> Array[Chain] {
  let result = Array::new()
  for model in self.models {
    for chain in model.chains {
      result.push(chain)
    }
  }
  result
}

///|
/// Get all residues in this structure.
pub fn Structure::get_residues(self : Structure) -> Array[Residue] {
  let result = Array::new()
  for model in self.models {
    for chain in model.chains {
      for residue in chain.residues {
        result.push(residue)
      }
    }
  }
  result
}

///|
/// Get all atoms in this structure.
pub fn Structure::get_atoms(self : Structure) -> Array[Atom] {
  let result = Array::new()
  for model in self.models {
    for chain in model.chains {
      for residue in chain.residues {
        for atom in residue.atoms {
          result.push(atom)
        }
      }
    }
  }
  result
}

///|
/// Calculate the distance between two atoms by their identifiers.
pub fn Structure::distance(
  self : Structure,
  atom1_name : String,
  atom2_name : String,
) -> Double? {
  for model in self.models {
    for chain in model.chains {
      for residue in chain.residues {
        for atom in residue.atoms {
          if atom.name == atom1_name {
            for model2 in self.models {
              for chain2 in model2.chains {
                for residue2 in chain2.residues {
                  for atom2 in residue2.atoms {
                    if atom2.name == atom2_name {
                      return Some(atom.distance(atom2))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  None
}
