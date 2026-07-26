///|
/// Bio.PDB.Dice - PDB structure extraction and selection tools.
/// Provides functions to extract chains, residues, and atoms from PDB structures,
/// similar to Biopython's Bio.PDB.Dice module.
/// Also includes Selection utilities for advanced structure selection.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// A residue identifier (chain, resseq, icode).
pub struct ResidueId {
  chain_id : Char
  resseq : Int
  icode : Char
} derive(Eq)

///|
/// Selection level enum for structure traversal.
pub enum SelectionLevel {
  AtomLevel
  ResidueLevel
  ChainLevel
  ModelLevel
  StructureLevel
}

///|
/// Dice selection result containing extracted structure info.
pub struct DiceResult {
  structure_id : String
  num_atoms : Int
  num_residues : Int
  num_chains : Int
  num_models : Int
} derive(Eq)

// ============================================================================
// Constructors
// ============================================================================

///|
pub fn ResidueId::new(
  chain_id : Char,
  resseq : Int,
  icode? : Char = ' ',
) -> ResidueId {
  ResidueId::{ chain_id, resseq, icode }
}

///|
pub fn DiceResult::new(
  structure_id : String,
  num_atoms : Int,
  num_residues : Int,
  num_chains : Int,
  num_models : Int,
) -> DiceResult {
  DiceResult::{ structure_id, num_atoms, num_residues, num_chains, num_models }
}

// ============================================================================
// Chain extraction
// ============================================================================

///|
/// Extract specific chains from a structure.
/// Returns a new Structure containing only the specified chains.
pub fn extract_chains(
  structure : Structure,
  chain_ids : Array[Char],
) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let mut selected = false
      for k = 0; k < chain_ids.length() && !selected; k = k + 1 {
        if chain.id == chain_ids[k] {
          selected = true
        }
      }
      if selected {
        new_chains.push(chain)
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

///|
/// Extract a single chain from a structure.
pub fn extract_chain(
  structure : Structure,
  chain_id : Char,
) -> Structure {
  let chain_arr = [chain_id]
  extract_chains(structure, chain_arr)
}

// ============================================================================
// Residue range extraction
// ============================================================================

///|
/// Extract a range of residues from a specific chain.
/// Returns a new Structure with only residues in [start, end] range.
pub fn extract_residue_range(
  structure : Structure,
  chain_id : Char,
  start_res : Int,
  end_res : Int,
) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      if chain.id == chain_id {
        let new_residues : Array[Residue] = []
        for k = 0; k < chain.residues.length(); k = k + 1 {
          let residue = chain.residues[k]
          if residue.resseq >= start_res && residue.resseq <= end_res {
            new_residues.push(residue)
          }
        }
        if new_residues.length() > 0 {
          let new_chain = Chain::new(id=chain.id, residues=new_residues)
          new_chains.push(new_chain)
        }
      } else {
        // Keep other chains as-is? No, dice extracts specific parts
        // Skip other chains
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

///|
/// Extract specific residues by their IDs from a structure.
pub fn extract_residues(
  structure : Structure,
  residue_ids : Array[ResidueId],
) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    let chain_residues_map = Map([], capacity=8) // chain_id -> Array[Residue]
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let selected_residues : Array[Residue] = []
      
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        let mut is_selected = false
        
        for r = 0; r < residue_ids.length() && !is_selected; r = r + 1 {
          let rid = residue_ids[r]
          if chain.id == rid.chain_id &&
             residue.resseq == rid.resseq &&
             residue.icode == rid.icode {
            is_selected = true
          }
        }
        
        if is_selected {
          selected_residues.push(residue)
        }
      }
      
      if selected_residues.length() > 0 {
        let new_chain = Chain::new(id=chain.id, residues=selected_residues)
        new_chains.push(new_chain)
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

// ============================================================================
// Atom selection
// ============================================================================

///|
/// Extract only atoms with specified names from a structure.
pub fn extract_atoms_by_name(
  structure : Structure,
  atom_names : Array[String],
) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let new_residues : Array[Residue] = []
      
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        let new_atoms : Array[Atom] = []
        
        for a = 0; a < residue.atoms.length(); a = a + 1 {
          let atom = residue.atoms[a]
          let mut selected = false
          for n = 0; n < atom_names.length() && !selected; n = n + 1 {
            if atom.name == atom_names[n] {
              selected = true
            }
          }
          if selected {
            new_atoms.push(atom)
          }
        }
        
        if new_atoms.length() > 0 {
          let new_residue = Residue::new(
            resname=residue.resname,
            chainid=residue.chainid,
            resseq=residue.resseq,
            atoms=new_atoms,
            hetfield=residue.hetfield,
          )
          new_residues.push(new_residue)
        }
      }
      
      if new_residues.length() > 0 {
        let new_chain = Chain::new(id=chain.id, residues=new_residues)
        new_chains.push(new_chain)
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

///|
/// Extract only CA (alpha carbon) atoms from a structure.
pub fn extract_ca_atoms(structure : Structure) -> Structure {
  let ca_names = ["CA"]
  extract_atoms_by_name(structure, ca_names)
}

///|
/// Extract only backbone atoms (N, CA, C, O) from a structure.
pub fn extract_backbone_atoms(structure : Structure) -> Structure {
  let backbone_names = ["N", "CA", "C", "O"]
  extract_atoms_by_name(structure, backbone_names)
}

// ============================================================================
// Model extraction
// ============================================================================

///|
/// Extract specific models from a structure.
pub fn extract_models(
  structure : Structure,
  model_ids : Array[Int],
) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let mut selected = false
    for j = 0; j < model_ids.length() && !selected; j = j + 1 {
      if model.id == model_ids[j] {
        selected = true
      }
    }
    if selected {
      new_models.push(model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

///|
/// Extract the first model from a structure.
pub fn extract_first_model(structure : Structure) -> Structure {
  if structure.models.length() == 0 {
    return structure
  }
  let first_id = structure.models[0].id
  extract_models(structure, [first_id])
}

// ============================================================================
// Het atom filtering
// ============================================================================

///|
/// Remove hetero residues (HETATM) from a structure.
/// Keeps only standard amino acid residues.
pub fn remove_hetero_residues(structure : Structure) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let new_residues : Array[Residue] = []
      
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        if !residue.is_het() {
          new_residues.push(residue)
        }
      }
      
      if new_residues.length() > 0 {
        let new_chain = Chain::new(id=chain.id, residues=new_residues)
        new_chains.push(new_chain)
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

///|
/// Keep only hetero residues (HETATM) in a structure.
pub fn keep_hetero_residues(structure : Structure) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let new_residues : Array[Residue] = []
      
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        if residue.is_het() {
          new_residues.push(residue)
        }
      }
      
      if new_residues.length() > 0 {
        let new_chain = Chain::new(id=chain.id, residues=new_residues)
        new_chains.push(new_chain)
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

// ============================================================================
// Selection utilities (Selection module)
// ============================================================================

///|
/// Get all unique chain IDs from a structure.
pub fn get_unique_chain_ids(structure : Structure) -> Array[Char] {
  let ids : Array[Char] = []
  let seen = Map([], capacity=8)
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let key = "" + chain.id.to_string()
      if seen.get(key) == None {
        seen.set(key, "1")
        ids.push(chain.id)
      }
    }
  }
  
  ids
}

///|
/// Get all unique residue names from a structure.
pub fn get_unique_residue_names(structure : Structure) -> Array[String] {
  let names : Array[String] = []
  let seen = Map([], capacity=32)
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        if seen.get(residue.resname) == None {
          seen.set(residue.resname, "1")
          names.push(residue.resname)
        }
      }
    }
  }
  
  names
}

///|
/// Get all unique atom names from a structure.
pub fn get_unique_atom_names(structure : Structure) -> Array[String] {
  let names : Array[String] = []
  let seen = Map([], capacity=32)
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        for a = 0; a < residue.atoms.length(); a = a + 1 {
          let atom = residue.atoms[a]
          if seen.get(atom.name) == None {
            seen.set(atom.name, "1")
            names.push(atom.name)
          }
        }
      }
    }
  }
  
  names
}

///|
/// Count atoms in a structure by element type.
pub fn count_atoms_by_element(structure : Structure) -> Map[String, Int] {
  let counts = Map([], capacity=16)
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        for a = 0; a < residue.atoms.length(); a = a + 1 {
          let atom = residue.atoms[a]
          let element = if atom.element != "" { atom.element } else { "?" }
          let current = counts.get(element).unwrap_or(0)
          counts.set(element, current + 1)
        }
      }
    }
  }
  
  counts
}

// ============================================================================
// Structure statistics
// ============================================================================

///|
/// Get detailed statistics about a structure.
pub fn structure_stats(structure : Structure) -> DiceResult {
  let mut num_atoms = 0
  let mut num_residues = 0
  let mut num_chains = 0
  let num_models = structure.models.length()
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    num_chains = num_chains + model.chains.length()
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      num_residues = num_residues + chain.residues.length()
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        num_atoms = num_atoms + residue.atoms.length()
      }
    }
  }
  
  DiceResult::new(
    structure.id,
    num_atoms,
    num_residues,
    num_chains,
    num_models,
  )
}

///|
/// Get the sequence of a specific chain as a string (one-letter codes).
pub fn chain_sequence(structure : Structure, chain_id : Char) -> String {
  let mut seq = ""
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      if chain.id == chain_id {
        for k = 0; k < chain.residues.length(); k = k + 1 {
          let residue = chain.residues[k]
          if !residue.is_het() {
            seq = seq + residue_name_to_letter(residue.resname)
          }
        }
        return seq // Return first model's sequence
      }
    }
  }
  
  seq
}

///|
fn residue_name_to_letter(resname : String) -> String {
  match resname {
    "ALA" => "A"
    "ARG" => "R"
    "ASN" => "N"
    "ASP" => "D"
    "CYS" => "C"
    "GLN" => "Q"
    "GLU" => "E"
    "GLY" => "G"
    "HIS" => "H"
    "ILE" => "I"
    "LEU" => "L"
    "LYS" => "K"
    "MET" => "M"
    "PHE" => "F"
    "PRO" => "P"
    "SER" => "S"
    "THR" => "T"
    "TRP" => "W"
    "TYR" => "Y"
    "VAL" => "V"
    _ => "X"
  }
}

// ============================================================================
// Structure filtering by B-factor
// ============================================================================

///|
/// Filter atoms by B-factor (temperature factor) range.
pub fn filter_by_bfactor(
  structure : Structure,
  min_bfactor : Double,
  max_bfactor : Double,
) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let new_residues : Array[Residue] = []
      
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        let new_atoms : Array[Atom] = []
        
        for a = 0; a < residue.atoms.length(); a = a + 1 {
          let atom = residue.atoms[a]
          if atom.bfactor >= min_bfactor && atom.bfactor <= max_bfactor {
            new_atoms.push(atom)
          }
        }
        
        if new_atoms.length() > 0 {
          let new_residue = Residue::new(
            resname=residue.resname,
            chainid=residue.chainid,
            resseq=residue.resseq,
            atoms=new_atoms,
            hetfield=residue.hetfield,
          )
          new_residues.push(new_residue)
        }
      }
      
      if new_residues.length() > 0 {
        let new_chain = Chain::new(id=chain.id, residues=new_residues)
        new_chains.push(new_chain)
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

///|
/// Filter atoms by occupancy value.
pub fn filter_by_occupancy(
  structure : Structure,
  min_occupancy : Double,
) -> Structure {
  let new_models : Array[Model] = []
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    let new_chains : Array[Chain] = []
    
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      let new_residues : Array[Residue] = []
      
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        let new_atoms : Array[Atom] = []
        
        for a = 0; a < residue.atoms.length(); a = a + 1 {
          let atom = residue.atoms[a]
          if atom.occupancy >= min_occupancy {
            new_atoms.push(atom)
          }
        }
        
        if new_atoms.length() > 0 {
          let new_residue = Residue::new(
            resname=residue.resname,
            chainid=residue.chainid,
            resseq=residue.resseq,
            atoms=new_atoms,
            hetfield=residue.hetfield,
          )
          new_residues.push(new_residue)
        }
      }
      
      if new_residues.length() > 0 {
        let new_chain = Chain::new(id=chain.id, residues=new_residues)
        new_chains.push(new_chain)
      }
    }
    
    if new_chains.length() > 0 {
      let new_model = Model::new(id=model.id, chains=new_chains)
      new_models.push(new_model)
    }
  }
  
  Structure::new(id=structure.id, models=new_models)
}

// ============================================================================
// Geometric selection
// ============================================================================

///|
/// Select atoms within a certain radius of a point in 3D space.
pub fn select_atoms_near_point(
  structure : Structure,
  point : Vector3,
  radius : Double,
) -> Array[Atom] {
  let selected : Array[Atom] = []
  let radius_sq = radius * radius
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        for a = 0; a < residue.atoms.length(); a = a + 1 {
          let atom = residue.atoms[a]
          let dx = atom.coord.x - point.x
          let dy = atom.coord.y - point.y
          let dz = atom.coord.z - point.z
          let dist_sq = dx * dx + dy * dy + dz * dz
          if dist_sq <= radius_sq {
            selected.push(atom)
          }
        }
      }
    }
  }
  
  selected
}

///|
/// Select residues that have any atom within a certain radius of a point.
pub fn select_residues_near_point(
  structure : Structure,
  point : Vector3,
  radius : Double,
) -> Array[ResidueId] {
  let selected : Array[ResidueId] = []
  let seen = Map([], capacity=32)
  let radius_sq = radius * radius
  
  for i = 0; i < structure.models.length(); i = i + 1 {
    let model = structure.models[i]
    for j = 0; j < model.chains.length(); j = j + 1 {
      let chain = model.chains[j]
      for k = 0; k < chain.residues.length(); k = k + 1 {
        let residue = chain.residues[k]
        let mut has_close_atom = false
        
        for a = 0; a < residue.atoms.length() && !has_close_atom; a = a + 1 {
          let atom = residue.atoms[a]
          let dx = atom.coord.x - point.x
          let dy = atom.coord.y - point.y
          let dz = atom.coord.z - point.z
          let dist_sq = dx * dx + dy * dy + dz * dz
          if dist_sq <= radius_sq {
            has_close_atom = true
          }
        }
        
        if has_close_atom {
          let key = "" + chain.id.to_string() + ":" + residue.resseq.to_string() + ":" + residue.icode.to_string()
          if seen.get(key) == None {
            seen.set(key, "1")
            let rid = ResidueId::new(chain.id, residue.resseq, icode=residue.icode)
            selected.push(rid)
          }
        }
      }
    }
    break // Only check first model
  }
  
  selected
}
