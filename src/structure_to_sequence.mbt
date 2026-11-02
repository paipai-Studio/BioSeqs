///|
/// Structure-to-sequence conversion — extracts amino acid or nucleotide
/// sequences from 3D atomic coordinates.
///
/// Biotite provides `to_sequence()` on structure objects. This module
/// provides the same capability for BioSeqs' Chain/Structure types.

///|
/// Three-letter to one-letter amino acid mapping.
fn stoa_aa(three : String) -> Char {
  match three {
    "ALA" => 'A'
    "ARG" => 'R'
    "ASN" => 'N'
    "ASP" => 'D'
    "CYS" => 'C'
    "GLN" => 'Q'
    "GLU" => 'E'
    "GLY" => 'G'
    "HIS" => 'H'
    "ILE" => 'I'
    "LEU" => 'L'
    "LYS" => 'K'
    "MET" => 'M'
    "PHE" => 'F'
    "PRO" => 'P'
    "SER" => 'S'
    "THR" => 'T'
    "TRP" => 'W'
    "TYR" => 'Y'
    "VAL" => 'V'
    "SEC" => 'U'
    "PYL" => 'O'
    "MSE" => 'M'
    _ => 'X'
  }
}

///|
/// Three-letter to one-letter nucleotide mapping.
fn stoa_nt(three : String) -> Char {
  match three {
    "DA" => 'A'
    "DC" => 'C'
    "DG" => 'G'
    "DT" => 'T'
    "A" => 'A'
    "C" => 'C'
    "G" => 'G'
    "T" => 'T'
    "U" => 'U'
    "ADE" => 'A'
    "CYT" => 'C'
    "GUA" => 'G'
    "THY" => 'T'
    "URA" => 'U'
    _ => 'X'
  }
}

///|
/// Molecule type for sequence extraction.
pub enum SeqMoleculeType {
  Protein
  Nucleotide
  Auto // auto-detect
}

///|
/// Result of structure-to-sequence conversion.
pub struct StructureSequence {
  chain_id : Char
  sequence : String
  molecule_type : SeqMoleculeType
  residue_count : Int32
}

///|
/// Get the extracted sequence.
pub fn StructureSequence::get_sequence(self : StructureSequence) -> String {
  self.sequence
}

///|
/// Get the chain ID.
pub fn StructureSequence::get_chain_id(self : StructureSequence) -> Char {
  self.chain_id
}

///|
/// Get the residue count.
pub fn StructureSequence::get_residue_count(self : StructureSequence) -> Int32 {
  self.residue_count
}

///|
/// Detect molecule type from residue names.
fn detect_molecule_type(residues : Array[Residue]) -> SeqMoleculeType {
  let mut a = 0
  let mut n = 0
  for r in residues {
    let one = stoa_aa(r.resname)
    if one != 'X' {
      a = a + 1
    }
    let nt = stoa_nt(r.resname)
    if nt != 'X' {
      n = n + 1
    }
  }
  if a >= n {
    SeqMoleculeType::Protein
  } else {
    SeqMoleculeType::Nucleotide
  }
}

///|
/// Convert a single chain to a sequence.
///
/// Extracts the one-letter sequence from residue names in the chain.
/// For protein chains, returns amino acid sequence; for nucleotide
/// chains, returns nucleotide sequence. Auto-detects by default.
pub fn chain_to_sequence(
  chain : Chain,
  mol_type~ : SeqMoleculeType = SeqMoleculeType::Auto,
) -> StructureSequence {
  let residues = chain.get_residues()
  let actual_type = match mol_type {
    SeqMoleculeType::Auto => detect_molecule_type(residues)
    _ => mol_type
  }
  let mut seq = ""
  for r in residues {
    let c = match actual_type {
      SeqMoleculeType::Protein => stoa_aa(r.resname)
      SeqMoleculeType::Nucleotide => stoa_nt(r.resname)
      _ => stoa_aa(r.resname)
    }
    seq = seq + c.to_string()
  }
  StructureSequence::{
    chain_id: chain.get_id(),
    sequence: seq,
    molecule_type: actual_type,
    residue_count: residues.length(),
  }
}

///|
/// Convert all chains in a structure to sequences.
///
/// Returns one `StructureSequence` per chain.
pub fn structure_to_sequences(
  structure : Structure,
  mol_type~ : SeqMoleculeType = SeqMoleculeType::Auto,
) -> Array[StructureSequence] {
  let result : Array[StructureSequence] = []
  let models = structure.get_models()
  if models.length() == 0 {
    return result
  }
  let chains = models[0].get_chains()
  for chain in chains {
    result.push(chain_to_sequence(chain, mol_type~))
  }
  result
}

///|
/// Get all sequences from a structure as a FASTA string.
pub fn structure_to_fasta(structure : Structure) -> String {
  let seqs = structure_to_sequences(structure)
  let mut fasta = ""
  for s in seqs {
    fasta = fasta +
      ">" +
      structure.get_id() +
      "_chain_" +
      s.chain_id.to_string() +
      "\n"
    fasta = fasta + s.sequence + "\n"
  }
  fasta
}
