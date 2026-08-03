///|
/// Bio.SeqIO.PdbIO - extract amino acid sequences from PDB structures.
/// Mirrors Biopython's Bio.SeqIO.PdbIO module.
///
/// Two sources of sequence are supported:
///   - SEQRES records: the reference sequence deposited with the structure
///     (may include residues not resolved in the density).
///   - ATOM records:   the residues actually present in the structure.
///
/// Note: `parse_pdb` (in pdb_io.mbt) only retains ATOM/HETATM records, so
/// SEQRES data is only available when parsing raw PDB text via
/// `PdbSeqresParser` / `pdb_to_seqrecords(_, source="seqres")`. The
/// Structure-based `PdbSeqresRecordParser` falls back to the residues that
/// are present in the Structure (i.e. the resolved ATOM sequence).

///|
/// Standard 3-letter to 1-letter amino acid mapping used for PDB sequences.
/// Codes are upper case, matching SEQRES and ATOM record format. Includes
/// SEC (selenocysteine, U) and PYL (pyrrolysine, O).
fn pdb_aa_3_to_1_map() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=22)
  m["ALA"] = "A"
  m["ARG"] = "R"
  m["ASN"] = "N"
  m["ASP"] = "D"
  m["CYS"] = "C"
  m["GLN"] = "Q"
  m["GLU"] = "E"
  m["GLY"] = "G"
  m["HIS"] = "H"
  m["ILE"] = "I"
  m["LEU"] = "L"
  m["LYS"] = "K"
  m["MET"] = "M"
  m["PHE"] = "F"
  m["PRO"] = "P"
  m["SER"] = "S"
  m["THR"] = "T"
  m["TRP"] = "W"
  m["TYR"] = "Y"
  m["VAL"] = "V"
  m["SEC"] = "U"
  m["PYL"] = "O"
  m
}

///|
/// Reverse mapping: 1-letter amino acid code to upper-case 3-letter code.
/// Returns "UNK" for codes without a standard 3-letter equivalent.
fn pdb_aa_1_to_3_map() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=22)
  m["A"] = "ALA"
  m["R"] = "ARG"
  m["N"] = "ASN"
  m["D"] = "ASP"
  m["C"] = "CYS"
  m["Q"] = "GLN"
  m["E"] = "GLU"
  m["G"] = "GLY"
  m["H"] = "HIS"
  m["I"] = "ILE"
  m["L"] = "LEU"
  m["K"] = "LYS"
  m["M"] = "MET"
  m["F"] = "PHE"
  m["P"] = "PRO"
  m["S"] = "SER"
  m["T"] = "THR"
  m["W"] = "TRP"
  m["Y"] = "TYR"
  m["V"] = "VAL"
  m["U"] = "SEC"
  m["O"] = "PYL"
  m
}

///|
/// Convert a 3-letter amino acid code to its 1-letter equivalent.
/// Accepts any case (e.g. "ALA", "ala", "Ala"). Returns "X" for
/// non-standard or unrecognised residues (modified residues, ligands, ions).
pub fn three_to_one_aa(three_letter : String) -> String {
  let upper = three_letter.to_upper()
  match pdb_aa_3_to_1_map().get(upper) {
    Some(one) => one
    None => "X"
  }
}

///|
/// Convert a 1-letter amino acid code to its upper-case 3-letter equivalent.
/// Returns "UNK" for unrecognised codes.
fn one_to_three_aa(one_letter : String) -> String {
  match pdb_aa_1_to_3_map().get(one_letter) {
    Some(three) => three
    None => "UNK"
  }
}

///|
/// Left-pad `s` with `pad_char` to reach `width`.
fn pdb_pad_left(s : String, width : Int, pad_char : Char) -> String {
  let n = s.length()
  if n >= width {
    s
  } else {
    pad_char.to_string().repeat(width - n) + s
  }
}

///|
/// Extract the PDB ID code from the HEADER record (columns 63-66).
/// Returns an empty string when no HEADER record or no ID is present.
fn extract_pdb_id_from_content(content : String) -> String {
  let lines = content.split("\n").to_array()
  for line in lines {
    let n = line.length()
    if n >= 6 && line[0:6] == "HEADER" {
      if n >= 66 {
        let id = line[62:66].trim().to_owned()
        if id.length() > 0 {
          return id
        }
      }
      return ""
    }
  }
  ""
}

///|
/// Parse a single SEQRES line.
/// Returns `(chain_id, residue_names)` where `residue_names` is the list of
/// 3-letter residue codes on this line, in column order.
///
/// SEQRES column layout (1-indexed):
///   1-6    record name "SEQRES"
///   8-10   serial number of this line within the chain
///   12     chain identifier
///   14-17  total number of residues in the chain
///   20-22  residue name (1st)
///   24-26  residue name (2nd)
///   ...
///   68-70  residue name (13th)
pub fn parse_seqres_line(line : String) -> (String, Array[String]) {
  let n = line.length()
  let chain_id = if n > 11 {
    line.get_char(11).unwrap_or(' ').to_string()
  } else {
    " "
  }
  let residues = Array::new()
  let mut i = 19
  while i + 3 <= n {
    if i > 69 {
      break
    }
    let name = line[i:i + 3].trim().to_owned()
    if name.length() > 0 {
      residues.push(name)
    }
    i = i + 4
  }
  (chain_id, residues)
}

///|
/// Parse PDB SEQRES records and return one SeqRecord per chain.
///
/// Each chain's SEQRES residues are converted from 3-letter to 1-letter codes
/// (unknown residues become "X"). The record `id` is `<pdb_id>:<chain_id>`
/// (or just `<chain_id>` when no PDB ID is available), `name` is the chain
/// identifier.
pub fn pdb_seqres_parser(content : String) -> Array[SeqRecord] {
  let lines = content.split("\n").to_array()
  let pdb_id = extract_pdb_id_from_content(content)
  let chain_ids : Array[String] = Array::new()
  let chain_residues : Array[Array[String]] = Array::new()

  for line in lines {
    let n = line.length()
    if n < 6 {
      continue
    }
    if line[0:6] != "SEQRES" {
      continue
    }
    let (chain_id, residues) = parse_seqres_line(line.to_owned())
    if chain_id.length() == 0 {
      continue
    }
    let mut idx = -1
    for k = 0; k < chain_ids.length(); k = k + 1 {
      if chain_ids[k] == chain_id {
        idx = k
        break
      }
    }
    if idx < 0 {
      chain_ids.push(chain_id)
      let new_arr : Array[String] = Array::new()
      for r in residues {
        new_arr.push(r)
      }
      chain_residues.push(new_arr)
    } else {
      let existing = chain_residues[idx]
      for r in residues {
        existing.push(r)
      }
    }
  }

  let records = Array::new()
  for k = 0; k < chain_ids.length(); k = k + 1 {
    let chain_id = chain_ids[k]
    let residues = chain_residues[k]
    let seq_buf = StringBuilder::new()
    for r in residues {
      seq_buf.write_string(three_to_one_aa(r))
    }
    let id = if pdb_id.length() > 0 {
      pdb_id + ":" + chain_id
    } else {
      chain_id
    }
    records.push(SeqRecord::new(
      Seq::new(seq_buf.to_string()),
      id=id,
      name=chain_id,
      description="sequence from SEQRES, chain " + chain_id,
    ))
  }
  records
}

///|
/// Extract the amino acid sequence of a single chain as a 1-letter string.
/// Water molecules (HOH/WAT) are skipped; non-standard residues become "X".
pub fn extract_sequence_from_chain(chain : Chain) -> String {
  let seq_buf = StringBuilder::new()
  for residue in chain.residues {
    let name = residue.resname
    if name == "HOH" || name == "WAT" {
      continue
    }
    seq_buf.write_string(three_to_one_aa(name))
  }
  seq_buf.to_string()
}

///|
/// Build SeqRecords from a Structure's residues, one record per chain of the
/// first model. `source_label` is used in the description ("ATOM" or "SEQRES").
fn build_records_from_structure(
  structure : Structure,
  pdb_id : String,
  source_label : String,
) -> Array[SeqRecord] {
  let records = Array::new()
  if structure.models.length() == 0 {
    return records
  }
  let model = structure.models[0]
  for chain in model.chains {
    let seq_str = extract_sequence_from_chain(chain)
    let chain_id_str = chain.id.to_string()
    let id = if pdb_id.length() > 0 {
      pdb_id + ":" + chain_id_str
    } else {
      chain_id_str
    }
    records.push(SeqRecord::new(
      Seq::new(seq_str),
      id=id,
      name=chain_id_str,
      description="sequence from " + source_label + " records, chain " + chain_id_str,
    ))
  }
  records
}

///|
/// Extract amino acid sequences from a parsed Structure's ATOM records.
/// Returns one SeqRecord per chain (of the first model). Water molecules are
/// skipped and non-standard residues become "X".
pub fn pdb_atom_record_parser(structure : Structure) -> Array[SeqRecord] {
  build_records_from_structure(structure, structure.id, "ATOM")
}

///|
/// Extract sequences from a parsed Structure as if from SEQRES records.
///
/// Because `parse_pdb` does not retain SEQRES records, this falls back to the
/// residues actually present in the Structure (i.e. the resolved ATOM
/// sequence). Use `PdbSeqresParser` to read the deposited SEQRES sequence
/// directly from raw PDB text.
pub fn pdb_seqres_record_parser(structure : Structure) -> Array[SeqRecord] {
  build_records_from_structure(structure, structure.id, "SEQRES")
}

///|
/// Parse PDB ATOM/HETATM records and return one SeqRecord per chain.
///
/// Uses `parse_pdb` to build the Structure and then extracts each chain's
/// sequence. Water molecules (HOH/WAT) are skipped and non-standard residues
/// become "X". The record `id` is `<pdb_id>:<chain_id>`.
pub fn pdb_atom_parser(content : String) -> Array[SeqRecord] {
  let pdb_id = extract_pdb_id_from_content(content)
  let structure = parse_pdb(content)
  build_records_from_structure(structure, pdb_id, "ATOM")
}

///|
/// Extract the amino acid sequence of `chain_id` directly from PDB text lines.
///
/// Iterates ATOM/HETATM records belonging to `chain_id`, emitting one residue
/// per change of (resseq, icode). Residue names are converted from 3-letter
/// to 1-letter codes; water (HOH/WAT) is skipped and non-standard residues
/// become "X".
pub fn extract_sequence_from_pdb_lines(
  lines : Array[String],
  chain_id : String,
) -> String {
  let seq_buf = StringBuilder::new()
  let mut prev_resseq : Int? = None
  let mut prev_icode : Char = ' '
  for line in lines {
    let n = line.length()
    if n < 6 {
      continue
    }
    let is_atom = line[0:4] == "ATOM"
    let is_het = line[0:6] == "HETATM"
    if !(is_atom || is_het) {
      continue
    }
    if n <= 21 {
      continue
    }
    let cid = line.get_char(21).unwrap_or(' ').to_string()
    if cid != chain_id {
      continue
    }
    let resname = if n >= 20 {
      line[17:20].trim().to_owned()
    } else {
      "UNK"
    }
    if resname == "HOH" || resname == "WAT" {
      continue
    }
    let resseq = if n >= 26 {
      parse_int(line[22:26].trim().to_owned())
    } else {
      0
    }
    let icode = if n >= 27 {
      line.get_char(26).unwrap_or(' ')
    } else {
      ' '
    }
    let is_new = match prev_resseq {
      Some(p) => p != resseq || prev_icode != icode
      None => true
    }
    if is_new {
      seq_buf.write_string(three_to_one_aa(resname))
      prev_resseq = Some(resseq)
      prev_icode = icode
    }
  }
  seq_buf.to_string()
}

///|
/// Unified entry point for extracting sequences from PDB text.
///
/// `source="atom"` (default) extracts from ATOM/HETATM records via
/// `PdbAtomParser`. `source="seqres"` extracts from SEQRES records via
/// `PdbSeqresParser`.
pub fn pdb_to_seqrecords(
  content : String,
  source~ : String = "atom",
) -> Array[SeqRecord] {
  let lower = source.to_lower()
  if lower == "seqres" {
    pdb_seqres_parser(content)
  } else {
    pdb_atom_parser(content)
  }
}

///|
/// Write an array of SeqRecords back to PDB SEQRES format.
///
/// Each record's sequence is converted from 1-letter codes to upper-case
/// 3-letter codes (unknown codes become "UNK"). One chain is written per
/// record, using `record.name` as the chain identifier. Up to 13 residues
/// are placed on each SEQRES line.
pub fn write_pdb_seqrecords(records : Array[SeqRecord]) -> String {
  let buf = StringBuilder::new()
  for record in records {
    let seq_str = record.seq.to_string()
    let num_res = seq_str.length()
    if num_res == 0 {
      continue
    }
    let chain_id_str = record.name
    let chain_id = if chain_id_str.length() > 0 {
      chain_id_str.get_char(0).unwrap_or(' ')
    } else {
      ' '
    }
    let max_per_line = 13
    let mut ser = 1
    let mut idx = 0
    while idx < num_res {
      let line_buf = StringBuilder::new()
      line_buf.write_string("SEQRES ")
      line_buf.write_string(pdb_pad_left(ser.to_string(), 3, ' '))
      line_buf.write_char(' ')
      line_buf.write_char(chain_id)
      line_buf.write_char(' ')
      line_buf.write_string(pdb_pad_left(num_res.to_string(), 4, ' '))
      line_buf.write_string("  ")
      let mut count = 0
      while idx < num_res && count < max_per_line {
        let one = seq_str[idx:idx + 1].to_owned()
        line_buf.write_string(one_to_three_aa(one))
        line_buf.write_char(' ')
        idx = idx + 1
        count = count + 1
      }
      buf.write_string(line_buf.to_string())
      buf.write_char('\n')
      ser = ser + 1
    }
  }
  buf.to_string()
}
