///|
/// PDB format parser and writer for protein structures.

///|
/// Parse a PDB file content into a Structure.
pub fn parse_pdb(content : String) -> Structure {
  let lines = content.split("\n").to_array()
  let models = Array::new()
  let mut current_model_chains = Array::new()
  let mut current_chain_residues = Array::new()
  let mut current_residue_atoms = Array::new()
  let mut current_resname = ""
  let mut current_chainid = ' '
  let mut current_resseq = 1
  let mut current_hetfield = " "

  for line in lines {
    let n = line.length()
    if n < 6 {
      continue
    }

    let record_type = line[0:6]
    if record_type[0:4] != "ATOM" &&
      record_type[0:6] != "HETATM" &&
      record_type[0:5] != "MODEL" &&
      record_type[0:3] != "TER" &&
      record_type[0:6] != "ENDMDL" {
      continue
    }

    if record_type == "MODEL" {
      if current_residue_atoms.length() > 0 {
        let r = Residue::new(
          resname=current_resname,
          chainid=current_chainid,
          resseq=current_resseq,
          atoms=current_residue_atoms,
          hetfield=current_hetfield,
        )
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()
      }
      if current_chain_residues.length() > 0 {
        let c = Chain::new(id=current_chainid, residues=current_chain_residues)
        current_model_chains.push(c)
        current_chain_residues = Array::new()
      }
      if current_model_chains.length() > 0 {
        let m = Model::new(id=models.length() + 1, chains=current_model_chains)
        models.push(m)
        current_model_chains = Array::new()
      }
    } else if record_type[0:4] == "ATOM" || record_type == "HETATM" {
      let atom = parse_atom_line(line, record_type == "HETATM")

      if current_residue_atoms.length() > 0 &&
        (atom.chainid != current_chainid || atom.resseq != current_resseq) {
        let r = Residue::new(
          resname=current_resname,
          chainid=current_chainid,
          resseq=current_resseq,
          atoms=current_residue_atoms,
          hetfield=current_hetfield,
        )
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()

        if atom.chainid != current_chainid {
          let c = Chain::new(
            id=current_chainid,
            residues=current_chain_residues,
          )
          current_model_chains.push(c)
          current_chain_residues = Array::new()
        }
      }

      current_resname = atom.resname
      current_chainid = atom.chainid
      current_resseq = atom.resseq
      current_hetfield = if record_type == "HETATM" { "H" } else { " " }
      current_residue_atoms.push(atom)
    } else if record_type == "TER" {
      if current_residue_atoms.length() > 0 {
        let r = Residue::new(
          resname=current_resname,
          chainid=current_chainid,
          resseq=current_resseq,
          atoms=current_residue_atoms,
          hetfield=current_hetfield,
        )
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()
      }
      if current_chain_residues.length() > 0 {
        let c = Chain::new(id=current_chainid, residues=current_chain_residues)
        current_model_chains.push(c)
        current_chain_residues = Array::new()
      }
    } else if record_type == "ENDMDL" {
      if current_residue_atoms.length() > 0 {
        let r = Residue::new(
          resname=current_resname,
          chainid=current_chainid,
          resseq=current_resseq,
          atoms=current_residue_atoms,
          hetfield=current_hetfield,
        )
        current_chain_residues.push(r)
        current_residue_atoms = Array::new()
      }
      if current_chain_residues.length() > 0 {
        let c = Chain::new(id=current_chainid, residues=current_chain_residues)
        current_model_chains.push(c)
        current_chain_residues = Array::new()
      }
      if current_model_chains.length() > 0 {
        let m = Model::new(id=models.length() + 1, chains=current_model_chains)
        models.push(m)
        current_model_chains = Array::new()
      }
    }
  }

  if current_residue_atoms.length() > 0 {
    let r = Residue::new(
      resname=current_resname,
      chainid=current_chainid,
      resseq=current_resseq,
      atoms=current_residue_atoms,
      hetfield=current_hetfield,
    )
    current_chain_residues.push(r)
  }
  if current_chain_residues.length() > 0 {
    let c = Chain::new(id=current_chainid, residues=current_chain_residues)
    current_model_chains.push(c)
  }
  if current_model_chains.length() > 0 {
    let m = Model::new(id=models.length() + 1, chains=current_model_chains)
    models.push(m)
  }

  Structure::new(id="", models~)
}

///|
/// Parse an ATOM or HETATM line into an Atom.
fn parse_atom_line(line : StringView, _is_hetatm : Bool) -> Atom {
  let n = line.length()

  let name = if n >= 16 { pdb_trim_string(line[12:16]) } else { "" }
  let resname = if n >= 20 { pdb_trim_string(line[17:20]) } else { "UNK" }
  let chainid = if n > 21 { line.get_char(21).unwrap_or(' ') } else { ' ' }
  let resseq = if n >= 26 { pdb_parse_int_fast(line[22:26]) } else { 1 }

  let x = if n >= 38 { pdb_parse_double_fast(line[30:38]) } else { 0.0 }
  let y = if n >= 46 { pdb_parse_double_fast(line[38:46]) } else { 0.0 }
  let z = if n >= 54 { pdb_parse_double_fast(line[46:54]) } else { 0.0 }

  let occupancy = if n >= 60 { pdb_parse_double_fast(line[54:60]) } else { 1.0 }
  let bfactor = if n >= 66 { pdb_parse_double_fast(line[60:66]) } else { 0.0 }

  let element = if n >= 78 { pdb_trim_string(line[76:78]) } else { "" }

  Atom::new(
    name~,
    coord=Vector3::new(x, y, z),
    resname~,
    chainid~,
    resseq~,
    occupancy~,
    bfactor~,
    element~,
    charge=0,
  )
}

///|
/// Trim whitespace from both ends of a string view and return String.
fn pdb_trim_string(s : StringView) -> String {
  let n = s.length()
  if n == 0 {
    return ""
  }

  let mut start = 0
  while start < n {
    let c = s.unsafe_get(start).to_int()
    if c != 32 && c != 9 && c != 10 && c != 13 {
      break
    }
    start = start + 1
  }

  if start >= n {
    return ""
  }

  let mut end = n - 1
  while end >= start {
    let c = s.unsafe_get(end).to_int()
    if c != 32 && c != 9 && c != 10 && c != 13 {
      break
    }
    end = end - 1
  }

  s[start:end + 1].to_string()
}

///|
/// Fast integer parser for PDB format (no error checking).
fn pdb_parse_int_fast(s : StringView) -> Int {
  let n = s.length()
  if n == 0 {
    return 0
  }

  let mut result = 0
  let zero = '0'.to_int()
  let mut i = 0
  let mut is_negative = false

  if s.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }

  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c >= zero && c <= zero + 9 {
      result = result * 10 + (c - zero)
    }
    i = i + 1
  }

  if is_negative {
    -result
  } else {
    result
  }
}

///|
/// Parse an integer from a string.
fn pdb_parse_int(s : String) -> Int? {
  let n = s.length()
  if n == 0 {
    None
  } else {
    let mut result = 0
    let zero = '0'.to_int()
    let mut i = 0
    let mut is_negative = false

    if s.unsafe_get(0).to_int() == 45 {
      is_negative = true
      i = 1
    }

    while i < n {
      let c = s.unsafe_get(i).to_int()
      if c < zero || c > zero + 9 {
        return None
      }
      result = result * 10 + (c - zero)
      i = i + 1
    }

    if is_negative {
      Some(-result)
    } else {
      Some(result)
    }
  }
}

///|
/// Fast double parser for PDB format (no error checking).
fn pdb_parse_double_fast(s : StringView) -> Double {
  let n = s.length()
  if n == 0 {
    return 0.0
  }

  let mut result = 0.0
  let mut decimal = 0.0
  let mut decimal_places = 0
  let mut is_negative = false
  let zero = '0'.to_int()
  let mut i = 0

  if s.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }

  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c == 46 {
      i = i + 1
      while i < n {
        let dc = s.unsafe_get(i).to_int()
        if dc >= zero && dc <= zero + 9 {
          decimal = decimal * 10.0 + (dc - zero).to_double()
          decimal_places = decimal_places + 1
        }
        i = i + 1
      }
      break
    }
    if c >= zero && c <= zero + 9 {
      result = result * 10.0 + (c - zero).to_double()
    }
    i = i + 1
  }

  while decimal_places > 0 {
    decimal = decimal / 10.0
    decimal_places = decimal_places - 1
  }

  let final_result = result + decimal
  if is_negative {
    -final_result
  } else {
    final_result
  }
}

///|
/// Parse a double from a string.
fn pdb_parse_double(s : String) -> Double? {
  let n = s.length()
  if n == 0 {
    None
  } else {
    let mut result = 0.0
    let mut decimal = 0.0
    let mut decimal_places = 0
    let mut is_negative = false
    let zero = '0'.to_int()
    let mut i = 0

    if s.unsafe_get(0).to_int() == 45 {
      is_negative = true
      i = 1
    }

    while i < n {
      let c = s.unsafe_get(i).to_int()
      if c == 46 {
        i = i + 1
        while i < n {
          let dc = s.unsafe_get(i).to_int()
          if dc < zero || dc > zero + 9 {
            break
          }
          decimal = decimal * 10.0 + (dc - zero).to_double()
          decimal_places = decimal_places + 1
          i = i + 1
        }
        break
      }
      if c < zero || c > zero + 9 {
        return None
      }
      result = result * 10.0 + (c - zero).to_double()
      i = i + 1
    }

    while decimal_places > 0 {
      decimal = decimal / 10.0
      decimal_places = decimal_places - 1
    }

    let final_result = result + decimal
    if is_negative {
      Some(-final_result)
    } else {
      Some(final_result)
    }
  }
}

///|
/// Write a Structure to PDB format string.
pub fn write_pdb(structure : Structure) -> String {
  let buf = StringBuilder::new()

  let mut model_idx = 0
  for model in structure.models {
    model_idx = model_idx + 1
    buf.write_string("MODEL     ")
    buf.write_string(model_idx.to_string())
    buf.write_char('\n')

    for chain in model.chains {
      for residue in chain.residues {
        for atom in residue.atoms {
          buf.write_string(write_atom_line(atom, residue.hetfield == "H"))
          buf.write_char('\n')
        }
        buf.write_string("TER\n")
      }
    }

    buf.write_string("ENDMDL\n")
  }

  buf.write_string("END\n")
  buf.to_string()
}

///|
/// Write an atom to a PDB ATOM or HETATM line.
fn write_atom_line(atom : Atom, is_hetatm : Bool) -> String {
  let buf = StringBuilder::new()

  buf.write_string(if is_hetatm { "HETATM" } else { "ATOM  " })

  let serial_str = pad_left("1", 5, ' ')
  buf.write_string(serial_str)
  buf.write_char(' ')

  let name_str = pad_right(atom.name, 4, ' ')
  buf.write_string(name_str)

  buf.write_char(' ')
  buf.write_char(' ')

  let resname_str = pad_right(atom.resname, 3, ' ')
  buf.write_string(resname_str)
  buf.write_char(' ')

  buf.write_char(atom.chainid)
  buf.write_char(' ')

  let resseq_str = pad_left(atom.resseq.to_string(), 4, ' ')
  buf.write_string(resseq_str)

  buf.write_char(' ')
  buf.write_string("   ")

  let x_str = format_double(atom.coord.x, 8, 3)
  buf.write_string(x_str)

  let y_str = format_double(atom.coord.y, 8, 3)
  buf.write_string(y_str)

  let z_str = format_double(atom.coord.z, 8, 3)
  buf.write_string(z_str)

  let occ_str = format_double(atom.occupancy, 6, 2)
  buf.write_string(occ_str)

  let bfactor_str = format_double(atom.bfactor, 6, 2)
  buf.write_string(bfactor_str)

  buf.write_string("      ")

  let elem_str = pad_right(atom.element, 2, ' ')
  buf.write_string(elem_str)

  buf.write_string("  ")

  buf.to_string()
}

///|
/// Pad string to the left.
fn pad_left(s : String, width : Int, pad_char : Char) -> String {
  let n = s.length()
  if n >= width {
    s
  } else {
    pad_char.to_string().repeat(width - n) + s
  }
}

///|
/// Pad string to the right.
fn pad_right(s : String, width : Int, pad_char : Char) -> String {
  let n = s.length()
  if n >= width {
    s
  } else {
    s + pad_char.to_string().repeat(width - n)
  }
}

///|
/// Format a double with specified width and precision.
fn format_double(value : Double, width : Int, precision : Int) -> String {
  let s = value.to_string()
  let parts = s.split(".").to_array()
  let int_part = if parts.length() > 0 { parts[0].to_string() } else { "0" }
  let frac_part = if parts.length() > 1 { parts[1].to_string() } else { "" }

  let mut padded_frac = frac_part
  if padded_frac.length() < precision {
    padded_frac = padded_frac + "0".repeat(precision - padded_frac.length())
  } else if padded_frac.length() > precision {
    padded_frac = padded_frac[0:precision].to_string()
  }

  let result = int_part + "." + padded_frac
  pad_left(result, width, ' ')
}
