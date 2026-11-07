/// rdkit_io_ext.mbt — I/O completion (RDKit gap fill, item A8)
///
/// Molecule <-> V2000 molblock, SDF batch read/write, RXN read/write,
/// batch SMILES library (.smi) files, and a limited InChI -> SMILES
/// backbone reconstruction (formula + connections layer only; H layer and
/// charges are not reconstructed, failures return None).

///|
fn io_char(s : String, i : Int32) -> Int32 {
  s.unsafe_get(i).to_int()
}

///|
fn io_split_lines(s : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let n = s.length()
  let mut sb = StringBuilder()
  let mut i = 0
  while i < n {
    let ci = io_char(s, i)
    if ci == 10 {
      out.push(sb.to_string())
      sb = StringBuilder()
      i = i + 1
    } else if ci == 13 {
      i = i + 1
    } else {
      let _ = sb.write_char(s.unsafe_get(i).unsafe_to_char())
      i = i + 1
    }
  }
  out.push(sb.to_string())
  out
}

///|
fn io_pad_left_num(v : Int32, w : Int32) -> String {
  io_pad_left(v.to_string(), w)
}

///|
fn io_pad_left(s : String, w : Int32) -> String {
  let n = s.length()
  if n >= w {
    return s
  }
  let sb = StringBuilder()
  for _i = 0; _i < w - n; _i = _i + 1 {
    let _ = sb.write(" ")
  }
  let _ = sb.write(s)
  sb.to_string()
}

///|
fn io_pad_right(s : String, w : Int32) -> String {
  let n = s.length()
  if n >= w {
    return s
  }
  let sb = StringBuilder()
  let _ = sb.write(s)
  for _i = 0; _i < w - n; _i = _i + 1 {
    let _ = sb.write(" ")
  }
  sb.to_string()
}

///|
fn io_fmt_10_4(x : Double) -> String {
  let neg = x < 0.0
  let ax = if neg { -x } else { x }
  let scaled = (ax * 10000.0 + 0.5).floor() / 10000.0
  let int_part = scaled.floor()
  let frac = ((scaled - int_part) * 10000.0 + 0.5).floor().to_int()
  let sb = StringBuilder()
  if neg {
    let _ = sb.write("-")
  }
  let _ = sb.write(int_part.to_int().to_string())
  let _ = sb.write(".")
  let _ = sb.write(io_pad_left_num(frac, 4))
  io_pad_left(sb.to_string(), 10)
}

///|
fn io_atom_charge_code(charge : Int32) -> Int32 {
  match charge {
    3 => 1
    2 => 2
    1 => 3
    -1 => 5
    -2 => 6
    -3 => 7
    _ => 0
  }
}

///|
fn io_code_charge(code : Int32) -> Int32 {
  match code {
    1 => 3
    2 => 2
    3 => 1
    5 => -1
    6 => -2
    7 => -3
    _ => 0
  }
}

///|
/// Convert a Molecule to a V2000 molblock (coordinates are all zeros;
/// Molecule carries no 3D geometry).
pub fn rdkit_mol_to_molblock(mol : Molecule, name~ : String = "") -> String {
  let sb = StringBuilder()
  let _ = sb.write(name)
  let _ = sb.write("\n  BioSeqs          2D\n\n")
  let _ = sb.write(
    io_pad_left_num(mol.num_atoms(), 3) +
    io_pad_left_num(mol.num_bonds(), 3) +
    "  0  0  0  0  0  0  0  0999 V2000\n",
  )
  for atom in mol.atoms {
    let _ = sb.write(
      io_fmt_10_4(0.0) +
      io_fmt_10_4(0.0) +
      io_fmt_10_4(0.0) +
      " " +
      io_pad_right(atom.symbol, 3) +
      " 0" +
      io_pad_left_num(io_atom_charge_code(atom.charge), 3) +
      "  0  0  0  0  0  0  0  0  0  0\n",
    )
  }
  for bond in mol.bonds {
    let _ = sb.write(
      io_pad_left_num(bond.from + 1, 3) +
      io_pad_left_num(bond.to + 1, 3) +
      io_pad_left_num(bond.order, 3) +
      "  0\n",
    )
  }
  let _ = sb.write("M  END\n")
  sb.to_string()
}

///|
fn io_parse_int_field(line : String, start : Int32, stop : Int32) -> Int32 {
  let n = line.length()
  if n < stop {
    return 0
  }
  io_trim_parse(io_sub(line, start, stop))
}

///|
fn io_trim_parse(s : String) -> Int32 {
  let n = s.length()
  let mut start = 0
  let mut stop = n
  while start < n && io_char(s, start) == 32 {
    start = start + 1
  }
  while stop > start && io_char(s, stop - 1) == 32 {
    stop = stop - 1
  }
  if start >= stop {
    return 0
  }
  let mut neg = false
  let mut i = start
  if io_char(s, i) == 45 {
    neg = true
    i = i + 1
  }
  let mut v = 0
  while i < stop {
    let c = io_char(s, i)
    if c >= 48 && c <= 57 {
      v = v * 10 + (c - 48)
    } else {
      return 0
    }
    i = i + 1
  }
  if neg {
    -v
  } else {
    v
  }
}

///|
fn io_sub(line : String, start : Int32, stop : Int32) -> String {
  let n = line.length()
  if n <= start {
    return ""
  }
  let e = if stop > n { n } else { stop }
  line[start:e].to_owned()
}

///|
fn io_trim_str(s : String) -> String {
  let n = s.length()
  let mut start = 0
  let mut stop = n
  while start < n && io_char(s, start) == 32 {
    start = start + 1
  }
  while stop > start && io_char(s, stop - 1) == 32 {
    stop = stop - 1
  }
  s[start:stop].to_owned()
}

///|
fn io_starts_with(s : String, prefix : String) -> Bool {
  let m = prefix.length()
  if s.length() < m {
    return false
  }
  for j = 0; j < m; j = j + 1 {
    if io_char(s, j) != io_char(prefix, j) {
      return false
    }
  }
  true
}

///|
/// Parse a V2000 molblock back into a Molecule (heavy atoms + bonds).
pub fn rdkit_molblock_to_mol(text : String) -> Molecule {
  let lines = io_split_lines(text)
  if lines.length() < 4 {
    return Molecule::{ atoms: [], bonds: [], }
  }
  let n_atoms = io_parse_int_field(lines[3], 0, 3)
  let n_bonds = io_parse_int_field(lines[3], 3, 6)
  let atoms : Array[MolAtom] = Array(capacity=0)
  for i = 0; i < n_atoms; i = i + 1 {
    let idx = 4 + i
    if idx >= lines.length() {
      break
    }
    let line = lines[idx]
    let sym = io_trim_str(io_sub(line, 31, 34))
    let sym = if sym.length() == 0 { "C" } else { sym }
    let code = io_trim_parse(io_sub(line, 36, 39))
    atoms.push(MolAtom::{
      symbol: sym,
      charge: io_code_charge(code),
      aromatic: false,
    })
  }
  let bonds : Array[MolBond] = Array(capacity=0)
  for i = 0; i < n_bonds; i = i + 1 {
    let idx = 4 + n_atoms + i
    if idx >= lines.length() {
      break
    }
    let line = lines[idx]
    let from = io_parse_int_field(line, 0, 3) - 1
    let to = io_parse_int_field(line, 3, 6) - 1
    let order = io_parse_int_field(line, 6, 9)
    if from >= 0 && to >= 0 && from < atoms.length() && to < atoms.length() {
      let order = if order < 1 || order > 3 { 1 } else { order }
      bonds.push(MolBond::{ from, to, order, })
    }
  }
  mol_perceive_aromaticity(Molecule::{ atoms, bonds, })
}

///|
/// Split a text on a literal marker (records separated by "$$$$").
fn io_split_marker(text : String, marker : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let n = text.length()
  let m = marker.length()
  let mut sb = StringBuilder()
  let mut i = 0
  while i < n {
    let mut is_match = m > 0 && i + m <= n
    if is_match {
      for j = 0; j < m; j = j + 1 {
        if io_char(text, i + j) != io_char(marker, j) {
          is_match = false
        }
      }
    }
    if is_match {
      out.push(sb.to_string())
      sb = StringBuilder()
      i = i + m
    } else {
      let _ = sb.write_char(text.unsafe_get(i).unsafe_to_char())
      i = i + 1
    }
  }
  out.push(sb.to_string())
  out
}

///|
/// Read an SDF text into individual Molecules (records separated by a
/// "$$$$" line).
pub fn rdkit_read_sdf(text : String) -> Array[Molecule] {
  let out : Array[Molecule] = Array(capacity=0)
  let lines = io_split_lines(text)
  let mut sb = StringBuilder()
  for line in lines {
    if io_trim_str(line) == "$$$$" {
      let mol = rdkit_molblock_to_mol(sb.to_string())
      if mol.num_atoms() > 0 {
        out.push(mol)
      }
      sb = StringBuilder()
    } else {
      let _ = sb.write(line)
      let _ = sb.write("\n")
    }
  }
  let rest = sb.to_string()
  if io_trim_str(rest).length() > 0 {
    let mol = rdkit_molblock_to_mol(rest)
    if mol.num_atoms() > 0 {
      out.push(mol)
    }
  }
  out
}

///|
/// Write an SDF text from molecules and their names (parallel arrays).
pub fn rdkit_write_sdf(mols : Array[Molecule], names : Array[String]) -> String {
  let sb = StringBuilder()
  for i = 0; i < mols.length(); i = i + 1 {
    let name = if i < names.length() { names[i] } else { "" }
    let _ = sb.write(rdkit_mol_to_molblock(mols[i], name~))
    let _ = sb.write("$$$$\n")
  }
  sb.to_string()
}

///|
fn io_split_ws(s : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let n = s.length()
  let mut sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    let c = io_char(s, i)
    if c == 32 || c == 9 {
      if sb.to_string().length() > 0 {
        out.push(sb.to_string())
        sb = StringBuilder()
      }
    } else {
      let _ = sb.write_char(s.unsafe_get(i).unsafe_to_char())
    }
  }
  if sb.to_string().length() > 0 {
    out.push(sb.to_string())
  }
  out
}

///|
/// Parse an RXN V2000 file. Reactant/product templates are stored as
/// SMILES strings (consistent with RdkitReaction).
pub fn rdkit_parse_rxn_file(text : String) -> RdkitReaction {
  let lines = io_split_lines(text)
  let mut name = ""
  let mut i = 0
  while i < lines.length() {
    if io_trim_str(lines[i]) == "$RXN" {
      i = i + 1
      if i < lines.length() {
        name = io_trim_str(lines[i])
      }
      i = i + 3 // skip name/program/comment region
      break
    }
    i = i + 1
  }
  let mut n_reactants = 0
  let mut n_products = 0
  while i < lines.length() {
    let t = io_trim_str(lines[i])
    if t.length() > 0 && io_char(t, 0) >= 48 && io_char(t, 0) <= 57 {
      let fields = io_split_ws(t)
      if fields.length() >= 2 {
        n_reactants = io_trim_parse(fields[0])
        n_products = io_trim_parse(fields[1])
      }
      i = i + 1
      break
    }
    i = i + 1
  }
  let blocks : Array[String] = Array(capacity=0)
  while i < lines.length() {
    if io_trim_str(lines[i]) == "$MOL" {
      let mut sb = StringBuilder()
      i = i + 1
      while i < lines.length() && io_trim_str(lines[i]) != "$MOL" {
        let _ = sb.write(lines[i])
        let _ = sb.write("\n")
        i = i + 1
      }
      blocks.push(sb.to_string())
      sb = StringBuilder()
    } else {
      i = i + 1
    }
  }
  let reactants : Array[String] = Array(capacity=0)
  let products : Array[String] = Array(capacity=0)
  for k = 0; k < blocks.length(); k = k + 1 {
    let mol = rdkit_molblock_to_mol(blocks[k])
    let smi = mol_to_smiles(mol)
    if k < n_reactants {
      reactants.push(smi)
    } else if k < n_reactants + n_products {
      products.push(smi)
    }
  }
  RdkitReaction::{ reactants, products, name, }
}

///|
/// Write an RXN V2000 file from a RdkitReaction (SMILES templates are
/// parsed back into Molecules to emit $MOL blocks).
pub fn rdkit_write_rxn_file(rxn : RdkitReaction) -> String {
  let sb = StringBuilder()
  let _ = sb.write("$RXN\n")
  let _ = sb.write(rxn.name)
  let _ = sb.write("\n  BioSeqs\n\n")
  let _ = sb.write(
    "  " +
    io_pad_left_num(rxn.reactants.length(), 2) +
    " " +
    io_pad_left_num(rxn.products.length(), 2) +
    "\n",
  )
  for smi in rxn.reactants {
    let _ = sb.write("$MOL\n")
    let _ = sb.write(rdkit_mol_to_molblock(parse_smiles(smi)))
  }
  for smi in rxn.products {
    let _ = sb.write("$MOL\n")
    let _ = sb.write(rdkit_mol_to_molblock(parse_smiles(smi)))
  }
  sb.to_string()
}

///|
/// Read a batch SMILES library file: one record per line,
/// "SMILES<space/TAB>name" (name optional, '#' starts a comment line).
/// Returns (raw_smiles, name, parsed_mol) triples; blank and comment
/// lines are skipped.
pub fn rdkit_read_smi(text : String) -> Array[(String, String, Molecule)] {
  let out : Array[(String, String, Molecule)] = Array(capacity=0)
  for line in io_split_lines(text) {
    let t = io_trim_str(line)
    if t.length() == 0 || io_char(t, 0) == 35 {
      continue
    }
    let fields = io_split_ws(t)
    if fields.length() == 0 {
      continue
    }
    let smi = fields[0]
    let sb = StringBuilder()
    for i = 1; i < fields.length(); i = i + 1 {
      if i > 1 {
        let _ = sb.write(" ")
      }
      let _ = sb.write(fields[i])
    }
    let name = sb.to_string()
    out.push((smi, name, parse_smiles(smi)))
  }
  out
}

///|
/// Write a batch SMILES library file: "SMILES<TAB>name" per line.
pub fn rdkit_write_smi(records : Array[(String, String)]) -> String {
  let sb = StringBuilder()
  for (smi, name) in records {
    let _ = sb.write(smi)
    let _ = sb.write("\t")
    let _ = sb.write(name)
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
/// Limited InChI -> SMILES reconstruction.
///
/// Only the formula layer and the connections (/c) layer are interpreted:
/// the heavy-atom backbone is rebuilt, aromaticity is perceived and a
/// canonical SMILES is produced. The H layer, charges and stereo layers
/// are ignored. Returns None when the connection layer references an
/// atom outside the formula-derived backbone or cannot be tokenized.
pub fn rdkit_inchi_to_smiles(inchi : String) -> String? {
  let layers = io_split_marker(inchi, "/")
  if layers.length() == 0 {
    return None
  }
  if io_starts_with(layers[0], "InChI=1/") {
    // legacy prefix unsupported
    return None
  }
  let formula = if layers[0] == "InChI=1S" {
    if layers.length() > 1 {
      layers[1]
    } else {
      ""
    }
  } else if io_starts_with(layers[0], "InChI=1S") {
    layers[0]["InChI=1S".length():].to_owned()
  } else {
    return None
  }
  let elements = io_inchi_formula_backbone(formula)
  let n = elements.length()
  if n == 0 {
    return None
  }
  let mut conn : String? = None
  for i = 1; i < layers.length(); i = i + 1 {
    let layer = layers[i]
    if layer.length() > 1 && io_char(layer, 0) == 99 {
      conn = Some(layer[1:].to_owned())
    }
  }
  let conn_str = match conn {
    Some(s) => s
    None => ""
  }
  let atoms : Array[MolAtom] = Array(capacity=0)
  for sym in elements {
    atoms.push(MolAtom::{ symbol: sym, charge: 0, aromatic: false, })
  }
  let bonds : Array[MolBond] = Array(capacity=0)
  if conn_str.length() > 0 {
    match io_inchi_conn_tokens(conn_str) {
      None => return None
      Some(tok) =>
        match io_inchi_apply_tokens(tok, n) {
          None => return None
          Some(bs) =>
            for b in bs {
              bonds.push(b)
            }
        }
    }
  }
  let mol = mol_perceive_aromaticity(Molecule::{ atoms, bonds, })
  Some(mol_to_smiles(mol))
}

///|
/// Expand an InChI formula layer into the heavy-atom element sequence
/// (hydrogens are dropped; InChI numbers atoms in formula element order).
fn io_inchi_formula_backbone(formula : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let n = formula.length()
  let mut i = 0
  while i < n {
    let c = io_char(formula, i)
    if c >= 65 && c <= 90 {
      let sb = StringBuilder()
      let _ = sb.write_char(formula.unsafe_get(i).unsafe_to_char())
      i = i + 1
      while i < n {
        let c2 = io_char(formula, i)
        if c2 >= 97 && c2 <= 122 {
          let _ = sb.write_char(formula.unsafe_get(i).unsafe_to_char())
          i = i + 1
        } else {
          break
        }
      }
      let sym = sb.to_string()
      let mut count = 1
      let dsb = StringBuilder()
      while i < n {
        let c3 = io_char(formula, i)
        if c3 >= 48 && c3 <= 57 {
          let _ = dsb.write_char(formula.unsafe_get(i).unsafe_to_char())
          i = i + 1
        } else {
          break
        }
      }
      if dsb.to_string().length() > 0 {
        count = io_trim_parse(dsb.to_string())
      }
      if sym != "H" && sym != "D" && sym != "T" {
        for _k = 0; _k < count; _k = _k + 1 {
          out.push(sym)
        }
      }
    } else {
      i = i + 1
    }
  }
  out
}

///|
/// Tokenize an InChI connections layer body: integers separated by '-',
/// '(' opens a branch, ')' closes it.
fn io_inchi_conn_tokens(body : String) -> Array[Int32]? {
  let out : Array[Int32] = Array(capacity=0)
  let n = body.length()
  let mut i = 0
  while i < n {
    let c = io_char(body, i)
    if c == 40 || c == 41 {
      out.push(if c == 40 { 1000000 } else { 1000001 })
      i = i + 1
    } else if c == 45 {
      i = i + 1
    } else if c >= 48 && c <= 57 {
      let mut v = 0
      while i < n {
        let c2 = io_char(body, i)
        if c2 >= 48 && c2 <= 57 {
          v = v * 10 + (c2 - 48)
          i = i + 1
        } else {
          break
        }
      }
      out.push(v)
    } else {
      return None
    }
  }
  Some(out)
}

///|
/// Interpret tokenized connections: positive = new atom (bonded to the
/// current atom, current moves there), negative = ring closure to |t|
/// (current stays), 1000000/1000001 = open/close a branch.
fn io_inchi_apply_tokens(tokens : Array[Int32], n : Int32) -> Array[MolBond]? {
  let bonds : Array[MolBond] = Array(capacity=0)
  let mut cur = -1
  let stack : Array[Int32] = Array(capacity=0)
  for t in tokens {
    if t == 1000000 {
      stack.push(cur)
    } else if t == 1000001 {
      if stack.length() > 0 {
        cur = stack[stack.length() - 1]
        let _ = stack.pop()
      }
    } else if t < 0 {
      let target = -t
      if target < 1 || target > n || cur < 1 || cur > n {
        return None
      }
      bonds.push(MolBond::{ from: cur - 1, to: target - 1, order: 1, })
    } else {
      if t < 1 || t > n {
        return None
      }
      if cur >= 1 && cur <= n {
        bonds.push(MolBond::{ from: cur - 1, to: t - 1, order: 1, })
      }
      cur = t
    }
  }
  if stack.length() > 0 {
    return None
  }
  Some(bonds)
}
