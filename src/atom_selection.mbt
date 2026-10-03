///|
/// ProDy-style Atom Selection Language.
///
/// Supports queries like:
///   "calpha"                      — all C-alpha atoms
///   "backbone"                    — N, CA, C, O atoms
///   "resnum 1:50"                 — residues 1-50
///   "calpha and resnum 1:50"      — intersection
///   "calpha or backbone"          — union
///   "not water"                   — complement
///   "name CA CB"                  — atoms named CA or CB
///   "chain A and resname ALA"     — chain A alanines
///   "within 5.0 of calpha"        — spatial selection
///   "(ca or cb) and not water"    — grouped expression

// ===== AST =====

///|
pub(all) enum SelExpr {
  All
  None_
  ByName(Array[String])
  ByResName(Array[String])
  ByResNum(Array[(Int, Int)])  // list of (start, end) inclusive ranges
  ByChain(Array[Char])
  ByElement(Array[String])
  Calpha
  Backbone
  Sidechain
  Protein
  Water
  Nucleic
  Hetero
  Hydrogen
  Heavy
  And(SelExpr, SelExpr)
  Or(SelExpr, SelExpr)
  Not(SelExpr)
  Within(Double, SelExpr)
}

// ===== Tokenizer =====

///|
pub(all) enum Token {
  Word(String)
  Number(Double)
  Range(Int, Int)
  LParen
  RParen
  EOF
}

///|
fn tokenize(s : String) -> Array[Token] {
  let tokens : Array[Token] = Array(capacity=0)
  let chars = s.to_array()
  let n = chars.length()
  let mut i = 0
  while i < n {
    let ch = chars[i]
    // Skip whitespace
    if ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' {
      i = i + 1
      continue
    }
    // Parentheses
    if ch == '(' {
      tokens.push(Token::LParen)
      i = i + 1
      continue
    }
    if ch == ')' {
      tokens.push(Token::RParen)
      i = i + 1
      continue
    }
    // Read a word/number/range
    let buf = StringBuilder()
    while i < n {
      let c = chars[i]
      if c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '(' || c == ')' {
        break
      }
      buf.write_char(c)
      i = i + 1
    }
    let word = buf.to_string()
    // Check if it's a range (contains ':')
    let colon_idx = sel_find_char(word, ':')
    if colon_idx >= 0 {
      let left = substring(word, 0, colon_idx)
      let right = substring(word, colon_idx + 1, word.length() - colon_idx - 1)
      let start = safe_parse_int(left)
      let end = safe_parse_int(right)
      tokens.push(Token::Range(start, end))
    } else if is_number(word) {
      tokens.push(Token::Number(safe_parse_double(word)))
    } else {
      tokens.push(Token::Word(word))
    }
  }
  tokens.push(Token::EOF)
  tokens
}

///|
fn sel_find_char(s : String, ch : Char) -> Int {
  let chars = s.to_array()
  for i = 0; i < chars.length(); i = i + 1 {
    if chars[i] == ch {
      return i
    }
  }
  -1
}

///|
fn is_number(s : String) -> Bool {
  if s.length() == 0 {
    return false
  }
  let chars = s.to_array()
  let mut has_digit = false
  for i = 0; i < chars.length(); i = i + 1 {
    let c = chars[i]
    if c >= '0' && c <= '9' {
      has_digit = true
    } else if c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E' {
      ()
    } else {
      return false
    }
  }
  has_digit
}

// ===== Parser (recursive descent) =====
// Precedence: not > and > or
// Grammar:
//   expr   := or_expr
//   or_expr := and_expr ('or' and_expr)*
//   and_expr := not_expr ('and' not_expr)*
//   not_expr := 'not' not_expr | primary
//   primary := '(' expr ')' | atom_sel | spatial_sel

///|
struct SelParser {
  tokens : Array[Token]
  pos : Ref[Int]
}

///|
fn parse_selection(s : String) -> SelExpr? {
  let tokens = tokenize(s)
  let p = SelParser::{ tokens, pos: Ref(0) }
  let result = parse_or(p)
  // Must consume all tokens (except EOF)
  let pos = p.pos.val
  if pos < p.tokens.length() - 1 {
    return None
  }
  Some(result)
}

///|
fn parse_or(p : SelParser) -> SelExpr {
  let left = parse_and(p)
  let mut result = left
  while true {
    let tok = peek(p)
    match tok {
      Token::Word(w) => {
        if w == "or" {
          advance(p)
          let right = parse_and(p)
          result = SelExpr::Or(result, right)
        } else {
          break
        }
      }
      _ => break
    }
  }
  result
}

///|
fn parse_and(p : SelParser) -> SelExpr {
  let left = parse_not(p)
  let mut result = left
  while true {
    let tok = peek(p)
    match tok {
      Token::Word(w) => {
        if w == "and" {
          advance(p)
          let right = parse_not(p)
          result = SelExpr::And(result, right)
        } else {
          break
        }
      }
      _ => break
    }
  }
  result
}

///|
fn parse_not(p : SelParser) -> SelExpr {
  let tok = peek(p)
  match tok {
    Token::Word(w) => {
      if w == "not" {
        advance(p)
        let inner = parse_not(p)
        SelExpr::Not(inner)
      } else {
        parse_primary(p)
      }
    }
    _ => parse_primary(p)
  }
}

///|
fn parse_primary(p : SelParser) -> SelExpr {
  let tok = peek(p)
  match tok {
    Token::LParen => {
      advance(p)
      let expr = parse_or(p)
      let next = peek(p)
      match next {
        Token::RParen => advance(p)
        _ => ()
      }
      expr
    }
    Token::Word(w) => parse_keyword(p, w)
    _ => {
      advance(p)
      SelExpr::None_
    }
  }
}

///|
fn parse_keyword(p : SelParser, w : String) -> SelExpr {
  let lw = to_lowercase(w)
  match lw {
    "all" => {
      advance(p)
      SelExpr::All
    }
    "none" => {
      advance(p)
      SelExpr::None_
    }
    "calpha" | "ca" => {
      advance(p)
      SelExpr::Calpha
    }
    "backbone" | "bb" => {
      advance(p)
      SelExpr::Backbone
    }
    "sidechain" | "sc" => {
      advance(p)
      SelExpr::Sidechain
    }
    "protein" | "aa" => {
      advance(p)
      SelExpr::Protein
    }
    "water" | "hoh" => {
      advance(p)
      SelExpr::Water
    }
    "nucleic" | "na" => {
      advance(p)
      SelExpr::Nucleic
    }
    "hetero" => {
      advance(p)
      SelExpr::Hetero
    }
    "hydrogen" | "h" => {
      advance(p)
      SelExpr::Hydrogen
    }
    "heavy" | "noh" => {
      advance(p)
      SelExpr::Heavy
    }
    "name" => {
      advance(p)
      let names = read_words(p)
      SelExpr::ByName(names)
    }
    "resname" => {
      advance(p)
      let names = read_words(p)
      SelExpr::ByResName(names)
    }
    "resnum" => {
      advance(p)
      let ranges = read_ranges(p)
      SelExpr::ByResNum(ranges)
    }
    "chain" => {
      advance(p)
      let chains = read_chains(p)
      SelExpr::ByChain(chains)
    }
    "element" => {
      advance(p)
      let elems = read_words(p)
      SelExpr::ByElement(elems)
    }
    "within" => {
      advance(p)
      // Expect: <distance> of <selection>
      let dist = read_number(p)
      // Expect "of"
      let tok = peek(p)
      match tok {
        Token::Word(w2) => {
          if to_lowercase(w2) == "of" {
            advance(p)
          }
        }
        _ => ()
      }
      let inner = parse_primary(p)
      SelExpr::Within(dist, inner)
    }
    _ => {
      advance(p)
      SelExpr::None_
    }
  }
}

///|
fn read_words(p : SelParser) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  while true {
    let tok = peek(p)
    match tok {
      Token::Word(w) => {
        let lw = to_lowercase(w)
        // Stop at keywords
        if is_logical_op(lw) {
          break
        }
        advance(p)
        result.push(w)
      }
      _ => break
    }
  }
  if result.length() == 0 {
    result.push("")
  }
  result
}

///|
fn read_chains(p : SelParser) -> Array[Char] {
  let result : Array[Char] = Array(capacity=0)
  while true {
    let tok = peek(p)
    match tok {
      Token::Word(w) => {
        let lw = to_lowercase(w)
        if is_logical_op(lw) {
          break
        }
        advance(p)
        // Take first char of the word
        let chars = w.to_array()
        if chars.length() > 0 {
          result.push(chars[0])
        }
      }
      _ => break
    }
  }
  if result.length() == 0 {
    result.push(' ')
  }
  result
}

///|
fn read_ranges(p : SelParser) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = Array(capacity=0)
  while true {
    let tok = peek(p)
    match tok {
      Token::Range(start, end) => {
        advance(p)
        result.push((start, end))
      }
      Token::Number(n) => {
        advance(p)
        let v = double_to_int(n)
        result.push((v, v))
      }
      Token::Word(w) => {
        let lw = to_lowercase(w)
        if is_logical_op(lw) {
          break
        }
        // Try to parse as number
        let v = safe_parse_int(w)
        advance(p)
        result.push((v, v))
      }
      _ => break
    }
  }
  if result.length() == 0 {
    result.push((0, 0))
  }
  result
}

///|
fn read_number(p : SelParser) -> Double {
  let tok = peek(p)
  match tok {
    Token::Number(n) => {
      advance(p)
      n
    }
    Token::Word(w) => {
      advance(p)
      safe_parse_double(w)
    }
    _ => 0.0
  }
}

///|
fn is_keyword(w : String) -> Bool {
  match w {
    "and" | "or" | "not" | "all" | "none" |
    "calpha" | "ca" | "backbone" | "bb" | "sidechain" | "sc" |
    "protein" | "aa" | "water" | "hoh" | "nucleic" | "na" |
    "hetero" | "hydrogen" | "h" | "heavy" | "noh" |
    "name" | "resname" | "resnum" | "chain" | "element" |
    "within" | "of" | "same" | "residue" | "as" => true
    _ => false
  }
}

///|
fn peek(p : SelParser) -> Token {
  let pos = p.pos.val
  if pos < p.tokens.length() {
    p.tokens[pos]
  } else {
    Token::EOF
  }
}

///|
fn advance(p : SelParser) -> Unit {
  let pos = p.pos.val
  p.pos.val = pos + 1
}

// ===== Evaluator =====

/// Amino acid residue names (3-letter codes).
fn is_aa(resname : String) -> Bool {
  match resname {
    "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "GLN" | "GLU" |
    "GLY" | "HIS" | "ILE" | "LEU" | "LYS" | "MET" | "PHE" |
    "PRO" | "SER" | "THR" | "TRP" | "TYR" | "VAL" |
    "MSE" | "SEC" | "PYL" => true
    _ => false
  }
}

///|
fn is_water(resname : String) -> Bool {
  match resname {
    "HOH" | "WAT" | "TIP" | "TIP3" | "TIP4" | "DOD" | "SOL" => true
    _ => false
  }
}

///|
fn is_nucleic(resname : String) -> Bool {
  match resname {
    "DA" | "DC" | "DG" | "DT" | "DU" |
    "A" | "C" | "G" | "T" | "U" |
    "ADE" | "CYT" | "GUA" | "THY" | "URA" => true
    _ => false
  }
}

///|
fn is_backbone_name(name : String) -> Bool {
  match name {
    "N" | "CA" | "C" | "O" | "OXT" | "H" | "HA" => true
    _ => false
  }
}

///|
/// Evaluate a selection expression against a list of atoms.
/// Returns a boolean array (true = selected).
pub fn eval_selection(atoms : Array[Atom], expr : SelExpr) -> Array[Bool] {
  let n = atoms.length()
  let result = Array::make(n, false)
  eval_expr(atoms, expr, result)
  result
}

///|
fn eval_expr(atoms : Array[Atom], expr : SelExpr, result : Array[Bool]) -> Unit {
  let n = atoms.length()
  match expr {
    SelExpr::All => {
      for i = 0; i < n; i = i + 1 {
        result[i] = true
      }
    }
    SelExpr::None_ => ()
    SelExpr::ByName(names) => {
      for i = 0; i < n; i = i + 1 {
        for name in names {
          if atoms[i].name == name {
            result[i] = true
          }
        }
      }
    }
    SelExpr::ByResName(names) => {
      for i = 0; i < n; i = i + 1 {
        for name in names {
          if atoms[i].resname == name {
            result[i] = true
          }
        }
      }
    }
    SelExpr::ByResNum(ranges) => {
      for i = 0; i < n; i = i + 1 {
        let rn = atoms[i].resseq
        for (start, end) in ranges {
          if rn >= start && rn <= end {
            result[i] = true
          }
        }
      }
    }
    SelExpr::ByChain(chains) => {
      for i = 0; i < n; i = i + 1 {
        for ch in chains {
          if atoms[i].chainid == ch {
            result[i] = true
          }
        }
      }
    }
    SelExpr::ByElement(elems) => {
      for i = 0; i < n; i = i + 1 {
        for elem in elems {
          if atoms[i].element == elem {
            result[i] = true
          }
        }
      }
    }
    SelExpr::Calpha => {
      for i = 0; i < n; i = i + 1 {
        if atoms[i].name == "CA" && is_aa(atoms[i].resname) {
          result[i] = true
        }
      }
    }
    SelExpr::Backbone => {
      for i = 0; i < n; i = i + 1 {
        if is_backbone_name(atoms[i].name) && is_aa(atoms[i].resname) {
          result[i] = true
        }
      }
    }
    SelExpr::Sidechain => {
      for i = 0; i < n; i = i + 1 {
        if is_aa(atoms[i].resname) && !is_backbone_name(atoms[i].name) {
          result[i] = true
        }
      }
    }
    SelExpr::Protein => {
      for i = 0; i < n; i = i + 1 {
        if is_aa(atoms[i].resname) {
          result[i] = true
        }
      }
    }
    SelExpr::Water => {
      for i = 0; i < n; i = i + 1 {
        if is_water(atoms[i].resname) {
          result[i] = true
        }
      }
    }
    SelExpr::Nucleic => {
      for i = 0; i < n; i = i + 1 {
        if is_nucleic(atoms[i].resname) {
          result[i] = true
        }
      }
    }
    SelExpr::Hetero => {
      for i = 0; i < n; i = i + 1 {
        if !is_aa(atoms[i].resname) && !is_water(atoms[i].resname) {
          result[i] = true
        }
      }
    }
    SelExpr::Hydrogen => {
      for i = 0; i < n; i = i + 1 {
        if atoms[i].element == "H" || char_at(atoms[i].name, 0) == "H" {
          result[i] = true
        }
      }
    }
    SelExpr::Heavy => {
      for i = 0; i < n; i = i + 1 {
        if !(atoms[i].element == "H" || char_at(atoms[i].name, 0) == "H") {
          result[i] = true
        }
      }
    }
    SelExpr::And(left, right) => {
      let l = Array::make(n, false)
      let r = Array::make(n, false)
      eval_expr(atoms, left, l)
      eval_expr(atoms, right, r)
      for i = 0; i < n; i = i + 1 {
        result[i] = l[i] && r[i]
      }
    }
    SelExpr::Or(left, right) => {
      let l = Array::make(n, false)
      let r = Array::make(n, false)
      eval_expr(atoms, left, l)
      eval_expr(atoms, right, r)
      for i = 0; i < n; i = i + 1 {
        result[i] = l[i] || r[i]
      }
    }
    SelExpr::Not(inner) => {
      let tmp = Array::make(n, false)
      eval_expr(atoms, inner, tmp)
      for i = 0; i < n; i = i + 1 {
        result[i] = !tmp[i]
      }
    }
    SelExpr::Within(dist, inner) => {
      // First evaluate inner selection
      let inner_sel = Array::make(n, false)
      eval_expr(atoms, inner, inner_sel)
      // For each atom not in inner_sel, check if within dist of any inner atom
      for i = 0; i < n; i = i + 1 {
        if inner_sel[i] {
          result[i] = true
        } else {
          for j = 0; j < n; j = j + 1 {
            if inner_sel[j] {
              let dx = atoms[i].coord.x - atoms[j].coord.x
              let dy = atoms[i].coord.y - atoms[j].coord.y
              let dz = atoms[i].coord.z - atoms[j].coord.z
              let d2 = dx * dx + dy * dy + dz * dz
              if d2 <= dist * dist {
                result[i] = true
              }
            }
          }
        }
      }
    }
  }
}

// ===== Public API =====

///|
/// Parse a selection string into an AST.
/// Returns None if parsing fails.
pub fn parse_sel(s : String) -> SelExpr? {
  parse_selection(s)
}

///|
/// Select atoms from a chain using a selection string.
/// Returns the selected atoms.
pub fn Chain::select(self : Chain, sel : String) -> Array[Atom] {
  let atoms = self.get_atoms()
  select_atoms(atoms, sel)
}

///|
/// Select atoms from an atom array using a selection string.
/// Returns the selected atoms.
pub fn select_atoms(atoms : Array[Atom], sel : String) -> Array[Atom] {
  match parse_selection(sel) {
    Some(expr) => {
      let mask = eval_selection(atoms, expr)
      let result : Array[Atom] = Array(capacity=0)
      for i = 0; i < atoms.length(); i = i + 1 {
        if mask[i] {
          result.push(atoms[i])
        }
      }
      result
    }
    None => Array(capacity=0)
  }
}

///|
/// Get a boolean mask for atoms using a selection string.
pub fn select_mask(atoms : Array[Atom], sel : String) -> Array[Bool] {
  match parse_selection(sel) {
    Some(expr) => eval_selection(atoms, expr)
    None => Array::make(atoms.length(), false)
  }
}

///|
/// Count selected atoms.
pub fn count_selected(atoms : Array[Atom], sel : String) -> Int {
  match parse_selection(sel) {
    Some(expr) => {
      let mask = eval_selection(atoms, expr)
      let mut count = 0
      for i = 0; i < mask.length(); i = i + 1 {
        if mask[i] {
          count = count + 1
        }
      }
      count
    }
    None => 0
  }
}

///|
/// Get indices of selected atoms (0-based).
pub fn select_indices(atoms : Array[Atom], sel : String) -> Array[Int] {
  match parse_selection(sel) {
    Some(expr) => {
      let mask = eval_selection(atoms, expr)
      let result : Array[Int] = Array(capacity=0)
      for i = 0; i < mask.length(); i = i + 1 {
        if mask[i] {
          result.push(i)
        }
      }
      result
    }
    None => Array(capacity=0)
  }
}

///|
/// Select residues from a chain using a selection string.
/// A residue is selected if any of its atoms are selected.
pub fn Chain::select_residues(self : Chain, sel : String) -> Array[Residue] {
  let atoms = self.get_atoms()
  match parse_selection(sel) {
    Some(expr) => {
      let mask = eval_selection(atoms, expr)
      let result : Array[Residue] = Array(capacity=0)
      let residues = self.get_residues()
      let mut atom_idx = 0
      for r = 0; r < residues.length(); r = r + 1 {
        let num_atoms = residues[r].get_num_atoms()
        let mut selected = false
        for a = 0; a < num_atoms; a = a + 1 {
          if atom_idx < mask.length() && mask[atom_idx] {
            selected = true
          }
          atom_idx = atom_idx + 1
        }
        if selected {
          result.push(residues[r])
        }
      }
      result
    }
    None => Array(capacity=0)
  }
}

///|
fn double_to_int(d : Double) -> Int {
  let mut result = 0
  let mut x = d
  if x < 0.0 {
    x = -x
  }
  while x >= 1.0 {
    result = result + 1
    x = x - 1.0
  }
  if d < 0.0 {
    -result
  } else {
    result
  }
}

///|
fn safe_parse_double(s : String) -> Double {
  let chars = s.to_array()
  let n = chars.length()
  if n == 0 {
    return 0.0
  }
  let mut i = 0
  let mut neg = false
  if chars[0] == '-' {
    neg = true
    i = 1
  } else if chars[0] == '+' {
    i = 1
  }
  let mut int_part = 0.0
  while i < n && chars[i] >= '0' && chars[i] <= '9' {
    int_part = int_part * 10.0 + (chars[i].to_int() - '0'.to_int()).to_double()
    i = i + 1
  }
  let mut frac_part = 0.0
  if i < n && chars[i] == '.' {
    i = i + 1
    let mut scale = 0.1
    while i < n && chars[i] >= '0' && chars[i] <= '9' {
      frac_part = frac_part + (chars[i].to_int() - '0'.to_int()).to_double() * scale
      scale = scale * 0.1
      i = i + 1
    }
  }
  let result = int_part + frac_part
  if neg {
    -result
  } else {
    result
  }
}

///|
fn safe_parse_int(s : String) -> Int {
  let chars = s.to_array()
  let n = chars.length()
  if n == 0 {
    return 0
  }
  let mut i = 0
  let mut neg = false
  if chars[0] == '-' {
    neg = true
    i = 1
  } else if chars[0] == '+' {
    i = 1
  }
  let mut result = 0
  while i < n && chars[i] >= '0' && chars[i] <= '9' {
    result = result * 10 + (chars[i].to_int() - '0'.to_int())
    i = i + 1
  }
  if neg {
    -result
  } else {
    result
  }
}

///|
fn is_logical_op(w : String) -> Bool {
  match w {
    "and" | "or" | "not" => true
    _ => false
  }
}
