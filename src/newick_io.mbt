///|
/// Newick format parser and writer for phylogenetic trees.

///|
/// Parse a Newick string into a Tree.
pub fn parse_newick(s : String) -> Tree {
  let root = parse_newick_helper(s, 0).0
  Tree::new(root)
}

///|
/// Write a Tree to Newick format string.
pub fn write_newick(tree : Tree) -> String {
  let buf = StringBuilder::new()
  write_clade_to_buffer(tree.root, buf)
  buf.write_string(";")
  buf.to_string()
}

///|
fn parse_newick_helper(s : String, pos : Int) -> (Clade, Int) {
  let mut p = skip_whitespace(s, pos)

  if p < s.length() && s.get_char(p) == Some('(') {
    p = p + 1
    let clades = Array::new()

    let mut pos = skip_whitespace(s, p)
    while pos < s.length() && s.get_char(pos) != Some(')') {
      if clades.length() > 0 && s.get_char(pos) == Some(',') {
        pos = pos + 1
      }
      let (clade, new_pos) = parse_newick_helper(s, pos)
      clades.push(clade)
      pos = skip_whitespace(s, new_pos)
    }

    p = pos + 1
    let (name, new_pos1) = parse_name(s, p)
    p = new_pos1
    let (branch_length, new_pos2) = parse_branch_length(s, p)
    p = new_pos2

    (
      Clade::new(branch_length~, name~, confidence=None, clades~, parent=None),
      p,
    )
  } else {
    let (name, new_pos1) = parse_name(s, p)
    p = new_pos1
    let (branch_length, new_pos2) = parse_branch_length(s, p)
    p = new_pos2
    let empty_clades = Array::new()
    (
      Clade::new(
        branch_length~,
        name~,
        confidence=None,
        clades=empty_clades,
        parent=None,
      ),
      p,
    )
  }
}

///|
fn skip_whitespace(s : String, pos : Int) -> Int {
  let mut p = pos
  while p < s.length() {
    match s.get_char(p) {
      Some(c) => {
        let code = c.to_int()
        if code == 32 || code == 9 || code == 10 || code == 13 {
          p = p + 1
        } else {
          break
        }
      }
      None => break
    }
  }
  p
}

///|
fn parse_name(s : String, pos : Int) -> (String?, Int) {
  let mut p = skip_whitespace(s, pos)

  if p >= s.length() {
    return (None, p)
  }

  match s.get_char(p) {
    Some(c) => {
      let code = c.to_int()
      if code == 40 || code == 41 || code == 44 || code == 58 || code == 59 {
        return (None, p)
      }

      let buf = StringBuilder::new()
      while p < s.length() {
        match s.get_char(p) {
          Some(ch) => {
            let ch_code = ch.to_int()
            if ch_code == 40 ||
              ch_code == 41 ||
              ch_code == 44 ||
              ch_code == 58 ||
              ch_code == 59 {
              break
            }
            buf.write_char(ch)
            p = p + 1
          }
          None => break
        }
      }

      let name = buf.to_string()
      if name.length() == 0 {
        (None, p)
      } else {
        (Some(name), p)
      }
    }
    None => (None, p)
  }
}

///|
fn parse_branch_length(s : String, pos : Int) -> (Double?, Int) {
  let mut p = skip_whitespace(s, pos)

  if p >= s.length() {
    return (None, p)
  }

  match s.get_char(p) {
    Some(':') => {
      p = p + 1
      p = skip_whitespace(s, p)

      let buf = StringBuilder::new()
      let mut has_dot = false
      while p < s.length() {
        match s.get_char(p) {
          Some(ch) => {
            let code = ch.to_int()
            if code >= 48 && code <= 57 {
              buf.write_char(ch)
              p = p + 1
            } else if code == 46 && !has_dot {
              buf.write_char(ch)
              has_dot = true
              p = p + 1
            } else if code == 45 && buf.to_string().length() == 0 {
              buf.write_char(ch)
              p = p + 1
            } else {
              break
            }
          }
          None => break
        }
      }

      let num_str = buf.to_string()
      if num_str.length() == 0 {
        (None, p)
      } else {
        (parse_double(num_str), p)
      }
    }
    _ => (None, p)
  }
}

///|
fn write_clade_to_buffer(clade : Clade, buf : StringBuilder) -> Unit {
  if !clade.is_terminal() {
    buf.write_char('(')
    let mut first = true
    for c in clade.clades {
      if first {
        first = false
      } else {
        buf.write_char(',')
      }
      write_clade_to_buffer(c, buf)
    }
    buf.write_char(')')
  }

  match clade.name {
    Some(n) => buf.write_string(n)
    None => ()
  }

  match clade.branch_length {
    Some(l) => {
      buf.write_char(':')
      buf.write_string(l.to_string())
    }
    None => ()
  }
}
