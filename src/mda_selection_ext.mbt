///|
/// Extended MDAnalysis-style atom selection language for `Universe` /
/// `AtomGroup` (`mda_universe.mbt`).
///
/// In addition to the basic keywords of `mda_selection.mbt`
/// (protein / water / backbone / name / resname / chain / element / resid /
/// all / none, combined with and / or / not), this evaluator supports:
/// - `segid X` / `type X` — segment identifier / force-field atom type
/// - `bynum 1:10` / `bynum 1 3 5` — 1-based atom numbers (ranges allowed)
/// - `byres <sel>` — expand the sub-selection to whole residues
/// - `around <cutoff> <sel>` — atoms within cutoff distance of the
///   sub-selection (needs `positions`; orthogonal minimum-image PBC when
///   `box` is given); excludes the sub-selection atoms themselves
/// - `point x y z <cutoff>` — atoms within cutoff of a fixed point
/// - `prop mass <op> <v>` / `prop charge <op> <v>` — numeric property
///   comparison with op in {>, <, >=, <=, ==}
///
/// Boolean combinations use standard MDAnalysis precedence:
/// `not` > `and` > `or`, with parentheses for grouping, e.g.
/// `(name N or name O) and not water`.

///|
/// Evaluate an extended selection expression; returns selected atom indices
/// (0-based, ascending) into `atoms`.
pub fn mda_select_universe(
  expr : String,
  atoms : Array[UnivAtom],
  positions~ : Array[Array[Double]]? = None,
  box~ : Array[Double]? = None,
) -> Array[Int32] {
  let tokens = mda_sel_ext_tokens(expr)
  let n = atoms.length()
  if tokens.length() == 0 || n == 0 {
    return []
  }
  let (mask, _) = mda_sel_ext_eval_expr(tokens, 0, atoms, positions, box)
  let selected : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if mask[i] {
      selected.push(i)
    }
  }
  selected
}

///|
/// Tokenize an expression, keeping `and` / `or` / `not` / `(` / `)` as
/// control tokens; all words are lowercased. Parentheses attached to words
/// (e.g. "(protein" or "backbone)") are split off into separate tokens.
fn mda_sel_ext_tokens(expr : String) -> Array[String] {
  let tokens : Array[String] = Array(capacity=0)
  let lower = expr.to_lower()
  for w in lower.split(" ") {
    let word = w.to_owned().trim().to_owned()
    if word.length() > 0 {
      mda_sel_ext_push_word(tokens, word)
    }
  }
  tokens
}

///|
fn mda_sel_ext_push_word(tokens : Array[String], word : String) -> Unit {
  let open = 40 // '('
  let close = 41 // ')'
  let n = word.length()
  let mut lead = 0
  while lead < n && word.unsafe_get(lead).to_int() == open {
    tokens.push("(")
    lead = lead + 1
  }
  let mut end = n
  while end > lead && word.unsafe_get(end - 1).to_int() == close {
    end = end - 1
  }
  if end > lead {
    let sb = StringBuilder()
    let mut i = lead
    while i < end {
      sb.write_char(word.unsafe_get(i).unsafe_to_char())
      i = i + 1
    }
    tokens.push(sb.to_string())
  }
  let mut tail = n
  while tail > end {
    tokens.push(")")
    tail = tail - 1
  }
}

///|
/// expr := term ("or" term)*
fn mda_sel_ext_eval_expr(
  tokens : Array[String],
  start : Int32,
  atoms : Array[UnivAtom],
  positions : Array[Array[Double]]?,
  box : Array[Double]?,
) -> (Array[Bool], Int32) {
  let r0 = mda_sel_ext_eval_term(tokens, start, atoms, positions, box)
  let mask = r0.0
  let mut p = r0.1
  while p < tokens.length() && tokens[p] == "or" {
    let r = mda_sel_ext_eval_term(tokens, p + 1, atoms, positions, box)
    let m2 = r.0
    let p2 = r.1
    for i = 0; i < mask.length(); i = i + 1 {
      mask[i] = mask[i] || m2[i]
    }
    p = p2
  }
  (mask, p)
}

///|
/// term := factor ("and" factor)*
fn mda_sel_ext_eval_term(
  tokens : Array[String],
  start : Int32,
  atoms : Array[UnivAtom],
  positions : Array[Array[Double]]?,
  box : Array[Double]?,
) -> (Array[Bool], Int32) {
  let r0 = mda_sel_ext_eval_factor(tokens, start, atoms, positions, box)
  let mask = r0.0
  let mut p = r0.1
  while p < tokens.length() && tokens[p] == "and" {
    let r = mda_sel_ext_eval_factor(tokens, p + 1, atoms, positions, box)
    let m2 = r.0
    let p2 = r.1
    for i = 0; i < mask.length(); i = i + 1 {
      mask[i] = mask[i] && m2[i]
    }
    p = p2
  }
  (mask, p)
}

///|
/// factor := "not" factor | "(" expr ")" | keyword-section
fn mda_sel_ext_eval_factor(
  tokens : Array[String],
  start : Int32,
  atoms : Array[UnivAtom],
  positions : Array[Array[Double]]?,
  box : Array[Double]?,
) -> (Array[Bool], Int32) {
  let n = atoms.length()
  if start >= tokens.length() {
    let m : Array[Bool] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      m.push(false)
    }
    return (m, start)
  }
  if tokens[start] == "not" {
    let (m, p) = mda_sel_ext_eval_factor(
      tokens,
      start + 1,
      atoms,
      positions,
      box,
    )
    for i = 0; i < m.length(); i = i + 1 {
      m[i] = !m[i]
    }
    return (m, p)
  }
  if tokens[start] == "(" {
    let (m, p) = mda_sel_ext_eval_expr(tokens, start + 1, atoms, positions, box)
    let pn = if p < tokens.length() && tokens[p] == ")" { p + 1 } else { p }
    return (m, pn)
  }
  // keyword section: consume tokens up to the next operator / parenthesis
  let mut end = start
  while end < tokens.length() &&
        tokens[end] != "and" &&
        tokens[end] != "or" &&
        tokens[end] != "not" &&
        tokens[end] != ")" {
    end = end + 1
  }
  let seg : Array[String] = Array(capacity=0)
  for i = start; i < end; i = i + 1 {
    seg.push(tokens[i])
  }
  let mask = mda_sel_ext_mask(seg, atoms, positions, box)
  (mask, end)
}

///|
/// Slice `words[from..]` (used by keyword sections with sub-selections).
fn mda_sel_ext_tail(words : Array[String], from : Int32) -> Array[String] {
  let tail : Array[String] = Array(capacity=0)
  for i = from; i < words.length(); i = i + 1 {
    tail.push(words[i])
  }
  tail
}

///|
fn mda_sel_ext_mask(
  words : Array[String],
  atoms : Array[UnivAtom],
  positions : Array[Array[Double]]?,
  box : Array[Double]?,
) -> Array[Bool] {
  let n = atoms.length()
  let mask : Array[Bool] = Array(capacity=n)
  if words.length() == 0 {
    for i = 0; i < n; i = i + 1 {
      mask.push(false)
    }
    return mask
  }
  let kw = words[0]
  if kw == "protein" {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].is_protein)
    }
  } else if kw == "water" {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].is_water)
    }
  } else if kw == "backbone" {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].is_backbone)
    }
  } else if kw == "all" {
    for i = 0; i < n; i = i + 1 {
      mask.push(true)
    }
  } else if kw == "none" {
    for i = 0; i < n; i = i + 1 {
      mask.push(false)
    }
  } else if kw == "name" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].name.to_lower() == words[1])
    }
  } else if kw == "resname" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].resname.to_lower() == words[1])
    }
  } else if kw == "chain" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].chain.to_lower() == words[1])
    }
  } else if kw == "segid" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].segid.to_lower() == words[1])
    }
  } else if kw == "element" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].element.to_lower() == words[1])
    }
  } else if kw == "type" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      mask.push(atoms[i].atom_type.to_lower() == words[1])
    }
  } else if kw == "resid" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      mask.push(mda_sel_ext_resid_match(words[1], atoms[i].resid))
    }
  } else if kw == "bynum" && words.length() > 1 {
    for i = 0; i < n; i = i + 1 {
      let mut hit = false
      for k = 1; k < words.length(); k = k + 1 {
        if mda_sel_ext_resid_match(words[k], i + 1) {
          hit = true
        }
      }
      mask.push(hit)
    }
  } else if kw == "byres" && words.length() > 1 {
    let seeds = mda_sel_ext_mask(
      mda_sel_ext_tail(words, 1),
      atoms,
      positions,
      box,
    )
    for i = 0; i < n; i = i + 1 {
      mask.push(false)
    }
    for i = 0; i < n; i = i + 1 {
      if seeds[i] {
        let key = mda_sel_ext_reskey(atoms[i])
        for j = 0; j < n; j = j + 1 {
          if mda_sel_ext_reskey(atoms[j]) == key {
            mask[j] = true
          }
        }
      }
    }
  } else if kw == "around" && words.length() > 2 {
    let cut = match parse_double(words[1]) {
      Some(v) => v
      None => 0.0
    }
    let seeds = mda_sel_ext_mask(
      mda_sel_ext_tail(words, 2),
      atoms,
      positions,
      box,
    )
    for i = 0; i < n; i = i + 1 {
      mask.push(false)
    }
    match positions {
      Some(pos) =>
        if pos.length() >= n {
          for i = 0; i < n; i = i + 1 {
            if !seeds[i] {
              let mut hit = false
              for j = 0; j < n; j = j + 1 {
                if seeds[j] && mda_sel_ext_dist(pos[i], pos[j], box) <= cut {
                  hit = true
                }
              }
              mask[i] = hit
            }
          }
        }
      None => ()
    }
  } else if kw == "point" && words.length() > 4 {
    let px = match parse_double(words[1]) {
      Some(v) => v
      None => 0.0
    }
    let py = match parse_double(words[2]) {
      Some(v) => v
      None => 0.0
    }
    let pz = match parse_double(words[3]) {
      Some(v) => v
      None => 0.0
    }
    let cut = match parse_double(words[4]) {
      Some(v) => v
      None => 0.0
    }
    match positions {
      Some(pos) => {
        let ref_pt = [px, py, pz]
        for i = 0; i < n; i = i + 1 {
          if i < pos.length() {
            mask.push(mda_sel_ext_dist(pos[i], ref_pt, box) <= cut)
          } else {
            mask.push(false)
          }
        }
      }
      None =>
        for i = 0; i < n; i = i + 1 {
          mask.push(false)
        }
    }
  } else if kw == "prop" && words.length() > 3 {
    let field = words[1]
    let op = words[2]
    let val = match parse_double(words[3]) {
      Some(v) => v
      None => 0.0
    }
    for i = 0; i < n; i = i + 1 {
      let x = if field == "mass" {
        atoms[i].mass
      } else if field == "charge" {
        atoms[i].charge
      } else if field == "resid" {
        atoms[i].resid.to_double()
      } else {
        0.0
      }
      mask.push(mda_sel_ext_compare(x, op, val))
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      mask.push(false)
    }
  }
  mask
}

///|
fn mda_sel_ext_reskey(a : UnivAtom) -> String {
  a.chain + "|" + a.segid + "|" + a.resid.to_string()
}

///|
fn mda_sel_ext_resid_match(word : String, resid : Int32) -> Bool {
  if word.contains(":") {
    let parts = word.split(":").to_array()
    if parts.length() == 2 {
      let lo = parse_int(parts[0].to_owned())
      let hi = parse_int(parts[1].to_owned())
      return resid >= lo && resid <= hi
    }
    return false
  }
  resid == parse_int(word)
}

///|
fn mda_sel_ext_compare(x : Double, op : String, v : Double) -> Bool {
  if op == ">" {
    x > v
  } else if op == "<" {
    x < v
  } else if op == ">=" {
    x >= v
  } else if op == "<=" {
    x <= v
  } else if op == "==" {
    x == v
  } else {
    false
  }
}

///|
fn mda_sel_ext_dist(
  a : Array[Double],
  b : Array[Double],
  box : Array[Double]?,
) -> Double {
  let mut dx = a[0] - b[0]
  let mut dy = a[1] - b[1]
  let mut dz = a[2] - b[2]
  match box {
    Some(bx) =>
      if bx.length() >= 3 {
        dx = mda_sel_ext_min_image(dx, bx[0])
        dy = mda_sel_ext_min_image(dy, bx[1])
        dz = mda_sel_ext_min_image(dz, bx[2])
      }
    None => ()
  }
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
fn mda_sel_ext_min_image(d : Double, l : Double) -> Double {
  if l <= 0.0 {
    return d
  }
  let k = (d / l + 0.5).floor()
  d - k * l
}
