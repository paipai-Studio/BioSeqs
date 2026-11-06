/// biotite_dotbracket.mbt — Dot-bracket notation for RNA secondary structure
/// Inspired by biotite.sequence.dotbracket

///|
/// RNA secondary structure in dot-bracket notation.
/// Supports pseudoknots using multiple bracket levels: (), [], {}, <>.
pub struct BiotiteDotBracket {
  base_pairs : Array[(Int32, Int32)] // (i, j) pairs, i < j
  length : Int32 // sequence length
}

///|
/// Parse dot-bracket notation into base pairs.
/// Supports multi-bracket pseudoknot notation: ()[]{}<>
pub fn BiotiteDotBracket::from_string(dotbracket : String) -> BiotiteDotBracket {
  let n = dotbracket.length()
  let pairs : Array[(Int32, Int32)] = Array(capacity=0)
  // Stack for each bracket type
  let open_paren : Array[Int32] = Array(capacity=0)
  let open_bracket : Array[Int32] = Array(capacity=0)
  let open_brace : Array[Int32] = Array(capacity=0)
  let open_angle : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let ch = dotbracket[i]
    match ch {
      '(' => open_paren.push(i)
      ')' =>
        if open_paren.length() > 0 {
          let open_idx = open_paren[open_paren.length() - 1]
          let _ = open_paren.remove(open_paren.length() - 1)
          pairs.push((open_idx, i))
        }
      '[' => open_bracket.push(i)
      ']' =>
        if open_bracket.length() > 0 {
          let open_idx = open_bracket[open_bracket.length() - 1]
          let _ = open_bracket.remove(open_bracket.length() - 1)
          pairs.push((open_idx, i))
        }
      '{' => open_brace.push(i)
      '}' =>
        if open_brace.length() > 0 {
          let open_idx = open_brace[open_brace.length() - 1]
          let _ = open_brace.remove(open_brace.length() - 1)
          pairs.push((open_idx, i))
        }
      '<' => open_angle.push(i)
      '>' =>
        if open_angle.length() > 0 {
          let open_idx = open_angle[open_angle.length() - 1]
          let _ = open_angle.remove(open_angle.length() - 1)
          pairs.push((open_idx, i))
        }
      _ => ()
    }
  }
  // Sort pairs by first index
  pairs.sort_by(fn(a, b) {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })
  { base_pairs: pairs, length: n, }
}

///|
/// Convert base pairs to dot-bracket notation.
/// Uses multi-bracket notation for pseudoknots.
pub fn BiotiteDotBracket::to_string(self : BiotiteDotBracket) -> String {
  let chars : Array[Int32] = Array::make(self.length, 46) // '.' = 46
  let brackets = ["()", "[]", "{}", "<>"]
  // Assign bracket levels greedily
  // A pair (i,j) can use level L if no other pair at level L crosses it
  let levels : Array[Array[(Int32, Int32)]] = Array::make(4, [])
  for pair_idx = 0; pair_idx < self.base_pairs.length(); pair_idx = pair_idx + 1 {
    let (i, j) = self.base_pairs[pair_idx]
    let mut assigned = false
    for level = 0; level < 4; level = level + 1 {
      let mut can_use = true
      for k = 0; k < levels[level].length(); k = k + 1 {
        let (oi, oj) = levels[level][k]
        // Check for crossing: i < oi < j < oj or oi < i < oj < j
        if (i < oi && oi < j && j < oj) || (oi < i && i < oj && oj < j) {
          can_use = false
          break
        }
      }
      if can_use {
        levels[level].push((i, j))
        chars[i] = brackets[level][0].to_int()
        chars[j] = brackets[level][1].to_int()
        assigned = true
        break
      }
    }
    // If no level available, leave as dots (shouldn't happen with 4 levels for typical structures)
    ignore(assigned)
  }
  let sb = StringBuilder()
  for i = 0; i < self.length; i = i + 1 {
    let _ = sb.write_char(chars[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Get the number of base pairs.
pub fn BiotiteDotBracket::num_base_pairs(self : BiotiteDotBracket) -> Int32 {
  self.base_pairs.length()
}

///|
/// Check if the structure contains pseudoknots.
/// A pseudoknot exists if two pairs (i,j) and (k,l) cross: i < k < j < l.
pub fn BiotiteDotBracket::has_pseudoknot(self : BiotiteDotBracket) -> Bool {
  let n = self.base_pairs.length()
  for a = 0; a < n; a = a + 1 {
    let (i, j) = self.base_pairs[a]
    for b = a + 1; b < n; b = b + 1 {
      let (k, l) = self.base_pairs[b]
      // Check crossing
      if (i < k && k < j && j < l) || (k < i && i < l && l < j) {
        return true
      }
    }
  }
  false
}

///|
/// Get the partner of a position, if any.
pub fn BiotiteDotBracket::get_partner(
  self : BiotiteDotBracket,
  pos : Int32,
) -> Int32? {
  for i = 0; i < self.base_pairs.length(); i = i + 1 {
    let (a, b) = self.base_pairs[i]
    if a == pos {
      return Some(b)
    }
    if b == pos {
      return Some(a)
    }
  }
  None
}

///|
/// Check if a position is paired.
pub fn BiotiteDotBracket::is_paired(
  self : BiotiteDotBracket,
  pos : Int32,
) -> Bool {
  self.get_partner(pos) is Some(_)
}

///|
/// Get all unpaired positions.
pub fn BiotiteDotBracket::unpaired_positions(
  self : BiotiteDotBracket,
) -> Array[Int32] {
  let paired : Array[Bool] = Array::make(self.length, false)
  for i = 0; i < self.base_pairs.length(); i = i + 1 {
    let (a, b) = self.base_pairs[i]
    paired[a] = true
    paired[b] = true
  }
  let result : Array[Int32] = Array(capacity=0)
  for i = 0; i < self.length; i = i + 1 {
    if biotite_db_not(paired[i]) {
      result.push(i)
    }
  }
  result
}

///|
/// Helper: negate boolean
fn biotite_db_not(b : Bool) -> Bool {
  !b
}
