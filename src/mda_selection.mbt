/// mda_selection.mbt — Atom selection language parser and evaluator
/// Item 21 (selection language: "protein and not resname ALA", "name CA", "resid 1:50", etc.)

///|
/// Selection token types
pub enum MdaSelToken {
  MdaSelWord(String) // keyword or name
  MdaSelNumber(Double)
  MdaSelOp(String) // and, or, not
  MdaSelLParen
  MdaSelRParen
  MdaSelColon // range operator
}

///|
/// Atom info for selection evaluation
pub struct MdaAtomInfo {
  name : String // atom name: CA, N, O, etc.
  resname : String // residue name: ALA, GLY, etc.
  resid : Int32 // residue id
  chain : String // chain identifier
  element : String // element symbol
  is_protein : Bool // is protein atom
  is_water : Bool // is water
  is_backbone : Bool // is backbone atom
  x : Double
  y : Double
  z : Double
}

///|
/// Parse selection expression into tokens
pub fn mda_sel_tokenize(expr : String) -> Array[MdaSelToken] {
  let tokens : Array[MdaSelToken] = Array(capacity=0)
  let words = expr.split(" ")
  for w in words {
    let word = w.to_owned().trim().to_owned()
    if word.length() == 0 {
      continue
    }
    if word == "and" || word == "or" || word == "not" {
      tokens.push(MdaSelOp(word))
    } else if word == "(" {
      tokens.push(MdaSelLParen)
    } else if word == ")" {
      tokens.push(MdaSelRParen)
    } else if word == ":" {
      tokens.push(MdaSelColon)
    } else {
      let parsed = parse_int(word)
      if parsed != 0 {
        tokens.push(MdaSelNumber(parsed.to_double()))
      } else {
        tokens.push(MdaSelWord(word))
      }
    }
  }
  tokens
}

///|
/// Evaluate a simple selection expression against a list of atoms.
pub fn mda_sel_evaluate(
  expr : String,
  atoms : Array[MdaAtomInfo],
) -> Array[Int32] {
  let expr_lower = expr.to_lower()
  let words = expr_lower.split(" ")
  let word_list : Array[String] = Array(capacity=0)
  for w in words {
    let word = w.to_owned().trim().to_owned()
    if word.length() > 0 {
      word_list.push(word)
    }
  }
  let selected : Array[Int32] = Array(capacity=0)
  let n_atoms = atoms.length()
  let has_and = mda_sel_contains(word_list, "and")
  let has_or = mda_sel_contains(word_list, "or")
  let has_not = mda_sel_contains(word_list, "not")
  if has_and || has_or || has_not {
    let parts = mda_sel_split_ops(word_list)
    let part_masks : Array[Array[Bool]] = Array(capacity=parts.length())
    for p = 0; p < parts.length(); p = p + 1 {
      let mask : Array[Bool] = Array(capacity=n_atoms)
      for i = 0; i < n_atoms; i = i + 1 {
        mask.push(mda_sel_match_atom(parts[p], atoms[i]))
      }
      part_masks.push(mask)
    }
    let result_mask : Array[Bool] = Array(capacity=n_atoms)
    for i = 0; i < n_atoms; i = i + 1 {
      result_mask.push(part_masks[0][i])
    }
    for p = 1; p < parts.length(); p = p + 1 {
      for i = 0; i < n_atoms; i = i + 1 {
        if has_and && !has_or {
          result_mask[i] = result_mask[i] && part_masks[p][i]
        } else if has_or {
          result_mask[i] = result_mask[i] || part_masks[p][i]
        }
      }
    }
    if has_not && !has_and && !has_or {
      for i = 0; i < n_atoms; i = i + 1 {
        result_mask[i] = !result_mask[i]
      }
    }
    for i = 0; i < n_atoms; i = i + 1 {
      if result_mask[i] {
        selected.push(i)
      }
    }
    return selected
  }
  for i = 0; i < n_atoms; i = i + 1 {
    if mda_sel_match_atom(word_list, atoms[i]) {
      selected.push(i)
    }
  }
  selected
}

///|
fn mda_sel_contains(arr : Array[String], val : String) -> Bool {
  for i = 0; i < arr.length(); i = i + 1 {
    if arr[i] == val {
      return true
    }
  }
  false
}

///|
fn mda_sel_split_ops(words : Array[String]) -> Array[Array[String]] {
  let parts : Array[Array[String]] = Array(capacity=0)
  let current : Array[String] = Array(capacity=0)
  for i = 0; i < words.length(); i = i + 1 {
    let w = words[i]
    if w == "and" || w == "or" || w == "not" {
      if current.length() > 0 {
        parts.push(current.copy())
        current.clear()
      }
    } else {
      current.push(w)
    }
  }
  if current.length() > 0 {
    parts.push(current.copy())
  }
  parts
}

///|
fn mda_sel_match_atom(words : Array[String], atom : MdaAtomInfo) -> Bool {
  if words.length() == 0 {
    return false
  }
  let kw = words[0]
  if kw == "protein" {
    return atom.is_protein
  }
  if kw == "water" {
    return atom.is_water
  }
  if kw == "backbone" {
    return atom.is_backbone
  }
  if kw == "name" && words.length() > 1 {
    return atom.name.to_lower() == words[1]
  }
  if kw == "resname" && words.length() > 1 {
    return atom.resname.to_lower() == words[1]
  }
  if kw == "chain" && words.length() > 1 {
    return atom.chain.to_lower() == words[1]
  }
  if kw == "element" && words.length() > 1 {
    return atom.element.to_lower() == words[1]
  }
  if kw == "resid" && words.length() > 1 {
    let range_str = words[1]
    if mda_sel_str_contains(range_str, ":") {
      let parts = range_str.split(":")
      let part_arr = parts.to_array()
      if part_arr.length() == 2 {
        let lo = parse_int(part_arr[0].to_owned())
        let hi = parse_int(part_arr[1].to_owned())
        return atom.resid >= lo && atom.resid <= hi
      }
      return false
    }
    let parsed = parse_int(range_str)
    return atom.resid == parsed
  }
  if kw == "all" {
    return true
  }
  if kw == "none" {
    return false
  }
  false
}

///|
fn mda_sel_str_contains(s : String, ch : String) -> Bool {
  s.contains(ch)
}

///|
pub fn mda_selection_dummy() -> Unit {
  ignore(MdaSelWord("test"))
  ignore(MdaSelNumber(0.0))
  ignore(MdaSelOp("and"))
  ignore(MdaSelLParen)
  ignore(MdaSelRParen)
  ignore(MdaSelColon)
  ignore(MdaAtomInfo::{
    name: "",
    resname: "",
    resid: 0,
    chain: "",
    element: "",
    is_protein: false,
    is_water: false,
    is_backbone: false,
    x: 0.0,
    y: 0.0,
    z: 0.0,
  })
}
