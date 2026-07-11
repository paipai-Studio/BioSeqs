///|
/// Bio.Phylo - Phylogenetic tree data structures.

///|
/// A clade (node) in a phylogenetic tree.
/// Contains branch length, name, confidence, child clades, and parent reference.
pub struct Clade {
  branch_length : Double?
  name : String?
  confidence : Double?
  clades : Array[Clade]
  parent : Clade?
}

///|
/// Construct a Clade with optional attributes.
pub fn Clade::new(
  branch_length? : Double? = None,
  name? : String? = None,
  confidence? : Double? = None,
  clades? : Array[Clade] = [],
  parent? : Clade? = None,
) -> Clade {
  let mut result = Clade::{
    branch_length,
    name,
    confidence,
    clades: [],
    parent,
  }
  for clade in clades {
    result = result.add_child(clade)
  }
  result
}

///|
/// Check if this clade is the root of the tree.
pub fn Clade::is_root(self : Clade) -> Bool {
  self.parent is None
}

///|
/// Add a child clade to this clade and return new Clade.
pub fn Clade::add_child(self : Clade, child : Clade) -> Clade {
  let child_with_parent = Clade::{
    branch_length: child.branch_length,
    name: child.name,
    confidence: child.confidence,
    clades: child.clades,
    parent: Some(self),
  }
  let new_clades = self.clades
  new_clades.push(child_with_parent)
  Clade::{
    branch_length: self.branch_length,
    name: self.name,
    confidence: self.confidence,
    clades: new_clades,
    parent: self.parent,
  }
}

///|
/// Check if this clade is a terminal (leaf) node.
pub fn Clade::is_terminal(self : Clade) -> Bool {
  self.clades.length() == 0
}

///|
/// Return the number of terminal (leaf) nodes in this clade's subtree.
pub fn Clade::count_terminals(self : Clade) -> Int {
  if self.is_terminal() {
    1
  } else {
    let mut count = 0
    for clade in self.clades {
      count = count + clade.count_terminals()
    }
    count
  }
}

///|
/// Get all terminal (leaf) clades in this subtree.
pub fn Clade::get_terminals(self : Clade) -> Array[Clade] {
  let result = Array::new()
  self.collect_terminals_into(result)
  result
}

///|
/// Helper to collect terminals into an existing array.
fn Clade::collect_terminals_into(self : Clade, result : Array[Clade]) -> Unit {
  if self.is_terminal() {
    result.push(self)
  } else {
    for clade in self.clades {
      clade.collect_terminals_into(result)
    }
  }
}

///|
/// Get all non-terminal (internal) clades in this subtree.
pub fn Clade::get_nonterminals(self : Clade) -> Array[Clade] {
  let result = Array::new()
  self.collect_nonterminals_into(result)
  result
}

///|
/// Helper to collect non-terminals into an existing array.
fn Clade::collect_nonterminals_into(
  self : Clade,
  result : Array[Clade],
) -> Unit {
  if !self.is_terminal() {
    result.push(self)
    for clade in self.clades {
      clade.collect_nonterminals_into(result)
    }
  }
}

///|
/// Find all clades matching given attributes.
pub fn Clade::find_clades(
  self : Clade,
  name? : String? = None,
  terminal? : Bool? = None,
) -> Array[Clade] {
  let result = Array::new()
  self.find_clades_into(name, terminal, result)
  result
}

///|
/// Helper to find clades and collect into existing array.
fn Clade::find_clades_into(
  self : Clade,
  name : String?,
  terminal : Bool?,
  result : Array[Clade],
) -> Unit {
  let match_name = if name is Some(_) { self.name == name } else { true }
  let match_terminal = if terminal is Some(_) {
    self.is_terminal() == terminal.unwrap()
  } else {
    true
  }

  if match_name && match_terminal {
    result.push(self)
  }

  for clade in self.clades {
    clade.find_clades_into(name, terminal, result)
  }
}

///|
/// Find the first clade matching given attributes, or None.
pub fn Clade::find_any(
  self : Clade,
  name? : String? = None,
  terminal? : Bool? = None,
) -> Clade? {
  let match_name = if name is Some(_) { self.name == name } else { true }
  let match_terminal = if terminal is Some(_) {
    self.is_terminal() == terminal.unwrap()
  } else {
    true
  }

  if match_name && match_terminal {
    return Some(self)
  }

  for clade in self.clades {
    match clade.find_any(name~, terminal~) {
      Some(found) => return Some(found)
      None => ()
    }
  }

  None
}

///|
/// Check if the tree downstream of this node is strictly bifurcating.
pub fn Clade::is_bifurcating(self : Clade) -> Bool {
  let n = self.clades.length()
  if n != 0 && n != 2 {
    return false
  }
  for clade in self.clades {
    if !clade.is_bifurcating() {
      return false
    }
  }
  true
}

///|
/// Find a path from this clade to a descendant with the given name.
pub fn Clade::get_path(self : Clade, name : String) -> Array[Clade]? {
  if self.name == Some(name) {
    let path = Array::new()
    path.push(self)
    Some(path)
  } else {
    for clade in self.clades {
      match clade.get_path(name) {
        Some(child_path) => {
          let path = Array::new()
          path.push(self)
          for c in child_path {
            path.push(c)
          }
          return Some(path)
        }
        None => ()
      }
    }
    None
  }
}

///|
/// Calculate depths of all clades from this clade.
pub fn Clade::depths(
  self : Clade,
  current_depth? : Double? = None,
) -> Array[(String, Double)] {
  let result = Array::new()
  self.depths_helper(current_depth, result)
  result
}

///|
/// Calculate depths and return as a Map for faster lookup.
pub fn Clade::depths_map(
  self : Clade,
  current_depth? : Double? = None,
) -> Map[String, Double] {
  let result : Map[String, Double] = Map([], capacity=16)
  self.depths_map_helper(current_depth, result)
  result
}

///|
/// Helper to collect depths into an existing array.
fn Clade::depths_helper(
  self : Clade,
  current_depth : Double?,
  result : Array[(String, Double)],
) -> Unit {
  let depth = match current_depth {
    Some(d) => d + self.branch_length.unwrap_or(0.0)
    None => self.branch_length.unwrap_or(0.0)
  }

  match self.name {
    Some(n) => result.push((n, depth))
    None => ()
  }

  for clade in self.clades {
    clade.depths_helper(Some(depth), result)
  }
}

///|
/// Helper to collect depths into an existing Map.
fn Clade::depths_map_helper(
  self : Clade,
  current_depth : Double?,
  result : Map[String, Double],
) -> Unit {
  let depth = match current_depth {
    Some(d) => d + self.branch_length.unwrap_or(0.0)
    None => self.branch_length.unwrap_or(0.0)
  }

  match self.name {
    Some(n) => result[n] = depth
    None => ()
  }

  for clade in self.clades {
    clade.depths_map_helper(Some(depth), result)
  }
}

///|
/// A phylogenetic tree.
pub struct Tree {
  root : Clade
  rooted : Bool
  name : String?
  weight : Double
}

///|
/// Construct a Tree with root clade.
pub fn Tree::new(
  root : Clade,
  rooted? : Bool = false,
  name? : String? = None,
  weight? : Double = 1.0,
) -> Tree {
  Tree::{ root, rooted, name, weight }
}

///|
/// Check if this tree is rooted.
pub fn Tree::is_rooted(self : Tree) -> Bool {
  self.rooted
}

///|
/// Get the root clade.
pub fn Tree::get_root(self : Tree) -> Clade {
  self.root
}

///|
/// Return the number of terminal (leaf) nodes in this tree.
pub fn Tree::count_terminals(self : Tree) -> Int {
  self.root.count_terminals()
}

///|
/// Get all terminal (leaf) clades in this tree.
pub fn Tree::get_terminals(self : Tree) -> Array[Clade] {
  self.root.get_terminals()
}

///|
/// Get all non-terminal (internal) clades in this tree.
pub fn Tree::get_nonterminals(self : Tree) -> Array[Clade] {
  self.root.get_nonterminals()
}

///|
/// Find all clades matching given attributes.
pub fn Tree::find_clades(
  self : Tree,
  name? : String? = None,
  terminal? : Bool? = None,
) -> Array[Clade] {
  self.root.find_clades(name~, terminal~)
}

///|
/// Find the first clade matching given attributes, or None.
pub fn Tree::find_any(
  self : Tree,
  name? : String? = None,
  terminal? : Bool? = None,
) -> Clade? {
  self.root.find_any(name~, terminal~)
}

///|
/// Check if the tree is strictly bifurcating.
pub fn Tree::is_bifurcating(self : Tree) -> Bool {
  self.root.is_bifurcating()
}

///|
/// Find a path from root to a clade with the given name.
pub fn Tree::get_path(self : Tree, name : String) -> Array[Clade]? {
  self.root.get_path(name)
}

///|
/// Calculate depths of all clades from root.
pub fn Tree::depths(self : Tree) -> Array[(String, Double)] {
  self.root.depths()
}

///|
/// Calculate the distance between two clades by name.
pub fn Tree::distance(self : Tree, name1 : String, name2 : String) -> Double? {
  match (self.root.get_path(name1), self.root.get_path(name2)) {
    (Some(path1), Some(path2)) => {
      let mut i = 0
      let min_len = if path1.length() < path2.length() {
        path1.length()
      } else {
        path2.length()
      }
      while i < min_len && path1[i].name == path2[i].name {
        i = i + 1
      }

      let mut distance = 0.0
      for j in i..<path1.length() {
        distance = distance + path1[j].branch_length.unwrap_or(0.0)
      }
      for j in i..<path2.length() {
        distance = distance + path2[j].branch_length.unwrap_or(0.0)
      }

      Some(distance)
    }
    _ => None
  }
}

///|
/// Find the most recent common ancestor of clades with given names.
pub fn Tree::common_ancestor(self : Tree, names : Array[String]) -> Clade? {
  if names.length() == 0 {
    return None
  }

  let paths : Array[Array[Clade]] = Array::new()
  for name in names {
    match self.root.get_path(name) {
      Some(path) => paths.push(path)
      None => return None
    }
  }

  if paths.length() == 0 {
    return None
  }

  let mut min_len = paths[0].length()
  for path in paths {
    if path.length() < min_len {
      min_len = path.length()
    }
  }

  let mut i = 0
  while i < min_len {
    let current = paths[0][i]
    let mut all_match = true
    for j in 1..<paths.length() {
      if paths[j][i].name != current.name {
        all_match = false
        break
      }
    }
    if !all_match {
      break
    }
    i = i + 1
  }

  if i > 0 {
    Some(paths[0][i - 1])
  } else {
    None
  }
}

///|
/// Generate ASCII art representation of the tree.
pub fn Tree::draw_ascii(self : Tree) -> String {
  let buf = StringBuilder::new()
  Clade::draw_ascii_helper(self.root, "", true, buf)
  buf.to_string()
}

///|
fn Clade::draw_ascii_helper(
  self : Clade,
  prefix : String,
  is_last : Bool,
  buf : StringBuilder,
) -> Unit {
  let name = match self.name {
    Some(n) => n
    None => ""
  }

  buf.write_string(prefix)
  if is_last {
    buf.write_string("└── ")
  } else {
    buf.write_string("├── ")
  }
  buf.write_string(name)

  buf.write_string("\n")

  let child_count = self.clades.length()
  for i in 0..<child_count {
    let child = self.clades[i]
    let child_is_last = i == child_count - 1
    let child_prefix = if is_last { prefix + "    " } else { prefix + "│   " }
    child.draw_ascii_helper(child_prefix, child_is_last, buf)
  }
}

///|
/// Return a string representation of the tree.
pub fn Tree::repr(self : Tree) -> String {
  let buf = StringBuilder::new()
  buf.write_string("Tree(rooted=")
  buf.write_string(if self.rooted { "true" } else { "false" })
  buf.write_string(", weight=")
  buf.write_string(self.weight.to_string())
  buf.write_string(")")
  buf.to_string()
}
