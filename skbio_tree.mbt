///|
/// scikit-bio Tree module - Phylogenetic tree data structures.

///|
/// A node within a phylogenetic tree.
/// Stores links to parent and child nodes, branch length, and support value.
pub struct TreeNode {
  name : String?
  length : Double?
  support : Double?
  parent : TreeNode?
  children : Array[TreeNode]
}

///|
/// Construct a TreeNode with optional attributes.
pub fn TreeNode::new(
  name? : String? = None,
  length? : Double? = None,
  support? : Double? = None,
  parent? : TreeNode? = None,
  children? : Array[TreeNode] = [],
) -> TreeNode {
  let mut result = TreeNode::{ name, length, support, parent, children: [] }
  for child in children {
    result = result.add_child(child)
  }
  result
}

///|
/// Check if this node is the root of the tree.
pub fn TreeNode::is_root(self : TreeNode) -> Bool {
  self.parent is None
}

///|
/// Check if this node is a tip (leaf) of the tree.
pub fn TreeNode::is_tip(self : TreeNode) -> Bool {
  self.children.length() == 0
}

///|
/// Check if this node has any children.
pub fn TreeNode::has_children(self : TreeNode) -> Bool {
  self.children.length() > 0
}

///|
/// Return the number of children this node has.
pub fn TreeNode::num_children(self : TreeNode) -> Int {
  self.children.length()
}

///|
/// Add a child node to this node and return new TreeNode.
pub fn TreeNode::add_child(self : TreeNode, child : TreeNode) -> TreeNode {
  let child_with_parent = TreeNode::{
    name: child.name,
    length: child.length,
    support: child.support,
    parent: Some(self),
    children: child.children,
  }
  let new_children = self.children
  new_children.push(child_with_parent)
  TreeNode::{
    name: self.name,
    length: self.length,
    support: self.support,
    parent: self.parent,
    children: new_children,
  }
}

///|
/// Remove a child node from this node and return new TreeNode.
pub fn TreeNode::remove_child(
  self : TreeNode,
  child : TreeNode,
) -> (TreeNode, Bool) {
  let mut removed = false
  let new_children = Array::new()
  for c in self.children {
    if !(c.name == child.name && c.length == child.length) {
      new_children.push(c)
    } else {
      removed = true
    }
  }
  (
    TreeNode::{
      name: self.name,
      length: self.length,
      support: self.support,
      parent: self.parent,
      children: new_children,
    },
    removed,
  )
}

///|
/// Return root of the tree which contains this node.
pub fn TreeNode::root(self : TreeNode) -> TreeNode {
  match self.parent {
    Some(p) => p.root()
    None => self
  }
}

///|
/// Return all nodes in a path from this node to another node.
pub fn TreeNode::path(self : TreeNode, other : TreeNode) -> Array[TreeNode]? {
  let self_root = self.root()
  let other_root = other.root()

  if !(self_root.name == other_root.name &&
    self_root.children.length() == other_root.children.length()) {
    return None
  }

  let ancestors_self = self.ancestors()
  let ancestors_other = other.ancestors()

  let path : Array[TreeNode] = Array::new()
  path.push(self)

  let mut found = false
  for node in ancestors_self {
    path.push(node)
    if node.name == other.name {
      found = true
      break
    }
  }

  if !found {
    let other_path : Array[TreeNode] = Array::new()
    for node in ancestors_other {
      other_path.push(node)
      if node.name == self.name {
        break
      }
    }
    let n = other_path.length()
    for i = 0; i < n; i = i + 1 {
      path.push(other_path[n - 1 - i])
    }
  }

  Some(path)
}

///|
/// Return all ancestral nodes from this node back to the root.
pub fn TreeNode::ancestors(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  let mut current = self.parent
  while current is Some(_) {
    let p = current.unwrap()
    result.push(p)
    current = p.parent
  }
  result
}

///|
/// Return all nodes that are siblings of this node.
pub fn TreeNode::siblings(self : TreeNode) -> Array[TreeNode] {
  match self.parent {
    None => Array::new()
    Some(p) => {
      let result = Array::new()
      for child in p.children {
        if !(child.name == self.name && child.length == self.length) {
          result.push(child)
        }
      }
      result
    }
  }
}

///|
/// Return all nodes that are neighbors of this node.
pub fn TreeNode::neighbors(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  match self.parent {
    Some(p) => result.push(p)
    None => ()
  }
  for child in self.children {
    result.push(child)
  }
  result
}

///|
/// Find the lowest common ancestor of a list of nodes.
pub fn TreeNode::lca(self : TreeNode, nodes : Array[TreeNode]) -> TreeNode? {
  if nodes.length() == 0 {
    return None
  }

  let paths : Array[Array[TreeNode]] = Array::new()
  for node in nodes {
    let path = node.ancestors()
    path.push(node)
    paths.push(path)
  }

  let mut min_len = paths[0].length()
  for path in paths {
    if path.length() < min_len {
      min_len = path.length()
    }
  }

  let mut i = 0
  while i < min_len {
    let current = paths[0][paths[0].length() - 1 - i]
    let mut all_match = true
    let m = paths.length()
    for j = 1; j < m; j = j + 1 {
      let other = paths[j][paths[j].length() - 1 - i]
      if !(current.name == other.name && current.length == other.length) {
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
    Some(paths[0][paths[0].length() - i])
  } else {
    None
  }
}

///|
/// Perform preorder traversal over tree starting from this node.
pub fn TreeNode::preorder(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  self.preorder_helper(result)
  result
}

///|
fn TreeNode::preorder_helper(self : TreeNode, result : Array[TreeNode]) -> Unit {
  result.push(self)
  for child in self.children {
    child.preorder_helper(result)
  }
}

///|
/// Perform postorder traversal over tree starting from this node.
pub fn TreeNode::postorder(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  self.postorder_helper(result)
  result
}

///|
fn TreeNode::postorder_helper(
  self : TreeNode,
  result : Array[TreeNode],
) -> Unit {
  for child in self.children {
    child.postorder_helper(result)
  }
  result.push(self)
}

///|
/// Perform level order traversal over tree starting from this node.
pub fn TreeNode::levelorder(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  let queue = Array::new()
  queue.push(self)

  while queue.length() > 0 {
    let node = queue.remove(0)
    result.push(node)
    for child in node.children {
      queue.push(child)
    }
  }

  result
}

///|
/// Return all tips (leaves) descended from this node.
pub fn TreeNode::tips(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  self.collect_tips(result)
  result
}

///|
fn TreeNode::collect_tips(self : TreeNode, result : Array[TreeNode]) -> Unit {
  if self.is_tip() {
    result.push(self)
  } else {
    for child in self.children {
      child.collect_tips(result)
    }
  }
}

///|
/// Return all non-tip nodes descended from this node.
pub fn TreeNode::non_tips(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  self.collect_non_tips(result)
  result
}

///|
fn TreeNode::collect_non_tips(
  self : TreeNode,
  result : Array[TreeNode],
) -> Unit {
  if !self.is_tip() {
    result.push(self)
    for child in self.children {
      child.collect_non_tips(result)
    }
  }
}

///|
/// Calculate the distance from this node to another node.
pub fn TreeNode::distance(self : TreeNode, other : TreeNode) -> Double? {
  match self.path(other) {
    None => None
    Some(path) => {
      let mut dist = 0.0
      let n = path.length() - 1
      for i = 0; i < n; i = i + 1 {
        dist = dist + path[i].length.unwrap_or(0.0)
      }
      Some(dist)
    }
  }
}

///|
/// Calculate the total length of the tree (sum of all branch lengths).
pub fn TreeNode::total_length(self : TreeNode) -> Double {
  let mut total = self.length.unwrap_or(0.0)
  for child in self.children {
    total = total + child.total_length()
  }
  total
}

///|
/// Find a node by name in the subtree rooted at this node.
pub fn TreeNode::find(self : TreeNode, name : String) -> TreeNode? {
  if self.name == Some(name) {
    Some(self)
  } else {
    for child in self.children {
      match child.find(name) {
        Some(found) => return Some(found)
        None => ()
      }
    }
    None
  }
}

///|
/// Find all nodes matching a name.
pub fn TreeNode::find_by_name(
  self : TreeNode,
  name : String,
) -> Array[TreeNode] {
  let result = Array::new()
  self.find_by_name_helper(name, result)
  result
}

///|
fn TreeNode::find_by_name_helper(
  self : TreeNode,
  name : String,
  result : Array[TreeNode],
) -> Unit {
  if self.name == Some(name) {
    result.push(self)
  }
  for child in self.children {
    child.find_by_name_helper(name, result)
  }
}

///|
/// Return the number of tips in the subtree.
pub fn TreeNode::count_tips(self : TreeNode) -> Int {
  if self.is_tip() {
    1
  } else {
    let mut count = 0
    for child in self.children {
      count = count + child.count_tips()
    }
    count
  }
}

///|
/// Return the number of nodes in the subtree.
pub fn TreeNode::count(self : TreeNode) -> Int {
  let mut count = 1
  for child in self.children {
    count = count + child.count()
  }
  count
}

///|
/// Newick format parser.
pub fn TreeNode::from_newick(newick : String) -> TreeNode raise ParseError {
  skbio_parse_newick(newick)
}

///|
/// Convert tree to Newick format string.
pub fn TreeNode::to_newick(self : TreeNode) -> String {
  let buf = StringBuilder::new()
  self.to_newick_helper(buf)
  buf.write_string(";")
  buf.to_string()
}

///|
fn TreeNode::to_newick_helper(self : TreeNode, buf : StringBuilder) -> Unit {
  if self.children.length() > 0 {
    buf.write_char('(')
    let mut first = true
    for child in self.children {
      if !first {
        buf.write_char(',')
      }
      first = false
      child.to_newick_helper(buf)
    }
    buf.write_char(')')
  }

  match self.name {
    Some(n) => buf.write_string(n)
    None => ()
  }

  match self.length {
    Some(l) => {
      buf.write_char(':')
      buf.write_string(l.to_string())
    }
    None => ()
  }
}

///|
/// ParseError for Newick parsing.
pub suberror ParseError {
  InvalidFormat(String)
  MissingSemicolon
  UnbalancedParentheses
  InvalidNumber(String)
}

///|
/// Parse Newick format string into TreeNode.
fn skbio_parse_newick(input : String) -> TreeNode raise ParseError {
  let trimmed = input.trim()

  if !trimmed.has_suffix(";") {
    raise ParseError::MissingSemicolon
  }

  let s = trimmed[0:trimmed.length() - 1].to_owned()
  let n = s.length()
  let (root, pos) = parse_node_rec(s, 0, n)

  if pos < n {
    raise ParseError::InvalidFormat("Extra characters after root".to_string())
  }

  root
}

///|
fn parse_node_rec(
  s : String,
  pos : Int,
  len : Int,
) -> (TreeNode, Int) raise ParseError {
  let mut pos = skip_ws(s, pos, len)

  let children : Array[TreeNode] = Array::new()
  let mut name : String? = None
  let mut length : Double? = None

  if pos < len && s.unsafe_get(pos).unsafe_to_char() == '(' {
    pos = pos + 1
    pos = skip_ws(s, pos, len)

    if pos < len && s.unsafe_get(pos).unsafe_to_char() != ')' {
      let (child, new_pos) = parse_node_rec(s, pos, len)
      pos = new_pos
      children.push(child)

      while pos < len && s.unsafe_get(pos).unsafe_to_char() == ',' {
        pos = pos + 1
        pos = skip_ws(s, pos, len)
        let (child, new_pos) = parse_node_rec(s, pos, len)
        pos = new_pos
        children.push(child)
        pos = skip_ws(s, pos, len)
      }
    }

    if pos >= len || s.unsafe_get(pos).unsafe_to_char() != ')' {
      raise ParseError::UnbalancedParentheses
    }
    pos = pos + 1
    pos = skip_ws(s, pos, len)
  }

  if pos < len {
    let c = s.unsafe_get(pos).unsafe_to_char()
    if c != ':' && c != ',' && c != ')' && c != ';' {
      let (label, new_pos) = parse_label_rec(s, pos, len)
      name = Some(label)
      pos = new_pos
      pos = skip_ws(s, pos, len)
    }
  }

  if pos < len && s.unsafe_get(pos).unsafe_to_char() == ':' {
    pos = pos + 1
    pos = skip_ws(s, pos, len)
    let (num, new_pos) = parse_number_rec(s, pos, len)
    length = Some(num)
    pos = new_pos
    pos = skip_ws(s, pos, len)
  }

  let node = TreeNode::new(name~, length~, children~)

  (node, pos)
}

///|
fn parse_label_rec(s : String, pos : Int, len : Int) -> (String, Int) {
  let buf = StringBuilder::new()
  let mut pos = pos

  if pos < len && s.unsafe_get(pos).unsafe_to_char() == '\'' {
    pos = pos + 1
    while pos < len && s.unsafe_get(pos).unsafe_to_char() != '\'' {
      let c = s.unsafe_get(pos).unsafe_to_char()
      if c == '\'' &&
        pos + 1 < len &&
        s.unsafe_get(pos + 1).unsafe_to_char() == '\'' {
        buf.write_char('\'')
        pos = pos + 1
      } else {
        buf.write_char(c)
      }
      pos = pos + 1
    }
    if pos < len && s.unsafe_get(pos).unsafe_to_char() == '\'' {
      pos = pos + 1
    }
  } else {
    while pos < len {
      let c = s.unsafe_get(pos).unsafe_to_char()
      if c == ':' || c == ',' || c == ')' || c == '(' || is_whitespace(c) {
        break
      }
      buf.write_char(c)
      pos = pos + 1
    }
  }

  let s = buf.to_string()
  (s.replace(old="_".to_string_view(), new=" ".to_string_view()), pos)
}

///|
fn parse_number_rec(
  s : String,
  pos : Int,
  len : Int,
) -> (Double, Int) raise ParseError {
  let mut pos = pos

  let mut result = 0.0
  let mut decimal = 0.0
  let mut decimal_places = 0
  let mut is_negative = false
  let mut is_scientific = false
  let mut sci_sign = 1
  let mut sci_value = 0
  let zero = '0'.to_int()

  if pos < len && s.unsafe_get(pos).to_int() == 45 {
    is_negative = true
    pos = pos + 1
  }

  let mut has_digit = false

  while pos < len {
    let c = s.unsafe_get(pos).to_int()
    if c == 46 {
      pos = pos + 1
      while pos < len {
        let dc = s.unsafe_get(pos).to_int()
        if dc < zero || dc > zero + 9 {
          break
        }
        decimal = decimal * 10.0 + (dc - zero).to_double()
        decimal_places = decimal_places + 1
        has_digit = true
        pos = pos + 1
      }
      break
    }
    if c == 101 || c == 69 {
      is_scientific = true
      pos = pos + 1
      if pos < len {
        let sc = s.unsafe_get(pos).to_int()
        if sc == 45 {
          sci_sign = -1
          pos = pos + 1
        } else if sc == 43 {
          pos = pos + 1
        }
      }
      while pos < len {
        let dc = s.unsafe_get(pos).to_int()
        if dc < zero || dc > zero + 9 {
          break
        }
        sci_value = sci_value * 10 + (dc - zero)
        pos = pos + 1
      }
      break
    }
    if c < zero || c > zero + 9 {
      break
    }
    result = result * 10.0 + (c - zero).to_double()
    has_digit = true
    pos = pos + 1
  }

  if !has_digit {
    raise ParseError::InvalidNumber("".to_string())
  }

  while decimal_places > 0 {
    decimal = decimal / 10.0
    decimal_places = decimal_places - 1
  }

  let mut final_result = result + decimal
  if is_scientific {
    let mut pow10 = 1.0
    let mut j = 0
    while j < sci_value {
      pow10 = pow10 * 10.0
      j = j + 1
    }
    if sci_sign == 1 {
      final_result = final_result * pow10
    } else {
      final_result = final_result / pow10
    }
  }

  if is_negative {
    final_result = -final_result
  }

  (final_result, pos)
}

///|
fn skip_ws(s : String, pos : Int, len : Int) -> Int {
  let mut pos = pos
  while pos < len && is_whitespace(s.unsafe_get(pos).unsafe_to_char()) {
    pos = pos + 1
  }
  pos
}

///|
fn is_whitespace(c : Char) -> Bool {
  c == ' ' || c == '\t' || c == '\n' || c == '\r'
}

///|
fn is_digit(c : Char) -> Bool {
  c >= '0' && c <= '9'
}
