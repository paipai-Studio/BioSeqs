///|
/// scikit-bio Tree module - Phylogenetic tree data structures.

///|
/// A node within a phylogenetic tree.
/// Stores links to parent and child nodes, branch length, and support value.
pub struct TreeNode {
  name : String?
  length : Double?
  support : Double?
  children : Array[TreeNode]
}

///|
/// Construct a TreeNode with optional attributes.
pub fn TreeNode::new(
  name? : String? = None,
  length? : Double? = None,
  support? : Double? = None,
  children? : Array[TreeNode] = [],
) -> TreeNode {
  TreeNode::{ name, length, support, children }
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
/// Find the lowest common ancestor of nodes by names.
/// Uses backtracking path search to avoid intermediate array allocations.
pub fn TreeNode::lca(self : TreeNode, names : Array[String]) -> TreeNode? {
  if names.length() == 0 {
    return None
  }

  // Find paths from root to each target node using backtracking
  let paths : Array[Array[TreeNode]] = Array::new()
  for name in names {
    let path = Array::new()
    if find_path_backtrack(self, name, path) {
      paths.push(path)
    } else {
      return None
    }
  }

  // Find the last common node across all paths
  let min_len = {
    let mut m = paths[0].length()
    for i in 1..<paths.length() {
      if paths[i].length() < m {
        m = paths[i].length()
      }
    }
    m
  }

  let mut lca_idx = 0
  for i in 0..<min_len {
    let candidate = paths[0][i]
    let mut all_match = true
    for j in 1..<paths.length() {
      if paths[j][i].name != candidate.name {
        all_match = false
        break
      }
    }
    if all_match {
      lca_idx = i
    } else {
      break
    }
  }

  Some(paths[0][lca_idx])
}

///|
/// Backtracking path search: appends nodes to path as it descends,
/// removes them if the branch doesn't contain the target.
fn find_path_backtrack(
  node : TreeNode,
  name : String,
  path : Array[TreeNode],
) -> Bool {
  path.push(node)
  if node.name == Some(name) {
    return true
  }
  for child in node.children {
    if find_path_backtrack(child, name, path) {
      return true
    }
  }
  ignore(path.remove(path.length() - 1))
  false
}

///|
/// Find the path from this node to a node with the given name.
/// Returns the path as an array of TreeNodes from root to target.
pub fn TreeNode::find_path(self : TreeNode, name : String) -> Array[TreeNode]? {
  let path = Array::new()
  if find_path_backtrack(self, name, path) {
    Some(path)
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
/// Uses index-based queue for O(n) performance.
pub fn TreeNode::levelorder(self : TreeNode) -> Array[TreeNode] {
  let result = Array::new()
  let queue = Array::new()
  queue.push(self)
  let mut head = 0
  while head < queue.length() {
    let node = queue[head]
    head = head + 1
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
/// Calculate the distance between two nodes by their names.
/// Uses path-based approach: finds paths to both nodes, identifies LCA,
/// then sums branch lengths from LCA to each node.
pub fn TreeNode::distance(
  self : TreeNode,
  name1 : String,
  name2 : String,
) -> Double? {
  match (self.find_path(name1), self.find_path(name2)) {
    (Some(path1), Some(path2)) => {
      // Find LCA index (last common node on both paths)
      let mut lca_idx = 0
      let min_len = if path1.length() < path2.length() {
        path1.length()
      } else {
        path2.length()
      }
      for i in 0..<min_len {
        if path1[i].name == path2[i].name {
          lca_idx = i
        } else {
          break
        }
      }
      // Sum distances from LCA to each node
      let mut dist = 0.0
      for i in (lca_idx + 1)..<path1.length() {
        dist = dist + path1[i].length.unwrap_or(0.0)
      }
      for i in (lca_idx + 1)..<path2.length() {
        dist = dist + path2[i].length.unwrap_or(0.0)
      }
      Some(dist)
    }
    _ => None
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
  parse_tree_node_from_newick(newick)
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
fn parse_tree_node_from_newick(input : String) -> TreeNode raise ParseError {
  let trimmed = input.trim().to_owned()

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

  (buf.to_string(), pos)
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

  // Apply decimal places using a single division
  if decimal_places > 0 {
    let mut divisor = 1.0
    for i = 0; i < decimal_places; i = i + 1 {
      divisor = divisor * 10.0
    }
    decimal = decimal / divisor
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
