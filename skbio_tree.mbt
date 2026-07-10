///|
/// scikit-bio Tree module - Phylogenetic tree data structures.

///|
/// A node within a phylogenetic tree.
/// Stores links to parent and child nodes, branch length, and support value.
pub struct TreeNode {
  name : String?
  length : Double?
  support : Double?
  pub children : Array[TreeNode]
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
pub fn TreeNode::lca(self : TreeNode, names : Array[String]) -> TreeNode? {
  if names.length() == 0 {
    return None
  }
  
  let nodes : Array[TreeNode] = Array::new()
  for name in names {
    match self.find(name) {
      Some(n) => nodes.push(n)
      None => return None
    }
  }
  
  find_lca(self, nodes)
}

///|
fn find_lca(root : TreeNode, nodes : Array[TreeNode]) -> TreeNode? {
  if nodes.length() == 0 {
    return None
  }
  
  if nodes.length() == 1 {
    return Some(nodes[0])
  }
  
  find_lca_helper(root, nodes)
}

///|
fn find_lca_helper(current : TreeNode, nodes : Array[TreeNode]) -> TreeNode? {
  let mut found_all = true
  for node in nodes {
    if !is_descendant(current, node) {
      found_all = false
      break
    }
  }
  
  if !found_all {
    return None
  }
  
  for child in current.children {
    let mut child_has_all = true
    for node in nodes {
      if !is_descendant(child, node) {
        child_has_all = false
        break
      }
    }
    if child_has_all {
      return find_lca_helper(child, nodes)
    }
  }
  
  Some(current)
}

///|
fn is_descendant(ancestor : TreeNode, descendant : TreeNode) -> Bool {
  if ancestor.name == descendant.name && 
     ancestor.children.length() == descendant.children.length() {
    return true
  }
  for child in ancestor.children {
    if is_descendant(child, descendant) {
      return true
    }
  }
  false
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
fn TreeNode::postorder_helper(self : TreeNode, result : Array[TreeNode]) -> Unit {
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
fn TreeNode::collect_non_tips(self : TreeNode, result : Array[TreeNode]) -> Unit {
  if !self.is_tip() {
    result.push(self)
    for child in self.children {
      child.collect_non_tips(result)
    }
  }
}

///|
/// Calculate the distance between two nodes by their names.
pub fn TreeNode::distance(self : TreeNode, name1 : String, name2 : String) -> Double? {
  match (self.find(name1), self.find(name2)) {
    (Some(node1), Some(node2)) => {
      let names = Array::new()
      names.push(name1)
      names.push(name2)
      let lca = self.lca(names)
      match lca {
        Some(ancestor) => {
          let dist1 = distance_to_ancestor(node1, ancestor, 0.0)
          let dist2 = distance_to_ancestor(node2, ancestor, 0.0)
          Some(dist1 + dist2)
        }
        None => None
      }
    }
    _ => None
  }
}

///|
fn distance_to_ancestor(node : TreeNode, ancestor : TreeNode, current : Double) -> Double {
  if node.name == ancestor.name && 
     node.children.length() == ancestor.children.length() {
    return current
  }
  for child in ancestor.children {
    if is_descendant(child, node) {
      return distance_to_ancestor(node, child, current + child.length.unwrap_or(0.0))
    }
  }
  current
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
pub fn TreeNode::find_by_name(self : TreeNode, name : String) -> Array[TreeNode] {
  let result = Array::new()
  self.find_by_name_helper(name, result)
  result
}

///|
fn TreeNode::find_by_name_helper(self : TreeNode, name : String, result : Array[TreeNode]) -> Unit {
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
fn parse_node_rec(s : String, pos : Int, len : Int) -> (TreeNode, Int) raise ParseError {
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
  
  let node = TreeNode::new(name=name, length=length, children=children)
  
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
      if c == '\'' && pos + 1 < len && s.unsafe_get(pos + 1).unsafe_to_char() == '\'' {
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
fn parse_number_rec(s : String, pos : Int, len : Int) -> (Double, Int) raise ParseError {
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


