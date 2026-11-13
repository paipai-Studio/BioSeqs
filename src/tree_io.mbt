///|
/// TreeIO - Phylogenetic tree parsing inspired by Biopython Bio.TreeIO.

///|
pub struct TreeNode {
  name : String
  branch_length : Double
  children : Array[TreeNode]
}

///|
pub fn tree_node(
  name : String,
  branch_length : Double,
  children : Array[TreeNode],
) -> TreeNode {
  TreeNode::{ name, branch_length, children, }
}

///|
pub fn parse_newick(content : String) -> TreeNode {
  let trimmed = content.trim().to_string()
  // Strip optional trailing semicolon.
  let without_semi = if trimmed.length() > 0 && char_at(trimmed, trimmed.length() - 1) == ";" {
    trimmed.substring(start=0, end=trimmed.length() - 1).to_string()
  } else {
    trimmed
  }

  if without_semi.length() == 0 {
    return TreeNode::{ name: "", branch_length: 0.0, children: [], }
  }

  let (root, p) = parse_subtree(without_semi, 0)
  // Consume anything trailing; tolerate whitespace only.
  let _ = skip_ws(without_semi, p)
  root
}

///|
/// Parse one subtree starting at position `p`. Returns the parsed TreeNode
/// and the position just past the consumed segment. Handles both leaves and
/// internal nodes: `((child,...)[label][:len])` or `label[:len]`.
fn parse_subtree(s : String, p : Int32) -> (TreeNode, Int32) {
  let q = skip_ws(s, p)
  if q >= s.length() {
    return (TreeNode::{ name: "", branch_length: 0.0, children: [], }, q)
  }
  let first = char_at(s, q)
  if first == "(" {
    // Internal node — parse children until matching ")".
    let after_open = q + 1
    let mut pos = after_open
    let children : Array[TreeNode] = []
    let mut first_child = true
    while true {
      let c = skip_ws(s, pos)
      if c >= s.length() || char_at(s, c) == ")" {
        pos = c
        break
      }
      if !first_child {
        // Expect a comma between siblings (strict Newick).
        let comma = skip_ws(s, c)
        if comma < s.length() && char_at(s, comma) == "," {
          pos = skip_ws(s, comma + 1)
        } else {
          pos = c
        }
      }
      first_child = false
      let (child, next) = parse_subtree(s, pos)
      children.push(child)
      pos = next
    }
    // Consume closing ")".
    let close = skip_ws(s, pos)
    if close < s.length() && char_at(s, close) == ")" {
      pos = close + 1
    }
    // Optional internal-node label + branch length.
    let (name, p1) = parse_label(s, pos)
    let (length, p2) = parse_branch_length(s, p1)
    (TreeNode::{ name, branch_length: length, children, }, p2)
  } else {
    // Leaf / unparenthesized single root.
    let (name, p1) = parse_label(s, q)
    let (length, p2) = parse_branch_length(s, p1)
    (TreeNode::{ name, branch_length: length, children: [], }, p2)
  }
}

///|
/// Return the first non-whitespace position ≥ p (or s.length() if none).
fn skip_ws(s : String, p : Int32) -> Int32 {
  let n = s.length()
  let mut q = p
  while q < n {
    let c = char_at(s, q)
    if c == " " || c == "\t" || c == "\n" || c == "\r" {
      q = q + 1
    } else {
      break
    }
  }
  q
}

///|
/// Parse an optional node label. Handles quoted ("..." / '...') and unquoted
/// identifiers. Returns (label, position_past_label). If no label is present
/// (next char is `( ) , : ; [` or whitespace), returns ("", p).
fn parse_label(s : String, p : Int32) -> (String, Int32) {
  let n = s.length()
  let q = skip_ws(s, p)
  if q >= n {
    return ("", q)
  }
  let first = char_at(s, q)
  if first == "\"" || first == "'" {
    let quote = first
    let mut pos = q + 1
    let mut label = ""
    while pos < n {
      let c = char_at(s, pos)
      if c == quote {
        pos = pos + 1
        return (label, pos)
      }
      label = label + c
      pos = pos + 1
    }
    return (label, pos)
  }
  // Unquoted label.
  let mut pos = q
  let mut label = ""
  while pos < n {
    let c = char_at(s, pos)
    if c == "(" || c == ")" || c == "," || c == ":" || c == ";" || c == "[" {
      break
    }
    if c == " " || c == "\t" || c == "\n" || c == "\r" {
      break
    }
    label = label + c
    pos = pos + 1
  }
  (label, pos)
}

///|
/// Parse an optional branch length. If `:` is present, consumes it plus the
/// numeric value and returns (length, position_past_number). Otherwise returns
/// (0.0, p). Returns 0.0 on parse failure.
fn parse_branch_length(s : String, p : Int32) -> (Double, Int32) {
  let n = s.length()
  let q = skip_ws(s, p)
  if q >= n || char_at(s, q) != ":" {
    return (0.0, p)
  }
  let mut pos = q + 1
  // Optional whitespace after colon.
  while pos < n {
    let c = char_at(s, pos)
    if c == " " || c == "\t" {
      pos = pos + 1
    } else {
      break
    }
  }
  let start = pos
  while pos < n {
    let c = char_at(s, pos)
    if (c >= "0" && c <= "9") || c == "." || c == "-" || c == "+" ||
       c == "e" || c == "E" {
      pos = pos + 1
    } else {
      break
    }
  }
  let num_str = s.substring(start=start, end=pos).to_string()
  let v = match parse_double(num_str) {
    Some(d) => d
    None => 0.0
  }
  (v, pos)
}

///|
pub fn parse_nhx(content : String) -> TreeNode {
  let n = content.length()

  let mut cleaned = ""
  let mut i = 0
  while i < n {
    if starts_with(content, i, "[&&NHX") {
      while i < n && char_at(content, i) != "]" {
        i = i + 1
      }
      i = i + 1
    } else {
      cleaned = cleaned + char_at(content, i)
      i = i + 1
    }
  }

  parse_newick(cleaned)
}

///|
pub fn write_newick(node : TreeNode, show_lengths : Bool) -> String {
  if node.children.length() == 0 {
    let mut result = node.name
    if show_lengths && node.branch_length > 0.0 {
      result = result + ":" + node.branch_length.to_string()
    }
    result
  } else {
    let mut result = "("

    for i = 0; i < node.children.length(); i = i + 1 {
      result = result + write_newick(node.children[i], show_lengths)
      if i < node.children.length() - 1 {
        result = result + ","
      }
    }

    result = result + ")" + node.name

    if show_lengths && node.branch_length > 0.0 {
      result = result + ":" + node.branch_length.to_string()
    }

    result
  }
}

///|
pub fn get_leaves(node : TreeNode) -> Array[String] {
  let leaves : Array[String] = Array(capacity=0)
  let stack : Array[TreeNode] = Array(capacity=0)

  stack.push(node)

  while stack.length() > 0 {
    let current = match stack.pop() {
      Some(v) => v
      None => node
    }

    if current.children.length() == 0 {
      leaves.push(current.name)
    } else {
      for i = current.children.length() - 1; i >= 0; i = i - 1 {
        stack.push(current.children[i])
      }
    }
  }

  leaves
}

///|
pub fn get_leaf_count(node : TreeNode) -> Int32 {
  let mut count = 0
  let stack : Array[TreeNode] = Array(capacity=0)

  stack.push(node)

  while stack.length() > 0 {
    let current = match stack.pop() {
      Some(v) => v
      None => node
    }

    if current.children.length() == 0 {
      count = count + 1
    } else {
      for i = 0; i < current.children.length(); i = i + 1 {
        stack.push(current.children[i])
      }
    }
  }

  count
}

///|
pub fn get_node_count(node : TreeNode) -> Int32 {
  let mut count = 0
  let stack : Array[TreeNode] = Array(capacity=0)

  stack.push(node)

  while stack.length() > 0 {
    let current = match stack.pop() {
      Some(v) => v
      None => node
    }
    count = count + 1

    for i = 0; i < current.children.length(); i = i + 1 {
      stack.push(current.children[i])
    }
  }

  count
}

///|
pub fn get_tree_depth(node : TreeNode) -> Double {
  let mut max_depth = 0.0
  let stack : Array[(TreeNode, Double)] = Array(capacity=0)

  stack.push((node, 0.0))

  while stack.length() > 0 {
    let (current, current_depth) = match stack.pop() {
      Some(v) => v
      None => (node, 0.0)
    }
    let total_depth = current_depth + current.branch_length

    if current.children.length() == 0 {
      if total_depth > max_depth {
        max_depth = total_depth
      }
    } else {
      for i = 0; i < current.children.length(); i = i + 1 {
        stack.push((current.children[i], total_depth))
      }
    }
  }

  max_depth
}

///|
pub fn total_branch_length(node : TreeNode) -> Double {
  let mut total = 0.0
  let stack : Array[TreeNode] = Array(capacity=0)

  stack.push(node)

  while stack.length() > 0 {
    let current = match stack.pop() {
      Some(v) => v
      None => node
    }
    total = total + current.branch_length

    for i = 0; i < current.children.length(); i = i + 1 {
      stack.push(current.children[i])
    }
  }

  total
}

///|
pub fn find_node(node : TreeNode, name : String) -> TreeNode? {
  let stack : Array[TreeNode] = Array(capacity=0)

  stack.push(node)

  while stack.length() > 0 {
    let current = match stack.pop() {
      Some(v) => v
      None => node
    }

    if current.name == name {
      return Some(current)
    }

    for i = 0; i < current.children.length(); i = i + 1 {
      stack.push(current.children[i])
    }
  }

  None
}

///|
pub fn prune_tree(node : TreeNode, keep_names : Array[String]) -> TreeNode? {
  let keep_map : Map[String, Bool] = Map([], capacity=keep_names.length())
  for i = 0; i < keep_names.length(); i = i + 1 {
    keep_map[keep_names[i]] = true
  }

  prune_node(node, keep_map)
}

///|
fn prune_node(n : TreeNode, keep : Map[String, Bool]) -> TreeNode? {
  if n.children.length() == 0 {
    if keep.get(n.name) is Some(_) {
      Some(n)
    } else {
      None
    }
  } else {
    let new_children : Array[TreeNode] = Array(capacity=0)

    for i = 0; i < n.children.length(); i = i + 1 {
      match prune_node(n.children[i], keep) {
        Some(child) => new_children.push(child)
        None => ignore(())
      }
    }

    if new_children.length() > 0 {
      Some(TreeNode::{
        name: n.name,
        branch_length: n.branch_length,
        children: new_children,
      })
    } else {
      None
    }
  }
}
