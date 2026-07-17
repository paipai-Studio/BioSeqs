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
  TreeNode::{ name, branch_length, children }
}

///|
pub fn parse_newick(content : String) -> TreeNode {
  let trimmed = content.trim().replace(old=";", new="").to_owned()
  let n = trimmed.length()

  if n == 0 {
    return TreeNode::{ name: "", branch_length: 0.0, children: [] }
  }

  let leaves : Array[TreeNode] = Array::new()
  let mut i = 0

  while i < n {
    let c = char_at(trimmed, i)

    if c == "(" || c == ")" || c == "," {
      i = i + 1
    } else if c == ":" {
      i = i + 1
      while i < n {
        let nc = char_at(trimmed, i)
        if nc == ":" || nc == "," || nc == ")" || nc == "(" {
          break
        }
        i = i + 1
      }
    } else {
      let mut token = ""
      while i < n {
        let nc = char_at(trimmed, i)
        if nc == ":" || nc == "," || nc == ")" || nc == "(" {
          break
        }
        token = token + nc
        i = i + 1
      }

      let mut length = 0.0
      if i < n && char_at(trimmed, i) == ":" {
        i = i + 1
        let mut num = ""
        while i < n {
          let nc = char_at(trimmed, i)
          if nc == ":" || nc == "," || nc == ")" || nc == "(" {
            break
          }
          num = num + nc
          i = i + 1
        }
        length = match parse_double(num) {
          Some(v) => v
          None => 0.0
        }
      }

      leaves.push(TreeNode::{ name: token, branch_length: length, children: [] })
    }
  }

  TreeNode::{ name: "", branch_length: 0.0, children: leaves }
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
  let leaves : Array[String] = Array::new()
  let stack : Array[TreeNode] = Array::new()

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
pub fn get_leaf_count(node : TreeNode) -> Int {
  let mut count = 0
  let stack : Array[TreeNode] = Array::new()

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
pub fn get_node_count(node : TreeNode) -> Int {
  let mut count = 0
  let stack : Array[TreeNode] = Array::new()

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
  let stack : Array[(TreeNode, Double)] = Array::new()

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
  let stack : Array[TreeNode] = Array::new()

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
  let stack : Array[TreeNode] = Array::new()

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
    let new_children : Array[TreeNode] = Array::new()

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
