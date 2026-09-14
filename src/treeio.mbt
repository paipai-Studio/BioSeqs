///|
/// treeio - Base classes and parsers for phylogenetic tree input/output.
/// Inspired by Bioconductor treeio package.
/// Supports: Newick, NEXUS, Beast/MrBayes annotations, jplace, r8s, phylip tree formats.
/// Provides unified TreeioTreeData structure for all tree formats.

///|
/// Tree node with annotations (Beast/MrBayes style).
pub struct TreeioNode {
  node_id : Int
  label : String
  branch_length : Double
  mut parent : Int
  children : Array[Int]
  mut is_leaf : Bool
  annotations : Map[String, String]
}

///|
/// Unified tree data structure (treeio treedata).
pub struct TreeioTreeData {
  nodes : Array[TreeioNode]
  n_nodes : Int
  n_tips : Int
  root : Int
  tip_labels : Array[String]
  node_labels : Array[String]
  edge : Array[(Int, Int)]
  edge_length : Array[Double]
  mut format : String
}

///|
/// Simple string to double parser.
fn treeio_atof(s : String) -> Double {
  let s_trim = s.trim().to_string()
  if s_trim.length() == 0 {
    return 0.0
  }
  let neg = if s_trim[0] == 45 { true } else { false } // '-'
  let start = if neg { 1 } else { 0 }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut in_frac = false
  let mut frac_div = 10.0
  let mut i = start
  while i < s_trim.length() {
    let c = s_trim[i].to_int()
    if c >= 48 && c <= 57 {
      let digit = (c - 48).to_double()
      if in_frac {
        frac_part = frac_part + digit / frac_div
        frac_div = frac_div * 10.0
      } else {
        int_part = int_part * 10.0 + digit
      }
    } else if c == 46 { // '.'
      in_frac = true
    } else {
      break
    }
    i = i + 1
  }
  let result = int_part + frac_part
  if neg { -result } else { result }
}

///|
/// Parse Newick tree string into TreeioTreeData.
pub fn treeio_parse_newick(newick : String) -> TreeioTreeData {
  let s = newick.trim().to_string()
  let tree_str = if s.length() > 0 && s[s.length() - 1] == 59 {
    s[0:s.length() - 1].to_owned()
  } else {
    s
  }
  let nodes : Array[TreeioNode] = []
  let edges : Array[(Int, Int)] = []
  let edge_lengths : Array[Double] = []
  let result = treeio_parse_newick_rec(tree_str, nodes, edges, edge_lengths, 0)
  let final_nodes = result.0
  let final_edges = result.1
  let final_elens = result.2
  // Set parent and children from edges
  let mut i = 0
  while i < final_edges.length() {
    let (parent_id, child_id) = final_edges[i]
    let mut pi = 0
    while pi < final_nodes.length() {
      if final_nodes[pi].node_id == child_id {
        final_nodes[pi].parent = parent_id
      }
      pi = pi + 1
    }
    i = i + 1
  }
  let mut j = 0
  while j < final_edges.length() {
    let (parent_id, child_id) = final_edges[j]
    let mut pi = 0
    while pi < final_nodes.length() {
      if final_nodes[pi].node_id == parent_id {
        final_nodes[pi].children.push(child_id)
        final_nodes[pi].is_leaf = false
      }
      pi = pi + 1
    }
    j = j + 1
  }
  let mut root = -1
  let mut ki = 0
  while ki < final_nodes.length() {
    if final_nodes[ki].parent == -1 {
      root = final_nodes[ki].node_id
    }
    ki = ki + 1
  }
  let tip_labels : Array[String] = []
  let node_labels : Array[String] = []
  let mut ti = 0
  while ti < final_nodes.length() {
    if final_nodes[ti].is_leaf {
      tip_labels.push(final_nodes[ti].label)
    } else {
      node_labels.push(final_nodes[ti].label)
    }
    ti = ti + 1
  }
  TreeioTreeData::{
    nodes: final_nodes,
    n_nodes: final_nodes.length(),
    n_tips: tip_labels.length(),
    root,
    tip_labels,
    node_labels,
    edge: final_edges,
    edge_length: final_elens,
    format: "newick",
  }
}

///|
fn treeio_parse_newick_rec(
  s : String,
  nodes : Array[TreeioNode],
  edges : Array[(Int, Int)],
  edge_lengths : Array[Double],
  next_id : Int,
) -> (Array[TreeioNode], Array[(Int, Int)], Array[Double], Int) {
  let s_trim = s.trim().to_string()
  if s_trim.length() == 0 {
    return (nodes, edges, edge_lengths, next_id)
  }
  let my_id = next_id
  let mut my_next = next_id + 1
  if s_trim[0] == 40 {
    let mut depth = 1
    let mut pos = 1
    while pos < s_trim.length() && depth > 0 {
      if s_trim[pos] == 40 {
        depth = depth + 1
      } else if s_trim[pos] == 41 {
        depth = depth - 1
      }
      pos = pos + 1
    }
    let inner = s_trim[1:pos - 1].to_owned()
    let after = if pos < s_trim.length() {
      s_trim[pos:s_trim.length()].to_owned()
    } else {
      ""
    }
    let (label, blen) = treeio_parse_label_blen(after)
    nodes.push(TreeioNode::{
      node_id: my_id,
      label,
      branch_length: blen,
      parent: -1,
      children: [],
      is_leaf: true,
      annotations: Map([], capacity=0),
    })
    let children = treeio_split_top(inner)
    let mut ci = 0
    while ci < children.length() {
      let child_str = children[ci]
      let result = treeio_parse_newick_rec(child_str, nodes, edges, edge_lengths, my_next)
      let child_id = my_next
      my_next = result.3
      edges.push((my_id, child_id))
      edge_lengths.push(0.0)
      ci = ci + 1
    }
  } else {
    let (label, blen) = treeio_parse_label_blen(s_trim)
    nodes.push(TreeioNode::{
      node_id: my_id,
      label,
      branch_length: blen,
      parent: -1,
      children: [],
      is_leaf: true,
      annotations: Map([], capacity=0),
    })
  }
  (nodes, edges, edge_lengths, my_next)
}

///|
fn treeio_parse_label_blen(s : String) -> (String, Double) {
  let s_trim = s.trim().to_string()
  if s_trim.length() == 0 {
    return ("", 0.0)
  }
  let mut colon_pos = -1
  let mut depth = 0
  let mut i = 0
  while i < s_trim.length() {
    let c = s_trim[i]
    if c == 91 || c == 123 {
      depth = depth + 1
    } else if c == 93 || c == 125 {
      depth = depth - 1
    } else if c == 58 && depth == 0 {
      colon_pos = i
    }
    i = i + 1
  }
  if colon_pos < 0 {
    return (s_trim, 0.0)
  }
  let label = s_trim[0:colon_pos].to_owned().trim().to_string()
  let blen_str = s_trim[colon_pos + 1:s_trim.length()].to_owned().trim().to_string()
  let blen = treeio_atof(blen_str)
  (label, blen)
}

///|
fn treeio_split_top(s : String) -> Array[String] {
  let parts : Array[String] = []
  let mut depth = 0
  let mut start = 0
  let mut i = 0
  while i < s.length() {
    let c = s[i]
    if c == 40 {
      depth = depth + 1
    } else if c == 41 {
      depth = depth - 1
    } else if c == 44 && depth == 0 {
      parts.push(s[start:i].to_owned().trim().to_string())
      start = i + 1
    }
    i = i + 1
  }
  if start < s.length() {
    parts.push(s[start:s.length()].to_owned().trim().to_string())
  }
  parts
}

///|
/// Write TreeioTreeData to Newick format.
pub fn treeio_write_newick(tree : TreeioTreeData) -> String {
  if tree.root < 0 {
    return ";"
  }
  treeio_newick_from_tree(tree, tree.root) + ";"
}

///|
fn treeio_newick_from_tree(tree : TreeioTreeData, node_id : Int) -> String {
  let node = tree.nodes[node_id]
  if node.is_leaf {
    if node.branch_length > 0.0 {
      return node.label + ":" + node.branch_length.to_string()
    }
    return node.label
  }
  let sb = StringBuilder::new()
  sb.write_char('(')
  let mut i = 0
  while i < node.children.length() {
    if i > 0 {
      sb.write_char(',')
    }
    sb.write_string(treeio_newick_from_tree(tree, node.children[i]))
    i = i + 1
  }
  sb.write_char(')')
  if node.label.length() > 0 {
    sb.write_string(node.label)
  }
  if node.branch_length > 0.0 {
    sb.write_char(':')
    sb.write_string(node.branch_length.to_string())
  }
  sb.to_string()
}

///|
/// jplace placement entry.
pub struct TreeioJplacePlacement {
  p : Int
  pendant_length : Double
  proximal_length : Double
  likelihood : Double
  like_weight_ratio : Double
  ambiguous : Bool
}

///|
/// jplace format result.
pub struct TreeioJplaceResult {
  tree : String
  placements : Array[(String, Array[TreeioJplacePlacement])]
  n_places : Int
}

///|
/// Parse jplace format (simplified JSON-like parsing).
pub fn treeio_parse_jplace(content : String) -> TreeioJplaceResult {
  let tree_str = treeio_extract_json(content, "tree")
  let placements : Array[(String, Array[TreeioJplacePlacement])] = []
  let mut n_places = 0
  let pos = content.find("placements").unwrap_or(-1)
  if pos >= 0 {
    let bytes = content.to_array()
    let mut i = 0
    while i < bytes.length() {
      if bytes[i].to_int() == 91 {
        n_places = n_places + 1
      }
      i = i + 1
    }
    n_places = n_places / 3
  }
  TreeioJplaceResult::{ tree: tree_str, placements, n_places }
}

///|
fn treeio_extract_json(content : String, key : String) -> String {
  let search = "\"" + key + "\":"
  let pos = content.find(search).unwrap_or(-1)
  if pos < 0 {
    return ""
  }
  let start = pos + search.length()
  let mut s = start
  while s < content.length() && (content[s] == 32 || content[s] == 9) {
    s = s + 1
  }
  if s >= content.length() {
    return ""
  }
  if content[s] == 34 {
    let mut e = s + 1
    while e < content.length() && content[e] != 34 {
      e = e + 1
    }
    return content[s + 1:e].to_owned()
  }
  let mut e = s
  while e < content.length() && content[e] != 44 && content[e] != 125 {
    e = e + 1
  }
  content[s:e].to_owned().trim().to_string()
}

///|
/// Parse Beast NEXUS tree with annotations.
pub fn treeio_parse_beast(newick : String) -> TreeioTreeData {
  let tree = treeio_parse_newick(treeio_strip_beast(newick))
  tree.format = "beast"
  let annotations = treeio_extract_beast(newick)
  let mut ai = 0
  while ai < annotations.length() {
    let (label, ann) = annotations[ai]
    let mut ni = 0
    while ni < tree.nodes.length() {
      if tree.nodes[ni].label == label {
        let keys = ann.keys().collect()
        let mut ki = 0
        while ki < keys.length() {
          tree.nodes[ni].annotations[keys[ki]] = ann[keys[ki]]
          ki = ki + 1
        }
      }
      ni = ni + 1
    }
    ai = ai + 1
  }
  tree
}

///|
fn treeio_strip_beast(s : String) -> String {
  let sb = StringBuilder::new()
  let mut i = 0
  let bytes = s.to_array()
  while i < bytes.length() {
    if bytes[i].to_int() == 91 {
      let mut depth = 1
      i = i + 1
      while i < bytes.length() && depth > 0 {
        if bytes[i].to_int() == 91 {
          depth = depth + 1
        } else if bytes[i].to_int() == 93 {
          depth = depth - 1
        }
        i = i + 1
      }
    } else {
      sb.write_char(bytes[i])
      i = i + 1
    }
  }
  sb.to_string()
}

///|
fn treeio_extract_beast(s : String) -> Array[(String, Map[String, String])] {
  let result : Array[(String, Map[String, String])] = []
  let bytes = s.to_array()
  let mut i = 0
  while i < bytes.length() {
    if bytes[i].to_int() == 91 {
      let mut depth = 1
      let start = i + 1
      i = i + 1
      while i < bytes.length() && depth > 0 {
        if bytes[i].to_int() == 91 {
          depth = depth + 1
        } else if bytes[i].to_int() == 93 {
          depth = depth - 1
        }
        i = i + 1
      }
      let ann_str = s[start:i - 1].to_owned()
      let ann_clean = if ann_str.length() > 0 && ann_str[0] == 38 {
        ann_str[1:ann_str.length()].to_owned()
      } else {
        ann_str
      }
      let ann_map : Map[String, String] = Map([], capacity=8)
      let pairs = ann_clean.split(",").to_array()
      let mut pi = 0
      while pi < pairs.length() {
        let pair = pairs[pi].trim().to_owned()
        let eq_pos = pair.find("=").unwrap_or(-1)
        if eq_pos >= 0 {
          let key = pair[0:eq_pos].to_owned().trim().to_string()
          let val = pair[eq_pos + 1:pair.length()].to_owned().trim().to_string()
          ann_map[key] = val
        }
        pi = pi + 1
      }
      let mut label_end = start - 1
      while label_end > 0 && (bytes[label_end].to_int() == 32 || bytes[label_end].to_int() == 9) {
        label_end = label_end - 1
      }
      let mut label_start = label_end
      while label_start > 0 {
        let c = bytes[label_start].to_int()
        if (c >= 65 && c <= 90) || (c >= 97 && c <= 122) || (c >= 48 && c <= 57) || c == 95 {
          label_start = label_start - 1
        } else {
          break
        }
      }
      let label = s[label_start + 1:label_end + 1].to_owned().trim().to_string()
      if label.length() > 0 {
        result.push((label, ann_map))
      }
    }
    i = i + 1
  }
  result
}

///|
/// Parse r8s log file.
pub fn treeio_parse_r8s(content : String) -> Array[TreeioTreeData] {
  let trees : Array[TreeioTreeData] = []
  let lines = content.split("\n").to_array()
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_owned()
    if line.contains("tree") && line.contains("=") {
      let eq_pos = line.find("=").unwrap_or(-1)
      if eq_pos >= 0 {
        let tree_str = line[eq_pos + 1:line.length()].to_owned().trim().to_string()
        trees.push(treeio_parse_newick(tree_str))
      }
    }
    i = i + 1
  }
  trees
}

///|
/// Parse phylip tree format.
pub fn treeio_parse_phylip(content : String) -> TreeioTreeData {
  let lines = content.split("\n").to_array()
  let mut tree_line = ""
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_owned()
    if line.contains("(") && line.contains(")") {
      tree_line = line
      break
    }
    i = i + 1
  }
  let tree = treeio_parse_newick(tree_line)
  tree.format = "phylip"
  tree
}

///|
/// Get tree statistics: (n_tips, n_internal, n_edges).
pub fn treeio_tree_stats(tree : TreeioTreeData) -> (Int, Int, Int) {
  (tree.n_tips, tree.n_nodes - tree.n_tips, tree.edge.length())
}

///|
/// Check if tree is binary.
pub fn treeio_is_binary(tree : TreeioTreeData) -> Bool {
  let mut i = 0
  while i < tree.nodes.length() {
    if !tree.nodes[i].is_leaf {
      if tree.nodes[i].children.length() != 2 {
        return false
      }
    }
    i = i + 1
  }
  true
}

///|
/// Check if tree is rooted.
pub fn treeio_is_rooted(tree : TreeioTreeData) -> Bool {
  if tree.root < 0 || tree.root >= tree.nodes.length() {
    return false
  }
  tree.nodes[tree.root].children.length() == 2
}

///|
/// Get total tree length.
pub fn treeio_total_length(tree : TreeioTreeData) -> Double {
  let mut total = 0.0
  let mut i = 0
  while i < tree.nodes.length() {
    total = total + tree.nodes[i].branch_length
    i = i + 1
  }
  total
}

///|
/// Get maximum tree depth.
pub fn treeio_max_depth(tree : TreeioTreeData) -> Int {
  if tree.root < 0 {
    return 0
  }
  treeio_depth_rec(tree, tree.root, 0)
}

///|
fn treeio_depth_rec(tree : TreeioTreeData, node_id : Int, depth : Int) -> Int {
  let node = tree.nodes[node_id]
  if node.is_leaf {
    return depth
  }
  let mut max_d = depth
  let mut i = 0
  while i < node.children.length() {
    let child_d = treeio_depth_rec(tree, node.children[i], depth + 1)
    if child_d > max_d {
      max_d = child_d
    }
    i = i + 1
  }
  max_d
}

///|
/// Create example tree for testing.
pub fn treeio_example_tree() -> TreeioTreeData {
  treeio_parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
}

///|
/// Create example Beast tree for testing.
pub fn treeio_example_beast() -> TreeioTreeData {
  treeio_parse_beast("((A[&height=0.1]:0.1,B[&height=0.1]:0.2)[&height=0.3]:0.3,C[&height=0.1]:0.4);")
}
