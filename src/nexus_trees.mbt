///|
/// nexus_trees - NEXUS Trees/Nodes block: tree collection management.
/// Inspired by Biopython Bio.Nexus.Trees and Bio.Nexus.Nodes modules.
/// Parses the TREES block of a NEXUS file: BEGIN TREES; ... TRANSLATE ... TREE name = newick; END;

///|
/// A single node in a NEXUS tree (Bio.Nexus.Nodes Node analog).
pub(all) struct NexusNode {
  node_id : Int
  mut label : String
  mut branch_length : Double
  mut taxon : String
  mut children : Array[Int]
  mut parent : Int
}

///|
/// A NEXUS tree: root node + node registry (Bio.Nexus.Trees Tree analog).
pub struct NexusTree {
  name : String
  mut rooted : Bool
  mut weight : Double
  nodes : Map[Int, NexusNode]
  mut root : Int
  mut n_taxa : Int
}

///|
/// Collection of trees from a NEXUS TREES block.
pub struct NexusTreeCollection {
  trees : Array[NexusTree]
  translate : Map[String, String] // taxon number -> label
}

///|
/// Create an empty node.
pub fn NexusNode::new(node_id : Int) -> NexusNode {
  NexusNode::{ node_id, label: "", branch_length: 0.0, taxon: "", children: [], parent: -1 }
}

///|
/// Check if a node is a leaf (no children).
pub fn NexusNode::is_leaf(self : NexusNode) -> Bool {
  self.children.length() == 0
}

///|
/// Create an empty tree.
pub fn NexusTree::new(name : String) -> NexusTree {
  NexusTree::{ name, rooted: true, weight: 1.0, nodes: Map([], capacity=16), root: -1, n_taxa: 0 }
}

///|
/// Add a node to the tree.
pub fn NexusTree::add_node(self : NexusTree, node : NexusNode) -> Unit {
  self.nodes[node.node_id] = node
}

///|
/// Get a node by id.
pub fn NexusTree::get_node(self : NexusTree, node_id : Int) -> NexusNode? {
  self.nodes.get(node_id)
}

///|
/// Get all leaf (tip) nodes.
pub fn NexusTree::get_leaves(self : NexusTree) -> Array[NexusNode] {
  let leaves : Array[NexusNode] = []
  for node in self.nodes.values() {
    if node.is_leaf() {
      leaves.push(node)
    }
  }
  leaves
}

///|
/// Count nodes in the tree.
pub fn NexusTree::node_count(self : NexusTree) -> Int {
  self.nodes.length()
}

///|
/// Serialize tree back to Newick format.
pub fn NexusTree::to_newick(self : NexusTree) -> String {
  if self.root < 0 {
    return ";"
  }
  let sb = StringBuilder()
  fn write_node(nid : Int) -> Unit {
    match self.nodes.get(nid) {
      Some(node) => {
        if node.is_leaf() {
          sb.write_string(if node.taxon.length() > 0 { node.taxon } else { node.label })
        } else {
          sb.write_string("(")
          let mut c = 0
          while c < node.children.length() {
            if c > 0 {
              sb.write_string(",")
            }
            write_node(node.children[c])
            c = c + 1
          }
          sb.write_string(")")
          sb.write_string(node.label)
        }
        if node.branch_length > 0.0 {
          sb.write_string(":")
          sb.write_string(node.branch_length.to_string())
        }
      }
      None => ()
    }
  }

  write_node(self.root)
  sb.write_string(";")
  sb.to_string()
}

///|
/// Create an empty tree collection.
pub fn NexusTreeCollection::new() -> NexusTreeCollection {
  NexusTreeCollection::{ trees: [], translate: Map([], capacity=16) }
}

///|
/// Split a Newick string into top-level tokens (respecting parentheses).
fn nexus_trees_split_commas(s : String) -> Array[String] {
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
      parts.push(s[start:i].to_string())
      start = i + 1
    }
    i = i + 1
  }
  if start < s.length() {
    parts.push(s[start:s.length()].to_string())
  }
  parts
}

///|
/// Simple Newick parser into a NexusTree (recursive descent).
fn nexus_trees_parse_newick_inner(
  newick : String,
  tree : NexusTree,
  next_id : Int
) -> (Int, Int) {
  // Returns (node_id, next_available_id)
  let trimmed = newick.trim(chars=" \t\r;").to_string()
  if trimmed.length() == 0 {
    return (-1, next_id)
  }
  if trimmed[0] == 40 {
    // '(' internal node
    // find matching ')'
    let mut depth = 0
    let mut close_pos = -1
    let mut i = 0
    while i < trimmed.length() {
      if trimmed[i] == 40 {
        depth = depth + 1
      } else if trimmed[i] == 41 {
        depth = depth - 1
        if depth == 0 {
          close_pos = i
          i = trimmed.length()
          continue
        }
      }
      i = i + 1
    }
    if close_pos < 0 {
      return (-1, next_id)
    }
    let inner_str = trimmed[1:close_pos].to_string()
    // remainder after ')': label and/or :branch_length
    let rest = trimmed[close_pos + 1:trimmed.length()].to_string()
    // parse label + branch length from rest
    let (label, blen) = nexus_trees_parse_label_branch(rest)
    // create internal node
    let node_id = next_id
    let node = NexusNode::new(node_id)
    node.label = label
    node.branch_length = blen
    tree.add_node(node)
    // parse children
    let child_parts = nexus_trees_split_commas(inner_str)
    let child_ids : Array[Int] = []
    let mut next = node_id + 1
    let mut ci = 0
    while ci < child_parts.length() {
      let (child_id, next2) = nexus_trees_parse_newick_inner(child_parts[ci], tree, next)
      if child_id >= 0 {
        child_ids.push(child_id)
        // set parent
        match tree.nodes.get(child_id) {
          Some(child) => {
            child.parent = node_id
            tree.nodes[child_id] = child
          }
          None => ()
        }
        next = next2
      }
      ci = ci + 1
    }
    // update children
    match tree.nodes.get(node_id) {
      Some(nd) => {
        nd.children = child_ids
        tree.nodes[node_id] = nd
      }
      None => ()
    }
    return (node_id, next)
  } else {
    // leaf: taxon name + optional :branch_length
    let (taxon, blen) = nexus_trees_parse_label_branch(trimmed)
    let node_id = next_id
    let node = NexusNode::new(node_id)
    node.taxon = taxon
    node.branch_length = blen
    tree.add_node(node)
    return (node_id, node_id + 1)
  }
}

///|
/// Parse "label:branch" suffix into (label, branch_length).
fn nexus_trees_parse_label_branch(s : String) -> (String, Double) {
  let trimmed = s.trim(chars=" \t\r").to_string()
  // find last ':' not inside quotes (simplified: first ':')
  let colon = trimmed.find(":")
  match colon {
    Some(pos) => {
      let label = trimmed[0:pos].to_string()
      let branch_str = trimmed[pos + 1:trimmed.length()].to_string()
      let blen = nexus_trees_atof(branch_str)
      (label, blen)
    }
    None => (trimmed, 0.0)
  }
}

///|
/// Simple string to double.
fn nexus_trees_atof(s : String) -> Double {
  let t = s.trim(chars=" \t\r").to_string()
  if t.length() == 0 {
    return 0.0
  }
  let neg = t[0] == 45
  let start = if neg { 1 } else { 0 }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut in_frac = false
  let mut frac_div = 10.0
  let mut k = start
  while k < t.length() {
    let c = t[k]
    if c >= 48 && c <= 57 {
      let digit = (c - 48).to_int().to_double()
      if in_frac {
        frac_part = frac_part + digit / frac_div
        frac_div = frac_div * 10.0
      } else {
        int_part = int_part * 10.0 + digit
      }
    } else if c == 46 {
      in_frac = true
    }
    k = k + 1
  }
  if neg {
    -(int_part + frac_part)
  } else {
    int_part + frac_part
  }
}

///|
/// Parse a single TREE line: "tree name = newick" or "TREE name = [&R] newick".
fn nexus_trees_parse_tree_line(line : String) -> NexusTree? {
  let trimmed = line.trim(chars=" \t\r;").to_string()
  if trimmed.length() == 0 {
    return None
  }
  // find '='
  match trimmed.find("=") {
    Some(eq_pos) => {
      let left = trimmed[0:eq_pos].trim(chars=" \t").to_string()
      let right = trimmed[eq_pos + 1:trimmed.length()].trim(chars=" \t").to_string()
      // left: "tree name" or "translate"
      let left_fields = left.split(" ").to_array()
      let kw = if left_fields.length() > 0 { left_fields[0].to_string() } else { "" }
      if kw == "tree" || kw == "TREE" {
        let name = if left_fields.length() > 1 { left_fields[1].to_string() } else { "tree" }
        let tree = NexusTree::new(name)
        // check rooted flag [&R] or [&U]
        let newick_part = right
        if right.length() >= 4 && right[0] == 91 && right[1] == 38 {
          // [&R] or [&U] prefix
          let rb = right.find("]")
          match rb {
            Some(pos) => {
              let flag = right[0:pos].to_string()
              tree.rooted = flag.contains("R") && not(flag.contains("U"))
              tree.weight = 1.0
              let rest_start = pos + 1
              let newick_part = right[rest_start:right.length()].trim(chars=" \t").to_string()
              let (root_id, _) = nexus_trees_parse_newick_inner(newick_part, tree, 0)
              tree.root = root_id
              let leaves = tree.get_leaves()
              tree.n_taxa = leaves.length()
              return Some(tree)
            }
            None => return None
          }
        }
        let (root_id, _) = nexus_trees_parse_newick_inner(newick_part, tree, 0)
        tree.root = root_id
        let leaves = tree.get_leaves()
        tree.n_taxa = leaves.length()
        return Some(tree)
      }
      None
    }
    None => None
  }
}

///|
/// Parse a NEXUS TREES block into a tree collection.
pub fn nexus_trees_parse(content : String) -> NexusTreeCollection {
  let collection = NexusTreeCollection::new()
  let lines = content.split("\n").to_array()
  let in_trees = false
  let mut in_trees_block = in_trees
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim(chars=" \t\r").to_string()
    let upper = line.to_lower()
    if upper.contains("begin trees") {
      in_trees_block = true
    } else if upper.contains("end;") && in_trees_block {
      in_trees_block = false
    } else if in_trees_block {
      if upper.contains("translate") {
        // TRANSLATE lines: "1 Homo_sapiens," etc (until ';')
        let mut j = i
        while j < lines.length() {
          let tline = lines[j].trim(chars=" \t\r").to_string()
          if j > i {
            let tupper = tline.to_lower()
            if tupper.contains("end;") {
              j = lines.length()
              break
            }
          }
          // parse "num label," entries
          let clean = tline.trim(chars=" \t\r,;").to_string()
          if clean.length() > 0 && not(clean.to_lower().contains("translate")) {
            let fields = clean.split(" ").to_array()
            if fields.length() >= 2 {
              let num = fields[0].to_string()
              let label = fields[1].to_string()
              collection.translate[num] = label
            }
          }
          if tline.contains(";") {
            i = j
            break
          }
          j = j + 1
        }
      } else if upper.contains("tree") {
        match nexus_trees_parse_tree_line(line) {
          Some(tree) => collection.trees.push(tree)
          None => ()
        }
      }
    }
    i = i + 1
  }
  collection
}

///|
/// Parse multiple Newick strings as a collection (no NEXUS wrapper).
pub fn nexus_trees_from_newicks(newicks : Array[(String, String)]) -> NexusTreeCollection {
  let collection = NexusTreeCollection::new()
  let mut i = 0
  while i < newicks.length() {
    let (name, newick) = newicks[i]
    let line = "tree " + name + " = " + newick
    match nexus_trees_parse_tree_line(line) {
      Some(tree) => collection.trees.push(tree)
      None => ()
    }
    i = i + 1
  }
  collection
}

///|
/// Find a tree by name in the collection.
pub fn NexusTreeCollection::get_tree(self : NexusTreeCollection, name : String) -> NexusTree? {
  let mut i = 0
  while i < self.trees.length() {
    if self.trees[i].name == name {
      return Some(self.trees[i])
    }
    i = i + 1
  }
  None
}

///|
/// Resolve translated taxon labels (e.g. "1" -> "Homo_sapiens").
pub fn NexusTreeCollection::resolve_taxon(self : NexusTreeCollection, taxon : String) -> String {
  match self.translate.get(taxon) {
    Some(label) => label
    None => taxon
  }
}

///|
/// Compute patristic distance between two leaves in a tree.
pub fn NexusTree::patristic_distance(self : NexusTree, leaf_a : String, leaf_b : String) -> Double {
  // Collect ancestor paths with branch lengths
  fn path_to_root(nid : Int) -> Array[(Int, Double)] {
    let path : Array[(Int, Double)] = []
    let mut cur = nid
    while cur >= 0 {
      match self.nodes.get(cur) {
        Some(node) => {
          path.push((cur, node.branch_length))
          cur = node.parent
        }
        None => {
          cur = -1
        }
      }
    }
    path
  }

  // Find leaf node ids by taxon
  let mut id_a = -1
  let mut id_b = -1
  for node in self.nodes.values() {
    if node.taxon == leaf_a {
      id_a = node.node_id
    }
    if node.taxon == leaf_b {
      id_b = node.node_id
    }
  }
  if id_a < 0 || id_b < 0 {
    return 0.0
  }
  let path_a = path_to_root(id_a)
  let path_b = path_to_root(id_b)
  // Find MRCA
  let mut mrca = -1
  let mut ai = 0
  while ai < path_a.length() {
    let mut bi = 0
    while bi < path_b.length() {
      if path_a[ai].0 == path_b[bi].0 {
        mrca = path_a[ai].0
        ai = path_a.length()
        break
      }
      bi = bi + 1
    }
    ai = ai + 1
  }
  if mrca < 0 {
    return 0.0
  }
  // Sum branch lengths from each leaf up to (excluding) MRCA
  let mut dist = 0.0
  let mut ai2 = 0
  while ai2 < path_a.length() {
    if path_a[ai2].0 != mrca {
      dist = dist + path_a[ai2].1
    } else {
      break
    }
    ai2 = ai2 + 1
  }
  let mut bi2 = 0
  while bi2 < path_b.length() {
    if path_b[bi2].0 != mrca {
      dist = dist + path_b[bi2].1
    } else {
      break
    }
    bi2 = bi2 + 1
  }
  dist
}

///|
/// Compute total tree length (sum of all branch lengths).
pub fn NexusTree::total_length(self : NexusTree) -> Double {
  let mut total = 0.0
  for node in self.nodes.values() {
    total = total + node.branch_length
  }
  total
}

///|
/// Example NEXUS content with a TREES block.
pub fn nexus_trees_example_content() -> String {
  "#NEXUS\nBEGIN TREES;\n  TREE tree1 = ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\nEND;\n"
}