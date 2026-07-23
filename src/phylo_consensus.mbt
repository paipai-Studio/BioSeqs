///|
/// Bio.Phylo.Consensus - Consensus tree building algorithms.
/// Inspired by Biopython Bio.Phylo.Consensus.
/// Provides majority rule consensus and strict consensus tree building.

///|
pub struct ConsensusTree {
  nodes : Array[ConsensusNode]
  root : ConsensusNode?
}

///|
pub struct ConsensusNode {
  name : String
  children : Array[ConsensusNode]
  support : Double
  is_internal : Bool
}

///|
pub fn ConsensusNode::new(name : String, is_internal : Bool) -> ConsensusNode {
  ConsensusNode::{
    name,
    children: Array::new(),
    support: 0.0,
    is_internal: is_internal
  }
}

///|
pub fn ConsensusNode::add_child(self : ConsensusNode, child : ConsensusNode) -> ConsensusNode {
  let new_children = self.children
  new_children.push(child)
  ConsensusNode::{
    name: self.name,
    children: new_children,
    support: self.support,
    is_internal: self.is_internal
  }
}

///|
pub fn ConsensusNode::set_support(self : ConsensusNode, support : Double) -> ConsensusNode {
  ConsensusNode::{
    name: self.name,
    children: self.children,
    support: support,
    is_internal: self.is_internal
  }
}

///|
pub fn ConsensusNode::is_leaf(self : ConsensusNode) -> Bool {
  !self.is_internal && self.children.length() == 0
}

///|
pub fn ConsensusNode::get_leaves(self : ConsensusNode) -> Array[String] {
  let result : Array[String] = Array::new()
  ignore(get_leaves_recursive(self, result))
  result
}

///|
fn get_leaves_recursive(node : ConsensusNode, result : Array[String]) -> Unit {
  if node.is_leaf() {
    result.push(node.name)
  } else {
    let mut i = 0
    while i < node.children.length() {
      ignore(get_leaves_recursive(node.children[i], result))
      i = i + 1
    }
  }
}

///|
pub struct Split {
  taxa : Array[String]
  support : Double
}

///|
pub fn Split::new(taxa : Array[String]) -> Split {
  Split::{
    taxa,
    support: 0.0
  }
}

///|
pub fn Split::set_support(self : Split, support : Double) -> Split {
  Split::{
    taxa: self.taxa,
    support: support
  }
}

///|
pub fn Split::normalize(self : Split) -> Split {
  let sorted = sort_strings(self.taxa)
  Split::{
    taxa: sorted,
    support: self.support
  }
}

///|
fn sort_taxa_strings(arr : Array[String]) -> Array[String] {
  let n = arr.length()
  let result = Array::make(n, "")
  
  let mut i = 0
  while i < n {
    result[i] = arr[i]
    i = i + 1
  }
  
  i = 0
  while i < n - 1 {
    let mut j = i + 1
    while j < n {
      if string_gt(result[i], result[j]) {
        let tmp = result[i]
        result[i] = result[j]
        result[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  
  result
}

///|
fn string_gt(a : String, b : String) -> Bool {
  let n = a.length()
  let m = b.length()
  let min_len = if n < m { n } else { m }
  let mut i = 0
  while i < min_len {
    let ca = a.unsafe_get(i).to_int()
    let cb = b.unsafe_get(i).to_int()
    if ca > cb {
      return true
    } else if ca < cb {
      return false
    }
    i = i + 1
  }
  n > m
}

///|
pub fn Split::hash(self : Split) -> String {
  let mut result = self.taxa[0]
  let mut i = 1
  while i < self.taxa.length() {
    result = result + "," + self.taxa[i]
    i = i + 1
  }
  result
}

///|
pub fn get_all_splits(node : ConsensusNode, all_taxa : Array[String]) -> Array[Split] {
  let result : Array[Split] = Array::new()
  ignore(get_splits_recursive(node, all_taxa, result))
  result
}

///|
fn get_splits_recursive(node : ConsensusNode, all_taxa : Array[String], result : Array[Split]) -> Unit {
  if !node.is_leaf() {
    let node_leaves = node.get_leaves()
    
    let other_leaves : Array[String] = Array::new()
    let mut i = 0
    while i < all_taxa.length() {
      let found = contains_string(node_leaves, all_taxa[i])
      if !found {
        other_leaves.push(all_taxa[i])
      }
      i = i + 1
    }
    
    if node_leaves.length() > 0 && other_leaves.length() > 0 {
      let split1 = Split::new(node_leaves).normalize()
      result.push(split1)
    }
    
    let mut j = 0
    while j < node.children.length() {
      ignore(get_splits_recursive(node.children[j], all_taxa, result))
      j = j + 1
    }
  }
}

///|
fn contains_string(arr : Array[String], s : String) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == s {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn build_consensus(trees : Array[ConsensusTree], threshold : Double) -> ConsensusTree {
  if trees.length() == 0 {
    return ConsensusTree::{ nodes: Array::new(), root: None }
  }
  
  let all_taxa = get_all_taxa(trees[0])
  
  let split_counts : @hashmap.HashMap[String, Int] = @hashmap.HashMap([])
  
  let mut i = 0
  while i < trees.length() {
    match trees[i].root {
      Some(root) => {
        let splits = get_all_splits(root, all_taxa)
        let mut j = 0
        while j < splits.length() {
          let h = splits[j].hash()
          let current = if split_counts.contains(h) { split_counts.get(h).unwrap() } else { 0 }
          split_counts.set(h, current + 1)
          j = j + 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  
  let n_trees = trees.length().to_double()
  let consensus_splits : Array[Split] = Array::new()
  
  let entries = split_counts.iter().collect()
  let mut k = 0
  while k < entries.length() {
    let (hash, count) = entries[k]
    let support = count.to_double() / n_trees
    if support >= threshold {
      let taxa = split_by_char(hash, ','.to_int())
      let split = Split::new(taxa).set_support(support)
      consensus_splits.push(split)
    }
    k = k + 1
  }
  
  let consensus_root = build_tree_from_splits(all_taxa, consensus_splits)
  
  ConsensusTree::{
    nodes: Array::new(),
    root: Some(consensus_root)
  }
}

///|
pub fn build_majority_consensus(trees : Array[ConsensusTree]) -> ConsensusTree {
  build_consensus(trees, 0.5)
}

///|
pub fn build_strict_consensus(trees : Array[ConsensusTree]) -> ConsensusTree {
  build_consensus(trees, 1.0)
}

///|
pub fn build_half_consensus(trees : Array[ConsensusTree]) -> ConsensusTree {
  build_consensus(trees, 0.5)
}

///|
fn get_all_taxa(tree : ConsensusTree) -> Array[String] {
  match tree.root {
    Some(root) => root.get_leaves()
    None => Array::new()
  }
}

///|
fn build_tree_from_splits(taxa : Array[String], splits : Array[Split]) -> ConsensusNode {
  if splits.length() == 0 || taxa.length() == 1 {
    let root = ConsensusNode::new("root", true)
    let mut i = 0
    while i < taxa.length() {
      let leaf = ConsensusNode::new(taxa[i], false)
      root.children.push(leaf)
      i = i + 1
    }
    return root
  }
  
  let sorted_splits = sort_splits_by_support(splits)
  
  let root = ConsensusNode::new("root", true)
  let leaves : @hashmap.HashMap[String, ConsensusNode] = @hashmap.HashMap([])
  
  let mut i = 0
  while i < taxa.length() {
    leaves.set(taxa[i], ConsensusNode::new(taxa[i], false))
    i = i + 1
  }
  
  let mut j = 0
  while j < sorted_splits.length() {
    let split = sorted_splits[j]
    let split_taxa = split.taxa
    
    let group = Array::new()
    let remaining = Array::new()
    
    let keys = leaves.keys().collect()
    let mut k = 0
    while k < keys.length() {
      let taxon = keys[k]
      let node = leaves.get(taxon).unwrap()
      if contains_string(split_taxa, taxon) {
        group.push(node)
      } else {
        remaining.push(node)
      }
      k = k + 1
    }
    
    if group.length() > 1 && remaining.length() > 0 {
      let internal_node = ConsensusNode::new("", true).set_support(split.support)
      
      let mut l = 0
      while l < group.length() {
        internal_node.children.push(group[l])
        l = l + 1
      }
      
      leaves.clear()
      
      leaves.set("internal_" + j.to_string(), internal_node)
      
      let mut m = 0
      while m < remaining.length() {
        leaves.set(remaining[m].name, remaining[m])
        m = m + 1
      }
    }
    
    j = j + 1
  }
  
  let final_nodes = leaves.values().collect()
  let mut n = 0
  while n < final_nodes.length() {
    root.children.push(final_nodes[n])
    n = n + 1
  }
  
  root
}

///|
fn sort_splits_by_support(splits : Array[Split]) -> Array[Split] {
  let n = splits.length()
  let result = Array::make(n, Split::new(Array::new()))
  
  let mut i = 0
  while i < n {
    result[i] = splits[i]
    i = i + 1
  }
  
  i = 0
  while i < n - 1 {
    let mut j = i + 1
    while j < n {
      if result[i].support < result[j].support {
        let tmp = result[i]
        result[i] = result[j]
        result[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  
  result
}

///|
pub fn parse_consensus_newick(newick : String) -> ConsensusNode? {
  let trimmed = trim_newick(newick)
  if trimmed == "" {
    return None
  }
  
  let (root, _) = parse_newick_recursive(trimmed, 0)
  
  Some(root)
}

///|
fn trim_newick(s : String) -> String {
  let n = s.length()
  let mut start = 0
  while start < n && (s.unsafe_get(start).to_int() == ' '.to_int() || s.unsafe_get(start).to_int() == '\t'.to_int()) {
    start = start + 1
  }
  let mut end = n
  while end > start && (s.unsafe_get(end - 1).to_int() == ' '.to_int() || s.unsafe_get(end - 1).to_int() == '\t'.to_int() || s.unsafe_get(end - 1).to_int() == ';'.to_int()) {
    end = end - 1
  }
  substring(s, start, end - start)
}

///|
fn parse_newick_recursive(s : String, i : Int) -> (ConsensusNode, Int) {
  let children = Array::new()
  
  let mut idx = i + 1
  
  while idx < s.length() && s.unsafe_get(idx).to_int() != ')'.to_int() {
    if s.unsafe_get(idx).to_int() == '('.to_int() {
      let (child, new_idx) = parse_newick_recursive(s, idx)
      children.push(child)
      idx = new_idx
    } else if s.unsafe_get(idx).to_int() != ','.to_int() && s.unsafe_get(idx).to_int() != ')'.to_int() {
      let mut name = ""
      while idx < s.length() && s.unsafe_get(idx).to_int() != ','.to_int() && s.unsafe_get(idx).to_int() != ')'.to_int() && s.unsafe_get(idx).to_int() != ':'.to_int() {
        name = name + char_at(s, idx)
        idx = idx + 1
      }
      children.push(ConsensusNode::new(name, false))
      
      while idx < s.length() && (s.unsafe_get(idx).to_int() == ':'.to_int() || s.unsafe_get(idx).to_int() == ','.to_int() || s.unsafe_get(idx).to_int() == ')'.to_int()) {
        if s.unsafe_get(idx).to_int() == ')'.to_int() {
          break
        }
        idx = idx + 1
      }
    } else {
      idx = idx + 1
    }
  }
  
  idx = idx + 1
  
  let mut name = ""
  while idx < s.length() && s.unsafe_get(idx).to_int() != ';'.to_int() && s.unsafe_get(idx).to_int() != ','.to_int() && s.unsafe_get(idx).to_int() != ':'.to_int() {
    name = name + char_at(s, idx)
    idx = idx + 1
  }
  
  while idx < s.length() && (s.unsafe_get(idx).to_int() == ':'.to_int() || s.unsafe_get(idx).to_int() == ','.to_int() || s.unsafe_get(idx).to_int() == ')'.to_int()) {
    idx = idx + 1
  }
  
  let node = if children.length() > 0 {
    ConsensusNode::new(name, true)
  } else {
    ConsensusNode::new(name, false)
  }
  
  let mut j = 0
  while j < children.length() {
    node.children.push(children[j])
    j = j + 1
  }
  
  (node, idx)
}

///|
pub fn newick_to_tree(newick : String) -> ConsensusTree {
  match parse_consensus_newick(newick) {
    Some(root) => ConsensusTree::{ nodes: Array::new(), root: Some(root) }
    None => ConsensusTree::{ nodes: Array::new(), root: None }
  }
}

///|
pub fn ConsensusTree::to_newick(self : ConsensusTree) -> String {
  match self.root {
    Some(root) => node_to_newick(root) + ";"
    None => ""
  }
}

///|
fn node_to_newick(node : ConsensusNode) -> String {
  if node.is_leaf() {
    node.name
  } else {
    let mut result = "("
    
    let mut i = 0
    while i < node.children.length() {
      if i > 0 {
        result = result + ","
      }
      result = result + node_to_newick(node.children[i])
      i = i + 1
    }
    
    result = result + ")"
    
    if node.name != "" {
      result = result + node.name
    }
    
    if node.support > 0.0 {
      result = result + ":" + node.support.to_string()
    }
    
    result
  }
}

///|
pub fn create_example_trees() -> Array[ConsensusTree] {
  let tree1 = newick_to_tree("((A,B),(C,D));")
  let tree2 = newick_to_tree("((A,C),(B,D));")
  let tree3 = newick_to_tree("((A,B),(C,D));")
  let tree4 = newick_to_tree("((A,B,C),D);")
  let tree5 = newick_to_tree("((A,B),(C,D));")
  
  [tree1, tree2, tree3, tree4, tree5]
}

///|
pub fn create_simple_tree() -> ConsensusTree {
  newick_to_tree("((A,B),(C,D));")
}

///|
pub fn ConsensusNode::to_string(self : ConsensusNode) -> String {
  node_to_string(self, 0)
}

///|
fn node_to_string(node : ConsensusNode, indent : Int) -> String {
  let mut prefix = ""
  let mut i = 0
  while i < indent {
    prefix = prefix + "  "
    i = i + 1
  }
  
  let mut result = prefix
  
  if node.is_internal {
    result = result + "[internal"
    if node.support > 0.0 {
      result = result + " support=" + node.support.to_string()
    }
    if node.name != "" {
      result = result + " name=" + node.name
    }
    result = result + "]\n"
  } else {
    result = result + "[leaf] " + node.name + "\n"
  }
  
  let mut j = 0
  while j < node.children.length() {
    result = result + node_to_string(node.children[j], indent + 1)
    j = j + 1
  }
  
  result
}

///|
pub fn ConsensusTree::to_string(self : ConsensusTree) -> String {
  match self.root {
    Some(root) => root.to_string()
    None => "Empty tree\n"
  }
}

///|
pub fn calculate_consensus_support(trees : Array[ConsensusTree], split : Split) -> Double {
  if trees.length() == 0 {
    return 0.0
  }
  
  let all_taxa = get_all_taxa(trees[0])
  let normalized_split = split.normalize()
  let split_hash = normalized_split.hash()
  
  let mut count = 0
  let mut i = 0
  while i < trees.length() {
    match trees[i].root {
      Some(root) => {
        let splits = get_all_splits(root, all_taxa)
        let mut j = 0
        while j < splits.length() {
          if splits[j].normalize().hash() == split_hash {
            count = count + 1
            break
          }
          j = j + 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  
  count.to_double() / trees.length().to_double()
}

///|
pub fn get_all_splits_from_trees(trees : Array[ConsensusTree]) -> Array[Split] {
  if trees.length() == 0 {
    return Array::new()
  }
  
  let all_taxa = get_all_taxa(trees[0])
  let seen : @hashmap.HashMap[String, Split] = @hashmap.HashMap([])
  
  let mut i = 0
  while i < trees.length() {
    match trees[i].root {
      Some(root) => {
        let splits = get_all_splits(root, all_taxa)
        let mut j = 0
        while j < splits.length() {
          let normalized = splits[j].normalize()
          let h = normalized.hash()
          if !seen.contains(h) {
            seen.set(h, normalized)
          }
          j = j + 1
        }
      }
      None => ()
    }
    i = i + 1
  }
  
  seen.values().collect()
}