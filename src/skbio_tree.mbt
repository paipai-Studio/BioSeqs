///|
/// scikit-bio `skbio.tree`（TreeNode 核心）+ 基于树的 UniFrac 距离。
///
/// 这里实现的是 UniFrac 所需的 TreeNode 最小子集：带分支长度的有根树 +
/// Newick 解析 + tips 访问，以及 unweighted/weighted UniFrac。语义对照
/// scikit-bio 0.7.3（`tree/_tree.py` `to_array` 的分支长度约定 + `diversity/
/// beta/_unifrac.py` 的算法）。

///|
/// A node in a rooted tree with branch lengths. `length` is the length of the
/// edge connecting this node to its **parent** (the `name:length` value in
/// Newick); the root's length is 0. `name` is `None` for unnamed internal
/// nodes.
pub struct SkbioNode {
  name : String?
  length : Double
  children : Array[SkbioNode]
}

///|
pub fn SkbioNode::new(name : String?, length : Double, children : Array[SkbioNode]) -> SkbioNode {
  SkbioNode::{ name, length, children }
}

///|
/// A tip (leaf) is a node with no children.
pub fn SkbioNode::is_tip(self : SkbioNode) -> Bool {
  self.children.length() == 0
}

///|
fn skbio_collect_tips(node : SkbioNode, out : Array[String]) -> Unit {
  if node.is_tip() {
    match node.name {
      Some(n) => out.push(n)
      None => ()
    }
  } else {
    for child in node.children {
      skbio_collect_tips(child, out)
    }
  }
}

///|
/// Names of all tip (leaf) nodes, in preorder.
pub fn SkbioNode::tips(self : SkbioNode) -> Array[String] {
  let out : Array[String] = Array::new()
  skbio_collect_tips(self, out)
  out
}

///|
/// Parse `name[:length]` starting at `i`; returns (name, length, next_index).
fn skbio_parse_label(s : String, i : Int) -> (String?, Double, Int) {
  let n = s.length()
  let mut j = i
  let mut name = ""
  while j < n {
    let c = char_at(s, j)
    if c == ":" || c == "," || c == ")" || c == "(" {
      break
    }
    name = name + c
    j = j + 1
  }
  let name_opt : String? = if name == "" { None } else { Some(name) }
  let mut length = 0.0
  if j < n && char_at(s, j) == ":" {
    j = j + 1
    let mut num = ""
    while j < n {
      let c = char_at(s, j)
      if c == ":" || c == "," || c == ")" || c == "(" {
        break
      }
      num = num + c
      j = j + 1
    }
    length = match parse_double(num) {
      Some(v) => v
      None => 0.0
    }
  }
  (name_opt, length, j)
}

///|
/// Recursive-descent Newick node parser; returns (node, next_index).
fn skbio_parse_node(s : String, i : Int) -> (SkbioNode, Int) {
  if char_at(s, i) == "(" {
    let mut j = i + 1
    let children : Array[SkbioNode] = Array::new()
    let (c, nj) = skbio_parse_node(s, j)
    children.push(c)
    j = nj
    while j < s.length() && char_at(s, j) == "," {
      let (c2, nj2) = skbio_parse_node(s, j + 1)
      children.push(c2)
      j = nj2
    }
    // s[j] is now ')'
    let after_close = j + 1
    let (name, length, nj3) = skbio_parse_label(s, after_close)
    (SkbioNode::{ name, length, children }, nj3)
  } else {
    let (name, length, nj) = skbio_parse_label(s, i)
    (SkbioNode::{ name, length, children: [] }, nj)
  }
}

///|
/// Read a rooted Newick string into a `SkbioNode` tree (the `TreeNode.read`
/// subset). Supports nested parentheses, comma-separated children, tip and
/// internal-node labels, and `:branch-length` on every node.
pub fn skbio_tree_read(newick : String) -> SkbioNode {
  let s = newick.trim().replace(old=";", new="").to_owned()
  if s.length() == 0 {
    return SkbioNode::{ name: None, length: 0.0, children: [] }
  }
  let (node, _) = skbio_parse_node(s, 0)
  node
}

///|
fn skbio_taxon_index(taxa : Array[String], name : String) -> Int {
  for i = 0; i < taxa.length(); i = i + 1 {
    if taxa[i] == name {
      return i
    }
  }
  -1
}

///|
/// Postorder accumulation: for every node, count of `u`/`v` reads descending
/// from (i.e. in the subtree below) that node. Appends (length, cu, cv) per
/// node into the three out arrays.
fn skbio_unifrac_accumulate(
  node : SkbioNode,
  taxa : Array[String],
  u_counts : Array[Double],
  v_counts : Array[Double],
  lengths : Array[Double],
  cu : Array[Double],
  cv : Array[Double],
) -> (Double, Double) {
  let mut u = 0.0
  let mut v = 0.0
  if node.is_tip() {
    match node.name {
      None => ()
      Some(n) => {
        let idx = skbio_taxon_index(taxa, n)
        if idx >= 0 {
          u = u_counts[idx]
          v = v_counts[idx]
        }
      }
    }
  }
  for child in node.children {
    let (cu_sub, cv_sub) = skbio_unifrac_accumulate(child, taxa, u_counts, v_counts, lengths, cu, cv)
    u = u + cu_sub
    v = v + cv_sub
  }
  lengths.push(node.length)
  cu.push(u)
  cv.push(v)
  (u, v)
}

///|
fn skbio_unifrac_arrays(
  tree : SkbioNode,
  taxa : Array[String],
  u_counts : Array[Double],
  v_counts : Array[Double],
) -> (Array[Double], Array[Double], Array[Double]) {
  let lengths : Array[Double] = Array::new()
  let cu : Array[Double] = Array::new()
  let cv : Array[Double] = Array::new()
  let (_, _) = skbio_unifrac_accumulate(tree, taxa, u_counts, v_counts, lengths, cu, cv)
  (lengths, cu, cv)
}

///|
/// Unweighted UniFrac distance between two samples (`unweighted_unifrac`):
/// Σ len·XOR(present_u, present_v) / Σ len·OR(present_u, present_v), where
/// present = node count > 0, summed over all nodes (root contributes len 0).
pub fn skbio_unweighted_unifrac(
  u_counts : Array[Double],
  v_counts : Array[Double],
  taxa : Array[String],
  tree : SkbioNode,
) -> Double {
  let (lengths, cu, cv) = skbio_unifrac_arrays(tree, taxa, u_counts, v_counts)
  let mut unique = 0.0
  let mut observed = 0.0
  for i = 0; i < lengths.length(); i = i + 1 {
    let pu = cu[i] > 0.0
    let pv = cv[i] > 0.0
    if pu != pv {
      unique = unique + lengths[i]
    }
    if pu || pv {
      observed = observed + lengths[i]
    }
  }
  if observed == 0.0 {
    0.0
  } else {
    unique / observed
  }
}

///|
/// Weighted UniFrac distance (unnormalized, `weighted_unifrac` default):
/// Σ len · |prop_u − prop_v|, where prop = node count / total reads.
pub fn skbio_weighted_unifrac(
  u_counts : Array[Double],
  v_counts : Array[Double],
  taxa : Array[String],
  tree : SkbioNode,
) -> Double {
  let (lengths, cu, cv) = skbio_unifrac_arrays(tree, taxa, u_counts, v_counts)
  let mut tu = 0.0
  for x in u_counts {
    tu = tu + x
  }
  let mut tv = 0.0
  for x in v_counts {
    tv = tv + x
  }
  let mut wu = 0.0
  for i = 0; i < lengths.length(); i = i + 1 {
    let pu = if tu > 0.0 { cu[i] / tu } else { cu[i] }
    let pv = if tv > 0.0 { cv[i] / tv } else { cv[i] }
    wu = wu + lengths[i] * (pu - pv).abs()
  }
  wu
}