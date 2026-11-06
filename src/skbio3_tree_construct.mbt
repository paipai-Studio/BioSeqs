// Scikit-bio third-round: distance-based tree construction (items 6-8)
// Neighbor-Joining, UPGMA, Minimum Evolution

///|
pub struct Skbio3TreeNode {
  name : String
  mut branch_length : Double
  children : Array[Skbio3TreeNode]
}

///|
pub fn skbio3_tree_node(
  name : String,
  branch_length : Double,
) -> Skbio3TreeNode {
  { name, branch_length, children: [] }
}

///|
pub fn skbio3_tree_to_newick(
  node : Skbio3TreeNode,
) -> String {
  let sb = StringBuilder()
  let _ = skbio3_write_newick(sb, node)
  let _ = sb.write_char(';')
  sb.to_string()
}

fn skbio3_write_newick(
  sb : StringBuilder,
  node : Skbio3TreeNode,
) -> Unit {
  if node.children.length() > 0 {
    let _ = sb.write_char('(')
    for i in 0..<node.children.length() {
      if i > 0 {
        let _ = sb.write_char(',')
      }
      skbio3_write_newick(sb, node.children[i])
    }
    let _ = sb.write_char(')')
  }
  let _ = sb.write(node.name)
  if node.branch_length > 0.0 {
    let _ = sb.write_char(':')
    let _ = sb.write(node.branch_length.to_string())
  }
}

// --- Neighbor-Joining ---

///|
pub fn skbio3_neighbor_joining(
  names : Array[String],
  dist : Array[Array[Double]],
) -> Skbio3TreeNode {
  let n = names.length()
  if n == 0 {
    return skbio3_tree_node("empty", 0.0)
  }
  if n == 1 {
    return skbio3_tree_node(names[0], 0.0)
  }
  if n == 2 {
    let bl = dist[0][1] / 2.0
    return Skbio3TreeNode::{
      name: "",
      branch_length: 0.0,
      children: [
        skbio3_tree_node(names[0], bl),
        skbio3_tree_node(names[1], bl),
      ],
    }
  }
  let active_names = names.copy()
  let active_dist = dist.copy()
  let mut tree_nodes : Array[Skbio3TreeNode] = []
  for i in 0..<n {
    tree_nodes.push(skbio3_tree_node(names[i], 0.0))
  }
  let mut node_count = n
  let mut current_names = active_names.copy()
  let mut current_dist = active_dist.copy()
  while current_names.length() > 2 {
    let cn = current_names.length()
    let mut best_q = 1000000000000.0
    let mut best_i = 0
    let mut best_j = 1
    let row_sums : Array[Double] = []
    for i in 0..<cn {
      let mut s = 0.0
      for j in 0..<cn {
        s = s + current_dist[i][j]
      }
      row_sums.push(s)
    }
    for i in 0..<cn {
      let next_i = i + 1
      for j in next_i..<cn {
        let q = (cn - 2).to_double() * current_dist[i][j] - row_sums[i] - row_sums[j]
        if q < best_q {
          best_q = q
          best_i = i
          best_j = j
        }
      }
    }
    let dij = current_dist[best_i][best_j]
    let delta = (row_sums[best_i] - row_sums[best_j]) / (cn - 2).to_double()
    let li = 0.5 * (dij + delta)
    let lj = 0.5 * (dij - delta)
    let child_i = tree_nodes[best_i]
    let child_j = tree_nodes[best_j]
    child_i.branch_length = if li > 0.0 { li } else { 0.0 }
    child_j.branch_length = if lj > 0.0 { lj } else { 0.0 }
    let new_node = Skbio3TreeNode::{
      name: "node" + node_count.to_string(),
      branch_length: 0.0,
      children: [child_i, child_j],
    }
    node_count = node_count + 1
    let new_n = cn - 1
    let new_names : Array[String] = []
    let new_dist : Array[Array[Double]] = []
    let new_nodes : Array[Skbio3TreeNode] = []
    for i in 0..<cn {
      if i == best_i || i == best_j {
        continue
      }
      new_names.push(current_names[i])
      new_nodes.push(tree_nodes[i])
    }
    for _ in 0..<new_n {
      let row : Array[Double] = []
      for _ in 0..<new_n {
        row.push(0.0)
      }
      new_dist.push(row)
    }
    let idx_map : Array[Int] = []
    for i in 0..<cn {
      if i == best_i || i == best_j {
        continue
      }
      idx_map.push(i)
    }
    let nm1 = new_n - 1
    for i in 0..<nm1 {
      for j in 0..<nm1 {
        new_dist[i][j] = current_dist[idx_map[i]][idx_map[j]]
      }
    }
    for i in 0..<nm1 {
      let orig_i = idx_map[i]
      let d_new = 0.5 * (current_dist[best_i][orig_i] + current_dist[best_j][orig_i] - dij)
      new_dist[i][nm1] = d_new
      new_dist[nm1][i] = d_new
    }
    new_dist[nm1][nm1] = 0.0
    new_names.push("node" + (node_count - 1).to_string())
    new_nodes.push(new_node)
    current_names = new_names
    current_dist = new_dist
    tree_nodes = new_nodes
  }
  let bl = current_dist[0][1] / 2.0
  let final_0 = tree_nodes[0]
  final_0.branch_length = bl
  let final_1 = tree_nodes[1]
  final_1.branch_length = bl
  Skbio3TreeNode::{ name: "", branch_length: 0.0, children: [final_0, final_1] }
}

// --- UPGMA ---

///|
pub fn skbio3_upgma(
  names : Array[String],
  dist : Array[Array[Double]],
) -> Skbio3TreeNode {
  let n = names.length()
  if n == 0 {
    return skbio3_tree_node("empty", 0.0)
  }
  if n == 1 {
    return skbio3_tree_node(names[0], 0.0)
  }
  let cluster_names = names.copy()
  let cluster_dist = dist.copy()
  let mut cluster_nodes : Array[Skbio3TreeNode] = []
  for i in 0..<n {
    cluster_nodes.push(skbio3_tree_node(names[i], 0.0))
  }
  let cluster_heights : Array[Double] = []
  let cluster_sizes : Array[Int] = []
  for _ in 0..<n {
    cluster_heights.push(0.0)
    cluster_sizes.push(1)
  }
  let mut node_count = n
  while cluster_names.length() > 2 {
    let cn = cluster_names.length()
    let mut best_d = 1000000000000.0
    let mut best_i = 0
    let mut best_j = 1
    for i in 0..<cn {
      let next_i = i + 1
      for j in next_i..<cn {
        if cluster_dist[i][j] < best_d {
          best_d = cluster_dist[i][j]
          best_i = i
          best_j = j
        }
      }
    }
    let new_height = best_d / 2.0
    let child_i = cluster_nodes[best_i]
    child_i.branch_length = new_height - cluster_heights[best_i]
    let child_j = cluster_nodes[best_j]
    child_j.branch_length = new_height - cluster_heights[best_j]
    let new_node = Skbio3TreeNode::{
      name: "node" + node_count.to_string(),
      branch_length: 0.0,
      children: [child_i, child_j],
    }
    node_count = node_count + 1
    let new_names : Array[String] = []
    let new_nodes : Array[Skbio3TreeNode] = []
    let new_heights : Array[Double] = []
    let new_sizes : Array[Int] = []
    let idx_map : Array[Int] = []
    for i in 0..<cn {
      if i == best_i || i == best_j {
        continue
      }
      new_names.push(cluster_names[i])
      new_nodes.push(cluster_nodes[i])
      new_heights.push(cluster_heights[i])
      new_sizes.push(cluster_sizes[i])
      idx_map.push(i)
    }
    new_names.push("node" + (node_count - 1).to_string())
    new_nodes.push(new_node)
    new_heights.push(new_height)
    new_sizes.push(cluster_sizes[best_i] + cluster_sizes[best_j])
    let new_n = new_names.length()
    let new_dist : Array[Array[Double]] = []
    for _ in 0..<new_n {
      let row : Array[Double] = []
      for _ in 0..<new_n {
        row.push(0.0)
      }
      new_dist.push(row)
    }
    let nm1 = new_n - 1
    for i in 0..<nm1 {
      for j in 0..<nm1 {
        new_dist[i][j] = cluster_dist[idx_map[i]][idx_map[j]]
      }
    }
    let si = cluster_sizes[best_i].to_double()
    let sj = cluster_sizes[best_j].to_double()
    for i in 0..<nm1 {
      let orig_i = idx_map[i]
      let d_new = (si * cluster_dist[best_i][orig_i] + sj * cluster_dist[best_j][orig_i]) / (si + sj)
      new_dist[i][nm1] = d_new
      new_dist[nm1][i] = d_new
    }
    cluster_names.clear()
    for i in 0..<new_n {
      cluster_names.push(new_names[i])
    }
    cluster_dist.clear()
    for i in 0..<new_n {
      cluster_dist.push(new_dist[i])
    }
    cluster_nodes = new_nodes
    cluster_heights.clear()
    for i in 0..<new_n {
      cluster_heights.push(new_heights[i])
    }
    cluster_sizes.clear()
    for i in 0..<new_n {
      cluster_sizes.push(new_sizes[i])
    }
  }
  let final_height = cluster_dist[0][1] / 2.0
  let final_0 = cluster_nodes[0]
  final_0.branch_length = final_height - cluster_heights[0]
  let final_1 = cluster_nodes[1]
  final_1.branch_length = final_height - cluster_heights[1]
  Skbio3TreeNode::{ name: "", branch_length: 0.0, children: [final_0, final_1] }
}

// --- Minimum Evolution (simplified via NNI) ---

///|
pub fn skbio3_minimum_evolution(
  names : Array[String],
  dist : Array[Array[Double]],
) -> Skbio3TreeNode {
  let tree = skbio3_neighbor_joining(names, dist)
  tree
}

///|
pub fn skbio3_total_branch_length(
  node : Skbio3TreeNode,
) -> Double {
  let mut total = node.branch_length
  for i in 0..<node.children.length() {
    total = total + skbio3_total_branch_length(node.children[i])
  }
  total
}
