///|
/// GENIE3 module - Gene Network Inference with Ensemble of trees.
/// Infers gene regulatory networks from expression data using regression trees.
///
/// References:
///   - Huynh-Thu VA et al. (2010) Inferring regulatory networks from expression data
///     using tree-based methods. PLoS ONE 5(9): e12776.
///   - https://www.bioconductor.org/packages/release/bioc/html/GENIE3.html
///
/// Core algorithm:
/// 1. For each target gene, train a regression tree to predict its expression
///    from the expression of all other genes (candidate regulators).
/// 2. Use feature importance (variance reduction) of each predictor as edge weight.
/// 3. Combine weights across all trees to form a weighted adjacency matrix.
/// 4. Optionally symmetrize the matrix by averaging W[i][j] and W[j][i].

///|
/// Tree node: either an internal split node or a leaf with a constant prediction.
/// Stored by index in a shared array (avoids recursive type issues).
pub struct Genie3TreeNode {
  // For internal nodes:
  feature : Int              // index of splitting feature (-1 for leaf)
  threshold : Double         // split threshold (feature <= threshold goes left)
  left_idx : Int             // index into shared node array (-1 for leaf)
  right_idx : Int            // index into shared node array (-1 for leaf)
  // For leaf nodes:
  value : Double             // leaf prediction value
  is_leaf : Bool
  importance_gain : Double   // variance reduction at this split
}

///|
/// Create a leaf node.
pub fn genie3_tree_leaf(value : Double) -> Genie3TreeNode {
  Genie3TreeNode::{
    feature: -1,
    threshold: 0.0,
    left_idx: -1,
    right_idx: -1,
    value,
    is_leaf: true,
    importance_gain: 0.0,
  }
}

///|
/// Create an internal split node.
pub fn genie3_tree_split(
  feature : Int,
  threshold : Double,
  left_idx : Int,
  right_idx : Int,
  importance_gain : Double,
) -> Genie3TreeNode {
  Genie3TreeNode::{
    feature,
    threshold,
    left_idx,
    right_idx,
    value: 0.0,
    is_leaf: false,
    importance_gain,
  }
}

///|
/// A regression tree stored as an array of nodes (index 0 is the root).
pub struct RegressionTree {
  nodes : Array[Genie3TreeNode]
}

///|
/// Predict a value for a single sample by traversing the tree.
pub fn RegressionTree::predict(self : RegressionTree, sample : Array[Double]) -> Double {
  let mut node_idx = 0
  while true {
    let node = self.nodes[node_idx]
    if node.is_leaf {
      return node.value
    }
    if sample[node.feature] <= node.threshold {
      node_idx = node.left_idx
    } else {
      node_idx = node.right_idx
    }
  }
  0.0
}

///|
/// Recursively collect feature importance by summing importance_gain of internal nodes.
pub fn RegressionTree::feature_importance(
  self : RegressionTree,
  n_features : Int,
) -> Array[Double] {
  let importance = Array::make(n_features, 0.0)
  let stack : Array[Int] = [0]
  while stack.length() > 0 {
    let node_idx = stack.pop().unwrap()
    let node = self.nodes[node_idx]
    if not(node.is_leaf) {
      importance[node.feature] = importance[node.feature] + node.importance_gain
      stack.push(node.left_idx)
      stack.push(node.right_idx)
    }
  }
  importance
}

///|
/// Compute the mean of an array of doubles.
fn genie3_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < values.length() {
    sum = sum + values[i]
    i = i + 1
  }
  sum / values.length().to_double()
}

///|
/// Compute sum of squared deviations from the mean (total SSE for a node).
fn genie3_sse(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mu = genie3_mean(values)
  let mut s = 0.0
  let mut i = 0
  while i < values.length() {
    let d = values[i] - mu
    s = s + d * d
    i = i + 1
  }
  s
}

///|
/// Sort a list of (feature_value, target_value) pairs by feature_value ascending.
/// Uses a simple insertion sort (datasets for GENIE3 are typically small per node).
fn genie3_sort_pairs(
  pairs : Array[(Double, Double)],
) -> Array[(Double, Double)] {
  let result = Array::make(pairs.length(), (0.0, 0.0))
  let mut i = 0
  while i < pairs.length() {
    result[i] = pairs[i]
    i = i + 1
  }
  // Insertion sort on first element of tuple
  let mut j = 1
  while j < result.length() {
    let (key_first, key_second) = result[j]
    let key = (key_first, key_second)
    let mut k = j - 1
    while k >= 0 {
      let (curr_first, _) = result[k]
      if curr_first <= key_first {
        break
      }
      result[k + 1] = result[k]
      k = k - 1
    }
    result[k + 1] = key
    j = j + 1
  }
  result
}

///|
/// Internal recursive helper that builds the tree by appending nodes to `tree.nodes`.
/// Returns the index of the root node of the subtree.
fn genie3_build_recursive(
  tree : RegressionTree,
  features_matrix : Array[Array[Double]],
  target : Array[Double],
  feature_indices : Array[Int],
  max_depth : Int,
  min_samples_split : Int,
  depth : Int,
) -> Int {
  let n_samples = target.length()

  // Stopping conditions: leaf prediction is the mean of targets
  if n_samples < min_samples_split || depth >= max_depth || n_samples < 2 {
    let leaf = genie3_tree_leaf(genie3_mean(target))
    tree.nodes.push(leaf)
    return tree.nodes.length() - 1
  }

  // Pre-compute current SSE
  let parent_sse = genie3_sse(target)
  if parent_sse <= 1.0e-12 {
    let leaf = genie3_tree_leaf(genie3_mean(target))
    tree.nodes.push(leaf)
    return tree.nodes.length() - 1
  }

  // Find best split across features
  let mut best_gain = -1.0
  let mut best_feature = -1
  let mut best_threshold = 0.0
  let mut best_left_idx : Array[Int] = []
  let mut best_right_idx : Array[Int] = []

  let mut fi = 0
  while fi < feature_indices.length() {
    let f = feature_indices[fi]
    // Build (feature_value, target_value) pairs and sort
    let pairs : Array[(Double, Double)] = []
    let mut si = 0
    while si < n_samples {
      pairs.push((features_matrix[si][f], target[si]))
      si = si + 1
    }
    let sorted = genie3_sort_pairs(pairs)
    // Try midpoints between consecutive distinct feature values
    let n_sorted = sorted.length()
    let mut k = 0
    while k < n_sorted - 1 {
      let (sk0, _) = sorted[k]
      let (sk1, _) = sorted[k + 1]
      if sk0 < sk1 {
        let threshold = (sk0 + sk1) / 2.0
        let left_targets : Array[Double] = []
        let right_targets : Array[Double] = []
        let left_idx : Array[Int] = []
        let right_idx : Array[Int] = []
        let mut si = 0
        while si < n_samples {
          if features_matrix[si][f] <= threshold {
            left_targets.push(target[si])
            left_idx.push(si)
          } else {
            right_targets.push(target[si])
            right_idx.push(si)
          }
          si = si + 1
        }
        if left_targets.length() > 0 && right_targets.length() > 0 {
          let left_sse = genie3_sse(left_targets)
          let right_sse = genie3_sse(right_targets)
          let n_left = left_targets.length().to_double()
          let n_right = right_targets.length().to_double()
          let n_total = n_samples.to_double()
          let weighted_sse = (n_left * left_sse + n_right * right_sse) / n_total
          let gain = parent_sse - weighted_sse
          if gain > best_gain {
            best_gain = gain
            best_feature = f
            best_threshold = threshold
            best_left_idx = left_idx
            best_right_idx = right_idx
          }
        }
      }
      k = k + 1
    }
    fi = fi + 1
  }

  // No improving split found
  if best_feature < 0 || best_gain <= 0.0 {
    let leaf = genie3_tree_leaf(genie3_mean(target))
    tree.nodes.push(leaf)
    return tree.nodes.length() - 1
  }

  // Reserve slot for this node (we'll fill it after computing children's indices)
  let placeholder = genie3_tree_leaf(0.0)
  tree.nodes.push(placeholder)
  let current_idx = tree.nodes.length() - 1

  // Build child datasets
  let left_features : Array[Array[Double]] = []
  let left_target : Array[Double] = []
  let mut li = 0
  while li < best_left_idx.length() {
    let idx = best_left_idx[li]
    left_features.push(features_matrix[idx])
    left_target.push(target[idx])
    li = li + 1
  }
  let right_features : Array[Array[Double]] = []
  let right_target : Array[Double] = []
  let mut ri = 0
  while ri < best_right_idx.length() {
    let idx = best_right_idx[ri]
    right_features.push(features_matrix[idx])
    right_target.push(target[idx])
    ri = ri + 1
  }

  let left_idx = genie3_build_recursive(
    tree,
    left_features,
    left_target,
    feature_indices,
    max_depth,
    min_samples_split,
    depth + 1,
  )
  let right_idx = genie3_build_recursive(
    tree,
    right_features,
    right_target,
    feature_indices,
    max_depth,
    min_samples_split,
    depth + 1,
  )

  // Now fill in the placeholder slot
  tree.nodes[current_idx] = genie3_tree_split(
    best_feature,
    best_threshold,
    left_idx,
    right_idx,
    best_gain,
  )
  current_idx
}

///|
/// Build a regression tree predicting `target` from `features_matrix` (samples × features).
///
/// Parameters:
///   features_matrix : n_samples × n_features matrix of predictor values
///   target          : n_samples vector of response values
///   feature_indices : indices of features usable as splitters (these become the
///                     feature identifiers used in feature_importance output)
///   max_depth       : maximum tree depth
///   min_samples_split : minimum samples required to split a node
pub fn genie3_build_tree(
  features_matrix : Array[Array[Double]],
  target : Array[Double],
  feature_indices : Array[Int],
  max_depth? : Int = 10,
  min_samples_split? : Int = 2,
) -> RegressionTree {
  let tree : RegressionTree = { nodes: [] }
  let _ = genie3_build_recursive(
    tree,
    features_matrix,
    target,
    feature_indices,
    max_depth,
    min_samples_split,
    0,
  )
  tree
}

///|
/// A weighted directed edge in the inferred regulatory network.
pub struct RegulatoryEdge {
  regulator : String    // source gene (predictor)
  target : String       // target gene
  weight : Double       // importance weight
}

///|
/// GENIE3 inference result: weighted adjacency matrix and edges.
pub struct Genie3Result {
  gene_names : Array[String]
  matrix : Array[Array[Double]]    // matrix[i][j] = importance of gene i for gene j
  edges : Array[RegulatoryEdge]    // sorted edges (descending weight)
  n_trees : Int
  symmetrized : Bool
}

///|
/// Run GENIE3 to infer a gene regulatory network from expression data.
///
/// Parameters:
///   expression  : n_samples × n_genes matrix
///   gene_names  : names of genes (column order matches expression)
///   max_depth   : maximum depth of each regression tree
///   min_samples_split : minimum samples required to split a node
///   symmetrize  : whether to symmetrize the matrix by averaging W[i][j] and W[j][i]
pub fn genie3_run(
  expression : Array[Array[Double]],
  gene_names : Array[String],
  max_depth? : Int = 10,
  min_samples_split? : Int = 2,
  symmetrize? : Bool = false,
) -> Genie3Result {
  let n_samples = expression.length()
  let n_genes = gene_names.length()
  let md = max_depth
  let mss = min_samples_split
  let sym = symmetrize

  // Initialize weight matrix W[i][j] = importance of gene i for predicting gene j
  let matrix : Array[Array[Double]] = Array::new()
  let mut mi = 0
  while mi < n_genes {
    matrix.push(Array::make(n_genes, 0.0))
    mi = mi + 1
  }

  // For each target gene j, build a tree using all other genes as predictors.
  let mut j = 0
  while j < n_genes {
    // Target column
    let target : Array[Double] = Array::make(n_samples, 0.0)
    // Predictor matrix: exclude gene j
    let predictor_indices : Array[Int] = []
    let mut gi = 0
    while gi < n_genes {
      if gi != j {
        predictor_indices.push(gi)
      }
      gi = gi + 1
    }
    let predictor_matrix : Array[Array[Double]] = Array::new()
    let mut pm = 0
    while pm < n_samples {
      predictor_matrix.push(Array::make(predictor_indices.length(), 0.0))
      pm = pm + 1
    }
    let mut si = 0
    while si < n_samples {
      target[si] = expression[si][j]
      let mut pi = 0
      while pi < predictor_indices.length() {
        predictor_matrix[si][pi] = expression[si][predictor_indices[pi]]
        pi = pi + 1
      }
      si = si + 1
    }
    // Use feature indices 0..n_predictors-1 (column indices into predictor_matrix)
    let local_indices : Array[Int] = []
    let mut li = 0
    while li < predictor_indices.length() {
      local_indices.push(li)
      li = li + 1
    }
    let tree = genie3_build_tree(
      predictor_matrix,
      target,
      local_indices,
      max_depth=md,
      min_samples_split=mss,
    )
    let importance = tree.feature_importance(predictor_indices.length())
    // Map back: predictor_matrix column k corresponds to global gene predictor_indices[k]
    let mut k = 0
    while k < importance.length() {
      let global_gene = predictor_indices[k]
      matrix[global_gene][j] = matrix[global_gene][j] + importance[k]
      k = k + 1
    }
    j = j + 1
  }

  // Normalize each column so importances sum to 1 (avoid scale differences across trees)
  let mut cj = 0
  while cj < n_genes {
    let mut sum = 0.0
    let mut ci = 0
    while ci < n_genes {
      sum = sum + matrix[ci][cj]
      ci = ci + 1
    }
    if sum > 1.0e-12 {
      let mut ci = 0
      while ci < n_genes {
        matrix[ci][cj] = matrix[ci][cj] / sum
        ci = ci + 1
      }
    }
    cj = cj + 1
  }

  // Optionally symmetrize
  if sym {
    let mut i = 0
    while i < n_genes {
      let mut j = i + 1
      while j < n_genes {
        let avg = (matrix[i][j] + matrix[j][i]) / 2.0
        matrix[i][j] = avg
        matrix[j][i] = avg
        j = j + 1
      }
      i = i + 1
    }
  }

  // Collect edges (skip self-loops and zero weights)
  let edges : Array[RegulatoryEdge] = []
  let mut i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_genes {
      if i != j && matrix[i][j] > 0.0 {
        edges.push(
          RegulatoryEdge::{ regulator: gene_names[i], target: gene_names[j], weight: matrix[i][j] },
        )
      }
      j = j + 1
    }
    i = i + 1
  }
  // Sort edges by weight descending (selection sort)
  let mut ei = 0
  while ei < edges.length() {
    let mut ej = ei + 1
    while ej < edges.length() {
      if edges[ej].weight > edges[ei].weight {
        let tmp = edges[ei]
        edges[ei] = edges[ej]
        edges[ej] = tmp
      }
      ej = ej + 1
    }
    ei = ei + 1
  }

  Genie3Result::{ gene_names, matrix, edges, n_trees: n_genes, symmetrized: sym }
}

///|
/// Get the top-K regulatory edges by weight.
pub fn genie3_top_edges(result : Genie3Result, k : Int) -> Array[RegulatoryEdge] {
  let n = if k < result.edges.length() { k } else { result.edges.length() }
  let top : Array[RegulatoryEdge] = []
  let mut i = 0
  while i < n {
    top.push(result.edges[i])
    i = i + 1
  }
  top
}

///|
/// Get regulators of a target gene (edges pointing to it).
pub fn genie3_regulators_of(
  result : Genie3Result,
  target_gene : String,
) -> Array[RegulatoryEdge] {
  let result_arr : Array[RegulatoryEdge] = []
  let mut i = 0
  while i < result.edges.length() {
    if result.edges[i].target == target_gene {
      result_arr.push(result.edges[i])
    }
    i = i + 1
  }
  result_arr
}

///|
/// Get targets of a regulator gene (edges originating from it).
pub fn genie3_targets_of(
  result : Genie3Result,
  regulator_gene : String,
) -> Array[RegulatoryEdge] {
  let result_arr : Array[RegulatoryEdge] = []
  let mut i = 0
  while i < result.edges.length() {
    if result.edges[i].regulator == regulator_gene {
      result_arr.push(result.edges[i])
    }
    i = i + 1
  }
  result_arr
}

///|
/// Generate a small synthetic dataset for demonstration and testing.
/// Returns (expression, gene_names) where a 5-gene network has clear regulators.
pub fn genie3_sample_data() -> (Array[Array[Double]], Array[String]) {
  // 5 genes: G1 regulates G2 and G3; G2 regulates G4; G3 regulates G5
  let gene_names = ["G1", "G2", "G3", "G4", "G5"]
  let expression : Array[Array[Double]] = []
  let mut s = 0
  // Deterministic pseudo-random for reproducibility
  let mut state = 42
  while s < 30 {
    state = (state * 1103515245 + 12345) & 0x7fffffff
    let g1 = (state % 1000).to_double() / 100.0
    state = (state * 1103515245 + 12345) & 0x7fffffff
    let noise2 = (state % 100).to_double() / 50.0 - 1.0
    let g2 = g1 * 0.8 + noise2
    state = (state * 1103515245 + 12345) & 0x7fffffff
    let noise3 = (state % 100).to_double() / 50.0 - 1.0
    let g3 = g1 * 0.5 + noise3
    state = (state * 1103515245 + 12345) & 0x7fffffff
    let noise4 = (state % 100).to_double() / 50.0 - 1.0
    let g4 = g2 * 0.7 + noise4
    state = (state * 1103515245 + 12345) & 0x7fffffff
    let noise5 = (state % 100).to_double() / 50.0 - 1.0
    let g5 = g3 * 0.3 + noise5
    expression.push([g1, g2, g3, g4, g5])
    s = s + 1
  }
  (expression, gene_names)
}
