///|
/// TreeConstruction - Distance-based phylogenetic tree construction.
/// Inspired by Biopython Bio.Phylo.TreeConstruction.
/// Implements UPGMA, WPGMA, and Neighbor-Joining algorithms.

///|
/// Represents a node in the tree construction process.
struct BuildNode {
  id : Int
  name : String
  height : Double
  children : Array[BuildNode]
  branch_lengths : Array[Double]
}

///|
/// Create a leaf build node.
fn build_node_leaf(id : Int, name : String) -> BuildNode {
  BuildNode::{ id, name, height: 0.0, children: [], branch_lengths: [] }
}

///|
/// Create an internal build node from two children.
fn build_node_internal(
  id : Int,
  left : BuildNode,
  right : BuildNode,
  left_branch : Double,
  right_branch : Double,
  height : Double,
) -> BuildNode {
  BuildNode::{
    id,
    name: "",
    height,
    children: [left, right],
    branch_lengths: [left_branch, right_branch],
  }
}

///|
/// Convert a BuildNode to a Clade (phylogenetic tree node).
fn build_node_to_clade(node : BuildNode) -> Clade {
  if node.children.length() == 0 {
    Clade::new(name=Some(node.name))
  } else {
    let mut clade = Clade::new(
      name=if node.name.length() > 0 { Some(node.name) } else { None },
    )
    let mut i = 0
    while i < node.children.length() {
      let child_clade = build_node_to_clade(node.children[i])
      let bl = node.branch_lengths[i]
      let child_with_bl = Clade::{
        branch_length: Some(bl),
        name: child_clade.name,
        confidence: child_clade.confidence,
        clades: child_clade.clades,
        parent: None,
      }
      clade = clade.add_child(child_with_bl)
      i = i + 1
    }
    clade
  }
}

///|
/// Internal distance model kinds.
enum DistanceModelKind {
  /// Biopython "identity" model: plain character identity.
  Identity
  /// Biopython substitution scoring matrix models (blosum62, blastn, ...).
  Matrix
  /// Jukes-Cantor correction (MoonBit extension, not in Biopython 1.88).
  Jc69
  /// Kimura 2-parameter correction (MoonBit extension, not in Biopython 1.88).
  K80
}

///|
/// Distance calculator for building distance matrices from sequences.
///
/// Port of Biopython `Bio.Phylo.TreeConstruction.DistanceCalculator`.
/// Supports `"identity"` plus all 30 built-in substitution matrices
/// (see `distance_models()`), e.g. `"blastn"`, `"hoxd70"`, `"blosum62"`,
/// `"pam250"`, `"dayhoff"`, `"gonnet1992"`. The `"jc69"`/`"jukes-cantor"`
/// and `"k80"`/`"kimura"` corrected-distance models are MoonBit extensions.
pub struct DistanceCalculator {
  model : String
  kind : DistanceModelKind
  scoring_matrix : SubsMatrix?
  skip_letters : Array[String]
}

///|
/// Resolve a model name to a registered scoring matrix.
/// "blastn" maps to the canonical "NUC.4.4" IUPAC matrix; other names are
/// matched case-insensitively against the registry.
fn find_scoring_matrix(model : String) -> SubsMatrix? {
  let target = to_lowercase(model)
  let canonical = if target == "blastn" {
    "nuc.4.4"
  } else {
    target
  }
  let keys = registry_keys.val
  for i in 0..<keys.length() {
    if to_lowercase(keys[i]) == canonical {
      return Some(registry_matrices.val[i])
    }
  }
  None
}

///|
/// Join strings with a separator.
fn join_model_names(names : Array[String], sep : String) -> String {
  let mut result = ""
  for i in 0..<names.length() {
    if i > 0 {
      result = result + sep
    }
    result = result + names[i]
  }
  result
}

///|
/// Resolve a model name to (kind, scoring matrix option).
fn resolve_distance_model(
  model : String,
) -> (DistanceModelKind, SubsMatrix?) raise SeqError {
  let lower_model = to_lowercase(model)
  if model == "identity" {
    (Identity, None)
  } else if lower_model == "jc69" || lower_model == "jukes-cantor" ||
    lower_model == "jukes_cantor" {
    (Jc69, None)
  } else if lower_model == "k80" || lower_model == "k2p" ||
    lower_model == "kimura" {
    (K80, None)
  } else {
    subs_ensure_registry()
    match find_scoring_matrix(model) {
      Some(m) => (Matrix, Some(m))
      None => {
        let available = join_model_names(bio_distance_model_names(), ", ")
        raise SeqError("Model not supported. Available models: " + available)
      }
    }
  }
}

///|
/// Create a distance calculator with a substitution model.
///
/// Raises `SeqError` for unknown models. By default no letters are skipped
/// for the identity model and `"-"`, `"*"` are skipped for matrix models
/// (Biopython behavior); pass `skip_letters` to override.
pub fn DistanceCalculator::new(
  model? : String = "identity",
  skip_letters? : Array[String]? = None,
) -> DistanceCalculator raise SeqError {
  let (kind, scoring_matrix) = resolve_distance_model(model)
  let default_skip : Array[String] = match kind {
    Identity => []
    Matrix => ["-", "*"]
    Jc69 => ["-", "?"]
    K80 => ["-", "?"]
  }
  let resolved_skip = match skip_letters {
    Some(custom) => if custom.length() > 0 { custom } else { default_skip }
    None => default_skip
  }
  DistanceCalculator::{ model, kind, scoring_matrix, skip_letters: resolved_skip }
}

///|
/// Name of the distance model.
pub fn DistanceCalculator::model_name(self : DistanceCalculator) -> String {
  self.model
}

///|
/// Convert a byte/UTF-16 code unit from string indexing to a one-character
/// string. Sequences here are ASCII DNA/protein, so the direct code point
/// conversion is exact.
fn codeunit_to_string(c : UInt16) -> String {
  let buf = StringBuilder()
  buf.write_char(Char::from_int(c.to_int()))
  buf.to_string()
}

///|
/// Check whether a letter is in the skip list.
fn is_skipped(skip_letters : Array[String], letter : String) -> Bool {
  for s in skip_letters {
    if s == letter {
      return true
    }
  }
  false
}

///|
/// Calculate pairwise distance between two sequences.
///
/// For identity: `1 - matches / len(seq1)` (all columns counted, no skipping).
/// For matrix models: skipped columns (`-`, `*`) are ignored; letters absent
/// from the scoring matrix raise `SeqError`; the maximum attainable score is
/// the larger summed self-score of the two sequences; if no column remains
/// the distance is 1. Returns a value in [0, 1] (negative distances possible
/// for asymmetrical matrices).
pub fn DistanceCalculator::pairwise(
  self : DistanceCalculator,
  seq1 : String,
  seq2 : String,
) -> Double raise SeqError {
  self.pairwise_named(seq1, "", seq2, "")
}

///|
/// Pairwise distance with sequence identifiers for error messages.
fn DistanceCalculator::pairwise_named(
  self : DistanceCalculator,
  seq1 : String,
  id1 : String,
  seq2 : String,
  id2 : String,
) -> Double raise SeqError {
  let len1 = seq1.length()
  let len2 = seq2.length()
  let min_len = if len1 < len2 { len1 } else { len2 }

  match self.kind {
    Identity => {
      // Score by character identity. Default skip list is empty, so every
      // column (including gaps) counts; user-supplied skip letters are
      // honored. max_score is always len(seq1).
      if len1 == 0 {
        return 1.0
      }
      let mut score = 0
      for i in 0..<min_len {
        let l1 = codeunit_to_string(seq1.unsafe_get(i))
        let l2 = codeunit_to_string(seq2.unsafe_get(i))
        if !is_skipped(self.skip_letters, l1) &&
          !is_skipped(self.skip_letters, l2) &&
          seq1.unsafe_get(i) == seq2.unsafe_get(i) {
          score = score + 1
        }
      }
      return 1.0 - score.to_double() / len1.to_double()
    }
    Matrix => {
      let m = match self.scoring_matrix {
        Some(m) => m
        None => abort("DistanceCalculator matrix model without scoring matrix")
      }
      let mut score = 0.0
      let mut max_score1 = 0.0
      let mut max_score2 = 0.0
      for i in 0..<min_len {
        let l1 = codeunit_to_string(seq1.unsafe_get(i))
        let l2 = codeunit_to_string(seq2.unsafe_get(i))
        if is_skipped(self.skip_letters, l1) ||
          is_skipped(self.skip_letters, l2) {
          continue
        }
        if !m.has_letter(l1) {
          raise SeqError(
            "Bad letter '" + l1 + "' in sequence '" + id1 +
              "' at position '" + i.to_string() + "'",
          )
        }
        if !m.has_letter(l2) {
          raise SeqError(
            "Bad letter '" + l2 + "' in sequence '" + id2 +
              "' at position '" + i.to_string() + "'",
          )
        }
        max_score1 = max_score1 + m.get_score(l1, l1)
        max_score2 = max_score2 + m.get_score(l2, l2)
        score = score + m.get_score(l1, l2)
      }
      // Take the higher score if the matrix is asymmetrical
      let max_score = if max_score1 > max_score2 {
        max_score1
      } else {
        max_score2
      }
      if max_score == 0.0 {
        return 1.0
      }
      return 1.0 - score / max_score
    }
    Jc69 | K80 => {
      // MoonBit extensions: p-distance based corrections, skipping gaps.
      if min_len == 0 {
        return 1.0
      }
      let mut mismatches = 0
      let mut valid = 0
      for i in 0..<min_len {
        let c1 = seq1.unsafe_get(i)
        let c2 = seq2.unsafe_get(i)
        if !is_skipped(self.skip_letters, codeunit_to_string(c1)) &&
          !is_skipped(self.skip_letters, codeunit_to_string(c2)) {
          valid = valid + 1
          if c1 != c2 {
            mismatches = mismatches + 1
          }
        }
      }
      if valid == 0 {
        return 1.0
      }
      let p = mismatches.to_double() / valid.to_double()
      match self.kind {
        Jc69 => {
          let d = 1.0 - 4.0 / 3.0 * p
          if d <= 0.0 {
            return 0.0
          } else {
            return -0.75 * @math.ln(d)
          }
        }
        K80 => {
          let d = 1.0 - 2.0 * p
          if d <= 0.0 {
            return 0.0
          } else {
            return -0.5 * @math.ln(d)
          }
        }
        _ => abort("unreachable")
      }
    }
  }
}

///|
/// Build a distance matrix from a set of sequences.
pub fn DistanceCalculator::get_distance(
  self : DistanceCalculator,
  sequences : Array[(String, String)],
) -> Array[Array[Double]] raise SeqError {
  let n = sequences.length()
  let matrix : Array[Array[Double]] = Array::new()

  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row.push(0.0)
      } else if j < i {
        // Distances are symmetric; reuse the computed lower triangle.
        row.push(matrix[j][i])
      } else {
        let d = self.pairwise_named(
          sequences[i].1,
          sequences[i].0,
          sequences[j].1,
          sequences[j].0,
        )
        row.push(d)
      }
    }
    matrix.push(row)
  }

  matrix
}

// ---------------------------------------------------------------------------
// Distance model name lists (Biopython DistanceCalculator.models attributes)
// ---------------------------------------------------------------------------

///|
/// Names of DNA distance models (including matrix models whose alphabet is
/// not a superset of the protein alphabet; matches Biopython 1.88).
pub fn distance_dna_models() -> Array[String] {
  subs_ensure_registry()
  bio_dna_model_names()
}

///|
/// Names of protein distance models (matrices whose alphabet contains the
/// full 23-letter Biopython protein alphabet).
pub fn distance_protein_models() -> Array[String] {
  subs_ensure_registry()
  bio_protein_model_names()
}

///|
/// All supported distance model names: "identity" + DNA + protein models.
pub fn distance_models() -> Array[String] {
  subs_ensure_registry()
  bio_distance_model_names()
}

///|
/// Tree constructor using distance-based methods.
pub struct TreeConstructor {
  tree_method : String
}

///|
/// Create a tree constructor with a specific method.
pub fn TreeConstructor::new(tree_method? : String = "nj") -> TreeConstructor {
  TreeConstructor::{ tree_method, }
}

///|
/// Build a phylogenetic tree from a distance matrix.
pub fn TreeConstructor::build_tree(
  self : TreeConstructor,
  distance_matrix : Array[Array[Double]],
  names : Array[String],
) -> Clade {
  if self.tree_method == "upgma" {
    build_upgma(distance_matrix, names)
  } else if self.tree_method == "wpgma" {
    build_wpgma(distance_matrix, names)
  } else {
    build_nj(distance_matrix, names)
  }
}

///|
/// UPGMA (Unweighted Pair Group Method with Arithmetic mean) tree construction.
fn build_upgma(
  distances : Array[Array[Double]],
  names : Array[String],
) -> Clade {
  let n = names.length()

  if n == 0 {
    return Clade::new()
  }
  if n == 1 {
    return Clade::new(name=Some(names[0]))
  }

  let mut clusters : Array[BuildNode] = Array::new()
  for i = 0; i < n; i = i + 1 {
    clusters.push(build_node_leaf(i, names[i]))
  }

  let mut dist = copy_matrix(distances)
  let mut next_id = n

  while clusters.length() > 1 {
    let (i, j) = find_min_pair(dist)
    let dij = dist[i][j]
    let height = dij / 2.0

    let ni = clusters.length()

    let new_cluster = build_node_internal(
      next_id,
      clusters[i],
      clusters[j],
      height - clusters[i].height,
      height - clusters[j].height,
      height,
    )
    next_id = next_id + 1

    let new_dists_row : Array[Double] = Array::new()

    let mut k = 0
    while k < ni {
      if k != i && k != j {
        let size_i = count_leaves(clusters[i])
        let size_j = count_leaves(clusters[j])
        let d = (
            size_i.to_double() * dist[i][k] + size_j.to_double() * dist[j][k]
          ) /
          (size_i.to_double() + size_j.to_double())
        new_dists_row.push(d)
      }
      k = k + 1
    }
    new_dists_row.push(0.0)

    let new_clusters : Array[BuildNode] = Array::new()
    k = 0
    while k < ni {
      if k != i && k != j {
        new_clusters.push(clusters[k])
      }
      k = k + 1
    }
    new_clusters.push(new_cluster)

    let new_dist_matrix : Array[Array[Double]] = Array::new()
    let mut ridx = 0
    let mut r = 0
    while r < ni {
      if r != i && r != j {
        let row : Array[Double] = Array::new()
        let mut c = 0
        while c < ni {
          if c != i && c != j {
            row.push(dist[r][c])
          }
          c = c + 1
        }
        row.push(new_dists_row[ridx])
        new_dist_matrix.push(row)
        ridx = ridx + 1
      }
      r = r + 1
    }
    new_dist_matrix.push(new_dists_row)

    clusters = new_clusters
    dist = new_dist_matrix
  }

  build_node_to_clade(clusters[0])
}

///|
/// WPGMA (Weighted Pair Group Method with Arithmetic mean) tree construction.
fn build_wpgma(
  distances : Array[Array[Double]],
  names : Array[String],
) -> Clade {
  let n = names.length()

  if n == 0 {
    return Clade::new()
  }
  if n == 1 {
    return Clade::new(name=Some(names[0]))
  }

  let mut clusters : Array[BuildNode] = Array::new()
  for i = 0; i < n; i = i + 1 {
    clusters.push(build_node_leaf(i, names[i]))
  }

  let mut dist = copy_matrix(distances)
  let mut next_id = n

  while clusters.length() > 1 {
    let (i, j) = find_min_pair(dist)
    let dij = dist[i][j]
    let height = dij / 2.0

    let ni = clusters.length()

    let new_cluster = build_node_internal(
      next_id,
      clusters[i],
      clusters[j],
      height - clusters[i].height,
      height - clusters[j].height,
      height,
    )
    next_id = next_id + 1

    let new_dists_row : Array[Double] = Array::new()
    let mut k = 0
    while k < ni {
      if k != i && k != j {
        let d = (dist[i][k] + dist[j][k]) / 2.0
        new_dists_row.push(d)
      }
      k = k + 1
    }
    new_dists_row.push(0.0)

    let new_clusters : Array[BuildNode] = Array::new()
    k = 0
    while k < ni {
      if k != i && k != j {
        new_clusters.push(clusters[k])
      }
      k = k + 1
    }
    new_clusters.push(new_cluster)

    let new_dist_matrix : Array[Array[Double]] = Array::new()
    let mut ridx = 0
    let mut r = 0
    while r < ni {
      if r != i && r != j {
        let row : Array[Double] = Array::new()
        let mut c = 0
        while c < ni {
          if c != i && c != j {
            row.push(dist[r][c])
          }
          c = c + 1
        }
        row.push(new_dists_row[ridx])
        new_dist_matrix.push(row)
        ridx = ridx + 1
      }
      r = r + 1
    }
    new_dist_matrix.push(new_dists_row)

    clusters = new_clusters
    dist = new_dist_matrix
  }

  build_node_to_clade(clusters[0])
}

///|
/// Neighbor-Joining tree construction (Saitou & Nei, 1987).
fn build_nj(distances : Array[Array[Double]], names : Array[String]) -> Clade {
  let n = names.length()

  if n == 0 {
    return Clade::new()
  }
  if n == 1 {
    return Clade::new(name=Some(names[0]))
  }
  if n == 2 {
    let d = distances[0][1]
    let c1 = Clade::new(name=Some(names[0]), branch_length=Some(d / 2.0))
    let c2 = Clade::new(name=Some(names[1]), branch_length=Some(d / 2.0))
    return Clade::new(clades=[c1, c2])
  }

  let mut active : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    active.push(i)
  }

  let mut dist = copy_matrix(distances)
  let nodes : Array[BuildNode] = Array::new()
  for i = 0; i < n; i = i + 1 {
    nodes.push(build_node_leaf(i, names[i]))
  }
  let mut next_id = n

  while active.length() > 2 {
    let m = active.length()

    let row_sums : Array[Double] = Array::new()
    for i = 0; i < m; i = i + 1 {
      let mut s = 0.0
      for j = 0; j < m; j = j + 1 {
        s = s + dist[i][j]
      }
      row_sums.push(s)
    }

    let mut min_q = 100000000000000000000000000000.0
    let mut min_i = 0
    let mut min_j = 1

    for i = 0; i < m; i = i + 1 {
      for j = i + 1; j < m; j = j + 1 {
        let factor = if m > 2 { (m - 2).to_double() } else { 1.0 }
        let q = dist[i][j] - (row_sums[i] + row_sums[j]) / factor
        if q < min_q {
          min_q = q
          min_i = i
          min_j = j
        }
      }
    }

    let dij = dist[min_i][min_j]
    let diu = dij / 2.0 +
      (row_sums[min_i] - row_sums[min_j]) / (2.0 * (m - 2).to_double())
    let dju = dij - diu

    let new_node = build_node_internal(
      next_id,
      nodes[active[min_i]],
      nodes[active[min_j]],
      diu,
      dju,
      0.0,
    )
    next_id = next_id + 1

    let idx_new_node = next_id - 1
    nodes.push(new_node)

    let keep_indices : Array[Int] = Array::new()
    for k = 0; k < m; k = k + 1 {
      if k != min_i && k != min_j {
        keep_indices.push(k)
      }
    }

    let new_dist_matrix : Array[Array[Double]] = Array::new()

    for r = 0; r < keep_indices.length(); r = r + 1 {
      let row : Array[Double] = Array::new()
      for c = 0; c < keep_indices.length(); c = c + 1 {
        row.push(dist[keep_indices[r]][keep_indices[c]])
      }
      let ku = keep_indices[r]
      let duk = (dist[min_i][ku] + dist[min_j][ku] - dij) / 2.0
      row.push(duk)
      new_dist_matrix.push(row)
    }

    let last_row : Array[Double] = Array::new()
    for c = 0; c < keep_indices.length(); c = c + 1 {
      let ku = keep_indices[c]
      let duk = (dist[min_i][ku] + dist[min_j][ku] - dij) / 2.0
      last_row.push(duk)
    }
    last_row.push(0.0)
    new_dist_matrix.push(last_row)

    let new_active_list : Array[Int] = Array::new()
    for k in keep_indices {
      new_active_list.push(active[k])
    }
    new_active_list.push(idx_new_node)

    active = new_active_list
    dist = new_dist_matrix
  }

  if active.length() == 2 {
    let d = dist[0][1]
    let n1 = nodes[active[0]]
    let n2 = nodes[active[1]]

    let root = build_node_internal(next_id, n1, n2, d / 2.0, d / 2.0, 0.0)
    build_node_to_clade(root)
  } else {
    build_node_to_clade(nodes[active[0]])
  }
}

///|
/// Copy a distance matrix.
fn copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for row in matrix {
    let new_row : Array[Double] = Array::new()
    for v in row {
      new_row.push(v)
    }
    result.push(new_row)
  }
  result
}

///|
/// Find the pair with minimum distance.
fn find_min_pair(dist : Array[Array[Double]]) -> (Int, Int) {
  let n = dist.length()
  let mut min_d = 100000000000000000000000000000.0
  let mut min_i = 0
  let mut min_j = 1

  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if dist[i][j] < min_d {
        min_d = dist[i][j]
        min_i = i
        min_j = j
      }
    }
  }

  (min_i, min_j)
}

///|
/// Count leaf nodes in a tree.
fn count_leaves(node : BuildNode) -> Int {
  if node.children.length() == 0 {
    1
  } else {
    let mut count = 0
    for child in node.children {
      count = count + count_leaves(child)
    }
    count
  }
}
