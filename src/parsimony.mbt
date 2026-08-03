///|
/// Parsimony scoring module - Fitch and Sankoff parsimony algorithms.
/// Inspired by Biopython's Bio.Phylo.TreeConstruction.ParsimonyScorer.
///
/// Parsimony methods find the tree that requires the fewest evolutionary
/// changes (character state transitions) to explain the observed data.
///
/// Two algorithms are implemented:
///   1. Fitch parsimony: equal cost for all changes (cost = 1 per transition)
///   2. Sankoff parsimony: variable cost via a cost matrix (dynamic programming)
///
/// The alignment is provided as a Map from sequence name to sequence string.
/// Leaf clades are matched to sequences by their `name` field.

///|
/// A parsimony cost matrix for Sankoff scoring.
/// Maps (from_state, to_state) -> cost.
pub struct ParsimonyMatrix {
  states : Array[String]
  costs : Array[Array[Double]]
}

///|
/// Create a parsimony cost matrix.
pub fn ParsimonyMatrix::new(
  states : Array[String],
  costs : Array[Array[Double]]
) -> ParsimonyMatrix {
  ParsimonyMatrix::{ states, costs }
}

///|
/// Get the cost of transitioning from state i to state j.
pub fn parsimony_matrix_cost(
  matrix : ParsimonyMatrix,
  from_idx : Int,
  to_idx : Int
) -> Double {
  matrix.costs[from_idx][to_idx]
}

///|
/// Get the index of a state in the matrix. Returns -1 if not found.
pub fn parsimony_matrix_index(matrix : ParsimonyMatrix, state : String) -> Int {
  let mut i = 0
  while i < matrix.states.length() {
    if matrix.states[i] == state {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Create a standard Fitch cost matrix for DNA.
/// All changes cost 1, same state costs 0.
pub fn parsimony_fitch_dna_matrix() -> ParsimonyMatrix {
  let states = ["A", "C", "G", "T"]
  let costs = [
    [0.0, 1.0, 1.0, 1.0],
    [1.0, 0.0, 1.0, 1.0],
    [1.0, 1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0, 0.0],
  ]
  ParsimonyMatrix::new(states, costs)
}

///|
/// Create a Sankoff cost matrix for DNA with transition/transversion bias.
/// Transitions (A<->G, C<->T) cost 1, transversions cost 2.
pub fn parsimony_sankoff_dna_matrix(
  transition_cost? : Double = 1.0,
  transversion_cost? : Double = 2.0
) -> ParsimonyMatrix {
  let states = ["A", "C", "G", "T"]
  let costs = [
    [0.0, transversion_cost, transition_cost, transversion_cost],
    [transversion_cost, 0.0, transversion_cost, transition_cost],
    [transition_cost, transversion_cost, 0.0, transversion_cost],
    [transversion_cost, transition_cost, transversion_cost, 0.0],
  ]
  ParsimonyMatrix::new(states, costs)
}

///|
/// Calculate the Fitch parsimony score of a tree given an alignment.
///
/// Uses the Fitch algorithm:
///   - Bottom-up: at each internal node, intersect children's state sets.
///     If intersection is non-empty, use it (no cost). If empty, use union (cost +1).
///   - The total score is the sum of costs across all positions.
pub fn fitch_parsimony_score(
  tree : Tree,
  alignment : Map[String, String]
) -> Double {
  if alignment.size() == 0 {
    return -1.0
  }
  let mut aln_len = 0
  let keys = alignment.keys()
  let key_arr : Array[String] = keys.collect()
  if key_arr.length() > 0 {
    aln_len = alignment[key_arr[0]].length()
  }
  if aln_len == 0 {
    return 0.0
  }
  let mut total_score = 0.0
  let mut pos = 0
  while pos < aln_len {
    total_score = total_score + fitch_score_position(tree.root, alignment, pos)
    pos = pos + 1
  }
  total_score
}

///|
/// Calculate the Fitch parsimony score for a single position.
fn fitch_score_position(
  clade : Clade,
  alignment : Map[String, String],
  pos : Int
) -> Double {
  let (states, score) = fitch_compute_states(clade, alignment, pos)
  ignore(states)
  score
}

///|
/// Recursively compute the set of possible states and accumulated score
/// for a clade at a given alignment position using the Fitch algorithm.
/// Returns (state_set, accumulated_score).
fn fitch_compute_states(
  clade : Clade,
  alignment : Map[String, String],
  pos : Int
) -> (Array[String], Double) {
  if clade.is_terminal() {
    let name = match clade.name {
      Some(n) => n
      None => return ([], 0.0)
    }
    let seq_opt = alignment.get(name)
    match seq_opt {
      Some(seq) =>
        if pos < seq.length() {
          let state = seq[pos:pos + 1].to_owned().to_upper()
          ([state], 0.0)
        } else {
          ([], 0.0)
        }
      None => ([], 0.0)
    }
  } else {
    let mut combined : Array[String] = []
    let mut total_score = 0.0
    let mut first = true
    let mut ci = 0
    while ci < clade.clades.length() {
      let (child_states, child_score) = fitch_compute_states(clade.clades[ci], alignment, pos)
      total_score = total_score + child_score
      if first {
        combined = child_states
        first = false
      } else {
        let intersection = fitch_set_intersection(combined, child_states)
        if intersection.length() > 0 {
          combined = intersection
        } else {
          combined = fitch_set_union(combined, child_states)
          total_score = total_score + 1.0
        }
      }
      ci = ci + 1
    }
    (combined, total_score)
  }
}

///|
/// Compute the intersection of two state sets.
fn fitch_set_intersection(a : Array[String], b : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < a.length() {
    let mut j = 0
    while j < b.length() {
      if a[i] == b[j] {
        result.push(a[i])
        break
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Compute the union of two state sets (preserving order, no duplicates).
fn fitch_set_union(a : Array[String], b : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < a.length() {
    result.push(a[i])
    i = i + 1
  }
  let mut j = 0
  while j < b.length() {
    let mut found = false
    let mut k = 0
    while k < result.length() {
      if result[k] == b[j] {
        found = true
        break
      }
      k = k + 1
    }
    if !found {
      result.push(b[j])
    }
    j = j + 1
  }
  result
}

///|
/// Calculate the Sankoff parsimony score of a tree given an alignment and cost matrix.
///
/// Uses the Sankoff dynamic programming algorithm:
///   - For each node v and state t, compute s_t(v) = minimum cost of the
///     subtree rooted at v if v has state t.
///   - For leaves: s_t(leaf) = 0 if leaf has state t, infinity otherwise.
///   - For internal nodes: s_t(v) = sum over children c of
///       min over i of (cost(t, i) + s_i(c))
///   - Total score = min over t of s_t(root).
pub fn sankoff_parsimony_score(
  tree : Tree,
  alignment : Map[String, String],
  matrix : ParsimonyMatrix
) -> Double {
  if alignment.size() == 0 {
    return -1.0
  }
  let mut aln_len = 0
  let keys = alignment.keys()
  let key_arr : Array[String] = keys.collect()
  if key_arr.length() > 0 {
    aln_len = alignment[key_arr[0]].length()
  }
  if aln_len == 0 {
    return 0.0
  }
  let mut total_score = 0.0
  let mut pos = 0
  while pos < aln_len {
    total_score = total_score + sankoff_score_position(tree.root, alignment, pos, matrix)
    pos = pos + 1
  }
  total_score
}

///|
/// Calculate the Sankoff parsimony score for a single position.
fn sankoff_score_position(
  clade : Clade,
  alignment : Map[String, String],
  pos : Int,
  matrix : ParsimonyMatrix
) -> Double {
  let scores = sankoff_compute_scores(clade, alignment, pos, matrix)
  let mut min_score = 1.0e30
  let mut i = 0
  while i < scores.length() {
    if scores[i] < min_score {
      min_score = scores[i]
    }
    i = i + 1
  }
  if min_score >= 1.0e30 { 0.0 } else { min_score }
}

///|
/// Recursively compute the Sankoff score vector for a clade.
/// Returns an array where scores[i] = minimum cost if this clade has state i.
fn sankoff_compute_scores(
  clade : Clade,
  alignment : Map[String, String],
  pos : Int,
  matrix : ParsimonyMatrix
) -> Array[Double] {
  let n_states = matrix.states.length()
  if clade.is_terminal() {
    let scores = Array::make(n_states, 1.0e30)
    let name = match clade.name {
      Some(n) => n
      None => return scores
    }
    let seq_opt = alignment.get(name)
    match seq_opt {
      Some(seq) =>
        if pos < seq.length() {
          let state = seq[pos:pos + 1].to_owned().to_upper()
          let idx = parsimony_matrix_index(matrix, state)
          if idx >= 0 {
            scores[idx] = 0.0
          }
        }
      None => ()
    }
    scores
  } else {
    let scores = Array::make(n_states, 0.0)
    let mut ci = 0
    while ci < clade.clades.length() {
      let child = clade.clades[ci]
      let child_scores = sankoff_compute_scores(child, alignment, pos, matrix)
      let mut t = 0
      while t < n_states {
        let mut min_child = 1.0e30
        let mut i = 0
        while i < n_states {
          let cost = matrix.costs[t][i] + child_scores[i]
          if cost < min_child {
            min_child = cost
          }
          i = i + 1
        }
        scores[t] = scores[t] + min_child
        t = t + 1
      }
      ci = ci + 1
    }
    scores
  }
}

///|
/// Calculate parsimony score using either Fitch (no matrix) or Sankoff (with matrix).
pub fn parsimony_get_score(
  tree : Tree,
  alignment : Map[String, String],
  matrix? : ParsimonyMatrix? = None
) -> Double {
  match matrix {
    Some(m) => sankoff_parsimony_score(tree, alignment, m)
    None => fitch_parsimony_score(tree, alignment)
  }
}

///|
/// Create a sample alignment for testing parsimony.
/// Returns a 4-taxon, 5-character alignment.
pub fn parsimony_sample_alignment() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=4)
  m["A"] = "ACGTA"
  m["B"] = "ACGTA"
  m["C"] = "ATCTA"
  m["D"] = "ATCTA"
  m
}

///|
/// Create a sample tree for testing parsimony.
/// Topology: ((A, B), (C, D))
pub fn parsimony_sample_tree() -> Tree {
  let leaf_a = Clade::new(name=Some("A"))
  let leaf_b = Clade::new(name=Some("B"))
  let leaf_c = Clade::new(name=Some("C"))
  let leaf_d = Clade::new(name=Some("D"))
  let inner_left = Clade::new(clades=[leaf_a, leaf_b])
  let inner_right = Clade::new(clades=[leaf_c, leaf_d])
  let root = Clade::new(clades=[inner_left, inner_right])
  Tree::new(root, rooted=true)
}

///|
/// Create a sample alignment with more variation.
pub fn parsimony_sample_alignment_varied() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=5)
  m["A"] = "ACGTACGT"
  m["B"] = "ACGTTCGT"
  m["C"] = "ATCTAGCT"
  m["D"] = "ATCTAGCT"
  m["E"] = "ATCTAGCC"
  m
}
