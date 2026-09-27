///|
/// Maximum-likelihood phylogenetic tree inference (Jukes-Cantor model).
///
/// Complements the parsimony constructor with a likelihood-based NNI search:
/// given a multiple-sequence alignment, find the tree topology (reachable from
/// a starting tree by nearest-neighbor-interchange) that maximizes the
/// Felsenstein pruning log-likelihood under the Jukes-Cantor (JC69) model.

///|
/// Map a nucleotide code unit to a state index (0=A,1=C,2=G,3=T); -1 for
/// gaps/unknown characters. Lower-case is folded to upper-case.
fn ml_base_to_state(c : UInt16) -> Int {
  let upper = if c >= 97 && c <= 122 { c - 32 } else { c }
  if upper == 65 {
    0
  } else if upper == 67 {
    1
  } else if upper == 71 {
    2
  } else if upper == 84 {
    3
  } else {
    -1
  }
}

///|
/// Jukes-Cantor transition probability: same state.
fn ml_jc_same(t : Double) -> Double {
  0.25 + 0.75 * @math.exp(-4.0 * t / 3.0)
}

///|
/// Jukes-Cantor transition probability: different state.
fn ml_jc_diff(t : Double) -> Double {
  0.25 - 0.25 * @math.exp(-4.0 * t / 3.0)
}

///|
/// Build a per-site map (taxon name -> nucleotide code unit) from the aligned
/// sequences for one column.
fn ml_build_site(alignment : Map[String, String], col : Int) -> Map[String, UInt16] {
  let site : Map[String, UInt16] = Map([], capacity=0)
  for (name, seq) in alignment {
    if col < seq.length() {
      site.set(name, seq.unsafe_get(col))
    } else {
      site.set(name, 45)
    }
  }
  site
}

///|
/// Felsenstein pruning: conditional likelihood vector (4 states) of the subtree
/// rooted at `clade`, for one site. Leaf = indicator vector (or uniform for
/// gaps); internal node = product over children of P·L.
fn ml_node_likelihood(
  clade : Clade,
  site : Map[String, UInt16],
  default_len : Double,
) -> Array[Double] {
  if clade.clades.length() == 0 {
    let name = match clade.name {
      Some(n) => n
      None => ""
    }
    let state = match site.get(name) {
      Some(c) => ml_base_to_state(c)
      None => -1
    }
    if state >= 0 {
      let lik = [0.0, 0.0, 0.0, 0.0]
      lik[state] = 1.0
      lik
    } else {
      [1.0, 1.0, 1.0, 1.0]
    }
  } else {
    let mut result = [1.0, 1.0, 1.0, 1.0]
    for child in clade.clades {
      let bl = match child.branch_length {
        Some(l) => l
        None => default_len
      }
      let child_lik = ml_node_likelihood(child, site, default_len)
      let updated = [0.0, 0.0, 0.0, 0.0]
      for s = 0; s < 4; s = s + 1 {
        let mut sum = 0.0
        for sp = 0; sp < 4; sp = sp + 1 {
          let p = if s == sp { ml_jc_same(bl) } else { ml_jc_diff(bl) }
          sum = sum + p * child_lik[sp]
        }
        updated[s] = result[s] * sum
      }
      result = updated
    }
    result
  }
}

///|
/// Log-likelihood of a tree under the Jukes-Cantor model, via Felsenstein
/// pruning. `alignment` maps taxon name -> aligned sequence (all equal length).
/// `branch_length` is used for any branch that lacks its own length.
pub fn ml_tree_log_likelihood(
  tree : Tree,
  alignment : Map[String, String],
  branch_length? : Double = 0.1,
) -> Double {
  let seqs = alignment.values().collect()
  if seqs.length() == 0 {
    return 0.0
  }
  let n_cols = seqs[0].length()
  let mut total = 0.0
  for col = 0; col < n_cols; col = col + 1 {
    let site = ml_build_site(alignment, col)
    let root_lik = ml_node_likelihood(tree.root, site, branch_length)
    // stationary distribution pi = 1/4 each
    let site_lik = 0.25 * (root_lik[0] + root_lik[1] + root_lik[2] + root_lik[3])
    if site_lik > 0.0 {
      total = total + @math.ln(site_lik)
    }
  }
  total
}

///|
/// NNI hill-climbing that maximizes the log-likelihood, starting from
/// `starting_tree`. Stops when no neighbor improves the score (or `max_iter`
/// sweeps have run). Missing branch lengths fall back to `branch_length`, so
/// topology-only trees can be searched directly.
pub fn ml_nni_search(
  starting_tree : Tree,
  alignment : Map[String, String],
  branch_length? : Double = 0.1,
  max_iter? : Int = 100,
) -> Tree {
  let mut best_tree = starting_tree
  let mut best_score = ml_tree_log_likelihood(best_tree, alignment, branch_length~)
  let mut iter = 0
  while iter < max_iter {
    let mut improved = false
    let mut local_best = best_tree
    let mut local_score = best_score
    for neighbor in nni_get_neighbors(best_tree) {
      let score = ml_tree_log_likelihood(neighbor, alignment, branch_length~)
      if score > local_score {
        local_score = score
        local_best = neighbor
      }
    }
    if local_score > best_score {
      best_score = local_score
      best_tree = local_best
      improved = true
    }
    if not(improved) {
      break
    }
    iter = iter + 1
  }
  best_tree
}
