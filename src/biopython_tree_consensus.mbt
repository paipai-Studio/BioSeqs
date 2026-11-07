/// biopython_tree_consensus.mbt — Consensus tree construction
/// Inspired by Biopython Bio.Phylo.Consensus: majority-rule consensus

///|
/// Frequency of each bipartition across multiple trees.
pub struct BiopythonSplitFrequency {
  split : Array[Bool]
  frequency : Double // 0.0 to 1.0
}

///|
/// Compute split frequencies across multiple trees.
pub fn biopython_split_frequencies(
  trees : Array[BiopythonTreeBipartitions],
) -> Array[BiopythonSplitFrequency] {
  let n_trees = trees.length()
  if n_trees == 0 {
    return []
  }
  let n_taxa = trees[0].taxa.length()
  // Collect all unique splits
  let all_splits : Array[Array[Bool]] = Array(capacity=0)
  let counts : Array[Int32] = Array(capacity=0)
  for t = 0; t < n_trees; t = t + 1 {
    for s = 0; s < trees[t].splits.length(); s = s + 1 {
      let split = trees[t].splits[s]
      let mut found = false
      for u = 0; u < all_splits.length(); u = u + 1 {
        if biopython_splits_equal(all_splits[u], split, n_taxa) {
          counts[u] = counts[u] + 1
          found = true
          break
        }
      }
      if !found {
        all_splits.push(split)
        counts.push(1)
      }
    }
  }
  // Convert to frequencies
  let result : Array[BiopythonSplitFrequency] = Array(capacity=0)
  for i = 0; i < all_splits.length(); i = i + 1 {
    result.push({
      split: all_splits[i],
      frequency: counts[i].to_double() / n_trees.to_double(),
    })
  }
  // Sort by frequency descending
  result.sort_by(fn(a, b) {
    if b.frequency < a.frequency {
      -1
    } else if b.frequency > a.frequency {
      1
    } else {
      0
    }
  })
  result
}

///|
/// Check if two splits are equal (same or complement).
fn biopython_splits_equal(
  s1 : Array[Bool],
  s2 : Array[Bool],
  n : Int32,
) -> Bool {
  let mut same = true
  let mut comp = true
  for j = 0; j < n; j = j + 1 {
    if s1[j] != s2[j] {
      same = false
    }
    if s1[j] == s2[j] {
      comp = false
    }
  }
  same || comp
}

///|
/// Majority-rule consensus tree.
/// Includes splits that appear in > 50% of trees.
pub fn biopython_majority_rule_consensus(
  trees : Array[BiopythonTreeBipartitions],
) -> Array[BiopythonSplitFrequency] {
  let freqs = biopython_split_frequencies(trees)
  let result : Array[BiopythonSplitFrequency] = Array(capacity=0)
  for i = 0; i < freqs.length(); i = i + 1 {
    if freqs[i].frequency > 0.5 {
      result.push(freqs[i])
    }
  }
  result
}

///|
/// Strict consensus tree.
/// Includes only splits that appear in ALL trees.
pub fn biopython_strict_consensus(
  trees : Array[BiopythonTreeBipartitions],
) -> Array[BiopythonSplitFrequency] {
  let freqs = biopython_split_frequencies(trees)
  let result : Array[BiopythonSplitFrequency] = Array(capacity=0)
  for i = 0; i < freqs.length(); i = i + 1 {
    if freqs[i].frequency >= 1.0 {
      result.push(freqs[i])
    }
  }
  result
}

///|
/// Adams consensus tree (simplified).
/// Includes splits compatible with all trees, frequency-weighted.
pub fn biopython_adams_consensus(
  trees : Array[BiopythonTreeBipartitions],
) -> Array[BiopythonSplitFrequency] {
  // Simplified: same as majority rule but with frequency >= 1/3
  let freqs = biopython_split_frequencies(trees)
  let result : Array[BiopythonSplitFrequency] = Array(capacity=0)
  for i = 0; i < freqs.length(); i = i + 1 {
    if freqs[i].frequency >= 0.334 {
      result.push(freqs[i])
    }
  }
  result
}
