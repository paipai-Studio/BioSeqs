// Scikit-bio third-round: tree topology operations (items 16-19)
// prune_taxa, triplet_distance, quartet_distance, transfer_bootstrap_expectation

///|
pub fn skbio3_tree_prune_taxa(
  tip_names : Array[String],
  keep : Array[String],
) -> Array[String] {
  // Prune taxa: return only tip names that are in the keep set
  let keep_set : Array[String] = keep.copy()
  let result : Array[String] = []
  for i in 0..<tip_names.length() {
    let mut found = false
    for j in 0..<keep_set.length() {
      if keep_set[j] == tip_names[i] {
        found = true
        break
      }
    }
    if found {
      result.push(tip_names[i])
    }
  }
  result
}

///|
pub fn skbio3_triplet_distance(
  tips_a : Array[String],
  bipartitions_a : Array[(Array[String], Array[String])],
  tips_b : Array[String],
  bipartitions_b : Array[(Array[String], Array[String])],
) -> Int32 {
  // Triplet distance: count triplets of taxa with different topologies in two trees
  // A triplet is defined by choosing 3 taxa; the topology is determined by
  // which pair is more closely related (sister)
  let n = tips_a.length()
  if n < 3 {
    return 0
  }
  let mut diff_count = 0
  // For each triplet, check if the topology differs
  for i in 0..<n {
    let _next = i + 1
    for j in _next..<n {
      let _next = j + 1
      for k in _next..<n {
        // Get the topology in tree A: which pair is sister?
        let topo_a = skbio3_triplet_topology(
          tips_a[i],
          tips_a[j],
          tips_a[k],
          bipartitions_a,
        )
        let topo_b = skbio3_triplet_topology(
          tips_b[i],
          tips_b[j],
          tips_b[k],
          bipartitions_b,
        )
        if topo_a != topo_b {
          diff_count = diff_count + 1
        }
      }
    }
  }
  diff_count
}

///|
fn skbio3_triplet_topology(
  t1 : String,
  t2 : String,
  t3 : String,
  bipartitions : Array[(Array[String], Array[String])],
) -> Int32 {
  // Determine which pair is sister: returns 0 (t1,t2), 1 (t1,t3), 2 (t2,t3)
  // Check bipartitions: if a bipartition separates {t1,t2} from {t3}, then t1,t2 are sister
  for i in 0..<bipartitions.length() {
    let (side1, _) = bipartitions[i]
    let t1_in_1 = skbio3_in_array(t1, side1)
    let t2_in_1 = skbio3_in_array(t2, side1)
    let t3_in_1 = skbio3_in_array(t3, side1)
    // t1 and t2 on same side, t3 on other → t1,t2 sister
    if t1_in_1 == t2_in_1 && t1_in_1 != t3_in_1 {
      return 0
    }
    // t1 and t3 on same side, t2 on other → t1,t3 sister
    if t1_in_1 == t3_in_1 && t1_in_1 != t2_in_1 {
      return 1
    }
    // t2 and t3 on same side, t1 on other → t2,t3 sister
    if t2_in_1 == t3_in_1 && t2_in_1 != t1_in_1 {
      return 2
    }
  }
  0 // default
}

///|
fn skbio3_in_array(s : String, arr : Array[String]) -> Bool {
  for i in 0..<arr.length() {
    if arr[i] == s {
      return true
    }
  }
  false
}

///|
pub fn skbio3_quartet_distance(
  tips_a : Array[String],
  bipartitions_a : Array[(Array[String], Array[String])],
  tips_b : Array[String],
  bipartitions_b : Array[(Array[String], Array[String])],
) -> Int32 {
  // Quartet distance: count quartets of taxa with different topologies
  let n = tips_a.length()
  if n < 4 {
    return 0
  }
  let mut diff_count = 0
  for i in 0..<n {
    let _next = i + 1
    for j in _next..<n {
      let _next = j + 1
      for k in _next..<n {
        let _next = k + 1
        for l in _next..<n {
          let topo_a = skbio3_quartet_topology(
            tips_a[i],
            tips_a[j],
            tips_a[k],
            tips_a[l],
            bipartitions_a,
          )
          let topo_b = skbio3_quartet_topology(
            tips_b[i],
            tips_b[j],
            tips_b[k],
            tips_b[l],
            bipartitions_b,
          )
          if topo_a != topo_b {
            diff_count = diff_count + 1
          }
        }
      }
    }
  }
  diff_count
}

///|
fn skbio3_quartet_topology(
  t1 : String,
  t2 : String,
  t3 : String,
  t4 : String,
  bipartitions : Array[(Array[String], Array[String])],
) -> Int32 {
  // Quartet topology: which pair of pairs are grouped together
  // Returns: 0 = (t1,t2)|(t3,t4), 1 = (t1,t3)|(t2,t4), 2 = (t1,t4)|(t2,t3)
  for i in 0..<bipartitions.length() {
    let (side1, _) = bipartitions[i]
    let t1_1 = skbio3_in_array(t1, side1)
    let t2_1 = skbio3_in_array(t2, side1)
    let t3_1 = skbio3_in_array(t3, side1)
    let t4_1 = skbio3_in_array(t4, side1)
    // (t1,t2)|(t3,t4)
    if t1_1 == t2_1 && t3_1 == t4_1 && t1_1 != t3_1 {
      return 0
    }
    // (t1,t3)|(t2,t4)
    if t1_1 == t3_1 && t2_1 == t4_1 && t1_1 != t2_1 {
      return 1
    }
    // (t1,t4)|(t2,t3)
    if t1_1 == t4_1 && t2_1 == t3_1 && t1_1 != t2_1 {
      return 2
    }
  }
  0
}

///|
pub fn skbio3_transfer_bootstrap_expectation(
  site_supports : Array[Int32],
  n_taxa : Int32,
) -> Array[Double] {
  // Transfer Bootstrap Expectation (TBE): for each site, compute TBE support
  // TBE = 1 - (min_transfer_distance / (n_taxa - 3))
  // site_supports[i] = number of bootstrap trees that match the reference at site i
  let n_trees = site_supports.length()
  if n_trees == 0 || n_taxa < 4 {
    return []
  }
  let tbe_values : Array[Double] = []
  let denom = (n_taxa - 3).to_double()
  for i in 0..<n_trees {
    // Transfer distance = n_trees - support (simplified)
    let transfer_dist = (n_trees - site_supports[i]).to_double()
    let tbe = 1.0 - transfer_dist / (denom * n_trees.to_double())
    tbe_values.push(if tbe < 0.0 { 0.0 } else if tbe > 1.0 { 1.0 } else { tbe })
  }
  tbe_values
}
