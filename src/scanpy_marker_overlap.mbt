/// Scanpy-inspired: tl.marker_gene_overlap — compare marker genes with reference sets

///|
/// Compute overlap between computed marker genes and reference marker sets.
/// Uses Jaccard index: |A ∩ B| / |A ∪ B|
///
/// Parameters:
/// - `computed_markers`: marker genes per group (from rank_genes_groups)
/// - `reference_markers`: reference marker gene sets (map: name → genes)
/// - `top_n`: number of top markers to consider (default 100)
pub fn scanpy_marker_gene_overlap(
  computed_markers : Array[Array[String]],
  reference_markers : Map[String, Array[String]],
  top_n~ : Int32 = 100,
) -> Map[String, Map[String, Double]] {
  let result : Map[String, Map[String, Double]] = Map([], capacity=0)
  let n_groups = computed_markers.length()
  let ref_sets : Map[String, Set[String]] = Map([], capacity=0)
  for ref_key in reference_markers.keys() {
    match reference_markers.get(ref_key) {
      Some(genes) => {
        let s : Set[String] = Set::from_iter(genes.iter())
        ref_sets[ref_key] = s
      }
      None => ()
    }
  }
  for g in 0..<n_groups {
    let n_take = if top_n < computed_markers[g].length() {
      top_n
    } else {
      computed_markers[g].length()
    }
    let computed_set : Set[String] = Set([], capacity=n_take)
    for j in 0..<n_take {
      computed_set.add(computed_markers[g][j])
    }
    let group_name = "group_" + g.to_string()
    let overlaps : Map[String, Double] = Map([], capacity=0)
    for ref_name in ref_sets.keys() {
      match ref_sets.get(ref_name) {
        Some(ref_set) => {
          let mut intersection = 0
          for gene in computed_set {
            if ref_set.contains(gene) {
              intersection += 1
            }
          }
          let union_size = computed_set.length() +
            ref_set.length() -
            intersection
          let jaccard = if union_size > 0 {
            intersection.to_double() / union_size.to_double()
          } else {
            0.0
          }
          overlaps[ref_name] = jaccard
        }
        None => ()
      }
    }
    result[group_name] = overlaps
  }
  result
}

///|
/// Compute overlap using overlap coefficient (|A ∩ B| / min(|A|, |B|)).
pub fn scanpy_marker_overlap_coefficient(
  computed_markers : Array[Array[String]],
  reference_markers : Map[String, Array[String]],
  top_n~ : Int32 = 100,
) -> Map[String, Map[String, Double]] {
  let result : Map[String, Map[String, Double]] = Map([], capacity=0)
  let n_groups = computed_markers.length()
  let ref_sets : Map[String, Set[String]] = Map([], capacity=0)
  for ref_key in reference_markers.keys() {
    match reference_markers.get(ref_key) {
      Some(genes) => {
        let s : Set[String] = Set::from_iter(genes.iter())
        ref_sets[ref_key] = s
      }
      None => ()
    }
  }
  for g in 0..<n_groups {
    let n_take = if top_n < computed_markers[g].length() {
      top_n
    } else {
      computed_markers[g].length()
    }
    let computed_set : Set[String] = Set([], capacity=n_take)
    for j in 0..<n_take {
      computed_set.add(computed_markers[g][j])
    }
    let group_name = "group_" + g.to_string()
    let overlaps : Map[String, Double] = Map([], capacity=0)
    for ref_name in ref_sets.keys() {
      match ref_sets.get(ref_name) {
        Some(ref_set) => {
          let mut intersection = 0
          for gene in computed_set {
            if ref_set.contains(gene) {
              intersection += 1
            }
          }
          let min_size = if computed_set.length() < ref_set.length() {
            computed_set.length()
          } else {
            ref_set.length()
          }
          let coeff = if min_size > 0 {
            intersection.to_double() / min_size.to_double()
          } else {
            0.0
          }
          overlaps[ref_name] = coeff
        }
        None => ()
      }
    }
    result[group_name] = overlaps
  }
  result
}
