// ===========================================================================
// contact_eval.mbt — DCA Contact Map vs 3D Structure Comparison
//
// Evaluates DCA-predicted residue contacts against true contacts extracted
// from a 3D protein structure. Provides precision/recall/F1, stratified
// evaluation by sequence separation, coupling-distance correlation, and
// top-K precision curves.
//
// References:
//   - Morcos et al., PNAS 2011 (DCA contact evaluation)
//   - Ekeberg et al., Phys Rev E 2013 (plmDCA contact evaluation)
//   - Marks et al., PLoS ONE 2011 (contact prediction benchmarks)
// ===========================================================================

///|
/// Result of evaluating predicted contacts against true contacts.
pub(all) struct ContactEvaluation {
  /// Precision = TP / (TP + FP)
  precision : Double
  /// Recall = TP / (TP + FN)
  recall : Double
  /// F1 = 2 * P * R / (P + R)
  f1 : Double
  /// True positives count
  tp : Int
  /// False positives count
  fp : Int
  /// False negatives count
  fn_count : Int
  /// Precision at top-K: [(k, precision), ...]
  top_k_precisions : Array[(Int, Double)]
  /// Mean 3D distance of true positive predictions (Å)
  mean_tp_dist : Double
  /// Mean 3D distance of false positive predictions (Å)
  mean_fp_dist : Double
  /// Number of predicted contacts
  num_predicted : Int
  /// Number of true contacts
  num_true : Int
}

///|
/// Extract true residue contacts from a 3D structure.
/// Two residues i, j are in contact if their CA-CA distance < threshold.
/// Only pairs with |i - j| >= min_seq_dist are included.
///
/// Returns sorted array of (i, j) with i < j.
pub fn extract_true_contacts(
  chain : Chain,
  threshold? : Double = 8.0,
  min_seq_dist? : Int = 5,
) -> Array[(Int, Int)] {
  let residues = chain.get_residues()
  let ca_coords : Array[Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => ca_coords.push(ca.coord)
      None => ()
    }
  }
  let n = ca_coords.length()
  let contacts = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j in (i + 1)..<n {
      if j - i >= min_seq_dist {
        let d = ca_coords[i].distance(ca_coords[j])
        if d < threshold {
          contacts.push((i, j))
        }
      }
    }
  }
  contacts
}

///|
/// Extract true contacts with distance information.
/// Returns array of (i, j, distance) for all qualifying pairs.
pub fn extract_contacts_with_dist(
  chain : Chain,
  threshold? : Double = 8.0,
  min_seq_dist? : Int = 5,
) -> Array[(Int, Int, Double)] {
  let residues = chain.get_residues()
  let ca_coords : Array[Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => ca_coords.push(ca.coord)
      None => ()
    }
  }
  let n = ca_coords.length()
  let contacts = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j in (i + 1)..<n {
      if j - i >= min_seq_dist {
        let d = ca_coords[i].distance(ca_coords[j])
        if d < threshold {
          contacts.push((i, j, d))
        }
      }
    }
  }
  contacts
}

///|
/// Get CA-CA distance for a specific residue pair (i, j).
/// Returns -1.0 if either residue lacks a CA atom.
pub fn residue_ca_distance(
  chain : Chain,
  i : Int,
  j : Int,
) -> Double {
  let residues = chain.get_residues()
  if i < 0 || i >= residues.length() || j < 0 || j >= residues.length() {
    return -1.0
  }
  match (residues[i].get_ca(), residues[j].get_ca()) {
    (Some(ca_i), Some(ca_j)) => ca_i.coord.distance(ca_j.coord)
    _ => -1.0
  }
}

// ---------------------------------------------------------------------------
// Contact set operations
// ---------------------------------------------------------------------------

///|
/// Build a hash set of contact pairs for fast lookup.
/// Key = i * 100000 + j (assuming i, j < 100000)
fn contact_key(i : Int, j : Int) -> Int {
  i * 100000 + j
}

///|
/// Check if a pair (i, j) is in the true contact set.
fn is_true_contact(
  true_map : Map[Int, Unit],
  i : Int,
  j : Int,
) -> Bool {
  let k1 = contact_key(i, j)
  let k2 = contact_key(j, i)
  true_map.contains(k1) || true_map.contains(k2)
}

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

///|
/// Evaluate predicted contacts against true contacts.
/// predicted and true_contacts are arrays of (i, j) pairs with i < j.
pub fn evaluate_contact_prediction(
  predicted : Array[(Int, Int)],
  true_contacts : Array[(Int, Int)],
) -> ContactEvaluation {
  // Build true contact set
  let true_map : Map[Int, Unit] = Map([], capacity=0)
  for pair in true_contacts {
    let (i, j) = pair
    true_map.set(contact_key(i, j), ())
  }
  let mut tp = 0
  let mut fp = 0
  for pair in predicted {
    let (i, j) = pair
    if is_true_contact(true_map, i, j) {
      tp = tp + 1
    } else {
      fp = fp + 1
    }
  }
  let fn_count = true_contacts.length() - tp
  let precision = if predicted.length() > 0 {
    tp.to_double() / predicted.length().to_double()
  } else {
    0.0
  }
  let recall = if true_contacts.length() > 0 {
    tp.to_double() / true_contacts.length().to_double()
  } else {
    0.0
  }
  let f1 = if precision + recall > 0.0 {
    2.0 * precision * recall / (precision + recall)
  } else {
    0.0
  }
  ContactEvaluation::{
    precision,
    recall,
    f1,
    tp,
    fp,
    fn_count: fn_count,
    top_k_precisions: Array(capacity=0),
    mean_tp_dist: 0.0,
    mean_fp_dist: 0.0,
    num_predicted: predicted.length(),
    num_true: true_contacts.length(),
  }
}

///|
/// Evaluate predicted contacts with 3D distance information.
/// Computes mean distances for TP and FP predictions.
pub fn evaluate_with_distances(
  predicted : Array[(Int, Int)],
  true_contacts : Array[(Int, Int)],
  chain : Chain,
) -> ContactEvaluation {
  let base = evaluate_contact_prediction(predicted, true_contacts)
  let true_map : Map[Int, Unit] = Map([], capacity=0)
  for pair in true_contacts {
    let (i, j) = pair
    true_map.set(contact_key(i, j), ())
  }
  let mut tp_dist_sum = 0.0
  let mut tp_count = 0
  let mut fp_dist_sum = 0.0
  let mut fp_count = 0
  for pair in predicted {
    let (i, j) = pair
    let d = residue_ca_distance(chain, i, j)
    if d >= 0.0 {
      if is_true_contact(true_map, i, j) {
        tp_dist_sum = tp_dist_sum + d
        tp_count = tp_count + 1
      } else {
        fp_dist_sum = fp_dist_sum + d
        fp_count = fp_count + 1
      }
    }
  }
  ContactEvaluation::{
    ..base,
    mean_tp_dist: if tp_count > 0 { tp_dist_sum / tp_count.to_double() } else { 0.0 },
    mean_fp_dist: if fp_count > 0 { fp_dist_sum / fp_count.to_double() } else { 0.0 },
  }
}

// ---------------------------------------------------------------------------
// DCA score → contact prediction
// ---------------------------------------------------------------------------

///|
/// Rank residue pairs by score (descending) and return top-K as (i, j) pairs.
/// Only includes pairs with i < j and |i - j| >= min_seq_dist.
pub fn rank_dca_pairs(
  scores : Array[Array[Double]],
  top_k : Int,
  min_seq_dist : Int,
) -> Array[(Int, Double)] {
  // Returns array of (pair_key, score) sorted descending
  let l = scores.length()
  let pairs = Array(capacity=0)
  for i = 0; i < l; i = i + 1 {
    for j in (i + 1)..<l {
      if j - i >= min_seq_dist {
        pairs.push((contact_key(i, j), scores[i][j]))
      }
    }
  }
  // Sort by score descending (selection sort)
  let n = pairs.length()
  for i = 0; i < n; i = i + 1 {
    let mut max_idx = i
    for j in (i + 1)..<n {
      if pairs[j].1 > pairs[max_idx].1 {
        max_idx = j
      }
    }
    if max_idx != i {
      let tmp = pairs[i]
      pairs[i] = pairs[max_idx]
      pairs[max_idx] = tmp
    }
  }
  let k = if top_k < n { top_k } else { n }
  { let r = Array(capacity=k); for i in 0..<k { r.push(pairs[i]) }; r }
}

///|
/// Decode a contact key back to (i, j).
fn decode_key(key : Int) -> (Int, Int) {
  (key / 100000, key % 100000)
}

///|
/// Full DCA vs structure comparison.
/// Takes DCA scores, extracts top-K predicted contacts, compares against
/// true contacts from the 3D structure.
pub fn dca_structure_comparison(
  scores : Array[Array[Double]],
  chain : Chain,
  top_k : Int,
  threshold? : Double = 8.0,
  min_seq_dist? : Int = 5,
) -> ContactEvaluation {
  // Rank DCA scores
  let ranked = rank_dca_pairs(scores, top_k, min_seq_dist)
  let predicted = Array(capacity=ranked.length())
  for pair in ranked {
    let (key, _) = pair
    let (i, j) = decode_key(key)
    predicted.push((i, j))
  }
  // Extract true contacts
  let true_contacts = extract_true_contacts(chain, threshold~, min_seq_dist~)
  // Evaluate with distances
  evaluate_with_distances(predicted, true_contacts, chain)
}

// ---------------------------------------------------------------------------
// Top-K precision curve
// ---------------------------------------------------------------------------

///|
/// Compute precision at multiple top-K values.
/// Returns array of (k, precision) for each k in k_values.
pub fn top_k_precision_curve(
  scores : Array[Array[Double]],
  true_contacts : Array[(Int, Int)],
  k_values : Array[Int],
  min_seq_dist : Int,
) -> Array[(Int, Double)] {
  let l = scores.length()
  // Rank all pairs
  let pairs = Array(capacity=0)
  for i = 0; i < l; i = i + 1 {
    for j in (i + 1)..<l {
      if j - i >= min_seq_dist {
        pairs.push((contact_key(i, j), scores[i][j]))
      }
    }
  }
  // Sort descending
  let n = pairs.length()
  for i = 0; i < n; i = i + 1 {
    let mut max_idx = i
    for j in (i + 1)..<n {
      if pairs[j].1 > pairs[max_idx].1 {
        max_idx = j
      }
    }
    if max_idx != i {
      let tmp = pairs[i]
      pairs[i] = pairs[max_idx]
      pairs[max_idx] = tmp
    }
  }
  // Build true contact set
  let true_map : Map[Int, Unit] = Map([], capacity=0)
  for pair in true_contacts {
    let (i, j) = pair
    true_map.set(contact_key(i, j), ())
  }
  // Compute cumulative TP count
  let result = Array(capacity=k_values.length())
  for kv in k_values {
    let k = if kv < n { kv } else { n }
    let mut tp = 0
    for idx = 0; idx < k; idx = idx + 1 {
      let (key, _) = pairs[idx]
      let (i, j) = decode_key(key)
      if is_true_contact(true_map, i, j) {
        tp = tp + 1
      }
    }
    let prec = if k > 0 { tp.to_double() / k.to_double() } else { 0.0 }
    result.push((kv, prec))
  }
  result
}

// ---------------------------------------------------------------------------
// Stratified evaluation by sequence separation
// ---------------------------------------------------------------------------

///|
/// Evaluate DCA contacts stratified by minimum sequence separation.
/// For each separation s in separations, only considers pairs with |i-j| >= s.
/// Returns array of (separation, evaluation).
pub fn stratified_evaluation(
  scores : Array[Array[Double]],
  chain : Chain,
  separations : Array[Int],
  top_k : Int,
  threshold? : Double = 8.0,
) -> Array[(Int, ContactEvaluation)] {
  let result = Array(capacity=separations.length())
  for sep in separations {
    let eval = dca_structure_comparison(
      scores,
      chain,
      top_k,
      threshold~,
      min_seq_dist=sep,
    )
    result.push((sep, eval))
  }
  result
}

// ---------------------------------------------------------------------------
// Contact map overlap metrics
// ---------------------------------------------------------------------------

///|
/// Jaccard index between two contact maps (boolean matrices).
/// J = |A ∩ B| / |A ∪ B|
pub fn contact_map_jaccard(
  map1 : Array[Array[Bool]],
  map2 : Array[Array[Bool]],
) -> Double {
  let n = map1.length()
  let mut intersection = 0
  let mut union_count = 0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let a = map1[i][j]
      let b = map2[i][j]
      if a && b {
        intersection = intersection + 1
        union_count = union_count + 1
      } else if a || b {
        union_count = union_count + 1
      }
    }
  }
  if union_count == 0 {
    return 0.0
  }
  intersection.to_double() / union_count.to_double()
}

///|
/// Overlap coefficient (Szymkiewicz–Simpson) between two contact maps.
/// O = |A ∩ B| / min(|A|, |B|)
pub fn contact_map_overlap(
  map1 : Array[Array[Bool]],
  map2 : Array[Array[Bool]],
) -> Double {
  let n = map1.length()
  let mut intersection = 0
  let mut count1 = 0
  let mut count2 = 0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let a = map1[i][j]
      let b = map2[i][j]
      if a {
        count1 = count1 + 1
      }
      if b {
        count2 = count2 + 1
      }
      if a && b {
        intersection = intersection + 1
      }
    }
  }
  let min_count = if count1 < count2 { count1 } else { count2 }
  if min_count == 0 {
    return 0.0
  }
  intersection.to_double() / min_count.to_double()
}

// ---------------------------------------------------------------------------
// Coupling score vs 3D distance correlation
// ---------------------------------------------------------------------------

///|
/// Compute Pearson correlation between coupling scores and inverse distances.
/// High-quality DCA scores should anti-correlate with distance (high score →
/// short distance), i.e., correlate with 1/distance.
///
/// Returns (pearson_r, spearman_r) where r is correlation between
/// DCA score and 1/distance.
pub fn coupling_distance_correlation(
  scores : Array[Array[Double]],
  dist_matrix : Array[Array[Double]],
  min_seq_dist : Int,
) -> (Double, Double) {
  let l = scores.length()
  let score_vals = Array(capacity=0)
  let inv_dist_vals = Array(capacity=0)
  for i = 0; i < l; i = i + 1 {
    for j in (i + 1)..<l {
      if j - i >= min_seq_dist {
        let d = dist_matrix[i][j]
        if d > 0.0 {
          score_vals.push(scores[i][j])
          inv_dist_vals.push(1.0 / d)
        }
      }
    }
  }
  let n = score_vals.length()
  if n < 2 {
    return (0.0, 0.0)
  }
  // Pearson correlation
  let mut sx = 0.0
  let mut sy = 0.0
  for i = 0; i < n; i = i + 1 {
    sx = sx + score_vals[i]
    sy = sy + inv_dist_vals[i]
  }
  let mx = sx / n.to_double()
  let my = sy / n.to_double()
  let mut sxy = 0.0
  let mut sxx = 0.0
  let mut syy = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = score_vals[i] - mx
    let dy = inv_dist_vals[i] - my
    sxy = sxy + dx * dy
    sxx = sxx + dx * dx
    syy = syy + dy * dy
  }
  let denom = ce_sqrt(sxx * syy)
  let pearson = if denom > 0.0 { sxy / denom } else { 0.0 }
  // Spearman rank correlation
  let spearman = spearman_correlation(score_vals, inv_dist_vals)
  (pearson, spearman)
}

///|
/// Compute Spearman rank correlation between two arrays.
fn spearman_correlation(
  x : Array[Double],
  y : Array[Double],
) -> Double {
  let n = x.length()
  if n < 2 {
    return 0.0
  }
  let rx = compute_ranks(x)
  let ry = compute_ranks(y)
  let mut sx = 0.0
  let mut sy = 0.0
  for i = 0; i < n; i = i + 1 {
    sx = sx + rx[i]
    sy = sy + ry[i]
  }
  let mx = sx / n.to_double()
  let my = sy / n.to_double()
  let mut sxy = 0.0
  let mut sxx = 0.0
  let mut syy = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = rx[i] - mx
    let dy = ry[i] - my
    sxy = sxy + dx * dy
    sxx = sxx + dx * dx
    syy = syy + dy * dy
  }
  let denom = ce_sqrt(sxx * syy)
  if denom > 0.0 { sxy / denom } else { 0.0 }
}

///|
/// Compute ranks of values (1-based, average for ties).
fn compute_ranks(values : Array[Double]) -> Array[Double] {
  let n = values.length()
  let ranks = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut rank = 0.0
    let mut count = 0
    for j = 0; j < n; j = j + 1 {
      if values[j] < values[i] {
        rank = rank + 1.0
      } else if values[j] == values[i] {
        count = count + 1
      }
    }
    // Average rank for ties: rank + (count + 1) / 2
    ranks[i] = rank + (count.to_double() + 1.0) / 2.0
  }
  ranks
}

///|
/// Square root helper.
fn ce_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = if x > 1.0 { x } else { 1.0 }
  let mut prev = 0.0
  for i = 0; i < 50; i = i + 1 {
    prev = guess
    guess = (guess + x / guess) * 0.5
    if (guess - prev).abs() < 1.0e-15 * guess {
      break
    }
  }
  guess
}

// ---------------------------------------------------------------------------
// Summary report
// ---------------------------------------------------------------------------

///|
/// Comprehensive DCA contact evaluation report.
/// Combines all metrics into a single struct for easy inspection.
pub(all) struct ContactEvalReport {
  /// Overall evaluation (all separations)
  overall : ContactEvaluation
  /// Stratified by sequence separation
  stratified : Array[(Int, ContactEvaluation)]
  /// Top-K precision curve
  precision_curve : Array[(Int, Double)]
  /// Pearson correlation (score vs 1/distance)
  pearson_r : Double
  /// Spearman correlation (score vs 1/distance)
  spearman_r : Double
  /// Jaccard index between predicted and true contact maps
  jaccard : Double
}

///|
/// Generate a comprehensive contact evaluation report.
pub fn contact_eval_report(
  scores : Array[Array[Double]],
  chain : Chain,
  top_k : Int,
  threshold? : Double = 8.0,
  min_seq_dist? : Int = 5,
) -> ContactEvalReport {
  let overall = dca_structure_comparison(
    scores,
    chain,
    top_k,
    threshold~,
    min_seq_dist~,
  )
  let stratified = stratified_evaluation(
    scores,
    chain,
    [5, 12, 24],
    top_k,
    threshold~,
  )
  let true_contacts = extract_true_contacts(chain, threshold~, min_seq_dist~)
  let precision_curve = top_k_precision_curve(
    scores,
    true_contacts,
    [1, 5, 10, 20, 50, 100],
    min_seq_dist,
  )
  // Distance matrix for correlation
  let residues = chain.get_residues()
  let ca_coords : Array[Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => ca_coords.push(ca.coord)
      None => ()
    }
  }
  let n = ca_coords.length()
  let dist_matrix = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      row.push(ca_coords[i].distance(ca_coords[j]))
    }
    dist_matrix.push(row)
  }
  let (pearson_r, spearman_r) = coupling_distance_correlation(
    scores,
    dist_matrix,
    min_seq_dist,
  )
  // Build contact maps for Jaccard
  let l = scores.length()
  let ranked = rank_dca_pairs(scores, top_k, min_seq_dist)
  let pred_map = Array::make(l, Array::make(l, false))
  for pair in ranked {
    let (key, _) = pair
    let (i, j) = decode_key(key)
    if i < l && j < l {
      pred_map[i][j] = true
      pred_map[j][i] = true
    }
  }
  let true_map_arr = Array::make(l, Array::make(l, false))
  for tc in true_contacts {
    let (i, j) = tc
    if i < l && j < l {
      true_map_arr[i][j] = true
      true_map_arr[j][i] = true
    }
  }
  let jaccard = contact_map_jaccard(pred_map, true_map_arr)
  ContactEvalReport::{
    overall,
    stratified,
    precision_curve,
    pearson_r,
    spearman_r,
    jaccard,
  }
}
