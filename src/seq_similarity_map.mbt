/// seq_similarity_map.mbt — sliding-window similarity maps
/// (SimilarityMaps analogue, item B6)
///
/// For every window of the query sequence, the k-mer set overlap with
/// the reference is computed, giving a per-position contribution profile.
/// Peaks mark sub-regions that drive the global similarity. Output is a
/// numeric vector plus an optional ASCII rendering (matching the
/// visualization.mbt ASCII-first style).

///|
/// Sliding-window k-mer Tanimoto profile: entry i is the k-mer set
/// similarity between query[i:i+window] and the whole reference.
pub fn seq_similarity_map(
  query : Seq,
  reference : Seq,
  k~ : Int32 = 3,
  window~ : Int32 = 12,
) -> Array[Double] {
  let q = query.to_string()
  let r = reference.to_string()
  let n = q.length()
  let out : Array[Double] = Array(capacity=0)
  if window <= 0 || n < window {
    return out
  }
  // reference k-mer set (set semantics for the small windows)
  let ref_set : Set[String] = Set([])
  let rn = r.length()
  if rn >= k {
    let mut i = 0
    while i <= rn - k {
      let _ = ref_set.add(r[i:i + k].to_owned())
      i = i + 1
    }
  }
  let mut i = 0
  while i + window <= n {
    // window k-mer set
    let win : Set[String] = Set([])
    if window >= k {
      let mut j = 0
      while j <= window - k {
        let _ = win.add(q[i + j:i + j + k].to_owned())
        j = j + 1
      }
    }
    // Jaccard between window set and reference set
    let mut inter = 0
    for km in win {
      if ref_set.contains(km) {
        inter = inter + 1
      }
    }
    let union = win.length() + ref_set.length() - inter
    if union == 0 {
      out.push(0.0)
    } else {
      out.push(inter.to_double() / union.to_double())
    }
    i = i + 1
  }
  out
}

///|
/// Local maxima of the profile above `threshold` (peak positions).
pub fn seq_similarity_peaks(
  map : Array[Double],
  threshold~ : Double = 0.5,
) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  let n = map.length()
  for i = 0; i < n; i = i + 1 {
    if map[i] < threshold {
      continue
    }
    let left_ok = i == 0 || map[i - 1] <= map[i]
    let right_ok = i == n - 1 || map[i + 1] <= map[i]
    if left_ok && right_ok {
      out.push(i)
    }
  }
  out
}

///|
/// ASCII rendering of a similarity profile: values are bucketed into
/// `width` columns (averaged) and rendered as an 8-level character bar
/// (space .. '#').
pub fn seq_similarity_map_ascii(
  map : Array[Double],
  width~ : Int32 = 60,
) -> String {
  if map.length() == 0 || width <= 0 {
    return ""
  }
  let levels = " .:-=+*#"
  let sb = StringBuilder()
  let nf = map.length().to_double()
  let wf = width.to_double()
  for c = 0; c < width; c = c + 1 {
    let lo = (c.to_double() * nf / wf).floor().to_int()
    let mut hi = ((c.to_double() + 1.0) * nf / wf).floor().to_int()
    if hi <= lo {
      hi = lo + 1
    }
    if hi > map.length() {
      hi = map.length()
    }
    let mut sum = 0.0
    for i = lo; i < hi; i = i + 1 {
      sum = sum + map[i]
    }
    let avg = sum / (hi - lo).to_double()
    let mut idx = (avg * 8.0).floor().to_int()
    if idx > 7 {
      idx = 7
    }
    if idx < 0 {
      idx = 0
    }
    let _ = sb.write_char(levels.unsafe_get(idx).unsafe_to_char())
  }
  sb.to_string()
}
