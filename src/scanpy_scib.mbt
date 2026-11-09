///|
/// scIB-style batch-integration suite.  `scanpy_scib_correct` dispatches to
/// the batch-correction backends that already exist in the repo (Harmony,
/// Scanorama; `none` returns an untouched copy), and `scanpy_scib_metrics`
/// scores an embedding with the standard scIB evaluation metrics: graph
/// connectivity, iLISI (batch mixing), kBET (neighborhood batch composition)
/// and silhouette scores by batch (should be low) and by cell label (should
/// be high).

///|
/// Result of `scanpy_scib_correct`.
pub struct ScanpyScibCorrection {
  corrected : Array[Array[Double]]
  method_ : String
  n_batches : Int32
}

///|
/// Batch-correct an embedding.  `method_` is one of `"harmony"` (default),
/// `"scanorama"` (mutual-NN alignment per batch) or `"none"` (identity copy).
pub fn scanpy_scib_correct(
  embedding : Array[Array[Double]],
  batch : Array[Int32],
  method_~ : String = "harmony",
  n_clusters~ : Int32 = 10,
  max_iter~ : Int32 = 30,
  lambda~ : Double = 1.0,
) -> ScanpyScibCorrection {
  let n = embedding.length()
  let batches = Set([])
  for b in batch {
    batches.add(b)
  }
  let n_batches = batches.length()
  match method_ {
    "none" => {
      let copy : Array[Array[Double]] = Array(capacity=n)
      for row in embedding {
        copy.push(row.copy())
      }
      ScanpyScibCorrection::{ corrected: copy, method_: "none", n_batches, }
    }
    "scanorama" => {
      // split per batch, align, then reassemble in the original order
      let per_batch : Map[Int32, Array[Int32]] = Map([], capacity=n_batches)
      for i = 0; i < n; i = i + 1 {
        match per_batch.get(batch[i]) {
          None => {
            let idx : Array[Int32] = Array(capacity=0)
            idx.push(i)
            per_batch.set(batch[i], idx)
          }
          Some(idx) => idx.push(i)
        }
      }
      let grouped : Array[Array[Array[Double]]] = Array(capacity=n_batches)
      let group_order : Array[Int32] = Array(capacity=n_batches)
      for b in batches {
        let idx = per_batch.get(b).unwrap()
        let rows : Array[Array[Double]] = Array(capacity=idx.length())
        for i in idx {
          rows.push(embedding[i])
        }
        grouped.push(rows)
        group_order.push(b)
      }
      let aligned = scanpy_scanorama(grouped)
      // scanpy_scanorama returns all cells flattened in `grouped` order
      let corrected : Array[Array[Double]] = Array(capacity=n)
      for _i = 0; _i < n; _i = _i + 1 {
        corrected.push([])
      }
      let mut offset = 0
      for b in group_order {
        let idx = per_batch.get(b).unwrap()
        for t = 0; t < idx.length(); t = t + 1 {
          corrected[idx[t]] = aligned[offset + t]
        }
        offset = offset + idx.length()
      }
      ScanpyScibCorrection::{ corrected, method_: "scanorama", n_batches, }
    }
    _ => {
      let res = scanpy_harmonize(
        embedding,
        batch,
        n_clusters~,
        max_iter~,
        lambda~,
      )
      ScanpyScibCorrection::{
        corrected: res.corrected,
        method_: "harmony",
        n_batches,
      }
    }
  }
}

///|
/// Mean silhouette score over a grouping (ids per row).  For each point the
/// score compares the mean distance to its own group (excluding itself) with
/// the minimum mean distance to any other group; single-group / single-member
/// edge cases score 0.
pub fn scanpy_scib_silhouette(
  embedding : Array[Array[Double]],
  group : Array[Int32],
) -> Double {
  let n = embedding.length()
  if n < 2 {
    return 0.0
  }
  let groups = Set([])
  for g in group {
    groups.add(g)
  }
  let n_groups = groups.length()
  if n_groups < 2 {
    return 0.0
  }
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    let own = group[i]
    // mean distance to each group
    let sums : Map[Int32, Double] = Map([], capacity=n_groups)
    let counts : Map[Int32, Int32] = Map([], capacity=n_groups)
    for j = 0; j < n; j = j + 1 {
      if j == i {
        continue
      }
      let mut d2 = 0.0
      let dim = if embedding[i].length() < embedding[j].length() {
        embedding[i].length()
      } else {
        embedding[j].length()
      }
      for d = 0; d < dim; d = d + 1 {
        let diff = embedding[i][d] - embedding[j][d]
        d2 = d2 + diff * diff
      }
      let g = group[j]
      match sums.get(g) {
        None => sums.set(g, d2.sqrt())
        Some(x) => sums.set(g, x + d2.sqrt())
      }
      match counts.get(g) {
        None => counts.set(g, 1)
        Some(c) => counts.set(g, c + 1)
      }
    }
    let a = match sums.get(own) {
      None => 0.0
      Some(s) => s / counts.get(own).unwrap().to_double()
    }
    let mut b = 1.0e300
    for g in sums {
      if g.0 != own {
        let m = g.1 / counts.get(g.0).unwrap().to_double()
        if m < b {
          b = m
        }
      }
    }
    let s = if a > 0.0 || b < 1.0e299 {
      let denom = if a > b { a } else { b }
      if denom > 0.0 {
        (b - a) / denom
      } else {
        0.0
      }
    } else {
      0.0
    }
    total = total + s
  }
  total / n.to_double()
}

///|
/// kNN graph connectivity: the fraction `sum |c|^2 / n^2` of connected-
/// component sizes (1.0 when the graph is a single component).  Uses the
/// mutual-union of directed kNN edges.
pub fn scanpy_scib_graph_connectivity(
  embedding : Array[Array[Double]],
  k~ : Int32 = 15,
) -> Double {
  let n = embedding.length()
  if n == 0 {
    return 0.0
  }
  let kk = if k < n { k } else { n }
  let (knn_idx, _d) = scanpy_knn_distances(embedding, kk)
  // union-find over kNN edges (skip self edges)
  let parent = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    parent[i] = i
  }
  let find = fn(v0 : Int32) -> Int32 {
    let mut v = v0
    while parent[v] != v {
      parent[v] = parent[parent[v]]
      v = parent[v]
    }
    v
  }
  for i = 0; i < n; i = i + 1 {
    for t = 1; t < kk; t = t + 1 {
      let j = knn_idx[i][t]
      if j != i {
        let ri = find(i)
        let rj = find(j)
        if ri != rj {
          parent[ri] = rj
        }
      }
    }
  }
  let sizes : Map[Int32, Int32] = Map([], capacity=16)
  for i = 0; i < n; i = i + 1 {
    let r = find(i)
    match sizes.get(r) {
      None => sizes.set(r, 1)
      Some(c) => sizes.set(r, c + 1)
    }
  }
  let mut acc = 0.0
  for kv in sizes {
    let s = kv.1.to_double()
    acc = acc + s * s
  }
  acc / (n.to_double() * n.to_double())
}

///|
/// scIB metric bundle over an integrated embedding:
///
/// - `graph_connectivity`: see `scanpy_scib_graph_connectivity` (higher is
///   better),
/// - `ilisi`: mean exponential Shannon entropy of the batch distribution in
///   each cell's kNN neighborhood (higher is better),
/// - `kbet`: mean overlap `1 - L1(p_neighborhood, p_global) / 2` of the
///   neighborhood batch distribution with the global one (higher is better),
/// - `silhouette_batch` / `silhouette_label`: `scanpy_scib_silhouette` by
///   batch (lower is better) and by label (higher is better).
pub struct ScanpyScibMetrics {
  graph_connectivity : Double
  ilisi : Double
  kbet : Double
  silhouette_batch : Double
  silhouette_label : Double
}

///|
/// Compute the scIB metrics for an embedding with per-cell batch ids and
/// (optional) cell-label ids.  Pass an empty `labels` to skip the label
/// silhouette (reported as 0.0).
pub fn scanpy_scib_metrics(
  embedding : Array[Array[Double]],
  batch : Array[Int32],
  labels : Array[Int32],
  k~ : Int32 = 15,
) -> ScanpyScibMetrics {
  let n = embedding.length()
  if n == 0 {
    return ScanpyScibMetrics::{
      graph_connectivity: 0.0,
      ilisi: 0.0,
      kbet: 0.0,
      silhouette_batch: 0.0,
      silhouette_label: 0.0,
    }
  }
  let conn = scanpy_scib_graph_connectivity(embedding, k~)
  let kk = if k < n { k } else { n }
  let (knn_idx, _d) = scanpy_knn_distances(embedding, kk)
  // global batch distribution
  let batch_ids = Set([])
  for b in batch {
    batch_ids.add(b)
  }
  let global_counts : Map[Int32, Int32] = Map([], capacity=batch_ids.length())
  for b in batch {
    match global_counts.get(b) {
      None => global_counts.set(b, 1)
      Some(c) => global_counts.set(b, c + 1)
    }
  }
  let global_p : Map[Int32, Double] = Map([], capacity=batch_ids.length())
  for kv in global_counts {
    global_p.set(kv.0, kv.1.to_double() / n.to_double())
  }
  // iLISI + kBET over kNN neighborhoods (including the cell itself)
  let mut ilisi_sum = 0.0
  let mut kbet_sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let local_counts : Map[Int32, Int32] = Map([], capacity=batch_ids.length())
    for t = 0; t < kk; t = t + 1 {
      let j = knn_idx[i][t]
      let b = batch[j]
      match local_counts.get(b) {
        None => local_counts.set(b, 1)
        Some(c) => local_counts.set(b, c + 1)
      }
    }
    // iLISI = exp(H(p))
    let mut h = 0.0
    for kv in local_counts {
      let p = kv.1.to_double() / kk.to_double()
      if p > 0.0 {
        h = h - p * @math.ln(p)
      }
    }
    ilisi_sum = ilisi_sum + @math.exp(h)
    // kBET = 1 - L1(p_local, p_global) / 2
    let mut l1 = 0.0
    for kv in global_p {
      let pl = match local_counts.get(kv.0) {
        None => 0.0
        Some(c) => c.to_double() / kk.to_double()
      }
      l1 = l1 + (pl - kv.1).abs()
    }
    kbet_sum = kbet_sum + (1.0 - l1 / 2.0).max(0.0)
  }
  let label_sil = if labels.length() == n {
    scanpy_scib_silhouette(embedding, labels)
  } else {
    0.0
  }
  ScanpyScibMetrics::{
    graph_connectivity: conn,
    ilisi: ilisi_sum / n.to_double(),
    kbet: kbet_sum / n.to_double(),
    silhouette_batch: scanpy_scib_silhouette(embedding, batch),
    silhouette_label: label_sil,
  }
}
