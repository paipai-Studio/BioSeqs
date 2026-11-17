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
/// Adjusted Rand index between two labelings (truth vs. predicted clusters):
/// `RI = sum_ij C(n_ij,2)`, chance-corrected against the hypergeometric
/// expectation; 1.0 = perfect agreement, 0.0 = random overlap, negative =
/// worse than chance.  Degenerate (single-class) identical partitions score
/// 1.0, matching sklearn.
pub fn scanpy_scib_ari(truth : Array[Int32], pred : Array[Int32]) -> Double {
  let n = truth.length()
  if n < 2 || pred.length() != n {
    return 0.0
  }
  let pair : Map[(Int32, Int32), Int32] = Map([], capacity=n)
  let row : Map[Int32, Int32] = Map([], capacity=n)
  let col : Map[Int32, Int32] = Map([], capacity=n)
  for i in 0..<n {
    let t = truth[i]
    let p = pred[i]
    match pair.get((t, p)) {
      None => pair.set((t, p), 1)
      Some(c) => pair.set((t, p), c + 1)
    }
    match row.get(t) {
      None => row.set(t, 1)
      Some(c) => row.set(t, c + 1)
    }
    match col.get(p) {
      None => col.set(p, 1)
      Some(c) => col.set(p, c + 1)
    }
  }
  let comb2 = fn(m : Int32) -> Double {
    let v = m.to_double()
    v * (v - 1.0) / 2.0
  }
  let mut sum_ij = 0.0
  for kv in pair {
    sum_ij = sum_ij + comb2(kv.1)
  }
  let mut sum_t = 0.0
  for kv in row {
    sum_t = sum_t + comb2(kv.1)
  }
  let mut sum_p = 0.0
  for kv in col {
    sum_p = sum_p + comb2(kv.1)
  }
  let expected = sum_t * sum_p / comb2(n)
  let denom = 0.5 * (sum_t + sum_p) - expected
  if denom <= 0.0 {
    return 1.0
  }
  (sum_ij - expected) / denom
}

///|
/// Normalized mutual information between two labelings, normalized by the
/// arithmetic mean of the marginal entropies (sklearn default); 1.0 when the
/// two labelings carry identical information.  Two constant labelings agree
/// trivially and score 1.0.
pub fn scanpy_scib_nmi(truth : Array[Int32], pred : Array[Int32]) -> Double {
  let n = truth.length()
  if n < 2 || pred.length() != n {
    return 0.0
  }
  let pair : Map[(Int32, Int32), Int32] = Map([], capacity=n)
  let row : Map[Int32, Int32] = Map([], capacity=n)
  let col : Map[Int32, Int32] = Map([], capacity=n)
  for i in 0..<n {
    let t = truth[i]
    let p = pred[i]
    match pair.get((t, p)) {
      None => pair.set((t, p), 1)
      Some(c) => pair.set((t, p), c + 1)
    }
    match row.get(t) {
      None => row.set(t, 1)
      Some(c) => row.set(t, c + 1)
    }
    match col.get(p) {
      None => col.set(p, 1)
      Some(c) => col.set(p, c + 1)
    }
  }
  let nf = n.to_double()
  let entropy = fn(counts : Map[Int32, Int32]) -> Double {
    let mut h = 0.0
    for kv in counts {
      let p = kv.1.to_double() / nf
      if p > 0.0 {
        h = h - p * @math.ln(p)
      }
    }
    h
  }
  let h_t = entropy(row)
  let h_p = entropy(col)
  let mut mi = 0.0
  for kv in pair {
    let p_ij = kv.1.to_double() / nf
    let p_t = row.get(kv.0.0).unwrap().to_double() / nf
    let p_p = col.get(kv.0.1).unwrap().to_double() / nf
    mi = mi + p_ij * @math.ln(p_ij / (p_t * p_p))
  }
  let norm = 0.5 * (h_t + h_p)
  if norm <= 0.0 {
    return if (h_t - h_p).abs() < 1.0e-12 { 1.0 } else { 0.0 }
  }
  (mi / norm).max(0.0)
}

///|
/// Shared kNN-neighborhood batch/label mixing kernel behind
/// `scanpy_scib_ilisi` / `scanpy_scib_kbet`.  Returns
/// `(mean_iLISI, mean_kBET)` over each cell's k nearest neighbours
/// (including the cell itself).
fn scib_knn_mixing(
  embedding : Array[Array[Double]],
  group : Array[Int32],
  k : Int32,
) -> (Double, Double) {
  let n = embedding.length()
  if n == 0 {
    return (0.0, 0.0)
  }
  let kk = if k < n { k } else { n }
  if kk == 0 {
    return (0.0, 0.0)
  }
  let (knn_idx, _d) = scanpy_knn_distances(embedding, kk)
  let group_ids = Set([])
  for g in group {
    group_ids.add(g)
  }
  let global_counts : Map[Int32, Int32] = Map([], capacity=group_ids.length())
  for g in group {
    match global_counts.get(g) {
      None => global_counts.set(g, 1)
      Some(c) => global_counts.set(g, c + 1)
    }
  }
  let global_p : Map[Int32, Double] = Map([], capacity=group_ids.length())
  for kv in global_counts {
    global_p.set(kv.0, kv.1.to_double() / n.to_double())
  }
  let mut ilisi_sum = 0.0
  let mut kbet_sum = 0.0
  for i in 0..<n {
    let local_counts : Map[Int32, Int32] = Map([], capacity=group_ids.length())
    for t in 0..<kk {
      let j = knn_idx[i][t]
      let g = group[j]
      match local_counts.get(g) {
        None => local_counts.set(g, 1)
        Some(c) => local_counts.set(g, c + 1)
      }
    }
    // iLISI = exp(H(p_local))
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
  (ilisi_sum / n.to_double(), kbet_sum / n.to_double())
}

///|
/// Mean iLISI (inverse Simpson-style exponential Shannon entropy) of the
/// `group` distribution in each cell's kNN neighborhood — higher values mean
/// the groups are better mixed.  Works for batch ids (batch mixing) as well
/// as for label ids (label mixing / iCLISI-style).
pub fn scanpy_scib_ilisi(
  embedding : Array[Array[Double]],
  group : Array[Int32],
  k~ : Int32 = 15,
) -> Double {
  scib_knn_mixing(embedding, group, k).0
}

///|
/// Mean kBET-style agreement between each cell's kNN neighborhood
/// composition and the global `group` distribution:
/// `1 - L1(p_neighborhood, p_global) / 2` averaged over cells (1.0 = every
/// neighborhood mirrors the global mix).
pub fn scanpy_scib_kbet(
  embedding : Array[Array[Double]],
  group : Array[Int32],
  k~ : Int32 = 15,
) -> Double {
  scib_knn_mixing(embedding, group, k).1
}

///|
/// Batch PCA entropy: PCA-project the embedding (`n_components` PCs), split
/// each PC's scores into `n_bins` equal-width bins and average the batch
/// distribution entropy inside each bin weighted by bin occupancy, with the
/// per-PC scores combined by explained-variance share, normalized by
/// `ln(n_batches)`.  Returns a value in [0, 1]; 1.0 means the dominant
/// variance directions carry no batch structure (batches perfectly mixed
/// along every informative PC), 0.0 means every PC bin is batch-pure.  A
/// single batch trivially scores 1.0.
pub fn scanpy_scib_pca_entropy(
  embedding : Array[Array[Double]],
  batch : Array[Int32],
  n_components~ : Int32 = 10,
  n_bins~ : Int32 = 4,
) -> Double {
  let n = embedding.length()
  if n == 0 {
    return 0.0
  }
  let batch_ids = Set([])
  for b in batch {
    batch_ids.add(b)
  }
  let n_batches = batch_ids.length()
  if n_batches <= 1 {
    return 1.0
  }
  let n_comps = if n_components < 1 { 1 } else { n_components }
  let pca = scanpy_pca(embedding, n_comps)
  let n_pc = if pca.coords.length() > 0 { pca.coords[0].length() } else { 0 }
  if n_pc == 0 {
    return 0.0
  }
  let ln_b = @math.ln(n_batches.to_double())
  let n_bins_eff = if n_bins < 1 { 1 } else { n_bins }
  // PCs are weighted by their explained-variance share (clipped at 0 so
  // numerically-negative eigenvalues from rank-deficient data contribute
  // nothing), so degenerate noise directions do not dilute the score.
  let weights : Array[Double] = Array(capacity=n_pc)
  let mut w_sum = 0.0
  for c in 0..<n_pc {
    let w = if c < pca.variance_ratio.length() {
      pca.variance_ratio[c].max(0.0)
    } else {
      0.0
    }
    weights.push(w)
    w_sum = w_sum + w
  }
  if w_sum <= 0.0 {
    return 0.0
  }
  let mut acc = 0.0
  for c in 0..<n_pc {
    let mut lo = pca.coords[0][c]
    let mut hi = lo
    for i in 0..<n {
      let s = pca.coords[i][c]
      if s < lo {
        lo = s
      }
      if s > hi {
        hi = s
      }
    }
    let width = (hi - lo) / n_bins_eff.to_double()
    // per-bin batch occupancy
    let bins : Map[Int32, Map[Int32, Int32]] = Map([], capacity=n_bins_eff)
    for i in 0..<n {
      let b = if width <= 0.0 {
        0
      } else {
        let t = ((pca.coords[i][c] - lo) / width).floor().to_int()
        if t < 0 {
          0
        } else if t >= n_bins_eff {
          n_bins_eff - 1
        } else {
          t
        }
      }
      match bins.get(b) {
        None => {
          let m : Map[Int32, Int32] = Map([])
          m.set(batch[i], 1)
          bins.set(b, m)
        }
        Some(m) =>
          match m.get(batch[i]) {
            None => m.set(batch[i], 1)
            Some(cnt) => m.set(batch[i], cnt + 1)
          }
      }
    }
    let mut pc_acc = 0.0
    for kv in bins {
      let mut bin_total = 0
      for kv2 in kv.1 {
        bin_total = bin_total + kv2.1
      }
      if bin_total > 0 {
        let mut h = 0.0
        for kv2 in kv.1 {
          let p = kv2.1.to_double() / bin_total.to_double()
          if p > 0.0 {
            h = h - p * @math.ln(p)
          }
        }
        pc_acc = pc_acc + bin_total.to_double() / n.to_double() * h
      }
    }
    acc = acc + weights[c] / w_sum * pc_acc / ln_b
  }
  acc
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
///   batch (lower is better) and by label (higher is better),
/// - `ari` / `nmi`: agreement of `clusters` (partition of the integrated
///   embedding) with the reference `labels`; 0.0 when `clusters` is empty,
/// - `pca_entropy`: batch mixing entropy along the embedding's principal
///   components (higher is better).
pub struct ScanpyScibMetrics {
  graph_connectivity : Double
  ilisi : Double
  kbet : Double
  silhouette_batch : Double
  silhouette_label : Double
  ari : Double
  nmi : Double
  pca_entropy : Double
}

///|
/// Compute the scIB metrics for an embedding with per-cell batch ids,
/// (optional) cell-label ids and (optional) cluster assignments of the
/// integrated embedding for the ARI/NMI label-agreement scores.  Pass an
/// empty `labels` to skip the label silhouette (reported as 0.0) and an
/// empty `clusters` to skip ARI/NMI (reported as 0.0).
pub fn scanpy_scib_metrics(
  embedding : Array[Array[Double]],
  batch : Array[Int32],
  labels : Array[Int32],
  clusters : Array[Int32],
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
      ari: 0.0,
      nmi: 0.0,
      pca_entropy: 0.0,
    }
  }
  let conn = scanpy_scib_graph_connectivity(embedding, k~)
  let (ilisi, kbet) = scib_knn_mixing(embedding, batch, k)
  let label_sil = if labels.length() == n {
    scanpy_scib_silhouette(embedding, labels)
  } else {
    0.0
  }
  let (ari, nmi) = if clusters.length() == n && labels.length() == n {
    (scanpy_scib_ari(labels, clusters), scanpy_scib_nmi(labels, clusters))
  } else {
    (0.0, 0.0)
  }
  ScanpyScibMetrics::{
    graph_connectivity: conn,
    ilisi,
    kbet,
    silhouette_batch: scanpy_scib_silhouette(embedding, batch),
    silhouette_label: label_sil,
    ari,
    nmi,
    pca_entropy: scanpy_scib_pca_entropy(embedding, batch),
  }
}
