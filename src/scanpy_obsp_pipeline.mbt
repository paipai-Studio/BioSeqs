///|
/// scanpy neighbors -> `adata.obsp` -> clustering pipeline.  In scanpy,
/// `sc.pp.neighbors` writes `obsp["distances"]` / `obsp["connectivities"]` and
/// the `tl.*` graph tools read them back.  This file wires the existing
/// matrix-level scanpy functions into the `AnnData` container:
/// `scanpy_neighbors_anndata` computes the kNN graphs and stores them as
/// dense pairwise matrices in `obsp`; `scanpy_graph_from_obsp` reads a graph
/// back as a `SparseGraph`; `scanpy_louvain_anndata` / `scanpy_leiden_anndata`
/// cluster directly from `obsp` and write the labels to `obs`.

///|
/// Result of `scanpy_neighbors_anndata`: the `obsp` keys that were written
/// and the number of observations.
pub struct ScanpyObspNeighbors {
  n_obs : Int32
  distances_key : String
  connectivities_key : String
}

///|
/// Compute neighbors on the AnnData expression matrix and store the two
/// sparse graphs as dense pairwise matrices under `adata.obsp[distances_key]`
/// (symmetrized kNN distances, 0 diagonal) and
/// `adata.obsp[connectivities_key]` (UMAP fuzzy simplicial set with
/// `method="umap"`, or Gaussian kernel with `method="gauss"`; diagonal 1).
/// Mirrors `sc.pp.neighbors(adata, n_neighbors, metric, method)`.
pub fn scanpy_neighbors_anndata(
  adata : AnnData,
  n_neighbors : Int32,
  metric~ : String = "euclidean",
  method_~ : String = "umap",
  distances_key~ : String = "distances",
  connectivities_key~ : String = "connectivities",
) -> ScanpyObspNeighbors {
  let x = adata.x_matrix()
  let res = scanpy_neighbors_metric(x, n_neighbors, metric~, method_~)
  adata.obsp[distances_key] = res.distances.to_dense()
  adata.obsp[connectivities_key] = res.connectivities.to_dense()
  ScanpyObspNeighbors::{
    n_obs: adata.n_obs(),
    distances_key,
    connectivities_key,
  }
}

///|
/// Read a pairwise matrix from `adata.obsp` as a `SparseGraph` (scanpy's
/// `adata.obsp[key]` adjacency).  The diagonal becomes the graph self-weight
/// and every strictly positive off-diagonal entry becomes an undirected edge.
/// Returns an empty graph when the key is absent.
pub fn scanpy_graph_from_obsp(
  adata : AnnData,
  key~ : String = "connectivities",
) -> SparseGraph {
  let m = match adata.obsp.get(key) {
    Some(m) => m
    None => return SparseGraph::{ n: 0, diag: [], adj: [], }
  }
  let n = m.length()
  let diag : Array[Double] = Array(capacity=n)
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    diag.push(if i < m[i].length() { m[i][i] } else { 0.0 })
    let row : Array[(Int32, Double)] = Array(capacity=0)
    for j = 0; j < m[i].length(); j = j + 1 {
      if j != i && m[i][j] > 0.0 {
        row.push((j, m[i][j]))
      }
    }
    adj.push(row)
  }
  SparseGraph::{ n, diag, adj, }
}

///|
/// Run Louvain community detection on `adata.obsp[obsp_key]` and store the
/// labels (as strings, dense `0..k-1`) under `adata.obs[obs_key]`.  Returns
/// the number of communities.
pub fn scanpy_louvain_anndata(
  adata : AnnData,
  obsp_key~ : String = "connectivities",
  obs_key~ : String = "louvain",
) -> Int32 {
  let graph = scanpy_graph_from_obsp(adata, key=obsp_key)
  let labels = scanpy_louvain(graph)
  let strs : Array[String] = Array(capacity=labels.length())
  for l in labels {
    strs.push(l.to_string())
  }
  adata.obs[obs_key] = strs
  let mut mx = 0
  for l in labels {
    if l > mx {
      mx = l
    }
  }
  mx + 1
}

///|
/// Run Leiden-style clustering on `adata.obsp[obsp_key]` and store the labels
/// (as strings, dense `0..k-1`) under `adata.obs[obs_key]`.  Returns the
/// number of communities.
pub fn scanpy_leiden_anndata(
  adata : AnnData,
  obsp_key~ : String = "connectivities",
  obs_key~ : String = "leiden",
) -> Int32 {
  let graph = scanpy_graph_from_obsp(adata, key=obsp_key)
  let labels = scanpy_leiden(graph)
  let strs : Array[String] = Array(capacity=labels.length())
  for l in labels {
    strs.push(l.to_string())
  }
  adata.obs[obs_key] = strs
  let mut mx = 0
  for l in labels {
    if l > mx {
      mx = l
    }
  }
  mx + 1
}
