/// anndata_convert.mbt — bridges between AnnData and the R-flavoured
/// containers already in BioSeqs: SingleCellExperiment, SummarizedExperiment
/// and SeuratObject (the anndata2ri / Seurat::CreateSeuratObject analogues).
///
/// Orientation conventions: AnnData X is cells × genes; SCE/SE/Seurat
/// assay matrices are genes × cells, so conversions transpose.

///|
/// Transpose a dense matrix (cells×genes -> genes×cells and back).
fn convert_transpose(x : Array[Array[Double]]) -> Array[Array[Double]] {
  let n_rows = x.length()
  let n_cols = if n_rows > 0 { x[0].length() } else { 0 }
  let out : Array[Array[Double]] = Array(capacity=n_cols)
  for j = 0; j < n_cols; j = j + 1 {
    let row : Array[Double] = Array(capacity=n_rows)
    for i = 0; i < n_rows; i = i + 1 {
      if j < x[i].length() {
        row.push(x[i][j])
      } else {
        row.push(0.0)
      }
    }
    out.push(row)
  }
  out
}

///|
/// Stringify an uns JsonValue for string-map style metadata.
fn convert_uns_strs(uns : Map[String, JsonValue]) -> Map[String, String] {
  let out : Map[String, String] = Map([], capacity=4)
  for k in uns.keys().collect() {
    out[k] = json_to_string(uns[k])
  }
  out
}

///|
/// Parse string-map style metadata back into uns (JSON when possible).
fn convert_strs_uns(m : Map[String, String]) -> Map[String, JsonValue] {
  let out : Map[String, JsonValue] = Map([], capacity=4)
  for k in m.keys().collect() {
    match parse_json(m[k]) {
      Some(v) => out[k] = v
      None => out[k] = JsonStr(m[k])
    }
  }
  out
}

///|
/// Convert an AnnData to a SingleCellExperiment: X (transposed) becomes the
/// "counts" assay, var -> row_data, obs -> col_data, obsm -> reduced_dims,
/// uns -> metadata.
pub fn anndata_to_sce(adata : AnnData) -> SingleCellExperiment {
  let sce = SingleCellExperiment::new(
    convert_transpose(adata.x_matrix()),
    adata.var_names.copy(),
    adata.obs_names.copy(),
  )
  for k in adata.variables.keys().collect() {
    sce.row_data[k] = adata.variables[k]
  }
  for k in adata.obs.keys().collect() {
    sce.col_data[k] = adata.obs[k]
  }
  for k in adata.obsm.keys().collect() {
    sce.reduced_dims[k] = adata.obsm[k]
  }
  for k in adata.uns.keys().collect() {
    sce.metadata[k] = json_to_string(adata.uns[k])
  }
  sce
}

///|
/// Convert a SingleCellExperiment back to an AnnData: the "counts" assay
/// (transposed) becomes X, row_data -> var, col_data -> obs,
/// reduced_dims -> obsm, metadata -> uns.
pub fn sce_to_anndata(sce : SingleCellExperiment) -> AnnData {
  let counts = match sce.assays.get("counts") {
    Some(m) => m
    None => Array(capacity=0)
  }
  let adata = AnnData::new(
    convert_transpose(counts),
    sce.col_names.copy(),
    sce.row_names.copy(),
  )
  for k in sce.row_data.keys().collect() {
    adata.variables[k] = sce.row_data[k]
  }
  for k in sce.col_data.keys().collect() {
    adata.obs[k] = sce.col_data[k]
  }
  for k in sce.reduced_dims.keys().collect() {
    adata.obsm[k] = sce.reduced_dims[k]
  }
  let meta = convert_strs_uns(sce.metadata)
  for k in meta.keys().collect() {
    adata.uns[k] = meta[k]
  }
  adata
}

///|
/// Convert an AnnData to a SeuratObject: X (transposed, rounded) becomes
/// counts, obs -> meta_data, obsm["X_pca"] becomes pca (cell × dim),
/// obsm["X_umap"] becomes umap. The log-normalized `data` slot is filled
/// with the same matrix as counts (as doubles) row-by-row.
pub fn anndata_to_seurat(adata : AnnData) -> SeuratObject {
  let genes_x_cells = convert_transpose(adata.x_matrix())
  let counts : Array[Array[Int32]] = Array(capacity=genes_x_cells.length())
  for row in genes_x_cells {
    let crow : Array[Int32] = Array(capacity=row.length())
    for v in row {
      let iv = if v >= 0.0 { v.to_int() } else { 0 }
      crow.push(iv)
    }
    counts.push(crow)
  }
  let obj = SeuratObject::new(
    counts,
    adata.var_names.copy(),
    adata.obs_names.copy(),
  )
  for i = 0; i < genes_x_cells.length() && i < obj.data.length(); i = i + 1 {
    let row : Array[Double] = Array(capacity=genes_x_cells[i].length())
    for v in genes_x_cells[i] {
      row.push(v)
    }
    obj.data[i] = row
  }
  for k in adata.obs.keys().collect() {
    obj.meta_data[k] = adata.obs[k]
  }
  match adata.obsm.get("X_pca") {
    Some(pca) => obj.pca = pca
    None => ()
  }
  match adata.obsm.get("X_umap") {
    Some(umap) => obj.umap = umap
    None => ()
  }
  obj
}

///|
/// Convert a SeuratObject back to an AnnData: `data` (transposed) becomes X,
/// falling back to counts when data is empty/all-zero (as left by
/// SeuratObject::new), meta_data -> obs, clusters -> an obs column, pca
/// (cell × dim) -> obsm["X_pca"], umap -> obsm["X_umap"], var_features
/// recorded in uns.
pub fn seurat_to_anndata(obj : SeuratObject) -> AnnData {
  let n_genes = obj.row_names.length()
  let n_cells = obj.col_names.length()
  // data is usable when it matches the shape and holds any non-zero entry;
  // an all-zero data block is SeuratObject::new's placeholder
  let mut data_has_values = false
  if obj.data.length() == n_genes &&
    n_genes > 0 &&
    obj.data[0].length() == n_cells &&
    n_cells > 0 {
    for row in obj.data {
      for v in row {
        if v != 0.0 {
          data_has_values = true
          break
        }
      }
      if data_has_values {
        break
      }
    }
  }
  let src_is_data = data_has_values
  let src = if src_is_data { obj.data } else { [] }
  let x : Array[Array[Double]] = Array(capacity=n_cells)
  for c = 0; c < n_cells; c = c + 1 {
    let row : Array[Double] = Array(capacity=n_genes)
    for g = 0; g < n_genes; g = g + 1 {
      if src_is_data && g < src.length() && c < src[g].length() {
        row.push(src[g][c])
      } else if not(src_is_data) &&
        g < obj.counts.length() &&
        c < obj.counts[g].length() {
        row.push(obj.counts[g][c].to_double())
      } else {
        row.push(0.0)
      }
    }
    x.push(row)
  }
  let adata = AnnData::new(x, obj.col_names.copy(), obj.row_names.copy())
  for k in obj.meta_data.keys().collect() {
    adata.obs[k] = obj.meta_data[k]
  }
  // clusters -> obs column (as strings)
  if obj.clusters.length() == n_cells && n_cells > 0 {
    let cl : Array[String] = Array(capacity=n_cells)
    for c in obj.clusters {
      cl.push(c.to_string())
    }
    adata.obs["cluster"] = cl
  }
  if obj.pca.length() == n_cells && n_cells > 0 {
    adata.obsm["X_pca"] = obj.pca
  }
  if obj.umap.length() == n_cells && n_cells > 0 {
    adata.obsm["X_umap"] = obj.umap
  }
  if obj.var_features.length() > 0 {
    let vf : Array[JsonValue] = Array(capacity=obj.var_features.length())
    for f in obj.var_features {
      vf.push(JsonStr(f))
    }
    adata.uns["var_features"] = JsonArr(vf)
  }
  adata
}

///|
/// Convert an AnnData to a SummarizedExperiment: X (transposed) becomes the
/// "counts" assay, var_names -> row_ranges ((name, 0, 0) triples),
/// obs -> col_data (one map per cell), uns -> metadata (stringified).
pub fn anndata_to_summarized(adata : AnnData) -> SummarizedExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([], capacity=2)
  assays["counts"] = convert_transpose(adata.x_matrix())
  let row_ranges : Array[(String, Int32, Int32)] = Array(
    capacity=adata.var_names.length(),
  )
  for n in adata.var_names {
    row_ranges.push((n, 0, 0))
  }
  let col_data : Array[Map[String, String]] = Array(
    capacity=adata.obs_names.length(),
  )
  for _i = 0; _i < adata.obs_names.length(); _i = _i + 1 {
    col_data.push(Map([], capacity=4))
  }
  for k in adata.obs.keys().collect() {
    let vals = adata.obs[k]
    for i = 0; i < vals.length() && i < col_data.length(); i = i + 1 {
      col_data[i][k] = vals[i]
    }
  }
  summarized_experiment(
    assays,
    row_ranges,
    col_data,
    convert_uns_strs(adata.uns),
  )
}

///|
/// Convert a SummarizedExperiment back to an AnnData: the first assay
/// (transposed) becomes X, row_ranges names -> var_names, col_data maps ->
/// obs columns, metadata -> uns. Returns `None` when the SE has no assay.
pub fn summarized_to_anndata(se : SummarizedExperiment) -> AnnData? {
  let assay_names = se.assays.keys().collect()
  if assay_names.length() == 0 {
    return None
  }
  let counts = se.assays[assay_names[0]]
  let n_genes = se.row_ranges.length()
  let n_cells = se.col_data.length()
  let var_names : Array[String] = Array(capacity=n_genes)
  for r in se.row_ranges {
    var_names.push(r.0)
  }
  let x = convert_transpose(counts)
  let obs_names : Array[String] = Array(capacity=n_cells)
  for i = 0; i < n_cells; i = i + 1 {
    obs_names.push("cell_" + i.to_string())
  }
  let adata = AnnData::new(x, obs_names, var_names)
  // col_data maps -> obs columns (keys from the first cell's map, falling
  // back to any key seen)
  let keys : Array[String] = Array(capacity=0)
  let seen : Set[String] = Set([])
  for m in se.col_data {
    for k in m.keys().collect() {
      if not(seen.contains(k)) {
        seen.add(k)
        keys.push(k)
      }
    }
  }
  for k in keys {
    let col : Array[String] = Array(capacity=n_cells)
    for m in se.col_data {
      match m.get(k) {
        Some(v) => col.push(v)
        None => col.push("")
      }
    }
    adata.obs[k] = col
  }
  for k in se.metadata.keys().collect() {
    match parse_json(se.metadata[k]) {
      Some(v) => adata.uns[k] = v
      None => adata.uns[k] = JsonStr(se.metadata[k])
    }
  }
  Some(adata)
}
