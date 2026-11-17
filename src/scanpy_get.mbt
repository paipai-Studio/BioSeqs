///|
/// scanpy `get` 层 — 与官方 `sc.get.obs_df / sc.get.var_df /
/// sc.get.rank_genes_groups_df` 对齐的查询工具。给定 AnnData 与一组名字，
/// 返回"表达列 + 注释列"对齐的 tidy 数据帧（`ScanpyGetDf`）。

///|
/// Query result: `df` holds the expression columns (row/column names aligned
/// with the AnnData axis), `meta` holds the annotation columns of the opposite
/// axis (`obs` for `obs_df`, `var` for `var_df`).
pub struct ScanpyGetDf {
  df : AnnDataFrame
  meta : Map[String, Array[String]]
}

///|
/// Copy of an obs/var annotation map (shared arrays are read-only in practice).
fn scanpy_get_copy_meta(
  m : Map[String, Array[String]],
) -> Map[String, Array[String]] {
  let keys = m.keys().collect()
  let out : Map[String, Array[String]] = Map([], capacity=keys.length())
  for k in keys {
    out[k] = m[k]
  }
  out
}

///|
/// `sc.get.obs_df` — extract expression columns for the named genes (`genes`;
/// empty means every var) together with all `obs` annotation columns. Rows are
/// `obs_names`, columns are the requested gene names (in order). `layer` selects
/// an alternative matrix layer; genes not present produce zero columns.
pub fn scanpy_get_obs_df(
  adata : AnnData,
  genes~ : Array[String] = [],
  layer~ : String = "",
) -> ScanpyGetDf {
  let wanted : Array[String] = if genes.length() == 0 {
    adata.var_names
  } else {
    genes
  }
  let columns : Map[String, Array[Double]] = Map([], capacity=wanted.length())
  for g in wanted {
    let col : Array[Double] = Array(capacity=adata.n_obs())
    match adata.obs_vector(g, layer~) {
      Some(v) =>
        for i = 0; i < adata.n_obs(); i = i + 1 {
          col.push(if i < v.length() { v[i] } else { 0.0 })
        }
      None =>
        for _i in 0..<adata.n_obs() {
          col.push(0.0)
        }
    }
    columns[g] = col
  }
  ScanpyGetDf::{
    df: AnnDataFrame::{
      columns,
      row_names: adata.obs_names,
      col_names: wanted,
    },
    meta: scanpy_get_copy_meta(adata.obs),
  }
}

///|
/// `sc.get.var_df` — extract expression columns for the named cells (`cells`;
/// empty means every obs) together with all `var` annotation columns. Rows are
/// `var_names`, columns are the requested cell names (in order). `layer` selects
/// an alternative matrix layer; cells not present produce zero columns.
pub fn scanpy_get_var_df(
  adata : AnnData,
  cells~ : Array[String] = [],
  layer~ : String = "",
) -> ScanpyGetDf {
  let wanted : Array[String] = if cells.length() == 0 {
    adata.obs_names
  } else {
    cells
  }
  let columns : Map[String, Array[Double]] = Map([], capacity=wanted.length())
  for c in wanted {
    let col : Array[Double] = Array(capacity=adata.n_vars())
    match adata.var_vector(c, layer~) {
      Some(v) =>
        for j = 0; j < adata.n_vars(); j = j + 1 {
          col.push(if j < v.length() { v[j] } else { 0.0 })
        }
      None =>
        for _j in 0..<adata.n_vars() {
          col.push(0.0)
        }
    }
    columns[c] = col
  }
  ScanpyGetDf::{
    df: AnnDataFrame::{
      columns,
      row_names: adata.var_names,
      col_names: wanted,
    },
    meta: scanpy_get_copy_meta(adata.variables),
  }
}

///|
/// Build the canonical `uns["rank_genes_groups"]` JSON payload from per-group
/// rank results (as produced by `scanpy_rank_genes_groups` over each group
/// pair). Layout: `{"groups": [...], "names": [[...], ...], "scores": ...,
/// "pvals": ..., "pvals_adj": ..., "logfoldchanges": ...}` where every matrix
/// is `[group][gene]`.
pub fn scanpy_rank_genes_groups_json(
  group_names : Array[String],
  gene_names : Array[String],
  scores : Array[Array[Double]],
  pvals : Array[Array[Double]],
  pvals_adj : Array[Array[Double]],
  logfoldchanges : Array[Array[Double]],
) -> JsonValue {
  let strs = fn(a : Array[String]) -> JsonValue {
    JsonValue::JsonArr(a.map(fn(s) { JsonValue::JsonStr(s) }))
  }
  let nums = fn(a : Array[Array[Double]]) -> JsonValue {
    JsonValue::JsonArr(
      a.map(fn(row) {
        JsonValue::JsonArr(row.map(fn(v) { JsonValue::JsonNum(v) }))
      }),
    )
  }
  let obj : Map[String, JsonValue] = Map([], capacity=6)
  obj["groups"] = strs(group_names)
  obj["names"] = strs(gene_names)
  obj["scores"] = nums(scores)
  obj["pvals"] = nums(pvals)
  obj["pvals_adj"] = nums(pvals_adj)
  obj["logfoldchanges"] = nums(logfoldchanges)
  JsonValue::JsonObj(obj)
}

///|
/// `sc.get.rank_genes_groups_df` — read the `uns[uns_key]` payload (built by
/// `scanpy_rank_genes_groups_json` / `set_uns`) back as a tidy frame for the
/// requested `group` (default: first group), keeping genes whose adjusted
/// p-value passes `pval_cutoff` and |log fold change| passes `logfc_cutoff`.
/// Columns are `score`, `p_value`, `p_value_adj`, `logfoldchange`; rows are the
/// gene names. Returns `None` when the key is absent or malformed.
pub fn scanpy_get_rank_genes_groups_df(
  adata : AnnData,
  group~ : String = "",
  uns_key~ : String = "rank_genes_groups",
  pval_cutoff~ : Double = 1.0,
  logfc_cutoff~ : Double = 0.0,
) -> ScanpyGetDf? {
  let payload = match adata.uns.get(uns_key) {
    Some(JsonValue::JsonObj(o)) => o
    _ => return None
  }
  let names = match payload.get("names") {
    Some(JsonValue::JsonArr(a)) => {
      let out : Array[String] = Array(capacity=0)
      for v in a {
        match v {
          JsonValue::JsonStr(s) => out.push(s)
          _ => return None
        }
      }
      out
    }
    _ => return None
  }
  let groups = match payload.get("groups") {
    Some(JsonValue::JsonArr(a)) => {
      let out : Array[String] = Array(capacity=0)
      for v in a {
        match v {
          JsonValue::JsonStr(s) => out.push(s)
          _ => out.push("")
        }
      }
      out
    }
    _ => Array::make(1, "")
  }
  let mut gi = 0
  if group.length() > 0 {
    let mut found = -1
    for i = 0; i < groups.length(); i = i + 1 {
      if groups[i] == group {
        found = i
        break
      }
    }
    if found < 0 {
      return None
    }
    gi = found
  }
  let get_row = fn(key : String) -> Array[Double]? {
    match payload.get(key) {
      Some(JsonValue::JsonArr(a)) =>
        match a.get(gi) {
          Some(JsonValue::JsonArr(row)) => {
            let out : Array[Double] = Array(capacity=0)
            for v in row {
              match v {
                JsonValue::JsonNum(n) => out.push(n)
                _ => return None
              }
            }
            Some(out)
          }
          _ => None
        }
      _ => None
    }
  }
  let scores = match get_row("scores") {
    Some(v) => v
    None => return None
  }
  let pvals = match get_row("pvals") {
    Some(v) => v
    None => return None
  }
  let pvals_adj = match get_row("pvals_adj") {
    Some(v) => v
    None => Array::make(scores.length(), 1.0)
  }
  let lfcs = match get_row("logfoldchanges") {
    Some(v) => v
    None => Array::make(scores.length(), 0.0)
  }
  let columns : Map[String, Array[Double]] = Map([], capacity=4)
  let score_col : Array[Double] = Array(capacity=0)
  let pval_col : Array[Double] = Array(capacity=0)
  let padj_col : Array[Double] = Array(capacity=0)
  let lfc_col : Array[Double] = Array(capacity=0)
  let rows : Array[String] = Array(capacity=0)
  let n = if names.length() < scores.length() {
    names.length()
  } else {
    scores.length()
  }
  for g = 0; g < n; g = g + 1 {
    let lfc = if g < lfcs.length() { lfcs[g] } else { 0.0 }
    let abs_lfc = if lfc < 0.0 { -lfc } else { lfc }
    let pa = if g < pvals_adj.length() { pvals_adj[g] } else { 1.0 }
    if pa <= pval_cutoff && abs_lfc >= logfc_cutoff {
      rows.push(names[g])
      score_col.push(scores[g])
      pval_col.push(if g < pvals.length() { pvals[g] } else { 1.0 })
      padj_col.push(pa)
      lfc_col.push(lfc)
    }
  }
  columns["score"] = score_col
  columns["p_value"] = pval_col
  columns["p_value_adj"] = padj_col
  columns["logfoldchange"] = lfc_col
  Some(ScanpyGetDf::{
    df: AnnDataFrame::{
      columns,
      row_names: rows,
      col_names: ["score", "p_value", "p_value_adj", "logfoldchange"],
    },
    meta: Map([], capacity=0),
  })
}
