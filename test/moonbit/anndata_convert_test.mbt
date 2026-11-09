///|
/// Tests for AnnData <-> SingleCellExperiment / Seurat / SummarizedExperiment
/// conversion bridges.

///|
fn cv_build() -> @src.AnnData {
  @src.AnnData::new([[1.0, 2.0], [3.0, 4.0]], ["c1", "c2"], ["g1", "g2"])
  .set_obs("batch", ["A", "B"])
  .set_var("symbol", ["X1", "Y2"])
  .set_obsm("X_pca", [[1.0, 0.0], [0.0, 1.0]])
  .set_uns("title", @src.JsonStr("demo"))
}

///|
test "convert_anndata_to_sce_and_back" {
  let ad = cv_build()
  let sce = @src.anndata_to_sce(ad)
  // SCE is gene x cells: row_names genes, col_names cells
  assert_eq(sce.row_names[0], "g1")
  assert_eq(sce.col_names[1], "c2")
  match sce.assays.get("counts") {
    Some(counts) => {
      // transposed: counts[gene][cell]
      assert_eq(counts[0][1], 3.0)
      assert_eq(counts[1][0], 2.0)
    }
    None => fail("counts assay missing")
  }
  // obs -> col_data, var -> row_data, obsm -> reduced_dims, uns -> metadata
  match sce.col_data.get("batch") {
    Some(col) => assert_eq(col[0], "A")
    None => fail("col_data missing")
  }
  match sce.row_data.get("symbol") {
    Some(row) => assert_eq(row[1], "Y2")
    None => fail("row_data missing")
  }
  match sce.reduced_dims.get("X_pca") {
    Some(m) => assert_eq(m[1][1], 1.0)
    None => fail("reduced_dims missing")
  }
  match sce.metadata.get("title") {
    Some(s) => assert_eq(s.contains("demo"), true)
    None => fail("metadata missing")
  }
  // round-trip back to AnnData
  let back = @src.sce_to_anndata(sce)
  assert_eq(back.n_obs(), 2)
  assert_eq(back.n_vars(), 2)
  assert_eq(back.get_x(0, 0), 1.0)
  assert_eq(back.get_x(0, 1), 2.0)
  assert_eq(back.get_x(1, 0), 3.0)
  assert_eq(back.obs_names[0], "c1")
  assert_eq(back.var_names[1], "g2")
  match back.obs_column("batch") {
    Some(col) => assert_eq(col[1], "B")
    None => fail("obs missing after roundtrip")
  }
  match back.get_obsm("X_pca") {
    Some(m) => assert_eq(m[0][0], 1.0)
    None => fail("obsm missing after roundtrip")
  }
}

///|
test "convert_anndata_to_seurat_and_back" {
  let ad = cv_build()
  let obj = @src.anndata_to_seurat(ad)
  // Seurat counts are gene x cells (Int32, rounded)
  assert_eq(obj.row_names[0], "g1")
  assert_eq(obj.col_names[1], "c2")
  assert_eq(obj.counts[1][0], 2)
  assert_eq(obj.counts[0][1], 3)
  // data slot filled with the same (double) matrix
  assert_eq(obj.data[1][0], 2.0)
  // obs -> meta_data, X_pca -> pca
  match obj.meta_data.get("batch") {
    Some(col) => assert_eq(col[0], "A")
    None => fail("meta_data missing")
  }
  assert_eq(obj.pca[0][0], 1.0)
  // round-trip back
  let back = @src.seurat_to_anndata(obj)
  assert_eq(back.n_obs(), 2)
  assert_eq(back.n_vars(), 2)
  assert_eq(back.get_x(0, 0), 1.0)
  assert_eq(back.get_x(1, 1), 4.0)
  assert_eq(back.obs_names[0], "c1")
  assert_eq(back.var_names[1], "g2")
  match back.obs_column("batch") {
    Some(col) => assert_eq(col[0], "A")
    None => fail("obs missing after seurat roundtrip")
  }
  match back.get_obsm("X_pca") {
    Some(m) => assert_eq(m[0][1], 0.0)
    None => fail("pca -> obsm missing")
  }
}

///|
test "convert_seurat_counts_fallback_when_data_empty" {
  // a freshly-constructed object has zero data; counts are used instead
  let obj = @src.SeuratObject::new([[1, 2], [3, 4]], ["g1", "g2"], ["c1", "c2"])
  let back = @src.seurat_to_anndata(obj)
  assert_eq(back.get_x(0, 0), 1.0)
  assert_eq(back.get_x(1, 1), 4.0)
  // clusters become an obs column
  match back.obs_column("cluster") {
    Some(col) => {
      assert_eq(col[0], "-1")
      assert_eq(col[1], "-1")
    }
    None => fail("cluster column missing")
  }
}

///|
test "convert_anndata_to_summarized_and_back" {
  let ad = cv_build()
  let se = @src.anndata_to_summarized(ad)
  match @src.se_assay(se, "counts") {
    Some(counts) => {
      assert_eq(counts[0][1], 3.0)
      assert_eq(counts[1][1], 4.0)
    }
    None => fail("counts assay missing")
  }
  assert_eq(@src.se_row_ranges(se)[1].0, "g2")
  assert_eq(@src.se_col_data(se)[0].get("batch"), Some("A"))
  match @src.se_metadata(se).get("title") {
    Some(s) => assert_eq(s.contains("demo"), true)
    None => fail("metadata missing")
  }
  // round-trip back
  match @src.summarized_to_anndata(se) {
    Some(back) => {
      assert_eq(back.n_obs(), 2)
      assert_eq(back.n_vars(), 2)
      assert_eq(back.get_x(1, 1), 4.0)
      assert_eq(back.var_names[0], "g1")
      match back.obs_column("batch") {
        Some(col) => assert_eq(col[1], "B")
        None => fail("obs missing after SE roundtrip")
      }
      match back.get_uns("title") {
        Some(@src.JsonStr(s)) => assert_eq(s, "demo")
        _ => fail("uns missing after SE roundtrip")
      }
    }
    None => fail("SE -> AnnData failed")
  }
}

///|
test "convert_empty_se_returns_none" {
  let empty = @src.summarized_experiment(
    Map([], capacity=0),
    [],
    [],
    Map([], capacity=0),
  )
  match @src.summarized_to_anndata(empty) {
    None => ()
    Some(_) => fail("empty SE should not convert")
  }
}
