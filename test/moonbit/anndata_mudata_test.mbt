///|
/// Tests for the multimodal MuData container.

///|
fn md_rna() -> @src.AnnData {
  @src.AnnData::new([[1.0, 2.0], [3.0, 4.0]], ["c1", "c2"], ["g1", "g2"]).set_obs(
    "batch",
    ["A", "A"],
  )
}

///|
fn md_adt() -> @src.AnnData {
  @src.AnnData::new([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]], ["c1", "c2"], [
    "p1", "p2", "p3",
  ]).set_obs("batch", ["A", "A"])
}

///|
fn md_mods() -> Map[String, @src.AnnData] {
  let mods : Map[String, @src.AnnData] = Map([], capacity=2)
  mods["rna"] = md_rna()
  mods["adt"] = md_adt()
  mods
}

///|
test "mudata_new_inner_join" {
  let md = @src.mudata_new(md_mods())
  assert_eq(md.n_mods(), 2)
  assert_eq(md.n_obs(), 2)
  assert_eq(md.obs_names[0], "c1")
  match md.mod("rna") {
    Some(rna) => assert_eq(rna.n_vars(), 2)
    None => fail("rna mod missing")
  }
  assert_eq(md.mod_n_vars("adt"), 3)
  assert_eq(md.mod_var_names("adt")[2], "p3")
  match md.mod("nope") {
    None => ()
    Some(_) => fail("unknown mod should be None")
  }
  assert_eq(md.mod_var_names("nope").length(), 0)
}

///|
test "mudata_inner_join_filters_shared_obs" {
  // a mod missing c2 shrinks the global obs set under inner join
  let partial = @src.AnnData::new([[7.0]], ["c1"], ["q1"])
  let mods = md_mods()
  mods["prot"] = partial
  let md = @src.mudata_new(mods)
  assert_eq(md.n_obs(), 1)
  assert_eq(md.obs_names[0], "c1")
}

///|
test "mudata_outer_join" {
  let partial = @src.AnnData::new([[7.0, 8.0]], ["c9"], ["q1", "q2"])
  let mods = md_mods()
  mods["prot"] = partial
  let md = @src.mudata_new(mods, join="outer")
  assert_eq(md.n_obs(), 3)
  // first-seen order: rna c1,c2 then prot c9
  assert_eq(md.obs_names[2], "c9")
}

///|
test "mudata_global_slots" {
  let md = @src.mudata_new(md_mods())
  let pca : Array[Array[Double]] = [[1.0, 2.0], [3.0, 4.0]]
  ignore(md.set_global_obsm("X_pca", pca))
  ignore(md.set_global_varm("PC_loadings", [[0.5]]))
  ignore(md.set_global_uns("title", @src.JsonStr("multi")))
  match md.get_global_obsm("X_pca") {
    Some(m) => assert_eq(m[1][0], 3.0)
    None => fail("global obsm missing")
  }
  match md.get_global_varm("PC_loadings") {
    Some(m) => assert_eq(m[0][0], 0.5)
    None => fail("global varm missing")
  }
  match md.get_global_uns("title") {
    Some(@src.JsonStr(s)) => assert_eq(s, "multi")
    _ => fail("global uns missing")
  }
  assert_eq(md.global_obsm_keys().length(), 1)
  assert_eq(md.uns_keys().length(), 1)
}

///|
test "mudata_to_adata_column_blocks" {
  let md = @src.mudata_new(md_mods())
  let ad = md.to_adata()
  // rna (2 vars) + adt (3 vars) = 5 columns, mod-prefixed
  assert_eq(ad.n_obs(), 2)
  assert_eq(ad.n_vars(), 5)
  assert_eq(ad.var_names[0], "mod:rna:g1")
  assert_eq(ad.var_names[2], "mod:adt:p1")
  assert_eq(ad.var_names[4], "mod:adt:p3")
  assert_eq(ad.get_x(0, 0), 1.0) // rna g1 c1
  assert_eq(ad.get_x(1, 4), 60.0) // adt p3 c2
  // global slots carried over
  ignore(md.set_global_uns("title", @src.JsonStr("multi")))
  let ad2 = @src.mudata_new(md_mods())
    .set_global_uns("title", @src.JsonStr("t"))
    .to_adata()
  match ad2.get_uns("title") {
    Some(@src.JsonStr(s)) => assert_eq(s, "t")
    _ => fail("uns carried over to to_adata failed")
  }
}

///|
test "mudata_json_roundtrip" {
  let md = @src.mudata_new(md_mods()).set_global_uns(
    "title",
    @src.JsonStr("rt"),
  )
  let content = @src.mudata_write_json(md)
  match @src.mudata_read_json(content) {
    Some(back) => {
      assert_eq(back.n_mods(), 2)
      assert_eq(back.n_obs(), 2)
      match back.mod("rna") {
        Some(rna) => {
          assert_eq(rna.n_vars(), 2)
          assert_eq(rna.get_x(0, 0), 1.0)
          match rna.obs_column("batch") {
            Some(col) => assert_eq(col[0], "A")
            None => fail("rna obs missing")
          }
        }
        None => fail("rna mod missing after roundtrip")
      }
      match back.mod("adt") {
        Some(adt) => assert_eq(adt.get_x(1, 2), 60.0)
        None => fail("adt mod missing after roundtrip")
      }
      match back.get_global_uns("title") {
        Some(@src.JsonStr(s)) => assert_eq(s, "rt")
        _ => fail("global uns missing after roundtrip")
      }
    }
    None => fail("mudata json roundtrip failed")
  }
  // malformed input rejected
  match @src.mudata_read_json("garbage") {
    None => ()
    Some(_) => fail("garbage should not parse")
  }
}
