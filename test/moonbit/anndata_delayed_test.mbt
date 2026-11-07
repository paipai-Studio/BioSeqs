///|
/// Tests for chunked lazy evaluation over AnnData X (dask analogue).

///|
fn dl_build() -> @src.AnnData {
  // 4 obs x 2 vars, chunked by 2 rows
  @src.AnnData::new(
    [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]],
    ["c1", "c2", "c3", "c4"],
    ["g1", "g2"],
  )
}

///|
test "delayed_new_and_chunks" {
  let d = @src.anndata_delayed_new(dl_build(), chunk_size=2)
  assert_eq(@src.anndata_delayed_n_chunks(d), 2)
  // identity: chunk 0 matches the source rows
  let c0 = @src.anndata_delayed_chunk(d, 0)
  assert_eq(c0.length(), 2)
  assert_eq(c0[0][0], 1.0)
  assert_eq(c0[1][1], 4.0)
  let c1 = @src.anndata_delayed_chunk(d, 1)
  assert_eq(c1[0][0], 5.0)
  // ragged final chunk: chunk size 3 -> chunks of 3 + 1 rows
  let d3 = @src.anndata_delayed_new(dl_build(), chunk_size=3)
  assert_eq(@src.anndata_delayed_n_chunks(d3), 2)
  assert_eq(@src.anndata_delayed_chunk(d3, 1).length(), 1)
  // out-of-range chunk is empty
  assert_eq(@src.anndata_delayed_chunk(d, 5).length(), 0)
}

///|
test "delayed_ops_chain" {
  let d = @src.anndata_delayed_new(dl_build(), chunk_size=2)
  let dlog = @src.anndata_delayed_map(d, @src.DLog1p)
  // log1p(1) = ln(2)
  let c = @src.anndata_delayed_chunk(dlog, 0)
  let diff = c[0][0] - @math.ln(2.0)
  assert_eq(diff < 1.0e-12 && diff > -1.0e-12, true)
  // chaining: log1p then +1 then *2 then square
  let chain = @src.anndata_delayed_map(
    @src.anndata_delayed_map(
      @src.anndata_delayed_map(dlog, @src.DAddScalar(1.0)),
      @src.DMulScalar(2.0),
    ),
    @src.DSquare,
  )
  let cc = @src.anndata_delayed_elem(chain, 0, 0)
  let expect = {
    let l = @math.ln(2.0)
    let s = (l + 1.0) * 2.0
    s * s
  }
  let delta = cc - expect
  assert_eq(delta < 1.0e-9 && delta > -1.0e-9, true)
  // original handle unaffected by the chain
  assert_eq(@src.anndata_delayed_elem(d, 0, 0), 1.0)
}

///|
test "delayed_collect_and_to_adata" {
  let d = @src.anndata_delayed_map(
    @src.anndata_delayed_new(dl_build(), chunk_size=2),
    @src.DMulScalar(10.0),
  )
  let full = @src.anndata_delayed_collect(d)
  assert_eq(full.length(), 4)
  assert_eq(full[3][1], 80.0)
  let ad = @src.anndata_delayed_to_adata(d)
  assert_eq(ad.n_obs(), 4)
  assert_eq(ad.get_x(3, 1), 80.0)
  assert_eq(ad.obs_names[0], "c1")
  assert_eq(ad.var_names[1], "g2")
}

///|
test "delayed_streaming_sums" {
  let d = @src.anndata_delayed_map(
    @src.anndata_delayed_new(dl_build(), chunk_size=2),
    @src.DAddScalar(1.0),
  )
  // row sums with +1 per element: 5, 9, 13, 17
  let rs = @src.anndata_delayed_row_sums(d)
  assert_eq(rs.length(), 4)
  assert_eq(rs[0], 5.0)
  assert_eq(rs[3], 17.0)
  // col sums: (1+3+5+7)+4=20, (2+4+6+8)+4=24
  let cs = @src.anndata_delayed_col_sums(d)
  assert_eq(cs[0], 20.0)
  assert_eq(cs[1], 24.0)
}

///|
test "delayed_read_elem" {
  let elem = @src.JsonArr([
    @src.JsonArr([@src.JsonNum(1.0), @src.JsonNum(2.0)]),
    @src.JsonArr([@src.JsonNum(3.0), @src.JsonNum(4.0)]),
  ])
  match @src.anndata_read_elem_as_delayed(elem, chunk_size=1) {
    Some(d) => {
      assert_eq(@src.anndata_delayed_n_chunks(d), 2)
      assert_eq(@src.anndata_delayed_elem(d, 1, 1), 4.0)
      // with an op applied
      let dl = @src.anndata_delayed_map(d, @src.DMulScalar(3.0))
      assert_eq(@src.anndata_delayed_elem(dl, 0, 1), 6.0)
    }
    None => fail("read_elem_as_delayed failed")
  }
  // non-array element rejected
  match @src.anndata_read_elem_as_delayed(@src.JsonStr("nope")) {
    None => ()
    Some(_) => fail("should reject non-array element")
  }
}
