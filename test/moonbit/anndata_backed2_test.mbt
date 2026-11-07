///|
/// Tests for backed-mode enhancements: r+ write-back, pad-aligned reads
/// and backed subsetting.

///|
fn bk_content() -> String {
  let ad = @src.AnnData::new(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]],
    ["c1", "c2", "c3", "c4"],
    ["g1", "g2", "g3"],
  )
  @src.scanpy_write_h5ad(ad)
}

///|
test "backed_read_only_vs_rplus" {
  let ro = @src.anndata_backed_open("f.h5ad", bk_content())
  assert_eq(ro.is_writable(), false)
  assert_eq(ro.mode(), "r")
  // writes are silently ignored in read-only mode
  ro.set_chunk(0, [[99.0, 99.0, 99.0]])
  assert_eq(ro.get_x(0, 0), 1.0)
  assert_eq(ro.flush(), None)

  let rw = @src.anndata_backed_open("f.h5ad", bk_content(), mode="r+")
  assert_eq(rw.is_writable(), true)
  assert_eq(rw.mode(), "r+")
}

///|
test "backed_rplus_writeback_and_flush" {
  let rw = @src.anndata_backed_open("f.h5ad", bk_content(), mode="r+")
  // overwrite rows 1..2
  rw.set_chunk(1, [[40.0, 50.0, 60.0], [70.0, 80.0, 90.0]])
  assert_eq(rw.get_x(1, 0), 40.0)
  assert_eq(rw.get_x(2, 2), 90.0)
  // untouched row unchanged
  assert_eq(rw.get_x(0, 0), 1.0)
  // flush serializes the mutated store
  match rw.flush() {
    Some(content) =>
      match @src.scanpy_read_h5ad(content) {
        Some(reloaded) => {
          assert_eq(reloaded.get_x(1, 0), 40.0)
          assert_eq(reloaded.get_x(0, 0), 1.0)
        }
        None => fail("flushed content unreadable")
      }
    None => fail("flush should return content in r+ mode")
  }
}

///|
test "backed_pad_read_alignment" {
  let bk = @src.anndata_backed_open("f.h5ad", bk_content(), chunk_size=3)
  assert_eq(bk.chunk_size(), 3)
  // request rows [1, 2) — padded to the full first chunk [0, 3)
  let padded = bk.pad_read(1, 2)
  assert_eq(padded.length(), 3)
  assert_eq(padded[0][0], 1.0)
  assert_eq(padded[2][0], 7.0)
  // last chunk padded but clamped to n_obs
  let last = bk.pad_read(3, 4)
  assert_eq(last.length(), 1)
  assert_eq(last[0][2], 12.0)
  // degenerate ranges
  assert_eq(bk.pad_read(0, 0).length(), 0)
  assert_eq(bk.pad_read(5, 9).length(), 0)
}

///|
test "backed_subset_yields_new_handle" {
  let bk = @src.anndata_backed_open("f.h5ad", bk_content())
  let sub = bk.subset([2, 0])
  assert_eq(sub.n_obs(), 2)
  assert_eq(sub.obs_names()[0], "c3")
  assert_eq(sub.obs_names()[1], "c1")
  assert_eq(sub.n_vars(), 3)
  // X resolves through the subset handle
  assert_eq(sub.get_x(0, 0), 7.0)
  assert_eq(sub.get_x(1, 0), 1.0)
  assert_eq(sub.get_row(1)[2], 3.0)
  // filename preserved
  assert_eq(sub.filename(), "f.h5ad")
}

///|
test "backed_mark_loaded" {
  let bk = @src.anndata_backed_open("f.h5ad", bk_content())
  assert_eq(bk.is_loaded(), false)
  ignore(bk.to_memory())
  bk.mark_loaded()
  assert_eq(bk.is_loaded(), true)
}
