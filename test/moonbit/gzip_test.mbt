///|
/// Tests for gzip decompression and the gzipped 10x-mtx path.

///|
fn bytes_str(bytes : Array[Int]) -> String {
  let sb = StringBuilder()
  for b in bytes {
    sb.write_char(b.to_uint16().unsafe_to_char())
  }
  sb.to_string()
}

///|
test "gzip_decompress" {
  let blob = bytes_str([
    31, 139, 8, 0, 84, 68, 179, 106, 2, 255, 115, 116, 118, 15, 81, 200, 72, 205,
    201, 201, 87, 40, 78, 78, 204, 43, 168, 84, 48, 52, 50, 54, 49, 5, 0, 86, 65,
    141, 87, 23, 0, 0, 0,
  ])
  let out = @src.gzip_decompress(blob)
  assert_eq(out, "ACGT hello scanpy 12345")
  // gunzip_if_needed leaves plain text unchanged
  assert_eq(@src.gunzip_if_needed("plain text"), "plain text")
  assert_eq(@src.gunzip_if_needed(blob), "ACGT hello scanpy 12345")
}

///|
test "scanpy_read_10x_mtx_gz" {
  // gzipped "%%MatrixMarket matrix coordinate integer general\n3 2 4\n1 1 10\n1 2 5\n2 1 3\n3 2 7\n"
  let mtx_gz = bytes_str([
    31, 139, 8, 0, 165, 68, 179, 106, 2, 255, 83, 85, 245, 77, 44, 41, 202, 172,
    240, 77, 44, 202, 78, 45, 81, 200, 5, 115, 20, 146, 243, 243, 139, 82, 50, 243,
    18, 75, 82, 21, 50, 243, 74, 82, 211, 83, 139, 20, 210, 83, 243, 82, 139, 18,
    115, 184, 140, 21, 140, 20, 76, 184, 12, 21, 128, 208, 0, 72, 25, 41, 152, 114,
    25, 1, 57, 198, 96, 9, 115, 46, 0, 16, 98, 85, 94, 80, 0, 0, 0,
  ])
  let features = "g1\tsym1\tGene Expression\ng2\tsym2\tGene Expression\ng3\tsym3\tAntibody Capture\n"
  let barcodes = "cell1\ncell2\n"
  let adata = @src.scanpy_read_10x_mtx(mtx_gz, features, barcodes, "gene_symbols", true, true)
  assert_eq(adata.n_obs(), 2)
  assert_eq(adata.n_vars(), 2)
  assert_true(assert_close(adata.x_matrix()[0][0], 10.0, 1.0e-9))
  assert_true(assert_close(adata.x_matrix()[1][0], 5.0, 1.0e-9))
}