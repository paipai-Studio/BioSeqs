///|
/// scanpy single-cell end-to-end demo: a tiny synthetic 10x dataset (2 clusters
/// of 3 cells, 2 marker genes each) → read_10x_mtx → normalize_total → log1p →
/// pca → neighbors → louvain / leiden → diffmap → dpt → umap → paga. Every step
/// maps to the corresponding `sc.pp.*` / `sc.tl.*` call.

///|
fn fmt_d(row : Array[Double]) -> String {
  let sb = StringBuilder::new()
  let mut first = true
  for x in row {
    if not(first) {
      sb.write_char(' ')
    }
    sb.write_string(x.to_string())
    first = false
  }
  sb.to_string()
}

///|
fn fmt_i(row : Array[Int]) -> String {
  let sb = StringBuilder::new()
  let mut first = true
  for x in row {
    if not(first) {
      sb.write_char(' ')
    }
    sb.write_string(x.to_string())
    first = false
  }
  sb.to_string()
}

///|
fn fmt_m(m : Array[Array[Double]]) -> String {
  let sb = StringBuilder::new()
  for row in m {
    sb.write_string("    ")
    sb.write_string(fmt_d(row))
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
fn main {
  // 4 genes x 6 cells (cells 1-3 = cluster A, cells 4-6 = cluster B)
  let mtx = "%%MatrixMarket matrix coordinate integer general\n%\n4 6 20\n1 1 10\n1 2 11\n1 3 9\n1 4 1\n1 6 1\n2 1 8\n2 2 9\n2 3 10\n2 4 2\n2 5 1\n3 2 1\n3 3 2\n3 4 10\n3 5 12\n3 6 11\n4 1 1\n4 3 1\n4 4 9\n4 5 10\n4 6 12\n"
  let features = "g1\tg1\tGene Expression\ng2\tg2\tGene Expression\ng3\tg3\tGene Expression\ng4\tg4\tGene Expression\n"
  let barcodes = "cell1\ncell2\ncell3\ncell4\ncell5\ncell6\n"

  // sc.read_10x_mtx (mtx is genes×cells here; transposed to cells×genes)
  let adata = @src.scanpy_read_10x_mtx(mtx, features, barcodes, "gene_symbols", true, true)
  println("read_10x_mtx -> AnnData \{adata.n_obs()} cells x \{adata.n_vars()} genes")

  // sc.pp.normalize_total + sc.pp.log1p
  let norm = @src.scanpy_normalize_total(adata.x_matrix(), 10000.0)
  let logx = @src.scanpy_log1p(norm)

  // sc.pp.pca (2 components)
  let pca = @src.scanpy_pca(logx, 2)
  println("\npp.pca variance_ratio = \{fmt_d(pca.variance_ratio)}")
  println("pp.pca coordinates (cells x PCs):\n\{fmt_m(pca.coords)}")

  // sc.pp.neighbors (on PCA coords, 3 neighbors) -> connectivities
  let conn = @src.scanpy_neighbors(pca.coords, 3)

  // sc.tl.louvain / sc.tl.leiden
  let louvain = @src.scanpy_louvain(conn)
  let leiden = @src.scanpy_leiden(conn)
  println("\ntl.louvain labels = \{fmt_i(louvain)}")
  println("tl.leiden  labels = \{fmt_i(leiden)}")

  // sc.tl.diffmap + sc.tl.dpt (root = first cell)
  let dm = @src.scanpy_diffmap(conn, 6)
  println("\ntl.diffmap eigenvalues = \{fmt_d(dm.eigen_values)}")
  let dpt = @src.scanpy_dpt(dm, 0)
  println("tl.dpt pseudotime (root cell1) = \{fmt_d(dpt)}")

  // sc.tl.umap (2D)
  let emb = @src.scanpy_umap(conn, 200)
  println("\ntl.umap 2D coordinates:\n\{fmt_m(emb)}")

  // sc.tl.paga (cluster connectivity, 2 x 2)
  let paga = @src.scanpy_paga(conn, louvain)
  println("tl.paga cluster-connectivity matrix (2x2):")
  println(fmt_m(paga))
}