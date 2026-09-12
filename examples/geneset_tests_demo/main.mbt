///|
fn fmt4(value : Double) -> String {
  let scaled = (value * 10000.0).round() / 10000.0
  scaled.to_string()
}

///|
fn print_rows(title : String, rows : Array[@src.CameraResult]) -> Unit {
  println("\n\{title}")
  println("  set        N  cor        dir     p-value    FDR")
  for r in rows {
    let cor = match r.correlation() {
      Option::None => "-"
      Some(c) => fmt4(c)
    }
    let fdr = match r.fdr() {
      Option::None => "-"
      Some(f) => fmt4(f)
    }
    println(
      "  \{r.name()}  \{r.n_genes()}  \{cor}  \{r.direction()}  \{fmt4(r.p_value())}  \{fdr}"
    )
  }
}

///|
fn main raise {
  println("=== limma competitive gene set tests Demo ===")

  // ------------------------------------------------------------------
  // 1. rankSumTestWithCorrelation
  //    Wilcoxon-Mann-Whitney rank-sum p-values for a precomputed
  //    gene-level statistic, optionally accounting for correlation
  //    between the genes inside the set.
  // ------------------------------------------------------------------
  let statistics = [
    2.4, 1.7, 1.1, -2.2, -1.8, 0.3, 0.8, -0.5, 0.1, -0.2, 0.9, -0.7,
  ]
  let index = [0, 1, 2] // strongly up-regulated genes (0-based)

  let (down0, up0) = @src.rank_sum_test_with_correlation(index, statistics)
  let (down_c, up_c) = @src.rank_sum_test_with_correlation(
    index,
    statistics,
    correlation=0.2,
  )
  println("\n1. rankSumTestWithCorrelation on a 12-gene statistic vector")
  println("   independent genes : down=\{fmt4(down0)}  up=\{fmt4(up0)}")
  println("   intra-set cor 0.2 : down=\{fmt4(down_c)}  up=\{fmt4(up_c)}")

  // ------------------------------------------------------------------
  // 2. cameraPR: competitive test from precomputed gene statistics.
  // ------------------------------------------------------------------
  let sets = [
    @src.CameraGeneSet::create("UP", [0, 1, 2]),
    @src.CameraGeneSet::create("DOWN", [3, 4]),
    @src.CameraGeneSet::create("NULLSET", [6, 9, 10]),
  ]

  let pr = @src.camera_pr(statistics, sets, sort=false)
  print_rows("2. cameraPR (parametric, fixed inter-gene correlation 0.01)", pr)

  let pr_ranks = @src.camera_pr(
    statistics,
    sets,
    use_ranks=true,
    sort=false,
  )
  print_rows("3. cameraPR with use_ranks=true (Wilcoxon ranks, normal ref.)", pr_ranks)

  // ------------------------------------------------------------------
  // 4. camera.default: fit a linear model per gene first, moderate the
  //    t-statistics with empirical Bayes, then run the competitive test.
  //    12 genes x 6 samples; two groups of three (columns 0-2, 3-5).
  //    Rows 0-2 are induced ~+2, rows 3-4 ~-2, the rest are background.
  // ------------------------------------------------------------------
  let y = @src.Matrix::new(
    12,
    6,
    [
      -0.2, 0.1, 0.0, 2.1, 1.9, 2.0,
      0.3, -0.1, 0.1, 2.2, 1.8, 2.1,
      -0.1, 0.2, -0.2, 1.9, 2.0, 2.2,
      0.0, -0.2, 0.1, -2.0, -1.8, -2.1,
      0.2, 0.0, -0.1, -1.9, -2.1, -2.0,
      0.1, -0.1, 0.1, -0.1, 0.1, -0.1,
      -0.2, 0.2, -0.2, 0.2, -0.2, 0.2,
      0.0, 0.1, -0.1, 0.0, -0.1, 0.1,
      0.15, -0.05, 0.05, -0.15, 0.05, -0.05,
      -0.1, 0.0, 0.1, 0.0, -0.1, 0.1,
      0.05, 0.1, -0.15, -0.05, 0.15, -0.1,
      -0.05, 0.05, 0.0, -0.1, 0.1, 0.0,
    ],
  )
  let design = @src.Matrix::new(
    6,
    2,
    [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
  )

  let cam = @src.camera_default(
    y,
    sets,
    design,
    inter_gene_cor=Some(0.01),
    sort=false,
  )
  print_rows("4. camera.default (fixed inter-gene correlation 0.01)", cam)

  let cam_est = @src.camera_default(
    y,
    sets,
    design,
    inter_gene_cor=Option::None,
    sort=false,
  )
  print_rows("5. camera.default (inter-gene correlation estimated per set)", cam_est)

  println("")
}
