///|
/// RankProd Module Demo.
///
/// Port of the Bioconductor RankProd package: rank-product / rank-sum
/// differential expression with exact analytical p-values and predicted
/// false-positive percentages (pfp, the rank-product FDR analogue).
///
/// A synthetic log2 expression matrix (12 genes x 6 samples, two classes of
/// three) carries four planted genes:
///   - genes 0-1: about 8-fold DOWN in class 1 (class1 < class2)
///   - genes 2-3: about 8-fold UP in class 1 (class1 > class2)
/// The other eight genes are balanced nulls.

/// Round to four significant digits, mirroring R prettyNum(digits = 4).
fn sig4(x : Double) -> String {
  if x.is_nan() {
    return "NaN"
  }
  if x == 0.0 {
    return "0"
  }
  let ax = x.abs()
  let mag = @math.floor(@math.log10(ax)).to_int()
  let factor = @math.pow(10.0, (4 - 1 - mag).to_double())
  let scaled = x * factor
  let rounded = if scaled >= 0.0 {
    @math.floor(scaled + 0.5)
  } else {
    -@math.floor(-scaled + 0.5)
  }
  (rounded / factor).to_string()
}

/// Left-pad a string to a fixed display width.
fn pad(s : String, width : Int) -> String {
  if s.length() >= width {
    s
  } else {
    " ".repeat(width - s.length()) + s
  }
}

fn main raise {
  println("=== RankProd Module Demo ===\n")

  // log2 expression: 12 genes (rows) x 6 samples (columns).
  // Columns 1-3 are class 0, columns 4-6 are class 1.
  let data = [
    [3.8, 4.2, 4.0, 6.1, 5.8, 6.3], // planted down (class1 < class2)
    [4.1, 3.7, 4.3, 5.9, 6.2, 5.7], // planted down
    [6.2, 5.9, 6.0, 4.1, 3.8, 4.2], // planted up
    [5.8, 6.3, 6.1, 4.3, 3.7, 4.0], // planted up
    [5.1, 4.8, 5.3, 4.9, 5.2, 4.7],
    [4.7, 5.2, 4.9, 5.3, 4.8, 5.1],
    [5.0, 4.6, 5.4, 4.8, 5.3, 4.9],
    [4.9, 5.3, 4.8, 5.1, 4.6, 5.4],
    [5.2, 4.9, 5.0, 4.7, 5.4, 4.8],
    [4.8, 5.1, 4.6, 5.4, 4.9, 5.2],
    [5.3, 4.7, 5.1, 5.0, 5.2, 4.6],
    [4.6, 5.4, 4.9, 5.2, 4.7, 5.3],
  ]
  let cl = [0, 0, 0, 1, 1, 1]
  println("Dataset: 12 genes x 6 samples, classes of 3 vs 3, log2 values.")

  // ---------------------------------------------------------------------------
  // 1. Rank product analysis
  // ---------------------------------------------------------------------------
  println("\n1. Rank product (geometric mean of per-pair ranks, all 9 pairings):")
  let rp = @src.rank_products(data=data, cl=cl, seed=2024)
  println(
    "  nrep = \{rp.nrep()}, n_features = \{rp.n_features()}, " +
      "pairings (k1 * k2 observed) = \{rp.random_pairs()}",
  )
  println("  gene    RPdown    RPup   AveFC(log2)   p.down    p.up   pfp.down  pfp.up")
  for g in 0..<rp.n_features() {
    println(
      "  \{pad(g.to_string(), 3)}  \{pad(sig4(rp.statistic_down()[g]), 7)} " +
        "\{pad(sig4(rp.statistic_up()[g]), 7)} " +
        "\{pad(sig4(rp.average_fold_change()[g]), 9)} " +
        "\{pad(sig4(rp.p_value_down()[g]), 9)} " +
        "\{pad(sig4(rp.p_value_up()[g]), 9)} " +
        "\{pad(sig4(rp.pfp_down()[g]), 9)} " +
        "\{pad(sig4(rp.pfp_up()[g]), 9)}",
    )
  }

  // ---------------------------------------------------------------------------
  // 2. Rank sum analysis
  // ---------------------------------------------------------------------------
  println("\n2. Rank sum (arithmetic mean of per-pair ranks, dice-sum p-values):")
  let rs = @src.rank_products(
    data=data,
    cl=cl,
    calculate_product=false,
    seed=2024,
  )
  for g in 0..<4 {
    println(
      "  gene \{g}: RSdown = \{sig4(rs.statistic_down()[g])}, " +
        "RSup = \{sig4(rs.statistic_up()[g])}, " +
        "p.down = \{sig4(rs.p_value_down()[g])}, p.up = \{sig4(rs.p_value_up()[g])}",
    )
  }

  // ---------------------------------------------------------------------------
  // 3. topGene-style selections
  // ---------------------------------------------------------------------------
  println("\n3. Top genes (topGene-style tables):")
  // Fixed count: two genes per direction.
  let (top_down_n, top_up_n) =
    @src.rankprod_top_genes(result=rp, num_gene=Some(2))
  let print_table =
    fn(title : String, table : @src.RankProdTopTable) {
      println("  \{title}")
      println("    gene  RP/RS    FC(class1/class2)    pfp    P.value")
      for i in 0..<table.indices().length() {
        println(
          "    \{pad(table.indices()[i].to_string(), 3)}  " +
            "\{pad(sig4(table.statistic()[i]), 7)} " +
            "\{pad(sig4(table.fold_change()[i]), 9)} " +
            "\{pad(sig4(table.pfp()[i]), 9)} " +
            "\{pad(sig4(table.p_value()[i]), 9)}",
        )
      }
    }
  print_table("Table 1: class1 < class2 (top 2)", top_down_n)
  print_table("Table 2: class1 > class2 (top 2)", top_up_n)

  // Threshold on pfp: the planted genes should be the only calls at 0.05.
  let (cut_down, cut_up) =
    @src.rankprod_top_genes(result=rp, cutoff=Some(0.05), method="pfp")
  println("\n  Genes called at pfp < 0.05:")
  // Tables are ordered by ascending statistic, so gene 3 precedes gene 2.
  println("    down: \{cut_down.indices()} (expected [0, 1])")
  println("    up:   \{cut_up.indices()} (expected [3, 2])")

  // Same selection with a raw p-value threshold.
  let (p_down, p_up) =
    @src.rankprod_top_genes(result=rp, cutoff=Some(0.01), method="pval")
  println("  Genes called at P.value < 0.01:")
  println("    down: \{p_down.indices()}, up: \{p_up.indices()}")

  // ---------------------------------------------------------------------------
  // 4. Exact rank-product / rank-sum distributions
  // ---------------------------------------------------------------------------
  println("\n4. Exact analytical p-values (no permutation needed):")
  // n = 12 genes, k = 3 pairings; smallest possible RP is 1.0 (rank 1 each).
  for rpval in [1.0, 1.26, 2.0, 3.0] {
    let p = @src.rank_product_pvalue(rpval, 12, 3)
    println("  RP = \{sig4(rpval)} with n = 12, k = 3 -> p = \{p}")
  }
  let (upper, lower, geometric) =
    @src.rank_product_pvalue_bounds(1.26, 12, 3)
  println(
    "  bounds at RP = 1.26: upper = \{upper}, lower = \{lower}, " +
      "geometric-mean interpolation = \{geometric}",
  )
  // Rank sum: sum s of k = 3 dice with n = 12 faces each (sums 3..36).
  // Argument order is (sum, k, n), mirroring the package dice.sum.cdf call.
  for s in [3.0, 6.0, 10.0, 18.0] {
    let p = @src.rank_sum_pvalue(s, 3, 12)
    let pmf = @src.rank_sum_pmf(s, 3, 12)
    println("  RS sum = \{s.to_int()} with n = 12, k = 3 -> P(S <= s) = \{p}, pmf = \{pmf}")
  }

  println("\n=== RankProd Module Demo Complete ===")
}
