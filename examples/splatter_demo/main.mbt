///|
/// Bioconductor splatter Demo - single-cell RNA-seq count simulation.
///
/// Demonstrates the port of the Bioconductor `splatter` package
/// (Zappia et al., 2017): the splat model with gamma-distributed gene
/// base means, log-normal library size factors, trended BCV dispersion,
/// negative-binomial counts, group-wise differential expression, batch
/// effects, logistic dropout and pseudotime paths. All randomness comes
/// from a deterministic LCG stream so every run is reproducible.

///|
fn mean(values : Array[Double]) -> Double {
  let mut total = 0.0
  for v in values {
    total = total + v
  }
  total / values.length().to_double()
}

///|
fn sum_counts(counts : Array[Array[Int]]) -> Int {
  let mut total = 0
  for row in counts {
    for c in row {
      total = total + c
    }
  }
  total
}

///|
fn nonzero_fraction(counts : Array[Array[Int]]) -> Double {
  let mut nonzero = 0
  let mut total = 0
  for row in counts {
    for c in row {
      if c > 0 {
        nonzero = nonzero + 1
      }
      total = total + 1
    }
  }
  nonzero.to_double() / total.to_double()
}

///|
fn main {
  println("=== splat: basic simulation (defaults: 100 genes x 100 cells) ===")
  let params = @src.SplatParams::default()
  let result = @src.splat_simulate(params, 12345)
  println("counts shape   : \{result.counts.length()} genes x \{result.counts[0].length()} cells")
  println("total counts   : \{sum_counts(result.counts)}")
  println("nonzero rate   : \{nonzero_fraction(result.counts)}")
  println("lib factors    : \{result.cell_lib_factors[0]}, \{result.cell_lib_factors[1]}, ... (mean=\{mean(result.cell_lib_factors)})")
  println("groups         : \{result.cell_group[0]} (single group)")
  println("batches        : \{result.cell_batch[0]} (single batch)")
  println("dropout events : \{result.n_dropout_events} (dropout disabled)")

  println("\n=== splat: group DE (2 groups x 50 cells, 60% DE) ===")
  let de_params = @src.SplatParams::create(
    200,
    0,
    group_cells=[50, 50],
    de_prob=0.6,
    de_down_prob=0.5,
    de_fac_loc=0.4,
    de_fac_scale=0.5,
  )
  let de_result = @src.splat_simulate(de_params, 67890)
  println("cell groups    : \{de_result.cell_group[0]}..=\{de_result.cell_group[49]}, \{de_result.cell_group[50]}..=\{de_result.cell_group[99]}")
  println("group log-means: \{de_result.group_log_means.length()} groups x \{de_result.group_log_means[0].length()} genes")
  // a strongly DE gene shows a large log2 fold change between groups
  let mut max_lfc = 0.0
  let mut max_gene = 0
  for g in 0..<de_result.group_log_means[0].length() {
    let lfc = (de_result.group_log_means[1][g] -
      de_result.group_log_means[0][g]).abs()
    if lfc > max_lfc {
      max_lfc = lfc
      max_gene = g
    }
  }
  println("max |log2FC|   : gene \{max_gene}, lfc=\{max_lfc}")

  println("\n=== splat: batch effects (2 batches, 40 + 60 cells) ===")
  let batch_params = @src.SplatParams::create(
    100,
    100,
    batch_cells=[40, 60],
    batch_fac_loc=0.2,
    batch_fac_scale=0.3,
  )
  let batch_result = @src.splat_simulate(batch_params, 2468)
  println("cell batches   : \{batch_result.cell_batch[0]} x40, \{batch_result.cell_batch[40]} x60")
  println("total counts   : \{sum_counts(batch_result.counts)}")

  println("\n=== splat: dropout (mid=2, shape=-1) ===")
  let drop_params = @src.SplatParams::create(
    100,
    100,
    enable_dropout=true,
    dropout_mid=2.0,
    dropout_shape=-1.0,
  )
  let drop_result = @src.splat_simulate(drop_params, 13579)
  println("dropout events : \{drop_result.n_dropout_events} zeros of 10000 entries")
  println("nonzero rate   : \{nonzero_fraction(drop_result.counts)} (vs \{nonzero_fraction(result.counts)} without dropout)")

  println("\n=== splat: pseudotime path (30 cells, 50% path genes) ===")
  let path_params = @src.SplatParams::create(
    100,
    100,
    enable_path=true,
    n_path_steps=30,
    path_probability=0.5,
    path_shift_scale=1.0,
  )
  let path_result = @src.splat_simulate(path_params, 97531)
  println("pseudotime     : \{path_result.cell_pseudotime[0]}, ..., \{path_result.cell_pseudotime[14]}, ..., \{path_result.cell_pseudotime[29]}")
  println("groups along path: \{path_result.cell_group[0]}..=\{path_result.cell_group[29]}")

  println("\n=== splat: deterministic RNG streams ===")
  let (s1, u1) = @src.splat_lcg(42)
  let (s2, n1) = @src.splat_normal(s1)
  let (s3, g1) = @src.splat_gamma(0.6, s2)
  let (s4, p1) = @src.splat_poisson(3.5, s3)
  println("lcg(42)    -> u=\{u1}")
  println("normal     -> \{n1}")
  println("gamma(0.6) -> \{g1}")
  println("poisson(3.5) -> \{p1}")
  let _ = s4
  // same seed reproduces the same simulation
  let replay = @src.splat_simulate(params, 12345)
  let mut identical = true
  for g in 0..<result.counts.length() {
    for c in 0..<result.counts[g].length() {
      if result.counts[g][c] != replay.counts[g][c] {
        identical = false
      }
    }
  }
  println("replay with seed 12345 identical: \{identical}")
}
