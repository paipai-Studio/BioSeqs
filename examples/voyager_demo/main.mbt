// Bioconductor Voyager-inspired spatial autocorrelation workflow.
//
// Exercises the full univariate pipeline on a synthetic SpatialExperiment:
//   1. Spatial weight construction (kNN, distance-band, inverse-distance).
//   2. Global Moran's I and Geary's c with analytical inference.
//   3. Local Moran's I (LISA) with permutation inference and quadrants.
//   4. Local Getis–Ord Gi* hotspot detection.
//   5. Bivariate Lee's L spatial association.
//   6. Empirical variogram and spherical model fitting.
//   7. Moran correlogram over distance bins.
//   8. Immutable SpatialExperiment write-back of local statistics.

///|
fn voyager_demo_coords(se : @src.SpatialExperiment) -> Array[Array[Double]] {
  let coords : Array[Array[Double]] = []
  for sc in se.spatial_coords {
    coords.push([sc.x, sc.y])
  }
  coords
}

///|
fn voyager_demo_round(value : Double, digits : Int) -> Double {
  let scale = @math.pow(10.0, digits.to_double())
  (value * scale).round() / scale
}

///|
fn voyager_demo_print_row(label : String, value : Double) -> Unit {
  println("  " + label + ": " + voyager_demo_round(value, 4).to_string())
}

///|
fn main {
  let se = @src.voyager_example_spatial_experiment()
  let coords = voyager_demo_coords(se)
  // assay[feature][spot]; the example ships 4 features on a 6x6 grid.
  let assay = se.assay["logcounts"]
  let gradient = assay[0]
  let hotspot = assay[1]
  let anti = assay[3]
  println("=== Bioconductor Voyager Demo ===")
  println(
    "SpatialExperiment: " +
    se.col_data.length().to_string() +
    " spots, " +
    assay.length().to_string() +
    " features, platform=" +
    se.metadata["platform"],
  )

  // -----------------------------------------------------------------------
  // 1. Spatial weight construction
  // -----------------------------------------------------------------------
  println("\n1. Spatial weight construction")
  let knn = @src.voyager_weights_knn(coords, 4) catch {
    VoyagerError(message) => abort("kNN weights failed: " + message)
  }
  println(
    "  kNN k=4: style=" +
    knn.style +
    ", S0=" +
    voyager_demo_round(@src.voyager_weights_s0(knn), 4).to_string(),
  )
  let dband = @src.voyager_weights_distance_band(coords, 1.0) catch {
    VoyagerError(message) => abort("distance-band weights failed: " + message)
  }
  println(
    "  distance-band 1.0: centre spot neighbours=" +
    dband.neighbors[18].length().to_string(),
  )
  let idw = @src.voyager_weights_inverse_distance(coords, 1.0) catch {
    VoyagerError(message) =>
      abort("inverse-distance weights failed: " + message)
  }
  println(
    "  inverse-distance power=1: centre spot neighbours=" +
    idw.neighbors[18].length().to_string(),
  )

  // -----------------------------------------------------------------------
  // 2. Global Moran's I and Geary's c
  // -----------------------------------------------------------------------
  println("\n2. Global spatial autocorrelation (gene_gradient)")
  let moran = @src.voyager_global_morans_i(gradient, knn) catch {
    VoyagerError(message) => abort("global Moran's I failed: " + message)
  }
  voyager_demo_print_row("Moran's I estimate", moran.estimate)
  voyager_demo_print_row("Moran's I expectation", moran.expectation)
  voyager_demo_print_row("Moran's I z-score", moran.z_score)
  voyager_demo_print_row("Moran's I p-value", moran.p_value)
  let geary = @src.voyager_global_gearys_c(gradient, knn) catch {
    VoyagerError(message) => abort("global Geary's c failed: " + message)
  }
  voyager_demo_print_row("Geary's c estimate", geary.estimate)
  voyager_demo_print_row("Geary's c p-value", geary.p_value)

  // -----------------------------------------------------------------------
  // 3. Local Moran's I (LISA) with permutation inference
  // -----------------------------------------------------------------------
  println("\n3. Local Moran's I (LISA) with permutation inference")
  let lisa = @src.voyager_local_morans_i(
    gradient,
    knn,
    permutations=199,
    seed=20240501,
    fdr_threshold=0.1,
  ) catch {
    VoyagerError(message) => abort("local Moran's I failed: " + message)
  }
  let mut hh = 0
  let mut ll = 0
  let mut sig = 0
  for i in 0..<lisa.n {
    let q = lisa.quadrants[i]
    if q == "HH" {
      hh = hh + 1
    }
    if q == "LL" {
      ll = ll + 1
    }
    if q != "not significant" {
      sig = sig + 1
    }
  }
  println(
    "  significant spots=" +
    sig.to_string() +
    ", HH=" +
    hh.to_string() +
    ", LL=" +
    ll.to_string(),
  )

  // -----------------------------------------------------------------------
  // 4. Local Getis–Ord Gi* hotspot detection
  // -----------------------------------------------------------------------
  println("\n4. Local Getis–Ord Gi* hotspot detection (gene_hotspot)")
  let getis = @src.voyager_local_getis_ord(
    hotspot,
    knn,
    star=true,
    fdr_threshold=0.1,
  ) catch {
    VoyagerError(message) => abort("Getis–Ord failed: " + message)
  }
  let mut max_z = -1.0e300
  let mut max_spot = 0
  for i in 0..<getis.n {
    if getis.z_scores[i] > max_z {
      max_z = getis.z_scores[i]
      max_spot = i
    }
  }
  println(
    "  peak Gi* z-score=" +
    voyager_demo_round(max_z, 4).to_string() +
    " at spot " +
    max_spot.to_string() +
    " (" +
    getis.classifications[max_spot] +
    ")",
  )

  // -----------------------------------------------------------------------
  // 5. Bivariate Lee's L
  // -----------------------------------------------------------------------
  println("\n5. Bivariate Lee's L (gene_gradient vs gene_anti)")
  let lees = @src.voyager_lees_l(
    gradient,
    anti,
    knn,
    permutations=199,
    seed=20240501,
    fdr_threshold=0.1,
  ) catch {
    VoyagerError(message) => abort("Lee's L failed: " + message)
  }
  voyager_demo_print_row("global Lee's L", lees.global_l)

  // -----------------------------------------------------------------------
  // 6. Empirical variogram and spherical model fit
  // -----------------------------------------------------------------------
  println("\n6. Empirical variogram and spherical model fit (gene_hotspot)")
  let empirical = @src.voyager_empirical_variogram(coords, hotspot, n_lags=8) catch {
    VoyagerError(message) => abort("empirical variogram failed: " + message)
  }
  println("  empirical points=" + empirical.length().to_string())
  let model = @src.voyager_fit_variogram(empirical, model_type="spherical") catch {
    VoyagerError(message) => abort("variogram fit failed: " + message)
  }
  println(
    "  spherical: nugget=" +
    voyager_demo_round(model.nugget, 4).to_string() +
    ", sill=" +
    voyager_demo_round(model.sill, 4).to_string() +
    ", range=" +
    voyager_demo_round(model.range, 4).to_string() +
    ", SSE=" +
    voyager_demo_round(model.fitted_sse, 4).to_string(),
  )
  voyager_demo_print_row(
    "variogram predicted semivariance at range",
    @src.voyager_variogram_predict(model, model.range),
  )

  // -----------------------------------------------------------------------
  // 7. Moran correlogram
  // -----------------------------------------------------------------------
  println("\n7. Moran correlogram over distance bins")
  let corr = @src.voyager_correlogram(coords, gradient, n_lags=5) catch {
    VoyagerError(message) => abort("correlogram failed: " + message)
  }
  println("  bins=" + corr.length().to_string())
  for point in corr {
    println(
      "  lag=" +
      voyager_demo_round(point.lag, 3).to_string() +
      ", Moran's I=" +
      voyager_demo_round(point.morans_i, 4).to_string() +
      ", npairs=" +
      point.npairs.to_string(),
    )
  }

  // -----------------------------------------------------------------------
  // 8. Immutable SpatialExperiment write-back
  // -----------------------------------------------------------------------
  println("\n8. Immutable SpatialExperiment write-back")
  let output = @src.voyager_run_univariate_sfe(
    se,
    [0, 1],
    assay_name="logcounts",
    stat_method="moran",
    permutations=99,
    seed=2024,
    fdr_threshold=0.1,
    output_prefix="voyager",
  ) catch {
    VoyagerError(message) =>
      abort("Voyager SpatialExperiment integration failed: " + message)
  }
  println(
    "  features analyzed=" +
    output.results.length().to_string() +
    ", method=" +
    output.results[0].stat_method,
  )
  println(
    "  input unchanged: " +
    (!se.col_data[0].contains("voyager.moran.local.gene_gradient")).to_string(),
  )
  println(
    "  output colData carries local Moran: " +
    output.experiment.col_data[0]
    .contains("voyager.moran.local.gene_gradient")
    .to_string(),
  )
  println(
    "  metadata method=" +
    output.experiment.metadata["voyager.method"] +
    ", n_features=" +
    output.experiment.metadata["voyager.n_features"],
  )
  println("=== Demo Complete ===")
}
