///|
fn main {
  let data = @src.droplet_utils_advanced_example() catch {
    _ => abort("failed to create the DropletUtils example")
  }
  println(
    "Droplet matrix: " +
    data.n_features().to_string() +
    " features x " +
    data.n_barcodes().to_string() +
    " barcodes",
  )

  let ambient = @src.droplet_utils_ambient_profile(data, lower=15) catch {
    _ => abort("failed to estimate the ambient profile")
  }
  let alpha = @src.droplet_utils_estimate_alpha(data, ambient) catch {
    _ => abort("failed to estimate Dirichlet-multinomial alpha")
  }
  println(
    "Ambient droplets: " +
    ambient.n_empty().to_string() +
    ", ambient molecules: " +
    ambient.total.to_string() +
    ", estimated alpha: " +
    alpha.to_string(),
  )

  let ranks = @src.droplet_utils_barcode_ranks(
    data,
    lower=15,
    exclude_from=0,
    window=0.2,
  ) catch {
    _ => abort("failed to trace the barcode-rank curve")
  }
  println(
    "Barcode-rank knee: " +
    ranks.knee.to_string() +
    ", inflection: " +
    ranks.inflection.to_string(),
  )

  let config = @src.DropletUtilsConfig::create(
    lower=15,
    iterations=499,
    fdr_threshold=0.05,
    retain=80.0,
    alpha~,
    rank_exclude=0,
    rank_window=0.2,
    seed=31,
  ) catch {
    _ => abort("failed to create the DropletUtils configuration")
  }
  let result = @src.droplet_utils_empty_drops(data, config~) catch {
    _ => abort("emptyDrops failed")
  }
  println(result.summary())

  println("\nCalled barcodes")
  for barcode in 0..<result.barcodes.length() {
    if result.calls[barcode] {
      let p_value = match result.p_values[barcode] {
        Some(value) => value.to_string()
        None => "NA"
      }
      let fdr = match result.fdr[barcode] {
        Some(value) => value.to_string()
        None => "NA"
      }
      println(
        result.barcodes[barcode] +
        "\ttotal=" +
        result.totals[barcode].to_string() +
        "\tp=" +
        p_value +
        "\tFDR=" +
        fdr +
        "\tretained=" +
        result.always_retained[barcode].to_string(),
      )
    }
  }

  let filtered = result.filter(data) catch {
    _ => abort("failed to filter called barcodes")
  }
  println(
    "\nFiltered matrix: " +
    filtered.n_features().to_string() +
    " features x " +
    filtered.n_barcodes().to_string() +
    " called barcodes",
  )

  let assay : Array[Array[Double]] = []
  for feature in 0..<data.n_features() {
    let row : Array[Double] = []
    for barcode in 0..<data.n_barcodes() {
      row.push(data.counts[feature][barcode].to_double())
    }
    assay.push(row)
  }
  let experiment = @src.SingleCellExperiment::new(
    assay,
    data.feature_names,
    data.barcodes,
  )
  let output = @src.droplet_utils_empty_drops_sce(
    experiment,
    output_prefix="droplet",
    config~,
  ) catch {
    _ => abort("DropletUtils SingleCellExperiment integration failed")
  }
  println(
    "SingleCellExperiment write-back: class=" +
    output.experiment.col_data["droplet.class"][40] +
    ", ambient=" +
    output.experiment.row_data["droplet.ambient"][0] +
    ", original_unchanged=" +
    (!experiment.col_data.contains("droplet.class")).to_string(),
  )
}
