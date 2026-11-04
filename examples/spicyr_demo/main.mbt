// Bioconductor spicyR-inspired differential spatial colocalization.

///|
fn main {
  let (cells, metadata, pairs) = @src.spicyr_example_data()
  let config = @src.SpicyConfig::create(
    radii=[0.5, 1.0, 2.0],
    edge_correction=true,
    use_weights=true,
    reference_condition="control",
    fdr_threshold=0.1,
    tolerance=0.01,
  ) catch {
    SpicyError(message) => abort("invalid spicyR configuration: " + message)
  }

  let result = @src.spicyr(
    cells,
    metadata,
    pairs~,
    covariate_names=["age"],
    config~,
  ) catch {
    SpicyError(message) => abort("spicyR analysis failed: " + message)
  }

  println("=== Direct spicyR analysis ===")
  println(result.summary())
  println("pair\tcondition\testimate\tSE\tp-value\tFDR\tmodel")
  for pair_result in result.pairs {
    for condition in pair_result.conditions {
      let model = if pair_result.model
        is @src.SpicyModelKind::WeightedRandomIntercept {
        "random-intercept"
      } else {
        "linear"
      }
      println(
        pair_result.pair.label() +
        "\t" +
        condition.condition_level +
        "\t" +
        condition.estimate.to_string() +
        "\t" +
        condition.standard_error.to_string() +
        "\t" +
        condition.p_value.to_string() +
        "\t" +
        condition.adjusted_p_value.to_string() +
        "\t" +
        model,
      )
    }
  }

  let first_pair = result.pairs[0]
  let first_image = first_pair.associations[0]
  println("\n=== Per-image cross-L curve ===")
  println(
    "image=" +
    first_image.image_id +
    ", pair=" +
    first_pair.pair.label() +
    ", summary=" +
    first_image.statistic.to_string() +
    ", precision_weight=" +
    first_image.weight.to_string(),
  )
  for index in 0..<first_image.radii.length() {
    println(
      "r=" +
      first_image.radii[index].to_string() +
      "\tL(r)=" +
      first_image.l_values[index].to_string(),
    )
  }

  let refitted = @src.spicy_fit_associations(
    first_pair.associations,
    metadata,
    covariate_names=["age"],
    config~,
  ) catch {
    SpicyError(message) =>
      abort("precomputed spicyR association fit failed: " + message)
  }
  println("\n=== Precomputed association input ===")
  println(refitted.summary())

  let experiment = @src.SpatialExperiment::new()
  for cell in cells {
    let mut image_index = -1
    for index in 0..<metadata.length() {
      if metadata[index].image_id == cell.image_id {
        image_index = index
        break
      }
    }
    if image_index < 0 {
      abort("example cell references unknown image")
    }
    let image = metadata[image_index]
    ignore(
      @src.se_add_col(
        experiment,
        Map([
          ("imageID", cell.image_id),
          ("cellType", cell.cell_type),
          ("condition", image.condition),
          ("subject", image.subject),
          ("age", image.covariates["age"].to_string()),
        ]),
      ),
    )
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(cell.x, cell.y),
      ),
    )
  }
  let integrated = @src.spicyr_spatial_experiment(
    experiment,
    subject_key="subject",
    covariate_keys=["age"],
    pairs~,
    config~,
  ) catch {
    SpicyError(message) =>
      abort("spicyR SpatialExperiment integration failed: " + message)
  }
  println("\n=== SpatialExperiment integration ===")
  println(
    "images=" +
    integrated.experiment.metadata["spicyr_images"] +
    ", pairs=" +
    integrated.experiment.metadata["spicyr_pairs"] +
    ", reference=" +
    integrated.experiment.metadata["spicyr_reference"],
  )
}
