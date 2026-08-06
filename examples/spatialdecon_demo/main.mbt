// Bioconductor SpatialDecon-inspired background-aware spatial deconvolution.

///|
fn main {
  let (data, profile, nuclei_counts) = @src.spatial_decon_example_data() catch {
    SpatialDeconError(message) => abort("example data failed: " + message)
  }
  let config = @src.SpatialDeconConfig::create(
    rescale_profile=false,
    tolerance=1.0e-9,
  ) catch {
    SpatialDeconError(message) => abort("invalid configuration: " + message)
  }
  let result = @src.spatial_decon(data, profile, nuclei_counts~, config~) catch {
    SpatialDeconError(message) => abort("deconvolution failed: " + message)
  }

  println("=== Background-aware spatial deconvolution ===")
  println(result.summary())
  println("spot\tcell type\tabundance\tproportion\tcells")
  for spot in 0..<result.n_spots() {
    let ranked = result.top_cell_types(
      result.spot_ids[spot],
      limit=result.n_cell_types(),
    ) catch {
      SpatialDeconError(message) => abort("ranking failed: " + message)
    }
    for entry in ranked {
      let cell_type = entry.cell_type
      let mut cell_index = 0
      while result.cell_types[cell_index] != cell_type {
        cell_index = cell_index + 1
      }
      println(
        entry.spot_id +
        "\t" +
        cell_type +
        "\t" +
        entry.abundance.to_string() +
        "\t" +
        entry.proportion.to_string() +
        "\t" +
        result.cell_counts[cell_index][spot].to_string(),
      )
    }
  }
  println("fit diagnostics:")
  for spot in 0..<result.n_spots() {
    println(
      result.spot_ids[spot] +
      "\tRMSE=" +
      result.rmse[spot].to_string() +
      "\tcorrelation=" +
      result.correlations[spot].to_string(),
    )
  }

  let lymphoid = @src.SpatialDeconCellMerge::create("Lymphoid", [
    "T_cell", "B_cell",
  ]) catch {
    SpatialDeconError(message) => abort("invalid cell merge: " + message)
  }
  let myeloid = @src.SpatialDeconCellMerge::create("Myeloid", ["Myeloid"]) catch {
    SpatialDeconError(message) => abort("invalid cell merge: " + message)
  }
  let collapsed = @src.collapse_spatial_decon(result, [lymphoid, myeloid]) catch {
    SpatialDeconError(message) => abort("cell-type collapse failed: " + message)
  }
  println("\n=== Collapsed cell types ===")
  for group in 0..<collapsed.n_cell_types() {
    let mut line = collapsed.cell_types[group]
    for abundance in collapsed.beta[group] {
      line = line + "\t" + abundance.to_string()
    }
    println(line)
  }

  let reverse = @src.reverse_spatial_decon(data, result) catch {
    SpatialDeconError(message) =>
      abort("reverse deconvolution failed: " + message)
  }
  println("\n=== Reverse deconvolution ===")
  println("first-gene intercept=" + reverse.coefficients[0][0].to_string())
  for cell_type in 0..<reverse.cell_types.length() {
    println(
      reverse.cell_types[cell_type] +
      "\tcoefficient=" +
      reverse.coefficients[0][cell_type + 1].to_string(),
    )
  }
  println("correlation=" + reverse.correlations[0].to_string())

  let background = @src.spatial_decon_background(
    ["Neg_A", "Marker_A", "Neg_B", "Marker_B"],
    ["spot_1", "spot_2"],
    [[1.0, 2.0], [10.0, 12.0], [3.0, 4.0], [20.0, 24.0]],
    ["pool_A", "pool_A", "pool_B", "pool_B"],
    ["Neg_A", "Neg_B"],
  ) catch {
    SpatialDeconError(message) =>
      abort("background estimation failed: " + message)
  }
  println("\n=== Probe-pool background ===")
  for row in background {
    println(row[0].to_string() + "\t" + row[1].to_string())
  }

  let learned_profile = @src.create_spatial_decon_profile(
    ["g1", "g2", "g3"],
    ["a1", "a2", "b1", "b2"],
    ["A", "A", "B", "B"],
    [[10.0, 8.0, 1.0, 2.0], [1.0, 2.0, 9.0, 11.0], [2.0, 2.0, 2.0, 2.0]],
    normalize=true,
    scaling_factor=5.0,
    min_cells=1,
    min_genes=0,
  ) catch {
    SpatialDeconError(message) =>
      abort("profile construction failed: " + message)
  }
  println("\n=== Single-cell-derived profile ===")
  println(
    "genes=" +
    learned_profile.n_genes().to_string() +
    ", cell types=" +
    learned_profile.n_cell_types().to_string(),
  )
  for cell_type in learned_profile.cell_types {
    println("profile cell type: " + cell_type)
  }

  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "counts", data.values))
  ignore(@src.se_add_assay(experiment, "background", data.background))
  for gene in data.gene_names {
    ignore(@src.se_add_row(experiment, Map([("gene_id", gene)])))
  }
  for spot in 0..<data.spot_ids.length() {
    ignore(
      @src.se_add_col(
        experiment,
        Map([("spot_id", data.spot_ids[spot]), ("sample", "synthetic")]),
      ),
    )
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(spot.to_double(), 0.0),
      ),
    )
  }
  let integrated = @src.spatial_decon_spatial_experiment(
    experiment,
    profile,
    background_assay_name="background",
    nuclei_counts~,
    config~,
  ) catch {
    SpatialDeconError(message) =>
      abort("SpatialExperiment integration failed: " + message)
  }
  println("\n=== SpatialExperiment integration ===")
  println(
    "genes=" +
    integrated.experiment.metadata["spatialdecon_genes"] +
    ", spots=" +
    integrated.experiment.metadata["spatialdecon_spots"] +
    ", original unchanged=" +
    (!experiment.col_data[0].contains("SpatialDecon:T_cell")).to_string(),
  )
}
