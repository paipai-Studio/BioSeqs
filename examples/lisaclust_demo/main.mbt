// Bioconductor lisaClust-inspired local spatial association and region discovery.

///|
fn main {
  let cells = @src.lisaclust_example_data()
  let config = @src.LisaConfig::create(
    radii=[2.5, 5.0],
    curve_kind=@src.lisa_k_curve(),
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=4.0,
    min_density=0.05,
    edge_correction=true,
    window_padding=0.1,
    edge_samples=96,
    n_clusters=2,
    n_starts=8,
    seed=11,
  ) catch {
    LisaError(message) => abort("invalid lisaClust configuration: " + message)
  }
  let result = @src.lisaclust(cells, config~) catch {
    LisaError(message) => abort("lisaClust analysis failed: " + message)
  }

  println("=== Direct lisaClust analysis ===")
  println(result.summary())
  println(
    "iterations=" +
    result.iterations.to_string() +
    ", converged=" +
    result.converged.to_string(),
  )
  for image_radii in result.curves.effective_radii {
    println(
      image_radii.image_id +
      " effective radii: " +
      image_radii.radii.to_string(),
    )
  }

  println("\n=== Region summaries ===")
  println("region\tsize\tdominant cell type\tmaximum enrichment")
  for summary in result.region_summaries {
    println(
      summary.region +
      "\t" +
      summary.size.to_string() +
      "\t" +
      summary.dominant_cell_type +
      "\t" +
      summary.maximum_enrichment.to_string(),
    )
  }

  let top = result.top_enrichments(limit=6, minimum_relative_frequency=1.0) catch {
    LisaError(message) => abort("enrichment ranking failed: " + message)
  }
  println("\n=== Top observed/expected enrichments ===")
  println("region\tcell type\tobserved\texpected\tratio")
  for entry in top {
    println(
      entry.region +
      "\t" +
      entry.cell_type +
      "\t" +
      entry.observed.to_string() +
      "\t" +
      entry.expected.to_string() +
      "\t" +
      entry.relative_frequency.to_string(),
    )
  }

  let l_config = @src.LisaConfig::create(
    radii=[2.5, 5.0],
    curve_kind=@src.lisa_l_curve(),
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=4.0,
    min_density=0.05,
    edge_correction=true,
    window_padding=0.1,
    edge_samples=96,
  ) catch {
    LisaError(message) => abort("invalid local-L configuration: " + message)
  }
  let l_curves = @src.lisa_curves(cells, config=l_config) catch {
    LisaError(message) => abort("local-L curve computation failed: " + message)
  }
  println("\n=== Centered local-L curve for the first cell ===")
  println("cell=" + l_curves.cell_ids[0])
  for feature in 0..<l_curves.n_features() {
    println(
      l_curves.feature_names[feature] +
      "\t" +
      l_curves.values[0][feature].to_string(),
    )
  }

  let experiment = @src.SpatialExperiment::new()
  for cell in cells {
    ignore(
      @src.se_add_col(
        experiment,
        Map([
          ("cellID", cell.cell_id),
          ("imageID", cell.image_id),
          ("cellType", cell.cell_type),
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
  let integrated = @src.lisaclust_spatial_experiment(
    experiment,
    cell_id_key="cellID",
    region_key="tissueRegion",
    config~,
  ) catch {
    LisaError(message) =>
      abort("lisaClust SpatialExperiment integration failed: " + message)
  }
  println("\n=== SpatialExperiment integration ===")
  println(
    "cells=" +
    integrated.experiment.metadata["lisaclust_cells"] +
    ", features=" +
    integrated.experiment.metadata["lisaclust_features"] +
    ", regions=" +
    integrated.experiment.metadata["lisaclust_regions"],
  )
  println(
    "first assigned tissueRegion=" +
    integrated.experiment.col_data[0]["tissueRegion"] +
    ", original unchanged=" +
    (!experiment.col_data[0].contains("tissueRegion")).to_string(),
  )
}
