///|
/// SummarizedExperiment - Multi-dimensional genomic data container inspired by Bioconductor SummarizedExperiment.

///|
pub struct SummarizedExperiment {
  assays : Map[String, Array[Array[Double]]]
  row_ranges : Array[(String, Int, Int)]
  col_data : Array[Map[String, String]]
  metadata : Map[String, String]
}

///|
pub fn summarized_experiment(
  assays : Map[String, Array[Array[Double]]],
  row_ranges : Array[(String, Int, Int)],
  col_data : Array[Map[String, String]],
  metadata : Map[String, String],
) -> SummarizedExperiment {
  SummarizedExperiment::{ assays, row_ranges, col_data, metadata }
}

///|
pub fn se_assay(
  se : SummarizedExperiment,
  name : String,
) -> Array[Array[Double]]? {
  se.assays.get(name)
}

///|
pub fn se_assay_names(se : SummarizedExperiment) -> Array[String] {
  se.assays.keys().collect()
}

///|
pub fn se_row_ranges(se : SummarizedExperiment) -> Array[(String, Int, Int)] {
  se.row_ranges
}

///|
pub fn se_col_data(se : SummarizedExperiment) -> Array[Map[String, String]] {
  se.col_data
}

///|
pub fn se_metadata(se : SummarizedExperiment) -> Map[String, String] {
  se.metadata
}

///|
pub fn se_nrow(se : SummarizedExperiment) -> Int {
  if se.row_ranges.length() > 0 {
    se.row_ranges.length()
  } else if se.assays.keys().collect().length() > 0 {
    let first_assay_name = se.assays.keys().collect()[0]
    se.assays[first_assay_name].length()
  } else {
    0
  }
}

///|
pub fn se_ncol(se : SummarizedExperiment) -> Int {
  if se.col_data.length() > 0 {
    se.col_data.length()
  } else if se.assays.keys().collect().length() > 0 {
    let first_assay_name = se.assays.keys().collect()[0]
    let assay = se.assays[first_assay_name]
    if assay.length() > 0 {
      assay[0].length()
    } else {
      0
    }
  } else {
    0
  }
}

///|
pub fn se_subset_rows(
  se : SummarizedExperiment,
  indices : Array[Int],
) -> SummarizedExperiment {
  let new_assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=se.assays.keys().collect().length(),
  )
  let assay_names = se.assays.keys().collect()

  for name_idx = 0; name_idx < assay_names.length(); name_idx = name_idx + 1 {
    let name = assay_names[name_idx]
    let old_assay = se.assays[name]
    let new_assay : Array[Array[Double]] = Array::make(
      indices.length(),
      Array::make(se_ncol(se), 0.0),
    )

    for i = 0; i < indices.length(); i = i + 1 {
      let old_idx = indices[i]
      for j = 0; j < se_ncol(se); j = j + 1 {
        new_assay[i][j] = old_assay[old_idx][j]
      }
    }

    new_assays[name] = new_assay
  }

  let new_row_ranges : Array[(String, Int, Int)] = Array::make(
    indices.length(),
    ("", 0, 0),
  )
  for i = 0; i < indices.length(); i = i + 1 {
    new_row_ranges[i] = se.row_ranges[indices[i]]
  }

  SummarizedExperiment::{
    assays: new_assays,
    row_ranges: new_row_ranges,
    col_data: se.col_data,
    metadata: se.metadata,
  }
}

///|
pub fn se_subset_cols(
  se : SummarizedExperiment,
  indices : Array[Int],
) -> SummarizedExperiment {
  let new_assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=se.assays.keys().collect().length(),
  )
  let assay_names = se.assays.keys().collect()

  for name_idx = 0; name_idx < assay_names.length(); name_idx = name_idx + 1 {
    let name = assay_names[name_idx]
    let old_assay = se.assays[name]
    let new_assay : Array[Array[Double]] = Array::make(
      se_nrow(se),
      Array::make(indices.length(), 0.0),
    )

    for i = 0; i < se_nrow(se); i = i + 1 {
      for j = 0; j < indices.length(); j = j + 1 {
        new_assay[i][j] = old_assay[i][indices[j]]
      }
    }

    new_assays[name] = new_assay
  }

  let new_col_data : Array[Map[String, String]] = Array::make(
    indices.length(),
    Map([], capacity=0),
  )
  for i = 0; i < indices.length(); i = i + 1 {
    new_col_data[i] = se.col_data[indices[i]]
  }

  SummarizedExperiment::{
    assays: new_assays,
    row_ranges: se.row_ranges,
    col_data: new_col_data,
    metadata: se.metadata,
  }
}

///|
pub fn se_row_means(
  se : SummarizedExperiment,
  assay_name : String,
) -> Array[Double] {
  match se.assays.get(assay_name) {
    Some(assay) => {
      let n_rows = assay.length()
      let n_cols = if assay.length() > 0 { assay[0].length() } else { 0 }
      let means : Array[Double] = Array::make(n_rows, 0.0)

      for i = 0; i < n_rows; i = i + 1 {
        let mut sum = 0.0
        for j = 0; j < n_cols; j = j + 1 {
          sum = sum + assay[i][j]
        }
        means[i] = sum / n_cols.to_double()
      }

      means
    }
    None => []
  }
}

///|
pub fn se_col_means(
  se : SummarizedExperiment,
  assay_name : String,
) -> Array[Double] {
  match se.assays.get(assay_name) {
    Some(assay) => {
      let n_rows = assay.length()
      let n_cols = if assay.length() > 0 { assay[0].length() } else { 0 }
      let means : Array[Double] = Array::make(n_cols, 0.0)

      for j = 0; j < n_cols; j = j + 1 {
        let mut sum = 0.0
        for i = 0; i < n_rows; i = i + 1 {
          sum = sum + assay[i][j]
        }
        means[j] = sum / n_rows.to_double()
      }

      means
    }
    None => []
  }
}

///|
pub fn se_merge(
  se1 : SummarizedExperiment,
  se2 : SummarizedExperiment,
) -> SummarizedExperiment {
  let new_assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=se1.assays.keys().collect().length(),
  )
  let assay_names = se1.assays.keys().collect()

  for name_idx = 0; name_idx < assay_names.length(); name_idx = name_idx + 1 {
    let name = assay_names[name_idx]
    let assay1 = se1.assays[name]
    let assay2 = match se2.assays.get(name) {
      Some(a) => a
      None => assay1
    }

    let n_rows1 = assay1.length()
    let n_rows2 = assay2.length()
    let n_cols = if assay1.length() > 0 { assay1[0].length() } else { 0 }

    let new_assay : Array[Array[Double]] = Array::make(
      n_rows1 + n_rows2,
      Array::make(n_cols, 0.0),
    )

    for i = 0; i < n_rows1; i = i + 1 {
      for j = 0; j < n_cols; j = j + 1 {
        new_assay[i][j] = assay1[i][j]
      }
    }

    for i = 0; i < n_rows2; i = i + 1 {
      for j = 0; j < n_cols; j = j + 1 {
        new_assay[n_rows1 + i][j] = assay2[i][j]
      }
    }

    new_assays[name] = new_assay
  }

  let new_row_ranges : Array[(String, Int, Int)] = Array::make(
    se1.row_ranges.length() + se2.row_ranges.length(),
    ("", 0, 0),
  )
  for i = 0; i < se1.row_ranges.length(); i = i + 1 {
    new_row_ranges[i] = se1.row_ranges[i]
  }
  for i = 0; i < se2.row_ranges.length(); i = i + 1 {
    new_row_ranges[se1.row_ranges.length() + i] = se2.row_ranges[i]
  }

  let new_col_data : Array[Map[String, String]] = Array::make(
    se1.col_data.length() + se2.col_data.length(),
    Map([], capacity=0),
  )
  for i = 0; i < se1.col_data.length(); i = i + 1 {
    new_col_data[i] = se1.col_data[i]
  }
  for i = 0; i < se2.col_data.length(); i = i + 1 {
    new_col_data[se1.col_data.length() + i] = se2.col_data[i]
  }

  let new_metadata : Map[String, String] = Map(
    [],
    capacity=se1.metadata.keys().collect().length(),
  )
  let meta_keys = se1.metadata.keys().collect()
  for k = 0; k < meta_keys.length(); k = k + 1 {
    new_metadata[meta_keys[k]] = se1.metadata[meta_keys[k]]
  }

  SummarizedExperiment::{
    assays: new_assays,
    row_ranges: new_row_ranges,
    col_data: new_col_data,
    metadata: new_metadata,
  }
}
