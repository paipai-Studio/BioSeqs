///|
pub struct ExpressionSet {
  assay_data : Array[Array[Double]]
  pheno_data : Map[String, Array[String]]
  feature_data : Map[String, Array[String]]
  annotation : String
  protocol_data : Map[String, Array[String]]
}

///|
pub struct AnnotatedDataFrame {
  data : Map[String, Array[String]]
  var_metadata : Map[String, String]
}

///|
pub fn ExpressionSet::new(
  assay_data : Array[Array[Double]],
  pheno_data : Map[String, Array[String]],
  feature_data : Map[String, Array[String]],
  annotation : String,
  protocol_data : Map[String, Array[String]],
) -> ExpressionSet {
  ExpressionSet::{
    assay_data,
    pheno_data,
    feature_data,
    annotation,
    protocol_data,
  }
}

///|
pub fn expression_set_nrow(eset : ExpressionSet) -> Int {
  eset.assay_data.length()
}

///|
pub fn expression_set_ncol(eset : ExpressionSet) -> Int {
  if eset.assay_data.length() > 0 {
    eset.assay_data[0].length()
  } else {
    0
  }
}

///|
pub fn expression_set_sample_names(eset : ExpressionSet) -> Array[String] {
  if eset.pheno_data.contains("sample_name") {
    eset.pheno_data["sample_name"]
  } else if eset.pheno_data.contains("samples") {
    eset.pheno_data["samples"]
  } else {
    let n = expression_set_ncol(eset)
    let result : Array[String] = Array::new()
    let mut i = 0
    while i < n {
      result.push("sample_" + i.to_string())
      i = i + 1
    }
    result
  }
}

///|
pub fn expression_set_feature_names(eset : ExpressionSet) -> Array[String] {
  if eset.feature_data.contains("feature_name") {
    eset.feature_data["feature_name"]
  } else if eset.feature_data.contains("genes") {
    eset.feature_data["genes"]
  } else {
    let n = expression_set_nrow(eset)
    let result : Array[String] = Array::new()
    let mut i = 0
    while i < n {
      result.push("gene_" + i.to_string())
      i = i + 1
    }
    result
  }
}

///|
pub fn expression_set_get_assay(eset : ExpressionSet) -> Array[Array[Double]] {
  eset.assay_data
}

///|
pub fn expression_set_get_pheno_data(
  eset : ExpressionSet,
) -> Map[String, Array[String]] {
  eset.pheno_data
}

///|
pub fn expression_set_get_feature_data(
  eset : ExpressionSet,
) -> Map[String, Array[String]] {
  eset.feature_data
}

///|
pub fn expression_set_get_annotation(eset : ExpressionSet) -> String {
  eset.annotation
}

///|
pub fn expression_set_subset(
  eset : ExpressionSet,
  rows : Array[Int],
  cols : Array[Int],
) -> ExpressionSet {
  let new_assay : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < rows.length() {
    let row_idx = rows[i]
    let new_row : Array[Double] = Array::new()
    let mut j = 0
    while j < cols.length() {
      let col_idx = cols[j]
      if row_idx < eset.assay_data.length() &&
        col_idx < eset.assay_data[row_idx].length() {
        new_row.push(eset.assay_data[row_idx][col_idx])
      } else {
        new_row.push(0.0)
      }
      j = j + 1
    }
    new_assay.push(new_row)
    i = i + 1
  }

  let new_pheno : Map[String, Array[String]] = Map(
    [],
    capacity=eset.pheno_data.length(),
  )
  let pheno_keys = eset.pheno_data.keys().collect()
  i = 0
  while i < pheno_keys.length() {
    let key = pheno_keys[i]
    let values = eset.pheno_data[key]
    let new_values : Array[String] = Array::new()
    let mut j = 0
    while j < cols.length() {
      let col_idx = cols[j]
      if col_idx < values.length() {
        new_values.push(values[col_idx])
      } else {
        new_values.push("")
      }
      j = j + 1
    }
    new_pheno[key] = new_values
    i = i + 1
  }

  let new_feature : Map[String, Array[String]] = Map(
    [],
    capacity=eset.feature_data.length(),
  )
  let feature_keys = eset.feature_data.keys().collect()
  i = 0
  while i < feature_keys.length() {
    let key = feature_keys[i]
    let values = eset.feature_data[key]
    let new_values : Array[String] = Array::new()
    let mut j = 0
    while j < rows.length() {
      let row_idx = rows[j]
      if row_idx < values.length() {
        new_values.push(values[row_idx])
      } else {
        new_values.push("")
      }
      j = j + 1
    }
    new_feature[key] = new_values
    i = i + 1
  }

  let annotation = eset.annotation
  let protocol_data = eset.protocol_data
  ExpressionSet::new(
    new_assay, new_pheno, new_feature, annotation, protocol_data,
  )
}

///|
pub fn expression_set_filter_by_row(
  eset : ExpressionSet,
  threshold : Double,
) -> ExpressionSet {
  let keep_rows : Array[Int] = Array::new()
  let mut i = 0
  while i < eset.assay_data.length() {
    let mut max_val = 0.0
    let mut j = 0
    while j < eset.assay_data[i].length() {
      if eset.assay_data[i][j] > max_val {
        max_val = eset.assay_data[i][j]
      }
      j = j + 1
    }
    if max_val >= threshold {
      keep_rows.push(i)
    }
    i = i + 1
  }
  let ncols = expression_set_ncol(eset)
  let cols : Array[Int] = Array::new()
  i = 0
  while i < ncols {
    cols.push(i)
    i = i + 1
  }
  expression_set_subset(eset, keep_rows, cols)
}

///|
pub fn expression_set_normalize(eset : ExpressionSet) -> ExpressionSet {
  let nrows = expression_set_nrow(eset)
  let ncols = expression_set_ncol(eset)

  let col_means : Array[Double] = Array::new()
  let mut i = 0
  while i < ncols {
    col_means.push(0.0)
    i = i + 1
  }

  i = 0
  while i < ncols {
    let mut sum = 0.0
    let mut j = 0
    while j < nrows {
      sum = sum + eset.assay_data[j][i]
      j = j + 1
    }
    col_means[i] = sum / nrows.to_double()
    i = i + 1
  }

  let col_sds : Array[Double] = Array::new()
  i = 0
  while i < ncols {
    col_sds.push(0.0)
    i = i + 1
  }

  i = 0
  while i < ncols {
    let mut sum_sq = 0.0
    let mut j = 0
    while j < nrows {
      let diff = eset.assay_data[j][i] - col_means[i]
      sum_sq = sum_sq + diff * diff
      j = j + 1
    }
    col_sds[i] = (sum_sq / nrows.to_double()).sqrt()
    i = i + 1
  }

  let new_assay : Array[Array[Double]] = Array::new()
  i = 0
  while i < nrows {
    let new_row : Array[Double] = Array::new()
    let mut j = 0
    while j < ncols {
      if col_sds[j] > 0.0 {
        new_row.push((eset.assay_data[i][j] - col_means[j]) / col_sds[j])
      } else {
        new_row.push(0.0)
      }
      j = j + 1
    }
    new_assay.push(new_row)
    i = i + 1
  }

  let pheno_data = eset.pheno_data
  let feature_data = eset.feature_data
  let annotation = eset.annotation
  let protocol_data = eset.protocol_data
  ExpressionSet::new(
    new_assay, pheno_data, feature_data, annotation, protocol_data,
  )
}

///|
pub fn expression_set_log2(eset : ExpressionSet) -> ExpressionSet {
  let nrows = expression_set_nrow(eset)
  let ncols = expression_set_ncol(eset)

  let new_assay : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < nrows {
    let new_row : Array[Double] = Array::new()
    let mut j = 0
    while j < ncols {
      let val = eset.assay_data[i][j]
      if val > 0.0 {
        new_row.push(@math.log2(val))
      } else {
        new_row.push(0.0)
      }
      j = j + 1
    }
    new_assay.push(new_row)
    i = i + 1
  }

  let pheno_data = eset.pheno_data
  let feature_data = eset.feature_data
  let annotation = eset.annotation
  let protocol_data = eset.protocol_data
  ExpressionSet::new(
    new_assay, pheno_data, feature_data, annotation, protocol_data,
  )
}

///|
pub fn expression_set_summary(eset : ExpressionSet) -> String {
  let nrows = expression_set_nrow(eset)
  let ncols = expression_set_ncol(eset)
  let sample_names = expression_set_sample_names(eset)
  let feature_names = expression_set_feature_names(eset)

  let mut result = "ExpressionSet Summary\n"
  result = result + "-------------------\n"
  result = result + "Features: " + nrows.to_string() + "\n"
  result = result + "Samples: " + ncols.to_string() + "\n"
  result = result + "Annotation: " + eset.annotation + "\n"
  result = result + "\nSample names:\n"

  let mut i = 0
  while i < sample_names.length() && i < 5 {
    result = result + "  - " + sample_names[i] + "\n"
    i = i + 1
  }
  if sample_names.length() > 5 {
    result = result +
      "  ... and " +
      (sample_names.length() - 5).to_string() +
      " more\n"
  }

  result = result + "\nFeature names:\n"
  i = 0
  while i < feature_names.length() && i < 5 {
    result = result + "  - " + feature_names[i] + "\n"
    i = i + 1
  }
  if feature_names.length() > 5 {
    result = result +
      "  ... and " +
      (feature_names.length() - 5).to_string() +
      " more\n"
  }

  result
}

///|
pub fn AnnotatedDataFrame::new(
  data : Map[String, Array[String]],
  var_metadata : Map[String, String],
) -> AnnotatedDataFrame {
  AnnotatedDataFrame::{ data, var_metadata }
}

///|
pub fn annotated_data_frame_nrow(adf : AnnotatedDataFrame) -> Int {
  let keys = adf.data.keys().collect()
  if keys.length() > 0 {
    adf.data[keys[0]].length()
  } else {
    0
  }
}

///|
pub fn annotated_data_frame_ncol(adf : AnnotatedDataFrame) -> Int {
  adf.data.keys().collect().length()
}

///|
pub fn annotated_data_frame_summary(adf : AnnotatedDataFrame) -> String {
  let nrows = annotated_data_frame_nrow(adf)
  let ncols = annotated_data_frame_ncol(adf)
  let keys = adf.data.keys().collect()

  let mut result = "AnnotatedDataFrame Summary\n"
  result = result + "--------------------------\n"
  result = result + "Rows: " + nrows.to_string() + "\n"
  result = result + "Columns: " + ncols.to_string() + "\n"

  if keys.length() > 0 {
    result = result + "\nColumn names:\n"
    let mut i = 0
    while i < keys.length() && i < 10 {
      result = result + "  - " + keys[i] + "\n"
      i = i + 1
    }
    if keys.length() > 10 {
      result = result +
        "  ... and " +
        (keys.length() - 10).to_string() +
        " more\n"
    }
  }

  result
}
