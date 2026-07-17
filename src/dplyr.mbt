///|
/// dplyr - Data manipulation library inspired by R/dplyr.

///|
pub enum DataType {
  IntType
  DoubleType
  StringType
  BoolType
}

///|
pub enum DataValue {
  IntValue(Int)
  DoubleValue(Double)
  StringValue(String)
  BoolValue(Bool)
}

///|
pub struct Column {
  name : String
  dtype : DataType
  values : Array[DataValue]
}

///|
pub struct DataFrame {
  columns : Array[Column]
  nrow : Int
  ncol : Int
}

///|
pub fn data_type_int() -> DataType {
  DataType::IntType
}

///|
pub fn data_type_double() -> DataType {
  DataType::DoubleType
}

///|
pub fn data_type_string() -> DataType {
  DataType::StringType
}

///|
pub fn data_type_bool() -> DataType {
  DataType::BoolType
}

///|
pub fn data_value_int(v : Int) -> DataValue {
  DataValue::IntValue(v)
}

///|
pub fn data_value_double(v : Double) -> DataValue {
  DataValue::DoubleValue(v)
}

///|
pub fn data_value_string(v : String) -> DataValue {
  DataValue::StringValue(v)
}

///|
pub fn data_value_bool(v : Bool) -> DataValue {
  DataValue::BoolValue(v)
}

///|
pub fn column(
  name : String,
  dtype : DataType,
  values : Array[DataValue],
) -> Column {
  Column::{ name, dtype, values }
}

///|
pub fn dataframe(columns : Array[Column]) -> DataFrame {
  let ncol = columns.length()
  let nrow = if ncol == 0 { 0 } else { columns[0].values.length() }
  DataFrame::{ columns, nrow, ncol }
}

///|
pub fn dataframe_nrow(df : DataFrame) -> Int {
  df.nrow
}

///|
pub fn dataframe_ncol(df : DataFrame) -> Int {
  df.ncol
}

///|
pub fn dataframe_colnames(df : DataFrame) -> Array[String] {
  let names : Array[String] = []
  let mut i = 0
  while i < df.ncol {
    names.push(df.columns[i].name)
    i = i + 1
  }
  names
}

///|
pub fn dataframe_get_column(df : DataFrame, name : String) -> Column? {
  let mut i = 0
  while i < df.ncol {
    if df.columns[i].name == name {
      return Option::Some(df.columns[i])
    }
    i = i + 1
  }
  Option::None
}

///|
pub fn dataframe_get_column_index(df : DataFrame, name : String) -> Int? {
  let mut i = 0
  while i < df.ncol {
    if df.columns[i].name == name {
      return Option::Some(i)
    }
    i = i + 1
  }
  Option::None
}

///|
pub fn dataframe_select(
  df : DataFrame,
  column_names : Array[String],
) -> DataFrame {
  let selected : Array[Column] = []
  let mut i = 0
  while i < column_names.length() {
    let col_opt = dataframe_get_column(df, column_names[i])
    if col_opt.is_some() {
      selected.push(col_opt.unwrap())
    }
    i = i + 1
  }
  dataframe(selected)
}

///|
pub fn dataframe_filter(
  df : DataFrame,
  predicate : (DataFrame, Int) -> Bool,
) -> DataFrame {
  let keep : Array[Int] = []
  let mut i = 0
  while i < df.nrow {
    if predicate(df, i) {
      keep.push(i)
    }
    i = i + 1
  }

  let new_cols : Array[Column] = []
  let mut j = 0
  while j < df.ncol {
    let col = df.columns[j]
    let new_values : Array[DataValue] = []
    let mut k = 0
    while k < keep.length() {
      new_values.push(col.values[keep[k]])
      k = k + 1
    }
    new_cols.push(Column::{
      name: col.name,
      dtype: col.dtype,
      values: new_values,
    })
    j = j + 1
  }
  dataframe(new_cols)
}

///|
pub fn dataframe_mutate(
  df : DataFrame,
  new_columns : Array[(String, DataType, (DataFrame, Int) -> DataValue)],
) -> DataFrame {
  let all_cols : Array[Column] = []
  let mut i = 0
  while i < df.ncol {
    all_cols.push(df.columns[i])
    i = i + 1
  }

  let mut j = 0
  while j < new_columns.length() {
    let (name, dtype, func) = new_columns[j]
    let values : Array[DataValue] = []
    let mut k = 0
    while k < df.nrow {
      values.push(func(df, k))
      k = k + 1
    }
    all_cols.push(Column::{ name, dtype, values })
    j = j + 1
  }
  dataframe(all_cols)
}

///|
fn data_value_to_string(v : DataValue) -> String {
  match v {
    DataValue::IntValue(av) => av.to_string()
    DataValue::DoubleValue(av) => av.to_string()
    DataValue::StringValue(av) => av
    DataValue::BoolValue(av) => if av { "true" } else { "false" }
  }
}

///|
fn data_value_less(a : DataValue, b : DataValue) -> Bool {
  match a {
    DataValue::IntValue(av) =>
      match b {
        DataValue::IntValue(bv) => av < bv
        _ => false
      }
    DataValue::DoubleValue(av) =>
      match b {
        DataValue::DoubleValue(bv) => av < bv
        _ => false
      }
    DataValue::StringValue(av) =>
      match b {
        DataValue::StringValue(bv) => av < bv
        _ => false
      }
    DataValue::BoolValue(av) =>
      match b {
        DataValue::BoolValue(bv) => av < bv
        _ => false
      }
  }
}

///|
fn data_value_eq(a : DataValue, b : DataValue) -> Bool {
  match a {
    DataValue::IntValue(av) =>
      match b {
        DataValue::IntValue(bv) => av == bv
        _ => false
      }
    DataValue::DoubleValue(av) =>
      match b {
        DataValue::DoubleValue(bv) => av == bv
        _ => false
      }
    DataValue::StringValue(av) =>
      match b {
        DataValue::StringValue(bv) => av == bv
        _ => false
      }
    DataValue::BoolValue(av) =>
      match b {
        DataValue::BoolValue(bv) => av == bv
        _ => false
      }
  }
}

///|
pub fn dataframe_arrange(
  df : DataFrame,
  column_name : String,
  descending : Bool,
) -> DataFrame {
  let col_idx_opt = dataframe_get_column_index(df, column_name)
  if col_idx_opt.is_none() {
    return df
  }
  let col_idx = col_idx_opt.unwrap()

  let indices : Array[Int] = []
  let mut i = 0
  while i < df.nrow {
    indices.push(i)
    i = i + 1
  }

  let mut swapped = true
  let mut n = df.nrow
  while swapped && n > 0 {
    swapped = false
    let mut j = 0
    while j < n - 1 {
      let val_j = df.columns[col_idx].values[indices[j]]
      let val_j1 = df.columns[col_idx].values[indices[j + 1]]
      let less = data_value_less(val_j, val_j1)
      let should_swap = if descending { !less } else { less }
      if !should_swap && !data_value_eq(val_j, val_j1) {
        let temp = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = temp
        swapped = true
      }
      j = j + 1
    }
    n = n - 1
  }

  let new_cols : Array[Column] = []
  let mut k = 0
  while k < df.ncol {
    let col = df.columns[k]
    let new_values : Array[DataValue] = []
    let mut l = 0
    while l < df.nrow {
      new_values.push(col.values[indices[l]])
      l = l + 1
    }
    new_cols.push(Column::{
      name: col.name,
      dtype: col.dtype,
      values: new_values,
    })
    k = k + 1
  }
  dataframe(new_cols)
}

///|
pub fn dataframe_group_by(
  df : DataFrame,
  group_column : String,
) -> Array[(String, DataFrame)] {
  let col_idx_opt = dataframe_get_column_index(df, group_column)
  if col_idx_opt.is_none() {
    return []
  }
  let col_idx = col_idx_opt.unwrap()

  let group_names : Array[String] = []
  let group_indices : Array[Array[Int]] = []

  let mut i = 0
  while i < df.nrow {
    let key = data_value_to_string(df.columns[col_idx].values[i])
    let mut found = false
    let mut g = 0
    while g < group_names.length() {
      if group_names[g] == key {
        group_indices[g].push(i)
        found = true
        break
      }
      g = g + 1
    }
    if !found {
      group_names.push(key)
      group_indices.push([i])
    }
    i = i + 1
  }

  let result : Array[(String, DataFrame)] = []
  let mut j = 0
  while j < group_names.length() {
    let key = group_names[j]
    let indices = group_indices[j]
    let group_cols : Array[Column] = []
    let mut k = 0
    while k < df.ncol {
      let col = df.columns[k]
      let new_values : Array[DataValue] = []
      let mut l = 0
      while l < indices.length() {
        new_values.push(col.values[indices[l]])
        l = l + 1
      }
      group_cols.push(Column::{
        name: col.name,
        dtype: col.dtype,
        values: new_values,
      })
      k = k + 1
    }
    result.push((key, dataframe(group_cols)))
    j = j + 1
  }
  result
}

///|
pub fn dataframe_summarize(
  df : DataFrame,
  group_column : String,
  agg_columns : Array[(String, (Array[DataValue]) -> DataValue)],
) -> DataFrame {
  let groups = dataframe_group_by(df, group_column)

  if groups.length() == 0 {
    return dataframe([])
  }

  let result_cols : Array[Column] = []

  let group_name_col = Column::{
    name: group_column,
    dtype: StringType,
    values: [],
  }
  result_cols.push(group_name_col)

  let mut i = 0
  while i < agg_columns.length() {
    let (name, agg_func) = agg_columns[i]
    let agg_col = Column::{ name, dtype: DoubleType, values: [] }
    result_cols.push(agg_col)
    i = i + 1
  }

  let mut j = 0
  while j < groups.length() {
    let (key, group_df) = groups[j]
    result_cols[0].values.push(DataValue::StringValue(key))

    let mut k = 0
    while k < agg_columns.length() {
      let (_, agg_func) = agg_columns[k]
      let col_idx_opt = dataframe_get_column_index(group_df, agg_columns[k].0)
      if col_idx_opt.is_some() {
        let col_idx = col_idx_opt.unwrap()
        result_cols[k + 1].values.push(
          agg_func(group_df.columns[col_idx].values),
        )
      } else {
        result_cols[k + 1].values.push(DataValue::DoubleValue(0.0))
      }
      k = k + 1
    }
    j = j + 1
  }

  dataframe(result_cols)
}

///|
pub fn dataframe_inner_join(
  df1 : DataFrame,
  df2 : DataFrame,
  by : String,
) -> DataFrame {
  let col1_idx_opt = dataframe_get_column_index(df1, by)
  let col2_idx_opt = dataframe_get_column_index(df2, by)
  if col1_idx_opt.is_none() {
    return df1
  }
  if col2_idx_opt.is_none() {
    return df2
  }
  let col1_idx = col1_idx_opt.unwrap()
  let col2_idx = col2_idx_opt.unwrap()

  let result_rows : Array[Int] = []
  let result_rows2 : Array[Int] = []

  let mut i = 0
  while i < df1.nrow {
    let val1 = df1.columns[col1_idx].values[i]
    let mut j = 0
    while j < df2.nrow {
      let val2 = df2.columns[col2_idx].values[j]
      if data_value_eq(val1, val2) {
        result_rows.push(i)
        result_rows2.push(j)
      }
      j = j + 1
    }
    i = i + 1
  }

  let new_cols : Array[Column] = []
  let mut k = 0
  while k < df1.ncol {
    let col = df1.columns[k]
    let new_values : Array[DataValue] = []
    let mut l = 0
    while l < result_rows.length() {
      new_values.push(col.values[result_rows[l]])
      l = l + 1
    }
    new_cols.push(Column::{
      name: col.name,
      dtype: col.dtype,
      values: new_values,
    })
    k = k + 1
  }

  let mut m = 0
  while m < df2.ncol {
    if df2.columns[m].name != by {
      let col = df2.columns[m]
      let new_values : Array[DataValue] = []
      let mut l = 0
      while l < result_rows.length() {
        new_values.push(col.values[result_rows2[l]])
        l = l + 1
      }
      new_cols.push(Column::{
        name: col.name,
        dtype: col.dtype,
        values: new_values,
      })
    }
    m = m + 1
  }
  dataframe(new_cols)
}

///|
pub fn dataframe_left_join(
  df1 : DataFrame,
  df2 : DataFrame,
  by : String,
) -> DataFrame {
  let col1_idx_opt = dataframe_get_column_index(df1, by)
  let col2_idx_opt = dataframe_get_column_index(df2, by)
  if col1_idx_opt.is_none() {
    return df1
  }
  if col2_idx_opt.is_none() {
    let all_cols : Array[Column] = []
    let mut i = 0
    while i < df1.ncol {
      all_cols.push(df1.columns[i])
      i = i + 1
    }
    let mut j = 0
    while j < df2.ncol {
      all_cols.push(Column::{
        name: df2.columns[j].name,
        dtype: df2.columns[j].dtype,
        values: [],
      })
      j = j + 1
    }
    return dataframe(all_cols)
  }
  let col1_idx = col1_idx_opt.unwrap()
  let col2_idx = col2_idx_opt.unwrap()

  let result_rows : Array[Int] = []
  let result_rows2 : Array[Int] = []

  let mut i = 0
  while i < df1.nrow {
    let val1 = df1.columns[col1_idx].values[i]
    let mut found = false
    let mut j = 0
    while j < df2.nrow {
      let val2 = df2.columns[col2_idx].values[j]
      if data_value_eq(val1, val2) {
        result_rows.push(i)
        result_rows2.push(j)
        found = true
      }
      j = j + 1
    }
    if !found {
      result_rows.push(i)
      result_rows2.push(-1)
    }
    i = i + 1
  }

  let new_cols : Array[Column] = []
  let mut k = 0
  while k < df1.ncol {
    let col = df1.columns[k]
    let new_values : Array[DataValue] = []
    let mut l = 0
    while l < result_rows.length() {
      new_values.push(col.values[result_rows[l]])
      l = l + 1
    }
    new_cols.push(Column::{
      name: col.name,
      dtype: col.dtype,
      values: new_values,
    })
    k = k + 1
  }

  let mut m = 0
  while m < df2.ncol {
    if df2.columns[m].name != by {
      let col = df2.columns[m]
      let new_values : Array[DataValue] = []
      let mut l = 0
      while l < result_rows.length() {
        if result_rows2[l] >= 0 {
          new_values.push(col.values[result_rows2[l]])
        } else {
          new_values.push(DataValue::IntValue(0))
        }
        l = l + 1
      }
      new_cols.push(Column::{
        name: col.name,
        dtype: col.dtype,
        values: new_values,
      })
    }
    m = m + 1
  }
  dataframe(new_cols)
}
