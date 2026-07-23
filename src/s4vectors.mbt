///|
/// S4Vectors module, inspired by Bioconductor S4Vectors package.
/// Provides foundational data structures for Bioconductor: Rle, DataFrame, and Hits.

///|
pub struct Rle {
  values : Array[String]
  lengths : Array[Int]
}

///|
pub fn Rle::new(values : Array[String], lengths : Array[Int]) -> Rle {
  Rle::{ values, lengths }
}

///|
pub fn Rle::from_vector(vec : Array[String]) -> Rle {
  let values : Array[String] = Array::new()
  let lengths : Array[Int] = Array::new()
  
  if vec.length() == 0 {
    return Rle::{ values, lengths }
  }
  
  let mut current_value = vec[0]
  let mut current_length = 1
  
  let mut i = 1
  while i < vec.length() {
    if vec[i] == current_value {
      current_length = current_length + 1
    } else {
      values.push(current_value)
      lengths.push(current_length)
      current_value = vec[i]
      current_length = 1
    }
    i = i + 1
  }
  
  values.push(current_value)
  lengths.push(current_length)
  
  Rle::{ values, lengths }
}

///|
pub fn Rle::length(self : Rle) -> Int {
  let mut total = 0
  let mut i = 0
  while i < self.lengths.length() {
    total = total + self.lengths[i]
    i = i + 1
  }
  total
}

///|
pub fn Rle::decode(self : Rle) -> Array[String] {
  let result : Array[String] = Array::new()
  
  let mut i = 0
  while i < self.values.length() {
    let mut j = 0
    while j < self.lengths[i] {
      result.push(self.values[i])
      j = j + 1
    }
    i = i + 1
  }
  
  result
}

///|
pub fn Rle::get(self : Rle, idx : Int) -> String {
  let mut pos = 0
  
  let mut i = 0
  while i < self.lengths.length() {
    if pos + self.lengths[i] > idx {
      return self.values[i]
    }
    pos = pos + self.lengths[i]
    i = i + 1
  }
  
  ""
}

///|
pub fn Rle::compress_ratio(self : Rle) -> Double {
  let original_size = self.length()
  let compressed_size = self.values.length() + self.lengths.length()
  
  if original_size == 0 {
    return 0.0
  }
  
  original_size.to_double() / compressed_size.to_double()
}

///|
pub struct S4DataFrameColumn {
  name : String
  values : Array[String]
}

///|
pub fn S4DataFrameColumn::new(name : String, values : Array[String]) -> S4DataFrameColumn {
  S4DataFrameColumn::{ name, values }
}

///|
pub struct S4DataFrame {
  columns : Array[S4DataFrameColumn]
  row_names : Array[String]
}

///|
pub fn S4DataFrame::new(columns : Array[S4DataFrameColumn], row_names : Array[String]) -> S4DataFrame {
  S4DataFrame::{ columns, row_names }
}

///|
pub fn S4DataFrame::nrow(self : S4DataFrame) -> Int {
  if self.columns.length() == 0 {
    return 0
  }
  self.columns[0].values.length()
}

///|
pub fn S4DataFrame::ncol(self : S4DataFrame) -> Int {
  self.columns.length()
}

///|
pub fn S4DataFrame::colnames(self : S4DataFrame) -> Array[String] {
  let names : Array[String] = Array::new()
  
  let mut i = 0
  while i < self.columns.length() {
    names.push(self.columns[i].name)
    i = i + 1
  }
  
  names
}

///|
pub fn S4DataFrame::get_column(self : S4DataFrame, name : String) -> S4DataFrameColumn? {
  let mut i = 0
  while i < self.columns.length() {
    if self.columns[i].name == name {
      return Some(self.columns[i])
    }
    i = i + 1
  }
  
  None
}

///|
pub fn S4DataFrame::get_row(self : S4DataFrame, idx : Int) -> Array[String] {
  let result : Array[String] = Array::new()
  
  let mut i = 0
  while i < self.columns.length() {
    if idx < self.columns[i].values.length() {
      result.push(self.columns[i].values[idx])
    } else {
      result.push("")
    }
    i = i + 1
  }
  
  result
}

///|
pub fn S4DataFrame::subset_rows(self : S4DataFrame, indices : Array[Int]) -> S4DataFrame {
  let new_columns : Array[S4DataFrameColumn] = Array::new()
  let new_row_names : Array[String] = Array::new()
  
  let mut i = 0
  while i < self.columns.length() {
    let col_values : Array[String] = Array::new()
    let mut j = 0
    while j < indices.length() {
      let idx = indices[j]
      if idx < self.columns[i].values.length() {
        col_values.push(self.columns[i].values[idx])
      }
      j = j + 1
    }
    new_columns.push(S4DataFrameColumn::new(self.columns[i].name, col_values))
    i = i + 1
  }
  
  let mut j = 0
  while j < indices.length() {
    let idx = indices[j]
    if idx < self.row_names.length() {
      new_row_names.push(self.row_names[idx])
    }
    j = j + 1
  }
  
  S4DataFrame::new(new_columns, new_row_names)
}

///|
pub struct Hits {
  query_hits : Array[Int]
  subject_hits : Array[Int]
  query_length : Int
  subject_length : Int
}

///|
pub fn Hits::new(query_hits : Array[Int], subject_hits : Array[Int], query_length : Int, subject_length : Int) -> Hits {
  Hits::{ query_hits, subject_hits, query_length, subject_length }
}

///|
pub fn Hits::n_hits(self : Hits) -> Int {
  self.query_hits.length()
}

///|
pub fn Hits::count_query_hits(self : Hits) -> Array[Int] {
  let counts : Array[Int] = Array::make(self.query_length, 0)
  
  let mut i = 0
  while i < self.query_hits.length() {
    let idx = self.query_hits[i]
    if idx >= 0 && idx < self.query_length {
      counts[idx] = counts[idx] + 1
    }
    i = i + 1
  }
  
  counts
}

///|
pub fn Hits::count_subject_hits(self : Hits) -> Array[Int] {
  let counts : Array[Int] = Array::make(self.subject_length, 0)
  
  let mut i = 0
  while i < self.subject_hits.length() {
    let idx = self.subject_hits[i]
    if idx >= 0 && idx < self.subject_length {
      counts[idx] = counts[idx] + 1
    }
    i = i + 1
  }
  
  counts
}

///|
pub fn Hits::find_hits(self : Hits, query_idx : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  
  let mut i = 0
  while i < self.query_hits.length() {
    if self.query_hits[i] == query_idx {
      result.push(self.subject_hits[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn Hits::as_matrix(self : Hits) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = Array::new()
  
  let mut i = 0
  while i < self.query_hits.length() {
    result.push((self.query_hits[i], self.subject_hits[i]))
    i = i + 1
  }
  
  result
}

///|
pub fn create_example_rle() -> Rle {
  let vec = ["A", "A", "A", "B", "B", "C", "C", "C", "C", "A", "A"]
  Rle::from_vector(vec)
}

///|
pub fn create_example_dataframe() -> S4DataFrame {
  let col1 = S4DataFrameColumn::new("gene_id", ["gene1", "gene2", "gene3", "gene4", "gene5"])
  let col2 = S4DataFrameColumn::new("expression", ["10.5", "25.3", "5.8", "100.2", "15.7"])
  let col3 = S4DataFrameColumn::new("p_value", ["0.001", "0.05", "0.8", "0.0001", "0.1"])
  
  S4DataFrame::new([col1, col2, col3], ["row1", "row2", "row3", "row4", "row5"])
}

///|
pub fn create_example_hits() -> Hits {
  let query_hits = [0, 0, 1, 2, 2, 2, 3]
  let subject_hits = [1, 3, 0, 1, 2, 4, 2]
  
  Hits::new(query_hits, subject_hits, 4, 5)
}