///|
/// BiocGenerics - Bioconductor generic functions for BioSeqs.
/// Inspired by Bioconductor BiocGenerics package.
/// Provides generic functions shared across Bioconductor packages.

///|
pub enum NAValue {
  NAInt
  NADouble
  NAString
  NABool
}

///|
pub fn na_int() -> Int {
  -2147483648
}

///|
pub fn na_double() -> Double {
  @double.not_a_number
}

///|
pub fn na_string() -> String {
  "NA"
}

///|
pub fn is_na_int(x : Int) -> Bool {
  x == -2147483648
}

///|
pub fn is_na_double(x : Double) -> Bool {
  x.is_nan()
}

///|
pub fn is_na_string(x : String) -> Bool {
  x == "NA" || x == "na" || x == ""
}

///|
pub fn any_na_int(arr : Array[Int]) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if is_na_int(arr[i]) {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn any_na_double(arr : Array[Double]) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if is_na_double(arr[i]) {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn any_na_string(arr : Array[String]) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if is_na_string(arr[i]) {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn na_omit_int(arr : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < arr.length() {
    if !is_na_int(arr[i]) {
      result.push(arr[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn na_omit_double(arr : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < arr.length() {
    if !is_na_double(arr[i]) {
      result.push(arr[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn na_omit_string(arr : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < arr.length() {
    if !is_na_string(arr[i]) {
      result.push(arr[i])
    }
    i = i + 1
  }
  result
}

///|
pub fn which_int(arr : Array[Int], value : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < arr.length() {
    if arr[i] == value {
      result.push(i)
    }
    i = i + 1
  }
  result
}

///|
pub fn which_double(arr : Array[Double], value : Double) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < arr.length() {
    if arr[i] == value {
      result.push(i)
    }
    i = i + 1
  }
  result
}

///|
pub fn which_bool(arr : Array[Bool]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < arr.length() {
    if arr[i] {
      result.push(i)
    }
    i = i + 1
  }
  result
}

///|
pub fn order_int(arr : Array[Int]) -> Array[Int] {
  let n = arr.length()
  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  
  let mut sorted = false
  while !sorted {
    sorted = true
    i = 0
    while i < n - 1 {
      let j = i + 1
      if arr[indices[i]] > arr[indices[j]] {
        let temp = indices[i]
        indices[i] = indices[j]
        indices[j] = temp
        sorted = false
      }
      i = i + 1
    }
  }
  
  indices
}

///|
pub fn order_double(arr : Array[Double]) -> Array[Int] {
  let n = arr.length()
  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  
  let mut sorted = false
  while !sorted {
    sorted = true
    i = 0
    while i < n - 1 {
      let j = i + 1
      if arr[indices[i]] > arr[indices[j]] {
        let temp = indices[i]
        indices[i] = indices[j]
        indices[j] = temp
        sorted = false
      }
      i = i + 1
    }
  }
  
  indices
}

///|
pub fn sort_int(arr : Array[Int]) -> Array[Int] {
  let indices = order_int(arr)
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < indices.length() {
    result.push(arr[indices[i]])
    i = i + 1
  }
  result
}

///|
pub fn sort_double(arr : Array[Double]) -> Array[Double] {
  let indices = order_double(arr)
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < indices.length() {
    result.push(arr[indices[i]])
    i = i + 1
  }
  result
}

///|
pub fn unique_int(arr : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let seen = Map([], capacity=arr.length())
  
  let mut i = 0
  while i < arr.length() {
    if !seen.contains(arr[i].to_string()) {
      seen.set(arr[i].to_string(), true)
      result.push(arr[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn unique_double(arr : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let seen = Map([], capacity=arr.length())
  
  let mut i = 0
  while i < arr.length() {
    if !seen.contains(arr[i].to_string()) {
      seen.set(arr[i].to_string(), true)
      result.push(arr[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn unique_string(arr : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen = Map([], capacity=arr.length())
  
  let mut i = 0
  while i < arr.length() {
    if !seen.contains(arr[i]) {
      seen.set(arr[i], true)
      result.push(arr[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn match_int(x : Array[Int], table : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let table_map = Map([], capacity=table.length())
  
  let mut i = 0
  while i < table.length() {
    table_map.set(table[i].to_string(), i)
    i = i + 1
  }
  
  i = 0
  while i < x.length() {
    if table_map.contains(x[i].to_string()) {
      result.push(table_map.get(x[i].to_string()).unwrap())
    } else {
      result.push(na_int())
    }
    i = i + 1
  }
  
  result
}

///|
pub fn match_string(x : Array[String], table : Array[String]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let table_map = Map([], capacity=table.length())
  
  let mut i = 0
  while i < table.length() {
    table_map.set(table[i], i)
    i = i + 1
  }
  
  i = 0
  while i < x.length() {
    if table_map.contains(x[i]) {
      result.push(table_map.get(x[i]).unwrap())
    } else {
      result.push(na_int())
    }
    i = i + 1
  }
  
  result
}

///|
pub fn intersect_int(x : Array[Int], y : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let y_set = Map([], capacity=y.length())
  
  let mut i = 0
  while i < y.length() {
    y_set.set(y[i].to_string(), true)
    i = i + 1
  }
  
  i = 0
  while i < x.length() {
    if y_set.contains(x[i].to_string()) {
      result.push(x[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn intersect_string(x : Array[String], y : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let y_set = Map([], capacity=y.length())
  
  let mut i = 0
  while i < y.length() {
    y_set.set(y[i], true)
    i = i + 1
  }
  
  i = 0
  while i < x.length() {
    if y_set.contains(x[i]) {
      result.push(x[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn union_int(x : Array[Int], y : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let seen = Map([], capacity=x.length() + y.length())
  
  let mut i = 0
  while i < x.length() {
    if !seen.contains(x[i].to_string()) {
      seen.set(x[i].to_string(), true)
      result.push(x[i])
    }
    i = i + 1
  }
  
  i = 0
  while i < y.length() {
    if !seen.contains(y[i].to_string()) {
      seen.set(y[i].to_string(), true)
      result.push(y[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn union_string(x : Array[String], y : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen = Map([], capacity=x.length() + y.length())
  
  let mut i = 0
  while i < x.length() {
    if !seen.contains(x[i]) {
      seen.set(x[i], true)
      result.push(x[i])
    }
    i = i + 1
  }
  
  i = 0
  while i < y.length() {
    if !seen.contains(y[i]) {
      seen.set(y[i], true)
      result.push(y[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn setdiff_int(x : Array[Int], y : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let y_set = Map([], capacity=y.length())
  
  let mut i = 0
  while i < y.length() {
    y_set.set(y[i].to_string(), true)
    i = i + 1
  }
  
  i = 0
  while i < x.length() {
    if !y_set.contains(x[i].to_string()) {
      result.push(x[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn setdiff_string(x : Array[String], y : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let y_set = Map([], capacity=y.length())
  
  let mut i = 0
  while i < y.length() {
    y_set.set(y[i], true)
    i = i + 1
  }
  
  i = 0
  while i < x.length() {
    if !y_set.contains(x[i]) {
      result.push(x[i])
    }
    i = i + 1
  }
  
  result
}

///|
pub fn table_int(arr : Array[Int]) -> Map[String, Int] {
  let result = Map([], capacity=arr.length())
  
  let mut i = 0
  while i < arr.length() {
    let key = arr[i].to_string()
    if result.contains(key) {
      result.set(key, result.get(key).unwrap() + 1)
    } else {
      result.set(key, 1)
    }
    i = i + 1
  }
  
  result
}

///|
pub fn table_string(arr : Array[String]) -> Map[String, Int] {
  let result = Map([], capacity=arr.length())
  
  let mut i = 0
  while i < arr.length() {
    let key = arr[i]
    if result.contains(key) {
      result.set(key, result.get(key).unwrap() + 1)
    } else {
      result.set(key, 1)
    }
    i = i + 1
  }
  
  result
}

///|
pub fn rev_int(arr : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = arr.length() - 1
  while i >= 0 {
    result.push(arr[i])
    i = i - 1
  }
  result
}

///|
pub fn rev_double(arr : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = arr.length() - 1
  while i >= 0 {
    result.push(arr[i])
    i = i - 1
  }
  result
}

///|
pub fn rev_string(arr : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut i = arr.length() - 1
  while i >= 0 {
    result.push(arr[i])
    i = i - 1
  }
  result
}

///|
pub fn rep_int(x : Int, times : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < times {
    result.push(x)
    i = i + 1
  }
  result
}

///|
pub fn rep_double(x : Double, times : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < times {
    result.push(x)
    i = i + 1
  }
  result
}

///|
pub fn rep_string(x : String, times : Int) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < times {
    result.push(x)
    i = i + 1
  }
  result
}

///|
pub fn seq_int(from : Int, to : Int, by? : Int) -> Array[Int] {
  let step = if by is Some(_) { by.unwrap() } else { 1 }
  let result : Array[Int] = Array::new()
  
  if step > 0 {
    let mut i = from
    while i <= to {
      result.push(i)
      i = i + step
    }
  } else if step < 0 {
    let mut i = from
    while i >= to {
      result.push(i)
      i = i + step
    }
  }
  
  result
}

///|
pub fn seq_double(from : Double, to : Double, by? : Double) -> Array[Double] {
  let step = if by is Some(_) { by.unwrap() } else { 1.0 }
  let result : Array[Double] = Array::new()
  
  if step > 0.0 {
    let mut i = from
    while i <= to {
      result.push(i)
      i = i + step
    }
  } else if step < 0.0 {
    let mut i = from
    while i >= to {
      result.push(i)
      i = i + step
    }
  }
  
  result
}

///|
pub fn head_int(arr : Array[Int], n : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let limit = if n < arr.length() { n } else { arr.length() }
  let mut i = 0
  while i < limit {
    result.push(arr[i])
    i = i + 1
  }
  result
}

///|
pub fn head_double(arr : Array[Double], n : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let limit = if n < arr.length() { n } else { arr.length() }
  let mut i = 0
  while i < limit {
    result.push(arr[i])
    i = i + 1
  }
  result
}

///|
pub fn head_string(arr : Array[String], n : Int) -> Array[String] {
  let result : Array[String] = Array::new()
  let limit = if n < arr.length() { n } else { arr.length() }
  let mut i = 0
  while i < limit {
    result.push(arr[i])
    i = i + 1
  }
  result
}

///|
pub fn tail_int(arr : Array[Int], n : Int) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let start = if n < arr.length() { arr.length() - n } else { 0 }
  let mut i = start
  while i < arr.length() {
    result.push(arr[i])
    i = i + 1
  }
  result
}

///|
pub fn tail_double(arr : Array[Double], n : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let start = if n < arr.length() { arr.length() - n } else { 0 }
  let mut i = start
  while i < arr.length() {
    result.push(arr[i])
    i = i + 1
  }
  result
}

///|
pub fn tail_string(arr : Array[String], n : Int) -> Array[String] {
  let result : Array[String] = Array::new()
  let start = if n < arr.length() { arr.length() - n } else { 0 }
  let mut i = start
  while i < arr.length() {
    result.push(arr[i])
    i = i + 1
  }
  result
}

///|
pub fn length_int(arr : Array[Int]) -> Int {
  arr.length()
}

///|
pub fn length_double(arr : Array[Double]) -> Int {
  arr.length()
}

///|
pub fn length_string(arr : Array[String]) -> Int {
  arr.length()
}

///|
pub fn names(arr : Array[String]) -> Array[String] {
  arr
}

///|
pub fn names_set(arr : Array[String], names : Array[String]) -> Array[String] {
  let result = arr.copy()
  let mut i = 0
  while i < names.length() && i < result.length() {
    result[i] = names[i]
    i = i + 1
  }
  result
}

///|
pub fn rownames(arr : Array[String]) -> Array[String] {
  arr
}

///|
pub fn colnames(arr : Array[String]) -> Array[String] {
  arr
}

///|
pub fn cbind(arrays : Array[Array[Double]]) -> Array[Array[Double]] {
  if arrays.length() == 0 {
    return []
  }
  
  let nrow = arrays[0].length()
  let result : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < nrow {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < arrays.length() {
      if i < arrays[j].length() {
        row.push(arrays[j][i])
      } else {
        row.push(na_double())
      }
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  
  result
}

///|
pub fn rbind(arrays : Array[Array[Double]]) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < arrays.length() {
    result.push(arrays[i].copy())
    i = i + 1
  }
  
  result
}
