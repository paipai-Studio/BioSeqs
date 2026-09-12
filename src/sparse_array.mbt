///|
/// Multidimensional sparse arrays inspired by Bioconductor SparseArray.
///
/// Values are stored in a canonical coordinate (COO) representation. Every
/// public constructor validates dimensions and coordinates, sorts entries,
/// merges duplicate coordinates, and removes explicit zero values.

///|
pub suberror SparseArrayError {
  SparseArrayError(String)
}

///|
pub struct SparseArrayEntry {
  coordinates : Array[Int]
  value : Double
}

///|
pub struct SparseArray {
  dimensions : Array[Int]
  entries : Array[SparseArrayEntry]
}

///|
pub fn SparseArrayEntry::coordinates(self : SparseArrayEntry) -> Array[Int] {
  self.coordinates.copy()
}

///|
pub fn SparseArrayEntry::value(self : SparseArrayEntry) -> Double {
  self.value
}

///|
fn sparse_array_compare_coordinates(
  left : Array[Int],
  right : Array[Int],
) -> Int {
  let limit = if left.length() < right.length() {
    left.length()
  } else {
    right.length()
  }
  let mut index = 0
  while index < limit {
    if left[index] < right[index] {
      return -1
    }
    if left[index] > right[index] {
      return 1
    }
    index = index + 1
  }
  if left.length() < right.length() {
    -1
  } else if left.length() > right.length() {
    1
  } else {
    0
  }
}

///|
fn sparse_array_validate_dimensions(
  dimensions : Array[Int],
) -> Unit raise SparseArrayError {
  if dimensions.length() == 0 {
    raise SparseArrayError("SparseArray requires at least one dimension")
  }
  let mut dimension = 0
  while dimension < dimensions.length() {
    if dimensions[dimension] < 0 {
      raise SparseArrayError("SparseArray dimensions must be non-negative")
    }
    dimension = dimension + 1
  }
}

///|
fn sparse_array_length(dimensions : Array[Int]) -> Int raise SparseArrayError {
  let mut result = 1
  let mut dimension = 0
  while dimension < dimensions.length() {
    let extent = dimensions[dimension]
    if extent == 0 {
      return 0
    }
    let next = result * extent
    if next / extent != result {
      raise SparseArrayError("SparseArray dimensions are too large")
    }
    result = next
    dimension = dimension + 1
  }
  result
}

///|
fn sparse_array_validate_coordinate(
  dimensions : Array[Int],
  coordinates : Array[Int],
) -> Unit raise SparseArrayError {
  if coordinates.length() != dimensions.length() {
    raise SparseArrayError(
      "SparseArray coordinate rank does not match its dimensions",
    )
  }
  let mut dimension = 0
  while dimension < dimensions.length() {
    let coordinate = coordinates[dimension]
    if coordinate < 0 || coordinate >= dimensions[dimension] {
      raise SparseArrayError(
        "SparseArray coordinate is outside dimension " + dimension.to_string(),
      )
    }
    dimension = dimension + 1
  }
}

///|
fn sparse_array_build(
  dimensions : Array[Int],
  entries : Array[SparseArrayEntry],
) -> SparseArray raise SparseArrayError {
  sparse_array_validate_dimensions(dimensions)
  ignore(sparse_array_length(dimensions))
  let sorted : Array[SparseArrayEntry] = []
  let mut index = 0
  while index < entries.length() {
    sparse_array_validate_coordinate(dimensions, entries[index].coordinates)
    sorted.push(SparseArrayEntry::{
      coordinates: entries[index].coordinates.copy(),
      value: entries[index].value,
    })
    index = index + 1
  }
  sorted.sort_by(fn(left : SparseArrayEntry, right : SparseArrayEntry) -> Int {
    sparse_array_compare_coordinates(left.coordinates, right.coordinates)
  })
  let canonical : Array[SparseArrayEntry] = []
  index = 0
  while index < sorted.length() {
    let coordinates = sorted[index].coordinates.copy()
    let mut value = sorted[index].value
    index = index + 1
    while index < sorted.length() &&
          sparse_array_compare_coordinates(
            coordinates,
            sorted[index].coordinates,
          ) ==
          0 {
      value = value + sorted[index].value
      index = index + 1
    }
    if value != 0.0 {
      canonical.push(SparseArrayEntry::{ coordinates, value, })
    }
  }
  SparseArray::{ dimensions: dimensions.copy(), entries: canonical, }
}

///|
/// Construct a sparse array from coordinates and values.
///
/// Coordinates are zero-based. Duplicate coordinates are added together and
/// entries whose resulting value is zero are omitted.
pub fn SparseArray::from_coo(
  dimensions : Array[Int],
  coordinates : Array[Array[Int]],
  values : Array[Double],
) -> SparseArray raise SparseArrayError {
  if coordinates.length() != values.length() {
    raise SparseArrayError(
      "SparseArray coordinates and values must have the same length",
    )
  }
  let entries : Array[SparseArrayEntry] = []
  let mut index = 0
  while index < values.length() {
    entries.push(SparseArrayEntry::{
      coordinates: coordinates[index].copy(),
      value: values[index],
    })
    index = index + 1
  }
  sparse_array_build(dimensions, entries)
}

///|
/// Construct an all-zero sparse array.
pub fn SparseArray::zeros(
  dimensions : Array[Int],
) -> SparseArray raise SparseArrayError {
  sparse_array_build(dimensions, [])
}

///|
fn sparse_array_linear_to_coordinates(
  linear_index : Int,
  dimensions : Array[Int],
) -> Array[Int] {
  let coordinates : Array[Int] = []
  let mut remaining = linear_index
  let mut dimension = 0
  while dimension < dimensions.length() {
    coordinates.push(remaining % dimensions[dimension])
    remaining = remaining / dimensions[dimension]
    dimension = dimension + 1
  }
  coordinates
}

///|
fn sparse_array_coordinates_to_linear(
  coordinates : Array[Int],
  dimensions : Array[Int],
) -> Int {
  let mut linear_index = 0
  let mut stride = 1
  let mut dimension = 0
  while dimension < dimensions.length() {
    linear_index = linear_index + coordinates[dimension] * stride
    stride = stride * dimensions[dimension]
    dimension = dimension + 1
  }
  linear_index
}

///|
/// Construct from column-major flat data, matching R array layout.
pub fn SparseArray::from_flat(
  dimensions : Array[Int],
  data : Array[Double],
) -> SparseArray raise SparseArrayError {
  sparse_array_validate_dimensions(dimensions)
  let expected = sparse_array_length(dimensions)
  if data.length() != expected {
    raise SparseArrayError(
      "Flat data length does not match SparseArray dimensions",
    )
  }
  let entries : Array[SparseArrayEntry] = []
  let mut index = 0
  while index < data.length() {
    if data[index] != 0.0 {
      entries.push(SparseArrayEntry::{
        coordinates: sparse_array_linear_to_coordinates(index, dimensions),
        value: data[index],
      })
    }
    index = index + 1
  }
  sparse_array_build(dimensions, entries)
}

///|
/// Construct a two-dimensional sparse array from a row-major MoonBit matrix.
pub fn SparseArray::from_matrix(
  data : Array[Array[Double]],
) -> SparseArray raise SparseArrayError {
  let rows = data.length()
  let columns = if rows == 0 { 0 } else { data[0].length() }
  let coordinates : Array[Array[Int]] = []
  let values : Array[Double] = []
  let mut row = 0
  while row < rows {
    if data[row].length() != columns {
      raise SparseArrayError("SparseArray input matrix must be rectangular")
    }
    let mut column = 0
    while column < columns {
      if data[row][column] != 0.0 {
        coordinates.push([row, column])
        values.push(data[row][column])
      }
      column = column + 1
    }
    row = row + 1
  }
  SparseArray::from_coo([rows, columns], coordinates, values)
}

///|
pub fn SparseArray::ndim(self : SparseArray) -> Int {
  self.dimensions.length()
}

///|
pub fn SparseArray::dim(self : SparseArray) -> Array[Int] {
  self.dimensions.copy()
}

///|
pub fn SparseArray::length(self : SparseArray) -> Int raise SparseArrayError {
  sparse_array_length(self.dimensions)
}

///|
pub fn SparseArray::nnzero(self : SparseArray) -> Int {
  self.entries.length()
}

///|
pub fn SparseArray::density(
  self : SparseArray,
) -> Double raise SparseArrayError {
  let total = self.length()
  if total == 0 {
    0.0
  } else {
    self.entries.length().to_double() / total.to_double()
  }
}

///|
pub fn SparseArray::is_sparse(self : SparseArray) -> Bool {
  true
}

///|
pub fn SparseArray::nzcoordinates(self : SparseArray) -> Array[Array[Int]] {
  let result : Array[Array[Int]] = []
  let mut index = 0
  while index < self.entries.length() {
    result.push(self.entries[index].coordinates.copy())
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::nzvalues(self : SparseArray) -> Array[Double] {
  let result : Array[Double] = []
  let mut index = 0
  while index < self.entries.length() {
    result.push(self.entries[index].value)
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::nonzero_entries(
  self : SparseArray,
) -> Array[SparseArrayEntry] {
  let result : Array[SparseArrayEntry] = []
  let mut index = 0
  while index < self.entries.length() {
    result.push(SparseArrayEntry::{
      coordinates: self.entries[index].coordinates.copy(),
      value: self.entries[index].value,
    })
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::get(
  self : SparseArray,
  coordinates : Array[Int],
) -> Double raise SparseArrayError {
  sparse_array_validate_coordinate(self.dimensions, coordinates)
  let mut low = 0
  let mut high = self.entries.length()
  while low < high {
    let middle = low + (high - low) / 2
    let comparison = sparse_array_compare_coordinates(
      self.entries[middle].coordinates,
      coordinates,
    )
    if comparison < 0 {
      low = middle + 1
    } else {
      high = middle
    }
  }
  if low < self.entries.length() &&
    sparse_array_compare_coordinates(self.entries[low].coordinates, coordinates) ==
    0 {
    self.entries[low].value
  } else {
    0.0
  }
}

///|
pub fn SparseArray::get_linear(
  self : SparseArray,
  linear_index : Int,
) -> Double raise SparseArrayError {
  let total = self.length()
  if linear_index < 0 || linear_index >= total {
    raise SparseArrayError("SparseArray linear index is out of bounds")
  }
  self.get(sparse_array_linear_to_coordinates(linear_index, self.dimensions))
}

///|
/// Return a new array with one coordinate replaced.
pub fn SparseArray::with_value(
  self : SparseArray,
  coordinates : Array[Int],
  value : Double,
) -> SparseArray raise SparseArrayError {
  sparse_array_validate_coordinate(self.dimensions, coordinates)
  let entries : Array[SparseArrayEntry] = []
  let mut index = 0
  while index < self.entries.length() {
    if sparse_array_compare_coordinates(
        self.entries[index].coordinates,
        coordinates,
      ) !=
      0 {
      entries.push(SparseArrayEntry::{
        coordinates: self.entries[index].coordinates.copy(),
        value: self.entries[index].value,
      })
    }
    index = index + 1
  }
  if value != 0.0 {
    entries.push(SparseArrayEntry::{ coordinates: coordinates.copy(), value, })
  }
  sparse_array_build(self.dimensions, entries)
}

///|
/// Return a new array after replacing several coordinates in order.
pub fn SparseArray::with_values(
  self : SparseArray,
  coordinates : Array[Array[Int]],
  values : Array[Double],
) -> SparseArray raise SparseArrayError {
  if coordinates.length() != values.length() {
    raise SparseArrayError(
      "Replacement coordinates and values must have the same length",
    )
  }
  let mut result = self
  let mut index = 0
  while index < values.length() {
    result = result.with_value(coordinates[index], values[index])
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::to_flat(
  self : SparseArray,
) -> Array[Double] raise SparseArrayError {
  let result = Array::make(self.length(), 0.0)
  let mut index = 0
  while index < self.entries.length() {
    let linear_index = sparse_array_coordinates_to_linear(
      self.entries[index].coordinates,
      self.dimensions,
    )
    result[linear_index] = self.entries[index].value
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::to_dense_matrix(
  self : SparseArray,
) -> Array[Array[Double]] raise SparseArrayError {
  if self.ndim() != 2 {
    raise SparseArrayError(
      "Dense matrix conversion requires a two-dimensional SparseArray",
    )
  }
  let result : Array[Array[Double]] = []
  let mut row = 0
  while row < self.dimensions[0] {
    result.push(Array::make(self.dimensions[1], 0.0))
    row = row + 1
  }
  let mut index = 0
  while index < self.entries.length() {
    let coordinates = self.entries[index].coordinates
    result[coordinates[0]][coordinates[1]] = self.entries[index].value
    index = index + 1
  }
  result
}

///|
/// Subset every dimension using zero-based index vectors.
///
/// The output dimension sizes and order follow the supplied vectors.
/// Repeated indices are preserved without materializing implicit zero values.
pub fn SparseArray::subset(
  self : SparseArray,
  indices : Array[Array[Int]],
) -> SparseArray raise SparseArrayError {
  if indices.length() != self.ndim() {
    raise SparseArrayError(
      "SparseArray subset requires one index vector per dimension",
    )
  }
  let output_dimensions : Array[Int] = []
  let mut dimension = 0
  while dimension < self.ndim() {
    output_dimensions.push(indices[dimension].length())
    let mut index = 0
    while index < indices[dimension].length() {
      let selected = indices[dimension][index]
      if selected < 0 || selected >= self.dimensions[dimension] {
        raise SparseArrayError(
          "SparseArray subset index is out of bounds in dimension " +
          dimension.to_string(),
        )
      }
      index = index + 1
    }
    dimension = dimension + 1
  }
  let output_entries : Array[SparseArrayEntry] = []
  let mut entry_index = 0
  while entry_index < self.entries.length() {
    let matches : Array[Array[Int]] = []
    let mut include = true
    dimension = 0
    while dimension < self.ndim() {
      let positions : Array[Int] = []
      let mut output_index = 0
      while output_index < indices[dimension].length() {
        if indices[dimension][output_index] ==
          self.entries[entry_index].coordinates[dimension] {
          positions.push(output_index)
        }
        output_index = output_index + 1
      }
      if positions.length() == 0 {
        include = false
      }
      matches.push(positions)
      dimension = dimension + 1
    }
    if include {
      let mut combination_count = 1
      dimension = 0
      while dimension < matches.length() {
        combination_count = combination_count * matches[dimension].length()
        dimension = dimension + 1
      }
      let mut combination = 0
      while combination < combination_count {
        let output_coordinates : Array[Int] = []
        let mut remaining = combination
        dimension = 0
        while dimension < matches.length() {
          let count = matches[dimension].length()
          output_coordinates.push(matches[dimension][remaining % count])
          remaining = remaining / count
          dimension = dimension + 1
        }
        output_entries.push(SparseArrayEntry::{
          coordinates: output_coordinates,
          value: self.entries[entry_index].value,
        })
        combination = combination + 1
      }
    }
    entry_index = entry_index + 1
  }
  sparse_array_build(output_dimensions, output_entries)
}

///|
/// Slice dimensions with zero-based, end-exclusive ranges.
pub fn SparseArray::slice(
  self : SparseArray,
  starts : Array[Int],
  ends : Array[Int],
) -> SparseArray raise SparseArrayError {
  if starts.length() != self.ndim() || ends.length() != self.ndim() {
    raise SparseArrayError("SparseArray slice bounds must match the array rank")
  }
  let indices : Array[Array[Int]] = []
  let mut dimension = 0
  while dimension < self.ndim() {
    if starts[dimension] < 0 ||
      ends[dimension] < starts[dimension] ||
      ends[dimension] > self.dimensions[dimension] {
      raise SparseArrayError(
        "SparseArray slice bounds are invalid in dimension " +
        dimension.to_string(),
      )
    }
    let selected : Array[Int] = []
    let mut index = starts[dimension]
    while index < ends[dimension] {
      selected.push(index)
      index = index + 1
    }
    indices.push(selected)
    dimension = dimension + 1
  }
  self.subset(indices)
}

///|
/// Permute dimensions. Output dimension i is input dimension permutation[i].
pub fn SparseArray::aperm(
  self : SparseArray,
  permutation : Array[Int],
) -> SparseArray raise SparseArrayError {
  if permutation.length() != self.ndim() {
    raise SparseArrayError("SparseArray permutation length must match its rank")
  }
  let seen = Array::make(self.ndim(), false)
  let output_dimensions : Array[Int] = []
  let mut dimension = 0
  while dimension < permutation.length() {
    let source = permutation[dimension]
    if source < 0 || source >= self.ndim() || seen[source] {
      raise SparseArrayError(
        "SparseArray permutation must contain each dimension exactly once",
      )
    }
    seen[source] = true
    output_dimensions.push(self.dimensions[source])
    dimension = dimension + 1
  }
  let output_entries : Array[SparseArrayEntry] = []
  let mut entry_index = 0
  while entry_index < self.entries.length() {
    let output_coordinates : Array[Int] = []
    dimension = 0
    while dimension < permutation.length() {
      output_coordinates.push(
        self.entries[entry_index].coordinates[permutation[dimension]],
      )
      dimension = dimension + 1
    }
    output_entries.push(SparseArrayEntry::{
      coordinates: output_coordinates,
      value: self.entries[entry_index].value,
    })
    entry_index = entry_index + 1
  }
  sparse_array_build(output_dimensions, output_entries)
}

///|
pub fn SparseArray::transpose(
  self : SparseArray,
) -> SparseArray raise SparseArrayError {
  if self.ndim() != 2 {
    raise SparseArrayError(
      "SparseArray transpose requires a two-dimensional array",
    )
  }
  self.aperm([1, 0])
}

///|
/// Bind two arrays along one dimension.
pub fn SparseArray::bind(
  self : SparseArray,
  other : SparseArray,
  along : Int,
) -> SparseArray raise SparseArrayError {
  if self.ndim() != other.ndim() {
    raise SparseArrayError("SparseArray bind requires arrays of equal rank")
  }
  if along < 0 || along >= self.ndim() {
    raise SparseArrayError("SparseArray bind dimension is out of bounds")
  }
  let output_dimensions = self.dimensions.copy()
  let mut dimension = 0
  while dimension < self.ndim() {
    if dimension != along &&
      self.dimensions[dimension] != other.dimensions[dimension] {
      raise SparseArrayError(
        "SparseArray bind dimensions differ outside the binding margin",
      )
    }
    dimension = dimension + 1
  }
  output_dimensions[along] = self.dimensions[along] + other.dimensions[along]
  let output_entries = self.nonzero_entries()
  let mut index = 0
  while index < other.entries.length() {
    let coordinates = other.entries[index].coordinates.copy()
    coordinates[along] = coordinates[along] + self.dimensions[along]
    output_entries.push(SparseArrayEntry::{
      coordinates,
      value: other.entries[index].value,
    })
    index = index + 1
  }
  sparse_array_build(output_dimensions, output_entries)
}

///|
fn sparse_array_require_same_dimensions(
  left : SparseArray,
  right : SparseArray,
) -> Unit raise SparseArrayError {
  if left.dimensions != right.dimensions {
    raise SparseArrayError(
      "SparseArray arithmetic requires identical dimensions",
    )
  }
}

///|
fn sparse_array_add_scaled(
  left : SparseArray,
  right : SparseArray,
  right_scale : Double,
) -> SparseArray raise SparseArrayError {
  sparse_array_require_same_dimensions(left, right)
  let output : Array[SparseArrayEntry] = []
  let mut left_index = 0
  let mut right_index = 0
  while left_index < left.entries.length() ||
        right_index < right.entries.length() {
    if left_index >= left.entries.length() {
      output.push(SparseArrayEntry::{
        coordinates: right.entries[right_index].coordinates.copy(),
        value: right.entries[right_index].value * right_scale,
      })
      right_index = right_index + 1
    } else if right_index >= right.entries.length() {
      output.push(SparseArrayEntry::{
        coordinates: left.entries[left_index].coordinates.copy(),
        value: left.entries[left_index].value,
      })
      left_index = left_index + 1
    } else {
      let comparison = sparse_array_compare_coordinates(
        left.entries[left_index].coordinates,
        right.entries[right_index].coordinates,
      )
      if comparison < 0 {
        output.push(SparseArrayEntry::{
          coordinates: left.entries[left_index].coordinates.copy(),
          value: left.entries[left_index].value,
        })
        left_index = left_index + 1
      } else if comparison > 0 {
        output.push(SparseArrayEntry::{
          coordinates: right.entries[right_index].coordinates.copy(),
          value: right.entries[right_index].value * right_scale,
        })
        right_index = right_index + 1
      } else {
        output.push(SparseArrayEntry::{
          coordinates: left.entries[left_index].coordinates.copy(),
          value: left.entries[left_index].value +
          right.entries[right_index].value * right_scale,
        })
        left_index = left_index + 1
        right_index = right_index + 1
      }
    }
  }
  sparse_array_build(left.dimensions, output)
}

///|
pub fn SparseArray::add(
  self : SparseArray,
  other : SparseArray,
) -> SparseArray raise SparseArrayError {
  sparse_array_add_scaled(self, other, 1.0)
}

///|
pub fn SparseArray::subtract(
  self : SparseArray,
  other : SparseArray,
) -> SparseArray raise SparseArrayError {
  sparse_array_add_scaled(self, other, -1.0)
}

///|
pub fn SparseArray::hadamard(
  self : SparseArray,
  other : SparseArray,
) -> SparseArray raise SparseArrayError {
  sparse_array_require_same_dimensions(self, other)
  let output : Array[SparseArrayEntry] = []
  let mut left_index = 0
  let mut right_index = 0
  while left_index < self.entries.length() &&
        right_index < other.entries.length() {
    let comparison = sparse_array_compare_coordinates(
      self.entries[left_index].coordinates,
      other.entries[right_index].coordinates,
    )
    if comparison < 0 {
      left_index = left_index + 1
    } else if comparison > 0 {
      right_index = right_index + 1
    } else {
      output.push(SparseArrayEntry::{
        coordinates: self.entries[left_index].coordinates.copy(),
        value: self.entries[left_index].value * other.entries[right_index].value,
      })
      left_index = left_index + 1
      right_index = right_index + 1
    }
  }
  sparse_array_build(self.dimensions, output)
}

///|
pub fn SparseArray::scale(
  self : SparseArray,
  factor : Double,
) -> SparseArray raise SparseArrayError {
  let output : Array[SparseArrayEntry] = []
  let mut index = 0
  while index < self.entries.length() {
    output.push(SparseArrayEntry::{
      coordinates: self.entries[index].coordinates.copy(),
      value: self.entries[index].value * factor,
    })
    index = index + 1
  }
  sparse_array_build(self.dimensions, output)
}

///|
/// Apply a zero-preserving operation to all stored nonzero values.
pub fn SparseArray::map_nonzero(
  self : SparseArray,
  transform : (Double) -> Double,
) -> SparseArray raise SparseArrayError {
  let output : Array[SparseArrayEntry] = []
  let mut index = 0
  while index < self.entries.length() {
    output.push(SparseArrayEntry::{
      coordinates: self.entries[index].coordinates.copy(),
      value: transform(self.entries[index].value),
    })
    index = index + 1
  }
  sparse_array_build(self.dimensions, output)
}

///|
pub fn SparseArray::sum(self : SparseArray) -> Double {
  let mut result = 0.0
  let mut index = 0
  while index < self.entries.length() {
    result = result + self.entries[index].value
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::mean(self : SparseArray) -> Double? raise SparseArrayError {
  let total = self.length()
  if total == 0 {
    None
  } else {
    Some(self.sum() / total.to_double())
  }
}

///|
pub fn SparseArray::minimum(
  self : SparseArray,
) -> Double? raise SparseArrayError {
  let total = self.length()
  if total == 0 {
    return None
  }
  let mut result = if self.entries.length() < total {
    0.0
  } else {
    self.entries[0].value
  }
  let mut index = 0
  while index < self.entries.length() {
    if self.entries[index].value < result {
      result = self.entries[index].value
    }
    index = index + 1
  }
  Some(result)
}

///|
pub fn SparseArray::maximum(
  self : SparseArray,
) -> Double? raise SparseArrayError {
  let total = self.length()
  if total == 0 {
    return None
  }
  let mut result = if self.entries.length() < total {
    0.0
  } else {
    self.entries[0].value
  }
  let mut index = 0
  while index < self.entries.length() {
    if self.entries[index].value > result {
      result = self.entries[index].value
    }
    index = index + 1
  }
  Some(result)
}

///|
fn sparse_array_require_matrix(
  array : SparseArray,
) -> Unit raise SparseArrayError {
  if array.ndim() != 2 {
    raise SparseArrayError(
      "SparseArray operation requires a two-dimensional array",
    )
  }
}

///|
pub fn SparseArray::row_sums(
  self : SparseArray,
) -> Array[Double] raise SparseArrayError {
  sparse_array_require_matrix(self)
  let result = Array::make(self.dimensions[0], 0.0)
  let mut index = 0
  while index < self.entries.length() {
    let row = self.entries[index].coordinates[0]
    result[row] = result[row] + self.entries[index].value
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::column_sums(
  self : SparseArray,
) -> Array[Double] raise SparseArrayError {
  sparse_array_require_matrix(self)
  let result = Array::make(self.dimensions[1], 0.0)
  let mut index = 0
  while index < self.entries.length() {
    let column = self.entries[index].coordinates[1]
    result[column] = result[column] + self.entries[index].value
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::row_means(
  self : SparseArray,
) -> Array[Double] raise SparseArrayError {
  sparse_array_require_matrix(self)
  if self.dimensions[1] == 0 {
    raise SparseArrayError(
      "Cannot compute row means for a matrix with zero columns",
    )
  }
  let result = self.row_sums()
  let mut row = 0
  while row < result.length() {
    result[row] = result[row] / self.dimensions[1].to_double()
    row = row + 1
  }
  result
}

///|
pub fn SparseArray::column_means(
  self : SparseArray,
) -> Array[Double] raise SparseArrayError {
  sparse_array_require_matrix(self)
  if self.dimensions[0] == 0 {
    raise SparseArrayError(
      "Cannot compute column means for a matrix with zero rows",
    )
  }
  let result = self.column_sums()
  let mut column = 0
  while column < result.length() {
    result[column] = result[column] / self.dimensions[0].to_double()
    column = column + 1
  }
  result
}

///|
pub fn SparseArray::row_nonzero_counts(
  self : SparseArray,
) -> Array[Int] raise SparseArrayError {
  sparse_array_require_matrix(self)
  let result = Array::make(self.dimensions[0], 0)
  let mut index = 0
  while index < self.entries.length() {
    let row = self.entries[index].coordinates[0]
    result[row] = result[row] + 1
    index = index + 1
  }
  result
}

///|
pub fn SparseArray::column_nonzero_counts(
  self : SparseArray,
) -> Array[Int] raise SparseArrayError {
  sparse_array_require_matrix(self)
  let result = Array::make(self.dimensions[1], 0)
  let mut index = 0
  while index < self.entries.length() {
    let column = self.entries[index].coordinates[1]
    result[column] = result[column] + 1
    index = index + 1
  }
  result
}

///|
/// Sparse matrix multiplication. The result is dense, matching the
/// Bioconductor SparseMatrix multiplication contract.
pub fn SparseArray::matmul(
  self : SparseArray,
  other : SparseArray,
) -> Array[Array[Double]] raise SparseArrayError {
  sparse_array_require_matrix(self)
  sparse_array_require_matrix(other)
  if self.dimensions[1] != other.dimensions[0] {
    raise SparseArrayError(
      "SparseArray matrix multiplication dimensions are incompatible",
    )
  }
  let result : Array[Array[Double]] = []
  let mut row = 0
  while row < self.dimensions[0] {
    result.push(Array::make(other.dimensions[1], 0.0))
    row = row + 1
  }
  let mut left_index = 0
  while left_index < self.entries.length() {
    let left_coordinates = self.entries[left_index].coordinates
    let mut right_index = 0
    while right_index < other.entries.length() {
      let right_coordinates = other.entries[right_index].coordinates
      if left_coordinates[1] == right_coordinates[0] {
        let output_row = left_coordinates[0]
        let output_column = right_coordinates[1]
        result[output_row][output_column] = result[output_row][output_column] +
          self.entries[left_index].value * other.entries[right_index].value
      }
      right_index = right_index + 1
    }
    left_index = left_index + 1
  }
  result
}

///|
pub fn SparseArray::crossprod(
  self : SparseArray,
) -> Array[Array[Double]] raise SparseArrayError {
  self.transpose().matmul(self)
}

///|
pub fn SparseArray::tcrossprod(
  self : SparseArray,
) -> Array[Array[Double]] raise SparseArrayError {
  self.matmul(self.transpose())
}

///|
pub fn SparseArray::summary(
  self : SparseArray,
) -> String raise SparseArrayError {
  let mut dimensions = ""
  let mut index = 0
  while index < self.dimensions.length() {
    if index > 0 {
      dimensions = dimensions + " x "
    }
    dimensions = dimensions + self.dimensions[index].to_string()
    index = index + 1
  }
  "SparseArray(" +
  dimensions +
  ", nnzero=" +
  self.nnzero().to_string() +
  ", density=" +
  self.density().to_string() +
  ")"
}

///|
pub fn sparse_array_sample() -> SparseArray raise SparseArrayError {
  SparseArray::from_coo(
    [4, 5, 2],
    [[0, 0, 0], [1, 0, 0], [3, 2, 0], [0, 4, 1], [2, 1, 1]],
    [10.0, 20.0, 30.0, 40.0, 50.0],
  )
}
