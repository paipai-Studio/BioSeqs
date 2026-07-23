///|
/// HilbertCurve - Hilbert curve coordinate mapping for genomic visualization.

///|
pub struct HilbertCurve {
  levels : Int
  dimension : Int
  max_coordinate : Int
}

///|
pub fn HilbertCurve::new(levels : Int, dimension : Int) -> HilbertCurve {
  let max_coord = (1 << levels) - 1
  HilbertCurve::{ levels, dimension, max_coordinate: max_coord }
}

///|
pub fn hilbert_encode(hc : HilbertCurve, coordinates : Array[Int]) -> Int {
  let mut result = 0
  let mut s = 1 << (hc.levels - 1)
  
  let x : Array[Int] = coordinates.copy()
  let n = hc.dimension
  
  while s > 0 {
    let mut region = 0
    let mut i = 0
    while i < n {
      if bit_is_set(x[i], s) {
        region = region | (1 << i)
      }
      i = i + 1
    }
    
    let shift_amount = n * (hc.levels - log2_int(s))
    result = result | (region << shift_amount)
    
    hilbert_inverse(hc, x, s)
    
    s = s >> 1
  }
  
  result
}

///|
fn bit_is_set(x : Int, mask : Int) -> Bool {
  (x & mask) != 0
}

///|
fn log2_int(x : Int) -> Int {
  let mut result = 0
  let mut value = x
  while value > 1 {
    value = value >> 1
    result = result + 1
  }
  result
}

///|
fn hilbert_inverse(hc : HilbertCurve, x : Array[Int], s : Int) -> Unit {
  let n = hc.dimension
  
  let shift_amount = log2_int(s)
  let bits = (x[0] >> shift_amount) & 1
  let mut i = 1
  while i < n {
    let bit_i = (x[i] >> shift_amount) & 1
    x[0] = x[0] ^ (bit_i << shift_amount)
    i = i + 1
  }
  
  let mut j = 1
  while j < n {
    let mut k = j + 1
    while k < n {
      let bit_j = (x[j] >> shift_amount) & 1
      let bit_k = (x[k] >> shift_amount) & 1
      x[j] = x[j] ^ ((bit_j ^ bit_k) << shift_amount)
      k = k + 1
    }
    j = j + 1
  }
}

///|
pub fn hilbert_decode(hc : HilbertCurve, index : Int) -> Array[Int] {
  let x : Array[Int] = Array::new()
  let mut i = 0
  while i < hc.dimension {
    x.push(0)
    i = i + 1
  }
  
  let mut s = 1
  let mut t = index
  
  let n = hc.dimension
  
  while s <= hc.max_coordinate {
    let region = t & ((1 << n) - 1)
    
    hilbert_inverse(hc, x, s)
    
    let mut i = 0
    while i < n {
      if bit_is_set(region, 1 << i) {
        x[i] = x[i] | s
      }
      i = i + 1
    }
    
    t = t >> n
    s = s << 1
  }
  
  x
}

///|
pub fn hilbert_point_to_segment(hc : HilbertCurve, start : Int, end : Int) -> Array[(Int, Int)] {
  let segments : Array[(Int, Int)] = Array::new()
  
  let mut i = start
  while i <= end {
    let coords = hilbert_decode(hc, i)
    if coords.length() >= 2 {
      segments.push((coords[0], coords[1]))
    }
    i = i + 1
  }
  
  segments
}

///|
pub fn hilbert_distance(hc : HilbertCurve, coord1 : Array[Int], coord2 : Array[Int]) -> Int {
  let idx1 = hilbert_encode(hc, coord1)
  let idx2 = hilbert_encode(hc, coord2)
  
  let diff = if idx1 > idx2 { idx1 - idx2 } else { idx2 - idx1 }
  diff
}

///|
pub fn hilbert_linearize_genome(genome_length : Int, levels : Int) -> Array[(Int, Int)] {
  let hc = HilbertCurve::new(levels, 2)
  let max_coord = (1 << levels) - 1
  
  let segments : Array[(Int, Int)] = Array::new()
  
  let mut i = 0
  while i < genome_length {
    let x = i % (max_coord + 1)
    let y = i / (max_coord + 1)
    
    let idx = hilbert_encode(hc, [x, y])
    segments.push((i, idx))
    i = i + 1
  }
  
  segments
}

///|
pub fn hilbert_map_to_grid(hc : HilbertCurve, values : Array[Double]) -> Array[Array[Double]] {
  let size = hc.max_coordinate + 1
  let grid : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < size {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < size {
      row.push(0.0)
      j = j + 1
    }
    grid.push(row)
    i = i + 1
  }
  
  let mut idx = 0
  while idx < values.length() {
    let coords = hilbert_decode(hc, idx)
    if coords.length() >= 2 {
      let x = coords[0]
      let y = coords[1]
      if x < size && y < size {
        grid[x][y] = grid[x][y] + values[idx]
      }
    }
    idx = idx + 1
  }
  
  grid
}

///|
pub fn create_example_hilbert_curve() -> HilbertCurve {
  HilbertCurve::new(4, 2)
}
