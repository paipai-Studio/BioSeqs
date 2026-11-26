/// Biotite `structure.density` — atomistic density as a 3D histogram.
///
/// Faithful port of `biotite.structure.density.density`: the
/// coordinates are binned into an even grid (spacing `delta`) spanning
/// the coordinate range; optional weights are accumulated per bin and
/// `density=true` normalizes the histogram into a probability density.
/// Returns the histogram together with the bin edges of all three axes.

///|
/// A 3D density histogram with the bin edges of each axis.
pub struct BiotiteDensity {
  // hist[x][y][z]
  hist : Array[Array[Array[Double]]]
  edges_x : Array[Double]
  edges_y : Array[Double]
  edges_z : Array[Double]
  delta : Double
}

///|
/// Compute the 3D density histogram of the given coordinates.
///
/// - `delta`: grid spacing in all three dimensions (default 1.0 A).
/// - `weights`: optional per-atom weights accumulated per bin
///   (default: each atom counts 1).
/// - `density`: when true, the histogram is normalized into a
///   probability density (count divided by total weight and bin
///   volume delta^3), mirroring `numpy.histogramdd(density=True)`.
pub fn biotite_density(
  coords : Array[Vector3],
  delta~ : Double = 1.0,
  weights~ : Array[Double]? = None,
  density~ : Bool = false,
) -> BiotiteDensity {
  let n = coords.length()
  // Grid spanned by the coordinate extremes (as histogramdd with
  // default bins based on arange(min, max + delta, delta)).
  let mut min_x = 0.0
  let mut min_y = 0.0
  let mut min_z = 0.0
  let mut max_x = 0.0
  let mut max_y = 0.0
  let mut max_z = 0.0
  if n > 0 {
    let c0 = coords[0]
    min_x = c0.x
    min_y = c0.y
    min_z = c0.z
    max_x = c0.x
    max_y = c0.y
    max_z = c0.z
    for i = 1; i < n; i = i + 1 {
      let c = coords[i]
      if c.x < min_x {
        min_x = c.x
      }
      if c.y < min_y {
        min_y = c.y
      }
      if c.z < min_z {
        min_z = c.z
      }
      if c.x > max_x {
        max_x = c.x
      }
      if c.y > max_y {
        max_y = c.y
      }
      if c.z > max_z {
        max_z = c.z
      }
    }
  }
  let (nx, ny, nz) = biotite_density_bin_count(
    max_x - min_x,
    max_y - min_y,
    max_z - min_z,
    delta,
  )
  let edges_x = biotite_density_edges(min_x, delta, nx)
  let edges_y = biotite_density_edges(min_y, delta, ny)
  let edges_z = biotite_density_edges(min_z, delta, nz)
  let hist = biotite_density_grid(nx, ny, nz)
  for i = 0; i < n; i = i + 1 {
    let w = match weights {
      Some(ws) => if i < ws.length() { ws[i] } else { 0.0 }
      None => 1.0
    }
    let c = coords[i]
    let ix = biotite_density_bin_index(c.x, min_x, delta, nx)
    let iy = biotite_density_bin_index(c.y, min_y, delta, ny)
    let iz = biotite_density_bin_index(c.z, min_z, delta, nz)
    hist[ix][iy][iz] = hist[ix][iy][iz] + w
  }
  let out = BiotiteDensity::{ hist, edges_x, edges_y, edges_z, delta, }
  if density {
    biotite_density_normalize(out)
  }
  out
}

///|
/// Number of bins per axis: ceil(span / delta), mirroring the number
/// of arange(min, max + delta, delta) edges minus one.
fn biotite_density_bin_count(
  span_x : Double,
  span_y : Double,
  span_z : Double,
  delta : Double,
) -> (Int32, Int32, Int32) {
  if delta <= 0.0 {
    return (1, 1, 1)
  }
  let count = fn(span : Double) -> Int32 {
    if span <= 0.0 {
      1
    } else {
      let ratio = span / delta
      let whole = ratio.to_int()
      if ratio > whole.to_double() {
        whole + 1
      } else {
        whole
      }
    }
  }
  (count(span_x), count(span_y), count(span_z))
}

///|
/// Bin edges along one axis: `n_bins + 1` values from `min` stepping
/// by `delta`.
fn biotite_density_edges(
  min : Double,
  delta : Double,
  n_bins : Int32,
) -> Array[Double] {
  let edges : Array[Double] = Array(capacity=n_bins + 1)
  for i = 0; i <= n_bins; i = i + 1 {
    edges.push(min + i.to_double() * delta)
  }
  edges
}

///|
/// All-zero grid; rows allocated independently to avoid shared-array
/// aliasing.
fn biotite_density_grid(
  nx : Int32,
  ny : Int32,
  nz : Int32,
) -> Array[Array[Array[Double]]] {
  let grid : Array[Array[Array[Double]]] = Array(capacity=nx)
  for _x = 0; _x < nx; _x = _x + 1 {
    let plane : Array[Array[Double]] = Array(capacity=ny)
    for _y = 0; _y < ny; _y = _y + 1 {
      plane.push(Array::make(nz, 0.0))
    }
    grid.push(plane)
  }
  grid
}

///|
/// Bin index of a coordinate (clamped into the grid; values at the
/// upper edge fall into the last bin like histogramdd's right closure).
fn biotite_density_bin_index(
  value : Double,
  min : Double,
  delta : Double,
  n_bins : Int32,
) -> Int32 {
  if delta <= 0.0 {
    return 0
  }
  let mut idx = ((value - min) / delta).to_int()
  if idx < 0 {
    idx = 0
  }
  if idx >= n_bins {
    idx = n_bins - 1
  }
  idx
}

///|
/// Normalize the histogram into a probability density:
/// hist / (total weight * bin volume).
fn biotite_density_normalize(d : BiotiteDensity) -> Unit {
  let mut total = 0.0
  for plane in d.hist {
    for row in plane {
      for v in row {
        total = total + v
      }
    }
  }
  let volume = d.delta * d.delta * d.delta
  if total > 0.0 && volume > 0.0 {
    let scale = 1.0 / (total * volume)
    for plane in d.hist {
      for row in plane {
        for k = 0; k < row.length(); k = k + 1 {
          row[k] = row[k] * scale
        }
      }
    }
  }
}

///|
/// Grid extents: (nx, ny, nz).
pub fn BiotiteDensity::shape(self : BiotiteDensity) -> (Int32, Int32, Int32) {
  let nx = self.hist.length()
  let ny = if nx > 0 { self.hist[0].length() } else { 0 }
  let nz = if ny > 0 && self.hist[0].length() > 0 {
    self.hist[0][0].length()
  } else {
    0
  }
  (nx, ny, nz)
}

///|
/// Total density over the whole grid.
pub fn BiotiteDensity::total(self : BiotiteDensity) -> Double {
  let mut sum = 0.0
  for plane in self.hist {
    for row in plane {
      for v in row {
        sum = sum + v
      }
    }
  }
  sum
}

///|
/// Maximum density value in the grid.
pub fn BiotiteDensity::max(self : BiotiteDensity) -> Double {
  let mut best = -1.7976931348623157e308
  for plane in self.hist {
    for row in plane {
      for v in row {
        if v > best {
          best = v
        }
      }
    }
  }
  best
}

///|
/// Bin edges of the x axis.
pub fn BiotiteDensity::get_edges_x(self : BiotiteDensity) -> Array[Double] {
  self.edges_x
}

///|
/// Bin edges of the y axis.
pub fn BiotiteDensity::get_edges_y(self : BiotiteDensity) -> Array[Double] {
  self.edges_y
}

///|
/// Bin edges of the z axis.
pub fn BiotiteDensity::get_edges_z(self : BiotiteDensity) -> Array[Double] {
  self.edges_z
}

///|
pub fn _biotite_structure_density_dummy() -> Unit {
  let d = biotite_density([], delta=1.0)
  let _ = d.shape()
  let _ = d.total()
  let _ = d.max()
  let _ = d.get_edges_x()
  let _ = d.get_edges_y()
  let _ = d.get_edges_z()
  ignore(d.delta)
}
