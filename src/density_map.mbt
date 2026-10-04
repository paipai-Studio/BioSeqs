///|
/// X-ray electron density maps: CCP4/MRC binary format, X-PLOR text format,
/// trilinear interpolation, contour extraction, simple electron-density
/// calculation from atomic coordinates (Gaussian model), and map-to-chain
/// conversion for visualisation.
///
/// ProDy densityMap subset: parseCCP4, parseXPLOR, sample, stats, contour,
/// add, subtract, map arithmetic.

///|
pub suberror DensityMapError {
  DensityMapError(String)
}

// ---------------------------------------------------------------------------
// Binary helpers (little-endian, same pattern as dcd.mbt)
// ---------------------------------------------------------------------------

///|
fn dm_read_i32_le(data : String, pos : Int) -> Int {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
}

///|
fn dm_read_f32_le(data : String, pos : Int) -> Double {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  let sign = if b3 >= 128 { -1.0 } else { 1.0 }
  let exponent = (b3 & 0x7F) * 2 + (b2 >> 7)
  let mantissa = (b2 & 0x7F) * 65536 + b1 * 256 + b0
  if exponent == 255 {
    0.0 / 0.0
  } else if exponent == 0 {
    sign * mantissa.to_double() * @math.pow(2.0, -149.0)
  } else {
    sign * (1.0 + mantissa.to_double() / 8388608.0) * @math.pow(
      2.0,
      (exponent - 127).to_double(),
    )
  }
}

///|
fn dm_write_i32_le(sb : StringBuilder, v : Int) -> Unit {
  sb.write_char(((v & 0xFF).to_uint16()).unsafe_to_char())
  sb.write_char((((v >> 8) & 0xFF).to_uint16()).unsafe_to_char())
  sb.write_char((((v >> 16) & 0xFF).to_uint16()).unsafe_to_char())
  sb.write_char((((v >> 24) & 0xFF).to_uint16()).unsafe_to_char())
}

///|
fn dm_write_f32_le(sb : StringBuilder, v : Double) -> Unit {
  if v == 0.0 {
    sb.write_char('\x00')
    sb.write_char('\x00')
    sb.write_char('\x00')
    sb.write_char('\x00')
    return
  }
  let sign : Int = if v < 0.0 { 1 } else { 0 }
  let absv = if v < 0.0 { -v } else { v }
  let log2v = @math.ln(absv) / @math.ln(2.0)
  let exponent = floor_d(log2v).to_int() + 127
  if exponent <= 0 {
    sb.write_char('\x00')
    sb.write_char('\x00')
    sb.write_char('\x00')
    sb.write_char('\x00')
    return
  }
  if exponent >= 255 {
    sb.write_char('\x00')
    sb.write_char('\x00')
    sb.write_char('\x80')
    sb.write_char((if sign == 1 { 0xFF } else { 0x7F }).to_uint16().unsafe_to_char())
    return
  }
  let pow2 = @math.pow(2.0, (exponent - 127).to_double())
  let mantissa_f = absv / pow2 - 1.0
  let mantissa = floor_d(mantissa_f * 8388608.0).to_int() & 0x7FFFFF
  let b3 = (sign << 7) | (exponent >> 1)
  let b2 = ((exponent & 1) << 7) | (mantissa >> 16)
  let b1 = (mantissa >> 8) & 0xFF
  let b0 = mantissa & 0xFF
  sb.write_char((b0.to_uint16()).unsafe_to_char())
  sb.write_char((b1.to_uint16()).unsafe_to_char())
  sb.write_char((b2.to_uint16()).unsafe_to_char())
  sb.write_char((b3.to_uint16()).unsafe_to_char())
}

///|
fn floor_d(x : Double) -> Double {
  let i = x.to_int()
  if x < 0.0 && i.to_double() != x {
    (i - 1).to_double()
  } else {
    i.to_double()
  }
}

///|
fn cos_d(deg : Double) -> Double {
  @math.cos(deg * 3.141592653589793 / 180.0)
}

///|
fn sin_d(deg : Double) -> Double {
  @math.sin(deg * 3.141592653589793 / 180.0)
}

// ---------------------------------------------------------------------------
// DensityMap struct
// ---------------------------------------------------------------------------

///|
pub(all) struct DensityMap {
  title : String
  nx : Int // grid columns (fastest)
  ny : Int // grid rows
  nz : Int // grid sections (slowest)
  mode : Int // data type (0=int8, 1=int16, 2=float32)
  nxstart : Int
  nystart : Int
  nzstart : Int
  mx : Int // unit-cell grid
  my : Int
  mz : Int
  xlen : Double // cell lengths (Å)
  ylen : Double
  zlen : Double
  alpha : Double // cell angles (degrees)
  beta : Double
  gamma : Double
  mapc : Int // axis mapping
  mapr : Int
  maps : Int
  dmin : Double // density stats (from header)
  dmax : Double
  dmean : Double
  ispg : Int // space-group number
  nsymbt : Int // symmetry-record bytes
  origin_x : Double // origin in Å (MRC2000)
  origin_y : Double
  origin_z : Double
  rms : Double
  data : Array[Double] // flattened density [nx*ny*nz], x fastest
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------

///|
pub fn DensityMap::get_title(self : DensityMap) -> String {
  self.title
}

///|
pub fn DensityMap::get_shape(self : DensityMap) -> (Int, Int, Int) {
  (self.nx, self.ny, self.nz)
}

///|
pub fn DensityMap::get_cell(
  self : DensityMap
) -> (Double, Double, Double, Double, Double, Double) {
  (self.xlen, self.ylen, self.zlen, self.alpha, self.beta, self.gamma)
}

///|
pub fn DensityMap::get_origin(self : DensityMap) -> (Double, Double, Double) {
  (self.origin_x, self.origin_y, self.origin_z)
}

///|
pub fn DensityMap::get_space_group(self : DensityMap) -> Int {
  self.ispg
}

///|
pub fn DensityMap::num_voxels(self : DensityMap) -> Int {
  self.nx * self.ny * self.nz
}

///|
pub fn DensityMap::get_data(self : DensityMap) -> Array[Double] {
  self.data
}

// ---------------------------------------------------------------------------
// Grid indexing
// ---------------------------------------------------------------------------

///|
/// Direct grid access: density at grid point (i, j, k).
/// Returns NaN for out-of-bounds indices.
pub fn DensityMap::get(self : DensityMap, i : Int, j : Int, k : Int) -> Double {
  if i < 0 || i >= self.nx || j < 0 || j >= self.ny || k < 0 || k >= self.nz {
    return 0.0 / 0.0
  }
  self.data[k * self.nx * self.ny + j * self.nx + i]
}

// ---------------------------------------------------------------------------
// Coordinate transforms
// ---------------------------------------------------------------------------

///|
/// Fractional → Cartesian for the map's unit cell.
pub fn DensityMap::frac_to_cart(
  self : DensityMap,
  u : Double,
  v : Double,
  w : Double,
) -> (Double, Double, Double) {
  let ca = cos_d(self.gamma)
  let cb = cos_d(self.beta)
  let cg = cos_d(self.alpha)
  let sg = sin_d(self.gamma)
  let vol_factor = (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
  let x = u * self.xlen + v * self.ylen * ca + w * self.zlen * cb
  let y = v * self.ylen * sg + w * self.zlen * (cg - cb * ca) / sg
  let z = w * self.zlen * vol_factor / sg
  (x, y, z)
}

///|
/// Cartesian → fractional for the map's unit cell.
pub fn DensityMap::cart_to_frac(
  self : DensityMap,
  x : Double,
  y : Double,
  z : Double,
) -> (Double, Double, Double) {
  let ca = cos_d(self.gamma)
  let cb = cos_d(self.beta)
  let cg = cos_d(self.alpha)
  let sg = sin_d(self.gamma)
  let vol_factor = (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
  let w = z * sg / (self.zlen * vol_factor)
  let v = (y - w * self.zlen * (cg - cb * ca) / sg) / (self.ylen * sg)
  let u = (x - v * self.ylen * ca - w * self.zlen * cb) / self.xlen
  (u, v, w)
}

// ---------------------------------------------------------------------------
// Trilinear interpolation
// ---------------------------------------------------------------------------

///|
/// Sample density at fractional coordinates (u, v, w) using trilinear
/// interpolation.  Fractional coords wrap around (periodic boundary).
pub fn DensityMap::sample_frac(
  self : DensityMap,
  u : Double,
  v : Double,
  w : Double,
) -> Double {
  let gx = u * self.nx.to_double() + self.nxstart.to_double()
  let gy = v * self.ny.to_double() + self.nystart.to_double()
  let gz = w * self.nz.to_double() + self.nzstart.to_double()
  self.sample_grid(gx, gy, gz)
}

///|
/// Sample density at grid coordinates (gx, gy, gz) with trilinear
/// interpolation and periodic boundary conditions.
pub fn DensityMap::sample_grid(
  self : DensityMap,
  gx : Double,
  gy : Double,
  gz : Double,
) -> Double {
  let i0 = floor_d(gx).to_int()
  let j0 = floor_d(gy).to_int()
  let k0 = floor_d(gz).to_int()
  let dx = gx - i0.to_double()
  let dy = gy - j0.to_double()
  let dz = gz - k0.to_double()
  let i0w = ((i0 % self.nx) + self.nx) % self.nx
  let j0w = ((j0 % self.ny) + self.ny) % self.ny
  let k0w = ((k0 % self.nz) + self.nz) % self.nz
  let i1w = (i0w + 1) % self.nx
  let j1w = (j0w + 1) % self.ny
  let k1w = (k0w + 1) % self.nz
  let c000 = self.data[k0w * self.nx * self.ny + j0w * self.nx + i0w]
  let c100 = self.data[k0w * self.nx * self.ny + j0w * self.nx + i1w]
  let c010 = self.data[k0w * self.nx * self.ny + j1w * self.nx + i0w]
  let c110 = self.data[k0w * self.nx * self.ny + j1w * self.nx + i1w]
  let c001 = self.data[k1w * self.nx * self.ny + j0w * self.nx + i0w]
  let c101 = self.data[k1w * self.nx * self.ny + j0w * self.nx + i1w]
  let c011 = self.data[k1w * self.nx * self.ny + j1w * self.nx + i0w]
  let c111 = self.data[k1w * self.nx * self.ny + j1w * self.nx + i1w]
  let c00 = c000 * (1.0 - dx) + c100 * dx
  let c10 = c010 * (1.0 - dx) + c110 * dx
  let c01 = c001 * (1.0 - dx) + c101 * dx
  let c11 = c011 * (1.0 - dx) + c111 * dx
  let c0 = c00 * (1.0 - dy) + c10 * dy
  let c1 = c01 * (1.0 - dy) + c11 * dy
  c0 * (1.0 - dz) + c1 * dz
}

///|
/// Sample density at Cartesian coordinates (x, y, z) in Å.
pub fn DensityMap::sample_cart(
  self : DensityMap,
  x : Double,
  y : Double,
  z : Double,
) -> Double {
  let (u, v, w) = self.cart_to_frac(x, y, z)
  self.sample_frac(u, v, w)
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

///|
/// Recompute density statistics: (min, max, mean, std, rms).
pub fn DensityMap::stats(
  self : DensityMap
) -> (Double, Double, Double, Double, Double) {
  let n = self.data.length()
  if n == 0 {
    return (0.0, 0.0, 0.0, 0.0, 0.0)
  }
  let mut mn = self.data[0]
  let mut mx = self.data[0]
  let mut sum = 0.0
  let mut sum2 = 0.0
  for i = 0; i < n; i = i + 1 {
    let v = self.data[i]
    if v < mn {
      mn = v
    }
    if v > mx {
      mx = v
    }
    sum = sum + v
    sum2 = sum2 + v * v
  }
  let mean = sum / n.to_double()
  let variance = sum2 / n.to_double() - mean * mean
  let std = if variance > 0.0 { variance.sqrt() } else { 0.0 }
  let rms = (sum2 / n.to_double()).sqrt()
  (mn, mx, mean, std, rms)
}

// ---------------------------------------------------------------------------
// Contour extraction
// ---------------------------------------------------------------------------

///|
/// Count the number of grid points with density above `level`.
pub fn DensityMap::count_above(self : DensityMap, level : Double) -> Int {
  let mut count = 0
  for i = 0; i < self.data.length(); i = i + 1 {
    if self.data[i] > level {
      count = count + 1
    }
  }
  count
}

///|
/// Extract grid points above `level` as (i, j, k, density) tuples.
pub fn DensityMap::contour_points(
  self : DensityMap,
  level : Double,
) -> Array[(Int, Int, Int, Double)] {
  let result : Array[(Int, Int, Int, Double)] = Array(capacity=0)
  for k = 0; k < self.nz; k = k + 1 {
    for j = 0; j < self.ny; j = j + 1 {
      for i = 0; i < self.nx; i = i + 1 {
        let v = self.data[k * self.nx * self.ny + j * self.nx + i]
        if v > level {
          result.push((i, j, k, v))
        }
      }
    }
  }
  result
}

///|
/// Find local maxima (peaks) in the density map.
/// A peak is a grid point whose density is higher than all 6 face neighbours.
pub fn DensityMap::find_peaks(
  self : DensityMap,
  min_density : Double,
) -> Array[(Int, Int, Int, Double)] {
  let result : Array[(Int, Int, Int, Double)] = Array(capacity=0)
  for k = 0; k < self.nz; k = k + 1 {
    for j = 0; j < self.ny; j = j + 1 {
      for i = 0; i < self.nx; i = i + 1 {
        let v = self.data[k * self.nx * self.ny + j * self.nx + i]
        if v < min_density {
          continue
        }
        let ip = (i + 1) % self.nx
        let im = (i - 1 + self.nx) % self.nx
        let jp = (j + 1) % self.ny
        let jm = (j - 1 + self.ny) % self.ny
        let kp = (k + 1) % self.nz
        let km = (k - 1 + self.nz) % self.nz
        let is_peak = v > self.get(im, j, k) &&
          v > self.get(ip, j, k) &&
          v > self.get(i, jm, k) &&
          v > self.get(i, jp, k) &&
          v > self.get(i, j, km) &&
          v > self.get(i, j, kp)
        if is_peak {
          result.push((i, j, k, v))
        }
      }
    }
  }
  result
}

// ---------------------------------------------------------------------------
// Map arithmetic
// ---------------------------------------------------------------------------

///|
/// Add two maps element-wise.  Maps must have the same shape.
pub fn DensityMap::add(
  self : DensityMap,
  other : DensityMap,
) -> DensityMap raise DensityMapError {
  if self.nx != other.nx || self.ny != other.ny || self.nz != other.nz {
    raise DensityMapError("Map shapes do not match")
  }
  let n = self.data.length()
  let data : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    data.push(self.data[i] + other.data[i])
  }
  DensityMap::{
    title: self.title + " + " + other.title,
    nx: self.nx,
    ny: self.ny,
    nz: self.nz,
    mode: self.mode,
    nxstart: self.nxstart,
    nystart: self.nystart,
    nzstart: self.nzstart,
    mx: self.mx,
    my: self.my,
    mz: self.mz,
    xlen: self.xlen,
    ylen: self.ylen,
    zlen: self.zlen,
    alpha: self.alpha,
    beta: self.beta,
    gamma: self.gamma,
    mapc: self.mapc,
    mapr: self.mapr,
    maps: self.maps,
    dmin: 0.0,
    dmax: 0.0,
    dmean: 0.0,
    ispg: self.ispg,
    nsymbt: 0,
    origin_x: self.origin_x,
    origin_y: self.origin_y,
    origin_z: self.origin_z,
    rms: 0.0,
    data,
  }
}

///|
/// Subtract `other` from `self` element-wise.
pub fn DensityMap::subtract(
  self : DensityMap,
  other : DensityMap,
) -> DensityMap raise DensityMapError {
  if self.nx != other.nx || self.ny != other.ny || self.nz != other.nz {
    raise DensityMapError("Map shapes do not match")
  }
  let n = self.data.length()
  let data : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    data.push(self.data[i] - other.data[i])
  }
  DensityMap::{
    title: self.title + " - " + other.title,
    nx: self.nx,
    ny: self.ny,
    nz: self.nz,
    mode: self.mode,
    nxstart: self.nxstart,
    nystart: self.nystart,
    nzstart: self.nzstart,
    mx: self.mx,
    my: self.my,
    mz: self.mz,
    xlen: self.xlen,
    ylen: self.ylen,
    zlen: self.zlen,
    alpha: self.alpha,
    beta: self.beta,
    gamma: self.gamma,
    mapc: self.mapc,
    mapr: self.mapr,
    maps: self.maps,
    dmin: 0.0,
    dmax: 0.0,
    dmean: 0.0,
    ispg: self.ispg,
    nsymbt: 0,
    origin_x: self.origin_x,
    origin_y: self.origin_y,
    origin_z: self.origin_z,
    rms: 0.0,
    data,
  }
}

///|
/// Scale all density values by `factor`.
pub fn DensityMap::scale(self : DensityMap, factor : Double) -> DensityMap {
  let n = self.data.length()
  let data : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    data.push(self.data[i] * factor)
  }
  DensityMap::{
    title: self.title,
    nx: self.nx,
    ny: self.ny,
    nz: self.nz,
    mode: self.mode,
    nxstart: self.nxstart,
    nystart: self.nystart,
    nzstart: self.nzstart,
    mx: self.mx,
    my: self.my,
    mz: self.mz,
    xlen: self.xlen,
    ylen: self.ylen,
    zlen: self.zlen,
    alpha: self.alpha,
    beta: self.beta,
    gamma: self.gamma,
    mapc: self.mapc,
    mapr: self.mapr,
    maps: self.maps,
    dmin: self.dmin * factor,
    dmax: self.dmax * factor,
    dmean: self.dmean * factor,
    ispg: self.ispg,
    nsymbt: self.nsymbt,
    origin_x: self.origin_x,
    origin_y: self.origin_y,
    origin_z: self.origin_z,
    rms: self.rms * factor,
    data,
  }
}

// ---------------------------------------------------------------------------
// CCP4/MRC binary format parser
// ---------------------------------------------------------------------------

///|
/// Parse a CCP4/MRC binary map (mode 2 = float32, the most common).
/// Returns a DensityMap.
pub fn parse_ccp4(content : String) -> DensityMap raise DensityMapError {
  if content.length() < 1024 {
    raise DensityMapError("CCP4 map too short: header is 1024 bytes")
  }
  // Header: 256 int32/float32 words
  let nx = dm_read_i32_le(content, 0)
  let ny = dm_read_i32_le(content, 4)
  let nz = dm_read_i32_le(content, 8)
  let mode = dm_read_i32_le(content, 12)
  let nxstart = dm_read_i32_le(content, 16)
  let nystart = dm_read_i32_le(content, 20)
  let nzstart = dm_read_i32_le(content, 24)
  let mx = dm_read_i32_le(content, 28)
  let my = dm_read_i32_le(content, 32)
  let mz = dm_read_i32_le(content, 36)
  let xlen = dm_read_f32_le(content, 40)
  let ylen = dm_read_f32_le(content, 44)
  let zlen = dm_read_f32_le(content, 48)
  let alpha = dm_read_f32_le(content, 52)
  let beta = dm_read_f32_le(content, 56)
  let gamma = dm_read_f32_le(content, 60)
  let mapc = dm_read_i32_le(content, 64)
  let mapr = dm_read_i32_le(content, 68)
  let maps = dm_read_i32_le(content, 72)
  let dmin = dm_read_f32_le(content, 76)
  let dmax = dm_read_f32_le(content, 80)
  let dmean = dm_read_f32_le(content, 84)
  let ispg = dm_read_i32_le(content, 88)
  let nsymbt = dm_read_i32_le(content, 92)
  // Words 25-27 (offset 96..108): extra
  // Word 28-31 (offset 108..124): exttyp (4 chars)
  // Word 32 (offset 124): nversion
  // Words 33-35 (offset 128..140): origin in Å (MRC2000)
  let origin_x = dm_read_f32_le(content, 128)
  let origin_y = dm_read_f32_le(content, 132)
  let origin_z = dm_read_f32_le(content, 136)
  // Word 50 (offset 196): rms
  let rms = dm_read_f32_le(content, 196)
  // Word 51 (offset 200): nlabl
  let nlabl = dm_read_i32_le(content, 200)
  // Words 52-256 (offset 204..1024): labels (80 chars each)
  let mut title = ""
  if nlabl > 0 {
    let sb = StringBuilder()
    let label_start = 204
    let label_bytes = 80 * nlabl
    let mut li = 0
    while li < label_bytes && label_start + li < 1024 {
      let ch = content.unsafe_get(label_start + li)
      if ch.to_int() != 0 {
        sb.write_char(ch.unsafe_to_char())
      }
      li = li + 1
    }
    title = sb.to_string()
  }
  // Validate
  if nx <= 0 || ny <= 0 || nz <= 0 {
    raise DensityMapError("Invalid map dimensions: " + nx.to_string() + "x" + ny.to_string() + "x" + nz.to_string())
  }
  if mode != 2 && mode != 0 && mode != 1 {
    raise DensityMapError("Unsupported map mode: " + mode.to_string() + " (only 0/1/2 supported)")
  }
  // Data starts after header + symmetry
  let data_start = 1024 + nsymbt
  let nvox = nx * ny * nz
  let bytes_per = if mode == 0 { 1 } else if mode == 1 { 2 } else { 4 }
  let expected = data_start + nvox * bytes_per
  if content.length() < expected {
    raise DensityMapError("CCP4 map truncated: need " + expected.to_string() + " bytes, got " + content.length().to_string())
  }
  // Read density data (x fastest, y middle, z slowest)
  let data : Array[Double] = Array(capacity=nvox)
  let mut pos = data_start
  if mode == 2 {
    // float32
    for _i = 0; _i < nvox; _i = _i + 1 {
      data.push(dm_read_f32_le(content, pos))
      pos = pos + 4
    }
  } else if mode == 1 {
    // int16
    for _i = 0; _i < nvox; _i = _i + 1 {
      let b0 = content.unsafe_get(pos).to_int() & 0xFF
      let b1 = content.unsafe_get(pos + 1).to_int() & 0xFF
      let val = b0 | (b1 << 8)
      let signed = if val >= 32768 { val - 65536 } else { val }
      data.push(signed.to_double())
      pos = pos + 2
    }
  } else {
    // int8 (mode 0)
    for _i = 0; _i < nvox; _i = _i + 1 {
      let b0 = content.unsafe_get(pos).to_int() & 0xFF
      let signed = if b0 >= 128 { b0 - 256 } else { b0 }
      data.push(signed.to_double())
      pos = pos + 1
    }
  }
  // If origin is zero but nxstart/nystart/nzstart are non-zero,
  // compute origin from start + cell
  let (ox, oy, oz) = if origin_x == 0.0 && origin_y == 0.0 && origin_z == 0.0 {
    let dx = xlen / mx.to_double()
    let dy = ylen / my.to_double()
    let dz = zlen / mz.to_double()
    (nxstart.to_double() * dx, nystart.to_double() * dy, nzstart.to_double() * dz)
  } else {
    (origin_x, origin_y, origin_z)
  }
  DensityMap::{
    title,
    nx,
    ny,
    nz,
    mode,
    nxstart,
    nystart,
    nzstart,
    mx,
    my,
    mz,
    xlen,
    ylen,
    zlen,
    alpha,
    beta,
    gamma,
    mapc,
    mapr,
    maps,
    dmin,
    dmax,
    dmean,
    ispg,
    nsymbt,
    origin_x: ox,
    origin_y: oy,
    origin_z: oz,
    rms,
    data,
  }
}

// ---------------------------------------------------------------------------
// CCP4/MRC binary writer
// ---------------------------------------------------------------------------

///|
/// Write a DensityMap to CCP4/MRC binary format (mode 2 = float32).
/// Returns the binary content as a String.
pub fn write_ccp4(map : DensityMap) -> String {
  let sb = StringBuilder()
  // Compute stats
  let (dmin, dmax, dmean, _std, rms) = map.stats()
  // Header: 256 words
  dm_write_i32_le(sb, map.nx) // 1
  dm_write_i32_le(sb, map.ny) // 2
  dm_write_i32_le(sb, map.nz) // 3
  dm_write_i32_le(sb, 2) // 4: mode = float32
  dm_write_i32_le(sb, map.nxstart) // 5
  dm_write_i32_le(sb, map.nystart) // 6
  dm_write_i32_le(sb, map.nzstart) // 7
  dm_write_i32_le(sb, map.mx) // 8
  dm_write_i32_le(sb, map.my) // 9
  dm_write_i32_le(sb, map.mz) // 10
  dm_write_f32_le(sb, map.xlen) // 11
  dm_write_f32_le(sb, map.ylen) // 12
  dm_write_f32_le(sb, map.zlen) // 13
  dm_write_f32_le(sb, map.alpha) // 14
  dm_write_f32_le(sb, map.beta) // 15
  dm_write_f32_le(sb, map.gamma) // 16
  dm_write_i32_le(sb, map.mapc) // 17
  dm_write_i32_le(sb, map.mapr) // 18
  dm_write_i32_le(sb, map.maps) // 19
  dm_write_f32_le(sb, dmin) // 20
  dm_write_f32_le(sb, dmax) // 21
  dm_write_f32_le(sb, dmean) // 22
  dm_write_i32_le(sb, map.ispg) // 23
  dm_write_i32_le(sb, map.nsymbt) // 24
  // Words 25-27: extra (zeros)
  dm_write_i32_le(sb, 0)
  dm_write_i32_le(sb, 0)
  dm_write_i32_le(sb, 0)
  // Words 28-31: exttyp
  sb.write_char('M')
  sb.write_char('A')
  sb.write_char('P')
  sb.write_char(' ')
  // Word 32: nversion = 20140
  dm_write_i32_le(sb, 20140)
  // Words 33-35: origin in Å
  dm_write_f32_le(sb, map.origin_x)
  dm_write_f32_le(sb, map.origin_y)
  dm_write_f32_le(sb, map.origin_z)
  // Words 36-49: extra (zeros)
  for _i = 0; _i < 14; _i = _i + 1 {
    dm_write_i32_le(sb, 0)
  }
  // Word 50: rms
  dm_write_f32_le(sb, rms)
  // Word 51: nlabl
  dm_write_i32_le(sb, 1)
  // Words 52-256: label (200 bytes = 50 words)
  let label = map.title
  let mut written = 0
  let label_chars = label.to_array()
  for ch in label_chars {
    if written < 200 {
      sb.write_char(ch)
      written = written + 1
    }
  }
  while written < 200 {
    sb.write_char('\x00')
    written = written + 1
  }
  // Symmetry records (if any)
  for _i = 0; _i < map.nsymbt; _i = _i + 1 {
    sb.write_char('\x00')
  }
  // Density data (float32, x fastest, z slowest)
  for i = 0; i < map.data.length(); i = i + 1 {
    dm_write_f32_le(sb, map.data[i])
  }
  sb.to_string()
}

// ---------------------------------------------------------------------------
// CCP4 convenience constructor
// ---------------------------------------------------------------------------

///|
/// Create a simple orthorhombic DensityMap from raw data.
/// Cell is orthorhombic (alpha=beta=gamma=90), origin at (0,0,0).
pub fn DensityMap::simple(
  title~ : String = "",
  nx~ : Int,
  ny~ : Int,
  nz~ : Int,
  xlen~ : Double,
  ylen~ : Double,
  zlen~ : Double,
  data~ : Array[Double],
) -> DensityMap {
  DensityMap::{
    title,
    nx,
    ny,
    nz,
    mode: 2,
    nxstart: 0,
    nystart: 0,
    nzstart: 0,
    mx: nx,
    my: ny,
    mz: nz,
    xlen,
    ylen,
    zlen,
    alpha: 90.0,
    beta: 90.0,
    gamma: 90.0,
    mapc: 1,
    mapr: 2,
    maps: 3,
    dmin: 0.0,
    dmax: 0.0,
    dmean: 0.0,
    ispg: 1,
    nsymbt: 0,
    origin_x: 0.0,
    origin_y: 0.0,
    origin_z: 0.0,
    rms: 0.0,
    data,
  }
}

// ---------------------------------------------------------------------------
// X-PLOR text format parser
// ---------------------------------------------------------------------------

///|
/// Parse an X-PLOR format density map (text).
///
/// Format sketch:
/// ```
///   title line
///   remark lines ...
///   !
///   na nagrid_start nagrid_end  nb nbgrid_start nbgrid_end  nc ncgrid_start ncgrid_end
///   xlen ylen zlen alpha beta gamma
///   zvalue
///   -nsec+1
///   density values (10 per line, fortran format)
///   ...
/// ```
pub fn parse_xplor(content : String) -> DensityMap raise DensityMapError {
  let lines = content.split("\n").to_array()
  let mut idx = 0
  // Skip blank lines at start
  while idx < lines.length() && lines[idx].trim().length() == 0 {
    idx = idx + 1
  }
  if idx >= lines.length() {
    raise DensityMapError("Empty X-PLOR map")
  }
  // Line 1: title
  let title = lines[idx].trim().to_string()
  idx = idx + 1
  // Remark lines until '!'
  while idx < lines.length() && !lines[idx].trim().starts_with("!") {
    idx = idx + 1
  }
  if idx >= lines.length() {
    raise DensityMapError("X-PLOR: missing '!' separator")
  }
  idx = idx + 1 // skip '!'
  // Grid line: na nagrid_start nagrid_end  nb ...  nc ...
  if idx >= lines.length() {
    raise DensityMapError("X-PLOR: missing grid specification")
  }
  let grid_parts = lines[idx].trim().split(" ").to_array().filter(
    fn(s) { s.length() > 0 },
  ).map(fn(s) { s.to_string() })
  if grid_parts.length() < 9 {
    raise DensityMapError("X-PLOR: malformed grid line")
  }
  let na = dm_safe_parse_int(grid_parts[0])
  let nagrid_start = dm_safe_parse_int(grid_parts[1])
  let nagrid_end = dm_safe_parse_int(grid_parts[2])
  let nb = dm_safe_parse_int(grid_parts[3])
  let nbgrid_start = dm_safe_parse_int(grid_parts[4])
  let nbgrid_end = dm_safe_parse_int(grid_parts[5])
  let nc = dm_safe_parse_int(grid_parts[6])
  let ncgrid_start = dm_safe_parse_int(grid_parts[7])
  let ncgrid_end = dm_safe_parse_int(grid_parts[8])
  idx = idx + 1
  // Cell line: xlen ylen zlen alpha beta gamma
  if idx >= lines.length() {
    raise DensityMapError("X-PLOR: missing cell parameters")
  }
  let cell_parts = lines[idx].trim().split(" ").to_array().filter(
    fn(s) { s.length() > 0 },
  ).map(fn(s) { s.to_string() })
  if cell_parts.length() < 6 {
    raise DensityMapError("X-PLOR: malformed cell line")
  }
  let xlen = dm_safe_parse_double(cell_parts[0])
  let ylen = dm_safe_parse_double(cell_parts[1])
  let zlen = dm_safe_parse_double(cell_parts[2])
  let alpha = dm_safe_parse_double(cell_parts[3])
  let beta = dm_safe_parse_double(cell_parts[4])
  let gamma = dm_safe_parse_double(cell_parts[5])
  idx = idx + 1
  // Density sections
  let nx = nagrid_end - nagrid_start + 1
  let ny = nbgrid_end - nbgrid_start + 1
  let nz = ncgrid_end - ncgrid_start + 1
  let nvox = nx * ny * nz
  let data : Array[Double] = Array(capacity=nvox)
  // Each section: zvalue line, then "-ny" line, then ny rows of nx values
  let mut sec = 0
  while sec < nz && idx < lines.length() {
    // Skip zvalue line
    idx = idx + 1
    if idx >= lines.length() {
      break
    }
    // Skip "-ny" line
    idx = idx + 1
    // Read ny rows
    let mut row = 0
    while row < ny && idx < lines.length() {
      let parts = lines[idx].trim().split(" ").to_array().filter(
        fn(s) { s.length() > 0 },
      ).map(fn(s) { s.to_string() })
      for p in parts {
        if data.length() < nvox {
          data.push(dm_safe_parse_double(p))
        }
      }
      idx = idx + 1
      row = row + 1
    }
    sec = sec + 1
  }
  if data.length() != nvox {
    raise DensityMapError("X-PLOR: expected " + nvox.to_string() + " values, got " + data.length().to_string())
  }
  DensityMap::{
    title,
    nx,
    ny,
    nz,
    mode: 2,
    nxstart: nagrid_start,
    nystart: nbgrid_start,
    nzstart: ncgrid_start,
    mx: na,
    my: nb,
    mz: nc,
    xlen,
    ylen,
    zlen,
    alpha,
    beta,
    gamma,
    mapc: 1,
    mapr: 2,
    maps: 3,
    dmin: 0.0,
    dmax: 0.0,
    dmean: 0.0,
    ispg: 1,
    nsymbt: 0,
    origin_x: 0.0,
    origin_y: 0.0,
    origin_z: 0.0,
    rms: 0.0,
    data,
  }
}

///|
fn dm_safe_parse_int(s : String) -> Int {
  let v : Int = @string.parse_int(s.trim()) catch {
    _ => 0
  }
  v
}

///|
fn dm_safe_parse_double(s : String) -> Double {
  let v : Double = @string.parse_double(s.trim()) catch {
    _ => 0.0
  }
  v
}

// ---------------------------------------------------------------------------
// Electron density from atomic coordinates (Gaussian model)
// ---------------------------------------------------------------------------

///|
/// Compute electron density on a grid from a set of atoms using a
/// Gaussian model:  rho(r) = sum_i Z_i * exp(-B_i * |r - r_i|^2 / 4) * w
///
/// Parameters:
///   atoms    – array of Atom
///   nx, ny, nz – grid dimensions
///   cell     – (xlen, ylen, zlen, alpha, beta, gamma)
///   resolution – resolution in Å (determines Gaussian width)
///   b_factor   – additional temperature factor (default 0)
///
/// Returns a DensityMap.
pub fn density_from_atoms(
  atoms : Array[Atom],
  nx~ : Int,
  ny~ : Int,
  nz~ : Int,
  xlen~ : Double,
  ylen~ : Double,
  zlen~ : Double,
  alpha~ : Double = 90.0,
  beta~ : Double = 90.0,
  gamma~ : Double = 90.0,
  resolution~ : Double = 2.0,
  b_factor~ : Double = 0.0,
) -> DensityMap {
  let nvox = nx * ny * nz
  let data : Array[Double] = Array(capacity=nvox)
  for _i = 0; _i < nvox; _i = _i + 1 {
    data.push(0.0)
  }
  // Gaussian width: sigma = resolution / (2 * sqrt(2 * ln(2)))
  // B_total = 8 * pi^2 * sigma^2 + b_factor
  let sigma = resolution / (2.0 * (2.0 * @math.ln(2.0)).sqrt())
  let b_total = 8.0 * 3.141592653589793 * 3.141592653589793 * sigma * sigma + b_factor
  // For each atom, spread Gaussian onto nearby grid points
  // Cutoff: 3 * sigma in fractional coords
  let dx = xlen / nx.to_double()
  let dy = ylen / ny.to_double()
  let dz = zlen / nz.to_double()
  let cutoff = 3.0 * sigma
  let icut = (cutoff / dx).to_int() + 1
  let jcut = (cutoff / dy).to_int() + 1
  let kcut = (cutoff / dz).to_int() + 1
  // For orthorhombic cell (most common case)
  let inv_four_b = 1.0 / (4.0 * b_total / (8.0 * 3.141592653589793 * 3.141592653589793))
  // Actually: Gaussian = exp(-|r-ri|^2 / (2*sigma^2))
  // With B = 8*pi^2*sigma^2, sigma^2 = B/(8*pi^2)
  // 1/(2*sigma^2) = 4*pi^2/B
  let inv_2sigma2 = 4.0 * 3.141592653589793 * 3.141592653589793 / b_total
  for ai = 0; ai < atoms.length(); ai = ai + 1 {
    let atom = atoms[ai]
    // Atomic number (approximate from element)
    let z = atomic_number(atom.element)
    if z == 0 {
      continue
    }
    let ax = atom.coord.x
    let ay = atom.coord.y
    let az = atom.coord.z
    // Atom's fractional position
    let fx = ax / xlen
    let fy = ay / ylen
    let fz = az / zlen
    let gi0 = (fx * nx.to_double()).to_int()
    let gj0 = (fy * ny.to_double()).to_int()
    let gk0 = (fz * nz.to_double()).to_int()
    // Spread onto nearby grid points
    for dk = -kcut; dk <= kcut; dk = dk + 1 {
      for dj = -jcut; dj <= jcut; dj = dj + 1 {
        for di = -icut; di <= icut; di = di + 1 {
          let gi = ((gi0 + di) % nx + nx) % nx
          let gj = ((gj0 + dj) % ny + ny) % ny
          let gk = ((gk0 + dk) % nz + nz) % nz
          // Cartesian distance
          let rx = gi.to_double() * dx - ax
          let ry = gj.to_double() * dy - ay
          let rz = gk.to_double() * dz - az
          let r2 = rx * rx + ry * ry + rz * rz
          if r2 > cutoff * cutoff {
            continue
          }
          let gauss = @math.exp(-r2 * inv_2sigma2)
          data[gk * nx * ny + gj * nx + gi] = data[
            gk * nx * ny + gj * nx + gi
          ] + z.to_double() * gauss
        }
      }
    }
  }
  DensityMap::{
    title: "density_from_atoms",
    nx,
    ny,
    nz,
    mode: 2,
    nxstart: 0,
    nystart: 0,
    nzstart: 0,
    mx: nx,
    my: ny,
    mz: nz,
    xlen,
    ylen,
    zlen,
    alpha,
    beta,
    gamma,
    mapc: 1,
    mapr: 2,
    maps: 3,
    dmin: 0.0,
    dmax: 0.0,
    dmean: 0.0,
    ispg: 1,
    nsymbt: 0,
    origin_x: 0.0,
    origin_y: 0.0,
    origin_z: 0.0,
    rms: 0.0,
    data,
  }
}

///|
/// Approximate atomic number from element symbol.
fn atomic_number(element : String) -> Int {
  let e = element.to_lower()
  match e {
    "h" => 1
    "he" => 2
    "li" => 3
    "be" => 4
    "b" => 5
    "c" => 6
    "n" => 7
    "o" => 8
    "f" => 9
    "ne" => 10
    "na" => 11
    "mg" => 12
    "al" => 13
    "si" => 14
    "p" => 15
    "s" => 16
    "cl" => 17
    "ar" => 18
    "k" => 19
    "ca" => 20
    "fe" => 26
    "cu" => 29
    "zn" => 30
    "se" => 34
    "br" => 35
    "i" => 53
    _ => 0
  }
}

// ---------------------------------------------------------------------------
// Mask and utility operations
// ---------------------------------------------------------------------------

///|
/// Create a binary mask: 1.0 where density > level, 0.0 elsewhere.
pub fn DensityMap::mask(self : DensityMap, level : Double) -> DensityMap {
  let n = self.data.length()
  let data : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    if self.data[i] > level {
      data.push(1.0)
    } else {
      data.push(0.0)
    }
  }
  DensityMap::{
    title: "mask(" + self.title + ")",
    nx: self.nx,
    ny: self.ny,
    nz: self.nz,
    mode: self.mode,
    nxstart: self.nxstart,
    nystart: self.nystart,
    nzstart: self.nzstart,
    mx: self.mx,
    my: self.my,
    mz: self.mz,
    xlen: self.xlen,
    ylen: self.ylen,
    zlen: self.zlen,
    alpha: self.alpha,
    beta: self.beta,
    gamma: self.gamma,
    mapc: self.mapc,
    mapr: self.mapr,
    maps: self.maps,
    dmin: 0.0,
    dmax: 1.0,
    dmean: 0.0,
    ispg: self.ispg,
    nsymbt: self.nsymbt,
    origin_x: self.origin_x,
    origin_y: self.origin_y,
    origin_z: self.origin_z,
    rms: 0.0,
    data,
  }
}

///|
/// Correlation coefficient between two maps.
/// CC = sum((a - <a>)(b - <b>)) / sqrt(sum((a-<a>)^2) * sum((b-<b>)^2))
pub fn DensityMap::correlation(
  self : DensityMap,
  other : DensityMap,
) -> Double raise DensityMapError {
  if self.nx != other.nx || self.ny != other.ny || self.nz != other.nz {
    raise DensityMapError("Map shapes do not match")
  }
  let n = self.data.length()
  if n == 0 {
    return 0.0
  }
  let mut sa = 0.0
  let mut sb = 0.0
  for i = 0; i < n; i = i + 1 {
    sa = sa + self.data[i]
    sb = sb + other.data[i]
  }
  let ma = sa / n.to_double()
  let mb = sb / n.to_double()
  let mut cov = 0.0
  let mut va = 0.0
  let mut vb = 0.0
  for i = 0; i < n; i = i + 1 {
    let da = self.data[i] - ma
    let db = other.data[i] - mb
    cov = cov + da * db
    va = va + da * da
    vb = vb + db * db
  }
  let denom = (va * vb).sqrt()
  if denom == 0.0 {
    return 0.0
  }
  cov / denom
}

///|
/// Grid point to Cartesian coordinates (Å).
pub fn DensityMap::grid_to_cart(
  self : DensityMap,
  i : Int,
  j : Int,
  k : Int,
) -> (Double, Double, Double) {
  let dx = self.xlen / self.mx.to_double()
  let dy = self.ylen / self.my.to_double()
  let dz = self.zlen / self.mz.to_double()
  let x = (i + self.nxstart).to_double() * dx + self.origin_x
  let y = (j + self.nystart).to_double() * dy + self.origin_y
  let z = (k + self.nzstart).to_double() * dz + self.origin_z
  (x, y, z)
}

///|
/// Volume per voxel in Å³.
pub fn DensityMap::voxel_volume(self : DensityMap) -> Double {
  let ca = cos_d(self.gamma)
  let cb = cos_d(self.beta)
  let cg = cos_d(self.alpha)
  let vol = self.xlen * self.ylen * self.zlen * (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
  vol / (self.mx * self.my * self.mz).to_double()
}

///|
/// Total volume above contour level in Å³.
pub fn DensityMap::volume_above(self : DensityMap, level : Double) -> Double {
  self.count_above(level).to_double() * self.voxel_volume()
}
