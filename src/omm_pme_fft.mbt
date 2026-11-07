/// omm_pme_fft.mbt — grid-based (FFT) reciprocal-space Ewald sum.
///
/// Item 9: replaces the direct k-space summation of `omm_pme_reciprocal_space`
/// with a true particle-mesh path:
///   1. charges are interpolated onto a regular grid (order-2 linear
///      weighting),
///   2. a 3D FFT transforms the charge density,
///   3. the reciprocal energy is accumulated as
///         E = (2*pi/V) * sum_{n != 0} exp(-k^2/(4*alpha^2)) / k^2
///                 * |Q_hat(n) / W(n)|^2 ,
///      where k = 2*pi*(n_x/L_x, n_y/L_y, n_z/L_z) and W(n) is the DFT of
///      the interpolation kernel; dividing by |W|^2 removes the grid
///      interpolation influence so the result matches the analytic sum.
///
/// Grid dimensions must be powers of two (radix-2 Cooley-Tukey FFT).
/// Units: kcal/mol (energy), A (length), e (charge).

///|
/// In-place radix-2 FFT of one line (n must be a power of two). `stride`
/// and `offset` address the line inside the flat 3D arrays. `inverse`
/// flips the twiddle sign and divides by n.
fn omm_fft_line(
  re : Array[Double],
  im : Array[Double],
  n : Int32,
  stride : Int32,
  offset : Int32,
  inverse : Bool,
) -> Unit {
  if n <= 1 {
    return
  }
  // Bit-reversal permutation into work buffers.
  let tre : Array[Double] = Array::make(n, 0.0)
  let tim : Array[Double] = Array::make(n, 0.0)
  let mut bits = 0
  let mut m = n
  while m > 1 {
    m = m / 2
    bits = bits + 1
  }
  for i = 0; i < n; i = i + 1 {
    let mut r = 0
    let mut v = i
    for _b = 0; _b < bits; _b = _b + 1 {
      r = r * 2 + v % 2
      v = v / 2
    }
    tre[i] = re[offset + r * stride]
    tim[i] = im[offset + r * stride]
  }
  // Butterfly passes.
  let mut len = 2
  while len <= n {
    let ang = (if inverse { 2.0 } else { -2.0 }) *
      3.141592653589793 /
      len.to_double()
    let wr0 = @math.cos(ang)
    let wi0 = @math.sin(ang)
    let half = len / 2
    let mut start = 0
    while start < n {
      let mut w_re = 1.0
      let mut w_im = 0.0
      for j = 0; j < half; j = j + 1 {
        let i1 = start + j
        let i2 = start + j + half
        let tr = w_re * tre[i2] - w_im * tim[i2]
        let ti = w_re * tim[i2] + w_im * tre[i2]
        let ur = tre[i1]
        let ui = tim[i1]
        tre[i1] = ur + tr
        tim[i1] = ui + ti
        tre[i2] = ur - tr
        tim[i2] = ui - ti
        let nwr = w_re * wr0 - w_im * wi0
        w_im = w_re * wi0 + w_im * wr0
        w_re = nwr
      }
      start = start + len
    }
    len = len * 2
  }
  let norm = if inverse { 1.0 / n.to_double() } else { 1.0 }
  for i = 0; i < n; i = i + 1 {
    re[offset + i * stride] = tre[i] * norm
    im[offset + i * stride] = tim[i] * norm
  }
}

///|
/// In-place 3D FFT over a flat nx*ny*nz grid (row-major, z fastest).
/// `inverse` normalizes by 1/N and flips the twiddle sign.
pub fn omm_fft_3d(
  re : Array[Double],
  im : Array[Double],
  nx : Int32,
  ny : Int32,
  nz : Int32,
  inverse~ : Bool = false,
) -> Unit {
  // z axis (contiguous)
  for iy = 0; iy < ny; iy = iy + 1 {
    for ix = 0; ix < nx; ix = ix + 1 {
      omm_fft_line(re, im, nz, 1, (ix * ny + iy) * nz, inverse)
    }
  }
  // y axis (stride nz)
  for ix = 0; ix < nx; ix = ix + 1 {
    for iz = 0; iz < nz; iz = iz + 1 {
      omm_fft_line(re, im, ny, nz, ix * ny * nz + iz, inverse)
    }
  }
  // x axis (stride ny*nz)
  for iy = 0; iy < ny; iy = iy + 1 {
    for iz = 0; iz < nz; iz = iz + 1 {
      omm_fft_line(re, im, nx, ny * nz, iy * nz + iz, inverse)
    }
  }
}

///|
/// Smallest power of two >= n (grid sizing helper).
pub fn omm_fft_next_pow2(n : Int32) -> Int32 {
  let mut p = 1
  while p < n {
    p = p * 2
  }
  p
}

///|
/// Assign charges to the grid with order-2 (linear) weighting and return
/// the flat charge-density grid.
fn omm_pme_assign_charges(
  positions : Array[Array[Double]],
  charges : Array[Double],
  params : OmmPMEParams,
  nx : Int32,
  ny : Int32,
  nz : Int32,
) -> Array[Double] {
  let grid : Array[Double] = Array::make(nx * ny * nz, 0.0)
  let box = params.box
  for i = 0; i < positions.length(); i = i + 1 {
    // Fractional coordinates in [0, 1).
    let fx = positions[i][0] / box[0] - (positions[i][0] / box[0]).floor()
    let fy = positions[i][1] / box[1] - (positions[i][1] / box[1]).floor()
    let fz = positions[i][2] / box[2] - (positions[i][2] / box[2]).floor()
    let gx = fx * nx.to_double()
    let gy = fy * ny.to_double()
    let gz = fz * nz.to_double()
    let i0 = gx.floor().to_int() % nx
    let j0 = gy.floor().to_int() % ny
    let k0 = gz.floor().to_int() % nz
    let tx = gx - gx.floor()
    let ty = gy - gy.floor()
    let tz = gz - gz.floor()
    let i1 = (i0 + 1) % nx
    let j1 = (j0 + 1) % ny
    let k1 = (k0 + 1) % nz
    let q = charges[i]
    // trilinear spread over the 8 surrounding nodes
    let wx0 = 1.0 - tx
    let wx1 = tx
    let wy0 = 1.0 - ty
    let wy1 = ty
    let wz0 = 1.0 - tz
    let wz1 = tz
    grid[(i0 * ny + j0) * nz + k0] += q * wx0 * wy0 * wz0
    grid[(i0 * ny + j0) * nz + k1] += q * wx0 * wy0 * wz1
    grid[(i0 * ny + j1) * nz + k0] += q * wx0 * wy1 * wz0
    grid[(i0 * ny + j1) * nz + k1] += q * wx0 * wy1 * wz1
    grid[(i1 * ny + j0) * nz + k0] += q * wx1 * wy0 * wz0
    grid[(i1 * ny + j0) * nz + k1] += q * wx1 * wy0 * wz1
    grid[(i1 * ny + j1) * nz + k0] += q * wx1 * wy1 * wz0
    grid[(i1 * ny + j1) * nz + k1] += q * wx1 * wy1 * wz1
  }
  grid
}

///|
/// sinc(x) = sin(x)/x with sinc(0) = 1.
fn omm_sinc(x : Double) -> Double {
  if x.abs() < 1.0e-12 {
    1.0
  } else {
    @math.sin(x) / x
  }
}

///|
/// Grid-based (FFT) reciprocal-space Ewald energy. Matches
/// `omm_pme_reciprocal_space` up to interpolation and grid-size error.
pub fn omm_pme_grid_reciprocal(
  positions : Array[Array[Double]],
  charges : Array[Double],
  params : OmmPMEParams,
) -> Double {
  let box = params.box
  if box[0] <= 0.0 || box[1] <= 0.0 || box[2] <= 0.0 {
    return 0.0
  }
  let nx = params.grid_x
  let ny = params.grid_y
  let nz = params.grid_z
  let grid_re = omm_pme_assign_charges(positions, charges, params, nx, ny, nz)
  let grid_im : Array[Double] = Array::make(nx * ny * nz, 0.0)
  omm_fft_3d(grid_re, grid_im, nx, ny, nz)
  let vol = box[0] * box[1] * box[2]
  let alpha_sq = params.alpha * params.alpha
  let pi = 3.141592653589793
  let mut energy = 0.0
  // Physical k vectors from the grid indices; use the signed convention
  // n in [-N/2, N/2) (the FFT bins at n > N/2 alias to negative n).
  let half_x = nx / 2
  let half_y = ny / 2
  let half_z = nz / 2
  for ix = 0; ix < nx; ix = ix + 1 {
    let sx = if ix <= half_x { ix } else { ix - nx }
    for iy = 0; iy < ny; iy = iy + 1 {
      let sy = if iy <= half_y { iy } else { iy - ny }
      for iz = 0; iz < nz; iz = iz + 1 {
        let sz = if iz <= half_z { iz } else { iz - nz }
        if sx == 0 && sy == 0 && sz == 0 {
          continue
        }
        let kx = 2.0 * pi * sx.to_double() / box[0]
        let ky = 2.0 * pi * sy.to_double() / box[1]
        let kz = 2.0 * pi * sz.to_double() / box[2]
        let k_sq = kx * kx + ky * ky + kz * kz
        let idx = (ix * ny + iy) * nz + iz
        let q_re = grid_re[idx]
        let q_im = grid_im[idx]
        // Undo the linear-interpolation window: |W|^2 = prod sinc^2(pi*s/N).
        let wx = omm_sinc(pi * sx.to_double() / nx.to_double())
        let wy = omm_sinc(pi * sy.to_double() / ny.to_double())
        let wz = omm_sinc(pi * sz.to_double() / nz.to_double())
        let w2 = wx * wx * wy * wy * wz * wz
        if w2 < 1.0e-14 {
          continue
        }
        let s_sq = (q_re * q_re + q_im * q_im) / w2
        energy = energy + s_sq * @math.exp(-k_sq / (4.0 * alpha_sq)) / k_sq
      }
    }
  }
  energy * 2.0 * pi / vol
}

///|
/// Full PME energy with the FFT-grid reciprocal part.
pub fn omm_pme_grid_energy(
  positions : Array[Array[Double]],
  charges : Array[Double],
  params : OmmPMEParams,
) -> Double {
  let e_real = omm_pme_real_space(positions, charges, params)
  let e_recip = omm_pme_grid_reciprocal(positions, charges, params)
  let e_self = omm_pme_self_energy(charges, params)
  e_real + e_recip + e_self
}
