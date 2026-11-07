///|
/// Site-specific radial distribution functions, mirroring MDAnalysis
/// `analysis.rdf.InterRDF_s` (site-wise RDF between atom groups) and a 2D
/// binned variant with quadrant folding and optional periodic minimum
/// images.

///|
/// Site-specific (pair-resolved) radial distribution counts
/// (MDAnalysis `InterRDF_s`).
///
/// For every frame, every pair `(sel1[i], sel2[j])` contributes 1 to
/// `counts[i * |sel2| + j][bin]` when their distance is `< range_max`
/// (pair order: `sel1` outer, `sel2` inner). The returned counts are raw
/// pair counts (not normalized); use `mda_rdf_sitewise_normalize` to turn
/// them into `g(r)`.
///
/// Returns `(bin_centers, counts)` where `bin_centers[b]` is the midpoint
/// of bin `b` and `counts` has one row of `nbins` entries per pair.
pub fn mda_rdf_sitewise(
  traj : Array[Array[Array[Double]]],
  sel1 : Array[Int32],
  sel2 : Array[Int32],
  nbins : Int32,
  range_max : Double,
) -> (Array[Double], Array[Array[Double]]) {
  let nb = nbins
  let dr = range_max / nbins.to_double()
  let np1 = sel1.length()
  let np2 = sel2.length()
  let npair = np1 * np2
  let counts : Array[Array[Double]] = Array(capacity=npair)
  for p = 0; p < npair; p = p + 1 {
    counts.push(Array::make(nb, 0.0))
  }
  for frame in traj {
    for i = 0; i < np1; i = i + 1 {
      let a = frame[sel1[i]]
      for j = 0; j < np2; j = j + 1 {
        let b = frame[sel2[j]]
        let d = mda_euclid(a, b)
        if d < range_max {
          let mut bin = (d / dr).to_int()
          if bin >= nb {
            bin = nb - 1
          }
          counts[i * np2 + j][bin] = counts[i * np2 + j][bin] + 1.0
        }
      }
    }
  }
  let centers : Array[Double] = Array(capacity=nb)
  for b = 0; b < nb; b = b + 1 {
    centers.push((b.to_double() + 0.5) * dr)
  }
  (centers, counts)
}

///|
/// Normalize raw site-wise pair counts into the radial distribution
/// function `g(r)`. For each pair row and bin:
///
/// `g[pair][bin] = counts[pair][bin] / (n_frames * shell_volume * rho2)`
///
/// where `shell_volume = (4/3) * pi * (r_out^3 - r_in^3)` is the volume of
/// the spherical shell of the bin and `rho2 = n2 / v_total` with
/// `v_total = (4/3) * pi * range_max^3` is the number density of the
/// `sel2` species inside the sampling sphere. This is the ideal-gas
/// reference for a single `sel1` site; pairs are normalized per pair, so
/// site-specific enrichment/depletion is directly visible.
pub fn mda_rdf_sitewise_normalize(
  counts : Array[Array[Double]],
  n1 : Int32,
  n2 : Int32,
  n_frames : Int32,
  nbins : Int32,
  range_max : Double,
) -> Array[Array[Double]] {
  let _ = n1
  let nb = nbins
  let dr = range_max / nbins.to_double()
  let v_total = 4.1887902047863905 * range_max * range_max * range_max
  let rho2 = n2.to_double() / v_total
  let nf = n_frames.to_double()
  let g : Array[Array[Double]] = Array(capacity=counts.length())
  for row in counts {
    let out_row : Array[Double] = Array(capacity=nb)
    for b = 0; b < nb; b = b + 1 {
      let r_in = b.to_double() * dr
      let r_out = (b + 1).to_double() * dr
      let shell = 4.1887902047863905 *
        (r_out * r_out * r_out - r_in * r_in * r_in)
      let denom = nf * shell * rho2
      out_row.push(if denom > 0.0 { row[b] / denom } else { 0.0 })
    }
    g.push(out_row)
  }
  g
}

///|
/// 2D binning of absolute pairwise displacements with quadrant folding:
/// for every frame and every pair `(sel1[i], sel2[j])` the displacement is
/// computed (minimum-image folded when `box` is given, orthorhombic edge
/// lengths `[Lx, Ly, Lz]`), reduced to `(|dx|, |dy|)`, and accumulated
/// into the histogram when both fall inside `[0, x_max) x [0, y_max)`.
///
/// Returns `(x_centers, y_centers, hist)` with `hist[bx][by]` a raw
/// (unnormalized) count; `x_centers`/`y_centers` are bin midpoints.
pub fn mda_rdf_2d(
  traj : Array[Array[Array[Double]]],
  sel1 : Array[Int32],
  sel2 : Array[Int32],
  x_max : Double,
  y_max : Double,
  nbins_x : Int32,
  nbins_y : Int32,
  box~ : Array[Double]? = None,
) -> (Array[Double], Array[Double], Array[Array[Double]]) {
  let nx = nbins_x
  let ny = nbins_y
  let dx = x_max / nbins_x.to_double()
  let dy = y_max / nbins_y.to_double()
  let hist : Array[Array[Double]] = Array(capacity=nx)
  for i = 0; i < nx; i = i + 1 {
    hist.push(Array::make(ny, 0.0))
  }
  for frame in traj {
    for i in sel1 {
      let a = frame[i]
      for j in sel2 {
        let b = frame[j]
        let disp = match box {
          Some(bx) => mda_guess_min_image(a, b, bx)
          None => [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
        }
        let ax = disp[0].abs()
        let ay = disp[1].abs()
        if ax < x_max && ay < y_max {
          let mut bx = (ax / dx).to_int()
          let mut by = (ay / dy).to_int()
          if bx >= nx {
            bx = nx - 1
          }
          if by >= ny {
            by = ny - 1
          }
          hist[bx][by] = hist[bx][by] + 1.0
        }
      }
    }
  }
  let xc : Array[Double] = Array(capacity=nx)
  let yc : Array[Double] = Array(capacity=ny)
  for b = 0; b < nx; b = b + 1 {
    xc.push((b.to_double() + 0.5) * dx)
  }
  for b = 0; b < ny; b = b + 1 {
    yc.push((b.to_double() + 0.5) * dy)
  }
  (xc, yc, hist)
}
