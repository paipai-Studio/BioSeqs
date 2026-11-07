/// Umbrella sampling for enhanced sampling along a reaction coordinate.
///
/// OpenMM equivalent: `CustomCVForce` with umbrella bias
///
/// Adds harmonic bias potential along a collective variable (reaction
/// coordinate), then reconstructs the PMF using WHAM.

///|
/// Umbrella sampling window
pub struct OmmUmbrellaWindow {
  center : Double // bias center
  k : Double // force constant
  samples : Array[Double] // collected CV values
}

///|
pub fn OmmUmbrellaWindow::new(
  center~ : Double,
  k~ : Double,
) -> OmmUmbrellaWindow {
  { center, k, samples: [], }
}

///|
pub fn OmmUmbrellaWindow::add_sample(
  self : OmmUmbrellaWindow,
  cv_value : Double,
) -> Unit {
  self.samples.push(cv_value)
}

///|
pub fn OmmUmbrellaWindow::center(self : OmmUmbrellaWindow) -> Double {
  self.center
}

///|
pub fn OmmUmbrellaWindow::k(self : OmmUmbrellaWindow) -> Double {
  self.k
}

///|
pub fn OmmUmbrellaWindow::n_samples(self : OmmUmbrellaWindow) -> Int32 {
  self.samples.length()
}

///|
/// Compute umbrella bias energy: U_bias = 0.5 * k * (cv - center)^2
pub fn omm_umbrella_energy(cv : Double, center : Double, k : Double) -> Double {
  0.5 * k * (cv - center) * (cv - center)
}

///|
/// Compute umbrella bias force on the CV: F = -k * (cv - center)
pub fn omm_umbrella_force(cv : Double, center : Double, k : Double) -> Double {
  -k * (cv - center)
}

///|
/// WHAM (Weighted Histogram Analysis Method) to reconstruct PMF.
///
/// windows: list of umbrella sampling windows with collected samples
/// n_bins: number of histogram bins
/// range_min, range_max: CV range for histogram
/// temperature: system temperature
/// max_iter: WHAM iterations
/// tolerance: convergence tolerance
///
/// Returns (bin_centers, pmf_values)
pub fn omm_wham(
  windows : Array[OmmUmbrellaWindow],
  n_bins : Int32,
  range_min : Double,
  range_max : Double,
  temperature~ : Double = 300.0,
  max_iter~ : Int32 = 10000,
  tolerance~ : Double = 0.0001,
) -> (Array[Double], Array[Double]) {
  let n_windows = windows.length()
  if n_windows == 0 || n_bins < 2 {
    return ([], [])
  }
  let kT = 0.001987204 * temperature
  let bin_width = (range_max - range_min) / n_bins.to_double()
  // Bin centers
  let bin_centers : Array[Double] = Array::make(n_bins, 0.0)
  for i = 0; i < n_bins; i = i + 1 {
    bin_centers[i] = range_min + (i.to_double() + 0.5) * bin_width
  }
  // Build histograms for each window
  let histograms : Array[Array[Double]] = Array::make(
    n_windows,
    Array::make(n_bins, 0.0),
  )
  for w = 0; w < n_windows; w = w + 1 {
    for s in windows[w].samples {
      let bin_idx = ((s - range_min) / bin_width).to_int()
      if bin_idx >= 0 && bin_idx < n_bins {
        histograms[w][bin_idx] = histograms[w][bin_idx] + 1.0
      }
    }
  }
  // WHAM iteration
  let f_const : Array[Double] = Array::make(n_windows, 0.0) // free energy constants
  let large_val = 10000000000.0
  let large_val2 = 1000000000.0
  let mut pmf : Array[Double] = Array::make(n_bins, 0.0)
  for iter = 0; iter < max_iter; iter = iter + 1 {
    // Compute probability at each bin
    let mut max_change = 0.0
    let new_pmf : Array[Double] = Array::make(n_bins, 0.0)
    for b = 0; b < n_bins; b = b + 1 {
      let cv = bin_centers[b]
      // Denominator: sum_w N_w * exp(-beta*(U_w(cv) - F_w))
      let mut denom = 0.0
      for w = 0; w < n_windows; w = w + 1 {
        let n_w = windows[w].n_samples().to_double()
        let u_w = omm_umbrella_energy(cv, windows[w].center, windows[w].k)
        denom = denom + n_w * @math.exp(-(u_w - f_const[w]) / kT)
      }
      // Numerator: sum_w h_w(b)
      let mut numer = 0.0
      for w = 0; w < n_windows; w = w + 1 {
        numer = numer + histograms[w][b]
      }
      if denom > 0.0 && numer > 0.0 {
        new_pmf[b] = -kT * @math.ln(numer / denom)
      } else {
        new_pmf[b] = large_val // inaccessible
      }
    }
    // Update F constants
    for w = 0; w < n_windows; w = w + 1 {
      let mut sum_exp = 0.0
      for b = 0; b < n_bins; b = b + 1 {
        let cv = bin_centers[b]
        let u_w = omm_umbrella_energy(cv, windows[w].center, windows[w].k)
        if new_pmf[b] < large_val2 {
          sum_exp = sum_exp + @math.exp(-(new_pmf[b] - u_w) / kT) * bin_width
        }
      }
      if sum_exp > 0.0 {
        let new_f = -kT * @math.ln(sum_exp)
        let diff = new_f - f_const[w]
        if diff < 0.0 {
          if -diff > max_change {
            max_change = -diff
          }
        } else if diff > max_change {
          max_change = diff
        }
        f_const[w] = new_f
      }
    }
    pmf = new_pmf
    if max_change < tolerance {
      break
    }
  }
  // Shift PMF so minimum is 0
  let mut min_pmf = large_val
  for b = 0; b < n_bins; b = b + 1 {
    if pmf[b] < min_pmf {
      min_pmf = pmf[b]
    }
  }
  for b = 0; b < n_bins; b = b + 1 {
    pmf[b] = pmf[b] - min_pmf
  }
  (bin_centers, pmf)
}
