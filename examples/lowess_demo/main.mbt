///|
/// LOWESS Demo.
///
/// Demonstrates Bio.Statistics.lowess module: locally weighted regression
/// for scatterplot smoothing. Inspired by Biopython's
/// Bio.Statistics.lowess module (Cleveland 1979).

/// Sine function via Taylor series with argument reduction to [-pi, pi].
fn local_sin(x : Double) -> Double {
  let pi = 3.141592653589793
  let two_pi = 6.283185307179586
  let mut y = x
  while y > pi {
    y = y - two_pi
  }
  while y < -pi {
    y = y + two_pi
  }
  let y2 = y * y
  let mut term = y
  let mut sum = y
  let mut k = 1
  while k <= 15 {
    term = -term * y2 / ((2 * k).to_double() * (2 * k + 1).to_double())
    sum = sum + term
    k = k + 1
  }
  sum
}

/// Deterministic pseudo-noise in [-0.2, 0.2] from an integer index.
fn pseudo_noise(i : Int) -> Double {
  let n = (i * 9301 + 49297) % 233280
  (n.to_double() / 233280.0 - 0.5) * 0.4
}

/// Format a Double to a fixed number of decimal places.
fn fmt(value : Double, decimals : Int) -> String {
  let s = value.to_string()
  let parts = s.split(".").to_array()
  let int_part = if parts.length() > 0 { parts[0].to_owned() } else { "0" }
  let frac_part = if parts.length() > 1 { parts[1].to_owned() } else { "" }
  let mut padded_frac = frac_part
  if padded_frac.length() < decimals {
    padded_frac = padded_frac + "0".repeat(decimals - padded_frac.length())
  } else if padded_frac.length() > decimals {
    padded_frac = padded_frac[0:decimals].to_owned()
  }
  int_part + "." + padded_frac
}

fn main {
  println("=== LOWESS Demo ===")
  println("Bio.Statistics.lowess: Locally Weighted Regression")
  println("Inspired by Biopython's Bio.Statistics.lowess module\n")

  // === 1. Generate noisy sine wave data ===
  println("=== 1. Generate Noisy Sine Wave Data ===")
  let n = 101
  let x : Array[Double] = Array::make(n, 0.0)
  let y : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    x[i] = i.to_double() * 0.1
    y[i] = local_sin(x[i]) + pseudo_noise(i)
    i = i + 1
  }
  println("  Generated \{n} points: x = 0, 0.1, ..., 10; y = sin(x) + noise")

  // === 2. LOWESS smoothing, print before/after ===
  println("\n=== 2. LOWESS Smoothing (every 5th point) ===")
  let smoothed = @bio.lowess(x, y, f=0.3, n_steps=0, delta=0.0)
  println("  x       y (noisy)    smoothed")
  i = 0
  while i < n {
    if i % 5 == 0 {
      println("  " + fmt(x[i], 2) + "     " + fmt(y[i], 6) + "    " + fmt(smoothed[i], 6))
    }
    i = i + 1
  }

  // === 3. Robust LOWESS with outliers ===
  println("\n=== 3. Robust LOWESS vs Standard LOWESS (with outliers) ===")
  let y_outlier = y.copy()
  y_outlier[20] = 3.0
  y_outlier[50] = -3.0
  y_outlier[80] = 4.0
  println("  Injected outliers at i=20 (y=3.0), i=50 (y=-3.0), i=80 (y=4.0)")
  let smoothed_std = @bio.lowess(x, y_outlier, f=0.3, n_steps=0, delta=0.0)
  let smoothed_rob = @bio.lowess_robust(x, y_outlier, f=0.3, n_steps=2, delta=0.0)
  println("  x       standard      robust")
  i = 0
  while i < n {
    if i % 5 == 0 {
      println(
        "  " +
        fmt(x[i], 2) +
        "     " +
        fmt(smoothed_std[i], 6) +
        "    " +
        fmt(smoothed_rob[i], 6),
      )
    }
    i = i + 1
  }
  println("  (robust LOWESS is less affected by the outlier points)")

  // === 4. lowess_smooth_series on time series ===
  println("\n=== 4. lowess_smooth_series on Time Series ===")
  let series = [
    1.0, 3.0, 2.0, 5.0, 4.0, 6.0, 5.0, 8.0, 7.0, 9.0, 8.0, 11.0, 10.0, 12.0,
    11.0, 14.0, 13.0, 15.0, 14.0, 17.0,
  ]
  let series_smoothed = @bio.lowess_smooth_series(series, window_frac=0.4)
  println("  idx   original   smoothed")
  i = 0
  while i < series.length() {
    println("  \{i}      " + fmt(series[i], 2) + "        " + fmt(series_smoothed[i], 4))
    i = i + 1
  }

  // === 5. lowess_predict on new points ===
  println("\n=== 5. lowess_predict on New Points ===")
  let pred_x = [0.05, 2.5, 5.0, 7.5, 9.95]
  println("  x_new   predicted")
  i = 0
  while i < pred_x.length() {
    let pred = @bio.lowess_predict(pred_x[i], x, y, smoothed, f=0.3)
    println("  " + fmt(pred_x[i], 2) + "      " + fmt(pred, 6))
    i = i + 1
  }

  // === 6. tricube_weight and bisquare_weight ===
  println("\n=== 6. Weight Functions ===")
  let u_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
  println("  u       tricube      bisquare")
  i = 0
  while i < u_vals.length() {
    let tw = @bio.lowess_tricube_weight(u_vals[i])
    let bw = @bio.lowess_bisquare_weight(u_vals[i])
    println("  " + fmt(u_vals[i], 2) + "     " + fmt(tw, 6) + "    " + fmt(bw, 6))
    i = i + 1
  }

  println("\n=== LOWESS Demo Complete ===")
}
