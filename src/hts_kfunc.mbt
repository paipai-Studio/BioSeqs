///|
/// htslib math functions (kfunc.c port).
/// Provides log-gamma, complementary error function, regularized incomplete
/// gamma and beta functions used in statistical tests.

///|
/// Log gamma function: log(Gamma(z)).
/// AS245, 2nd algorithm.
pub fn kf_lgamma(z : Double) -> Double {
  let x0 = 0.1659470187408462e-06 / (z + 7.0)
  let x1 = 0.9934937113930748e-05 / (z + 6.0)
  let x2 = 0.1385710331296526 / (z + 5.0)
  let x3 = 12.50734324009056 / (z + 4.0)
  let x4 = 176.6150291498386 / (z + 3.0)
  let x5 = 771.3234287757674 / (z + 2.0)
  let x6 = 1259.139216722289 / (z + 1.0)
  let x7 = 676.5203681218835 / z
  let x = x0 + x1 - x2 + x3 - x4 + x5 - x6 + x7 + 0.9999999999995183
  @math.ln(x) - 5.58106146679532777 - z + (z - 0.5) * @math.ln(z + 6.5)
}

///|
/// Complementary error function: erfc(x) = 2/sqrt(pi) * integral(x, inf, e^(-t^2)).
/// AS66, 2nd algorithm.
pub fn kf_erfc(x : Double) -> Double {
  let p0 = 220.2068679123761
  let p1 = 221.2135961699311
  let p2 = 112.0792914978709
  let p3 = 33.912866078383
  let p4 = 6.37396220353165
  let p5 = 0.7003830644436881
  let p6 = 0.03526249659989109
  let q0 = 440.4137358247522
  let q1 = 793.8265125199484
  let q2 = 637.3336333788311
  let q3 = 296.5642487796737
  let q4 = 86.78073220294608
  let q5 = 16.06417757920695
  let q6 = 1.755667163182642
  let q7 = 0.08838834764831844
  let sqrt2 = 1.4142135623730951
  let z = x.abs() * sqrt2
  if z > 37.0 {
    return if x > 0.0 { 0.0 } else { 2.0 }
  }
  let expntl = @math.exp(z * z * -0.5)
  let p = if z < 10.0 / sqrt2 {
    let num = (((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0
    let den = ((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0
    expntl * num / den
  } else {
    expntl / 2.506628274631001 / (z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 0.65)))))
  }
  if x > 0.0 { 2.0 * p } else { 2.0 * (1.0 - p) }
}

let kf_gamma_eps : Double = 1.0e-14

let kf_tiny : Double = 1.0e-290

///|
/// Regularized lower incomplete gamma function, by series expansion.
fn kf_gammap_series(s : Double, z : Double) -> Double {
  let mut sum = 1.0
  let mut x = 1.0
  let mut k = 1
  while k < 100 {
    x = x * z / (s + k.to_double())
    sum = sum + x
    if x / sum < kf_gamma_eps {
      break
    }
    k = k + 1
  }
  @math.exp(s * @math.ln(z) - z - kf_lgamma(s + 1.0) + @math.ln(sum))
}

///|
/// Regularized upper incomplete gamma function, by continued fraction.
/// Uses modified Lentz's algorithm.
fn kf_gammaq_contfrac(s : Double, z : Double) -> Double {
  let mut f = 1.0 + z - s
  let mut c = f
  let mut d = 0.0
  let mut j = 1
  while j < 100 {
    let a = j.to_double() * (s - j.to_double())
    let b = (2 * j).to_double() + 1.0 + z - s
    d = b + a * d
    if d < kf_tiny {
      d = kf_tiny
    }
    c = b + a / c
    if c < kf_tiny {
      c = kf_tiny
    }
    d = 1.0 / d
    let dd = c * d
    f = f * dd
    if (dd - 1.0).abs() < kf_gamma_eps {
      break
    }
    j = j + 1
  }
  @math.exp(s * @math.ln(z) - z - kf_lgamma(s) - @math.ln(f))
}

///|
/// Regularized lower incomplete gamma function P(s, z).
pub fn kf_gammap(s : Double, z : Double) -> Double {
  if z <= 1.0 || z < s {
    kf_gammap_series(s, z)
  } else {
    1.0 - kf_gammaq_contfrac(s, z)
  }
}

///|
/// Regularized upper incomplete gamma function Q(s, z).
pub fn kf_gammaq(s : Double, z : Double) -> Double {
  if z <= 1.0 || z < s {
    1.0 - kf_gammap_series(s, z)
  } else {
    kf_gammaq_contfrac(s, z)
  }
}

///|
/// Regularized incomplete beta function I_x(a, b).
/// Uses continued fraction expansion (Numerical Recipes in C, section 6.4).
fn kf_betai_cf(a : Double, b : Double, x : Double) -> Double {
  if x == 0.0 {
    return 0.0
  }
  if x == 1.0 {
    return 1.0
  }
  let mut f = 1.0
  let mut c = f
  let mut d = 0.0
  let mut j = 1
  while j < 200 {
    let m = j / 2
    let md = m.to_double()
    let aa = if j % 2 == 1 {
      -(a + md) * (a + b + md) * x / ((a + 2.0 * md) * (a + 2.0 * md + 1.0))
    } else {
      md * (b - md) * x / ((a + 2.0 * md - 1.0) * (a + 2.0 * md))
    }
    d = 1.0 + aa * d
    if d < kf_tiny {
      d = kf_tiny
    }
    c = 1.0 + aa / c
    if c < kf_tiny {
      c = kf_tiny
    }
    d = 1.0 / d
    let dd = c * d
    f = f * dd
    if (dd - 1.0).abs() < kf_gamma_eps {
      break
    }
    j = j + 1
  }
  @math.exp(kf_lgamma(a + b) - kf_lgamma(a) - kf_lgamma(b) + a * @math.ln(x) + b * @math.ln(1.0 - x)) / a / f
}

///|
/// Regularized incomplete beta function I_x(a, b).
pub fn kf_betai(a : Double, b : Double, x : Double) -> Double {
  if x < (a + 1.0) / (a + b + 2.0) {
    kf_betai_cf(a, b, x)
  } else {
    1.0 - kf_betai_cf(b, a, 1.0 - x)
  }
}
