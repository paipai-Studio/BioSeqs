///|
/// Shared numeric utilities (converged from per-module prefixed duplicates).

///|
/// NaN 与超出 ±1.0e300 的值视为非有限(与各模块原实现的多数派语义一致).
fn is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
/// 将数值限制在 [lower, upper] 区间内.
fn clamp(value : Double, lower : Double, upper : Double) -> Double {
  value.max(lower).min(upper)
}

///|
/// 权威 log-Gamma(x)(全库统一实现, 收敛各模块 16 份副本).
///
/// Lanczos g=7 近似(9 系数, 精度 ~1e-13); x < 0.5 走反射公式
/// lgamma(x) = ln(pi/sin(pi x)) - lgamma(1-x); x <= 0 或 NaN 返回 0.0
/// (多数派模块约定, 调用方均传正域参数).
fn lgamma(x : Double) -> Double {
  if x != x || x <= 0.0 {
    return 0.0
  }
  let pi = 3.14159265358979323846
  if x < 0.5 {
    return @math.ln(pi / @math.sin(pi * x)) - lgamma(1.0 - x)
  }
  let c = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7,
  ]
  let shifted = x - 1.0
  let mut series = c[0]
  for i in 1..<9 {
    series = series + c[i] / (shifted + i.to_double())
  }
  let t = shifted + 7.5
  0.9189385332046727 + (shifted + 0.5) * @math.ln(t) - t + @math.ln(series)
}
