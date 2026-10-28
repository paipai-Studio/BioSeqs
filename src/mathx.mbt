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
