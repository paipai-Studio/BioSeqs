///|
/// Shared test helpers. New test files should reference these instead of
/// defining per-file `*_close` helpers — which had previously collided as
/// "toplevel identifier declared twice" across files.

///|
/// Approximate equality within a tolerance.
pub fn assert_close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

///|
/// Sign-agnostic approximate equality (compares magnitudes).
pub fn assert_close_mag(a : Double, b : Double, tol : Double) -> Bool {
  (a.abs() - b.abs()).abs() < tol
}

///|
/// Approximate equality with a fixed 1e-6 tolerance.
pub fn assert_close_e6(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-6
}

///|
/// Approximate equality with a fixed 1e-9 tolerance.
pub fn assert_close_e9(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-9
}