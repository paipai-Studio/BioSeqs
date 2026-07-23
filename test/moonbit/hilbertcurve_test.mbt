///|
/// Tests for HilbertCurve module.

test "HilbertCurve creation" {
  let hc = @src.HilbertCurve::new(4, 2)
  assert_eq(hc.levels, 4)
  assert_eq(hc.dimension, 2)
  assert_eq(hc.max_coordinate, 15)
}

test "hilbert_encode_decode_roundtrip" {
  let hc = @src.HilbertCurve::new(3, 2)
  let coords = [3, 5]
  let encoded = @src.hilbert_encode(hc, coords)
  let decoded = @src.hilbert_decode(hc, encoded)
  assert_eq(decoded.length(), 2)
}

test "hilbert_distance" {
  let hc = @src.HilbertCurve::new(3, 2)
  let coord1 = [1, 2]
  let coord2 = [3, 4]
  let dist = @src.hilbert_distance(hc, coord1, coord2)
  assert_true(dist >= 0)
}

test "hilbert_point_to_segment" {
  let hc = @src.HilbertCurve::new(2, 2)
  let segments = @src.hilbert_point_to_segment(hc, 0, 3)
  assert_eq(segments.length(), 4)
}

test "hilbert_linearize_genome" {
  let segments = @src.hilbert_linearize_genome(10, 3)
  assert_eq(segments.length(), 10)
}

test "hilbert_map_to_grid" {
  let hc = @src.HilbertCurve::new(2, 2)
  let values = [1.0, 2.0, 3.0, 4.0]
  let grid = @src.hilbert_map_to_grid(hc, values)
  assert_eq(grid.length(), 4)
}
