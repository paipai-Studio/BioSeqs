///|
test "helix_geometry: best_fit_line collinear points" {
  // Points along x-axis
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
    @src.Vector3::new(2.0, 0.0, 0.0),
    @src.Vector3::new(3.0, 0.0, 0.0),
    @src.Vector3::new(4.0, 0.0, 0.0),
  ]
  let (origin, dir) = @src.calc_best_fit_line(coords)
  // Direction should be along x-axis (±1, 0, 0)
  assert_true(dir.x.abs() > 0.99)
  // Origin should be centroid (2, 0, 0)
  assert_true((origin.x - 2.0).abs() < 0.01)
}

///|
test "helix_geometry: best_fit_line z-axis" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 1.0),
    @src.Vector3::new(0.0, 0.0, 2.0),
    @src.Vector3::new(0.0, 0.0, 3.0),
  ]
  let (_, dir) = @src.calc_best_fit_line(coords)
  assert_true(dir.z.abs() > 0.99)
}

///|
test "helix_geometry: best_fit_line two points" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(3.0, 4.0, 0.0),
  ]
  let (origin, dir) = @src.calc_best_fit_line(coords)
  // Direction should be (3,4,0)/5 = (0.6, 0.8, 0)
  assert_true((dir.x.abs() - 0.6).abs() < 0.01)
  assert_true((dir.y.abs() - 0.8).abs() < 0.01)
  // Origin at midpoint
  assert_true((origin.x - 1.5).abs() < 0.01)
  assert_true((origin.y - 2.0).abs() < 0.01)
}

///|
test "helix_geometry: best_fit_line single point" {
  let coords = [@src.Vector3::new(1.0, 2.0, 3.0)]
  let (origin, _) = @src.calc_best_fit_line(coords)
  assert_true((origin.x - 1.0).abs() < 0.01)
  assert_true((origin.y - 2.0).abs() < 0.01)
  assert_true((origin.z - 3.0).abs() < 0.01)
}

///|
test "helix_geometry: best_fit_line empty" {
  let (origin, _) = @src.calc_best_fit_line([])
  assert_true(origin.x.abs() < 0.01)
  assert_true(origin.y.abs() < 0.01)
  assert_true(origin.z.abs() < 0.01)
}

///|
test "helix_geometry: best_fit_plane normal z" {
  // Points in xy-plane → normal should be z-axis
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 1.0, 0.0),
    @src.Vector3::new(1.0, 1.0, 0.0),
    @src.Vector3::new(2.0, 0.5, 0.0),
  ]
  let (_, normal) = @src.calc_best_fit_plane(coords)
  assert_true(normal.z.abs() > 0.99)
}

///|
test "helix_geometry: best_fit_plane normal x" {
  // Points in yz-plane → normal should be x-axis
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 1.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 1.0),
    @src.Vector3::new(0.0, 1.0, 1.0),
  ]
  let (_, normal) = @src.calc_best_fit_plane(coords)
  assert_true(normal.x.abs() > 0.99)
}

///|
test "helix_geometry: helix_axis ideal helix" {
  // Generate ideal alpha-helix CA coordinates (20 residues = ~5.6 turns)
  let coords = Array(capacity=20)
  let radius = 2.3
  let rise = 1.5
  let mut i = 0
  while i < 20 {
    let angle = i.to_double() * 100.0 * 3.14159265358979323846 / 180.0
    let x = radius * @math.cos(angle)
    let y = radius * @math.sin(angle)
    let z = i.to_double() * rise
    coords.push(@src.Vector3::new(x, y, z))
    i = i + 1
  }
  let (_, dir) = @src.calc_helix_axis(coords)
  // Helix axis should be approximately along z
  assert_true(dir.z.abs() > 0.95)
}

///|
test "helix_geometry: helix_kink straight helix" {
  // Long straight helix (24 residues) → small bend at any position
  let coords = Array(capacity=24)
  let mut i = 0
  while i < 24 {
    let angle = i.to_double() * 100.0 * 3.14159265358979323846 / 180.0
    let x = 2.3 * @math.cos(angle)
    let y = 2.3 * @math.sin(angle)
    let z = i.to_double() * 1.5
    coords.push(@src.Vector3::new(x, y, z))
    i = i + 1
  }
  let kink = @src.calc_helix_kink(coords, 12)
  match kink {
    Some(k) => assert_true(k.bend_angle < 20.0)
    None => assert_true(false) // should not be None for valid position
  }
}

///|
test "helix_geometry: helix_kink with bend" {
  // Two helix segments with different axes
  let coords = Array(capacity=0)
  // First segment: along z
  let mut i = 0
  while i < 8 {
    let angle = i.to_double() * 100.0 * 3.14159265358979323846 / 180.0
    coords.push(@src.Vector3::new(
      2.3 * @math.cos(angle),
      2.3 * @math.sin(angle),
      i.to_double() * 1.5,
    ))
    i = i + 1
  }
  // Second segment: tilted 30° from z
  let tilt = 30.0 * 3.14159265358979323846 / 180.0
  i = 0
  while i < 8 {
    let z_local = i.to_double() * 1.5
    let angle = (i + 8).to_double() * 100.0 * 3.14159265358979323846 / 180.0
    let x_local = 2.3 * @math.cos(angle)
    let y_local = 2.3 * @math.sin(angle)
    // Rotate around y by tilt
    let x = x_local * @math.cos(tilt) + z_local * @math.sin(tilt)
    let y = y_local
    let z = -x_local * @math.sin(tilt) + z_local * @math.cos(tilt) + 12.0
    coords.push(@src.Vector3::new(x, y, z))
    i = i + 1
  }
  let kink = @src.calc_helix_kink(coords, 8)
  match kink {
    Some(k) => assert_true(k.bend_angle > 15.0)
    None => assert_true(false)
  }
}

///|
test "helix_geometry: helix_kink edge cases" {
  // Too few residues
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 1.0),
    @src.Vector3::new(0.0, 0.0, 2.0),
    @src.Vector3::new(0.0, 0.0, 3.0),
  ]
  // kink at position 2: only 2 residues on each side
  let kink = @src.calc_helix_kink(coords, 2)
  // Need at least 3 on each side, so should be None
  match kink {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "helix_geometry: helix_kink position out of range" {
  let coords = Array(capacity=10)
  let mut i = 0
  while i < 10 {
    coords.push(@src.Vector3::new(0.0, 0.0, i.to_double()))
    i = i + 1
  }
  // Position too close to end
  let kink = @src.calc_helix_kink(coords, 8)
  match kink {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "helix_geometry: analyze_helix_kinks straight" {
  // Long straight helix — no kinks > 20°
  let coords = Array(capacity=24)
  let mut i = 0
  while i < 24 {
    let angle = i.to_double() * 100.0 * 3.14159265358979323846 / 180.0
    coords.push(@src.Vector3::new(
      2.3 * @math.cos(angle),
      2.3 * @math.sin(angle),
      i.to_double() * 1.5,
    ))
    i = i + 1
  }
  let kinks = @src.analyze_helix_kinks(coords)
  // Straight helix should have no kinks > 10°
  // With 24 residues (6.67 turns), axis fitting is accurate
  assert_eq(kinks.length(), 0)
}

///|
test "helix_geometry: analyze_helix_kinks with bend" {
  let coords = Array(capacity=0)
  let mut i = 0
  while i < 8 {
    let angle = i.to_double() * 100.0 * 3.14159265358979323846 / 180.0
    coords.push(@src.Vector3::new(
      2.3 * @math.cos(angle),
      2.3 * @math.sin(angle),
      i.to_double() * 1.5,
    ))
    i = i + 1
  }
  let tilt = 40.0 * 3.14159265358979323846 / 180.0
  i = 0
  while i < 8 {
    let z_local = i.to_double() * 1.5
    let angle = (i + 8).to_double() * 100.0 * 3.14159265358979323846 / 180.0
    let x_local = 2.3 * @math.cos(angle)
    let y_local = 2.3 * @math.sin(angle)
    let x = x_local * @math.cos(tilt) + z_local * @math.sin(tilt)
    let y = y_local
    let z = -x_local * @math.sin(tilt) + z_local * @math.cos(tilt) + 12.0
    coords.push(@src.Vector3::new(x, y, z))
    i = i + 1
  }
  let kinks = @src.analyze_helix_kinks(coords)
  assert_true(kinks.length() > 0)
}

///|
test "helix_geometry: best_fit_plane single point" {
  let (origin, _) = @src.calc_best_fit_plane([@src.Vector3::new(1.0, 2.0, 3.0)])
  assert_true((origin.x - 1.0).abs() < 0.01)
  assert_true((origin.y - 2.0).abs() < 0.01)
}

///|
test "helix_geometry: best_fit_line direction normalized" {
  let coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(10.0, 0.0, 0.0),
    @src.Vector3::new(20.0, 0.0, 0.0),
  ]
  let (_, dir) = @src.calc_best_fit_line(coords)
  let mag = (dir.x * dir.x + dir.y * dir.y + dir.z * dir.z).sqrt()
  assert_true((mag - 1.0).abs() < 0.01)
}
