///|
/// Tests for density_map.mbt — CCP4/MRC parsing, X-PLOR parsing,
/// trilinear interpolation, statistics, contour, peaks, arithmetic,
/// correlation, density from atoms, mask, volume.

test "dm_simple_create" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m = @src.DensityMap::simple(
    title="test",
    nx=2,
    ny=2,
    nz=2,
    xlen=10.0,
    ylen=10.0,
    zlen=10.0,
    data=data,
  )
  assert_eq(m.get_title(), "test")
  assert_eq(m.get_shape(), (2, 2, 2))
  let (xl, yl, zl, _, _, _) = m.get_cell()
  assert_eq(xl, 10.0)
  assert_eq(yl, 10.0)
  assert_eq(zl, 10.0)
  assert_eq(m.num_voxels(), 8)
}

test "dm_get_grid" {
  let data = [
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
  ]
  let m = @src.DensityMap::simple(
    title="t",
    nx=3,
    ny=2,
    nz=2,
    xlen=30.0,
    ylen=20.0,
    zlen=20.0,
    data=data,
  )
  assert_eq(m.get(0, 0, 0), 1.0)
  assert_eq(m.get(2, 0, 0), 3.0)
  assert_eq(m.get(0, 1, 0), 4.0)
  assert_eq(m.get(0, 0, 1), 7.0)
  assert_eq(m.get(2, 1, 1), 12.0)
  let oob = m.get(5, 5, 5)
  assert_true(oob != oob)
}

test "dm_stats" {
  let data = [1.0, 2.0, 3.0, 4.0]
  let m = @src.DensityMap::simple(
    title="t",
    nx=2,
    ny=2,
    nz=1,
    xlen=10.0,
    ylen=10.0,
    zlen=10.0,
    data=data,
  )
  let (mn, mx, mean, std, _rms) = m.stats()
  assert_eq(mn, 1.0)
  assert_eq(mx, 4.0)
  assert_eq(mean, 2.5)
  assert_true((std - 1.118).abs() < 0.001)
}

test "dm_sample_grid_corner" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m = @src.DensityMap::simple(
    title="t",
    nx=2,
    ny=2,
    nz=2,
    xlen=10.0,
    ylen=10.0,
    zlen=10.0,
    data=data,
  )
  let v = m.sample_grid(0.0, 0.0, 0.0)
  assert_eq(v, 1.0)
  let v2 = m.sample_grid(1.0, 0.0, 0.0)
  assert_eq(v2, 2.0)
}

test "dm_sample_grid_interp" {
  let data = [0.0, 10.0, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0]
  let m = @src.DensityMap::simple(
    title="t",
    nx=2,
    ny=2,
    nz=2,
    xlen=10.0,
    ylen=10.0,
    zlen=10.0,
    data=data,
  )
  let v = m.sample_grid(0.5, 0.0, 0.0)
  assert_true((v - 5.0).abs() < 1.0e-10)
}

test "dm_sample_grid_periodic" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m = @src.DensityMap::simple(
    title="t",
    nx=2,
    ny=2,
    nz=2,
    xlen=10.0,
    ylen=10.0,
    zlen=10.0,
    data=data,
  )
  let v = m.sample_grid(2.0, 0.0, 0.0)
  assert_eq(v, 1.0)
}

test "dm_count_above" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m = @src.DensityMap::simple(
    title="t",
    nx=2,
    ny=2,
    nz=2,
    xlen=10.0,
    ylen=10.0,
    zlen=10.0,
    data=data,
  )
  assert_eq(m.count_above(4.5), 4)
  assert_eq(m.count_above(0.0), 8)
  assert_eq(m.count_above(10.0), 0)
}

test "dm_contour_points" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m = @src.DensityMap::simple(
    title="t",
    nx=2,
    ny=2,
    nz=2,
    xlen=10.0,
    ylen=10.0,
    zlen=10.0,
    data=data,
  )
  let pts = m.contour_points(4.5)
  assert_eq(pts.length(), 4)
}

test "dm_find_peaks" {
  let data = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  ]
  let m = @src.DensityMap::simple(
    title="t",
    nx=3,
    ny=3,
    nz=3,
    xlen=30.0,
    ylen=30.0,
    zlen=30.0,
    data=data,
  )
  let peaks = m.find_peaks(5.0)
  assert_eq(peaks.length(), 1)
  let (pi, pj, pk, pv) = peaks[0]
  assert_eq(pi, 1)
  assert_eq(pj, 1)
  assert_eq(pk, 1)
  assert_eq(pv, 10.0)
}

test "dm_add" {
  let d1 = [1.0, 2.0, 3.0, 4.0]
  let d2 = [10.0, 20.0, 30.0, 40.0]
  let m1 = @src.DensityMap::simple(
    title="a", nx=2, ny=2, nz=1, xlen=10.0, ylen=10.0, zlen=10.0, data=d1,
  )
  let m2 = @src.DensityMap::simple(
    title="b", nx=2, ny=2, nz=1, xlen=10.0, ylen=10.0, zlen=10.0, data=d2,
  )
  let m3 = m1.add(m2)
  assert_eq(m3.get(0, 0, 0), 11.0)
  assert_eq(m3.get(1, 1, 0), 44.0)
}

test "dm_subtract" {
  let d1 = [10.0, 20.0, 30.0, 40.0]
  let d2 = [1.0, 2.0, 3.0, 4.0]
  let m1 = @src.DensityMap::simple(
    title="a", nx=2, ny=2, nz=1, xlen=10.0, ylen=10.0, zlen=10.0, data=d1,
  )
  let m2 = @src.DensityMap::simple(
    title="b", nx=2, ny=2, nz=1, xlen=10.0, ylen=10.0, zlen=10.0, data=d2,
  )
  let m3 = m1.subtract(m2)
  assert_eq(m3.get(0, 0, 0), 9.0)
  assert_eq(m3.get(1, 1, 0), 36.0)
}

test "dm_scale" {
  let data = [1.0, 2.0, 3.0, 4.0]
  let m = @src.DensityMap::simple(
    title="t", nx=2, ny=2, nz=1, xlen=10.0, ylen=10.0, zlen=10.0, data=data,
  )
  let m2 = m.scale(2.0)
  assert_eq(m2.get(0, 0, 0), 2.0)
  assert_eq(m2.get(1, 1, 0), 8.0)
}

test "dm_correlation_identical" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m1 = @src.DensityMap::simple(
    title="a", nx=2, ny=2, nz=2, xlen=10.0, ylen=10.0, zlen=10.0, data=data,
  )
  let m2 = @src.DensityMap::simple(
    title="b", nx=2, ny=2, nz=2, xlen=10.0, ylen=10.0, zlen=10.0, data=data,
  )
  let cc = m1.correlation(m2)
  assert_true((cc - 1.0).abs() < 1.0e-10)
}

test "dm_correlation_anticorrelated" {
  let d1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let d2 = [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
  let m1 = @src.DensityMap::simple(
    title="a", nx=2, ny=2, nz=2, xlen=10.0, ylen=10.0, zlen=10.0, data=d1,
  )
  let m2 = @src.DensityMap::simple(
    title="b", nx=2, ny=2, nz=2, xlen=10.0, ylen=10.0, zlen=10.0, data=d2,
  )
  let cc = m1.correlation(m2)
  assert_true((cc - (-1.0)).abs() < 1.0e-10)
}

test "dm_correlation_orthogonal" {
  let d1 = [1.0, -1.0, 1.0, -1.0]
  let d2 = [1.0, 1.0, -1.0, -1.0]
  let m1 = @src.DensityMap::simple(
    title="a", nx=2, ny=2, nz=1, xlen=10.0, ylen=10.0, zlen=10.0, data=d1,
  )
  let m2 = @src.DensityMap::simple(
    title="b", nx=2, ny=2, nz=1, xlen=10.0, ylen=10.0, zlen=10.0, data=d2,
  )
  let cc = m1.correlation(m2)
  assert_true(cc.abs() < 1.0e-10)
}

test "dm_mask" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m = @src.DensityMap::simple(
    title="t", nx=2, ny=2, nz=2, xlen=10.0, ylen=10.0, zlen=10.0, data=data,
  )
  let masked = m.mask(4.5)
  assert_eq(masked.get(0, 0, 0), 0.0)
  assert_eq(masked.get(1, 1, 1), 1.0)
  assert_eq(masked.count_above(0.5), 4)
}

test "dm_voxel_volume" {
  let m = @src.DensityMap::simple(
    title="t", nx=2, ny=2, nz=2, xlen=10.0, ylen=10.0, zlen=10.0,
    data=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
  )
  let vv = m.voxel_volume()
  assert_true((vv - 125.0).abs() < 1.0e-10)
}

test "dm_volume_above" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let m = @src.DensityMap::simple(
    title="t", nx=2, ny=2, nz=2, xlen=10.0, ylen=10.0, zlen=10.0, data=data,
  )
  let vol = m.volume_above(4.5)
  assert_true((vol - 500.0).abs() < 1.0e-10)
}

test "dm_grid_to_cart" {
  let m = @src.DensityMap::simple(
    title="t", nx=2, ny=2, nz=2, xlen=10.0, ylen=20.0, zlen=30.0,
    data=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
  )
  let (x, _, _) = m.grid_to_cart(1, 0, 0)
  assert_true((x - 5.0).abs() < 1.0e-10)
  let (_, y2, _) = m.grid_to_cart(0, 1, 0)
  assert_true((y2 - 10.0).abs() < 1.0e-10)
  let (_, _, z3) = m.grid_to_cart(0, 0, 1)
  assert_true((z3 - 15.0).abs() < 1.0e-10)
}

test "dm_frac_cart_roundtrip" {
  let m = @src.DensityMap::simple(
    title="t", nx=4, ny=4, nz=4, xlen=30.0, ylen=40.0, zlen=50.0,
    data=Array::make(64, 0.0),
  )
  let (u0, v0, w0) = (0.3, 0.7, 0.5)
  let (x, y, z) = m.frac_to_cart(u0, v0, w0)
  let (u, v, w) = m.cart_to_frac(x, y, z)
  assert_true((u - u0).abs() < 1.0e-8)
  assert_true((v - v0).abs() < 1.0e-8)
  assert_true((w - w0).abs() < 1.0e-8)
}
