///|
/// Tests for colorspiral module (Bio.Graphics.ColorSpiral).

test "spiral_hsv_to_rgb_red" {
  let c = @src.spiral_hsv_to_rgb(0.0, 1.0, 1.0)
  assert_true(c.r > 200)
  assert_true(c.g < 100)
}

test "spiral_hsv_to_rgb_green" {
  let c = @src.spiral_hsv_to_rgb(0.333333, 1.0, 1.0)
  assert_true(c.g > 200)
}

test "spiral_hsv_to_rgb_blue" {
  let c = @src.spiral_hsv_to_rgb(0.666667, 1.0, 1.0)
  assert_true(c.b > 200)
}

test "spiral_color_at" {
  let sp = @src.spiral_new()
  let c0 = @src.spiral_color_at(sp, 0)
  let c1 = @src.spiral_color_at(sp, 1)
  assert_true(c0.g > 200)
  assert_true(c1.r >= 0)
}

test "spiral_colors" {
  let sp = @src.spiral_new()
  let cs = @src.spiral_colors(sp, 5)
  assert_eq(cs.length(), 5)
  for c in cs {
    assert_true(c.r >= 0 && c.r <= 255)
    assert_true(c.g >= 0 && c.g <= 255)
    assert_true(c.b >= 0 && c.b <= 255)
  }
}

test "spiral_hex" {
  let c = @src.SpiralColor::new(255, 0, 0)
  let hex = @src.spiral_hex(c)
  assert_eq(hex, "#FF0000")
}

test "spiral_rgb" {
  let c = @src.SpiralColor::new(0, 255, 0)
  let rgb = @src.spiral_rgb(c)
  assert_eq(rgb, "rgb(0,255,0)")
}

test "spiral_luminance" {
  let white = @src.SpiralColor::new(255, 255, 255)
  let black = @src.SpiralColor::new(0, 0, 0)
  assert_true(@src.spiral_luminance(white) > @src.spiral_luminance(black))
}

test "spiral_contrast" {
  let white = @src.SpiralColor::new(255, 255, 255)
  let dark = @src.SpiralColor::new(30, 30, 30)
  let wc = @src.spiral_contrast(white)
  let dc = @src.spiral_contrast(dark)
  assert_eq(wc.r, 0)
  assert_eq(dc.r, 255)
}

test "spiral_svg_palette" {
  let cs = [@src.SpiralColor::new(255,0,0), @src.SpiralColor::new(0,255,0)]
  let svg = @src.spiral_svg_palette(cs)
  assert_true(svg.contains("<svg"))
  assert_true(svg.contains("#FF0000"))
}

test "spiral_scale_color" {
  let c = @src.spiral_scale_color(0.5, 0.0, 1.0)
  assert_true(c.r >= 0 && c.r <= 255)
}

test "spiral_new_with" {
  let sp = @src.spiral_new_with(0.0, 0.01, 0.8, 0.8, 0.02)
  let c = @src.spiral_color_at(sp, 0)
  assert_true(c.r >= 0)
}

test "spiral_hsv_at" {
  let sp = @src.spiral_new()
  let (h, s, v) = @src.spiral_hsv_at(sp, 0)
  assert_true(h >= 0.0)
  assert_true(s > 0.0)
}
