///|
/// Demo for colorspiral module: HSV color spiral for chart series.

fn main {
  println("=== colorspiral Demo: HSV Color Spiral ===")
  println("")

  // Example 1: Default color spiral
  println("--- Step 1: Default Color Spiral ---")
  let sp = @src.spiral_new()
  let cs = @src.spiral_colors(sp, 8)
  println("8 colors from default spiral:")
  for i in 0..<cs.length() {
    let hex = @src.spiral_hex(cs[i])
    let rgb = @src.spiral_rgb(cs[i])
    println("  " + i.to_string() + ": " + hex + " " + rgb)
  }
  println("")

  // Example 2: HSV to RGB conversion
  println("--- Step 2: HSV to RGB ---")
  let red = @src.spiral_hsv_to_rgb(0.0, 1.0, 1.0)
  let green = @src.spiral_hsv_to_rgb(0.333333, 1.0, 1.0)
  let blue = @src.spiral_hsv_to_rgb(0.666667, 1.0, 1.0)
  println("Red:   " + @src.spiral_hex(red))
  println("Green: " + @src.spiral_hex(green))
  println("Blue:  " + @src.spiral_hex(blue))
  println("")

  // Example 3: Custom spiral
  println("--- Step 3: Custom Spiral ---")
  let sp2 = @src.spiral_new_with(0.0, 0.05, 0.9, 0.9, 0.03)
  let cs2 = @src.spiral_colors(sp2, 6)
  println("6 colors from custom spiral (red start):")
  for i in 0..<cs2.length() {
    println("  " + i.to_string() + ": " + @src.spiral_hex(cs2[i]))
  }
  println("")

  // Example 4: HSV values
  println("--- Step 4: HSV Values ---")
  for i in 0..<5 {
    let (h, s, v) = @src.spiral_hsv_at(sp, i)
    println("  color " + i.to_string() + ": h=" + h.to_string() + " s=" + s.to_string() + " v=" + v.to_string())
  }
  println("")

  // Example 5: Luminance & contrast
  println("--- Step 5: Luminance & Contrast ---")
  let c = @src.spiral_hsv_to_rgb(0.5, 0.8, 0.5)
  let lum = @src.spiral_luminance(c)
  let contrast = @src.spiral_contrast(c)
  println("Color: " + @src.spiral_hex(c))
  println("Luminance: " + lum.to_string())
  println("Contrast (text) color: " + @src.spiral_hex(contrast))
  println("")

  // Example 6: SVG palette
  println("--- Step 6: SVG Palette ---")
  let palette = @src.spiral_svg_palette(cs, labels=["A", "B", "C", "D", "E", "F", "G", "H"])
  println("SVG palette length: " + palette.length().to_string() + " chars")
  println("")

  // Example 7: Value-to-color scale
  println("--- Step 7: Value-to-Color Scale ---")
  let vals = [-2.0, -1.0, 0.0, 1.0, 2.0]
  for v in vals {
    let c = @src.spiral_scale_color(v, -2.0, 2.0)
    println("  value " + v.to_string() + " -> " + @src.spiral_hex(c))
  }
  println("")

  println("=== colorspiral Demo Complete ===")
}
