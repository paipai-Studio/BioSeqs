///|
/// Tests for ebimage module (Image processing).

test "ebi_image_new" {
  let img = @src.EbiImage::new(10, 8)
  assert_eq(img.width, 10)
  assert_eq(img.height, 8)
  assert_eq(img.data.length(), 80)
}

test "ebi_image_from_rows" {
  let rows = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let img = @src.EbiImage::from_rows(rows)
  assert_eq(img.width, 3)
  assert_eq(img.height, 2)
  assert_true((img.data[0] - 1.0).abs() < 1.0e-9)
  assert_true((img.data[5] - 6.0).abs() < 1.0e-9)
}

test "ebi_image_to_rows_roundtrip" {
  let rows = [[1.0, 2.0], [3.0, 4.0]]
  let img = @src.EbiImage::from_rows(rows)
  let out = @src.EbiImage::to_rows(img)
  assert_eq(out.length(), 2)
  assert_true((out[0][0] - 1.0).abs() < 1.0e-9)
  assert_true((out[1][1] - 4.0).abs() < 1.0e-9)
}

test "ebi_image_get_set" {
  let img = @src.EbiImage::new(5, 5)
  @src.EbiImage::set(img, 2, 3, 7.5)
  let v = @src.EbiImage::get(img, 2, 3)
  assert_true((v - 7.5).abs() < 1.0e-9)
  assert_true(@src.EbiImage::get(img, 10, 10) == 0.0)
  @src.EbiImage::set(img, -1, 0, 9.0)
  @src.EbiImage::set(img, 0, -1, 9.0)
  assert_true(@src.EbiImage::get(img, 0, 0) == 0.0)
}

test "ebi_min_max" {
  let img = @src.ebi_example_image()
  let mn = @src.ebi_min(img)
  let mx = @src.ebi_max(img)
  assert_true(mn < mx)
  assert_true(mn >= 0.0)
  assert_true(mx <= 1.0)
}

test "ebi_mean_sd" {
  let img = @src.ebi_example_image()
  let m = @src.ebi_mean(img)
  let s = @src.ebi_sd(img)
  assert_true(m > 0.0 && m < 1.0)
  assert_true(s > 0.0)
}

test "ebi_normalize" {
  let img = @src.ebi_example_image()
  let norm = @src.ebi_normalize(img)
  let mn = @src.ebi_min(norm)
  let mx = @src.ebi_max(norm)
  assert_true(mn.abs() < 1.0e-9)
  assert_true((mx - 1.0).abs() < 1.0e-9)
}

test "ebi_otsu_threshold" {
  let img = @src.ebi_example_image()
  let t = @src.ebi_otsu_threshold(img)
  // The example image is two-level; Otsu must pick a split boundary in
  // [min, max] and it must separate the two levels.
  assert_true(t >= @src.ebi_min(img) && t <= @src.ebi_max(img))
  let binary = @src.ebi_threshold(img, t)
  let mut ones = 0
  for v in binary.data {
    if v == 1.0 {
      ones = ones + 1
    }
  }
  assert_true(ones > 0 && ones < 256)
}

test "ebi_threshold" {
  let img = @src.ebi_example_image()
  let binary = @src.ebi_threshold(img, 0.5)
  assert_eq(binary.width, img.width)
  assert_eq(binary.height, img.height)
  let mx = @src.ebi_max(binary)
  assert_true(mx <= 1.0)
  let mn = @src.ebi_min(binary)
  assert_true(mn >= 0.0)
}

test "ebi_threshold_otsu" {
  let img = @src.ebi_example_image()
  let binary = @src.ebi_threshold_otsu(img)
  let mut ones = 0
  for v in binary.data {
    if v == 1.0 {
      ones = ones + 1
    }
  }
  assert_true(ones > 0)
}

test "ebi_gaussian_blur" {
  let img = @src.ebi_example_image()
  let blurred = @src.ebi_gaussian_blur(img, 1.0)
  assert_eq(blurred.width, img.width)
  assert_eq(blurred.height, img.height)
  let mx = @src.ebi_max(blurred)
  assert_true(mx > 0.0 && mx <= @src.ebi_max(img) + 1.0e-9)
}

test "ebi_sobel_gradient" {
  let img = @src.ebi_example_image()
  let grad = @src.ebi_sobel_gradient(img)
  assert_eq(grad.width, img.width)
  assert_eq(grad.height, img.height)
  let mx = @src.ebi_max(grad)
  assert_true(mx > 0.0)
}

test "ebi_se_sizes" {
  assert_eq(@src.ebi_se_cross().length(), 4)
  assert_eq(@src.ebi_se_square().length(), 8)
  assert_true(@src.ebi_se_disk(1).length() >= 4)
}

test "ebi_erode_dilate_bounds" {
  let img = @src.ebi_example_image()
  let se = @src.ebi_se_cross()
  let eroded = @src.ebi_erode(img, se)
  let dilated = @src.ebi_dilate(img, se)
  assert_true(@src.ebi_min(eroded) <= @src.ebi_min(img))
  assert_true(@src.ebi_max(dilated) >= @src.ebi_max(img))
  assert_true(@src.ebi_max(eroded) <= @src.ebi_max(img))
  assert_true(@src.ebi_min(dilated) >= @src.ebi_min(img))
}

test "ebi_open_close" {
  let img = @src.ebi_example_image()
  let se = @src.ebi_se_cross()
  let opened = @src.ebi_open(img, se)
  let closed = @src.ebi_close(img, se)
  assert_eq(opened.width, img.width)
  assert_eq(closed.height, img.height)
}

test "ebi_label_two_objects" {
  let img = @src.ebi_example_image()
  let binary = @src.ebi_threshold_otsu(img)
  let (labels, count) = @src.ebi_label(binary, connectivity=8)
  assert_true(count >= 1 && count <= 2)
  assert_eq(labels.length(), img.width * img.height)
  for l in labels {
    assert_true(l >= 0 && l <= count)
  }
}

test "ebi_label_4_vs_8_connectivity" {
  let rows = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]]
  let img = @src.EbiImage::from_rows(rows)
  let (_, count8) = @src.ebi_label(img, connectivity=8)
  let (_, count4) = @src.ebi_label(img, connectivity=4)
  assert_eq(count8, 1)
  assert_eq(count4, 2)
}

test "ebi_object_features" {
  let img = @src.ebi_example_image()
  let binary = @src.ebi_threshold_otsu(img)
  let (labels, count) = @src.ebi_label(binary, connectivity=8)
  let features = @src.ebi_object_features(labels, img.width, img.height, count)
  assert_eq(features.length(), count)
  for obj in features {
    assert_true(obj.area > 0)
    assert_true(obj.cx >= 0.0 && obj.cx < img.width.to_double())
    assert_true(obj.cy >= 0.0 && obj.cy < img.height.to_double())
    assert_true(obj.max_x >= obj.min_x)
    assert_true(obj.max_y >= obj.min_y)
  }
}

test "ebi_watershed" {
  let img = @src.ebi_example_image()
  let binary = @src.ebi_threshold_otsu(img)
  let (labels, count) = @src.ebi_label(binary, connectivity=4)
  let ws = @src.ebi_watershed(img, labels)
  assert_eq(ws.length(), img.width * img.height)
  for l in ws {
    assert_true(l >= 0 && l <= count)
  }
}

test "ebi_to_ascii" {
  let img = @src.ebi_example_image()
  let ascii = @src.ebi_to_ascii(img)
  assert_true(ascii.length() > 0)
  assert_true(ascii.contains("\n"))
}

test "ebi_example_image" {
  let img = @src.ebi_example_image()
  assert_eq(img.width, 16)
  assert_eq(img.height, 16)
  assert_eq(img.data.length(), 256)
  assert_true(@src.ebi_max(img) <= 1.0)
}
