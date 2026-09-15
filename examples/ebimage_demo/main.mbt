///|
/// ebimage demo - Grayscale image processing for cellular imaging.
/// Inspired by Bioconductor {EBImage}.

fn main {
  println("=== ebimage Demo: Image Processing Pipeline ===")
  println("")

  // Step 1: synthetic image with two bright "cells"
  println("--- Step 1: Synthetic Image (two disks) ---")
  let img = @src.ebi_example_image()
  println(@src.ebi_to_ascii(img))
  println(
    "  stats: min=" + @src.ebi_min(img).to_string() +
    " max=" + @src.ebi_max(img).to_string() +
    " mean=" + @src.ebi_mean(img).to_string() +
    " sd=" + @src.ebi_sd(img).to_string(),
  )
  println("")

  // Step 2: smoothing and edge detection
  println("--- Step 2: Gaussian Blur + Sobel Edges ---")
  let blurred = @src.ebi_gaussian_blur(img, 1.0)
  let grad = @src.ebi_sobel_gradient(blurred)
  println("  max gradient magnitude: " + @src.ebi_max(grad).to_string())
  println("")

  // Step 3: Otsu thresholding
  println("--- Step 3: Otsu Threshold ---")
  let t = @src.ebi_otsu_threshold(img)
  let binary = @src.ebi_threshold(img, t)
  println("  otsu threshold = " + t.to_string())
  println(@src.ebi_to_ascii(binary))
  println("")

  // Step 4: morphological cleanup
  println("--- Step 4: Morphological Opening ---")
  let opened = @src.ebi_open(binary, @src.ebi_se_cross())
  let mut fg = 0
  for v in opened.data {
    if v > 0.0 {
      fg = fg + 1
    }
  }
  let mut fg0 = 0
  for v in binary.data {
    if v > 0.0 {
      fg0 = fg0 + 1
    }
  }
  println("  foreground pixels: " + fg0.to_string() + " -> after opening: " + fg.to_string())
  println("")

  // Step 5: connected components
  println("--- Step 5: Connected Components (8-connectivity) ---")
  let (labels, count) = @src.ebi_label(opened, connectivity=8)
  println("  objects found: " + count.to_string())
  let features = @src.ebi_object_features(labels, img.width, img.height, count)
  for obj in features {
    println(
      "  object " + obj.label.to_string() + ": area=" + obj.area.to_string() +
      " centroid=(" + obj.cx.to_string() + ", " + obj.cy.to_string() + ")" +
      " bbox=[" + obj.min_x.to_string() + "," + obj.min_y.to_string() +
      "]-[" + obj.max_x.to_string() + "," + obj.max_y.to_string() + "]",
    )
  }
  println("")

  // Step 6: marker-based watershed
  println("--- Step 6: Marker-based Watershed ---")
  let ws = @src.ebi_watershed(img, labels)
  let mut max_label = 0
  for l in ws {
    if l > max_label {
      max_label = l
    }
  }
  println("  watershed regions: " + max_label.to_string())
  println("")

  // Step 7: contrast normalisation
  println("--- Step 7: Contrast Normalisation ---")
  let norm = @src.ebi_normalize(grad)
  println("  normalised gradient range: [" + @src.ebi_min(norm).to_string() + ", " + @src.ebi_max(norm).to_string() + "]")

  println("")
  println("=== ebimage Demo Complete ===")
}
