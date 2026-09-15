///|
/// ebimage - Image processing for cellular / tissue imaging.
/// Inspired by Bioconductor {EBImage} (Pau G. et al., 2010).
///
/// Provides a compact grayscale image type together with the core operations
/// used in high-content screening pipelines: statistics and contrast
/// normalisation, Otsu thresholding, Gaussian smoothing, Sobel edge detection,
/// morphological filters, connected-component labelling, object feature
/// extraction and marker-based watershed segmentation.

///|
/// A single-channel (grayscale) image stored row-major with `data[y*width+x]`.
pub struct EbiImage {
  pub width : Int
  pub height : Int
  pub data : Array[Double]
}

///|
/// Construct a zero-filled image.
pub fn EbiImage::new(width : Int, height : Int) -> EbiImage {
  EbiImage::{ width, height, data: Array::make(width * height, 0.0) }
}

///|
/// Build an image from a list of rows (each row is a scanline).
pub fn EbiImage::from_rows(rows : Array[Array[Double]]) -> EbiImage {
  let h = rows.length()
  if h == 0 {
    return EbiImage::new(0, 0)
  }
  let w = rows[0].length()
  let data : Array[Double] = []
  for r in rows {
    for v in r {
      data.push(v)
    }
  }
  EbiImage::{ width: w, height: h, data }
}

///|
/// Row-major copy of the pixel data.
pub fn EbiImage::to_rows(self : EbiImage) -> Array[Array[Double]] {
  let rows : Array[Array[Double]] = []
  let mut y = 0
  while y < self.height {
    let row : Array[Double] = []
    let mut x = 0
    while x < self.width {
      row.push(self.data[y * self.width + x])
      x = x + 1
    }
    rows.push(row)
    y = y + 1
  }
  rows
}

///|
/// Read a pixel (0.0 outside the image).
pub fn EbiImage::get(self : EbiImage, x : Int, y : Int) -> Double {
  if x < 0 || x >= self.width || y < 0 || y >= self.height {
    0.0
  } else {
    self.data[y * self.width + x]
  }
}

///|
/// Write a pixel (out-of-range coordinates are ignored).
pub fn EbiImage::set(self : EbiImage, x : Int, y : Int, value : Double) -> Unit {
  if x >= 0 && x < self.width && y >= 0 && y < self.height {
    self.data[y * self.width + x] = value
  }
}

// ---- statistics -----------------------------------------------------------

///|
/// Minimum pixel value (0.0 for an empty image).
pub fn ebi_min(img : EbiImage) -> Double {
  if img.data.length() == 0 {
    return 0.0
  }
  let mut m = img.data[0]
  for v in img.data {
    if v < m {
      m = v
    }
  }
  m
}

///|
/// Maximum pixel value (0.0 for an empty image).
pub fn ebi_max(img : EbiImage) -> Double {
  if img.data.length() == 0 {
    return 0.0
  }
  let mut m = img.data[0]
  for v in img.data {
    if v > m {
      m = v
    }
  }
  m
}

///|
/// Mean pixel intensity.
pub fn ebi_mean(img : EbiImage) -> Double {
  let n = img.data.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in img.data {
    s = s + v
  }
  s / n.to_double()
}

///|
/// Sample standard deviation of pixel intensities.
pub fn ebi_sd(img : EbiImage) -> Double {
  let n = img.data.length()
  if n == 0 {
    return 0.0
  }
  let m = ebi_mean(img)
  let mut s = 0.0
  for v in img.data {
    let d = v - m
    s = s + d * d
  }
  (s / n.to_double()).sqrt()
}

///|
/// Linearly rescale intensities into [0, 1] (constant images become all zeros).
pub fn ebi_normalize(img : EbiImage) -> EbiImage {
  let mn = ebi_min(img)
  let mx = ebi_max(img)
  let out = EbiImage::new(img.width, img.height)
  let range = mx - mn
  let n = img.data.length()
  let mut i = 0
  while i < n {
    out.data[i] = if range < 1.0e-12 { 0.0 } else { (img.data[i] - mn) / range }
    i = i + 1
  }
  out
}

// ---- thresholding ---------------------------------------------------------

///|
/// Otsu's threshold: the intensity that maximises between-class variance.
pub fn ebi_otsu_threshold(img : EbiImage) -> Double {
  let mn = ebi_min(img)
  let mx = ebi_max(img)
  if mx - mn < 1.0e-12 {
    return mn
  }
  let bins = 256
  let hist = Array::make(bins, 0.0)
  let scale = (bins - 1).to_double() / (mx - mn)
  for v in img.data {
    let mut b = ((v - mn) * scale).to_int()
    if b < 0 {
      b = 0
    }
    if b >= bins {
      b = bins - 1
    }
    hist[b] = hist[b] + 1.0
  }
  let total = img.data.length().to_double()
  let mut sum_total = 0.0
  let mut i = 0
  while i < bins {
    sum_total = sum_total + i.to_double() * hist[i]
    i = i + 1
  }
  let mut w_b = 0.0
  let mut sum_b = 0.0
  let mut max_var = -1.0
  let mut best = 0
  i = 0
  while i < bins {
    w_b = w_b + hist[i]
    if w_b > 0.0 {
      let w_f = total - w_b
      if w_f <= 0.0 {
        break
      }
      sum_b = sum_b + i.to_double() * hist[i]
      let m_b = sum_b / w_b
      let m_f = (sum_total - sum_b) / w_f
      let diff = m_b - m_f
      let var_between = w_b * w_f * diff * diff
      if var_between > max_var {
        max_var = var_between
        best = i
      }
    }
    i = i + 1
  }
  mn + best.to_double() / scale
}

///|
/// Binarise the image: pixels strictly above `threshold` become 1.0, else 0.0.
pub fn ebi_threshold(img : EbiImage, threshold : Double) -> EbiImage {
  let out = EbiImage::new(img.width, img.height)
  let n = img.data.length()
  let mut i = 0
  while i < n {
    out.data[i] = if img.data[i] > threshold { 1.0 } else { 0.0 }
    i = i + 1
  }
  out
}

///|
/// Binarise using an automatically computed Otsu threshold.
pub fn ebi_threshold_otsu(img : EbiImage) -> EbiImage {
  ebi_threshold(img, ebi_otsu_threshold(img))
}

// ---- smoothing and edges --------------------------------------------------

///|
/// Separable Gaussian blur with clamped borders.
pub fn ebi_gaussian_blur(img : EbiImage, sigma : Double) -> EbiImage {
  if img.width == 0 || img.height == 0 || sigma <= 0.0 {
    return EbiImage::{ width: img.width, height: img.height, data: img.data.copy() }
  }
  let radius = {
    let r = @math.ceil(3.0 * sigma).to_int()
    if r < 1 {
      1
    } else {
      r
    }
  }
  let ksize = 2 * radius + 1
  let kernel = Array::make(ksize, 0.0)
  let mut ksum = 0.0
  let mut i = -radius
  while i <= radius {
    let v = @math.exp(-(i.to_double() * i.to_double()) / (2.0 * sigma * sigma))
    kernel[i + radius] = v
    ksum = ksum + v
    i = i + 1
  }
  let mut k = 0
  while k < ksize {
    kernel[k] = kernel[k] / ksum
    k = k + 1
  }

  // horizontal pass
  let tmp = EbiImage::new(img.width, img.height)
  let mut y = 0
  while y < img.height {
    let mut x = 0
    while x < img.width {
      let mut acc = 0.0
      let mut t = -radius
      while t <= radius {
        let mut sx = x + t
        if sx < 0 {
          sx = 0
        }
        if sx >= img.width {
          sx = img.width - 1
        }
        acc = acc + img.data[y * img.width + sx] * kernel[t + radius]
        t = t + 1
      }
      tmp.data[y * img.width + x] = acc
      x = x + 1
    }
    y = y + 1
  }

  // vertical pass
  let out = EbiImage::new(img.width, img.height)
  y = 0
  while y < img.height {
    let mut x = 0
    while x < img.width {
      let mut acc = 0.0
      let mut t = -radius
      while t <= radius {
        let mut sy = y + t
        if sy < 0 {
          sy = 0
        }
        if sy >= img.height {
          sy = img.height - 1
        }
        acc = acc + tmp.data[sy * img.width + x] * kernel[t + radius]
        t = t + 1
      }
      out.data[y * img.width + x] = acc
      x = x + 1
    }
    y = y + 1
  }
  out
}

///|
/// Sobel gradient magnitude (edge strength).
pub fn ebi_sobel_gradient(img : EbiImage) -> EbiImage {
  let out = EbiImage::new(img.width, img.height)
  let mut y = 0
  while y < img.height {
    let mut x = 0
    while x < img.width {
      let tl = img.get(x - 1, y - 1)
      let tc = img.get(x, y - 1)
      let tr = img.get(x + 1, y - 1)
      let ml = img.get(x - 1, y)
      let mr = img.get(x + 1, y)
      let bl = img.get(x - 1, y + 1)
      let bc = img.get(x, y + 1)
      let br = img.get(x + 1, y + 1)
      let gx = (tr + 2.0 * mr + br) - (tl + 2.0 * ml + bl)
      let gy = (bl + 2.0 * bc + br) - (tl + 2.0 * tc + tr)
      out.data[y * img.width + x] = (gx * gx + gy * gy).sqrt()
      x = x + 1
    }
    y = y + 1
  }
  out
}

// ---- morphology -----------------------------------------------------------

///|
/// Cross-shaped structuring element (4-neighbourhood offsets).
pub fn ebi_se_cross() -> Array[(Int, Int)] {
  [(1, 0), (-1, 0), (0, 1), (0, -1)]
}

///|
/// Full 3x3 square structuring element (8-neighbourhood offsets).
pub fn ebi_se_square() -> Array[(Int, Int)] {
  [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
}

///|
/// Disk structuring element of the given radius.
pub fn ebi_se_disk(radius : Int) -> Array[(Int, Int)] {
  let offsets : Array[(Int, Int)] = []
  let r = if radius < 1 { 1 } else { radius }
  let mut dy = -r
  while dy <= r {
    let mut dx = -r
    while dx <= r {
      if dx * dx + dy * dy <= r * r && !(dx == 0 && dy == 0) {
        offsets.push((dx, dy))
      }
      dx = dx + 1
    }
    dy = dy + 1
  }
  offsets
}

///|
/// Grayscale erosion (local minimum) with the given structuring element.
pub fn ebi_erode(img : EbiImage, se : Array[(Int, Int)]) -> EbiImage {
  let out = EbiImage::new(img.width, img.height)
  let mut y = 0
  while y < img.height {
    let mut x = 0
    while x < img.width {
      let mut m = img.data[y * img.width + x]
      for off in se {
        let v = img.get(x + off.0, y + off.1)
        if v < m {
          m = v
        }
      }
      out.data[y * img.width + x] = m
      x = x + 1
    }
    y = y + 1
  }
  out
}

///|
/// Grayscale dilation (local maximum) with the given structuring element.
pub fn ebi_dilate(img : EbiImage, se : Array[(Int, Int)]) -> EbiImage {
  let out = EbiImage::new(img.width, img.height)
  let mut y = 0
  while y < img.height {
    let mut x = 0
    while x < img.width {
      let mut m = img.data[y * img.width + x]
      for off in se {
        let v = img.get(x + off.0, y + off.1)
        if v > m {
          m = v
        }
      }
      out.data[y * img.width + x] = m
      x = x + 1
    }
    y = y + 1
  }
  out
}

///|
/// Morphological opening (erosion followed by dilation).
pub fn ebi_open(img : EbiImage, se : Array[(Int, Int)]) -> EbiImage {
  ebi_dilate(ebi_erode(img, se), se)
}

///|
/// Morphological closing (dilation followed by erosion).
pub fn ebi_close(img : EbiImage, se : Array[(Int, Int)]) -> EbiImage {
  ebi_erode(ebi_dilate(img, se), se)
}

// ---- connected components -------------------------------------------------

///|
/// Label the non-zero (foreground) pixels of a binary image.
/// `connectivity` must be 4 or 8. Background stays label 0 and component ids
/// start at 1. Returns (labels, number of components).
pub fn ebi_label(img : EbiImage, connectivity~ : Int = 8) -> (Array[Int], Int) {
  let w = img.width
  let h = img.height
  let labels = Array::make(w * h, 0)
  let offsets = if connectivity == 4 { ebi_se_cross() } else { ebi_se_square() }
  let mut count = 0
  let mut idx = 0
  let total = w * h
  while idx < total {
    if img.data[idx] != 0.0 && labels[idx] == 0 {
      count = count + 1
      let stack : Array[Int] = [idx]
      labels[idx] = count
      while stack.length() > 0 {
        let cur = stack.pop().unwrap()
        let cx = cur % w
        let cy = cur / w
        for off in offsets {
          let nx = cx + off.0
          let ny = cy + off.1
          if nx >= 0 && nx < w && ny >= 0 && ny < h {
            let nidx = ny * w + nx
            if img.data[nidx] != 0.0 && labels[nidx] == 0 {
              labels[nidx] = count
              stack.push(nidx)
            }
          }
        }
      }
    }
    idx = idx + 1
  }
  (labels, count)
}

///|
/// Per-object morphology measurements.
pub struct EbiObject {
  pub label : Int
  pub area : Int
  pub cx : Double
  pub cy : Double
  pub min_x : Int
  pub max_x : Int
  pub min_y : Int
  pub max_y : Int
}

///|
/// Compute area, centroid and bounding box for every labelled object.
pub fn ebi_object_features(
  labels : Array[Int],
  width : Int,
  height : Int,
  count : Int,
) -> Array[EbiObject] {
  let areas = Array::make(count + 1, 0)
  let sum_x = Array::make(count + 1, 0)
  let sum_y = Array::make(count + 1, 0)
  let min_x = Array::make(count + 1, width)
  let max_x = Array::make(count + 1, -1)
  let min_y = Array::make(count + 1, height)
  let max_y = Array::make(count + 1, -1)

  let mut y = 0
  while y < height {
    let mut x = 0
    while x < width {
      let l = labels[y * width + x]
      if l > 0 && l <= count {
        areas[l] = areas[l] + 1
        sum_x[l] = sum_x[l] + x
        sum_y[l] = sum_y[l] + y
        if x < min_x[l] {
          min_x[l] = x
        }
        if x > max_x[l] {
          max_x[l] = x
        }
        if y < min_y[l] {
          min_y[l] = y
        }
        if y > max_y[l] {
          max_y[l] = y
        }
      }
      x = x + 1
    }
    y = y + 1
  }

  let out : Array[EbiObject] = []
  let mut l = 1
  while l <= count {
    let a = areas[l]
    let cx = if a > 0 { sum_x[l].to_double() / a.to_double() } else { 0.0 }
    let cy = if a > 0 { sum_y[l].to_double() / a.to_double() } else { 0.0 }
    out.push(
      EbiObject::{
        label: l,
        area: a,
        cx,
        cy,
        min_x: min_x[l],
        max_x: max_x[l],
        min_y: min_y[l],
        max_y: max_y[l],
      },
    )
    l = l + 1
  }
  out
}

// ---- marker-based watershed ----------------------------------------------

///|
/// Marker-based watershed segmentation (Meyer's flooding algorithm) on the
/// gradient magnitude of the image. `markers` is a row-major label image where
/// 0 marks unseeded pixels and positive values are seed labels. Returns the
/// row-major label image after flooding.
pub fn ebi_watershed(img : EbiImage, markers : Array[Int]) -> Array[Int] {
  let w = img.width
  let h = img.height
  let total = w * h
  let labels = Array::make(total, 0)
  if markers.length() != total {
    return labels
  }
  let grad = ebi_sobel_gradient(img)
  let heap : Array[(Double, Int)] = []

  let mut i = 0
  while i < total {
    if markers[i] > 0 {
      labels[i] = markers[i]
      heap.push((grad.data[i], i))
    }
    i = i + 1
  }

  let nbrs = ebi_se_cross()
  let mut guard_count = 0
  while heap.length() > 0 && guard_count < total * 8 {
    let popped = heap.pop()
    if popped.is_none() {
      break
    }
    let cur_data = popped.unwrap()
    let priority = cur_data.0
    let cur = cur_data.1
    let cx = cur % w
    let cy = cur / w
    for off in nbrs {
      let nx = cx + off.0
      let ny = cy + off.1
      if nx >= 0 && nx < w && ny >= 0 && ny < h {
        let nidx = ny * w + nx
        if labels[nidx] == 0 {
          labels[nidx] = labels[cur]
          let g = grad.data[nidx]
          let np = if g > priority { g } else { priority }
          heap.push((np, nidx))
        }
      }
    }
    guard_count = guard_count + 1
  }
  labels
}

///|
/// Render a small image as ASCII art (useful for demos and debugging).
pub fn ebi_to_ascii(img : EbiImage) -> String {
  let ramp = " .:-=+*#%@"
  let mn = ebi_min(img)
  let mx = ebi_max(img)
  let range = mx - mn
  let buf = StringBuilder()
  let mut y = 0
  while y < img.height {
    let mut x = 0
    while x < img.width {
      let v = img.data[y * img.width + x]
      let frac = if range < 1.0e-12 { 0.0 } else { (v - mn) / range }
      let mut idx = (frac * 9.0).to_int()
      if idx < 0 {
        idx = 0
      }
      if idx > 9 {
        idx = 9
      }
      buf.write_string(ramp[idx:idx + 1].to_owned())
      x = x + 1
    }
    buf.write_string("\n")
    y = y + 1
  }
  buf.to_string()
}

///|
/// Small synthetic image used by tests and the demo: two bright disks on a
/// dark background.
pub fn ebi_example_image() -> EbiImage {
  let size = 16
  let img = EbiImage::new(size, size)
  let mut y = 0
  while y < size {
    let mut x = 0
    while x < size {
      let dx = (x - 5).to_double()
      let dy = (y - 5).to_double()
      let inside_a = dx * dx + dy * dy <= 9.0
      let ex = (x - 11).to_double()
      let ey = (y - 11).to_double()
      let inside_b = ex * ex + ey * ey <= 4.0
      img.data[y * size + x] = if inside_a || inside_b { 1.0 } else { 0.05 }
      x = x + 1
    }
    y = y + 1
  }
  img
}
