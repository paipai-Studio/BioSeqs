///|
/// xcms-style chromatographic peak detection for LC-MS metabolomics.
///
/// Inspired by the R/Bioconductor package xcms (Smith et al. 2006), this module
/// provides simplified centWave and matchedFilter peak detection for
/// chromatographic traces (retention time vs intensity).
///
/// Reference:
///   Tautenhahn et al. (2008) "Highly sensitive feature detection for high
///   resolution LC/MS" BMC Bioinformatics 9:504.

///|
/// A single chromatographic trace (e.g., extracted ion chromatogram).
pub struct XcpChromatogram {
  rt : Array[Double]
  intensity : Array[Double]
  label : String
}

///|
/// A detected chromatographic peak.
pub(all) struct XcpPeak {
  mz : Double
  rt : Double
  rtmin : Double
  rtmax : Double
  intensity : Double
  area : Double
  snr : Double
  fwhm : Double
  peak_method : String
}

///|
pub fn xcp_chrom_new(
  rt : Array[Double],
  intensity : Array[Double],
  label? : String = "",
) -> XcpChromatogram {
  XcpChromatogram::{
    rt: rt,
    intensity: intensity,
    label: label,
  }
}

///|
pub fn xcp_baseline_rolling_min(
  intensity : Array[Double],
  window_size : Int,
) -> Array[Double] {
  let n = intensity.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let half = window_size / 2
  for i in 0..<n {
    let lo = if i - half >= 0 { i - half } else { 0 }
    let hi = if i + half < n { i + half } else { n - 1 }
    let mut mn = intensity[lo]
    let j_start = lo + 1
    for j in j_start..=hi {
      if intensity[j] < mn {
        mn = intensity[j]
      }
    }
    result[i] = mn
  }
  result
}

///|
pub fn xcp_baseline_rolling_median(
  intensity : Array[Double],
  window_size : Int,
) -> Array[Double] {
  let n = intensity.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let half = window_size / 2
  for i in 0..<n {
    let lo = if i - half >= 0 { i - half } else { 0 }
    let hi = if i + half < n { i + half } else { n - 1 }
    let win_len = hi - lo + 1
    let win : Array[Double] = Array::make(win_len, 0.0)
    for k in 0..<win_len {
      win[k] = intensity[lo + k]
    }
    win.sort()
    result[i] = win[win_len / 2]
  }
  result
}

///|
pub fn xcp_noise_mad(intensity : Array[Double]) -> Double {
  let n = intensity.length()
  if n == 0 {
    return 0.0
  }
  let sorted = Array::make(n, 0.0)
  for i in 0..<n {
    sorted[i] = intensity[i]
  }
  sorted.sort()
  let med = sorted[n / 2]
  let devs : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    devs[i] = (intensity[i] - med).abs()
  }
  devs.sort()
  let mad = devs[n / 2]
  1.4826 * mad
}

///|
pub fn xcp_snr(intensity : Array[Double]) -> Double {
  let noise = xcp_noise_mad(intensity)
  if noise < 0.0000000001 {
    return 0.0
  }
  let mut mx = 0.0
  for v in intensity {
    if v > mx {
      mx = v
    }
  }
  mx / noise
}

///|
pub fn xcp_gaussian(x : Double, mu : Double, sigma : Double) -> Double {
  let z = (x - mu) / sigma
  @math.exp(-0.5 * z * z) / (sigma * 2.50662827463)
}

///|
fn xcp_fit_gaussian(
  rt : Array[Double],
  intensity : Array[Double],
  lo : Int,
  hi : Int,
) -> (Double, Double, Double, Double) {
  let n = hi - lo + 1
  if n < 3 {
    return (0.0, 0.0, 0.0, 0.0)
  }
  let mut sum_int = 0.0
  let mut sum_rt_int = 0.0
  let mut max_int = 0.0
  for i in lo..=hi {
    sum_int = sum_int + intensity[i]
    sum_rt_int = sum_rt_int + rt[i] * intensity[i]
    if intensity[i] > max_int {
      max_int = intensity[i]
    }
  }
  if sum_int < 0.000000000001 {
    return (0.0, 0.0, 0.0, 0.0)
  }
  let mu = sum_rt_int / sum_int
  let mut variance = 0.0
  for i in lo..=hi {
    let d = rt[i] - mu
    variance = variance + intensity[i] * d * d
  }
  variance = variance / sum_int
  let sigma = if variance > 0.0 { variance.sqrt() } else { (rt[hi] - rt[lo]) / 4.0 }
  let mut ss_res = 0.0
  let mut ss_tot = 0.0
  let mean_int = sum_int / n.to_double()
  for i in lo..=hi {
    let fitted = max_int * xcp_gaussian(rt[i], mu, sigma)
    let resid = intensity[i] - fitted
    ss_res = ss_res + resid * resid
    let dev = intensity[i] - mean_int
    ss_tot = ss_tot + dev * dev
  }
  let r2 = if ss_tot > 0.000000000001 { 1.0 - ss_res / ss_tot } else { 0.0 }
  (mu, sigma, max_int, r2)
}

///|
pub fn xcp_centwave(
  chrom : XcpChromatogram,
  snr_threshold? : Double = 3.0,
  min_peak_width? : Double = 0.0,
  max_peak_width? : Double = 1000.0,
  noise? : Double = 0.0,
) -> Array[XcpPeak] {
  let rt = chrom.rt
  let intensity = chrom.intensity
  let n = rt.length()
  let peaks : Array[XcpPeak] = Array::new()
  if n < 5 {
    return peaks
  }
  let noise_est = if noise > 0.0 { noise } else { xcp_noise_mad(intensity) }
  let threshold = noise_est * snr_threshold

  let mut i = 0
  while i < n {
    if intensity[i] > threshold {
      let start = i
      while i < n && intensity[i] > threshold {
        i = i + 1
      }
      let end = i - 1
      let region_len = end - start + 1
      if region_len >= 3 {

        let mut apex_int = intensity[start]
        for k in start..=end {
          if intensity[k] > apex_int {
            apex_int = intensity[k]

          }
        }
          let (mu, sigma, _amp, _r2) = xcp_fit_gaussian(rt, intensity, start, end)
        let fwhm = 2.3554820045 * sigma
        let rt_min = mu - 1.177741 * sigma
        let rt_max = mu + 1.177741 * sigma
        let rt_min_clamped = if rt_min >= rt[0] { rt_min } else { rt[0] }
        let rt_max_clamped = if rt_max <= rt[n - 1] { rt_max } else { rt[n - 1] }
        let width = rt_max_clamped - rt_min_clamped
        if width >= min_peak_width && width <= max_peak_width {
          let mut area = 0.0
          for k in start..<end {
            let dt = rt[k + 1] - rt[k]
            area = area + 0.5 * (intensity[k] + intensity[k + 1]) * dt
          }
          let local_noise = xcp_noise_mad(intensity)
          let snr_val = if local_noise > 0.000000000001 { apex_int / local_noise } else { 0.0 }
          if snr_val >= snr_threshold {
            peaks.push(XcpPeak::{
              mz: 0.0,
              rt: mu,
              rtmin: rt_min_clamped,
              rtmax: rt_max_clamped,
              intensity: apex_int,
              area: area,
              snr: snr_val,
              fwhm: fwhm,
              peak_method: "centWave",
            })
          }
        }
      }
    } else {
      i = i + 1
    }
  }

  peaks.sort_by(fn(a, b) { a.rt.compare(b.rt) })
  peaks
}

///|
pub fn xcp_matched_filter(
  chrom : XcpChromatogram,
  window_size? : Int = 5,
  snr_threshold? : Double = 3.0,
  noise? : Double = 0.0,
) -> Array[XcpPeak] {
  let rt = chrom.rt
  let intensity = chrom.intensity
  let n = rt.length()
  let peaks : Array[XcpPeak] = Array::new()
  if n < 5 {
    return peaks
  }
  let noise_est = if noise > 0.0 { noise } else { xcp_noise_mad(intensity) }
  let threshold = noise_est * snr_threshold
  let half = window_size / 2

  let smoothed : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let lo = if i - half >= 0 { i - half } else { 0 }
    let hi = if i + half < n { i + half } else { n - 1 }
    let mut sum = 0.0
    let mut cnt = 0
    for j in lo..=hi {
      sum = sum + intensity[j]
      cnt = cnt + 1
    }
    smoothed[i] = sum / cnt.to_double()
  }

  let mut i = 1
  while i < n - 1 {
    if smoothed[i] > threshold && smoothed[i] > smoothed[i - 1] && smoothed[i] >= smoothed[i + 1] {
      let apex_idx = i
      let apex_int = intensity[apex_idx]
      let mut left = apex_idx
      while left > 0 && intensity[left - 1] > threshold {
        left = left - 1
        if left > 0 && intensity[left] < intensity[left - 1] && intensity[left] < intensity[left + 1] {
          break
        }
      }
      let mut right = apex_idx
      while right < n - 1 && intensity[right + 1] > threshold {
        right = right + 1
        if right < n - 1 && intensity[right] < intensity[right - 1] && intensity[right] < intensity[right + 1] {
          break
        }
      }
      let mut area = 0.0
      for k in left..<right {
        let dt = rt[k + 1] - rt[k]
        area = area + 0.5 * (intensity[k] + intensity[k + 1]) * dt
      }
      let snr_val = if noise_est > 0.000000000001 { apex_int / noise_est } else { 0.0 }
      if snr_val >= snr_threshold {
        let fwhm = rt[right] - rt[left]
        peaks.push(XcpPeak::{
          mz: 0.0,
          rt: rt[apex_idx],
          rtmin: rt[left],
          rtmax: rt[right],
          intensity: apex_int,
          area: area,
          snr: snr_val,
          fwhm: fwhm,
          peak_method: "matchedFilter",
        })
      }
      i = right + 1
    } else {
      i = i + 1
    }
  }

  peaks.sort_by(fn(a, b) { a.rt.compare(b.rt) })
  peaks
}

///|
pub fn xcp_merge_peaks(
  peaks : Array[XcpPeak],
  min_overlap? : Double = 0.0,
) -> Array[XcpPeak] {
  let n = peaks.length()
  if n <= 1 {
    return peaks
  }
  let sorted = Array::make(n, peaks[0])
  for i in 0..<n {
    sorted[i] = peaks[i]
  }
  sorted.sort_by(fn(a, b) { a.rtmin.compare(b.rtmin) })

  let merged : Array[XcpPeak] = Array::new()
  let mut current = sorted[0]
  for i in 1..<n {
    let next = sorted[i]
    let overlap = current.rtmax - next.rtmin
    if overlap > min_overlap {
      let total_int = current.intensity + next.intensity
      let new_rt = (current.rt * current.intensity + next.rt * next.intensity) / total_int
      let new_apex = if current.intensity > next.intensity { current.intensity } else { next.intensity }
      current = XcpPeak::{
        mz: current.mz,
        rt: new_rt,
        rtmin: current.rtmin,
        rtmax: if next.rtmax > current.rtmax { next.rtmax } else { current.rtmax },
        intensity: new_apex,
        area: current.area + next.area,
        snr: if current.snr > next.snr { current.snr } else { next.snr },
        fwhm: current.rtmax - current.rtmin,
        peak_method: current.peak_method + "+" + next.peak_method,
      }
    } else {
      merged.push(current)
      current = next
    }
  }
  merged.push(current)
  merged
}

///|
pub fn xcp_export_peaks(peaks : Array[XcpPeak]) -> String {
  let sb = StringBuilder()
  sb.write_string("mz\trt\trtmin\trtmax\tintensity\tarea\tsnr\tfwhm\tmethod\n")
  for p in peaks {
    sb.write_string(p.mz.to_string())
    sb.write_string("\t")
    sb.write_string(p.rt.to_string())
    sb.write_string("\t")
    sb.write_string(p.rtmin.to_string())
    sb.write_string("\t")
    sb.write_string(p.rtmax.to_string())
    sb.write_string("\t")
    sb.write_string(p.intensity.to_string())
    sb.write_string("\t")
    sb.write_string(p.area.to_string())
    sb.write_string("\t")
    sb.write_string(p.snr.to_string())
    sb.write_string("\t")
    sb.write_string(p.fwhm.to_string())
    sb.write_string("\t")
    sb.write_string(p.peak_method)
    sb.write_string("\n")
  }
  sb.to_string()
}

///|
pub fn xcp_ascii_plot(
  chrom : XcpChromatogram,
  peaks : Array[XcpPeak],
  width? : Int = 60,
  height? : Int = 20,
) -> String {
  let n = chrom.intensity.length()
  if n == 0 {
    return ""
  }
  let mut mn = chrom.intensity[0]
  let mut mx = chrom.intensity[0]
  for v in chrom.intensity {
    if v < mn { mn = v }
    if v > mx { mx = v }
  }
  if mx - mn < 0.000000000001 {
    return ""
  }
  let plot : Array[Array[Char]] = Array::make(height, Array::make(width, ' '))
  for y in 0..<height {
    for x in 0..<width {
      plot[y][x] = ' '
    }
  }
  for x in 0..<width {
    let idx = (x * n / width).min(n - 1)
    let val = (chrom.intensity[idx] - mn) / (mx - mn)
    let y = (val * (height - 1).to_double()).to_int().min(height - 1)
    let row = height - 1 - y
    plot[row][x] = '*'
  }
  for p in peaks {
    let mut best_idx = 0
    let mut best_diff = (chrom.rt[0] - p.rt).abs()
    for k in 0..<n {
      let diff = (chrom.rt[k] - p.rt).abs()
      if diff < best_diff {
        best_diff = diff
        best_idx = k
      }
    }
    let x = (best_idx * width / n).min(width - 1)
    let val = (chrom.intensity[best_idx] - mn) / (mx - mn)
    let y = (val * (height - 1).to_double()).to_int().min(height - 1)
    let row = height - 1 - y
    plot[row][x] = 'P'
  }
  let mut result = ""
  for y in 0..<height {
    let mut line = "|"
    for x in 0..<width {
      line = line + plot[y][x].to_string()
    }
    result = result + line + "|\n"
  }
  result = result + "+" + "-".repeat(width) + "+\n"
  result
}

///|
pub fn xcp_tic_stats(intensity : Array[Double]) -> (Double, Double, Double) {
  let mut total = 0.0
  let mut mx = 0.0
  let n = intensity.length()
  for v in intensity {
    total = total + v
    if v > mx { mx = v }
  }
  let mean = if n > 0 { total / n.to_double() } else { 0.0 }
  (total, mx, mean)
}
