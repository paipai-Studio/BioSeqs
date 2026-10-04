///|
/// ProDy-style hinge site detection from slow normal modes.
///
/// Hinge residues are points where the protein bends — regions of minimal
/// displacement that separate domains moving in opposite directions.
///
/// Two detection strategies:
/// - **GNM (sign change)**: hinge between residues i, i+1 when mode[i]
///   and mode[i+1] have opposite signs
/// - **ANM (local minimum)**: per-residue displacement magnitude has a
///   local minimum at the hinge
///
/// Reference: Bahar, Atilgan, Erman (2001) "Direct evaluation of thermal
/// fluctuations in proteins using a single-parameter harmonic potential."

// ============================================================================
// Types
// ============================================================================

///|
/// A detected hinge site.
pub struct HingeSite {
  residue_idx : Int // index of the hinge residue
  mode_idx : Int // which mode detected this hinge (-1 = multi-mode)
  magnitude : Double // displacement magnitude at the hinge (small)
  flank_left : Double // displacement of left flanking region
  flank_right : Double // displacement of right flanking region
}

///|
/// Construct a HingeSite.
pub fn HingeSite::new(
  residue_idx : Int,
  mode_idx : Int,
  magnitude : Double,
  flank_left : Double,
  flank_right : Double,
) -> HingeSite {
  HingeSite::{ residue_idx, mode_idx, magnitude, flank_left, flank_right }
}

///|
pub fn HingeSite::get_residue_idx(self : HingeSite) -> Int {
  self.residue_idx
}

///|
pub fn HingeSite::get_mode_idx(self : HingeSite) -> Int {
  self.mode_idx
}

///|
pub fn HingeSite::get_magnitude(self : HingeSite) -> Double {
  self.magnitude
}

// ============================================================================
// GNM hinge detection — sign changes in scalar mode
// ============================================================================

///|
/// Detect hinge sites from a GNM mode (N scalar values per residue).
///
/// A hinge exists between residues i and i+1 when mode[i] and mode[i+1]
/// have opposite signs.  The hinge residue is the one with smaller |mode|.
pub fn detect_hinges_gnm(
  mode : Array[Double],
  mode_idx : Int,
) -> Array[HingeSite] {
  let n = mode.length()
  if n < 3 {
    return []
  }
  let hinges : Array[HingeSite] = Array(capacity=0)
  for i = 0; i < n - 1; i = i + 1 {
    let v0 = mode[i]
    let v1 = mode[i + 1]
    if v0 * v1 < 0.0 {
      let hinge_res = if v0.abs() <= v1.abs() { i } else { i + 1 }
      let mag = mode[hinge_res].abs()
      hinges.push(HingeSite::{ residue_idx: hinge_res, mode_idx, magnitude: mag, flank_left: v0.abs(), flank_right: v1.abs() })
    }
  }
  hinges
}

///|
/// Detect hinge sites from a GNM mode via local minima of |mode|.
///
/// A residue is a hinge if its |mode| is a local minimum and below
/// a threshold relative to the maximum.
pub fn detect_hinges_gnm_local_min(
  mode : Array[Double],
  mode_idx : Int,
  threshold? : Double = 0.1,
) -> Array[HingeSite] {
  let n = mode.length()
  if n < 3 {
    return []
  }
  let abs_mode = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    abs_mode[i] = mode[i].abs()
  }
  let mut max_val = 0.0
  for i = 0; i < n; i = i + 1 {
    if abs_mode[i] > max_val {
      max_val = abs_mode[i]
    }
  }
  if max_val < 1.0e-20 {
    return []
  }
  let thresh = threshold * max_val
  let hinges : Array[HingeSite] = Array(capacity=0)
  for i = 1; i < n - 1; i = i + 1 {
    if abs_mode[i] < abs_mode[i - 1] && abs_mode[i] < abs_mode[i + 1] && abs_mode[i] < thresh {
      hinges.push(HingeSite::{ residue_idx: i, mode_idx, magnitude: abs_mode[i], flank_left: abs_mode[i - 1], flank_right: abs_mode[i + 1] })
    }
  }
  hinges
}

// ============================================================================
// ANM hinge detection — local minima of per-residue displacement
// ============================================================================

///|
/// Compute per-residue displacement magnitude from an ANM mode (3N vector).
/// d[i] = sqrt(v[3i]² + v[3i+1]² + v[3i+2]²)
pub fn anm_mode_magnitudes(mode : Array[Double]) -> Array[Double] {
  let dim = mode.length()
  let n = dim / 3
  let mags = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let i3 = 3 * i
    mags[i] = (mode[i3] * mode[i3] + mode[i3 + 1] * mode[i3 + 1] + mode[i3 + 2] * mode[i3 + 2]).sqrt()
  }
  mags
}

///|
/// Detect hinge sites from an ANM mode (3N values).
///
/// Computes per-residue displacement magnitude and finds local minima
/// below a threshold relative to the maximum.
pub fn detect_hinges_anm(
  mode : Array[Double],
  mode_idx : Int,
  threshold? : Double = 0.3,
) -> Array[HingeSite] {
  let mags = anm_mode_magnitudes(mode)
  let n = mags.length()
  if n < 3 {
    return []
  }
  let mut max_val = 0.0
  for i = 0; i < n; i = i + 1 {
    if mags[i] > max_val {
      max_val = mags[i]
    }
  }
  if max_val < 1.0e-20 {
    return []
  }
  let thresh = threshold * max_val
  let hinges : Array[HingeSite] = Array(capacity=0)
  for i = 1; i < n - 1; i = i + 1 {
    if mags[i] < mags[i - 1] && mags[i] < mags[i + 1] && mags[i] < thresh {
      hinges.push(HingeSite::{ residue_idx: i, mode_idx, magnitude: mags[i], flank_left: mags[i - 1], flank_right: mags[i + 1] })
    }
  }
  hinges
}

// ============================================================================
// Multi-mode analysis
// ============================================================================

///|
/// Detect hinge sites common to multiple GNM slow modes.
///
/// A residue is a common hinge if it appears as a hinge in at least
/// `min_count` of the provided modes.
pub fn detect_common_hinges_gnm(
  modes : Array[Array[Double]],
  num_modes : Int,
  min_count? : Int = 2,
) -> Array[HingeSite] {
  let n_modes = if num_modes < modes.length() { num_modes } else { modes.length() }
  if n_modes == 0 {
    return []
  }
  let n = modes[0].length()
  let counts = Array::make(n, 0)
  let min_mag = Array::make(n, 1.0e30)
  for k = 0; k < n_modes; k = k + 1 {
    let hinges = detect_hinges_gnm(modes[k], k)
    for h in hinges {
      let idx = h.residue_idx
      counts[idx] = counts[idx] + 1
      if h.magnitude < min_mag[idx] {
        min_mag[idx] = h.magnitude
      }
    }
  }
  let result : Array[HingeSite] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if counts[i] >= min_count {
      result.push(HingeSite::{ residue_idx: i, mode_idx: -1, magnitude: min_mag[i], flank_left: counts[i].to_double(), flank_right: n_modes.to_double() })
    }
  }
  result
}

///|
/// Detect hinge sites common to multiple ANM slow modes.
pub fn detect_common_hinges_anm(
  modes : Array[Array[Double]],
  num_modes : Int,
  min_count? : Int = 2,
  threshold? : Double = 0.3,
) -> Array[HingeSite] {
  let n_modes = if num_modes < modes.length() { num_modes } else { modes.length() }
  if n_modes == 0 {
    return []
  }
  let dim = modes[0].length()
  let n = dim / 3
  let counts = Array::make(n, 0)
  let min_mag = Array::make(n, 1.0e30)
  for k = 0; k < n_modes; k = k + 1 {
    let hinges = detect_hinges_anm(modes[k], k, threshold~)
    for h in hinges {
      let idx = h.residue_idx
      counts[idx] = counts[idx] + 1
      if h.magnitude < min_mag[idx] {
        min_mag[idx] = h.magnitude
      }
    }
  }
  let result : Array[HingeSite] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if counts[i] >= min_count {
      result.push(HingeSite::{ residue_idx: i, mode_idx: -1, magnitude: min_mag[i], flank_left: counts[i].to_double(), flank_right: n_modes.to_double() })
    }
  }
  result
}

// ============================================================================
// Hinge scoring and classification
// ============================================================================

///|
/// Compute a hinge score for each residue across multiple GNM modes.
///
/// Score = fraction of modes in which the residue is a hinge.  Range [0, 1].
pub fn hinge_scores_gnm(
  modes : Array[Array[Double]],
  num_modes : Int,
) -> Array[Double] {
  let n_modes = if num_modes < modes.length() { num_modes } else { modes.length() }
  if n_modes == 0 {
    return []
  }
  let n = modes[0].length()
  let scores = Array::make(n, 0.0)
  for k = 0; k < n_modes; k = k + 1 {
    let hinges = detect_hinges_gnm(modes[k], k)
    for h in hinges {
      scores[h.residue_idx] = scores[h.residue_idx] + 1.0
    }
  }
  let nf = n_modes.to_double()
  for i = 0; i < n; i = i + 1 {
    scores[i] = scores[i] / nf
  }
  scores
}

///|
/// Compute a hinge score for each residue across multiple ANM modes.
pub fn hinge_scores_anm(
  modes : Array[Array[Double]],
  num_modes : Int,
  threshold? : Double = 0.3,
) -> Array[Double] {
  let n_modes = if num_modes < modes.length() { num_modes } else { modes.length() }
  if n_modes == 0 {
    return []
  }
  let dim = modes[0].length()
  let n = dim / 3
  let scores = Array::make(n, 0.0)
  for k = 0; k < n_modes; k = k + 1 {
    let hinges = detect_hinges_anm(modes[k], k, threshold~)
    for h in hinges {
      scores[h.residue_idx] = scores[h.residue_idx] + 1.0
    }
  }
  let nf = n_modes.to_double()
  for i = 0; i < n; i = i + 1 {
    scores[i] = scores[i] / nf
  }
  scores
}

///|
/// Classify residues as hinge (true) or non-hinge (false) based on
/// a score threshold.
pub fn classify_hinges(scores : Array[Double], threshold : Double) -> Array[Bool] {
  let n = scores.length()
  let result = Array::make(n, false)
  for i = 0; i < n; i = i + 1 {
    result[i] = scores[i] >= threshold
  }
  result
}

// ============================================================================
// Hinge geometry — bending angle and rotation axis
// ============================================================================

///|
/// Compute the bending angle at a hinge site from an ANM mode.
///
/// Takes the average displacement vectors of residues flanking the hinge
/// and computes the angle between them.  Returns angle in degrees.
pub fn hinge_angle(
  mode : Array[Double],
  hinge_idx : Int,
  window? : Int = 2,
) -> Double {
  let dim = mode.length()
  let n = dim / 3
  let w = if window > 0 { window } else { 1 }
  let mut left_v = Vector3::zero()
  let mut left_count = 0
  let lo_l = if hinge_idx - w >= 0 { hinge_idx - w } else { 0 }
  for i = lo_l; i < hinge_idx; i = i + 1 {
    let i3 = 3 * i
    left_v = Vector3::new(
      left_v.x + mode[i3],
      left_v.y + mode[i3 + 1],
      left_v.z + mode[i3 + 2],
    )
    left_count = left_count + 1
  }
  let mut right_v = Vector3::zero()
  let mut right_count = 0
  let hi_r = if hinge_idx + 1 + w < n { hinge_idx + 1 + w } else { n }
  for i = hinge_idx + 1; i < hi_r; i = i + 1 {
    let i3 = 3 * i
    right_v = Vector3::new(
      right_v.x + mode[i3],
      right_v.y + mode[i3 + 1],
      right_v.z + mode[i3 + 2],
    )
    right_count = right_count + 1
  }
  if left_count == 0 || right_count == 0 {
    return 0.0
  }
  let lc = left_count.to_double()
  let rc = right_count.to_double()
  left_v = Vector3::new(left_v.x / lc, left_v.y / lc, left_v.z / lc)
  right_v = Vector3::new(right_v.x / rc, right_v.y / rc, right_v.z / rc)
  let left_norm = left_v.norm()
  let right_norm = right_v.norm()
  if left_norm < 1.0e-12 || right_norm < 1.0e-12 {
    return 0.0
  }
  let dot = left_v.x * right_v.x + left_v.y * right_v.y + left_v.z * right_v.z
  let cos_theta = dot / (left_norm * right_norm)
  let cos_clamped = if cos_theta > 1.0 { 1.0 } else { if cos_theta < -1.0 { -1.0 } else { cos_theta } }
  let theta = @math.acos(cos_clamped)
  theta * 180.0 / 3.14159265358979323846
}

///|
/// Compute the hinge rotation axis from an ANM mode.
///
/// The axis is the cross product of the average displacement vectors
/// of the two flanking regions, normalized.
pub fn hinge_axis(
  mode : Array[Double],
  hinge_idx : Int,
  window? : Int = 2,
) -> Vector3 {
  let dim = mode.length()
  let n = dim / 3
  let w = if window > 0 { window } else { 1 }
  let mut left_v = Vector3::zero()
  let mut left_count = 0
  let lo_l = if hinge_idx - w >= 0 { hinge_idx - w } else { 0 }
  for i = lo_l; i < hinge_idx; i = i + 1 {
    let i3 = 3 * i
    left_v = Vector3::new(
      left_v.x + mode[i3],
      left_v.y + mode[i3 + 1],
      left_v.z + mode[i3 + 2],
    )
    left_count = left_count + 1
  }
  let mut right_v = Vector3::zero()
  let mut right_count = 0
  let hi_r = if hinge_idx + 1 + w < n { hinge_idx + 1 + w } else { n }
  for i = hinge_idx + 1; i < hi_r; i = i + 1 {
    let i3 = 3 * i
    right_v = Vector3::new(
      right_v.x + mode[i3],
      right_v.y + mode[i3 + 1],
      right_v.z + mode[i3 + 2],
    )
    right_count = right_count + 1
  }
  if left_count == 0 || right_count == 0 {
    return Vector3::zero()
  }
  let lc = left_count.to_double()
  let rc = right_count.to_double()
  left_v = Vector3::new(left_v.x / lc, left_v.y / lc, left_v.z / lc)
  right_v = Vector3::new(right_v.x / rc, right_v.y / rc, right_v.z / rc)
  // Cross product: left × right
  let axis = Vector3::new(
    left_v.y * right_v.z - left_v.z * right_v.y,
    left_v.z * right_v.x - left_v.x * right_v.z,
    left_v.x * right_v.y - left_v.y * right_v.x,
  )
  let norm = axis.norm()
  if norm < 1.0e-12 {
    return Vector3::zero()
  }
  Vector3::new(axis.x / norm, axis.y / norm, axis.z / norm)
}

///|
/// Compute bending angles for all detected hinges in an ANM mode.
pub fn hinge_angles(
  mode : Array[Double],
  hinges : Array[HingeSite],
  window? : Int = 2,
) -> Array[Double] {
  let result : Array[Double] = Array(capacity=hinges.length())
  for h in hinges {
    result.push(hinge_angle(mode, h.residue_idx, window~))
  }
  result
}

///|
/// Summary of hinge analysis: list of hinge residues with scores,
/// angles, and axes for each mode.
pub struct HingeAnalysis {
  n_residues : Int
  num_modes : Int
  scores : Array[Double] // per-residue hinge score [0, 1]
  hinge_indices : Array[Int] // residues classified as hinges
  angles : Array[Double] // bending angle per hinge (degrees)
}

///|
/// Run full hinge analysis on GNM slow modes.
pub fn analyze_hinges_gnm(
  modes : Array[Array[Double]],
  num_modes : Int,
  score_threshold? : Double = 0.5,
) -> HingeAnalysis {
  let n_modes = if num_modes < modes.length() { num_modes } else { modes.length() }
  let n = if modes.length() > 0 { modes[0].length() } else { 0 }
  let scores = hinge_scores_gnm(modes, n_modes)
  let is_hinge = classify_hinges(scores, score_threshold)
  let hinge_indices : Array[Int] = Array(capacity=0)
  let angles : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if is_hinge[i] {
      hinge_indices.push(i)
      // For GNM, angle is not applicable; use 0.0
      angles.push(0.0)
    }
  }
  HingeAnalysis::{ n_residues: n, num_modes: n_modes, scores, hinge_indices, angles }
}

///|
/// Run full hinge analysis on ANM slow modes.
pub fn analyze_hinges_anm(
  modes : Array[Array[Double]],
  num_modes : Int,
  score_threshold? : Double = 0.5,
  threshold? : Double = 0.3,
) -> HingeAnalysis {
  let n_modes = if num_modes < modes.length() { num_modes } else { modes.length() }
  let dim = if modes.length() > 0 { modes[0].length() } else { 0 }
  let n = dim / 3
  let scores = hinge_scores_anm(modes, n_modes, threshold~)
  let is_hinge = classify_hinges(scores, score_threshold)
  let hinge_indices : Array[Int] = Array(capacity=0)
  let angles : Array[Double] = Array(capacity=0)
  // Use the slowest mode for angle computation
  for i = 0; i < n; i = i + 1 {
    if is_hinge[i] {
      hinge_indices.push(i)
      if n_modes > 0 {
        angles.push(hinge_angle(modes[0], i))
      } else {
        angles.push(0.0)
      }
    }
  }
  HingeAnalysis::{ n_residues: n, num_modes: n_modes, scores, hinge_indices, angles }
}

///|
/// Get the number of hinge residues found.
pub fn HingeAnalysis::num_hinges(self : HingeAnalysis) -> Int {
  self.hinge_indices.length()
}

///|
/// Get the hinge residue indices.
pub fn HingeAnalysis::get_hinge_indices(self : HingeAnalysis) -> Array[Int] {
  self.hinge_indices
}

///|
/// Get the per-residue hinge scores.
pub fn HingeAnalysis::get_scores(self : HingeAnalysis) -> Array[Double] {
  self.scores
}

///|
/// Get the bending angles (degrees) for each hinge.
pub fn HingeAnalysis::get_angles(self : HingeAnalysis) -> Array[Double] {
  self.angles
}
