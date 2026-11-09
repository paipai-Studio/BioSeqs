///|
/// HELANAL helix axis ensemble analysis, mirroring MDAnalysis
/// `analysis.helix_analysis` (HELANAL). For every frame the global helix
/// axis is obtained as the principal direction of the covariance matrix of
/// all C-alpha points, while a sliding window of 4 consecutive residues
/// yields local helix axes. Reported per frame are the bending of each
/// window axis relative to the global axis and the twist (relative
/// rotation) between neighbouring window axes, both in degrees.

///|
/// Dot product of two 3-vectors.
fn mda_vec3_dot(a : Array[Double], b : Array[Double]) -> Double {
  a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

///|
/// Euclidean norm of a 3-vector.
fn mda_vec3_norm(a : Array[Double]) -> Double {
  mda_vec3_dot(a, a).sqrt()
}

///|
/// Angle between two 3-vectors in degrees (0.0 if either vector is zero).
fn mda_vec3_angle_deg(a : Array[Double], b : Array[Double]) -> Double {
  let na = mda_vec3_norm(a)
  let nb = mda_vec3_norm(b)
  if na == 0.0 || nb == 0.0 {
    return 0.0
  }
  let mut c = mda_vec3_dot(a, b) / (na * nb)
  if c > 1.0 {
    c = 1.0
  } else if c < -1.0 {
    c = -1.0
  }
  @math.acos(c) * 57.29577951308232
}

///|
/// Deterministic sign convention: flip the vector so that its
/// largest-magnitude component is positive.
fn mda_vec3_sign_fix(v : Array[Double]) -> Array[Double] {
  let mut mi = 0
  for i = 1; i < 3; i = i + 1 {
    if v[i] * v[i] > v[mi] * v[mi] {
      mi = i
    }
  }
  if v[mi] < 0.0 {
    [-v[0], -v[1], -v[2]]
  } else {
    [v[0], v[1], v[2]]
  }
}

///|
/// Centred 3x3 covariance matrix of a set of points
/// `S = sum_i (r_i - mean)(r_i - mean)^T`.
fn mda_helanal_covariance(pts : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = pts.length()
  if n == 0 {
    return [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  }
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for p in pts {
    cx = cx + p[0]
    cy = cy + p[1]
    cz = cz + p[2]
  }
  cx = cx / n.to_double()
  cy = cy / n.to_double()
  cz = cz / n.to_double()
  let mut sxx = 0.0
  let mut syy = 0.0
  let mut szz = 0.0
  let mut sxy = 0.0
  let mut sxz = 0.0
  let mut syz = 0.0
  for p in pts {
    let dx = p[0] - cx
    let dy = p[1] - cy
    let dz = p[2] - cz
    sxx = sxx + dx * dx
    syy = syy + dy * dy
    szz = szz + dz * dz
    sxy = sxy + dx * dy
    sxz = sxz + dx * dz
    syz = syz + dy * dz
  }
  [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
}

///|
/// Dominant eigenvector of a symmetric 3x3 matrix by power iteration.
/// Starts from the coordinate axes in turn to escape zero subspaces; for a
/// fully zero matrix (all points coincident) the fallback `[0, 0, 1]` is
/// returned. The result is normalised and sign-fixed deterministically.
fn mda_helanal_power_iter(m : Array[Array[Double]]) -> Array[Double] {
  let inits = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  for k = 0; k < 3; k = k + 1 {
    let mut vx = inits[k][0]
    let mut vy = inits[k][1]
    let mut vz = inits[k][2]
    let mut converged = false
    for it = 0; it < 200; it = it + 1 {
      let wx = m[0][0] * vx + m[0][1] * vy + m[0][2] * vz
      let wy = m[1][0] * vx + m[1][1] * vy + m[1][2] * vz
      let wz = m[2][0] * vx + m[2][1] * vy + m[2][2] * vz
      let n = (wx * wx + wy * wy + wz * wz).sqrt()
      if n < 1.0e-20 {
        break
      }
      vx = wx / n
      vy = wy / n
      vz = wz / n
      converged = true
    }
    if converged {
      return mda_vec3_sign_fix([vx, vy, vz])
    }
  }
  [0.0, 0.0, 1.0]
}

///|
/// Fit a local helix axis through 4 consecutive backbone points: the unit
/// principal direction (dominant eigenvector) of their centred covariance
/// matrix. The sign is fixed deterministically (largest-magnitude
/// component positive).
pub fn mda_helix_axis_window(coords4 : Array[Array[Double]]) -> Array[Double] {
  let cov = mda_helanal_covariance(coords4)
  mda_helanal_power_iter(cov)
}

///|
/// Per-frame HELANAL result.
///
/// global_axis: unit global helix axis of the whole fragment
/// window_axes: unit local axis per 4-residue window (n_res - 3 entries)
/// bend_angles: angle (degrees) of each window axis relative to the
///   global axis, one per window
/// twist_angles: angle (degrees) between neighbouring window axes
///   (window_axes.length() - 1 entries)
pub(all) struct HelanalFrame {
  global_axis : Array[Double]
  window_axes : Array[Array[Double]]
  bend_angles : Array[Double]
  twist_angles : Array[Double]
}

///|
/// HELANAL helix axis ensemble analysis over a trajectory
/// (MDAnalysis `analysis.helix_analysis`).
///
/// traj: frames, each a helix fragment of C-alpha coordinates
/// `[n_res][3]` with `n_res >= 4`; frames with fewer than 4 residues yield
/// empty window/angle arrays and a zero global axis.
///
/// For every frame the global axis is the principal direction of all
/// C-alpha points; each window axis (4 consecutive residues) is sign-
/// aligned with the global axis so that bend and twist angles fall in
/// `[0, 90]` degrees.
pub fn mda_helanal(traj : Array[Array[Array[Double]]]) -> Array[HelanalFrame] {
  let out : Array[HelanalFrame] = Array(capacity=0)
  for frame in traj {
    let n = frame.length()
    if n < 4 {
      out.push({
        global_axis: [0.0, 0.0, 0.0],
        window_axes: [],
        bend_angles: [],
        twist_angles: [],
      })
      continue
    }
    let global = mda_helanal_power_iter(mda_helanal_covariance(frame))
    let axes : Array[Array[Double]] = Array(capacity=0)
    for i = 0; i + 4 <= n; i = i + 1 {
      let w = mda_helix_axis_window([
        frame[i],
        frame[i + 1],
        frame[i + 2],
        frame[i + 3],
      ])
      if mda_vec3_dot(w, global) < 0.0 {
        axes.push([-w[0], -w[1], -w[2]])
      } else {
        axes.push(w)
      }
    }
    let bends : Array[Double] = Array(capacity=axes.length())
    for w in axes {
      bends.push(mda_vec3_angle_deg(w, global))
    }
    let twists : Array[Double] = Array(capacity=axes.length())
    for i = 0; i + 1 < axes.length(); i = i + 1 {
      twists.push(mda_vec3_angle_deg(axes[i], axes[i + 1]))
    }
    out.push({
      global_axis: global,
      window_axes: axes,
      bend_angles: bends,
      twist_angles: twists,
    })
  }
  out
}
