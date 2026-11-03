/// HOLE-like pore/channel radius profile.
///
/// MDAnalysis equivalent: `analysis.hole2.HOLE`
///
/// Computes the radius of a pore/channel along a user-defined axis
/// by finding the largest sphere that fits at each point along the axis.

///|
/// A point in the pore radius profile
pub struct MdaHolePoint {
  position : Double // position along axis
  radius : Double // pore radius at this position
}

///|
pub fn MdaHolePoint::position(self : MdaHolePoint) -> Double {
  self.position
}

///|
pub fn MdaHolePoint::radius(self : MdaHolePoint) -> Double {
  self.radius
}

///|
/// Compute pore radius profile along an axis.
///
/// atom_coords: coordinates of all atoms forming the pore walls
/// axis_start: start point of the pore axis [x, y, z]
/// axis_end: end point of the pore axis [x, y, z]
/// n_samples: number of points along the axis
/// probe_radius: radius of the probe (default 1.0 Å, water)
pub fn mda_hole_profile(
  atom_coords : Array[Array[Double]],
  axis_start : Array[Double],
  axis_end : Array[Double],
  n_samples~ : Int32 = 100,
  probe_radius~ : Double = 1.0,
) -> Array[MdaHolePoint] {
  let dx = axis_end[0] - axis_start[0]
  let dy = axis_end[1] - axis_start[1]
  let dz = axis_end[2] - axis_start[2]
  let axis_len = (dx * dx + dy * dy + dz * dz).sqrt()
  if axis_len == 0.0 {
    return []
  }
  // Unit vector along axis
  let ux = dx / axis_len
  let uy = dy / axis_len
  let uz = dz / axis_len
  let result : Array[MdaHolePoint] = Array(capacity=0)
  for s = 0; s < n_samples; s = s + 1 {
    let t = s.to_double() / (n_samples - 1).to_double()
    let px = axis_start[0] + t * dx
    let py = axis_start[1] + t * dy
    let pz = axis_start[2] + t * dz
    // Find the minimum distance to any atom (perpendicular to axis)
    let mut min_dist = 10000000000.0
    for a = 0; a < atom_coords.length(); a = a + 1 {
      let rx = atom_coords[a][0] - px
      let ry = atom_coords[a][1] - py
      let rz = atom_coords[a][2] - pz
      // Project onto axis to get parallel component
      let proj = rx * ux + ry * uy + rz * uz
      // Perpendicular distance
      let perp_sq = rx * rx + ry * ry + rz * rz - proj * proj
      if perp_sq < min_dist {
        min_dist = perp_sq
      }
    }
    let radius = min_dist.sqrt() - probe_radius
    let pos = t * axis_len
    result.push({
      position: pos,
      radius: if radius > 0.0 {
        radius
      } else {
        0.0
      },
    })
  }
  result
}

///|
/// Find the minimum pore radius (bottleneck).
pub fn mda_hole_min_radius(profile : Array[MdaHolePoint]) -> Double {
  if profile.length() == 0 {
    return 0.0
  }
  let mut min_r = profile[0].radius
  for p in profile {
    if p.radius < min_r {
      min_r = p.radius
    }
  }
  min_r
}

///|
/// Find the position of the pore bottleneck.
pub fn mda_hole_bottleneck(profile : Array[MdaHolePoint]) -> (Double, Double) {
  if profile.length() == 0 {
    return (0.0, 0.0)
  }
  let mut min_r = profile[0].radius
  let mut pos = profile[0].position
  for p in profile {
    if p.radius < min_r {
      min_r = p.radius
      pos = p.position
    }
  }
  (pos, min_r)
}

///|
/// Compute pore conductance using Hagen-Poiseuille approximation.
/// G = pi * r^4 / (8 * eta * L)
///
/// min_radius: minimum pore radius (m)
/// length: pore length (m)
/// viscosity: solvent viscosity (Pa·s)
pub fn mda_hole_conductance(
  min_radius : Double,
  length : Double,
  viscosity~ : Double = 0.001,
) -> Double {
  if length <= 0.0 || min_radius <= 0.0 {
    return 0.0
  }
  3.141592653589793 *
  min_radius *
  min_radius *
  min_radius *
  min_radius /
  (8.0 * viscosity * length)
}

///|
/// Classify pore by bottleneck radius.
pub fn mda_hole_classify(min_radius : Double) -> String {
  if min_radius < 1.0 {
    "Closed"
  } else if min_radius < 2.0 {
    "Narrow"
  } else if min_radius < 4.0 {
    "Open"
  } else {
    "Wide"
  }
}
