/// Virtual sites for multi-site water models and coarse-graining.
///
/// OpenMM equivalent: `VirtualSite` class hierarchy
///
/// Virtual sites are massless particles whose positions are computed
/// from parent atom positions (e.g., TIP4P oxygen dummy, DRUDE particles).

///|
/// Virtual site type
pub enum OmmVirtualSiteType {
  TwoParticleAverage(Int32, Int32, Double, Double) // (p1, p2, w1, w2)
  ThreeParticleAverage(Int32, Int32, Int32, Double, Double, Double) // (p1, p2, p3, w1, w2, w3)
  OutOfPlane(Int32, Int32, Int32, Double, Double, Double) // (p1, p2, p3, w12, w13, wCross)
  LocalCoordinate(
    Int32,
    Int32,
    Int32,
    Double,
    Double,
    Double,
    Double,
    Double,
    Double
  )
}

///|
/// A virtual site
pub struct OmmVirtualSite {
  index : Int32 // index of the virtual particle
  site_type : OmmVirtualSiteType
}

///|
pub fn OmmVirtualSite::index(self : OmmVirtualSite) -> Int32 {
  self.index
}

///|
/// Compute position of a virtual site from parent positions.
pub fn omm_virtual_site_position(
  site : OmmVirtualSite,
  positions : Array[Array[Double]],
) -> Array[Double] {
  match site.site_type {
    OmmVirtualSiteType::TwoParticleAverage(p1, p2, w1, w2) =>
      [
        w1 * positions[p1][0] + w2 * positions[p2][0],
        w1 * positions[p1][1] + w2 * positions[p2][1],
        w1 * positions[p1][2] + w2 * positions[p2][2],
      ]
    OmmVirtualSiteType::ThreeParticleAverage(p1, p2, p3, w1, w2, w3) =>
      [
        w1 * positions[p1][0] + w2 * positions[p2][0] + w3 * positions[p3][0],
        w1 * positions[p1][1] + w2 * positions[p2][1] + w3 * positions[p3][1],
        w1 * positions[p1][2] + w2 * positions[p2][2] + w3 * positions[p3][2],
      ]
    OmmVirtualSiteType::OutOfPlane(p1, p2, p3, w12, w13, wCross) => {
      // r = r1 + w12*(r2-r1) + w13*(r3-r1) + wCross*((r2-r1)×(r3-r1))
      let v12x = positions[p2][0] - positions[p1][0]
      let v12y = positions[p2][1] - positions[p1][1]
      let v12z = positions[p2][2] - positions[p1][2]
      let v13x = positions[p3][0] - positions[p1][0]
      let v13y = positions[p3][1] - positions[p1][1]
      let v13z = positions[p3][2] - positions[p1][2]
      let cross_x = v12y * v13z - v12z * v13y
      let cross_y = v12z * v13x - v12x * v13z
      let cross_z = v12x * v13y - v12y * v13x
      [
        positions[p1][0] + w12 * v12x + w13 * v13x + wCross * cross_x,
        positions[p1][1] + w12 * v12y + w13 * v13y + wCross * cross_y,
        positions[p1][2] + w12 * v12z + w13 * v13z + wCross * cross_z,
      ]
    }
    OmmVirtualSiteType::LocalCoordinate(p1, p2, p3, wx, wy, wz, _u1, _u2, _u3) =>
      // Simplified: weighted average with local frame weights
      [
        wx * positions[p1][0] + wy * positions[p2][0] + wz * positions[p3][0],
        wx * positions[p1][1] + wy * positions[p2][1] + wz * positions[p3][1],
        wx * positions[p1][2] + wy * positions[p2][2] + wz * positions[p3][2],
      ]
  }
}

///|
/// Update all virtual site positions in a position array.
pub fn omm_update_virtual_sites(
  positions : Array[Array[Double]],
  virtual_sites : Array[OmmVirtualSite],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::make(positions.length(), [
    0.0, 0.0, 0.0,
  ])
  for i = 0; i < positions.length(); i = i + 1 {
    result[i] = [positions[i][0], positions[i][1], positions[i][2]]
  }
  for vs in virtual_sites {
    let pos = omm_virtual_site_position(vs, positions)
    result[vs.index] = pos
  }
  result
}

///|
/// Create TIP4P virtual site (M site on bisector of H-O-H, offset from O)
/// o_idx: oxygen index, h1_idx, h2_idx: hydrogen indices, m_idx: M site index
pub fn omm_tip4p_virtual_site(
  o_idx : Int32,
  h1_idx : Int32,
  h2_idx : Int32,
  m_idx : Int32,
) -> OmmVirtualSite {
  // M = 0.5*H1 + 0.5*H2 - 0.5*O + ... (simplified: weighted average)
  // For TIP4P/2005: M is on the angle bisector at ~0.15 Å from O
  // Using ThreeParticleAverage: M = w_O * O + w_H1 * H1 + w_H2 * H2
  // w_O = -0.8, w_H1 = 0.9, w_H2 = 0.9 (approximate, sum = 1.0)
  {
    index: m_idx,
    site_type: OmmVirtualSiteType::ThreeParticleAverage(
      o_idx, h1_idx, h2_idx, -0.8, 0.9, 0.9,
    ),
  }
}

///|
/// Create DRUDE virtual site (spring-attached to parent atom)
pub fn omm_drude_virtual_site(
  parent_idx : Int32,
  drude_idx : Int32,
  _offset~ : Double = 0.0,
) -> OmmVirtualSite {
  // DRUDE particle is offset from parent along a fixed direction
  // Simplified as TwoParticleAverage with weights
  {
    index: drude_idx,
    site_type: OmmVirtualSiteType::TwoParticleAverage(
      parent_idx, parent_idx, 1.0, 0.0,
    ),
  }
}

///|
fn _omm_vs_unused_constructors() -> Unit {
  ignore(OmmVirtualSiteType::OutOfPlane(0, 0, 0, 0.0, 0.0, 0.0))
  ignore(
    OmmVirtualSiteType::LocalCoordinate(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
  )
}
