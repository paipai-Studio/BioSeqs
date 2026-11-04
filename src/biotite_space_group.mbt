/// biotite_space_group.mbt — Crystallographic space group symmetry operations
/// Inspired by biotite.structure.SpaceGroup

///|
/// A symmetry operation: rotation matrix (3x3) + translation vector (3).
pub struct BiotiteSymmetryOp {
  rotation : Array[Double] // 9 elements, row-major
  translation : Array[Double] // 3 elements
}

///|
/// Apply symmetry operation to a fractional coordinate.
pub fn BiotiteSymmetryOp::apply(
  self : BiotiteSymmetryOp,
  coord : Array[Double],
) -> Array[Double] {
  let r = self.rotation
  let t = self.translation
  [
    r[0] * coord[0] + r[1] * coord[1] + r[2] * coord[2] + t[0],
    r[3] * coord[0] + r[4] * coord[1] + r[5] * coord[2] + t[1],
    r[6] * coord[0] + r[7] * coord[1] + r[8] * coord[2] + t[2],
  ]
}

///|
/// Compose two symmetry operations: (self ∘ other).
pub fn BiotiteSymmetryOp::compose(
  self : BiotiteSymmetryOp,
  other : BiotiteSymmetryOp,
) -> BiotiteSymmetryOp {
  let a = self.rotation
  let b = other.rotation
  let at = self.translation
  let bt = other.translation
  // R = Ra * Rb, t = Ra * tb + ta
  let rotation : Array[Double] = Array::make(9, 0.0)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < 3; k = k + 1 {
        sum = sum + a[i * 3 + k] * b[k * 3 + j]
      }
      rotation[i * 3 + j] = sum
    }
  }
  let translation : Array[Double] = Array::make(3, 0.0)
  for i = 0; i < 3; i = i + 1 {
    let mut sum = at[i]
    for k = 0; k < 3; k = k + 1 {
      sum = sum + a[i * 3 + k] * bt[k]
    }
    translation[i] = sum
  }
  { rotation, translation, }
}

///|
/// Invert a symmetry operation.
pub fn BiotiteSymmetryOp::inverse(
  self : BiotiteSymmetryOp,
) -> BiotiteSymmetryOp {
  let r = self.rotation
  let t = self.translation
  // For orthogonal matrices, inverse = transpose
  let inv_r : Array[Double] = [
    r[0],
    r[3],
    r[6],
    r[1],
    r[4],
    r[7],
    r[2],
    r[5],
    r[8],
  ]
  // t_inv = -R^T * t
  let inv_t : Array[Double] = [
    -(inv_r[0] * t[0] + inv_r[1] * t[1] + inv_r[2] * t[2]),
    -(inv_r[3] * t[0] + inv_r[4] * t[1] + inv_r[5] * t[2]),
    -(inv_r[6] * t[0] + inv_r[7] * t[1] + inv_r[8] * t[2]),
  ]
  { rotation: inv_r, translation: inv_t, }
}

///|
/// Identity operation.
fn biotite_symmetry_identity() -> BiotiteSymmetryOp {
  {
    rotation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    translation: [0.0, 0.0, 0.0],
  }
}

///|
/// Get all symmetry operations for a space group.
/// Supports a subset of the 230 space groups.
pub fn biotite_space_group_transforms(
  space_group : String,
) -> Array[BiotiteSymmetryOp] {
  let id = biotite_symmetry_identity()
  match space_group {
    "P1" => [id]
    "P-1" =>
      [
        id,
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.0, 0.0, 0.0],
        },
      ]
    "P2" =>
      [
        id,
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.0, 0.0, 0.0],
        },
      ]
    "P21" =>
      [
        id,
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.0, 0.5, 0.0],
        },
      ]
    "P222" =>
      [
        id,
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0],
          translation: [0.0, 0.0, 0.0],
        },
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.0, 0.0, 0.0],
        },
        {
          rotation: [1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.0, 0.0, 0.0],
        },
      ]
    "P212121" =>
      [
        id,
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0],
          translation: [0.5, 0.0, 0.0],
        },
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.5, 0.5, 0.0],
        },
        {
          rotation: [1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.0, 0.5, 0.0],
        },
      ]
    "C2" =>
      [
        id,
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.0, 0.0, 0.0],
        },
        {
          rotation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
          translation: [0.5, 0.5, 0.0],
        },
        {
          rotation: [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
          translation: [0.5, 0.5, 0.0],
        },
      ]
    _ => [id]
  }
}

///|
/// Generate all equivalent positions for a given fractional coordinate.
pub fn biotite_space_group_equivalent_positions(
  space_group : String,
  coord : Array[Double],
) -> Array[Array[Double]] {
  let ops = biotite_space_group_transforms(space_group)
  let result : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < ops.length(); i = i + 1 {
    let c = ops[i].apply(coord)
    // Wrap to [0, 1)
    result.push([
      c[0] - @math.floor(c[0]),
      c[1] - @math.floor(c[1]),
      c[2] - @math.floor(c[2]),
    ])
  }
  result
}
