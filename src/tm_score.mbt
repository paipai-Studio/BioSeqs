///|
/// TM-score (Template Modeling score) — length-independent protein structure
/// similarity metric (Zhang & Skolnick 2004).
///
/// Complements RMSD superposition with the standard score: given two
/// coordinate sets whose residue correspondence is known, the optimal
/// (quaternion) superposition is found and the TM-score is
/// `(1/L) * sum_i 1/(1 + (d_i / d0)^2)` with
/// `d0 = 1.24 * (L-15)^{1/3} - 1.8` (and `d0 = 0.5` when `L < 22`).

///|
/// The TM-score length scale `d0`.
pub fn tm_score_d0(l : Int) -> Double {
  if l < 22 {
    0.5
  } else {
    1.24 * @math.pow(l.to_double() - 15.0, 1.0 / 3.0) - 1.8
  }
}

///|
/// TM-score for two already-aligned, already-superimposed, equal-length
/// coordinate sets. `l_target` is the normalizing (target) length, typically
/// the number of residues of the larger structure.
pub fn tm_score_aligned(
  coords1 : Array[QsPoint],
  coords2 : Array[QsPoint],
  l_target : Int,
) -> Double {
  let n = if coords1.length() < coords2.length() {
    coords1.length()
  } else {
    coords2.length()
  }
  if n == 0 || l_target == 0 {
    return 0.0
  }
  let d0 = tm_score_d0(l_target)
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    let dx = coords1[i].x - coords2[i].x
    let dy = coords1[i].y - coords2[i].y
    let dz = coords1[i].z - coords2[i].z
    let d = (dx * dx + dy * dy + dz * dz).sqrt()
    let q = d / d0
    sum = sum + 1.0 / (1.0 + q * q)
    i = i + 1
  }
  sum / l_target.to_double()
}

///|
/// TM-score after optimal superposition. `coords1` is rotated/translated onto
/// `coords2` (quaternion method), then the score is computed with `l_target`
/// equal to the larger structure's length. Returns 0.0 when either set is
/// empty.
pub fn tm_score(
  coords1 : Array[QsPoint],
  coords2 : Array[QsPoint],
) -> Double {
  let n1 = coords1.length()
  let n2 = coords2.length()
  if n1 == 0 || n2 == 0 {
    return 0.0
  }
  let l_target = if n1 > n2 { n1 } else { n2 }
  // qs_superimpose(a, b) 返回「b → a」的变换，所以变换 coords2 到 coords1 坐标系
  let r = qs_superimpose(coords1, coords2)
  let transformed = qs_apply_transform(coords2, r.rotation, r.translation)
  tm_score_aligned(coords1, transformed, l_target)
}
