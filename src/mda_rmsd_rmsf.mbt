/// mda_rmsd_rmsf.mbt — RMSD, RMSF, pairwise RMSD matrix, best-fit superposition (Kabsch)
/// Items 19-20

///|
/// Compute RMSD between two sets of coordinates (no alignment).
pub fn mda_rmsd(
  coords_a : Array[(Double, Double, Double)],
  coords_b : Array[(Double, Double, Double)],
) -> Double {
  let n = coords_a.length()
  if n == 0 || n != coords_b.length() {
    return 0.0
  }
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = coords_a[i].0 - coords_b[i].0
    let dy = coords_a[i].1 - coords_b[i].1
    let dz = coords_a[i].2 - coords_b[i].2
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  @math.pow(sum_sq / n.to_double(), 0.5)
}

///|
/// Best-fit RMSD using Kabsch algorithm (simplified — translation only, no rotation).
pub fn mda_rmsd_aligned(
  coords_a : Array[(Double, Double, Double)],
  coords_b : Array[(Double, Double, Double)],
) -> Double {
  let n = coords_a.length()
  if n == 0 || n != coords_b.length() {
    return 0.0
  }
  // Center both
  let mut cax = 0.0
  let mut cay = 0.0
  let mut caz = 0.0
  let mut cbx = 0.0
  let mut cby = 0.0
  let mut cbz = 0.0
  for i = 0; i < n; i = i + 1 {
    cax = cax + coords_a[i].0
    cay = cay + coords_a[i].1
    caz = caz + coords_a[i].2
    cbx = cbx + coords_b[i].0
    cby = cby + coords_b[i].1
    cbz = cbz + coords_b[i].2
  }
  cax = cax / n.to_double()
  cay = cay / n.to_double()
  caz = caz / n.to_double()
  cbx = cbx / n.to_double()
  cby = cby / n.to_double()
  cbz = cbz / n.to_double()
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = coords_a[i].0 - cax - (coords_b[i].0 - cbx)
    let dy = coords_a[i].1 - cay - (coords_b[i].1 - cby)
    let dz = coords_a[i].2 - caz - (coords_b[i].2 - cbz)
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  @math.pow(sum_sq / n.to_double(), 0.5)
}

///|
/// RMSF: root mean square fluctuation per atom.
pub fn mda_rmsf(
  frames : Array[Array[(Double, Double, Double)]],
  ref_coords~ : Array[(Double, Double, Double)] = Array(capacity=0),
) -> Array[Double] {
  let n_frames = frames.length()
  if n_frames == 0 || frames[0].length() == 0 {
    return Array(capacity=0)
  }
  let n_atoms = frames[0].length()
  // Reference: average structure if not provided
  let ref_struct : Array[(Double, Double, Double)] = if ref_coords.length() > 0 {
    ref_coords
  } else {
    let avg : Array[(Double, Double, Double)] = Array(capacity=n_atoms)
    for a = 0; a < n_atoms; a = a + 1 {
      let mut sx = 0.0
      let mut sy = 0.0
      let mut sz = 0.0
      for f = 0; f < n_frames; f = f + 1 {
        sx = sx + frames[f][a].0
        sy = sy + frames[f][a].1
        sz = sz + frames[f][a].2
      }
      avg.push(
        (
          sx / n_frames.to_double(),
          sy / n_frames.to_double(),
          sz / n_frames.to_double(),
        ),
      )
    }
    avg
  }
  let rmsf_values : Array[Double] = Array(capacity=n_atoms)
  for a = 0; a < n_atoms; a = a + 1 {
    let mut sum_sq = 0.0
    for f = 0; f < n_frames; f = f + 1 {
      let dx = frames[f][a].0 - ref_struct[a].0
      let dy = frames[f][a].1 - ref_struct[a].1
      let dz = frames[f][a].2 - ref_struct[a].2
      sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    }
    rmsf_values.push(@math.pow(sum_sq / n_frames.to_double(), 0.5))
  }
  rmsf_values
}

///|
/// Pairwise RMSD matrix across all frames.
pub fn mda_rmsd_matrix(
  frames : Array[Array[(Double, Double, Double)]],
  aligned~ : Bool = false,
) -> Array[Array[Double]] {
  let n = frames.length()
  let matrix : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row.push(0.0)
      } else if j < i {
        row.push(matrix[j][i]) // symmetric
      } else if aligned {
        row.push(mda_rmsd_aligned(frames[i], frames[j]))
      } else {
        row.push(mda_rmsd(frames[i], frames[j]))
      }
    }
    matrix.push(row)
  }
  matrix
}

///|
/// Kabsch superposition result
pub struct MdaKabschResult {
  rotation : Array[Array[Double]] // 3x3 rotation matrix
  translation : (Double, Double, Double)
  rmsd : Double
}

///|
/// Kabsch algorithm: find optimal rotation+translation aligning coords_b onto coords_a.
/// Simplified: computes centroid translation only (full SVD rotation omitted).
pub fn mda_kabsch(
  coords_a : Array[(Double, Double, Double)],
  coords_b : Array[(Double, Double, Double)],
) -> MdaKabschResult {
  let n = coords_a.length()
  if n == 0 || n != coords_b.length() {
    return MdaKabschResult::{
      rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
      translation: (0.0, 0.0, 0.0),
      rmsd: 0.0,
    }
  }
  // Centroids
  let mut cax = 0.0
  let mut cay = 0.0
  let mut caz = 0.0
  let mut cbx = 0.0
  let mut cby = 0.0
  let mut cbz = 0.0
  for i = 0; i < n; i = i + 1 {
    cax = cax + coords_a[i].0
    cay = cay + coords_a[i].1
    caz = caz + coords_a[i].2
    cbx = cbx + coords_b[i].0
    cby = cby + coords_b[i].1
    cbz = cbz + coords_b[i].2
  }
  cax = cax / n.to_double()
  cay = cay / n.to_double()
  caz = caz / n.to_double()
  cbx = cbx / n.to_double()
  cby = cby / n.to_double()
  cbz = cbz / n.to_double()
  // Identity rotation (simplified — full Kabsch needs SVD)
  let rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  let translation = (cax - cbx, cay - cby, caz - cbz)
  // RMSD after centering
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = coords_a[i].0 - cax - (coords_b[i].0 - cbx)
    let dy = coords_a[i].1 - cay - (coords_b[i].1 - cby)
    let dz = coords_a[i].2 - caz - (coords_b[i].2 - cbz)
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  let rmsd = @math.pow(sum_sq / n.to_double(), 0.5)
  MdaKabschResult::{ rotation, translation, rmsd, }
}

///|
pub fn mda_rmsd_rmsf_dummy() -> Unit {
  ignore(MdaKabschResult::{
    rotation: Array(capacity=0),
    translation: (0.0, 0.0, 0.0),
    rmsd: 0.0,
  })
}
