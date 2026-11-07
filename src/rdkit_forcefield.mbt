/// rdkit_forcefield.mbt — UFF-style force field + conformer sampling
/// (RDKit gap fill, item A2)
///
/// Simplified UFF-style potential on the bond-order graph: bond stretch,
/// angle bend (hybridization-based ideal angles), torsions and a
/// Lennard-Jones non-bonded term (1-2 / 1-3 pairs excluded). Geometry is
/// minimized with steepest descent on a numerically evaluated gradient,
/// which keeps the implementation small and dependency-free at the cost
/// of speed (fine for the small molecules this library targets).

///|
/// Build a conformer directly from atom symbols and coordinates (test /
/// programmatic construction helper).
pub fn RdkitConformer::from_coords(
  atoms : Array[String],
  coords : Array[(Double, Double, Double)],
) -> RdkitConformer {
  RdkitConformer::{ atoms, coords, energy: 0.0, }
}

///|
fn ff_bond_r0(a : String, b : String, order : Int32) -> Double {
  let key = a + "-" + b
  let base = match key {
    "C-C" => 1.52
    "C-N" => 1.47
    "C-O" => 1.43
    "C-H" => 1.09
    "C-S" => 1.82
    "N-H" => 1.01
    "O-H" => 0.96
    "N-N" => 1.45
    "N-O" => 1.40
    "O-O" => 1.48
    "S-O" => 1.43
    "S-H" => 1.34
    "C-F" => 1.35
    "C-Cl" => 1.77
    "C-Br" => 1.94
    "C-I" => 2.14
    _ => 1.50
  }
  match order {
    2 => base * 0.875
    3 => base * 0.78
    _ => base
  }
}

///|
fn ff_bond_kb(order : Int32) -> Double {
  match order {
    2 => 900.0
    3 => 1300.0
    _ => 500.0
  }
}

///|
/// Ideal angle (radians) for the central atom given its connectivity.
fn ff_angle_theta0(
  mol : Molecule,
  center : Int32,
  neighbors : Array[Array[Int32]],
) -> Double {
  let deg = neighbors[center].length()
  let mut double_bonds = 0
  for b in mol.bonds {
    if b.from == center || b.to == center {
      if b.order >= 2 {
        double_bonds = double_bonds + 1
      }
    }
  }
  if deg == 2 && double_bonds >= 2 {
    return 3.141592653589793 // sp, linear
  }
  if double_bonds >= 1 {
    return 2.0943951023931953 // sp2, 120 deg
  }
  1.9106332362490186 // sp3, 109.47 deg
}

///|
fn ff_vdw_radius(symbol : String) -> Double {
  match symbol {
    "H" => 1.20
    "C" => 1.90
    "N" => 1.80
    "O" => 1.75
    "F" => 1.70
    "S" => 2.00
    "Cl" => 1.95
    "Br" => 2.10
    "I" => 2.25
    _ => 1.85
  }
}

///|
fn ff_vdw_epsilon(symbol : String) -> Double {
  match symbol {
    "H" => 0.03
    "C" => 0.07
    "N" => 0.07
    "O" => 0.09
    "S" => 0.25
    "Cl" => 0.24
    "Br" => 0.28
    "I" => 0.32
    _ => 0.06
  }
}

///|
fn ff_adj(mol : Molecule) -> Array[Array[Int32]] {
  let n = mol.num_atoms()
  let neighbors : Array[Array[Int32]] = Array(capacity=0)
  for _k = 0; _k < n; _k = _k + 1 {
    neighbors.push([])
  }
  for b in mol.bonds {
    if b.from < n && b.to < n {
      neighbors[b.from].push(b.to)
      neighbors[b.to].push(b.from)
    }
  }
  neighbors
}

///|
/// Euclidean distance between two atoms in a coordinate array (public
/// helper for geometry inspection).
pub fn ff_bond_length(
  coords : Array[(Double, Double, Double)],
  i : Int32,
  j : Int32,
) -> Double {
  ff_dist(coords, i, j)
}

///|
fn ff_dist(
  coords : Array[(Double, Double, Double)],
  i : Int32,
  j : Int32,
) -> Double {
  let (x1, y1, z1) = coords[i]
  let (x2, y2, z2) = coords[j]
  let dx = x1 - x2
  let dy = y1 - y2
  let dz = z1 - z2
  rdkit_sqrt(dx * dx + dy * dy + dz * dz)
}

///|
/// Angle at atom b (in radians) formed by a-b-c.
fn ff_angle(
  coords : Array[(Double, Double, Double)],
  a : Int32,
  b : Int32,
  c : Int32,
) -> Double {
  let (xa, ya, za) = coords[a]
  let (xb, yb, zb) = coords[b]
  let (xc, yc, zc) = coords[c]
  let v1x = xa - xb
  let v1y = ya - yb
  let v1z = za - zb
  let v2x = xc - xb
  let v2y = yc - yb
  let v2z = zc - zb
  let d1 = rdkit_sqrt(v1x * v1x + v1y * v1y + v1z * v1z)
  let d2 = rdkit_sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
  if d1 < 0.0000001 || d2 < 0.0000001 {
    return 0.0
  }
  let mut cosv = (v1x * v2x + v1y * v2y + v1z * v2z) / (d1 * d2)
  if cosv > 1.0 {
    cosv = 1.0
  }
  if cosv < -1.0 {
    cosv = -1.0
  }
  @math.acos(cosv)
}

///|
/// Torsion angle a-b-c-d in radians (-pi..pi).
fn ff_torsion(
  coords : Array[(Double, Double, Double)],
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
) -> Double {
  let (xa, ya, za) = coords[a]
  let (xb, yb, zb) = coords[b]
  let (xc, yc, zc) = coords[c]
  let (xd, yd, zd) = coords[d]
  let b1x = xb - xa
  let b1y = yb - ya
  let b1z = zb - za
  let b2x = xc - xb
  let b2y = yc - yb
  let b2z = zc - zb
  let b3x = xd - xc
  let b3y = yd - yc
  let b3z = zd - zc
  // normals
  let n1x = b1y * b2z - b1z * b2y
  let n1y = b1z * b2x - b1x * b2z
  let n1z = b1x * b2y - b1y * b2x
  let n2x = b2y * b3z - b2z * b3y
  let n2y = b2z * b3x - b2x * b3z
  let n2z = b2x * b3y - b2y * b3x
  let n1l = rdkit_sqrt(n1x * n1x + n1y * n1y + n1z * n1z)
  let n2l = rdkit_sqrt(n2x * n2x + n2y * n2y + n2z * n2z)
  if n1l < 0.0000001 || n2l < 0.0000001 {
    return 0.0
  }
  let mut cosv = (n1x * n2x + n1y * n2y + n1z * n2z) / (n1l * n2l)
  if cosv > 1.0 {
    cosv = 1.0
  }
  if cosv < -1.0 {
    cosv = -1.0
  }
  @math.acos(cosv)
}

///|
/// Bond stretch + angle bend + torsion + LJ energy of the given geometry.
pub fn rdkit_ff_energy(
  mol : Molecule,
  coords : Array[(Double, Double, Double)],
) -> Double {
  let n = mol.num_atoms()
  if coords.length() < n {
    return 0.0
  }
  let neighbors = ff_adj(mol)
  let mut e = 0.0
  // bond stretch
  for b in mol.bonds {
    let r = ff_dist(coords, b.from, b.to)
    let r0 = ff_bond_r0(
      mol.atoms[b.from].symbol,
      mol.atoms[b.to].symbol,
      b.order,
    )
    let dr = r - r0
    e = e + ff_bond_kb(b.order) * dr * dr
  }
  // angle bend
  for c = 0; c < n; c = c + 1 {
    let nbs = neighbors[c]
    for i = 0; i < nbs.length(); i = i + 1 {
      for j = i + 1; j < nbs.length(); j = j + 1 {
        let theta = ff_angle(coords, nbs[i], c, nbs[j])
        let theta0 = ff_angle_theta0(mol, c, neighbors)
        let d = theta - theta0
        e = e + 0.12 * d * d
      }
    }
  }
  // torsions over non-terminal bonds
  for b in mol.bonds {
    for na in neighbors[b.from] {
      if na == b.to {
        continue
      }
      for nd in neighbors[b.to] {
        if nd == b.from || nd == na {
          continue
        }
        let phi = ff_torsion(coords, na, b.from, b.to, nd)
        // periodicity 3 for single bonds, barrier for double bonds
        let barrier = if b.order >= 2 { 5.0 } else { 0.12 }
        let period = if b.order >= 2 { 2.0 } else { 3.0 }
        e = e + barrier * (1.0 + @math.cos(period * phi))
      }
    }
  }
  // non-bonded LJ (exclude bonded pairs and 1-3 pairs)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let mut bonded = false
      for nb in neighbors[i] {
        if nb == j {
          bonded = true
        }
      }
      // 1-3 check: share a common neighbor
      let mut onethree = false
      if !bonded {
        for nb in neighbors[i] {
          for nb2 in neighbors[j] {
            if nb == nb2 {
              onethree = true
            }
          }
        }
      }
      if bonded || onethree {
        continue
      }
      let r = ff_dist(coords, i, j)
      if r > 0.1 {
        let sigma = ff_vdw_radius(mol.atoms[i].symbol) +
          ff_vdw_radius(mol.atoms[j].symbol)
        let eps = rdkit_sqrt(
          ff_vdw_epsilon(mol.atoms[i].symbol) *
          ff_vdw_epsilon(mol.atoms[j].symbol),
        )
        let sr6 = sigma / r * (sigma / r) * (sigma / r)
        let sr6 = sr6 * sr6
        e = e + eps * (sr6 * sr6 - 2.0 * sr6)
      }
    }
  }
  e
}

///|
/// Optimize the geometry with steepest descent on a central-difference
/// numerical gradient. Returns (final_energy, optimized_coords).
pub fn rdkit_ff_optimize(
  mol : Molecule,
  coords : Array[(Double, Double, Double)],
  max_iters~ : Int32 = 200,
  tol~ : Double = 0.000001,
) -> (Double, Array[(Double, Double, Double)]) {
  let n = mol.num_atoms()
  let cur : Array[(Double, Double, Double)] = Array(capacity=0)
  for c in coords {
    cur.push(c)
  }
  if n == 0 || coords.length() < n {
    return (rdkit_ff_energy(mol, cur), cur)
  }
  let eps = 0.0001
  let mut step = 0.01
  let mut e_cur = rdkit_ff_energy(mol, cur)
  for _iter = 0; _iter < max_iters; _iter = _iter + 1 {
    // numerical gradient
    let grad : Array[(Double, Double, Double)] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      let (x, y, z) = cur[i]
      // x
      let _ = cur.set(i, (x + eps, y, z))
      let ex = rdkit_ff_energy(mol, cur)
      let _ = cur.set(i, (x - eps, y, z))
      let ex2 = rdkit_ff_energy(mol, cur)
      let gx = (ex - ex2) / (2.0 * eps)
      // y
      let _ = cur.set(i, (x, y + eps, z))
      let ey = rdkit_ff_energy(mol, cur)
      let _ = cur.set(i, (x, y - eps, z))
      let ey2 = rdkit_ff_energy(mol, cur)
      let gy = (ey - ey2) / (2.0 * eps)
      // z
      let _ = cur.set(i, (x, y, z + eps))
      let ez = rdkit_ff_energy(mol, cur)
      let _ = cur.set(i, (x, y, z - eps))
      let ez2 = rdkit_ff_energy(mol, cur)
      let gz = (ez - ez2) / (2.0 * eps)
      let _ = cur.set(i, (x, y, z))
      grad.push((gx, gy, gz))
    }
    // take a step; halve if energy rises
    let mut accepted = false
    let mut local_step = step
    for _attempt = 0; _attempt < 12; _attempt = _attempt + 1 {
      let trial : Array[(Double, Double, Double)] = Array(capacity=n)
      for i = 0; i < n; i = i + 1 {
        let (x, y, z) = cur[i]
        let (gx, gy, gz) = grad[i]
        let nx = x - local_step * gx
        let ny = y - local_step * gy
        let nz = z - local_step * gz
        trial.push((nx, ny, nz))
      }
      let e_trial = rdkit_ff_energy(mol, trial)
      if e_trial < e_cur {
        for i = 0; i < n; i = i + 1 {
          let _ = cur.set(i, trial[i])
        }
        let drop = e_cur - e_trial
        e_cur = e_trial
        step = local_step * 1.2
        accepted = true
        if drop < tol {
          return (e_cur, cur)
        }
        break
      }
      local_step = local_step * 0.5
    }
    if !accepted {
      break
    }
  }
  (e_cur, cur)
}

///|
/// RMSD between two conformers with matching atom counts, after optimal
/// rigid-body superposition (self-contained Kabsch implementation).
pub fn rdkit_conformer_rmsd(a : RdkitConformer, b : RdkitConformer) -> Double {
  aligned_rms(a.coords, b.coords)
}

///|
/// Generate `n_confs` distance-geometry conformers, optionally optimize
/// each with the force field, deduplicate by RMSD threshold and return
/// them sorted by energy (lowest first).
pub fn rdkit_embed_multiple_conformers(
  smiles : String,
  n_confs~ : Int32 = 8,
  rms_threshold~ : Double = 0.5,
  optimize~ : Bool = true,
  seed~ : Int32 = 42,
) -> Array[RdkitConformer] {
  let mol = parse_smiles(smiles)
  let pool : Array[RdkitConformer] = Array(capacity=0)
  for k = 0; k < n_confs; k = k + 1 {
    let raw = rdkit_generate_conformer(smiles, seed=seed + k * 7919)
    let conf = if optimize {
      let (e, coords) = rdkit_ff_optimize(mol, raw.coords)
      RdkitConformer::{ atoms: raw.atoms, coords, energy: e, }
    } else {
      raw
    }
    // deduplicate against already accepted conformers
    let mut dup = false
    for kept in pool {
      if rdkit_conformer_rmsd(conf, kept) < rms_threshold {
        dup = true
      }
    }
    if !dup {
      pool.push(conf)
    }
  }
  // sort by energy ascending (insertion sort, pools are small)
  for i = 1; i < pool.length(); i = i + 1 {
    let key = pool[i]
    let mut j = i - 1
    while j >= 0 && pool[j].energy > key.energy {
      pool[j + 1] = pool[j]
      j = j - 1
    }
    pool[j + 1] = key
  }
  pool
}
