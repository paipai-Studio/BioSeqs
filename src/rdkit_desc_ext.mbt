/// rdkit_desc_ext.mbt — descriptor extensions (RDKit gap fill, item A4)
///
/// Kier-Kappa shape indices, Chi connectivity indices (degree-delta
/// approximation), BCUT2D descriptors (MW / Gasteiger charge extremes),
/// AUTOCORR2D (topological autocorrelation, property = atomic weight),
/// common fr_* fragment counts, and a structured batch API
/// rdkit_calc_descriptors returning Map[String, Double].

///|
/// Atomic weights for common elements (fallback 12.011).
fn rdkit_atomic_weight(symbol : String) -> Double {
  match symbol {
    "H" => 1.008
    "He" => 4.0026
    "Li" => 6.94
    "B" => 10.81
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "F" => 18.998
    "Na" => 22.99
    "Mg" => 24.305
    "Al" => 26.982
    "Si" => 28.085
    "P" => 30.974
    "S" => 32.06
    "Cl" => 35.45
    "K" => 39.098
    "Ca" => 40.078
    "Fe" => 55.845
    "Ni" => 58.693
    "Cu" => 63.546
    "Zn" => 65.38
    "Br" => 79.904
    "Se" => 78.97
    "Ag" => 107.87
    "I" => 126.9
    "Sn" => 118.71
    _ => 12.011
  }
}

///|
fn rdkit_desc_mol_weights(mol : Molecule) -> Array[Double] {
  let out : Array[Double] = Array(capacity=mol.num_atoms())
  for a in mol.atoms {
    out.push(rdkit_atomic_weight(a.symbol))
  }
  out
}

///|
/// Count simple paths of a given length (no repeated atoms), counted once
/// per undirected path.
fn rdkit_desc_count_paths(
  neighbors : Array[Array[Int32]],
  length : Int32,
) -> Int32 {
  let n = neighbors.length()
  if length < 1 || n == 0 {
    return 0
  }
  let mut total = 0
  for start = 0; start < n; start = start + 1 {
    total = total +
      rdkit_desc_path_dfs_impl(
        start,
        start,
        length,
        neighbors,
        Array::make(n, false),
      )
  }
  total / 2
}

///|
fn rdkit_desc_path_dfs_impl(
  start : Int32,
  cur : Int32,
  remaining : Int32,
  neighbors : Array[Array[Int32]],
  visited : Array[Bool],
) -> Int32 {
  if remaining == 0 {
    return 1
  }
  visited[cur] = true
  let mut count = 0
  for nb in neighbors[cur] {
    if !visited[nb] {
      count = count +
        rdkit_desc_path_dfs_impl(start, nb, remaining - 1, neighbors, visited)
    }
  }
  visited[cur] = false
  count
}

///|
/// Kier kappa shape indices (no alpha correction; degree-based path counts).
/// kappa1 = N(N-1)^2/P1^2, kappa2 = (N-1)(N-2)^2/P2^2,
/// kappa3 = (N-2)(N-3)^2/P3^2 (0 when N < required size).
pub fn rdkit_kappa_indices(mol : Molecule) -> (Double, Double, Double) {
  let n = mol.num_atoms()
  if n < 2 {
    return (0.0, 0.0, 0.0)
  }
  let neighbors = rdkit_desc_adj(mol)
  let p1 = mol.num_bonds().to_double()
  let p2 = rdkit_desc_count_paths(neighbors, 2).to_double()
  let p3 = rdkit_desc_count_paths(neighbors, 3).to_double()
  let nf = n.to_double()
  let k1 = if p1 > 0.0 { nf * (nf - 1.0) * (nf - 1.0) / (p1 * p1) } else { 0.0 }
  let k2 = if p2 > 0.0 && n >= 3 {
    (nf - 1.0) * (nf - 2.0) * (nf - 2.0) / (p2 * p2)
  } else {
    0.0
  }
  let k3 = if p3 > 0.0 && n >= 4 {
    (nf - 2.0) * (nf - 3.0) * (nf - 3.0) / (p3 * p3)
  } else {
    0.0
  }
  (k1, k2, k3)
}

///|
fn rdkit_desc_adj(mol : Molecule) -> Array[Array[Int32]] {
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
/// Chi connectivity indices using degree deltas (approximation of the
/// valence-delta values): Chi_m = sum over m-step simple paths of
/// 1/sqrt(prod of endpoint-and-inner degrees).
pub fn rdkit_chi_indices(
  mol : Molecule,
) -> (Double, Double, Double, Double, Double) {
  let n = mol.num_atoms()
  if n == 0 {
    return (0.0, 0.0, 0.0, 0.0, 0.0)
  }
  let neighbors = rdkit_desc_adj(mol)
  let deg : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    deg.push(neighbors[i].length().to_double())
  }
  let mut chi0 = 0.0
  for i = 0; i < n; i = i + 1 {
    chi0 = chi0 + 1.0 / rdkit_desc_sqrt(deg[i])
  }
  // Path-based terms for m = 1..4
  let sums : Array[Double] = Array::make(5, 0.0)
  for start = 0; start < n; start = start + 1 {
    let path : Array[Int32] = Array(capacity=0)
    path.push(start)
    rdkit_desc_chi_dfs(neighbors, deg, start, path, 4, sums)
  }
  // sums[m] counts each undirected path twice (once per endpoint) for m >= 1
  let chi1 = sums[1] / 2.0
  let chi2 = sums[2] / 2.0
  let chi3 = sums[3] / 2.0
  let chi4 = sums[4] / 2.0
  (chi0, chi1, chi2, chi3, chi4)
}

///|
fn rdkit_desc_chi_dfs(
  neighbors : Array[Array[Int32]],
  deg : Array[Double],
  cur : Int32,
  path : Array[Int32],
  remaining : Int32,
  sums : Array[Double],
) -> Unit {
  if remaining == 0 {
    return
  }
  for nb in neighbors[cur] {
    let mut seen = false
    for a in path {
      if a == nb {
        seen = true
      }
    }
    if seen {
      continue
    }
    let mut prod = 1.0
    for a in path {
      prod = prod * deg[a]
    }
    prod = prod * deg[nb]
    if prod > 0.0 {
      let m = path.length()
      if m <= 4 {
        sums[m] = sums[m] + 1.0 / rdkit_desc_sqrt(prod)
      }
    }
    path.push(nb)
    rdkit_desc_chi_dfs(neighbors, deg, nb, path, remaining - 1, sums)
    let _ = path.pop()
  }
}

///|
fn rdkit_desc_sqrt(x : Double) -> Double {
  @math.pow(x, 0.5)
}

///|
/// Jacobi eigenvalue solver for small symmetric matrices: returns all
/// eigenvalues (unsorted) after iterative Givens rotations.
fn rdkit_desc_jacobi_eigenvalues(
  mat : Array[Array[Double]],
  sweeps~ : Int32 = 60,
) -> Array[Double] {
  let n = mat.length()
  let a : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = mat[i][j]
    }
    a.push(row)
  }
  for _sweep = 0; _sweep < sweeps; _sweep = _sweep + 1 {
    // off-diagonal magnitude
    let mut off = 0.0
    for p = 0; p < n; p = p + 1 {
      for q = p + 1; q < n; q = q + 1 {
        off = off + rdkit_desc_abs(a[p][q] * a[p][q])
      }
    }
    if off < 0.000000000001 {
      break
    }
    for p = 0; p < n; p = p + 1 {
      for q = p + 1; q < n; q = q + 1 {
        let apq = a[p][q]
        if rdkit_desc_abs(apq * apq) > 0.0000000000001 {
          let theta = (a[q][q] - a[p][p]) / (2.0 * apq)
          let root = rdkit_desc_sqrt(1.0 + theta * theta)
          let t = if theta >= 0.0 {
            1.0 / (theta + root)
          } else {
            1.0 / (theta - root)
          }
          let c = 1.0 / rdkit_desc_sqrt(1.0 + t * t)
          let s = t * c
          for k = 0; k < n; k = k + 1 {
            let akp = a[k][p]
            let akq = a[k][q]
            a[k][p] = c * akp - s * akq
            a[k][q] = s * akp + c * akq
          }
          for k = 0; k < n; k = k + 1 {
            let apk = a[p][k]
            let aqk = a[q][k]
            a[p][k] = c * apk - s * aqk
            a[q][k] = s * apk + c * aqk
          }
        }
      }
    }
  }
  let out : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    out.push(a[i][i])
  }
  out
}

///|
fn rdkit_desc_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
/// Jacobi eigenvalue solver that also accumulates eigenvectors: returns
/// (eigenvalues, eigenvector matrix V) with mat = V * diag(vals) * V^T
/// (V's column k is the eigenvector for vals[k]).
fn rdkit_desc_jacobi_full(
  mat : Array[Array[Double]],
  sweeps~ : Int32 = 60,
) -> (Array[Double], Array[Array[Double]]) {
  let n = mat.length()
  let a : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = mat[i][j]
    }
    a.push(row)
  }
  let v : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    row[i] = 1.0
    v.push(row)
  }
  for _sweep = 0; _sweep < sweeps; _sweep = _sweep + 1 {
    let mut off = 0.0
    for p = 0; p < n; p = p + 1 {
      for q = p + 1; q < n; q = q + 1 {
        off = off + rdkit_desc_abs(a[p][q] * a[p][q])
      }
    }
    if off < 0.000000000001 {
      break
    }
    for p = 0; p < n; p = p + 1 {
      for q = p + 1; q < n; q = q + 1 {
        let apq = a[p][q]
        if rdkit_desc_abs(apq * apq) > 0.0000000000001 {
          let theta = (a[q][q] - a[p][p]) / (2.0 * apq)
          let root = rdkit_desc_sqrt(1.0 + theta * theta)
          let t = if theta >= 0.0 {
            1.0 / (theta + root)
          } else {
            1.0 / (theta - root)
          }
          let c = 1.0 / rdkit_desc_sqrt(1.0 + t * t)
          let s = t * c
          for k = 0; k < n; k = k + 1 {
            let akp = a[k][p]
            let akq = a[k][q]
            a[k][p] = c * akp - s * akq
            a[k][q] = s * akp + c * akq
          }
          for k = 0; k < n; k = k + 1 {
            let apk = a[p][k]
            let aqk = a[q][k]
            a[p][k] = c * apk - s * aqk
            a[q][k] = s * apk + c * aqk
          }
          for k = 0; k < n; k = k + 1 {
            let vkp = v[k][p]
            let vkq = v[k][q]
            v[k][p] = c * vkp - s * vkq
            v[k][q] = s * vkp + c * vkq
          }
        }
      }
    }
  }
  let vals : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    vals.push(a[i][i])
  }
  (vals, v)
}

///|
/// Build a Burden-style matrix: diagonal = per-atom property, bonded
/// off-diagonal entries = 1.0.
fn rdkit_desc_burden(
  mol : Molecule,
  props : Array[Double],
) -> Array[Array[Double]] {
  let n = mol.num_atoms()
  let mat : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    row[i] = props[i]
    mat[i] = row
  }
  for b in mol.bonds {
    if b.from < n && b.to < n {
      mat[b.from][b.to] = 1.0
      mat[b.to][b.from] = 1.0
    }
  }
  mat
}

///|
fn rdkit_desc_bcut_from_props(
  mol : Molecule,
  props : Array[Double],
) -> (Double, Double) {
  let n = mol.num_atoms()
  if n == 0 {
    return (0.0, 0.0)
  }
  let mat = rdkit_desc_burden(mol, props)
  let eigs = rdkit_desc_jacobi_eigenvalues(mat)
  let mut hi = eigs[0]
  let mut lo = eigs[0]
  for i = 1; i < eigs.length(); i = i + 1 {
    if eigs[i] > hi {
      hi = eigs[i]
    }
    if eigs[i] < lo {
      lo = eigs[i]
    }
  }
  (hi, lo)
}

///|
/// BCUT2D descriptors (4 of the 8 RDKit values): MW and Gasteiger charge
/// extremes of the Burden matrix.
pub fn rdkit_bcut2d(mol : Molecule) -> (Double, Double, Double, Double) {
  let (mw_hi, mw_lo) = rdkit_desc_bcut_from_props(
    mol,
    rdkit_desc_mol_weights(mol),
  )
  let charges = rdkit_gasteiger_charges(mol)
  let (ch_hi, ch_lo) = rdkit_desc_bcut_from_props(mol, charges)
  (mw_hi, mw_lo, ch_hi, ch_lo)
}

///|
/// AUTOCORR2D (Broto-Moreau style): ATSt = sum over atom pairs at
/// topological distance t of p_i * p_j, for t = 1..8 (property =
/// atomic weight, the ATSm-like variant).
pub fn rdkit_autocorr2d(mol : Molecule, max_lag~ : Int32 = 8) -> Array[Double] {
  let n = mol.num_atoms()
  let out : Array[Double] = Array::make(max_lag + 1, 0.0)
  if n == 0 {
    return out
  }
  let weights = rdkit_desc_mol_weights(mol)
  let neighbors = rdkit_desc_adj(mol)
  for start = 0; start < n; start = start + 1 {
    // BFS distances up to max_lag
    let dist = Array::make(n, -1)
    dist[start] = 0
    let queue : Array[Int32] = Array(capacity=0)
    queue.push(start)
    let mut head = 0
    while head < queue.length() {
      let cur = queue[head]
      head = head + 1
      if dist[cur] >= max_lag {
        continue
      }
      for nb in neighbors[cur] {
        if dist[nb] == -1 {
          dist[nb] = dist[cur] + 1
          queue.push(nb)
        }
      }
    }
    for other = start + 1; other < n; other = other + 1 {
      let d = dist[other]
      if d >= 1 && d <= max_lag {
        out[d] = out[d] + weights[start] * weights[other]
      }
    }
  }
  out
}

///|
/// Common fr_* fragment counts (approximations on the bond-order graph):
/// amide, ester, ether, ketone, nitrile, nitro, urea, guanidine,
/// aldehyde (C=O bonds), alkene (C=C bonds), halogen atoms, benzene
/// (aromatic six-rings).
pub fn rdkit_fragment_counts(mol : Molecule) -> Map[String, Int32] {
  let out : Map[String, Int32] = Map([], capacity=16)
  let smiles = mol_to_smiles(mol)
  out["fr_amide"] = rdkit_get_substruct_matches(smiles, "C(=O)N").length()
  out["fr_ester"] = rdkit_get_substruct_matches(smiles, "C(=O)O").length()
  out["fr_ether"] = rdkit_get_substruct_matches(smiles, "COC").length()
  out["fr_ketone"] = rdkit_get_substruct_matches(smiles, "CC(=O)C").length()
  out["fr_nitrile"] = rdkit_get_substruct_matches(smiles, "C#N").length()
  out["fr_nitro"] = rdkit_get_substruct_matches(smiles, "CN(=O)=O").length()
  out["fr_urea"] = rdkit_get_substruct_matches(smiles, "NC(=O)N").length()
  out["fr_guanidine"] = rdkit_get_substruct_matches(smiles, "NC(=N)N").length()
  out["fr_aldehyde"] = rdkit_get_substruct_matches(smiles, "C=O").length()
  out["fr_alkene"] = rdkit_get_substruct_matches(smiles, "C=C").length()
  let mut halogens = 0
  for a in mol.atoms {
    match a.symbol {
      "F" | "Cl" | "Br" | "I" => halogens = halogens + 1
      _ => ()
    }
  }
  out["fr_halogen"] = halogens
  let aromatic_mol = mol_perceive_aromaticity(mol)
  let mut benzene = 0
  for ring in mol_find_rings(aromatic_mol) {
    if ring.length() == 6 && mol_is_aromatic_ring(aromatic_mol, ring) {
      benzene = benzene + 1
    }
  }
  out["fr_benzene"] = benzene
  out
}

///|
/// Structured batch descriptor API (RDKit CalcMolDescriptors mode):
/// baseline descriptors + Kappa + Chi + BCUT2D + AUTOCORR2D + fr_* counts
/// in a single Map.
pub fn rdkit_calc_descriptors(mol : Molecule) -> Map[String, Double] {
  let out : Map[String, Double] = Map([], capacity=64)
  out["NumAtoms"] = mol.num_atoms().to_double()
  out["NumBonds"] = mol.num_bonds().to_double()
  out["HeavyAtoms"] = mol_heavy_atom_count(mol).to_double()
  out["MolWt"] = mol_molecular_weight(mol)
  out["LogP"] = mol_logp(mol)
  out["TPSA"] = mol_tpsa(mol)
  out["HBD"] = mol_hbd_count(mol).to_double()
  out["HBA"] = mol_hba_count(mol).to_double()
  out["RotatableBonds"] = mol_rotatable_bonds(mol).to_double()
  out["FormalCharge"] = mol_formal_charge(mol).to_double()
  out["FractionCSP3"] = mol_fraction_sp3(mol)
  out["RingCount"] = mol_ring_count(mol).to_double()
  out["AromaticRings"] = mol_aromatic_ring_count(mol).to_double()
  let (k1, k2, k3) = rdkit_kappa_indices(mol)
  out["Kappa1"] = k1
  out["Kappa2"] = k2
  out["Kappa3"] = k3
  let (c0, c1, c2, c3, c4) = rdkit_chi_indices(mol)
  out["Chi0"] = c0
  out["Chi1"] = c1
  out["Chi2"] = c2
  out["Chi3"] = c3
  out["Chi4"] = c4
  let (mw_hi, mw_lo, ch_hi, ch_lo) = rdkit_bcut2d(mol)
  out["BCUT2D_MWHI"] = mw_hi
  out["BCUT2D_MWLOW"] = mw_lo
  out["BCUT2D_CHGHI"] = ch_hi
  out["BCUT2D_CHGLOW"] = ch_lo
  let ats = rdkit_autocorr2d(mol)
  for t = 1; t <= 8; t = t + 1 {
    out["ATS" + t.to_string()] = ats[t]
  }
  for key, v in rdkit_fragment_counts(mol) {
    out[key] = v.to_double()
  }
  out
}
