/// Topological indices: Wiener, Balaban J, Randić connectivity.
///
/// RDKit equivalent: `rdkit.Chem.GraphDescriptors`
///
/// These are graph-theoretic invariants computed from the molecular graph.
/// They capture branching, compactness, and connectivity information.

///|
/// Compute the Wiener index: sum of all distances in the distance matrix
/// (upper triangle only).
/// W = sum_{i<j} d(i,j)
pub fn rdkit_wiener_index(mol : Molecule) -> Double {
  let dm = mol_distance_matrix(mol)
  let n = dm.length()
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      sum = sum + dm[i][j].to_double()
    }
  }
  sum
}

///|
/// Compute the Balaban J index.
/// J = (m / (m - n + 1)) * sum_{i<j, bond(i,j)} 1 / sqrt(d_i * d_j)
/// where d_i = sum of distances from atom i to all others (distance degree),
/// m = number of bonds, n = number of atoms.
pub fn rdkit_balaban_j(mol : Molecule) -> Double {
  let dm = mol_distance_matrix(mol)
  let n = mol.num_atoms()
  let m = mol.num_bonds()
  if n <= 1 || m == 0 {
    return 0.0
  }
  // Compute distance degrees
  let dist_deg : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n; j = j + 1 {
      if i != j {
        s = s + dm[i][j].to_double()
      }
    }
    dist_deg[i] = s
  }
  // Sum over bonds
  let mut edge_sum = 0.0
  for bond in mol.bonds {
    let a = bond.from
    let b = bond.to
    if dist_deg[a] > 0.0 && dist_deg[b] > 0.0 {
      edge_sum = edge_sum + 1.0 / (dist_deg[a] * dist_deg[b]).sqrt()
    }
  }
  let cycle_rank = m - n + 1
  if cycle_rank == 0 {
    return 0.0
  }
  m.to_double() / cycle_rank.to_double() * edge_sum
}

///|
/// Compute the Randić connectivity index (molecular connectivity index chi_1).
/// chi_1 = sum over bonds 1/sqrt(d_i * d_j)
/// where d_i = degree of atom i (number of bonded neighbors).
pub fn rdkit_randic_index(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  if n <= 1 {
    return 0.0
  }
  // Compute degrees
  let degree : Array[Int32] = Array::make(n, 0)
  for bond in mol.bonds {
    degree[bond.from] = degree[bond.from] + 1
    degree[bond.to] = degree[bond.to] + 1
  }
  // Sum over bonds
  let mut sum = 0.0
  for bond in mol.bonds {
    let di = degree[bond.from].to_double()
    let dj = degree[bond.to].to_double()
    if di > 0.0 && dj > 0.0 {
      sum = sum + 1.0 / (di * dj).sqrt()
    }
  }
  sum
}

///|
/// Compute the Platt index: sum of edge degrees (sum of all d_i for each bond).
/// = sum over bonds (degree[i] + degree[j])
pub fn rdkit_platt_index(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  if n <= 1 {
    return 0.0
  }
  let degree : Array[Int32] = Array::make(n, 0)
  for bond in mol.bonds {
    degree[bond.from] = degree[bond.from] + 1
    degree[bond.to] = degree[bond.to] + 1
  }
  let mut sum = 0.0
  for bond in mol.bonds {
    sum = sum + (degree[bond.from] + degree[bond.to]).to_double()
  }
  sum
}

///|
/// Compute the Zagreb M1 index: sum of squared vertex degrees.
/// M1 = sum_i d_i^2
pub fn rdkit_zagreb_m1(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  if n == 0 {
    return 0.0
  }
  let degree : Array[Int32] = Array::make(n, 0)
  for bond in mol.bonds {
    degree[bond.from] = degree[bond.from] + 1
    degree[bond.to] = degree[bond.to] + 1
  }
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let d = degree[i].to_double()
    sum = sum + d * d
  }
  sum
}

///|
/// Compute the Zagreb M2 index: sum over bonds of d_i * d_j.
/// M2 = sum_{bond(i,j)} d_i * d_j
pub fn rdkit_zagreb_m2(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  if n <= 1 {
    return 0.0
  }
  let degree : Array[Int32] = Array::make(n, 0)
  for bond in mol.bonds {
    degree[bond.from] = degree[bond.from] + 1
    degree[bond.to] = degree[bond.to] + 1
  }
  let mut sum = 0.0
  for bond in mol.bonds {
    sum = sum + (degree[bond.from] * degree[bond.to]).to_double()
  }
  sum
}

///|
/// All topological indices in one report.
pub struct RdkitTopoIndices {
  wiener : Double
  balaban_j : Double
  randic : Double
  platt : Double
  zagreb_m1 : Double
  zagreb_m2 : Double
}

///|
pub fn RdkitTopoIndices::wiener(self : RdkitTopoIndices) -> Double {
  self.wiener
}

///|
pub fn RdkitTopoIndices::balaban_j(self : RdkitTopoIndices) -> Double {
  self.balaban_j
}

///|
pub fn RdkitTopoIndices::randic(self : RdkitTopoIndices) -> Double {
  self.randic
}

///|
pub fn RdkitTopoIndices::platt(self : RdkitTopoIndices) -> Double {
  self.platt
}

///|
pub fn RdkitTopoIndices::zagreb_m1(self : RdkitTopoIndices) -> Double {
  self.zagreb_m1
}

///|
pub fn RdkitTopoIndices::zagreb_m2(self : RdkitTopoIndices) -> Double {
  self.zagreb_m2
}

///|
/// Compute all topological indices at once.
pub fn rdkit_topo_indices(mol : Molecule) -> RdkitTopoIndices {
  {
    wiener: rdkit_wiener_index(mol),
    balaban_j: rdkit_balaban_j(mol),
    randic: rdkit_randic_index(mol),
    platt: rdkit_platt_index(mol),
    zagreb_m1: rdkit_zagreb_m1(mol),
    zagreb_m2: rdkit_zagreb_m2(mol),
  }
}
