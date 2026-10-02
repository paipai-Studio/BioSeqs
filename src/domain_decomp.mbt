///|
/// Domain decomposition using slow normal modes.
///
/// Residues are clustered based on the sign patterns of the slowest
/// eigenvectors. Residues with similar sign patterns move together
/// as rigid domains.

///|
/// A structural domain.
pub struct Domain {
  id : Int
  residues : Array[Int]
  size : Int
}

///|
/// Domain decomposition result.
pub struct DomainDecomposition {
  n : Int
  num_domains : Int
  domains : Array[Domain]
  // Rigidity profile: how rigid each residue is [0, 1]
  rigidity : Array[Double]
}

///|
/// Decompose a structure into domains using GNM slow modes.
///
/// Strategy: use the slowest k eigenvectors. For each residue,
/// build a feature vector from the signs of the k slowest modes.
/// Cluster residues with identical sign patterns.
pub fn DomainDecomposition::from_gnm(
  gnm : GNM,
  num_slow_modes : Int,
) -> DomainDecomposition {
  let n = gnm.num_nodes()
  let k = if num_slow_modes > 0 && num_slow_modes < n {
    num_slow_modes
  } else {
    2
  }
  let modes = gnm.slow_modes(k)
  // Build sign patterns
  let sign_patterns : Array[Array[Int]] = Array::makei(n, { _ => Array::make(k, 0) })
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < k && j < modes.length(); j = j + 1 {
      let (_, _, evec) = modes[j]
      if evec[i] > 0.01 {
        sign_patterns[i][j] = 1
      } else if evec[i] < -0.01 {
        sign_patterns[i][j] = -1
      } else {
        sign_patterns[i][j] = 0
      }
    }
  }
  // Cluster by sign pattern
  let domain_id = Array::make(n, -1)
  let domains : Array[Domain] = []
  for i = 0; i < n; i = i + 1 {
    if domain_id[i] >= 0 {
      continue
    }
    // Find all residues with same pattern
    let residues : Array[Int] = []
    for j = i; j < n; j = j + 1 {
      if domain_id[j] >= 0 {
        continue
      }
      let mut same = true
      for m = 0; m < k; m = m + 1 {
        if sign_patterns[j][m] != sign_patterns[i][m] {
          same = false
          break
        }
      }
      if same {
        domain_id[j] = domains.length()
        residues.push(j)
      }
    }
    domains.push(Domain::{ id: domains.length(), residues, size: residues.length() })
  }
  // Rigidity profile from B-factors (inverse: low B-factor = rigid)
  let bfactors = gnm.predict_bfactors()
  let mut max_b = 0.0
  let mut min_b = 1.0e30
  for i = 0; i < n; i = i + 1 {
    if bfactors[i] > max_b {
      max_b = bfactors[i]
    }
    if bfactors[i] < min_b {
      min_b = bfactors[i]
    }
  }
  let rigidity = Array::make(n, 0.0)
  let range = max_b - min_b
  for i = 0; i < n; i = i + 1 {
    if range > 1.0e-10 {
      rigidity[i] = 1.0 - (bfactors[i] - min_b) / range
    } else {
      rigidity[i] = 0.5
    }
  }
  DomainDecomposition::{ n, num_domains: domains.length(), domains, rigidity }
}

///|
/// Decompose using ANM slow modes.
/// Projects 3N-dimensional eigenvectors to N-dimensional per-residue
/// amplitudes before sign clustering.
pub fn DomainDecomposition::from_anm(
  anm : ANM,
  num_slow_modes : Int,
) -> DomainDecomposition {
  let n = anm.num_nodes()
  let k = if num_slow_modes > 0 && num_slow_modes < n {
    num_slow_modes
  } else {
    2
  }
  let modes = anm.slow_modes(k)
  // Build per-residue amplitudes from 3N eigenvectors
  let amplitudes : Array[Array[Double]] = Array::makei(n, { _ => Array::make(k, 0.0) })
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < k && j < modes.length(); j = j + 1 {
      let (_, _, evec) = modes[j]
      let mut s = 0.0
      for a = 0; a < 3; a = a + 1 {
        s = s + evec[3 * i + a] * evec[3 * i + a]
      }
      amplitudes[i][j] = s.sqrt()
    }
  }
  // Sign patterns based on amplitude magnitude
  let sign_patterns : Array[Array[Int]] = Array::makei(n, { _ => Array::make(k, 0) })
  for j = 0; j < k && j < modes.length(); j = j + 1 {
    let (_, _, evec) = modes[j]
    for i = 0; i < n; i = i + 1 {
      let mut s = 0.0
      for a = 0; a < 3; a = a + 1 {
        s = s + evec[3 * i + a]
      }
      if s > 0.01 {
        sign_patterns[i][j] = 1
      } else if s < -0.01 {
        sign_patterns[i][j] = -1
      }
    }
  }
  // Cluster
  let domain_id = Array::make(n, -1)
  let domains : Array[Domain] = []
  for i = 0; i < n; i = i + 1 {
    if domain_id[i] >= 0 {
      continue
    }
    let residues : Array[Int] = []
    for j = i; j < n; j = j + 1 {
      if domain_id[j] >= 0 {
        continue
      }
      let mut same = true
      for m = 0; m < k; m = m + 1 {
        if sign_patterns[j][m] != sign_patterns[i][m] {
          same = false
          break
        }
      }
      if same {
        domain_id[j] = domains.length()
        residues.push(j)
      }
    }
    domains.push(Domain::{ id: domains.length(), residues, size: residues.length() })
  }
  // Rigidity from fluctuations
  let fluctuations = anm.fluctuations()
  let mut max_f = 0.0
  let mut min_f = 1.0e30
  for i = 0; i < n; i = i + 1 {
    if fluctuations[i] > max_f {
      max_f = fluctuations[i]
    }
    if fluctuations[i] < min_f {
      min_f = fluctuations[i]
    }
  }
  let rigidity = Array::make(n, 0.0)
  let range = max_f - min_f
  for i = 0; i < n; i = i + 1 {
    if range > 1.0e-10 {
      rigidity[i] = 1.0 - (fluctuations[i] - min_f) / range
    } else {
      rigidity[i] = 0.5
    }
  }
  DomainDecomposition::{ n, num_domains: domains.length(), domains, rigidity }
}

///|
/// Get all domains.
pub fn DomainDecomposition::domains(self : DomainDecomposition) -> Array[Domain] {
  self.domains
}

///|
/// Get domain boundaries: residue indices where domain changes.
pub fn DomainDecomposition::domain_boundaries(
  self : DomainDecomposition,
) -> Array[Int] {
  let boundaries : Array[Int] = []
  // Build residue → domain map
  let res_domain = Array::make(self.n, -1)
  for d = 0; d < self.domains.length(); d = d + 1 {
    for r in self.domains[d].residues {
      res_domain[r] = d
    }
  }
  for i = 1; i < self.n; i = i + 1 {
    if res_domain[i] != res_domain[i - 1] && res_domain[i] >= 0 && res_domain[i - 1] >= 0 {
      boundaries.push(i)
    }
  }
  boundaries
}

///|
/// Number of domains.
pub fn DomainDecomposition::num_domains(self : DomainDecomposition) -> Int {
  self.num_domains
}

///|
/// Rigidity profile.
pub fn DomainDecomposition::rigidity_profile(
  self : DomainDecomposition,
) -> Array[Double] {
  self.rigidity
}

///|
/// Get domain ID for a residue.
pub fn DomainDecomposition::domain_of(
  self : DomainDecomposition,
  residue_idx : Int,
) -> Int {
  for d = 0; d < self.domains.length(); d = d + 1 {
    for r in self.domains[d].residues {
      if r == residue_idx {
        return d
      }
    }
  }
  -1
}

///|
/// Largest domain size.
pub fn DomainDecomposition::largest_domain_size(
  self : DomainDecomposition,
) -> Int {
  let mut max_size = 0
  for d = 0; d < self.domains.length(); d = d + 1 {
    if self.domains[d].size > max_size {
      max_size = self.domains[d].size
    }
  }
  max_size
}

///|
/// Domain size distribution.
pub fn DomainDecomposition::domain_sizes(self : DomainDecomposition) -> Array[Int] {
  let sizes = Array::make(self.domains.length(), 0)
  for d = 0; d < self.domains.length(); d = d + 1 {
    sizes[d] = self.domains[d].size
  }
  sizes
}
