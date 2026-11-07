/// rdkit_cluster_butina.mbt — Butina leader clustering (RDKit gap fill, A4)
///
/// Butina 1999 unweighted leader algorithm: candidates sorted by the number
/// of neighbours within the cutoff, clusters are grown so that every member
/// is within `cutoff` of all other members of the same cluster.

///|
/// Butina clustering over a symmetric distance matrix.
/// Returns clusters as arrays of point indices; the first index of each
/// cluster is its leader (lowest average distance).
pub fn rdkit_butina_cluster(
  dist : Array[Array[Double]],
  cutoff : Double,
) -> Array[Array[Int32]] {
  let n = dist.length()
  let clusters : Array[Array[Int32]] = Array(capacity=0)
  if n == 0 {
    return clusters
  }
  // Neighbour counts within cutoff
  let n_neigh : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let mut c = 0
    for j = 0; j < n; j = j + 1 {
      if j != i && dist[i][j] <= cutoff {
        c = c + 1
      }
    }
    n_neigh[i] = c
  }
  // Candidate order: descending neighbour count, ties by index
  let order : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    order.push(i)
  }
  for a = 0; a < n; a = a + 1 {
    for b = a + 1; b < n; b = b + 1 {
      let ia = order[a]
      let ib = order[b]
      if n_neigh[ib] > n_neigh[ia] || (n_neigh[ib] == n_neigh[ia] && ib < ia) {
        let tmp = order[a]
        order[a] = order[b]
        order[b] = tmp
      }
    }
  }
  let assigned : Array[Int32] = Array::make(n, -1) // cluster id or -1
  let mut next_cluster = 0
  for c = 0; c < n; c = c + 1 {
    let leader = order[c]
    if assigned[leader] != -1 {
      continue
    }
    let members : Array[Int32] = Array(capacity=0)
    members.push(leader)
    assigned[leader] = next_cluster
    // Grow cluster: every member must be within cutoff of all others
    for k = c + 1; k < n; k = k + 1 {
      let cand = order[k]
      if assigned[cand] != -1 {
        continue
      }
      let mut ok = true
      for m in members {
        if dist[cand][m] > cutoff {
          ok = false
        }
      }
      if ok {
        assigned[cand] = next_cluster
        members.push(cand)
      }
    }
    clusters.push(members)
    next_cluster = next_cluster + 1
  }
  clusters
}

///|
/// Tanimoto similarity of two fingerprint hash multisets (min-count
/// intersection over max-count union). The existing rdkit_fp_tanimoto
/// expects bit vectors, which does not match Morgan hash multisets.
fn rdkit_fps_multiset_tanimoto(
  a : Array[Int32],
  b : Array[Int32],
) -> Double {
  let ca : Map[Int32, Int32] = Map([], capacity=0)
  for h in a {
    ca[h] = match ca.get(h) {
      Some(c) => c + 1
      None => 1
    }
  }
  let cb : Map[Int32, Int32] = Map([], capacity=0)
  for h in b {
    cb[h] = match cb.get(h) {
      Some(c) => c + 1
      None => 1
    }
  }
  let mut inter = 0
  for h, cnt in ca {
    match cb.get(h) {
      Some(c2) =>
        if cnt < c2 {
          inter = inter + cnt
        } else {
          inter = inter + c2
        }
      None => ()
    }
  }
  let union = a.length() + b.length() - inter
  if union == 0 {
    0.0
  } else {
    inter.to_double() / union.to_double()
  }
}

///|
/// Convenience wrapper: cluster molecules by fingerprint Tanimoto
/// distance (1 - Tanimoto). Fingerprints are hash multisets as produced
/// by the Morgan fingerprint API.
pub fn rdkit_butina_cluster_from_fps(
  fps : Array[Array[Int32]],
  cutoff : Double,
) -> Array[Array[Int32]] {
  let n = fps.length()
  let dist : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    dist.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let d = 1.0 - rdkit_fps_multiset_tanimoto(fps[i], fps[j])
      dist[i][j] = d
      dist[j][i] = d
    }
  }
  rdkit_butina_cluster(dist, cutoff)
}

///|
/// Cluster molecules directly from SMILES strings using Morgan
/// fingerprints (radius 2, the default cheminformatics choice).
pub fn rdkit_butina_cluster_smiles(
  smiles_list : Array[String],
  cutoff : Double,
  radius~ : Int32 = 2,
) -> Array[Array[Int32]] {
  let fps : Array[Array[Int32]] = Array(capacity=0)
  for smi in smiles_list {
    fps.push(morgan_fingerprint(parse_smiles(smi), radius))
  }
  rdkit_butina_cluster_from_fps(fps, cutoff)
}
