///|
/// Interaction graph / domain definition — PICIO + SCADIO layer.
///
/// PICIO (Protein Interaction Classifier I/O):
///   Chain-chain interface analysis — which residues form the
///   interface between two chains, interface area (ΔSASA),
///   and per-residue contact pairs across the interface.
///
/// SCADIO (Structure-Aware Domain Identification):
///   Decompose a structure into structural domains using
///   label-propagation community detection on the Cα contact graph.
///
/// Both layers build on top of existing primitives:
///   contact_map / residue_contact_map / inter_contact_pairs /
///   calc_sasa / get_residue_starts.

// ============================================================================
// Common helpers — extract per-chain Cα coordinates
// ============================================================================

fn ca_per_residue(
  arr : AtomArray,
) -> (Array[Vector3], Array[Char], Array[String]) {
  let starts = arr.get_residue_starts()
  let n_res = starts.length()
  let cas : Array[Vector3] = Array(capacity=n_res)
  let chains : Array[Char] = Array(capacity=n_res)
  let names : Array[String] = Array(capacity=n_res)
  let mut ri = 0
  while ri < n_res {
    let rs = starts[ri]
    let re = if ri + 1 < n_res { starts[ri + 1] } else { arr.len() }
    let mut ca_idx = -1
    let mut j = rs
    while j < re {
      if arr.names[j] == "CA" {
        ca_idx = j
        break
      }
      j = j + 1
    }
    if ca_idx >= 0 {
      cas.push(arr.coords[ca_idx])
    } else {
      cas.push(arr.coords[rs])
    }
    chains.push(arr.chain_ids[rs])
    names.push(arr.res_names[rs])
    ri = ri + 1
  }
  (cas, chains, names)
}

///|
/// Filter AtomArray to atoms belonging to a given chain id.
fn filter_by_chain(arr : AtomArray, chain : Char) -> AtomArray {
  let mask : Array[Bool] = Array(capacity=arr.len())
  let mut i = 0
  while i < arr.len() {
    mask.push(arr.chain_ids[i] == chain)
    i = i + 1
  }
  arr.filter(mask)
}

// ============================================================================
// PICIO — Chain-chain interface analysis
// ============================================================================

///|
/// One chain-chain interface summary.
pub struct ChainInterface {
  chain_a : Char
  chain_b : Char
  interface_area : Double
  n_residues_a : Int32
  n_residues_b : Int32
  n_contacts : Int32
  residues_a : Array[Int32]
  residues_b : Array[Int32]
}

///|
/// Find all residues on `chain_a` within `cutoff` Å of any residue on
/// `chain_b` (and vice versa). Returns two parallel arrays of global
/// residue indices.
pub fn AtomArray::chain_interface_residues(
  self : AtomArray,
  chain_a : Char,
  chain_b : Char,
  cutoff~ : Double = 5.0,
) -> (Array[Int32], Array[Int32]) {
  // cutoff~ is directly accessible as cutoff
  let (cas, chains, _) = ca_per_residue(self)
  let n = cas.length()
  let res_a : Array[Int32] = Array(capacity=n / 4)
  let res_b : Array[Int32] = Array(capacity=n / 4)
  let mut i = 0
  while i < n {
    if chains[i] != chain_a && chains[i] != chain_b {
      i = i + 1
      continue
    }
    let is_a = chains[i] == chain_a
    let mut in_interface = false
    let mut j = 0
    while j < n && !in_interface {
      if is_a {
        if chains[j] == chain_b {
          if cas[i].distance(cas[j]) <= cutoff {
            in_interface = true
          }
        }
      } else if chains[j] == chain_a {
        if cas[i].distance(cas[j]) <= cutoff {
          in_interface = true
        }
      }
      j = j + 1
    }
    if in_interface {
      if is_a {
        res_a.push(i)
      } else {
        res_b.push(i)
      }
    }
    i = i + 1
  }
  (res_a, res_b)
}

///|
/// Compute interface area between two chains using the ΔSASA method:
///   ΔSASA = SASA(chain_a) + SASA(chain_b) - SASA(chain_a + chain_b)
/// Standard buried-surface-area measure for interaction interfaces.
pub fn AtomArray::chain_interface_area(
  self : AtomArray,
  chain_a : Char,
  chain_b : Char,
) -> Double {
  let arr_a = filter_by_chain(self, chain_a)
  let arr_b = filter_by_chain(self, chain_b)
  if arr_a.len() == 0 || arr_b.len() == 0 {
    return 0.0
  }
  let sasa_a = sum_doubles(arr_a.calc_sasa())
  let sasa_b = sum_doubles(arr_b.calc_sasa())
  let sasa_complex = sum_doubles(self.calc_sasa())
  let delta = sasa_a + sasa_b - sasa_complex
  if delta < 0.0 {
    0.0
  } else {
    delta
  }
}

///|
/// Return all (res_i, res_j) pairs between chain_a and chain_b
/// whose Cα atoms are within `cutoff` Å.
pub fn AtomArray::chain_interaction_pairs(
  self : AtomArray,
  chain_a : Char,
  chain_b : Char,
  cutoff~ : Double = 5.0,
) -> Array[(Int32, Int32)] {
  // cutoff~ is directly accessible as cutoff
  let (cas, chains, _) = ca_per_residue(self)
  let n = cas.length()
  let pairs : Array[(Int32, Int32)] = Array(capacity=n / 2)
  let mut i = 0
  while i < n {
    if chains[i] != chain_a {
      i = i + 1
      continue
    }
    let mut j = 0
    while j < n {
      if chains[j] == chain_b {
        if cas[i].distance(cas[j]) <= cutoff {
          pairs.push((i, j))
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  pairs
}

///|
/// Analyze all chain-chain interfaces in this structure.
/// Returns one ChainInterface per unique (chain_a, chain_b) pair
/// that has at least one interface residue on each side.
pub fn AtomArray::all_chain_interfaces(
  self : AtomArray,
  cutoff~ : Double = 5.0,
) -> Array[ChainInterface] {
  // cutoff~ is directly accessible as cutoff
  let interfaces : Array[ChainInterface] = Array(capacity=4)
  let chain_set : Array[Char] = Array(capacity=8)
  let starts = self.get_residue_starts()
  let mut i = 0
  while i < starts.length() {
    let c = self.chain_ids[starts[i]]
    let mut found = false
    let mut j = 0
    while j < chain_set.length() {
      if chain_set[j] == c {
        found = true
        break
      }
      j = j + 1
    }
    if !found {
      chain_set.push(c)
    }
    i = i + 1
  }
  let mut a = 0
  while a < chain_set.length() {
    let mut b = a + 1
    while b < chain_set.length() {
      let ca = chain_set[a]
      let cb = chain_set[b]
      let (res_a, res_b) = self.chain_interface_residues(ca, cb, cutoff~)
      if res_a.length() >= 1 && res_b.length() >= 1 {
        let area = self.chain_interface_area(ca, cb)
        let pairs = self.chain_interaction_pairs(ca, cb, cutoff~)
        interfaces.push(ChainInterface::{
          chain_a: ca,
          chain_b: cb,
          interface_area: area,
          n_residues_a: res_a.length(),
          n_residues_b: res_b.length(),
          n_contacts: pairs.length(),
          residues_a: res_a,
          residues_b: res_b,
        })
      }
      b = b + 1
    }
    a = a + 1
  }
  interfaces
}

// ============================================================================
// SCADIO — Domain definition via label-propagation community detection
// ============================================================================

///|
/// One structural domain — a set of residue indices (global into
/// the AtomArray) with a dominant chain and size.
pub struct StructuralDomain {
  residues : Array[Int32]
  dominant_chain : Char
  n_residues : Int32
}

///|
/// Decompose this AtomArray into structural domains.
///
/// Algorithm: Label propagation on the residue-residue Cα contact graph.
///   1. Build CA-CA contact graph (undirected).
///   2. Initialize every residue with a unique label.
///   3. Iterate: each residue adopts the most common label among
///      its neighbors. Deterministic tie-break by smallest label id.
///   4. When no residue changes label → convergence.
///   5. Merge tiny communities (< min_domain_size residues) into
///      their largest adjacent neighbor.
pub fn AtomArray::find_domains(
  self : AtomArray,
  cutoff~ : Double = 8.0,
  min_domain_size~ : Int32 = 20,
  max_iter~ : Int32 = 50,
) -> Array[StructuralDomain] {
  // cutoff~ is directly accessible as cutoff
  // min_domain_size accessible directly
  // max_iter accessible directly
  let (cas, chains, _) = ca_per_residue(self)
  let n = cas.length()
  if n == 0 {
    return []
  }

  // --- Build undirected adjacency list ---
  let neighbors : Array[Array[Int32]] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let nb : Array[Int32] = Array(capacity=20)
    let mut j = 0
    while j < n {
      if i != j && cas[i].distance(cas[j]) <= cutoff {
        nb.push(j)
      }
      j = j + 1
    }
    neighbors.push(nb)
    i = i + 1
  }

  // --- Label propagation ---
  let mut labels : Array[Int32] = Array(capacity=n)
  let mut k = 0
  while k < n {
    labels.push(k)
    k = k + 1
  }

  let mut iteration = 0
  let mut changed = true
  while changed && iteration < max_iter {
    changed = false
    iteration = iteration + 1
    let new_labels : Array[Int32] = Array(capacity=n)
    let mut ci = 0
    while ci < n {
      // Count neighbor label frequencies
      let counts : Map[Int32, Int32] = Map([], capacity=8)
      for nb in neighbors[ci] {
        let lbl = labels[nb]
        let cur = match counts.get(lbl) {
          Some(v) => v
          None => 0
        }
        counts[lbl] = cur + 1
      }
      // Find label with max count; tie-break by smallest id
      let mut best_lbl = labels[ci]
      let mut best_cnt = -1
      let mut any = false
      for lbl in counts.keys() {
        let cnt = match counts.get(lbl) {
          Some(v) => v
          None => 0
        }
        if cnt > best_cnt {
          best_cnt = cnt
          best_lbl = lbl
          any = true
        } else if cnt == best_cnt && lbl < best_lbl {
          best_lbl = lbl
        }
      }
      if any && best_lbl != labels[ci] {
        changed = true
      }
      new_labels.push(best_lbl)
      ci = ci + 1
    }
    labels = new_labels
  }

  // --- Collect communities ---
  let domain_map : Map[Int32, Array[Int32]] = Map([], capacity=8)
  let mut ri = 0
  while ri < n {
    let lbl = labels[ri]
    let arr = match domain_map.get(lbl) {
      Some(a) => a
      None => Array(capacity=20)
    }
    arr.push(ri)
    domain_map[lbl] = arr
    ri = ri + 1
  }

  // Flatten to array
  let all_labels : Array[Int32] = Array(capacity=domain_map.length())
  let all_sizes : Array[Int32] = Array(capacity=domain_map.length())
  for lbl in domain_map.keys() {
    let sz = match domain_map.get(lbl) {
      Some(a) => a.length()
      None => 0
    }
    all_labels.push(lbl)
    all_sizes.push(sz)
  }

  // --- Filter tiny domains ---
  let big : Array[Bool] = Array(capacity=all_labels.length())
  let mut di = 0
  while di < all_labels.length() {
    big.push(all_sizes[di] >= min_domain_size)
    di = di + 1
  }

  // For each residue in a tiny domain, reassign to the most-frequent
  // big-domain label among its neighbors
  let final_labels : Array[Int32] = Array(capacity=n)
  let mut si = 0
  while si < n {
    let cur_lbl = labels[si]
    // Is cur_lbl a big domain?
    let mut cur_big = false
    let mut bi = 0
    while bi < all_labels.length() {
      if all_labels[bi] == cur_lbl {
        cur_big = big[bi]
        break
      }
      bi = bi + 1
    }
    if cur_big {
      final_labels.push(cur_lbl)
    } else {
      // Find most common big label among neighbors
      let counts2 : Map[Int32, Int32] = Map([], capacity=8)
      for nb in neighbors[si] {
        let nl = labels[nb]
        // Check if nl is a big domain
        let mut n_big = false
        let mut ni = 0
        while ni < all_labels.length() {
          if all_labels[ni] == nl {
            n_big = big[ni]
            break
          }
          ni = ni + 1
        }
        if n_big {
          let cur = match counts2.get(nl) {
            Some(v) => v
            None => 0
          }
          counts2[nl] = cur + 1
        }
      }
      let mut best_lbl = cur_lbl
      let mut best_cnt = 0
      let mut any = false
      for lbl in counts2.keys() {
        let cnt = match counts2.get(lbl) {
          Some(v) => v
          None => 0
        }
        if cnt > best_cnt {
          best_cnt = cnt
          best_lbl = lbl
          any = true
        } else if cnt == best_cnt && lbl < best_lbl {
          best_lbl = lbl
        }
      }
      if !any {
        final_labels.push(cur_lbl)
      } else {
        final_labels.push(best_lbl)
      }
    }
    si = si + 1
  }

  // --- Build final StructuralDomain array ---
  let final_map : Map[Int32, Array[Int32]] = Map([], capacity=8)
  let mut fi = 0
  while fi < n {
    let lbl = final_labels[fi]
    let arr = match final_map.get(lbl) {
      Some(a) => a
      None => Array(capacity=20)
    }
    arr.push(fi)
    final_map[lbl] = arr
    fi = fi + 1
  }

  let result : Array[StructuralDomain] = Array(capacity=final_map.length())
  for lbl in final_map.keys() {
    let residues = match final_map.get(lbl) {
      Some(a) => a
      None => Array(capacity=0)
    }
    // Find dominant chain
    let chain_counts : Map[Char, Int32] = Map([], capacity=4)
    for ri2 in residues {
      let c = chains[ri2]
      let cur = match chain_counts.get(c) {
        Some(v) => v
        None => 0
      }
      chain_counts[c] = cur + 1
    }
    let mut dom_chain : Char = if residues.length() > 0 {
      chains[residues[0]]
    } else {
      ' '
    }
    let mut dom_cnt = -1
    for c in chain_counts.keys() {
      let cnt = match chain_counts.get(c) {
        Some(v) => v
        None => 0
      }
      if cnt > dom_cnt {
        dom_cnt = cnt
        dom_chain = c
      }
    }
    result.push(StructuralDomain::{
      residues,
      dominant_chain: dom_chain,
      n_residues: residues.length(),
    })
  }
  result
}

///|
/// Quick text summary of domain decomposition.
pub fn summarize_domains(domains : Array[StructuralDomain]) -> String {
  let mut s = ""
  let mut i = 0
  while i < domains.length() {
    let d = domains[i]
    s = s +
      "Domain " +
      (i + 1).to_string() +
      " (chain " +
      d.dominant_chain.to_string() +
      "): " +
      d.n_residues.to_string() +
      " residues"
    if i + 1 < domains.length() {
      s = s + "\n"
    }
    i = i + 1
  }
  s
}

// ============================================================================
// Small helpers
// ============================================================================

///|
fn sum_doubles(arr : Array[Double]) -> Double {
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s
}
