///|
/// Model quality assessment: GDT_TS, GDT_HA, lDDT, DockQ, DRMSD.
///
/// These are standard scores used in CASP and structural biology
/// to compare a model structure against a reference (native) structure.

///|
/// Result of GDT (Global Distance Test) computation.
pub struct GDTResult {
  gdt_ts : Double
  gdt_ha : Double
  p05 : Double
  p1 : Double
  p2 : Double
  p4 : Double
  p8 : Double
}

///|
/// Result of lDDT (local Distance Difference Test) computation.
pub struct LDDTResult {
  lddt : Double
  per_residue : Array[Double]
  cutoff : Double
  thresholds : Array[Double]
}

///|
/// Result of DockQ computation.
pub struct DockQResult {
  dockq : Double
  fnat : Double
  fnonnat : Double
  irms : Double
  lrms : Double
}

// ─── helpers ───────────────────────────────────────────────

fn mq_sqrt(x : Double) -> Double {
  if x < 0.0 { 0.0 } else { x.sqrt() }
}

///|
fn to_qcp_coords(coords : Array[Vector3]) -> Array[QCPAtomCoordinate] {
  let n = coords.length()
  let result = Array(capacity=n)
  let mut i = 0
  while i < n {
    let v = coords[i]
    result.push(QCPAtomCoordinate::new(v.x, v.y, v.z))
    i = i + 1
  }
  result
}

///|
fn transform_point(
  rot : Array[Array[Double]],
  tran : Array[Double],
  v : Vector3,
) -> Vector3 {
  let rx = rot[0][0] * v.x + rot[0][1] * v.y + rot[0][2] * v.z + tran[0]
  let ry = rot[1][0] * v.x + rot[1][1] * v.y + rot[1][2] * v.z + tran[1]
  let rz = rot[2][0] * v.x + rot[2][1] * v.y + rot[2][2] * v.z + tran[2]
  Vector3::new(rx, ry, rz)
}

///|
/// Superpose model onto ref using QCP, return transformed model coords.
fn superpose_subset(
  ref_coords : Array[Vector3],
  model_coords : Array[Vector3],
  indices : Array[Int],
) -> Array[Vector3] {
  let ref_sub = Array(capacity=indices.length())
  let mod_sub = Array(capacity=indices.length())
  for idx in indices {
    ref_sub.push(ref_coords[idx])
    mod_sub.push(model_coords[idx])
  }
  // QCP: fixed=ref, moving=model → rotation transforms model onto ref
  let sup = QCPSuperimposer::new().set(
    to_qcp_coords(ref_sub),
    to_qcp_coords(mod_sub),
  )
  let rot = sup.get_rot()
  let tran = sup.get_tran()
  let n = model_coords.length()
  let result = Array(capacity=n)
  let mut i = 0
  while i < n {
    result.push(transform_point(rot, tran, model_coords[i]))
    i = i + 1
  }
  result
}

///|
fn gdt_maximize(
  ref_coords : Array[Vector3],
  model_coords : Array[Vector3],
  threshold : Double,
  seed_indices : Array[Int],
) -> Double {
  let n = ref_coords.length()
  if n == 0 { return 0.0 }
  let mut current = seed_indices
  let mut iter = 0
  while iter < 50 {
    let transformed = superpose_subset(ref_coords, model_coords, current)
    let new_set = Array(capacity=n)
    let mut i = 0
    while i < n {
      let d = ref_coords[i].distance(transformed[i])
      if d <= threshold { new_set.push(i) }
      i = i + 1
    }
    if new_set.length() > current.length() {
      current = new_set
    } else {
      break
    }
    iter = iter + 1
  }
  let transformed = superpose_subset(ref_coords, model_coords, current)
  let mut count = 0
  let mut i = 0
  while i < n {
    let d = ref_coords[i].distance(transformed[i])
    if d <= threshold { count = count + 1 }
    i = i + 1
  }
  count.to_double() / n.to_double() * 100.0
}

///|
fn gdt_for_threshold(
  ref_coords : Array[Vector3],
  model_coords : Array[Vector3],
  threshold : Double,
) -> Double {
  let n = ref_coords.length()
  if n < 3 {
    let all_indices = Array(capacity=n)
    let mut i = 0
    while i < n { all_indices.push(i); i = i + 1 }
    return gdt_maximize(ref_coords, model_coords, threshold, all_indices)
  }
  let mut best = 0.0
  let step = if n > 30 { n / 10 } else { 1 }
  let mut start = 0
  while start < n - 2 {
    let seed = [start, start + 1, start + 2]
    let score = gdt_maximize(ref_coords, model_coords, threshold, seed)
    if score > best { best = score }
    start = start + step
  }
  best
}

///|
/// Compute GDT_TS and GDT_HA.
/// GDT_TS = (P1 + P2 + P4 + P8) / 4
/// GDT_HA = (P0.5 + P1 + P2 + P4) / 4
pub fn calc_gdt_ts(
  ref_coords : Array[Vector3],
  model_coords : Array[Vector3],
) -> GDTResult {
  let n = ref_coords.length()
  if n == 0 || n != model_coords.length() {
    return GDTResult::{ gdt_ts: 0.0, gdt_ha: 0.0, p05: 0.0, p1: 0.0, p2: 0.0, p4: 0.0, p8: 0.0 }
  }
  let p05 = gdt_for_threshold(ref_coords, model_coords, 0.5)
  let p1 = gdt_for_threshold(ref_coords, model_coords, 1.0)
  let p2 = gdt_for_threshold(ref_coords, model_coords, 2.0)
  let p4 = gdt_for_threshold(ref_coords, model_coords, 4.0)
  let p8 = gdt_for_threshold(ref_coords, model_coords, 8.0)
  let gdt_ts = (p1 + p2 + p4 + p8) / 4.0
  let gdt_ha = (p05 + p1 + p2 + p4) / 4.0
  GDTResult::{ gdt_ts, gdt_ha, p05, p1, p2, p4, p8 }
}

// ─── lDDT ──────────────────────────────────────────────────

///|
/// Compute lDDT (local Distance Difference Test).
/// For each residue, compare distances to all others within cutoff
/// in the reference. Score = average step function over thresholds.
pub fn calc_lddt(
  ref_coords : Array[Vector3],
  model_coords : Array[Vector3],
  cutoff? : Double = 15.0,
) -> LDDTResult {
  let n = ref_coords.length()
  let thresholds = [0.5, 1.0, 2.0, 4.0]
  if n == 0 || n != model_coords.length() {
    return LDDTResult::{ lddt: 0.0, per_residue: [], cutoff, thresholds }
  }
  let per_residue = Array::make(n, 0.0)
  let mut total_score = 0.0
  let mut count_residues = 0
  let mut i = 0
  while i < n {
    let mut score_sum = 0.0
    let mut pair_count = 0
    let mut j = 0
    while j < n {
      if i == j { j = j + 1; continue }
      let d_ref = ref_coords[i].distance(ref_coords[j])
      if d_ref > cutoff { j = j + 1; continue }
      let d_model = model_coords[i].distance(model_coords[j])
      let diff = (d_model - d_ref).abs()
      let mut thresh_sum = 0.0
      for t in thresholds {
        if diff < t { thresh_sum = thresh_sum + 1.0 }
      }
      score_sum = score_sum + thresh_sum / thresholds.length().to_double()
      pair_count = pair_count + 1
      j = j + 1
    }
    if pair_count > 0 {
      per_residue[i] = score_sum / pair_count.to_double()
      total_score = total_score + per_residue[i]
      count_residues = count_residues + 1
    } else {
      per_residue[i] = 1.0
      total_score = total_score + 1.0
      count_residues = count_residues + 1
    }
    i = i + 1
  }
  let lddt = if count_residues > 0 {
    total_score / count_residues.to_double()
  } else { 0.0 }
  LDDTResult::{ lddt, per_residue, cutoff, thresholds }
}

// ─── DRMSD ─────────────────────────────────────────────────

///|
/// Compute distance RMSD (DRMSD).
/// DRMSD = sqrt( (1/N) * sum_{i<j} (d_model(i,j) - d_ref(i,j))^2 )
pub fn calc_drmsd(
  ref_coords : Array[Vector3],
  model_coords : Array[Vector3],
) -> Double {
  let n = ref_coords.length()
  if n < 2 || n != model_coords.length() { return 0.0 }
  let mut sum_sq = 0.0
  let mut count = 0
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let d_ref = ref_coords[i].distance(ref_coords[j])
      let d_model = model_coords[i].distance(model_coords[j])
      let diff = d_model - d_ref
      sum_sq = sum_sq + diff * diff
      count = count + 1
      j = j + 1
    }
    i = i + 1
  }
  if count == 0 { return 0.0 }
  mq_sqrt(sum_sq / count.to_double())
}

// ─── DockQ ─────────────────────────────────────────────────

///|
/// Compute DockQ score for protein-protein docking.
/// DockQ = 2 / (1 + exp(-z)), z = 0.05*(40*Fnat - iRMS - 0.5*LRMS)
pub fn calc_dockq(
  receptor_ref : Array[Vector3],
  ligand_ref : Array[Vector3],
  receptor_model : Array[Vector3],
  ligand_model : Array[Vector3],
  contact_cutoff? : Double = 5.0,
) -> DockQResult {
  let empty = DockQResult::{ dockq: 0.0, fnat: 0.0, fnonnat: 1.0, irms: 0.0, lrms: 0.0 }
  if receptor_ref.length() == 0 || ligand_ref.length() == 0 { return empty }
  if receptor_ref.length() != receptor_model.length() { return empty }
  if ligand_ref.length() != ligand_model.length() { return empty }
  // Native contacts
  let native_contacts = Array(capacity=0)
  let mut i = 0
  while i < receptor_ref.length() {
    let mut j = 0
    while j < ligand_ref.length() {
      let d = receptor_ref[i].distance(ligand_ref[j])
      if d <= contact_cutoff { native_contacts.push((i, j)) }
      j = j + 1
    }
    i = i + 1
  }
  // Model contacts
  let model_contacts = Array(capacity=0)
  i = 0
  while i < receptor_model.length() {
    let mut j = 0
    while j < ligand_model.length() {
      let d = receptor_model[i].distance(ligand_model[j])
      if d <= contact_cutoff { model_contacts.push((i, j)) }
      j = j + 1
    }
    i = i + 1
  }
  // Fnat
  let mut native_preserved = 0
  for (ri, lj) in native_contacts {
    let d_model = receptor_model[ri].distance(ligand_model[lj])
    if d_model <= contact_cutoff { native_preserved = native_preserved + 1 }
  }
  let fnat = if native_contacts.length() > 0 {
    native_preserved.to_double() / native_contacts.length().to_double()
  } else { 0.0 }
  // Fnonnat
  let mut non_native = 0
  for (ri, lj) in model_contacts {
    let mut is_native = false
    for (nri, nlj) in native_contacts {
      if nri == ri && nlj == lj { is_native = true; break }
    }
    if not(is_native) { non_native = non_native + 1 }
  }
  let fnonnat = if model_contacts.length() > 0 {
    non_native.to_double() / model_contacts.length().to_double()
  } else { 1.0 }
  // iRMS: interface RMSD
  let rec_iface = Map::new()
  let lig_iface = Map::new()
  for (ri, lj) in native_contacts { rec_iface[ri] = true; lig_iface[lj] = true }
  for (ri, lj) in model_contacts { rec_iface[ri] = true; lig_iface[lj] = true }
  let ref_interface = Array(capacity=0)
  let mod_interface = Array(capacity=0)
  for ri in rec_iface.keys() {
    ref_interface.push(receptor_ref[ri])
    mod_interface.push(receptor_model[ri])
  }
  for lj in lig_iface.keys() {
    ref_interface.push(ligand_ref[lj])
    mod_interface.push(ligand_model[lj])
  }
  let irms = if ref_interface.length() >= 3 {
    let sup = QCPSuperimposer::new().set(
      to_qcp_coords(ref_interface),
      to_qcp_coords(mod_interface),
    )
    sup.get_rmsd()
  } else { 0.0 }
  // LRMS: ligand RMSD after superposing receptor
  let rec_sup = QCPSuperimposer::new().set(
    to_qcp_coords(receptor_ref),
    to_qcp_coords(receptor_model),
  )
  let rot = rec_sup.get_rot()
  let tran = rec_sup.get_tran()
  let mut lrms_sum = 0.0
  let nl = ligand_model.length()
  let mut k = 0
  while k < nl {
    let transformed = transform_point(rot, tran, ligand_model[k])
    let dx = ligand_ref[k].x - transformed.x
    let dy = ligand_ref[k].y - transformed.y
    let dz = ligand_ref[k].z - transformed.z
    lrms_sum = lrms_sum + dx * dx + dy * dy + dz * dz
    k = k + 1
  }
  let lrms = mq_sqrt(lrms_sum / nl.to_double())
  // DockQ
  let z = 0.05 * (40.0 * fnat - irms - 0.5 * lrms)
  let dockq = 2.0 / (1.0 + @math.exp(-z))
  DockQResult::{ dockq, fnat, fnonnat, irms, lrms }
}
