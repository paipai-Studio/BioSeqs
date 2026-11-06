/// mda_interactions.mbt — salt bridges, hydrophobic contacts, pi-stacking, coordination number, contact map trajectory
/// Items 1-5

///|
/// Salt bridge result
pub struct MdaSaltBridge {
  res1 : Int32
  res2 : Int32
  distance : Double
  occupancy : Double // fraction of frames
}

///|
/// Detect salt bridges between charged residues (simplified: uses CA positions).
pub fn mda_salt_bridges(
  coords : Array[(Double, Double, Double)],
  res_ids : Array[Int32],
  res_names : Array[String],
  cutoff~ : Double = 4.0,
) -> Array[MdaSaltBridge] {
  let n = coords.length()
  let bridges : Array[MdaSaltBridge] = Array(capacity=0)
  // Identify charged residues: ASP/GLU (negative), LYS/ARG/HIS (positive)
  for i = 0; i < n; i = i + 1 {
    let r1 = res_names[i]
    let is_pos1 = r1 == "LYS" || r1 == "ARG" || r1 == "HIS"
    let is_neg1 = r1 == "ASP" || r1 == "GLU"
    if !is_pos1 && !is_neg1 {
      continue
    }
    for j = i + 1; j < n; j = j + 1 {
      let r2 = res_names[j]
      let is_pos2 = r2 == "LYS" || r2 == "ARG" || r2 == "HIS"
      let is_neg2 = r2 == "ASP" || r2 == "GLU"
      // Salt bridge: one positive, one negative
      if (is_pos1 && is_neg2) || (is_neg1 && is_pos2) {
        let d = mda_int_dist(coords[i], coords[j])
        if d <= cutoff {
          bridges.push(MdaSaltBridge::{
            res1: res_ids[i],
            res2: res_ids[j],
            distance: d,
            occupancy: 1.0,
          })
        }
      }
    }
  }
  bridges
}

///|
/// Detect salt bridges over a trajectory (time series).
pub fn mda_salt_bridges_traj(
  frames : Array[Array[(Double, Double, Double)]],
  res_ids : Array[Int32],
  res_names : Array[String],
  cutoff~ : Double = 4.0,
) -> Array[MdaSaltBridge] {
  let n_frames = frames.length()
  if n_frames == 0 {
    return Array(capacity=0)
  }
  // Accumulate per-frame results
  let bridge_map : Map[String, (MdaSaltBridge, Double)] = Map([], capacity=32)
  for f = 0; f < n_frames; f = f + 1 {
    let frame_bridges = mda_salt_bridges(frames[f], res_ids, res_names, cutoff~)
    for i = 0; i < frame_bridges.length(); i = i + 1 {
      let b = frame_bridges[i]
      let key = b.res1.to_string() + "-" + b.res2.to_string()
      let cur = bridge_map.get(key)
      match cur {
        Some(v) => bridge_map[key] = (v.0, v.1 + 1.0)
        None => bridge_map[key] = (b, 1.0)
      }
    }
  }
  let result : Array[MdaSaltBridge] = Array(capacity=0)
  bridge_map.each(fn(_k, v) {
    let occupancy = v.1 / n_frames.to_double()
    result.push(MdaSaltBridge::{
      res1: v.0.res1,
      res2: v.0.res2,
      distance: v.0.distance,
      occupancy,
    })
  })
  result
}

///|
/// Hydrophobic contact result
pub struct MdaHydrophobicContact {
  res1 : Int32
  res2 : Int32
  distance : Double
}

///|
/// Detect hydrophobic contacts between nonpolar residues.
pub fn mda_hydrophobic_contacts(
  coords : Array[(Double, Double, Double)],
  res_ids : Array[Int32],
  res_names : Array[String],
  cutoff~ : Double = 5.0,
) -> Array[MdaHydrophobicContact] {
  let hydrophobic = ["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "PRO"]
  let n = coords.length()
  let contacts : Array[MdaHydrophobicContact] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if !mda_int_in_array(res_names[i], hydrophobic) {
      continue
    }
    for j = i + 1; j < n; j = j + 1 {
      if !mda_int_in_array(res_names[j], hydrophobic) {
        continue
      }
      let d = mda_int_dist(coords[i], coords[j])
      if d <= cutoff {
        contacts.push(MdaHydrophobicContact::{
          res1: res_ids[i],
          res2: res_ids[j],
          distance: d,
        })
      }
    }
  }
  contacts
}

///|
/// Pi-stacking result
pub struct MdaPiStacking {
  ring1_res : Int32
  ring2_res : Int32
  distance : Double // centroid-centroid distance
  angle : Double // angle between ring normals (degrees)
  stacking_type : String // "face-to-face" or "edge-to-face"
}

///|
/// Detect pi-pi stacking between aromatic residues (PHE, TYR, TRP, HIS).
pub fn mda_pi_stacking(
  coords : Array[(Double, Double, Double)],
  res_ids : Array[Int32],
  res_names : Array[String],
  dist_cutoff~ : Double = 6.0,
  angle_cutoff~ : Double = 30.0,
) -> Array[MdaPiStacking] {
  let aromatic = ["PHE", "TYR", "TRP", "HIS"]
  let n = coords.length()
  let stackings : Array[MdaPiStacking] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if !mda_int_in_array(res_names[i], aromatic) {
      continue
    }
    for j = i + 1; j < n; j = j + 1 {
      if !mda_int_in_array(res_names[j], aromatic) {
        continue
      }
      let d = mda_int_dist(coords[i], coords[j])
      if d <= dist_cutoff {
        // Simplified: use random normal angles since we don't have ring atoms
        let angle = 0.0 // would compute from ring plane normals
        let stacking_type = if angle < angle_cutoff {
          "face-to-face"
        } else {
          "edge-to-face"
        }
        stackings.push(MdaPiStacking::{
          ring1_res: res_ids[i],
          ring2_res: res_ids[j],
          distance: d,
          angle,
          stacking_type,
        })
      }
    }
  }
  stackings
}

///|
/// Coordination number result
pub struct MdaCoordination {
  atom_idx : Int32
  n_neighbors : Int32
  neighbor_indices : Array[Int32]
}

///|
/// Compute coordination number for each atom.
pub fn mda_coordination_number(
  coords : Array[(Double, Double, Double)],
  cutoff~ : Double = 3.5,
) -> Array[MdaCoordination] {
  let n = coords.length()
  let result : Array[MdaCoordination] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let neighbors : Array[Int32] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      let d = mda_int_dist(coords[i], coords[j])
      if d <= cutoff {
        neighbors.push(j)
      }
    }
    result.push(MdaCoordination::{
      atom_idx: i,
      n_neighbors: neighbors.length(),
      neighbor_indices: neighbors,
    })
  }
  result
}

///|
/// Contact map trajectory: time-resolved contact frequency matrix.
pub fn mda_contact_map_trajectory(
  frames : Array[Array[(Double, Double, Double)]],
  cutoff~ : Double = 8.0,
) -> Array[Array[Double]] {
  let n_frames = frames.length()
  if n_frames == 0 {
    return Array(capacity=0)
  }
  let n_atoms = frames[0].length()
  let freq : Array[Array[Double]] = Array(capacity=n_atoms)
  for i = 0; i < n_atoms; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_atoms)
    for j = 0; j < n_atoms; j = j + 1 {
      row.push(0.0)
    }
    freq.push(row)
  }
  for f = 0; f < n_frames; f = f + 1 {
    let frame = frames[f]
    for i = 0; i < n_atoms; i = i + 1 {
      for j = i + 1; j < n_atoms; j = j + 1 {
        let d = mda_int_dist(frame[i], frame[j])
        if d <= cutoff {
          freq[i][j] = freq[i][j] + 1.0
          freq[j][i] = freq[j][i] + 1.0
        }
      }
    }
  }
  // Normalize by number of frames
  let nf = n_frames.to_double()
  for i = 0; i < n_atoms; i = i + 1 {
    for j = 0; j < n_atoms; j = j + 1 {
      freq[i][j] = freq[i][j] / nf
    }
  }
  freq
}

///|
/// Contact lifetime: average duration of persistent contacts.
pub fn mda_contact_lifetimes(
  frames : Array[Array[(Double, Double, Double)]],
  cutoff~ : Double = 8.0,
) -> Array[(Int32, Int32, Double)] {
  let n_frames = frames.length()
  if n_frames == 0 {
    return Array(capacity=0)
  }
  let n_atoms = frames[0].length()
  let lifetimes : Array[(Int32, Int32, Double)] = Array(capacity=0)
  for i = 0; i < n_atoms; i = i + 1 {
    for j = i + 1; j < n_atoms; j = j + 1 {
      let mut total_duration = 0
      let mut current_duration = 0
      for f = 0; f < n_frames; f = f + 1 {
        let d = mda_int_dist(frames[f][i], frames[f][j])
        if d <= cutoff {
          current_duration = current_duration + 1
        } else if current_duration > 0 {
          total_duration = total_duration + current_duration
          current_duration = 0
        }
      }
      if current_duration > 0 {
        total_duration = total_duration + current_duration
      }
      if total_duration > 0 {
        let avg_lifetime = total_duration.to_double() / n_frames.to_double()
        lifetimes.push((i, j, avg_lifetime))
      }
    }
  }
  lifetimes
}

///|
fn mda_int_dist(
  a : (Double, Double, Double),
  b : (Double, Double, Double),
) -> Double {
  let dx = a.0 - b.0
  let dy = a.1 - b.1
  let dz = a.2 - b.2
  @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
}

///|
fn mda_int_in_array(s : String, arr : Array[String]) -> Bool {
  for i = 0; i < arr.length(); i = i + 1 {
    if arr[i] == s {
      return true
    }
  }
  false
}

///|
pub fn mda_interactions_dummy() -> Unit {
  ignore(MdaSaltBridge::{ res1: 0, res2: 0, distance: 0.0, occupancy: 0.0, })
  ignore(MdaHydrophobicContact::{ res1: 0, res2: 0, distance: 0.0, })
  ignore(MdaPiStacking::{
    ring1_res: 0,
    ring2_res: 0,
    distance: 0.0,
    angle: 0.0,
    stacking_type: "",
  })
  ignore(MdaCoordination::{
    atom_idx: 0,
    n_neighbors: 0,
    neighbor_indices: Array(capacity=0),
  })
}
