/// Water bridge analysis and water residence time.
///
/// MDAnalysis equivalent: `HydrogenBondAnalysis` (water bridges),
/// `WaterResidenceTime`

///|
/// A water bridge between two residues
pub struct MdaWaterBridge {
  residue1 : Int32
  residue2 : Int32
  water_idx : Int32
  distance1 : Double
  distance2 : Double
}

///|
pub fn MdaWaterBridge::residue1(self : MdaWaterBridge) -> Int32 {
  self.residue1
}

///|
pub fn MdaWaterBridge::residue2(self : MdaWaterBridge) -> Int32 {
  self.residue2
}

///|
pub fn MdaWaterBridge::water_idx(self : MdaWaterBridge) -> Int32 {
  self.water_idx
}

///|
/// Find water bridges in a single frame.
///
/// protein_coords: coordinates of protein atoms (with residue assignment)
/// water_coords: coordinates of water oxygen atoms
/// residue_of_atom: residue index for each protein atom
/// cutoff: maximum distance for a hydrogen bond (default 3.0 Å)
pub fn mda_find_water_bridges(
  protein_coords : Array[Array[Double]],
  water_coords : Array[Array[Double]],
  residue_of_atom : Array[Int32],
  cutoff~ : Double = 3.0,
) -> Array[MdaWaterBridge] {
  let bridges : Array[MdaWaterBridge] = Array(capacity=0)
  let n_protein = protein_coords.length()
  let n_water = water_coords.length()
  for w = 0; w < n_water; w = w + 1 {
    // Find protein atoms near this water
    let nearby : Array[(Int32, Double)] = Array(capacity=0)
    for p = 0; p < n_protein; p = p + 1 {
      let dx = protein_coords[p][0] - water_coords[w][0]
      let dy = protein_coords[p][1] - water_coords[w][1]
      let dz = protein_coords[p][2] - water_coords[w][2]
      let d = (dx * dx + dy * dy + dz * dz).sqrt()
      if d <= cutoff {
        nearby.push((p, d))
      }
    }
    // Check for bridges between different residues
    for i = 0; i < nearby.length(); i = i + 1 {
      for j = i + 1; j < nearby.length(); j = j + 1 {
        let r1 = residue_of_atom[nearby[i].0]
        let r2 = residue_of_atom[nearby[j].0]
        if r1 != r2 {
          bridges.push({
            residue1: r1,
            residue2: r2,
            water_idx: w,
            distance1: nearby[i].1,
            distance2: nearby[j].1,
          })
        }
      }
    }
  }
  bridges
}

///|
/// Count unique residue pairs bridged by water.
pub fn mda_water_bridge_pairs(bridges : Array[MdaWaterBridge]) -> Int32 {
  let pairs : Array[(Int32, Int32)] = Array(capacity=0)
  for b in bridges {
    let (r1, r2) = if b.residue1 < b.residue2 {
      (b.residue1, b.residue2)
    } else {
      (b.residue2, b.residue1)
    }
    let mut found = false
    for p in pairs {
      if p.0 == r1 && p.1 == r2 {
        found = true
        break
      }
    }
    if !found {
      pairs.push((r1, r2))
    }
  }
  pairs.length()
}

///|
/// Water residence time near a residue.
///
/// traj_water: water oxygen coordinates per frame [frame][water][xyz]
/// residue_coords: residue center coordinates per frame [frame][xyz]
/// cutoff: distance cutoff for "near"
/// Returns average residence time (in frames).
pub fn mda_water_residence_time(
  traj_water : Array[Array[Array[Double]]],
  residue_coords : Array[Array[Double]],
  cutoff~ : Double = 4.0,
) -> Double {
  let n_frames = traj_water.length()
  if n_frames == 0 {
    return 0.0
  }
  let n_water = traj_water[0].length()
  // For each water, compute continuous residence durations
  let mut total_duration = 0.0
  let mut total_events = 0.0
  for w = 0; w < n_water; w = w + 1 {
    let mut in_residence = false
    let mut duration = 0
    for f = 0; f < n_frames; f = f + 1 {
      let dx = traj_water[f][w][0] - residue_coords[f][0]
      let dy = traj_water[f][w][1] - residue_coords[f][1]
      let dz = traj_water[f][w][2] - residue_coords[f][2]
      let d = (dx * dx + dy * dy + dz * dz).sqrt()
      if d <= cutoff {
        if !in_residence {
          in_residence = true
          duration = 1
        } else {
          duration = duration + 1
        }
      } else if in_residence {
        total_duration = total_duration + duration.to_double()
        total_events = total_events + 1.0
        in_residence = false
      }
    }
    // Handle residence at end of trajectory
    if in_residence {
      total_duration = total_duration + duration.to_double()
      total_events = total_events + 1.0
    }
  }
  if total_events > 0.0 {
    total_duration / total_events
  } else {
    0.0
  }
}

///|
/// Water occupancy: fraction of frames each water is near a residue.
pub fn mda_water_occupancy(
  traj_water : Array[Array[Array[Double]]],
  residue_coords : Array[Array[Double]],
  cutoff~ : Double = 4.0,
) -> Array[Double] {
  let n_frames = traj_water.length()
  let n_water = if n_frames > 0 { traj_water[0].length() } else { 0 }
  let occupancy : Array[Double] = Array::make(n_water, 0.0)
  if n_frames == 0 {
    return occupancy
  }
  for w = 0; w < n_water; w = w + 1 {
    let mut count = 0
    for f = 0; f < n_frames; f = f + 1 {
      let dx = traj_water[f][w][0] - residue_coords[f][0]
      let dy = traj_water[f][w][1] - residue_coords[f][1]
      let dz = traj_water[f][w][2] - residue_coords[f][2]
      let d = (dx * dx + dy * dy + dz * dz).sqrt()
      if d <= cutoff {
        count = count + 1
      }
    }
    occupancy[w] = count.to_double() / n_frames.to_double()
  }
  occupancy
}
