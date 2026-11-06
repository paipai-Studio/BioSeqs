///|
/// Hydrogen bond frequency analysis for multi-model structures.
///
/// Biotite's `hbond_frequency()` analyzes how often each hydrogen bond
/// occurs across multiple models (e.g., NMR ensembles, MD trajectories).
/// This module provides the same capability for BioSeqs.

///|
/// A hydrogen bond with its frequency across models.
pub struct HBondFrequency {
  donor_residue : Int32
  acceptor_residue : Int32
  donor_atom : String
  acceptor_atom : String
  frequency : Double // 0.0 to 1.0
  count : Int32 // number of models where HBond exists
  total_models : Int32
  avg_distance : Double // average D-A distance when present
  avg_angle : Double // average D-H-A angle when present
}

///|
/// Get the frequency value.
pub fn HBondFrequency::get_frequency(self : HBondFrequency) -> Double {
  self.frequency
}

///|
/// Get donor residue index.
pub fn HBondFrequency::get_donor_residue(self : HBondFrequency) -> Int32 {
  self.donor_residue
}

///|
/// Get acceptor residue index.
pub fn HBondFrequency::get_acceptor_residue(self : HBondFrequency) -> Int32 {
  self.acceptor_residue
}

///|
/// A single hydrogen bond observation in one model.
pub struct HBondObservation {
  donor_residue : Int32
  acceptor_residue : Int32
  donor_atom : String
  acceptor_atom : String
  distance : Double
  angle : Double
}

///|
/// Construct an HBondObservation.
pub fn HBondObservation::new(
  donor_residue~ : Int32,
  acceptor_residue~ : Int32,
  donor_atom~ : String,
  acceptor_atom~ : String,
  distance~ : Double,
  angle~ : Double,
) -> HBondObservation {
  HBondObservation::{
    donor_residue,
    acceptor_residue,
    donor_atom,
    acceptor_atom,
    distance,
    angle,
  }
}

///|
/// Compute hydrogen bond frequencies across multiple models.
///
/// Given a list of HBond observations per model, computes how often each
/// unique donor-acceptor pair forms a hydrogen bond.
///
/// Parameters:
/// - observations_per_model: array of arrays, one per model
/// - total_models: total number of models
pub fn hbond_frequency(
  observations_per_model : Array[Array[HBondObservation]],
  total_models~ : Int32,
) -> Array[HBondFrequency] {
  // Use a simple key: "donor_res:acceptor_res:donor_atom:acceptor_atom"
  let keys : Array[String] = []
  let counts : Array[Int32] = []
  let distances : Array[Double] = []
  let angles : Array[Double] = []
  for model_obs in observations_per_model {
    for obs in model_obs {
      let key = obs.donor_residue.to_string() +
        ":" +
        obs.acceptor_residue.to_string() +
        ":" +
        obs.donor_atom +
        ":" +
        obs.acceptor_atom
      let mut found = false
      for i in 0..<keys.length() {
        if keys[i] == key {
          counts[i] = counts[i] + 1
          distances[i] = distances[i] + obs.distance
          angles[i] = angles[i] + obs.angle
          found = true
          break
        }
      }
      if !found {
        keys.push(key)
        counts.push(1)
        distances.push(obs.distance)
        angles.push(obs.angle)
      }
    }
  }
  let result : Array[HBondFrequency] = []
  for i in 0..<keys.length() {
    let parts = keys[i].split(":").to_array()
    if parts.length() >= 4 {
      let d_res = parse_int_safe(parts[0].to_owned())
      let a_res = parse_int_safe(parts[1].to_owned())
      let d_atom = parts[2].to_owned()
      let a_atom = parts[3].to_owned()
      let cnt = counts[i]
      let freq = if total_models > 0 {
        cnt.to_double() / total_models.to_double()
      } else {
        0.0
      }
      let avg_d = if cnt > 0 { distances[i] / cnt.to_double() } else { 0.0 }
      let avg_a = if cnt > 0 { angles[i] / cnt.to_double() } else { 0.0 }
      result.push(HBondFrequency::{
        donor_residue: d_res,
        acceptor_residue: a_res,
        donor_atom: d_atom,
        acceptor_atom: a_atom,
        frequency: freq,
        count: cnt,
        total_models,
        avg_distance: avg_d,
        avg_angle: avg_a,
      })
    }
  }
  result
}

///|
/// Parse an integer from string, returning 0 on failure.
fn parse_int_safe(s : String) -> Int32 {
  let mut result = 0
  let mut sign = 1
  let mut i = 0
  if s.length() > 0 && s.unsafe_get(0).to_int() == 45 { // '-'
    sign = -1
    i = 1
  }
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c >= 48 && c <= 57 {
      result = result * 10 + (c - 48)
    } else {
      return 0
    }
    i = i + 1
  }
  sign * result
}

///|
/// Filter hydrogen bonds by minimum frequency.
pub fn filter_hbonds_by_frequency(
  hbonds : Array[HBondFrequency],
  min_frequency : Double,
) -> Array[HBondFrequency] {
  let result : Array[HBondFrequency] = []
  for hb in hbonds {
    if hb.frequency >= min_frequency {
      result.push(hb)
    }
  }
  result
}

///|
/// Sort hydrogen bonds by frequency (descending).
pub fn sort_hbonds_by_frequency(
  hbonds : Array[HBondFrequency],
) -> Array[HBondFrequency] {
  let result = Array(capacity=hbonds.length())
  for hb in hbonds {
    result.push(hb)
  }
  // Simple insertion sort
  for i in 1..<result.length() {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j].frequency < key.frequency {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
  }
  result
}
