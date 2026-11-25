///|
/// biotite_structure_traj.mbt — unified multi-model trajectory wrapper.

///|
pub(all) struct BiotiteTrajectoryFile {
  format : String
  models : Array[AtomArray]
}

///|
/// Build a trajectory from an explicit list of models.
pub fn biotite_trajectory_from_models(
  models : Array[AtomArray],
) -> BiotiteTrajectoryFile {
  { format: "models", models, }
}

///|
/// Parse a multi-MODEL PDB file into a trajectory.
/// A file without MODEL records becomes a single-frame trajectory.
pub fn biotite_trajectory_from_pdb(text : String) -> BiotiteTrajectoryFile {
  let lines = text.split("\n").to_array()
  let models : Array[AtomArray] = []
  let mut buf = StringBuilder()
  let mut in_model = false
  for line in lines {
    if line.has_prefix("MODEL") {
      in_model = true
      continue
    }
    if line.has_prefix("ENDMDL") {
      if in_model {
        models.push(AtomArray::from_pdb_text(buf.to_string()))
        buf = StringBuilder()
      }
      in_model = false
      continue
    }
    if in_model {
      let _ = buf.write_string(line.to_owned() + "\n")
    }
  }
  if models.length() == 0 {
    models.push(AtomArray::from_pdb_text(text))
  }
  { format: "pdb", models, }
}

///|
/// Number of frames.
pub fn biotite_trajectory_n_models(traj : BiotiteTrajectoryFile) -> Int32 {
  traj.models.length()
}

///|
/// One frame; out-of-range indices return an empty AtomArray.
pub fn biotite_trajectory_get_model(
  traj : BiotiteTrajectoryFile,
  index : Int32,
) -> AtomArray {
  if index >= 0 && index < traj.models.length() {
    traj.models[index]
  } else {
    AtomArray::new()
  }
}

///|
/// Concatenate trajectories frame by frame.
pub fn biotite_trajectory_concatenate(
  trajectories : Array[BiotiteTrajectoryFile],
) -> BiotiteTrajectoryFile {
  let models : Array[AtomArray] = []
  for traj in trajectories {
    for model in traj.models {
      models.push(model)
    }
  }
  { format: "pdb", models, }
}

///|
/// Serialize all frames as a multi-MODEL PDB string.
pub fn biotite_trajectory_to_pdb(traj : BiotiteTrajectoryFile) -> String {
  let sb = StringBuilder()
  for i in 0..<traj.models.length() {
    let _ = sb.write("MODEL     " + (i + 1).to_string() + "\n")
    let pdb = traj.models[i].to_pdb_string()
    let lines = pdb.split("\n").to_array()
    for line in lines {
      if line.to_owned() != "END" && line.length() > 0 {
        let _ = sb.write(line.to_owned() + "\n")
      }
    }
    let _ = sb.write("ENDMDL\n")
  }
  let _ = sb.write("END\n")
  sb.to_string()
}
