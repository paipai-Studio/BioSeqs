///|
/// biotite_application_autodock.mbt — AutoDock Vina configuration,
/// command construction and output parsing.

///|
/// Build a Vina configuration file body.
pub fn biotite_vina_config(
  receptor : String,
  ligand : String,
  center_x : Double,
  center_y : Double,
  center_z : Double,
  size_x : Int32,
  size_y : Int32,
  size_z : Int32,
  out : String,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("receptor = " + receptor + "\n")
  let _ = sb.write("ligand = " + ligand + "\n")
  let _ = sb.write("center_x = " + center_x.to_string() + "\n")
  let _ = sb.write("center_y = " + center_y.to_string() + "\n")
  let _ = sb.write("center_z = " + center_z.to_string() + "\n")
  let _ = sb.write("size_x = " + size_x.to_string() + "\n")
  let _ = sb.write("size_y = " + size_y.to_string() + "\n")
  let _ = sb.write("size_z = " + size_z.to_string() + "\n")
  let _ = sb.write("out = " + out + "\n")
  sb.to_string()
}

///|
/// Build the Vina command line for a config file.
pub fn biotite_vina_command(config_file : String, cpu~ : Int32 = 0) -> String {
  if cpu > 0 {
    "vina --config " + config_file + " --cpu " + cpu.to_string()
  } else {
    "vina --config " + config_file
  }
}

///|
/// Parse the Vina results table.
/// Returns (mode, affinity, rmsd_lower_bound, rmsd_upper_bound).
pub fn biotite_vina_parse_output(
  text : String,
) -> Array[(Int32, Double, Double, Double)] {
  let lines = text.split("\n").to_array()
  let result : Array[(Int32, Double, Double, Double)] = []
  for line in lines {
    let trimmed = line.trim().to_owned()
    let is_header = trimmed.length() == 0 ||
      trimmed.has_prefix("-") ||
      trimmed.has_prefix("mode")
    if is_header {
      continue
    }
    let fields = trimmed.split(" ").to_array()
    let parts : Array[String] = []
    for f in fields {
      if f.length() > 0 {
        parts.push(f.to_owned())
      }
    }
    if parts.length() < 4 {
      continue
    }
    let mode = parse_int(parts[0])
    let affinity = match parse_double(parts[1]) {
      Some(v) => v
      None => 0.0
    }
    let lb = match parse_double(parts[2]) {
      Some(v) => v
      None => 0.0
    }
    let ub = match parse_double(parts[3]) {
      Some(v) => v
      None => 0.0
    }
    result.push((mode, affinity, lb, ub))
  }
  result
}

///|
/// Pose with the best (most negative) affinity, if any.
pub fn biotite_vina_best_pose(
  poses : Array[(Int32, Double, Double, Double)],
) -> (Int32, Double, Double, Double)? {
  if poses.length() == 0 {
    return None
  }
  let mut best = poses[0]
  for pose in poses {
    if pose.1 < best.1 {
      best = pose
    }
  }
  Some(best)
}
