///|
/// biotite_application_dockq.mbt — DockQ quality score.

///|
/// DockQ = 1/3 * (1/(1+(RMSD/8.5)^2) + 1/(1+(iRMSD/1.5)^2) + f_nat)
pub fn biotite_dockq_score(
  rmsd : Double,
  irmsd : Double,
  fnat : Double,
) -> Double {
  let term1 = 1.0 / (1.0 + @math.pow(rmsd / 8.5, 2.0))
  let term2 = 1.0 / (1.0 + @math.pow(irmsd / 1.5, 2.0))
  (term1 + term2 + fnat) / 3.0
}

///|
/// Classify a DockQ score.
pub fn biotite_dockq_class(score : Double) -> String {
  if score >= 0.80 {
    "High"
  } else if score >= 0.49 {
    "Medium"
  } else if score >= 0.23 {
    "Acceptable"
  } else {
    "Incorrect"
  }
}

///|
/// Command-line construction for the DockQ executable.
pub fn biotite_dockq_command(model : String, native : String) -> String {
  "DockQ " + model + " " + native
}

///|
/// Parse DockQ stdout: `DockQ 0.742 ...`.
pub fn biotite_dockq_parse_output(text : String) -> Double? {
  let lines = text.split("\n").to_array()
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.contains("DockQ") {
      let fields = trimmed.split(" ").to_array()
      for i in 1..<fields.length() {
        match parse_double(fields[i].trim().to_owned()) {
          Some(v) => return Some(v)
          None => ()
        }
      }
    }
  }
  None
}
