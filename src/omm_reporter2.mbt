/// omm_reporter2.mbt — PDB trajectory reporter + XTC trajectory reporter
/// Items 18-19: Additional trajectory output formats

///|
// === PDB Reporter ===

///|
pub struct OmmPdbReporter {
  interval : Int32
  mut step_count : Int32
  frames : Array[String] // PDB text per frame
}

///|
/// Create a PDB reporter that writes every `interval` steps.
pub fn OmmPdbReporter::new(interval~ : Int32 = 100) -> OmmPdbReporter {
  OmmPdbReporter::{ interval, step_count: 0, frames: Array(capacity=0), }
}

///|
/// Check if reporter should write at current step.
pub fn OmmPdbReporter::should_report(self : OmmPdbReporter) -> Bool {
  self.step_count % self.interval == 0
}

///|
/// Generate PDB frame text from coordinates.
pub fn omm_pdb_frame(
  coords : Array[(Double, Double, Double)],
  atom_names : Array[String],
  res_names : Array[String],
  res_ids : Array[Int32],
  step : Int32,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("MODEL     ")
  let _ = sb.write(step.to_string())
  let _ = sb.write("\n")
  for i = 0; i < coords.length(); i = i + 1 {
    let aname = if i < atom_names.length() { atom_names[i] } else { "X" }
    let rname = if i < res_names.length() { res_names[i] } else { "X" }
    let rid = if i < res_ids.length() { res_ids[i] } else { 1 }
    // ATOM line: fixed-column PDB format
    let _ = sb.write("ATOM  ")
    let _ = sb.write(omm_r2_pad_int(i + 1, 5))
    let _ = sb.write(" ")
    let _ = sb.write(omm_r2_pad_str(aname, 4))
    let _ = sb.write(" ")
    let _ = sb.write(omm_r2_pad_str(rname, 3))
    let _ = sb.write(" A")
    let _ = sb.write(omm_r2_pad_int(rid, 4))
    let _ = sb.write("    ")
    let _ = sb.write(omm_r2_fmt_coord(coords[i].0))
    let _ = sb.write(omm_r2_fmt_coord(coords[i].1))
    let _ = sb.write(omm_r2_fmt_coord(coords[i].2))
    let _ = sb.write("  1.00  0.00\n")
  }
  let _ = sb.write("ENDMDL\n")
  sb.to_string()
}

///|
/// Report a frame to the PDB reporter.
pub fn OmmPdbReporter::report(
  self : OmmPdbReporter,
  coords : Array[(Double, Double, Double)],
  atom_names : Array[String],
  res_names : Array[String],
  res_ids : Array[Int32],
) -> Unit {
  let frame = omm_pdb_frame(
    coords,
    atom_names,
    res_names,
    res_ids,
    self.step_count,
  )
  self.frames.push(frame)
  self.step_count = self.step_count + 1
}

///|
// === XTC Reporter ===

///|
pub struct OmmXtcReporter {
  interval : Int32
  mut step_count : Int32
  // XTC stores compressed frames; we store raw coords for simplicity
  frames : Array[(Int32, Array[(Double, Double, Double)])] // (step, coords)
  mut box : (Double, Double, Double)
}

///|
/// Create an XTC reporter.
pub fn OmmXtcReporter::new(interval~ : Int32 = 100) -> OmmXtcReporter {
  OmmXtcReporter::{
    interval,
    step_count: 0,
    frames: Array(capacity=0),
    box: (0.0, 0.0, 0.0),
  }
}

///|
/// Check if XTC reporter should write.
pub fn OmmXtcReporter::should_report(self : OmmXtcReporter) -> Bool {
  self.step_count % self.interval == 0
}

///|
/// Report a frame to the XTC reporter.
pub fn OmmXtcReporter::report(
  self : OmmXtcReporter,
  coords : Array[(Double, Double, Double)],
  box~ : (Double, Double, Double) = (0.0, 0.0, 0.0),
) -> Unit {
  // Store frame (in real XTC, this would compress with GROMACS algorithm)
  let copied : Array[(Double, Double, Double)] = Array(capacity=coords.length())
  for c in coords {
    copied.push(c)
  }
  self.frames.push((self.step_count, copied))
  self.box = box
  self.step_count = self.step_count + 1
}

///|
/// Get number of frames stored.
pub fn OmmXtcReporter::n_frames(self : OmmXtcReporter) -> Int32 {
  self.frames.length()
}

///|
/// Get a specific frame.
pub fn OmmXtcReporter::get_frame(
  self : OmmXtcReporter,
  index : Int32,
) -> (Int32, Array[(Double, Double, Double)])? {
  if index < 0 || index >= self.frames.length() {
    None
  } else {
    Some(self.frames[index])
  }
}

///|
// === Helper functions ===

///|
fn omm_r2_pad_int(n : Int32, width : Int32) -> String {
  let s = n.to_string()
  let len = s.length()
  if len >= width {
    s
  } else {
    let sb = StringBuilder()
    for _i = 0; _i < width - len; _i = _i + 1 {
      let _ = sb.write(" ")
    }
    let _ = sb.write(s)
    sb.to_string()
  }
}

///|
fn omm_r2_pad_str(s : String, width : Int32) -> String {
  let len = s.length()
  if len >= width {
    s
  } else {
    let sb = StringBuilder()
    for _i = 0; _i < width - len; _i = _i + 1 {
      let _ = sb.write(" ")
    }
    let _ = sb.write(s)
    sb.to_string()
  }
}

///|
fn omm_r2_fmt_coord(x : Double) -> String {
  // Format as 8.3f (simplified)
  let scaled = x * 1000.0
  let int_part = scaled.floor().to_int()
  let whole = int_part / 1000
  let frac = omm_r2_abs(int_part % 1000)
  let frac_s = if frac < 10 {
    "00" + frac.to_string()
  } else if frac < 100 {
    "0" + frac.to_string()
  } else {
    frac.to_string()
  }
  let sign = if x < 0.0 { "-" } else { " " }
  let s = sign + omm_r2_pad_int(whole, 7) + "." + frac_s
  omm_r2_pad_str(s, 8)
}

///|
fn omm_r2_abs(x : Int32) -> Int32 {
  if x < 0 {
    -x
  } else {
    x
  }
}

///|
pub fn omm_reporter2_dummy() -> Unit {
  ignore(OmmPdbReporter::{
    interval: 0,
    step_count: 0,
    frames: Array(capacity=0),
  })
  ignore(OmmXtcReporter::{
    interval: 0,
    step_count: 0,
    frames: Array(capacity=0),
    box: (0.0, 0.0, 0.0),
  })
}
