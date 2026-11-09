/// omm_reporter_file.mbt — reporters that serialize to standard trajectory
/// file formats.
///
/// OpenMM equivalents: `DCDFile` / `PDBFile` writing from a `Simulation`.
///
/// - `OmmDCDWriter` accumulates simulation frames and serializes them to
///   X-PLOR DCD bytes through the shared `write_dcd` encoder, with a
///   Simulation-integrated `record_dcd` reporter loop.
/// - `omm_write_pdb_frames` renders multi-frame PDB text (MODEL/ENDMDL).
///
/// Units: A (coordinates), ps (timestep).

///|
/// Streaming DCD writer: collect frames one by one, serialize on demand.
pub struct OmmDCDWriter {
  n_atoms : Int32
  dt : Double
  frames : Array[Array[Array[Double]]]
}

///|
pub fn OmmDCDWriter::new(
  n_atoms~ : Int32,
  dt~ : Double = 0.001,
) -> OmmDCDWriter {
  { n_atoms, dt, frames: [], }
}

///|
/// Number of frames collected so far.
pub fn OmmDCDWriter::num_frames(self : OmmDCDWriter) -> Int32 {
  self.frames.length()
}

///|
/// Append one frame (coordinates are deep-copied).
pub fn OmmDCDWriter::add_frame(
  self : OmmDCDWriter,
  positions : Array[Array[Double]],
) -> Unit {
  let frame : Array[Array[Double]] = Array::make(positions.length(), [
    0.0, 0.0, 0.0,
  ])
  for i = 0; i < positions.length(); i = i + 1 {
    frame[i] = [positions[i][0], positions[i][1], positions[i][2]]
  }
  self.frames.push(frame)
}

///|
/// Serialize the collected frames to X-PLOR DCD bytes (file content).
pub fn OmmDCDWriter::to_dcd(self : OmmDCDWriter) -> String {
  write_dcd(Trajectory::{ frames: self.frames, })
}

///|
/// Simulation-integrated DCD reporter: run `n_steps` steps and record a
/// frame every `report_interval` steps (relative to the start of the run).
/// Returns the writer holding the recorded frames.
pub fn OmmSimulation::record_dcd(
  self : OmmSimulation,
  n_steps~ : Int32,
  report_interval~ : Int32,
  dt~ : Double = 0.001,
) -> OmmDCDWriter {
  let writer = OmmDCDWriter::new(n_atoms=self.state.positions.length(), dt~)
  let step0 = self.state.step
  for _s = 0; _s < n_steps; _s = _s + 1 {
    self.step(n_steps=1)
    if (self.state.step - step0) % report_interval == 0 {
      writer.add_frame(self.state.positions)
    }
  }
  writer
}

///|
/// Render frames as multi-model PDB text (MODEL/ENDMDL blocks). Atom
/// names/residue labels are synthetic (element-based); the output is a
/// valid multi-frame PDB for viewers that read MODEL records.
pub fn omm_write_pdb_frames(
  frames : Array[Array[Array[Double]]],
  atom_names~ : Array[String] = [],
  resname~ : String = "MOL",
) -> String {
  let sb : StringBuilder = StringBuilder()
  for f = 0; f < frames.length(); f = f + 1 {
    sb.write_string("MODEL     " + (f + 1).to_string() + "\n")
    let frame = frames[f]
    for i = 0; i < frame.length(); i = i + 1 {
      let name = if i < atom_names.length() { atom_names[i] } else { "X" }
      let serial = (i + 1).to_string()
      let mut padded = serial
      while padded.length() < 5 {
        padded = " " + padded
      }
      let mut nm = name
      while nm.length() < 4 {
        nm = " " + nm
      }
      let x = frame[i][0]
      let y = frame[i][1]
      let z = frame[i][2]
      let fx = format_coord(x)
      let fy = format_coord(y)
      let fz = format_coord(z)
      sb.write_string(
        "ATOM  " +
        padded +
        nm +
        " " +
        resname +
        "    1    " +
        fx +
        fy +
        fz +
        "\n",
      )
    }
    sb.write_string("TER\n")
    sb.write_string("ENDMDL\n")
  }
  sb.write_string("END\n")
  sb.to_string()
}

///|
/// Format a coordinate as PDB fixed-width %8.3f.
fn format_coord(v : Double) -> String {
  let scaled = (v * 1000.0).round().to_int()
  let neg = scaled < 0
  let abs_scaled = scaled.abs()
  let whole = abs_scaled / 1000
  let frac = abs_scaled % 1000
  let mut fr = frac.to_string()
  while fr.length() < 3 {
    fr = "0" + fr
  }
  let mut body = whole.to_string() + "." + fr
  if neg {
    body = "-" + body
  }
  while body.length() < 8 {
    body = " " + body
  }
  body
}
