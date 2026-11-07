///|
/// H5MD (HDF5 Molecular Data, https://nongnu.org/h5md/) object layer.
///
/// Real HDF5 binary I/O is out of scope; this module models the H5MD file
/// layout with self-contained lightweight structures (no dependency on the
/// rhdf5 types). The modeled layout is:
///
/// - `/h5md` group with a `version = [major, minor]` attribute
/// - `/particles/<group>/position|velocity|force` datasets, each carrying
///   `step` and `time` columns and a `value` of shape [frame][atom][3]
/// - `/particles/<group>/box` group with an `edges` dataset of shape
///   [frame][3][3] plus `dimension = 3`
/// - `/observables` (not stored here; validation ignores its absence)

///|
/// A time-dependent H5MD dataset: parallel `step`/`time` columns and a
/// `value` of shape [frame][atom][3].
pub(all) struct H5mdDataset {
  name : String
  step : Array[Int32]
  time : Array[Double]
  value : Array[Array[Array[Double]]]
}

///|
/// One `/particles/<name>` group: a position dataset (required), optional
/// velocity and force datasets, and optional box edges [frame][3][3].
pub(all) struct H5mdParticleGroup {
  name : String
  position : H5mdDataset
  velocity : H5mdDataset?
  force : H5mdDataset?
  box_edges : Array[Array[Array[Double]]]?
}

///|
/// In-memory H5MD file: version attribute plus particle groups and the
/// `/h5md` group metadata (author, creator).
pub(all) struct H5mdFile {
  h5md_version : (Int32, Int32)
  particles : Array[H5mdParticleGroup]
  author : String
  creator_name : String
  creator_version : String
}

///|
/// Error raised by H5MD object-layer operations.
pub suberror H5mdError {
  H5mdError(String)
}

///|
/// Create an empty H5MD file with version (1, 1).
pub fn h5md_file_new(
  author~ : String = "",
  creator_name~ : String = "BioSeqs",
  creator_version~ : String = "0.1.13",
) -> H5mdFile {
  let particles : Array[H5mdParticleGroup] = Array(capacity=0)
  H5mdFile::{
    h5md_version: (1, 1),
    particles,
    author,
    creator_name,
    creator_version,
  }
}

///|
/// Construct a dataset from name, step/time columns and [frame][atom][3]
/// values.
pub fn h5md_dataset_new(
  name : String,
  step : Array[Int32],
  time : Array[Double],
  value : Array[Array[Array[Double]]],
) -> H5mdDataset {
  H5mdDataset::{ name, step, time, value, }
}

///|
/// Add a particle group to the file and return its index. Raises on a
/// duplicate group name (H5MD names must be unique within /particles).
pub fn h5md_add_particles(
  f : H5mdFile,
  name : String,
  position : H5mdDataset,
  velocity~ : H5mdDataset? = None,
  force~ : H5mdDataset? = None,
  box_edges~ : Array[Array[Array[Double]]]? = None,
) -> Int32 raise H5mdError {
  for pg in f.particles {
    if pg.name == name {
      raise H5mdError("particle group '\{name}' already exists")
    }
  }
  let idx = f.particles.length()
  f.particles.push(H5mdParticleGroup::{
    name,
    position,
    velocity,
    force,
    box_edges,
  })
  idx
}

///|
/// Find a particle group by name; None when absent.
pub fn h5md_find_particles(f : H5mdFile, name : String) -> H5mdParticleGroup? {
  for pg in f.particles {
    if pg.name == name {
      return Some(pg)
    }
  }
  None
}

///|
/// Check one dataset: step/time lengths must be equal and match the number
/// of value frames. Problems are appended to `errors`.
fn h5md_check_dataset(
  prefix : String,
  ds : H5mdDataset,
  errors : Array[String],
) -> Unit {
  let n_frames = ds.value.length()
  if ds.step.length() != ds.time.length() {
    errors.push(
      "\{prefix}: step length \{ds.step.length()} != time length \{ds.time.length()}",
    )
  }
  if ds.step.length() != n_frames {
    errors.push(
      "\{prefix}: step length \{ds.step.length()} != frame count \{n_frames}",
    )
  }
}

///|
/// Whether two datasets agree on the per-frame atom count (when populated).
fn h5md_same_atom_shape(a : H5mdDataset, b : H5mdDataset) -> Bool {
  if a.value.length() == 0 || b.value.length() == 0 {
    return a.value.length() == b.value.length()
  }
  a.value[0].length() == b.value[0].length()
}

///|
/// Validate the file against the H5MD rules: version must be present, every
/// dataset's step/time lengths must be equal and match its frame count,
/// optional velocity/force must share the position atom shape, and box
/// edges must have as many frames as position. Returns the list of problems
/// (empty means valid). `/observables` is ignored when absent.
pub fn h5md_validate(f : H5mdFile) -> Array[String] {
  let errors : Array[String] = Array(capacity=0)
  let (major, minor) = f.h5md_version
  if major < 0 || minor < 0 {
    errors.push("h5md: version attribute missing or invalid")
  }
  for pg in f.particles {
    let prefix = "particles/\{pg.name}"
    h5md_check_dataset("\{prefix}/position", pg.position, errors)
    match pg.velocity {
      Some(v) => {
        h5md_check_dataset("\{prefix}/velocity", v, errors)
        if !h5md_same_atom_shape(pg.position, v) {
          errors.push("\{prefix}/velocity: atom count differs from position")
        }
      }
      None => ()
    }
    match pg.force {
      Some(fo) => {
        h5md_check_dataset("\{prefix}/force", fo, errors)
        if !h5md_same_atom_shape(pg.position, fo) {
          errors.push("\{prefix}/force: atom count differs from position")
        }
      }
      None => ()
    }
    match pg.box_edges {
      Some(box) =>
        if box.length() != pg.position.value.length() {
          errors.push(
            "\{prefix}/box: edge frame count \{box.length()} != position frame count \{pg.position.value.length()}",
          )
        }
      None => ()
    }
  }
  errors
}

///|
/// Convert a particle group's position dataset to a generic Trajectory.
pub fn h5md_to_trajectory(pg : H5mdParticleGroup) -> Trajectory {
  Trajectory::new(pg.position.value)
}

///|
/// Convert a Trajectory to an H5MD position dataset. Steps are 0..n-1 and
/// times are `i * time_step`.
pub fn trajectory_to_h5md_dataset(
  name : String,
  traj : Trajectory,
  time_step~ : Double = 0.0,
) -> H5mdDataset {
  let n = traj.num_frames()
  let step : Array[Int32] = Array(capacity=n)
  let time : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    step.push(i)
    time.push(i.to_double() * time_step)
  }
  H5mdDataset::{ name, step, time, value: traj.frames, }
}

///|
/// Whether an optional value is present (helper for JSON output).
fn[T] h5md_opt_is_some(o : T?) -> Bool {
  match o {
    Some(_) => true
    None => false
  }
}

///|
/// Escape a string for embedding in a JSON document.
fn h5md_json_escape(s : String) -> String {
  let buf = StringBuilder()
  for i = 0; i < s.length(); i = i + 1 {
    let c = s.unsafe_get(i).unsafe_to_char()
    if c == '"' {
      buf.write_string("\\\"")
    } else if c == '\\' {
      buf.write_string("\\\\")
    } else if c.to_int() < 32 {
      buf.write_char('?')
    } else {
      buf.write_char(c)
    }
  }
  buf.to_string()
}

///|
/// Serialize file metadata to a compact JSON summary: h5md version, author,
/// creator fields and one object per particle group (name, frame count,
/// atom count, presence flags for velocity/force/box).
pub fn h5md_to_json(f : H5mdFile) -> String {
  let buf = StringBuilder()
  buf.write_string(
    "{\"h5md_version\":[\{f.h5md_version.0},\{f.h5md_version.1}]",
  )
  buf.write_string(",\"author\":\"\{h5md_json_escape(f.author)}\"")
  buf.write_string(",\"creator_name\":\"\{h5md_json_escape(f.creator_name)}\"")
  buf.write_string(
    ",\"creator_version\":\"\{h5md_json_escape(f.creator_version)}\"",
  )
  buf.write_string(",\"particles\":[")
  for k = 0; k < f.particles.length(); k = k + 1 {
    if k > 0 {
      buf.write_char(',')
    }
    let pg = f.particles[k]
    let n_frames = pg.position.value.length()
    let n_atoms = if n_frames > 0 { pg.position.value[0].length() } else { 0 }
    buf.write_string(
      "{\"name\":\"\{h5md_json_escape(pg.name)}\",\"n_frames\":\{n_frames},\"n_atoms\":\{n_atoms}",
    )
    buf.write_string(
      ",\"has_velocity\":\{h5md_opt_is_some(pg.velocity)},\"has_force\":\{h5md_opt_is_some(pg.force)},\"has_box\":\{h5md_opt_is_some(pg.box_edges)}",
    )
    buf.write_char('}')
  }
  buf.write_string("]}")
  buf.to_string()
}
