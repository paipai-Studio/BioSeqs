///|
/// ProDy-style PDB trajectory I/O.
///
/// Treats a multi-model PDB file as a trajectory: each MODEL/ENDMDL
/// block is one frame (conformation snapshot).  Provides frame-level
/// operations: RMSD matrix, RMSF, mean structure, superposition,
/// slicing, and displacement vectors.
///
/// Reference: ProDy `parsePDB` / `writePDB` / `Trajectory` classes.

///|
/// A multi-model PDB treated as a trajectory.
///
/// `frames[frame_idx][atom_idx]` = `Vector3` coordinates.
/// Atom metadata (names, resnames, chainids, …) is taken from the
/// first model and assumed constant across frames.
pub struct PDBTrajectory {
  n_frames : Int
  n_atoms : Int
  atom_names : Array[String]
  resnames : Array[String]
  chainids : Array[Char]
  resseqs : Array[Int]
  elements : Array[String]
  bfactors : Array[Double]
  occupancies : Array[Double]
  hetflags : Array[Bool]
  frames : Array[Array[Vector3]]
}

///|
/// Build a PDBTrajectory from a parsed Structure.
///
/// All models must have the same number of atoms in the same order.
/// Returns `None` if the structure has no models or inconsistent atom counts.
pub fn PDBTrajectory::from_structure(structure : Structure) -> PDBTrajectory? {
  let models = structure.get_models()
  let n_models = models.length()
  if n_models == 0 {
    return None
  }
  let first_atoms = models[0].get_atoms()
  let n_atoms = first_atoms.length()
  if n_atoms == 0 {
    return None
  }
  let atom_names = Array::make(n_atoms, "")
  let resnames = Array::make(n_atoms, "")
  let chainids = Array::make(n_atoms, ' ')
  let resseqs = Array::make(n_atoms, 0)
  let elements = Array::make(n_atoms, "")
  let bfactors = Array::make(n_atoms, 0.0)
  let occupancies = Array::make(n_atoms, 1.0)
  let hetflags = Array::make(n_atoms, false)
  for i = 0; i < n_atoms; i = i + 1 {
    let a = first_atoms[i]
    atom_names[i] = a.name
    resnames[i] = a.resname
    chainids[i] = a.chainid
    resseqs[i] = a.resseq
    elements[i] = a.element
    bfactors[i] = a.bfactor
    occupancies[i] = a.occupancy
  }
  // Determine het flags from residues in first model
  let mut atom_idx = 0
  for chain in models[0].get_chains() {
    for res in chain.get_residues() {
      let is_het = res.is_het()
      for _ in 0..<res.get_num_atoms() {
        if atom_idx < n_atoms {
          hetflags[atom_idx] = is_het
        }
        atom_idx = atom_idx + 1
      }
    }
  }
  let frames : Array[Array[Vector3]] = Array(capacity=n_models)
  for m_idx = 0; m_idx < n_models; m_idx = m_idx + 1 {
    let atoms = models[m_idx].get_atoms()
    if atoms.length() != n_atoms {
      return None
    }
    let frame : Array[Vector3] = Array(capacity=n_atoms)
    for i = 0; i < n_atoms; i = i + 1 {
      frame.push(atoms[i].coord)
    }
    frames.push(frame)
  }
  Some(PDBTrajectory::{
    n_frames: n_models,
    n_atoms,
    atom_names,
    resnames,
    chainids,
    resseqs,
    elements,
    bfactors,
    occupancies,
    hetflags,
    frames,
  })
}

///|
/// Parse a PDB string directly into a PDBTrajectory.
pub fn parse_pdb_trajectory(content : String) -> PDBTrajectory? {
  let structure = parse_pdb(content)
  PDBTrajectory::from_structure(structure)
}

///|
/// Build a Model from a flat atom list by grouping consecutive atoms
/// with the same (chainid, resseq) into residues and chains.
fn build_model_from_atoms(model_id : Int, atoms : Array[Atom]) -> Model {
  let chains : Array[Chain] = Array(capacity=0)
  let mut chain_residues : Array[Residue] = Array(capacity=0)
  let mut res_atoms : Array[Atom] = Array(capacity=0)
  let mut curr_chainid = ' '
  let mut curr_resseq = -1
  let mut curr_resname = ""
  let mut curr_het = false
  let mut started = false
  for i = 0; i < atoms.length(); i = i + 1 {
    let a = atoms[i]
    if started && (a.chainid != curr_chainid || a.resseq != curr_resseq) {
      let r = Residue::new(
        resname=curr_resname,
        chainid=curr_chainid,
        resseq=curr_resseq,
        atoms=res_atoms,
        hetfield=if curr_het { "H" } else { " " },
      )
      chain_residues.push(r)
      res_atoms = Array(capacity=0)
      if a.chainid != curr_chainid {
        let c = Chain::new(id=curr_chainid, residues=chain_residues)
        chains.push(c)
        chain_residues = Array(capacity=0)
      }
    }
    curr_chainid = a.chainid
    curr_resseq = a.resseq
    curr_resname = a.resname
    curr_het = false
    started = true
    res_atoms.push(a)
  }
  if res_atoms.length() > 0 {
    let r = Residue::new(
      resname=curr_resname,
      chainid=curr_chainid,
      resseq=curr_resseq,
      atoms=res_atoms,
      hetfield=if curr_het { "H" } else { " " },
    )
    chain_residues.push(r)
  }
  if chain_residues.length() > 0 {
    let c = Chain::new(id=curr_chainid, residues=chain_residues)
    chains.push(c)
  }
  Model::new(id=model_id, chains~)
}

///|
/// Convert a PDBTrajectory back to a Structure.
pub fn PDBTrajectory::to_structure(self : PDBTrajectory) -> Structure {
  let models : Array[Model] = Array(capacity=self.n_frames)
  for f_idx = 0; f_idx < self.n_frames; f_idx = f_idx + 1 {
    let atoms : Array[Atom] = Array(capacity=self.n_atoms)
    for i = 0; i < self.n_atoms; i = i + 1 {
      let a = Atom::new(
        name=self.atom_names[i],
        coord=self.frames[f_idx][i],
        resname=self.resnames[i],
        chainid=self.chainids[i],
        resseq=self.resseqs[i],
        occupancy=self.occupancies[i],
        bfactor=self.bfactors[i],
        element=self.elements[i],
        charge=0,
      )
      atoms.push(a)
    }
    let model = build_model_from_atoms(f_idx + 1, atoms)
    models.push(model)
  }
  Structure::new(id="", models~)
}

///|
/// Number of frames (models).
pub fn PDBTrajectory::num_frames(self : PDBTrajectory) -> Int {
  self.n_frames
}

///|
/// Number of atoms per frame.
pub fn PDBTrajectory::num_atoms(self : PDBTrajectory) -> Int {
  self.n_atoms
}

///|
/// Get coordinates of one frame: Array[Vector3] of length n_atoms.
pub fn PDBTrajectory::get_frame(self : PDBTrajectory, i : Int) -> Array[Vector3] {
  self.frames[i]
}

///|
/// Get coordinates of one atom across all frames: Array[Vector3] of length n_frames.
pub fn PDBTrajectory::get_atom_coords(self : PDBTrajectory, atom_idx : Int) -> Array[Vector3] {
  let result : Array[Vector3] = Array(capacity=self.n_frames)
  for f = 0; f < self.n_frames; f = f + 1 {
    result.push(self.frames[f][atom_idx])
  }
  result
}

///|
/// Get all coordinates as [frame][atom] = Vector3.
pub fn PDBTrajectory::get_coords(self : PDBTrajectory) -> Array[Array[Vector3]] {
  self.frames
}

///|
/// Get all coordinates as flat double array [frame][atom][3].
pub fn PDBTrajectory::get_coords_flat(self : PDBTrajectory) -> Array[Array[Array[Double]]] {
  let result : Array[Array[Array[Double]]] = Array(capacity=self.n_frames)
  for f = 0; f < self.n_frames; f = f + 1 {
    let frame : Array[Array[Double]] = Array(capacity=self.n_atoms)
    for a = 0; a < self.n_atoms; a = a + 1 {
      let c = self.frames[f][a]
      frame.push([c.x, c.y, c.z])
    }
    result.push(frame)
  }
  result
}

///|
/// Extract a subset of frames (Python-style slice with step).
///
/// `slice(start, end, step)` — `end` is exclusive.
pub fn PDBTrajectory::slice(
  self : PDBTrajectory,
  start : Int,
  end : Int,
  step : Int,
) -> PDBTrajectory {
  let s = if start < 0 { 0 } else { start }
  let e = if end > self.n_frames { self.n_frames } else { end }
  let st = if step <= 0 { 1 } else { step }
  let selected : Array[Array[Vector3]] = Array(capacity=0)
  let mut i = s
  while i < e {
    selected.push(self.frames[i])
    i = i + st
  }
  PDBTrajectory::{
    n_frames: selected.length(),
    n_atoms: self.n_atoms,
    atom_names: self.atom_names,
    resnames: self.resnames,
    chainids: self.chainids,
    resseqs: self.resseqs,
    elements: self.elements,
    bfactors: self.bfactors,
    occupancies: self.occupancies,
    hetflags: self.hetflags,
    frames: selected,
  }
}

///|
/// Mean coordinates across all frames (reference structure).
pub fn PDBTrajectory::mean_coords(self : PDBTrajectory) -> Array[Vector3] {
  let result : Array[Vector3] = Array::make(self.n_atoms, Vector3::zero())
  for f = 0; f < self.n_frames; f = f + 1 {
    for i = 0; i < self.n_atoms; i = i + 1 {
      result[i] = Vector3::new(
        result[i].x + self.frames[f][i].x,
        result[i].y + self.frames[f][i].y,
        result[i].z + self.frames[f][i].z,
      )
    }
  }
  let nf = self.n_frames.to_double()
  for i = 0; i < self.n_atoms; i = i + 1 {
    result[i] = Vector3::new(result[i].x / nf, result[i].y / nf, result[i].z / nf)
  }
  result
}

///|
/// RMSD between two coordinate sets (same length).
fn traj_rmsd(a : Array[Vector3], b : Array[Vector3]) -> Double {
  let n = a.length()
  if n == 0 || b.length() < n {
    return 0.0
  }
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = a[i].x - b[i].x
    let dy = a[i].y - b[i].y
    let dz = a[i].z - b[i].z
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  (sum_sq / n.to_double()).sqrt()
}

///|
/// Pairwise RMSD matrix between all frames.
///
/// Returns `n_frames × n_frames` matrix.  Diagonal is zero.
pub fn PDBTrajectory::rmsd_matrix(self : PDBTrajectory) -> Array[Array[Double]] {
  let n = self.n_frames
  let mat : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let r = traj_rmsd(self.frames[i], self.frames[j])
      mat[i][j] = r
      mat[j][i] = r
    }
  }
  mat
}

///|
/// Per-atom RMSF across all frames (relative to mean structure).
///
/// RMSF_i = √( (1/N) Σ_f |r_i(f) - ⟨r_i⟩|² )
pub fn PDBTrajectory::rmsf(self : PDBTrajectory) -> Array[Double] {
  let mean = self.mean_coords()
  let result = Array::make(self.n_atoms, 0.0)
  for f = 0; f < self.n_frames; f = f + 1 {
    for i = 0; i < self.n_atoms; i = i + 1 {
      let dx = self.frames[f][i].x - mean[i].x
      let dy = self.frames[f][i].y - mean[i].y
      let dz = self.frames[f][i].z - mean[i].z
      result[i] = result[i] + dx * dx + dy * dy + dz * dz
    }
  }
  let nf = self.n_frames.to_double()
  for i = 0; i < self.n_atoms; i = i + 1 {
    result[i] = (result[i] / nf).sqrt()
  }
  result
}

///|
/// Displacement vectors from mean structure for each frame.
///
/// Returns `n_frames × n_atoms` array of `Vector3` displacements.
pub fn PDBTrajectory::deltag(self : PDBTrajectory) -> Array[Array[Vector3]] {
  let mean = self.mean_coords()
  let result : Array[Array[Vector3]] = Array(capacity=self.n_frames)
  for f = 0; f < self.n_frames; f = f + 1 {
    let frame : Array[Vector3] = Array(capacity=self.n_atoms)
    for i = 0; i < self.n_atoms; i = i + 1 {
      frame.push(Vector3::new(
        self.frames[f][i].x - mean[i].x,
        self.frames[f][i].y - mean[i].y,
        self.frames[f][i].z - mean[i].z,
      ))
    }
    result.push(frame)
  }
  result
}

///|
/// Superpose all frames to a reference using Kabsch superposition.
///
/// If `ref_coords` is `None`, uses the first frame as reference.
/// Returns a new PDBTrajectory with aligned coordinates.
pub fn PDBTrajectory::superpose_frames(
  self : PDBTrajectory,
  ref_coords? : Array[Vector3] = [],
) -> PDBTrajectory {
  let ref = if ref_coords.length() > 0 {
    ref_coords
  } else if self.n_frames > 0 {
    self.frames[0]
  } else {
    []
  }
  if ref.length() == 0 || ref.length() != self.n_atoms {
    return self
  }
  let aligned : Array[Array[Vector3]] = Array(capacity=self.n_frames)
  for f = 0; f < self.n_frames; f = f + 1 {
    let (rot, trans) = vector_superposition_transform(self.frames[f], ref)
    let frame : Array[Vector3] = Array(capacity=self.n_atoms)
    for i = 0; i < self.n_atoms; i = i + 1 {
      frame.push(rot.transform(self.frames[f][i]).add(trans))
    }
    aligned.push(frame)
  }
  PDBTrajectory::{
    n_frames: self.n_frames,
    n_atoms: self.n_atoms,
    atom_names: self.atom_names,
    resnames: self.resnames,
    chainids: self.chainids,
    resseqs: self.resseqs,
    elements: self.elements,
    bfactors: self.bfactors,
    occupancies: self.occupancies,
    hetflags: self.hetflags,
    frames: aligned,
  }
}

///|
/// Write the trajectory as a multi-model PDB string.
pub fn PDBTrajectory::write_pdb(self : PDBTrajectory) -> String {
  let buf = StringBuilder()
  for f_idx = 0; f_idx < self.n_frames; f_idx = f_idx + 1 {
    buf.write_string("MODEL     ")
    buf.write_string((f_idx + 1).to_string())
    buf.write_char('\n')
    let mut serial = 1
    for i = 0; i < self.n_atoms; i = i + 1 {
      let coord = self.frames[f_idx][i]
      let line = write_traj_atom_line(
        serial,
        self.atom_names[i],
        self.resnames[i],
        self.chainids[i],
        self.resseqs[i],
        coord,
        self.occupancies[i],
        self.bfactors[i],
        self.elements[i],
        self.hetflags[i],
      )
      buf.write_string(line)
      buf.write_char('\n')
      serial = serial + 1
    }
    buf.write_string("ENDMDL\n")
  }
  buf.write_string("END\n")
  buf.to_string()
}

///|
/// Format a single ATOM/HETATM line for trajectory output.
fn write_traj_atom_line(
  serial : Int,
  name : String,
  resname : String,
  chainid : Char,
  resseq : Int,
  coord : Vector3,
  occupancy : Double,
  bfactor : Double,
  element : String,
  is_het : Bool,
) -> String {
  let buf = StringBuilder()
  buf.write_string(if is_het { "HETATM" } else { "ATOM  " })
  // Serial: cols 7-11, right-justified
  let serial_str = if serial > 99999 { "99999" } else { serial.to_string() }
  buf.write_string(pad_left(serial_str, 5, ' '))
  buf.write_char(' ')
  // Atom name: cols 13-16
  buf.write_string(pad_right(name, 4, ' '))
  buf.write_char(' ')
  buf.write_char(' ')
  // Residue name: cols 18-20
  buf.write_string(pad_right(resname, 3, ' '))
  buf.write_char(' ')
  // Chain ID: col 22
  buf.write_char(chainid)
  // Residue seq: cols 23-26, right-justified
  let resseq_str = if resseq > 9999 { "9999" } else { if resseq < 0 { "0" } else { resseq.to_string() } }
  buf.write_string(pad_left(resseq_str, 4, ' '))
  buf.write_char(' ')
  // X: cols 31-38, %8.3f
  buf.write_string(format_float(coord.x, 8, 3))
  // Y: cols 39-46
  buf.write_string(format_float(coord.y, 8, 3))
  // Z: cols 47-54
  buf.write_string(format_float(coord.z, 8, 3))
  // Occupancy: cols 55-60, %6.2f
  buf.write_string(format_float(occupancy, 6, 2))
  // B-factor: cols 61-66, %6.2f
  buf.write_string(format_float(bfactor, 6, 2))
  buf.write_char(' ')
  buf.write_char(' ')
  buf.write_char(' ')
  buf.write_char(' ')
  // Element: cols 77-78
  buf.write_string(pad_left(element, 2, ' '))
  buf.to_string()
}

///|
/// Format a float as fixed-width with specified precision.
/// Mimics C printf "%w.pf".
fn format_float(val : Double, width : Int, precision : Int) -> String {
  // Handle sign
  let neg = val < 0.0
  let abs_val = if neg { -val } else { val }
  // Scale and round
  let mut scale = 1.0
  let mut pi = 0
  while pi < precision {
    scale = scale * 10.0
    pi = pi + 1
  }
  let scaled = (abs_val * scale + 0.5).floor()
  let int_part = (scaled / scale).floor()
  let frac_part = scaled - int_part * scale
  let int_str = int_part.to_int().to_string()
  // Build fractional part with leading zeros
  let frac_int = frac_part.to_int()
  let frac_str = if precision > 0 {
    let fs = frac_int.to_string()
    pad_left(fs, precision, '0')
  } else {
    ""
  }
  let num_str = if precision > 0 {
    int_str + "." + frac_str
  } else {
    int_str
  }
  let signed = if neg { "-" + num_str } else { num_str }
  pad_left(signed, width, ' ')
}
