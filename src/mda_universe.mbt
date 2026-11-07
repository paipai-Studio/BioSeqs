///|
/// MDAnalysis-style `Universe` / `AtomGroup` object model.
///
/// A `Universe` bundles a static topology (per-atom metadata) with a
/// `Trajectory` (coordinate frames) and a current-frame cursor, mirroring
/// `MDAnalysis.Universe(topology, trajectory)`. Builders cover the classic
/// workflows: PSF topology + DCD trajectory, multi-MODEL PDB, PQR structure
/// (+ optional trajectory). An `AtomGroup` is a lightweight selection view
/// over a `Universe` exposing positions, centers of geometry / mass, radius
/// of gyration, charges and sub-selections.

///|
/// Per-atom topology metadata (one entry per atom in the system).
pub(all) struct UnivAtom {
  name : String // atom name: CA, N, O, ...
  resname : String // residue name: ALA, HOH, ...
  resid : Int32 // residue id
  chain : String // chain identifier
  segid : String // segment identifier (PSF)
  element : String // element symbol: C, N, O, ...
  atom_type : String // force-field atom type (PSF)
  charge : Double // partial charge (e)
  mass : Double // atomic mass (amu); <= 0 means unknown
  is_protein : Bool // standard amino-acid residue
  is_water : Bool // water-like solvent residue
  is_backbone : Bool // protein backbone name in {N, CA, C, O}
}

///|
/// MDAnalysis-style Universe: topology + trajectory + current frame cursor.
pub struct Universe {
  atoms : Array[UnivAtom]
  traj : Trajectory
  mut current : Int32
}

///|
/// Create a Universe from atom metadata and a coordinate trajectory.
pub fn Universe::new(atoms : Array[UnivAtom], traj : Trajectory) -> Universe {
  Universe::{ atoms, traj, current: 0, }
}

///|
/// Number of atoms in the topology.
pub fn Universe::n_atoms(self : Universe) -> Int32 {
  self.atoms.length()
}

///|
/// Number of frames in the trajectory.
pub fn Universe::n_frames(self : Universe) -> Int32 {
  self.traj.num_frames()
}

///|
/// Current frame index (0-based).
pub fn Universe::frame_index(self : Universe) -> Int32 {
  self.current
}

///|
/// Set the current frame (clamped to the valid range).
pub fn Universe::set_frame(self : Universe, i : Int32) -> Unit {
  let n = self.traj.num_frames()
  if n == 0 {
    return
  }
  self.current = if i < 0 { 0 } else if i >= n { n - 1 } else { i }
}

///|
/// Atom metadata array (read-only view).
pub fn Universe::atom_info(self : Universe) -> Array[UnivAtom] {
  self.atoms
}

///|
/// Coordinates of all atoms at the current frame.
pub fn Universe::positions(self : Universe) -> Array[Array[Double]] {
  self.positions_at(self.current)
}

///|
/// Coordinates of all atoms at frame `i` (empty array when out of range).
pub fn Universe::positions_at(
  self : Universe,
  i : Int32,
) -> Array[Array[Double]] {
  let n = self.traj.num_frames()
  if i < 0 || i >= n {
    return []
  }
  self.traj.get_frame(i)
}

///|
/// Select atoms with the extended MDAnalysis-style selection language.
/// See `mda_select_universe` for the supported keywords.
pub fn Universe::select_atoms(
  self : Universe,
  expr : String,
  box~ : Array[Double]? = None,
) -> AtomGroup {
  let idx = mda_select_universe(
    expr,
    self.atoms,
    positions=Some(self.positions()),
    box~,
  )
  AtomGroup::new(self, idx)
}

///|
/// AtomGroup covering all atoms of the Universe.
pub fn Universe::all(self : Universe) -> AtomGroup {
  let idx : Array[Int32] = Array(capacity=0)
  for i = 0; i < self.atoms.length(); i = i + 1 {
    idx.push(i)
  }
  AtomGroup::new(self, idx)
}

///|
/// Charges of all atoms.
pub fn Universe::charges(self : Universe) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for a in self.atoms {
    out.push(a.charge)
  }
  out
}

///|
/// Masses of all atoms.
pub fn Universe::masses(self : Universe) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for a in self.atoms {
    out.push(a.mass)
  }
  out
}

///|
/// MDAnalysis-style AtomGroup: an index view over a Universe.
pub struct AtomGroup {
  universe : Universe
  indices : Array[Int32]
}

///|
/// Create an AtomGroup from a Universe and a list of atom indices.
pub fn AtomGroup::new(universe : Universe, indices : Array[Int32]) -> AtomGroup {
  AtomGroup::{ universe, indices, }
}

///|
/// Number of atoms in the group.
pub fn AtomGroup::n_atoms(self : AtomGroup) -> Int32 {
  self.indices.length()
}

///|
/// Atom indices (into the parent Universe, 0-based).
pub fn AtomGroup::indices(self : AtomGroup) -> Array[Int32] {
  self.indices
}

///|
/// Coordinates of the group atoms at the current frame.
pub fn AtomGroup::positions(self : AtomGroup) -> Array[Array[Double]] {
  let coords = self.universe.positions()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i in self.indices {
    if i >= 0 && i < coords.length() {
      out.push(coords[i])
    }
  }
  out
}

///|
/// Coordinates of the group atoms at frame `i`.
pub fn AtomGroup::positions_at(
  self : AtomGroup,
  i : Int32,
) -> Array[Array[Double]] {
  let coords = self.universe.positions_at(i)
  let out : Array[Array[Double]] = Array(capacity=0)
  for idx in self.indices {
    if idx >= 0 && idx < coords.length() {
      out.push(coords[idx])
    }
  }
  out
}

///|
/// Atom metadata of the group members.
pub fn AtomGroup::atom_info(self : AtomGroup) -> Array[UnivAtom] {
  let out : Array[UnivAtom] = Array(capacity=0)
  for i in self.indices {
    if i >= 0 && i < self.universe.atoms.length() {
      out.push(self.universe.atoms[i])
    }
  }
  out
}

///|
/// Atom names of the group.
pub fn AtomGroup::names(self : AtomGroup) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for a in self.atom_info() {
    out.push(a.name)
  }
  out
}

///|
/// Residue names of the group.
pub fn AtomGroup::resnames(self : AtomGroup) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for a in self.atom_info() {
    out.push(a.resname)
  }
  out
}

///|
/// Residue ids of the group.
pub fn AtomGroup::resids(self : AtomGroup) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  for a in self.atom_info() {
    out.push(a.resid)
  }
  out
}

///|
/// Charges of the group atoms.
pub fn AtomGroup::charges(self : AtomGroup) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for a in self.atom_info() {
    out.push(a.charge)
  }
  out
}

///|
/// Masses of the group atoms (non-positive mass counts as 1.0).
pub fn AtomGroup::masses(self : AtomGroup) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for a in self.atom_info() {
    out.push(if a.mass > 0.0 { a.mass } else { 1.0 })
  }
  out
}

///|
/// Sum of partial charges over the group.
pub fn AtomGroup::total_charge(self : AtomGroup) -> Double {
  let mut s = 0.0
  for c in self.charges() {
    s = s + c
  }
  s
}

///|
/// Sum of masses over the group (non-positive mass counts as 1.0).
pub fn AtomGroup::total_mass(self : AtomGroup) -> Double {
  let mut s = 0.0
  for m in self.masses() {
    s = s + m
  }
  s
}

///|
/// Geometric center of the group at the current frame.
pub fn AtomGroup::center_of_geometry(self : AtomGroup) -> Array[Double] {
  let pos = self.positions()
  let c = [0.0, 0.0, 0.0]
  for p in pos {
    c[0] = c[0] + p[0]
    c[1] = c[1] + p[1]
    c[2] = c[2] + p[2]
  }
  let n = pos.length().to_double()
  if n > 0.0 {
    c[0] = c[0] / n
    c[1] = c[1] / n
    c[2] = c[2] / n
  }
  c
}

///|
/// Mass-weighted center of the group at the current frame
/// (`AtomGroup.center_of_mass`). Non-positive masses count as 1.0.
pub fn AtomGroup::center_of_mass(self : AtomGroup) -> Array[Double] {
  let pos = self.positions()
  let ms = self.masses()
  let c = [0.0, 0.0, 0.0]
  let mut mt = 0.0
  for k = 0; k < pos.length(); k = k + 1 {
    let m = if k < ms.length() { ms[k] } else { 1.0 }
    c[0] = c[0] + m * pos[k][0]
    c[1] = c[1] + m * pos[k][1]
    c[2] = c[2] + m * pos[k][2]
    mt = mt + m
  }
  if mt > 0.0 {
    c[0] = c[0] / mt
    c[1] = c[1] / mt
    c[2] = c[2] / mt
  }
  c
}

///|
/// Mass-weighted radius of gyration about the center of mass
/// (`AtomGroup.radius_of_gyration`).
pub fn AtomGroup::radius_of_gyration(self : AtomGroup) -> Double {
  let pos = self.positions()
  let ms = self.masses()
  let com = self.center_of_mass()
  let mut acc = 0.0
  let mut mt = 0.0
  for k = 0; k < pos.length(); k = k + 1 {
    let m = if k < ms.length() { ms[k] } else { 1.0 }
    let dx = pos[k][0] - com[0]
    let dy = pos[k][1] - com[1]
    let dz = pos[k][2] - com[2]
    acc = acc + m * (dx * dx + dy * dy + dz * dz)
    mt = mt + m
  }
  if mt > 0.0 {
    (acc / mt).sqrt()
  } else {
    0.0
  }
}

///|
/// Sub-select within the group using the extended selection language.
/// The sub-expression is evaluated against the group's own atoms
/// (`bynum` indices are 1-based within the group).
pub fn AtomGroup::select(
  self : AtomGroup,
  expr : String,
  box~ : Array[Double]? = None,
) -> AtomGroup {
  let sub_atoms = self.atom_info()
  let idx = mda_select_universe(
    expr,
    sub_atoms,
    positions=Some(self.positions()),
    box~,
  )
  let mapped : Array[Int32] = Array(capacity=0)
  for i in idx {
    if i >= 0 && i < self.indices.length() {
      mapped.push(self.indices[i])
    }
  }
  AtomGroup::new(self.universe, mapped)
}

// ---------------------------------------------------------------------------
// Builders
// ---------------------------------------------------------------------------

///|
/// Classify a residue name: protein (standard amino acid), water-like
/// solvent (HOH/WAT/TIP/TIP3/TIP4P/SPC/SOL) or other.
pub fn mda_is_water_resname(resname : String) -> Bool {
  let r = resname.to_upper()
  r == "HOH" ||
  r == "WAT" ||
  r == "TIP" ||
  r == "TIP3" ||
  r == "TIP4P" ||
  r == "SPC" ||
  r == "SOL"
}

///|
/// Infer the element symbol from an atom name, e.g. "CA" -> "C" (alpha
/// carbon; "CA" is excluded from the two-letter set to avoid clashing with
/// calcium), "CL1" -> "Cl", "FE" -> "Fe", "1HG2" -> "H".
pub fn mda_infer_element(name : String) -> String {
  let upper_a = 'A'.to_int()
  let upper_z = 'Z'.to_int()
  let lower_a = 'a'.to_int()
  let lower_z = 'z'.to_int()
  let letters : Array[Int32] = Array(capacity=0)
  for i = 0; i < name.length(); i = i + 1 {
    let v = name.unsafe_get(i).to_int()
    if (v >= upper_a && v <= upper_z) || (v >= lower_a && v <= lower_z) {
      letters.push(v)
    } else if letters.length() > 0 {
      break
    }
  }
  if letters.length() == 0 {
    return ""
  }
  let c0 = letters[0]
  let c0_upper = if c0 >= lower_a && c0 <= lower_z { c0 - 32 } else { c0 }
  if letters.length() >= 2 {
    let c1 = letters[1]
    let c1_lower = if c1 >= upper_a && c1 <= upper_z { c1 + 32 } else { c1 }
    // normalized two-letter symbol: first upper, second lower
    let sb = StringBuilder()
    sb.write_char(c0_upper.to_uint16().unsafe_to_char())
    sb.write_char(c1_lower.to_uint16().unsafe_to_char())
    let two = sb.to_string()
    let two_upper = two.to_upper()
    if two_upper == "CL" ||
      two_upper == "BR" ||
      two_upper == "MG" ||
      two_upper == "MN" ||
      two_upper == "FE" ||
      two_upper == "ZN" ||
      two_upper == "SE" ||
      two_upper == "NA" ||
      two_upper == "SI" ||
      two_upper == "CU" {
      return two
    }
  }
  let sb1 = StringBuilder()
  sb1.write_char(c0_upper.to_uint16().unsafe_to_char())
  sb1.to_string()
}

///|
/// Approximate atomic mass by element symbol (amu); unknown elements
/// return 12.011 (carbon).
pub fn mda_element_mass(element : String) -> Double {
  let e = element.to_upper()
  if e == "H" {
    1.008
  } else if e == "C" {
    12.011
  } else if e == "N" {
    14.007
  } else if e == "O" {
    15.999
  } else if e == "P" {
    30.974
  } else if e == "S" {
    32.06
  } else if e == "NA" {
    22.99
  } else if e == "MG" {
    24.305
  } else if e == "CL" {
    35.45
  } else if e == "K" {
    39.098
  } else if e == "CA" {
    40.078
  } else if e == "FE" {
    55.845
  } else if e == "ZN" {
    65.38
  } else if e == "SE" {
    78.971
  } else if e == "F" {
    18.998
  } else if e == "BR" {
    79.904
  } else if e == "I" {
    126.904
  } else {
    12.011
  }
}

///|
fn mda_univ_classify_atom(
  name : String,
  resname : String,
  resid : Int32,
  chain : String,
  segid : String,
  element : String,
  atom_type : String,
  charge : Double,
  mass : Double,
) -> UnivAtom {
  let elem = if element.length() > 0 {
    element
  } else {
    mda_infer_element(name)
  }
  let is_protein = psf_is_protein_resname(resname)
  let bare = name.trim().to_owned().to_upper()
  UnivAtom::{
    name,
    resname,
    resid,
    chain,
    segid,
    element: elem,
    atom_type,
    charge,
    mass,
    is_protein,
    is_water: mda_is_water_resname(resname),
    is_backbone: is_protein &&
    (bare == "N" || bare == "CA" || bare == "C" || bare == "O"),
  }
}

///|
/// Build a Universe from a PSF topology text and a DCD trajectory byte
/// string (`MDAnalysis.Universe(psf, dcd)`).
pub fn universe_from_psf_dcd(
  psf_text : String,
  dcd_bytes : String,
) -> Universe raise PsfError {
  universe_from_psf_traj(psf_text, parse_dcd(dcd_bytes))
}

///|
/// Build a Universe from a PSF topology text and an in-memory Trajectory.
pub fn universe_from_psf_traj(
  psf_text : String,
  traj : Trajectory,
) -> Universe raise PsfError {
  let topo = psf_parse(psf_text)
  let atoms : Array[UnivAtom] = Array(capacity=0)
  for a in topo.atoms {
    atoms.push(
      mda_univ_classify_atom(
        a.atom_name,
        a.resname,
        a.resid,
        a.segid,
        a.segid,
        mda_infer_element(a.atom_name),
        a.atom_type,
        a.charge,
        a.mass,
      ),
    )
  }
  Universe::new(atoms, traj)
}

///|
/// Build a Universe from a multi-MODEL PDB text. Each `MODEL ... ENDMDL`
/// block becomes one frame; a file without MODEL records yields a single
/// frame. Atom metadata is taken from the first frame.
pub fn universe_from_pdb_text(text : String) -> Universe {
  let atoms : Array[UnivAtom] = Array(capacity=0)
  let frames : Array[Array[Array[Double]]] = Array(capacity=0)
  let mut cur : Array[Array[Double]]? = None
  let mut collecting = true
  let lines = text.split("\n")
  for line in lines {
    let l = line.to_owned()
    if l.length() >= 5 && mda_univ_slice(l, 0, 5) == "MODEL" {
      cur = Some([])
      continue
    }
    if l.length() >= 6 && mda_univ_slice(l, 0, 6) == "ENDMDL" {
      match cur {
        Some(f) => {
          frames.push(f)
          cur = None
        }
        None => ()
      }
      collecting = false
      continue
    }
    if l.length() >= 6 &&
      (
        mda_univ_slice(l, 0, 6) == "ATOM  " ||
        mda_univ_slice(l, 0, 6) == "HETATM"
      ) {
      if cur is None {
        cur = Some([])
      }
      let resname = mda_univ_slice(l, 17, 3).trim().to_owned()
      let name = mda_univ_slice(l, 12, 4).trim().to_owned()
      let chain = mda_univ_slice(l, 21, 1).trim().to_owned()
      let resid_w = mda_univ_slice(l, 22, 4).trim().to_owned()
      let resid_v = parse_int(resid_w)
      let x = mda_univ_col_double(l, 30, 8)
      let y = mda_univ_col_double(l, 38, 8)
      let z = mda_univ_col_double(l, 46, 8)
      let elem_col = mda_univ_slice(l, 76, 2).trim().to_owned()
      if collecting {
        atoms.push(
          mda_univ_classify_atom(
            name,
            resname,
            resid_v,
            chain,
            "",
            elem_col,
            "",
            0.0,
            mda_element_mass(
              if elem_col.length() > 0 {
                elem_col
              } else {
                mda_infer_element(name)
              },
            ),
          ),
        )
      }
      match cur {
        Some(f) => f.push([x, y, z])
        None => ()
      }
      continue
    }
  }
  match cur {
    Some(f) => if f.length() > 0 { frames.push(f) }
    None => ()
  }
  Universe::new(atoms, Trajectory::new(frames))
}

///|
/// Build a Universe from a PQR text (charges + radii). The coordinates of
/// the PQR become the single frame unless a trajectory is supplied.
pub fn universe_from_pqr_text(
  pqr_text : String,
  traj~ : Trajectory? = None,
) -> Universe raise PqrError {
  let s = pqr_parse(pqr_text)
  let atoms : Array[UnivAtom] = Array(capacity=0)
  for a in s.atoms {
    atoms.push(
      mda_univ_classify_atom(
        a.name,
        a.resname,
        a.resid,
        a.chain,
        "",
        mda_infer_element(a.name),
        "",
        a.charge,
        mda_element_mass(mda_infer_element(a.name)),
      ),
    )
  }
  let tr = match traj {
    Some(t) => t
    None => {
      let frames : Array[Array[Array[Double]]] = Array(capacity=0)
      frames.push(pqr_coords(s))
      Trajectory::new(frames)
    }
  }
  Universe::new(atoms, tr)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

///|
fn mda_univ_slice(s : String, start : Int32, len : Int32) -> String {
  let n = s.length()
  if start >= n || len <= 0 {
    return ""
  }
  let out = StringBuilder()
  let mut i = start
  let end = if start + len < n { start + len } else { n }
  while i < end {
    out.write_char(s.unsafe_get(i).unsafe_to_char())
    i = i + 1
  }
  out.to_string()
}

///|
fn mda_univ_col_double(s : String, start : Int32, len : Int32) -> Double {
  let w = mda_univ_slice(s, start, len).trim().to_owned()
  match parse_double(w) {
    Some(v) => v
    None => 0.0
  }
}
