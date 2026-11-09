/// omm_bridge.mbt — bridges between the OmmSystem world and the MDAnalysis
/// Universe world, plus a full AMBER prmtop/inpcrd parser.
///
/// Previously `omm_parse_prmtop` (omm_amber.mbt) returned empty arrays; this
/// file provides `omm_parse_prmtop_full`, a complete FLAG-based parser, and
/// conversion helpers in both directions:
/// - `omm_system_from_universe` / `omm_system_to_universe` (OmmSystem <-> Universe)
/// - `omm_universe_from_prmtop` (prmtop [+ inpcrd] -> Universe)
/// - `omm_system_from_prmtop` (prmtop -> OmmSystem with bonds/angles/dihedrals/LJ)

///|
/// Atomic number -> element symbol (1..54 plus common ones).
pub fn omm_element_from_number(z : Int32) -> String {
  match z {
    1 => "H"
    2 => "He"
    3 => "Li"
    4 => "Be"
    5 => "B"
    6 => "C"
    7 => "N"
    8 => "O"
    9 => "F"
    10 => "Ne"
    11 => "Na"
    12 => "Mg"
    13 => "Al"
    14 => "Si"
    15 => "P"
    16 => "S"
    17 => "Cl"
    18 => "Ar"
    19 => "K"
    20 => "Ca"
    26 => "Fe"
    27 => "Co"
    28 => "Ni"
    29 => "Cu"
    30 => "Zn"
    34 => "Se"
    35 => "Br"
    53 => "I"
    _ => "X"
  }
}

///|
/// Guess an element from an atomic mass (amu) using nearest reference mass.
pub fn omm_element_from_mass(mass : Double) -> String {
  let refs : Array[(String, Double)] = [
    ("H", 1.008),
    ("C", 12.011),
    ("N", 14.007),
    ("O", 15.999),
    ("F", 18.998),
    ("Na", 22.990),
    ("P", 30.974),
    ("S", 32.06),
    ("Cl", 35.453),
    ("K", 39.098),
    ("Fe", 55.845),
    ("Zn", 65.38),
    ("Br", 79.904),
    ("I", 126.904),
  ]
  let mut best = "X"
  let mut best_d = 1.0e18
  for (el, m) in refs {
    let d = (mass - m).abs()
    if d < best_d {
      best_d = d
      best = el
    }
  }
  if best_d < 0.65 {
    best
  } else {
    "X"
  }
}

///|
/// OmmSystem -> Universe: build per-atom topology metadata from particle
/// parameters. Elements are guessed from masses; residue/chain info is
/// synthetic (UNK/1..n). `positions` provides the trajectory frame
/// (a zero frame is created when omitted).
pub fn omm_system_to_universe(
  system : OmmSystem,
  positions~ : Array[Array[Double]]? = None,
) -> Universe {
  let n = system.num_particles()
  let atoms : Array[UnivAtom] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let p = system.particles[i]
    atoms.push({
      name: "A" + (i + 1).to_string(),
      resname: "UNK",
      resid: i + 1,
      chain: "A",
      segid: "SYS",
      element: omm_element_from_mass(p.mass),
      atom_type: "",
      charge: p.charge,
      mass: p.mass,
      is_protein: false,
      is_water: false,
      is_backbone: false,
    })
  }
  let frame : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  match positions {
    Some(p) =>
      for i = 0; i < n && i < p.length(); i = i + 1 {
        frame[i] = [p[i][0], p[i][1], p[i][2]]
      }
    None => ()
  }
  Universe::new(atoms, Trajectory::new([frame]))
}

///|
/// Universe -> OmmSystem: particle parameters from atom metadata.
/// LJ sigma/epsilon default to 0 (unknown); override with `sigma`/`epsilon`.
pub fn omm_system_from_universe(
  univ : Universe,
  sigma~ : Double = 0.0,
  epsilon~ : Double = 0.0,
) -> OmmSystem {
  let system = OmmSystem::new()
  for a in univ.atom_info() {
    let mass = if a.mass > 0.0 { a.mass } else { 1.0 }
    ignore(
      system.add_particle(
        particle=OmmParticle::new(mass~, charge=a.charge, sigma~, epsilon~),
      ),
    )
  }
  system
}

// ---------------------------------------------------------------------------
// Full AMBER prmtop parser
// ---------------------------------------------------------------------------

///|
/// Fully parsed AMBER prmtop (FLAG sections). Indices in `bonds`, `angles`,
/// `dihedrals` are 0-based atom indices; the trailing Int32 is a 0-based
/// parameter index into the corresponding *_FORCE_CONSTANT / *_VALUE arrays.
pub(all) struct OmmPrmtop {
  n_atoms : Int32
  n_types : Int32
  atom_names : Array[String]
  atom_types : Array[String]
  type_indices : Array[Int32]
  atomic_numbers : Array[Int32]
  charges : Array[Double] // elementary charge units
  masses : Array[Double]
  residue_labels : Array[String]
  residue_pointers : Array[Int32] // 0-based
  bonds : Array[(Int32, Int32, Int32)]
  angles : Array[(Int32, Int32, Int32, Int32)]
  dihedrals : Array[(Int32, Int32, Int32, Int32, Int32)]
  bond_force_constants : Array[Double]
  bond_equil_values : Array[Double]
  angle_force_constants : Array[Double]
  angle_equil_values : Array[Double]
  dihedral_force_constants : Array[Double]
  dihedral_phases : Array[Double] // degrees
  dihedral_periodicities : Array[Double]
  lj_a : Array[Double]
  lj_b : Array[Double]
}

///|
/// Split one prmtop data line into fixed-width fields.
fn omm_prmtop_fields(line : String, width : Int32) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let len = line.length()
  if width <= 0 {
    return out
  }
  let mut i = 0
  while i < len {
    let end = if i + width < len { i + width } else { len }
    out.push(line[i:end].to_owned().trim().to_owned())
    i = i + width
  }
  out
}

///|
/// Int-array conversion flag: values are integer fields.
fn omm_prmtop_to_ints(lines : Array[String], width : Int32) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  for l in lines {
    for f in omm_prmtop_fields(l, width) {
      if f.length() == 0 {
        continue
      }
      match parse_double(f) {
        Some(v) => out.push(v.to_int())
        None => ()
      }
    }
  }
  out
}

///|
/// Double-array conversion flag: values are floating point fields.
fn omm_prmtop_to_doubles(lines : Array[String], width : Int32) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for l in lines {
    for f in omm_prmtop_fields(l, width) {
      if f.length() == 0 {
        continue
      }
      match parse_double(f) {
        Some(v) => out.push(v)
        None => ()
      }
    }
  }
  out
}

///|
/// String-array (a4) flag: fixed 4-char names.
fn omm_prmtop_to_strings(lines : Array[String], width : Int32) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for l in lines {
    for f in omm_prmtop_fields(l, width) {
      if f.length() == 0 {
        continue
      }
      out.push(f)
    }
  }
  out
}

///|
/// Parse a complete AMBER prmtop file into an OmmPrmtop record.
pub fn omm_parse_prmtop_full(content : String) -> OmmPrmtop {
  let lines = content.split("\n").to_array()
  // Collect blocks: flag name -> raw data lines.
  let flags : Array[(String, Array[String])] = Array(capacity=0)
  let mut current : String = ""
  let mut buf : Array[String] = Array(capacity=0)
  for line in lines {
    let l = line.to_owned().trim().to_owned()
    if l.has_prefix("%FLAG") {
      if current.length() > 0 {
        flags.push((current, buf))
      }
      current = if l.length() > 6 {
        l[6:].to_owned().trim().to_owned()
      } else {
        ""
      }
      buf = Array(capacity=0)
      continue
    }
    if l.has_prefix("%FORMAT") {
      continue
    }
    if l.has_prefix("%") {
      continue
    }
    if l.length() == 0 {
      continue
    }
    buf.push(l)
  }
  if current.length() > 0 {
    flags.push((current, buf))
  }
  // Lookup helper.
  fn find_ints(name : String) -> Array[Int32] {
    for (f, ls) in flags {
      if f == name {
        // numeric formats use width 8 in prmtop
        return omm_prmtop_to_ints(ls, 8)
      }
    }
    []
  }
  fn find_doubles(name : String) -> Array[Double] {
    for (f, ls) in flags {
      if f == name {
        return omm_prmtop_to_doubles(ls, 16)
      }
    }
    []
  }
  fn find_strings(name : String) -> Array[String] {
    for (f, ls) in flags {
      if f == name {
        return omm_prmtop_to_strings(ls, 4)
      }
    }
    []
  }
  let pointers = find_ints("POINTERS")
  let n_atoms = if pointers.length() > 0 { pointers[0] } else { 0 }
  let n_types = if pointers.length() > 1 { pointers[1] } else { 0 }
  let raw_charges = find_doubles("CHARGE")
  let charges : Array[Double] = Array::make(raw_charges.length(), 0.0)
  for i = 0; i < raw_charges.length(); i = i + 1 {
    charges[i] = raw_charges[i] / 18.2223 // prmtop internal unit -> e
  }
  // Connectivity sections: triples / quadruples / quintuples of Int32.
  fn grouped_ints(name : String, group : Int32) -> Array[Array[Int32]] {
    let flat = find_ints(name)
    let out : Array[Array[Int32]] = Array(capacity=0)
    let n = flat.length() / group
    for g = 0; g < n; g = g + 1 {
      let row : Array[Int32] = Array::make(group, 0)
      for k = 0; k < group; k = k + 1 {
        row[k] = flat[g * group + k]
      }
      out.push(row)
    }
    out
  }
  let bond_rows = grouped_ints("BONDS_INC_HYDROGEN", 3)
  let bond_rows2 = grouped_ints("BONDS_WITHOUT_HYDROGEN", 3)
  let angle_rows = grouped_ints("ANGLES_INC_HYDROGEN", 4)
  let angle_rows2 = grouped_ints("ANGLES_WITHOUT_HYDROGEN", 4)
  let dih_rows = grouped_ints("DIHEDRALS_INC_HYDROGEN", 5)
  let dih_rows2 = grouped_ints("DIHEDRALS_WITHOUT_HYDROGEN", 5)
  let bonds : Array[(Int32, Int32, Int32)] = Array(capacity=0)
  for r in bond_rows {
    bonds.push((r[0] / 3, r[1] / 3, (r[2] - 1).abs()))
  }
  for r in bond_rows2 {
    bonds.push((r[0] / 3, r[1] / 3, (r[2] - 1).abs()))
  }
  let angles : Array[(Int32, Int32, Int32, Int32)] = Array(capacity=0)
  for r in angle_rows {
    angles.push((r[0] / 3, r[1] / 3, r[2] / 3, (r[3] - 1).abs()))
  }
  for r in angle_rows2 {
    angles.push((r[0] / 3, r[1] / 3, r[2] / 3, (r[3] - 1).abs()))
  }
  let dihedrals : Array[(Int32, Int32, Int32, Int32, Int32)] = Array(capacity=0)
  for r in dih_rows {
    dihedrals.push((r[0] / 3, r[1] / 3, r[2] / 3, r[3] / 3, (r[4] - 1).abs()))
  }
  for r in dih_rows2 {
    dihedrals.push((r[0] / 3, r[1] / 3, r[2] / 3, r[3] / 3, (r[4] - 1).abs()))
  }
  let residue_pointers_raw = find_ints("RESIDUE_POINTER")
  let residue_pointers : Array[Int32] = Array::make(
    residue_pointers_raw.length(),
    0,
  )
  for i = 0; i < residue_pointers_raw.length(); i = i + 1 {
    residue_pointers[i] = residue_pointers_raw[i] - 1
  }
  {
    n_atoms,
    n_types,
    atom_names: find_strings("ATOM_NAME"),
    atom_types: find_strings("AMBER_ATOM_TYPE"),
    type_indices: find_ints("ATOM_TYPE_INDEX"),
    atomic_numbers: find_ints("ATOMIC_NUMBER"),
    charges,
    masses: find_doubles("MASS"),
    residue_labels: find_strings("RESIDUE_LABEL"),
    residue_pointers,
    bonds,
    angles,
    dihedrals,
    bond_force_constants: find_doubles("BOND_FORCE_CONSTANT"),
    bond_equil_values: find_doubles("BOND_EQUIL_VALUE"),
    angle_force_constants: find_doubles("ANGLE_FORCE_CONSTANT"),
    angle_equil_values: find_doubles("ANGLE_EQUIL_VALUE"),
    dihedral_force_constants: find_doubles("DIHEDRAL_FORCE_CONSTANT"),
    dihedral_phases: find_doubles("DIHEDRAL_PHASE"),
    dihedral_periodicities: find_doubles("DIHEDRAL_PERIODICITY"),
    lj_a: find_doubles("LENNARD_JONES_ACOEF"),
    lj_b: find_doubles("LENNARD_JONES_BCOEF"),
  }
}

///|
/// Residue name of atom `i` (0-based) from residue pointers.
fn omm_prmtop_resname(top : OmmPrmtop, i : Int32) -> String {
  let nres = top.residue_labels.length()
  let mut res = "UNK"
  for r = 0; r < nres; r = r + 1 {
    if top.residue_pointers[r] <= i {
      res = top.residue_labels[r]
    }
  }
  res
}

///|
/// Standard protein residue names (for is_protein flags).
fn omm_is_protein_resname(resname : String) -> Bool {
  match resname {
    "ALA"
    | "ARG"
    | "ASN"
    | "ASP"
    | "CYS"
    | "GLN"
    | "GLU"
    | "GLY"
    | "HIS"
    | "ILE"
    | "LEU"
    | "LYS"
    | "MET"
    | "PHE"
    | "PRO"
    | "SER"
    | "THR"
    | "TRP"
    | "TYR"
    | "VAL" => true
    _ => false
  }
}

///|
/// Parse an AMBER inpcrd (ASCII coordinate) file into a coordinate frame.
fn omm_parse_inpcrd(content : String) -> Array[Array[Double]] {
  let lines = content.split("\n").to_array()
  let nonempty : Array[String] = Array(capacity=0)
  for l in lines {
    let t = l.to_owned().trim()
    if t.length() > 0 {
      nonempty.push(l.to_owned())
    }
  }
  let out : Array[Array[Double]] = Array(capacity=0)
  if nonempty.length() < 2 {
    return out
  }
  // Line 0 = title, line 1 = "n_atoms [time]".
  let header = nonempty[1].split(" ").to_array()
  let mut n = 0
  for p in header {
    let t = p.to_owned().trim()
    if t.length() > 0 {
      match parse_double(t.to_owned()) {
        Some(v) => {
          n = v.to_int()
          break
        }
        None => ()
      }
    }
  }
  let mut i = 0
  let mut li = 2
  while i < n && li < nonempty.length() {
    let parts = nonempty[li].split(" ").to_array()
    let vals : Array[Double] = Array(capacity=0)
    for p in parts {
      let t = p.to_owned().trim()
      if t.length() > 0 {
        match parse_double(t.to_owned()) {
          Some(v) => vals.push(v)
          None => ()
        }
      }
    }
    let mut k = 0
    while k + 2 < vals.length() + 1 && out.length() < n {
      if k + 2 < vals.length() {
        out.push([vals[k], vals[k + 1], vals[k + 2]])
        k = k + 3
        i = i + 1
      } else {
        break
      }
    }
    li = li + 1
  }
  out
}

///|
/// Build a Universe from a prmtop (topology) plus optional inpcrd
/// (coordinates). Without inpcrd the trajectory has zero frames.
pub fn omm_universe_from_prmtop(
  prmtop_content : String,
  inpcrd~ : String = "",
) -> Universe {
  let top = omm_parse_prmtop_full(prmtop_content)
  let n = top.n_atoms
  let atoms : Array[UnivAtom] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let resname = omm_prmtop_resname(top, i)
    let name = if i < top.atom_names.length() { top.atom_names[i] } else { "X" }
    let atom_type = if i < top.atom_types.length() {
      top.atom_types[i]
    } else {
      ""
    }
    let z = if i < top.atomic_numbers.length() {
      top.atomic_numbers[i]
    } else {
      0
    }
    let element = if z > 0 {
      omm_element_from_number(z)
    } else {
      omm_element_from_mass(
        if i < top.masses.length() {
          top.masses[i]
        } else {
          0.0
        },
      )
    }
    let charge = if i < top.charges.length() { top.charges[i] } else { 0.0 }
    let mass = if i < top.masses.length() { top.masses[i] } else { 0.0 }
    // Residue id from residue pointers (0-based atom index -> 1-based resid).
    let mut rid = 1
    for r = 0; r < top.residue_pointers.length(); r = r + 1 {
      if top.residue_pointers[r] <= i {
        rid = r + 1
      }
    }
    atoms.push({
      name,
      resname,
      resid: rid,
      chain: "A",
      segid: "SYS",
      element,
      atom_type,
      charge,
      mass,
      is_protein: omm_is_protein_resname(resname),
      is_water: resname == "WAT" || resname == "HOH",
      is_backbone: name == "N" || name == "CA" || name == "C" || name == "O",
    })
  }
  let frames : Array[Array[Array[Double]]] = Array(capacity=0)
  if inpcrd.length() > 0 {
    let coords = omm_parse_inpcrd(inpcrd)
    if coords.length() > 0 {
      frames.push(coords)
    }
  }
  Universe::new(atoms, Trajectory::new(frames))
}

///|
/// LJ (sigma, epsilon) for a type pair from ACOEF/BCOEF tables.
/// Type indices are 1-based; the symmetric pair (i <= j) is enforced.
pub fn omm_prmtop_lj(
  top : OmmPrmtop,
  type_i : Int32,
  type_j : Int32,
) -> (Double, Double) {
  let nt = top.n_types
  if nt <= 0 || top.lj_a.length() == 0 || type_i < 1 || type_j < 1 {
    return (0.0, 0.0)
  }
  let mut i = type_i
  let mut j = type_j
  if i > j {
    i = type_j
    j = type_i
  }
  if i > nt || j > nt {
    return (0.0, 0.0)
  }
  let idx = (2 * nt - i) * (i - 1) / 2 + j - 1
  if idx < 0 || idx >= top.lj_a.length() || idx >= top.lj_b.length() {
    return (0.0, 0.0)
  }
  let a = top.lj_a[idx]
  let b = top.lj_b[idx]
  if a <= 0.0 || b <= 0.0 {
    return (0.0, 0.0)
  }
  // A = 4 eps sigma^12, B = 4 eps sigma^6
  let sigma = @math.pow(a / b, 1.0 / 6.0)
  let epsilon = b * b / (4.0 * a)
  (sigma, epsilon)
}

///|
/// Build an OmmSystem (particles + valence terms + periodic box) from a
/// prmtop file. Box comes from IFLBOX-based POINTERS heuristics: use
/// `omm_system_from_prmtop_box` to set it explicitly.
pub fn omm_system_from_prmtop(prmtop_content : String) -> OmmSystem {
  let top = omm_parse_prmtop_full(prmtop_content)
  let system = OmmSystem::new()
  for i = 0; i < top.n_atoms; i = i + 1 {
    let mass = if i < top.masses.length() { top.masses[i] } else { 1.0 }
    let charge = if i < top.charges.length() { top.charges[i] } else { 0.0 }
    let ti = if i < top.type_indices.length() { top.type_indices[i] } else { 1 }
    let (sigma, epsilon) = omm_prmtop_lj(top, ti, ti)
    ignore(
      system.add_particle(
        particle=OmmParticle::new(mass~, charge~, sigma~, epsilon~),
      ),
    )
  }
  for b in top.bonds {
    let k = if b.2 < top.bond_force_constants.length() {
      top.bond_force_constants[b.2]
    } else {
      0.0
    }
    let r0 = if b.2 < top.bond_equil_values.length() {
      top.bond_equil_values[b.2]
    } else {
      0.0
    }
    system.add_bond(b.0, b.1, k, r0)
  }
  for a in top.angles {
    let k = if a.3 < top.angle_force_constants.length() {
      top.angle_force_constants[a.3]
    } else {
      0.0
    }
    // prmtop angle equilibrium is in degrees -> radians.
    let theta0 = if a.3 < top.angle_equil_values.length() {
      top.angle_equil_values[a.3] * 3.141592653589793 / 180.0
    } else {
      0.0
    }
    system.add_angle(a.0, a.1, a.2, k, theta0)
  }
  for d in top.dihedrals {
    let k = if d.4 < top.dihedral_force_constants.length() {
      top.dihedral_force_constants[d.4]
    } else {
      0.0
    }
    let phase_deg = if d.4 < top.dihedral_phases.length() {
      top.dihedral_phases[d.4]
    } else {
      0.0
    }
    let period_f = if d.4 < top.dihedral_periodicities.length() {
      top.dihedral_periodicities[d.4]
    } else {
      1.0
    }
    let n = if period_f < 1.0 { 1 } else { period_f.to_int().abs() }
    system.add_dihedral(
      d.0,
      d.1,
      d.2,
      d.3,
      k,
      phase_deg * 3.141592653589793 / 180.0,
      if n < 1 {
        1
      } else {
        n
      },
    )
  }
  system
}

///|
/// Same as `omm_system_from_prmtop` but also sets the periodic box.
pub fn omm_system_from_prmtop_box(
  prmtop_content : String,
  box~ : Array[Double],
) -> OmmSystem {
  let system = omm_system_from_prmtop(prmtop_content)
  system.set_periodic_box(box[0], box[1], box[2])
  system
}
