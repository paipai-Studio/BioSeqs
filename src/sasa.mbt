///|
/// Bio.PDB.SASA module, inspired by Biopython's Bio.PDB.SASA.
///
/// Calculates Solvent Accessible Surface Area (SASA) for protein
/// structures using the Shrake-Rupley algorithm (rolling ball method):
///
/// 1. For each atom, generate N test points on a sphere around the atom
///    (radius = van der Waals radius + probe radius).
/// 2. For each test point, check whether it is buried by any neighboring
///    atom (within that atom's van der Waals radius + probe radius).
/// 3. SASA = (accessible points / total points) * 4 * pi * r^2.
///
/// NOTE: The names `calc_sasa`, `calc_total_sasa` and `calc_residue_sasa`
//        already exist in `pdb_analysis.mbt` (operating on a single
//        `Chain` and returning tuples). To avoid duplicate-definition
//        errors (MoonBit has no function overloading) the structured
//        `Structure`-level entry points here use a `sasa_` prefix:
//        `sasa_calc`, `sasa_calc_residue`, `sasa_calc_total`.

///|
/// Per-atom SASA result.
pub struct SasaResult {
  atom_index : Int
  atom_name : String
  resname : String
  chainid : Char
  resseq : Int
  sasa : Double
}

///|
/// Construct a `SasaResult`.
pub fn SasaResult::new(
  atom_index~ : Int,
  atom_name~ : String,
  resname~ : String,
  chainid~ : Char,
  resseq~ : Int,
  sasa~ : Double,
) -> SasaResult {
  SasaResult::{ atom_index, atom_name, resname, chainid, resseq, sasa }
}

///|
pub fn SasaResult::get_atom_index(self : SasaResult) -> Int {
  self.atom_index
}

///|
pub fn SasaResult::get_atom_name(self : SasaResult) -> String {
  self.atom_name
}

///|
pub fn SasaResult::get_resname(self : SasaResult) -> String {
  self.resname
}

///|
pub fn SasaResult::get_chainid(self : SasaResult) -> Char {
  self.chainid
}

///|
pub fn SasaResult::get_resseq(self : SasaResult) -> Int {
  self.resseq
}

///|
pub fn SasaResult::get_sasa(self : SasaResult) -> Double {
  self.sasa
}

///|
/// Per-residue SASA aggregation, split into backbone and sidechain.
pub struct ResidueSasa {
  resname : String
  chainid : Char
  resseq : Int
  total_sasa : Double
  backbone_sasa : Double
  sidechain_sasa : Double
  num_atoms : Int
}

///|
/// Construct a `ResidueSasa`.
pub fn ResidueSasa::new(
  resname~ : String,
  chainid~ : Char,
  resseq~ : Int,
  total_sasa~ : Double,
  backbone_sasa~ : Double,
  sidechain_sasa~ : Double,
  num_atoms~ : Int,
) -> ResidueSasa {
  ResidueSasa::{
    resname,
    chainid,
    resseq,
    total_sasa,
    backbone_sasa,
    sidechain_sasa,
    num_atoms,
  }
}

///|
pub fn ResidueSasa::get_resname(self : ResidueSasa) -> String {
  self.resname
}

///|
pub fn ResidueSasa::get_chainid(self : ResidueSasa) -> Char {
  self.chainid
}

///|
pub fn ResidueSasa::get_resseq(self : ResidueSasa) -> Int {
  self.resseq
}

///|
pub fn ResidueSasa::get_total_sasa(self : ResidueSasa) -> Double {
  self.total_sasa
}

///|
pub fn ResidueSasa::get_backbone_sasa(self : ResidueSasa) -> Double {
  self.backbone_sasa
}

///|
pub fn ResidueSasa::get_sidechain_sasa(self : ResidueSasa) -> Double {
  self.sidechain_sasa
}

///|
pub fn ResidueSasa::get_num_atoms(self : ResidueSasa) -> Int {
  self.num_atoms
}

///|
/// Overall SASA summary for a structure.
pub struct SasaSummary {
  total_sasa : Double
  backbone_sasa : Double
  sidechain_sasa : Double
  num_atoms : Int
  num_residues : Int
}

///|
/// Construct a `SasaSummary`.
pub fn SasaSummary::new(
  total_sasa~ : Double,
  backbone_sasa~ : Double,
  sidechain_sasa~ : Double,
  num_atoms~ : Int,
  num_residues~ : Int,
) -> SasaSummary {
  SasaSummary::{
    total_sasa,
    backbone_sasa,
    sidechain_sasa,
    num_atoms,
    num_residues,
  }
}

///|
pub fn SasaSummary::get_total_sasa(self : SasaSummary) -> Double {
  self.total_sasa
}

///|
pub fn SasaSummary::get_backbone_sasa(self : SasaSummary) -> Double {
  self.backbone_sasa
}

///|
pub fn SasaSummary::get_sidechain_sasa(self : SasaSummary) -> Double {
  self.sidechain_sasa
}

///|
pub fn SasaSummary::get_num_atoms(self : SasaSummary) -> Int {
  self.num_atoms
}

///|
pub fn SasaSummary::get_num_residues(self : SasaSummary) -> Int {
  self.num_residues
}

///|
/// Standard van der Waals radii (Angstroms) for common elements.
/// Lookup is case-insensitive. Unknown elements default to 1.70 (carbon).
pub fn sasa_vdw_radius(element : String) -> Double {
  let e = element.to_upper()
  if e == "H" {
    1.20
  } else if e == "C" {
    1.70
  } else if e == "N" {
    1.55
  } else if e == "O" {
    1.52
  } else if e == "S" {
    1.80
  } else if e == "P" {
    1.80
  } else if e == "F" {
    1.47
  } else if e == "CL" {
    1.75
  } else if e == "BR" {
    1.85
  } else if e == "I" {
    1.98
  } else {
    1.70
  }
}

///|
/// Resolve the element string for an atom. Uses the `element` field when
/// present, otherwise derives the element from the first character of the
/// atom name (the common convention for protein atoms such as N, CA, C, O).
fn sasa_resolve_element(atom : Atom) -> String {
  let e = atom.element
  if e != "" {
    e.to_upper()
  } else {
    let name = atom.name
    if name.length() == 0 {
      "C"
    } else {
      let first = substring(name, 0, 1)
      if first == "" {
        "C"
      } else {
        first.to_upper()
      }
    }
  }
}

///|
/// Identify whether an atom name is a backbone atom.
/// Backbone atoms: N, CA, C, O, OXT.
fn sasa_is_backbone(name : String) -> Bool {
  let n = name.to_upper()
  n == "N" || n == "CA" || n == "C" || n == "O" || n == "OXT"
}

///|
/// Generate approximately uniform points on a unit sphere using the
/// Fibonacci spiral (golden angle) algorithm. Returns unit vectors.
/// `n` is clamped to at least 1.
pub fn sasa_generate_sphere_points(n : Int) -> Array[Vector3] {
  let points : Array[Vector3] = Array::new()
  let count = if n < 1 { 1 } else { n }
  if count == 1 {
    points.push(Vector3::new(0.0, 0.0, 1.0))
    return points
  }
  let golden_angle = 3.14159265358979323846 * (3.0 - 5.0.sqrt())
  let denom = (count - 1).to_double()
  for i in 0..<count {
    let y = 1.0 - (2.0 * i.to_double()) / denom
    let r = (1.0 - y * y).sqrt()
    let theta = golden_angle * i.to_double()
    let x = @math.cos(theta) * r
    let z = @math.sin(theta) * r
    points.push(Vector3::new(x, y, z))
  }
  points
}

///|
/// Calculate per-atom SASA for a structure using the Shrake-Rupley
/// rolling-ball algorithm.
///
/// `n_points`: number of test points per atom sphere.
/// `probe_radius`: solvent probe radius in Angstroms (1.4 for water).
pub fn sasa_calc(
  structure : Structure,
  n_points : Int,
  probe_radius : Double,
) -> Array[SasaResult] {
  let results : Array[SasaResult] = Array::new()
  let atoms = structure.get_atoms()
  let n_atoms = atoms.length()
  if n_atoms == 0 {
    return results
  }

  // Collect coordinates and per-atom probe-sphere radii.
  let coords : Array[Vector3] = Array::new()
  let radii : Array[Double] = Array::new()
  for i in 0..<n_atoms {
    let a = atoms[i]
    coords.push(a.coord)
    radii.push(sasa_vdw_radius(sasa_resolve_element(a)) + probe_radius)
  }

  // With non-positive point counts, every atom contributes zero SASA.
  if n_points <= 0 {
    for i in 0..<n_atoms {
      let a = atoms[i]
      results.push(
        SasaResult::new(
          atom_index=i,
          atom_name=a.name,
          resname=a.resname,
          chainid=a.chainid,
          resseq=a.resseq,
          sasa=0.0,
        ),
      )
    }
    return results
  }

  // Precompute neighbors for pruning: a point on atom i's sphere can be
  // buried by atom j only if dist(i, j) < radii[i] + radii[j].
  let neighbors : Array[Array[Int]] = Array::new()
  for i in 0..<n_atoms {
    let nb : Array[Int] = Array::new()
    let ci = coords[i]
    let ri = radii[i]
    for j in 0..<n_atoms {
      if i != j {
        let d = ci.distance(coords[j])
        if d < ri + radii[j] {
          nb.push(j)
        }
      }
    }
    neighbors.push(nb)
  }

  let sphere = sasa_generate_sphere_points(n_points)
  let n_sphere = sphere.length()
  let pi = 3.14159265358979323846

  for i in 0..<n_atoms {
    let ci = coords[i]
    let ri = radii[i]
    let r_sq = ri * ri
    let nb = neighbors[i]
    let mut n_accessible = 0
    for p in sphere {
      let px = ci.x + p.x * ri
      let py = ci.y + p.y * ri
      let pz = ci.z + p.z * ri
      let mut buried = false
      for j in nb {
        let cj = coords[j]
        let dx = px - cj.x
        let dy = py - cj.y
        let dz = pz - cj.z
        let rj = radii[j]
        if dx * dx + dy * dy + dz * dz < rj * rj {
          buried = true
          break
        }
      }
      if !buried {
        n_accessible = n_accessible + 1
      }
    }
    let sphere_area = 4.0 * pi * r_sq
    let sasa = (n_accessible.to_double() / n_sphere.to_double()) * sphere_area
    let a = atoms[i]
    results.push(
      SasaResult::new(
        atom_index=i,
        atom_name=a.name,
        resname=a.resname,
        chainid=a.chainid,
        resseq=a.resseq,
        sasa~,
      ),
    )
  }

  results
}

///|
/// Aggregate per-atom SASA results into per-residue SASA.
///
/// Atoms are matched to residues by index, assuming `atom_results` was
/// produced from `structure` (atoms are enumerated in model -> chain ->
/// residue -> atom order by `Structure::get_atoms`).
pub fn sasa_calc_residue(
  atom_results : Array[SasaResult],
  structure : Structure,
) -> Array[ResidueSasa] {
  let results : Array[ResidueSasa] = Array::new()

  // Map atom index -> sasa value for robust lookup.
  let sasa_map : Map[Int, Double] = Map([], capacity=atom_results.length())
  for r in atom_results {
    sasa_map[r.atom_index] = r.sasa
  }

  let mut idx = 0
  for residue in structure.get_residues() {
    let mut total = 0.0
    let mut backbone = 0.0
    let mut count = 0
    for atom in residue.atoms {
      let sasa = sasa_map.get(idx).unwrap_or(0.0)
      total = total + sasa
      if sasa_is_backbone(atom.name) {
        backbone = backbone + sasa
      }
      count = count + 1
      idx = idx + 1
    }
    let sidechain = total - backbone
    results.push(
      ResidueSasa::new(
        resname=residue.resname,
        chainid=residue.chainid,
        resseq=residue.resseq,
        total_sasa=total,
        backbone_sasa=backbone,
        sidechain_sasa=sidechain,
        num_atoms=count,
      ),
    )
  }

  results
}

///|
/// Summarize total SASA across all per-atom results.
///
/// Backbone/sidechain split is derived from atom names; the residue count
/// is the number of distinct (chainid, resseq) pairs.
pub fn sasa_calc_total(atom_results : Array[SasaResult]) -> SasaSummary {
  let mut total = 0.0
  let mut backbone = 0.0
  let mut count = 0
  let residue_keys : Map[String, Bool] = Map([], capacity=16)
  for r in atom_results {
    total = total + r.sasa
    if sasa_is_backbone(r.atom_name) {
      backbone = backbone + r.sasa
    }
    count = count + 1
    let key = r.chainid.to_string() + ":" + r.resseq.to_string()
    residue_keys[key] = true
  }
  let sidechain = total - backbone
  let num_residues = residue_keys.keys().collect().length()
  SasaSummary::new(
    total_sasa=total,
    backbone_sasa=backbone,
    sidechain_sasa=sidechain,
    num_atoms=count,
    num_residues~,
  )
}

///|
/// Build a small sample peptide structure (ALA-GLY-SER) for testing and
/// demonstration. All atoms carry explicit element labels.
pub fn sasa_sample_structure() -> Structure {
  let atoms1 : Array[Atom] = [
    Atom::new(
      name="N",
      coord=Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="N",
    ),
    Atom::new(
      name="CA",
      coord=Vector3::new(1.46, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
    Atom::new(
      name="C",
      coord=Vector3::new(2.01, 1.42, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
    Atom::new(
      name="O",
      coord=Vector3::new(1.23, 2.36, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="O",
    ),
    Atom::new(
      name="CB",
      coord=Vector3::new(1.05, -0.79, 1.21),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
  ]

  let atoms2 : Array[Atom] = [
    Atom::new(
      name="N",
      coord=Vector3::new(3.25, 1.56, 0.0),
      resname="GLY",
      chainid='A',
      resseq=2,
      element="N",
    ),
    Atom::new(
      name="CA",
      coord=Vector3::new(3.93, 2.86, 0.0),
      resname="GLY",
      chainid='A',
      resseq=2,
      element="C",
    ),
    Atom::new(
      name="C",
      coord=Vector3::new(5.43, 2.84, 0.0),
      resname="GLY",
      chainid='A',
      resseq=2,
      element="C",
    ),
    Atom::new(
      name="O",
      coord=Vector3::new(6.06, 3.88, 0.0),
      resname="GLY",
      chainid='A',
      resseq=2,
      element="O",
    ),
  ]

  let atoms3 : Array[Atom] = [
    Atom::new(
      name="N",
      coord=Vector3::new(6.13, 1.74, 0.0),
      resname="SER",
      chainid='A',
      resseq=3,
      element="N",
    ),
    Atom::new(
      name="CA",
      coord=Vector3::new(7.54, 1.42, 0.0),
      resname="SER",
      chainid='A',
      resseq=3,
      element="C",
    ),
    Atom::new(
      name="C",
      coord=Vector3::new(8.18, 2.07, 1.21),
      resname="SER",
      chainid='A',
      resseq=3,
      element="C",
    ),
    Atom::new(
      name="O",
      coord=Vector3::new(7.64, 2.18, 2.30),
      resname="SER",
      chainid='A',
      resseq=3,
      element="O",
    ),
    Atom::new(
      name="CB",
      coord=Vector3::new(7.66, -0.07, 0.0),
      resname="SER",
      chainid='A',
      resseq=3,
      element="C",
    ),
    Atom::new(
      name="OG",
      coord=Vector3::new(8.99, -0.46, 0.0),
      resname="SER",
      chainid='A',
      resseq=3,
      element="O",
    ),
  ]

  let res1 = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=atoms1)
  let res2 = Residue::new(resname="GLY", chainid='A', resseq=2, atoms=atoms2)
  let res3 = Residue::new(resname="SER", chainid='A', resseq=3, atoms=atoms3)
  let chain = Chain::new(id='A', residues=[res1, res2, res3])
  let model = Model::new(id=1, chains=[chain])
  Structure::new(id="sample", models=[model])
}
