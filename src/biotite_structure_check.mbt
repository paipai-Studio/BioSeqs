/// Biotite `structure.check` — steric clash detection and backbone
/// validation.
///
/// StericClashCheck flags atom pairs whose distance falls below a
/// fraction of the sum of their van-der-Waals radii;
/// check_backbone_structure validates the protein backbone geometry
/// (atom presence, adjacent C-N distance, N-CA-C bond angle).

///|
/// van-der-Waals radius (Angstrom, Bondi scale) per element.
pub fn biotite_vdw_radius(element : String) -> Double {
  match element {
    "H" => 1.10
    "D" => 1.10
    "C" => 1.70
    "N" => 1.55
    "O" => 1.52
    "S" => 1.80
    "P" => 1.80
    "F" => 1.47
    "Cl" => 1.75
    "Br" => 1.85
    "I" => 1.98
    "B" => 1.92
    "Si" => 2.10
    "Se" => 1.90
    "Na" => 2.27
    "K" => 2.75
    "Mg" => 1.73
    "Ca" => 2.31
    "Fe" => 2.04
    "Zn" => 2.10
    _ => 1.50
  }
}

///|
/// A single steric clash.
pub struct BiotiteClash {
  atom_i : Int32
  atom_j : Int32
  distance : Double
  mean_distance : Double
}

///|
/// Result of a steric clash check over an AtomArray.
pub struct BiotiteStericClashCheck {
  positions : Array[Vector3]
  elements : Array[String]
  cell_size : Double
  inner : CellList
  respect_stereochemistry : Bool
  bonded : Array[(Int32, Int32)]
}

///|
/// Build a clash check over the given structure.
///
/// With `respect_stereochemistry=true`, atom pairs that are directly
/// bonded (given via `bonds`) are not reported as clashes.
pub fn biotite_steric_clash_check(
  arr : AtomArray,
  bonds~ : BondList? = None,
  respect_stereochemistry~ : Bool = false,
) -> BiotiteStericClashCheck {
  let n = arr.len()
  let mut max_r = 0.0
  for i = 0; i < n; i = i + 1 {
    let r = biotite_vdw_radius(arr.elements[i])
    if r > max_r {
      max_r = r
    }
  }
  let cell_size = 2.0 * max_r
  let flat : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let c = arr.coords[i]
    flat.push([c.x, c.y, c.z])
  }
  let bonded : Array[(Int32, Int32)] = []
  match bonds {
    Some(b) =>
      for bd in b.get_bonds() {
        bonded.push((bd.atom1, bd.atom2))
      }
    None => ()
  }
  BiotiteStericClashCheck::{
    positions: arr.coords,
    elements: arr.elements,
    cell_size,
    inner: celllist_new(flat, cell_size),
    respect_stereochemistry,
    bonded,
  }
}

///|
/// All steric clashes in the structure.
///
/// A clash is reported when distance < (r_i + r_j) * 0.7.
pub fn BiotiteStericClashCheck::get_clashes(
  self : BiotiteStericClashCheck,
) -> Array[BiotiteClash] {
  let n = self.positions.length()
  let cutoff = self.cell_size
  let out : Array[BiotiteClash] = []
  for i = 0; i < n; i = i + 1 {
    let neigh = celllist_neighbors(self.inner, i, cutoff)
    for j in neigh {
      if i >= j {
        continue
      }
      let ri = biotite_vdw_radius(self.elements[i])
      let rj = biotite_vdw_radius(self.elements[j])
      let p1 = self.positions[i]
      let p2 = self.positions[j]
      let d = p1.distance(p2)
      if d < (ri + rj) * 0.7 {
        if self.respect_stereochemistry && self.is_bonded(i, j) {
          continue
        }
        out.push(BiotiteClash::{
          atom_i: i,
          atom_j: j,
          distance: d,
          mean_distance: (ri + rj) / 2.0,
        })
      }
    }
  }
  out
}

///|
/// Indices of all atoms that participate in at least one clash.
pub fn BiotiteStericClashCheck::get_clash_indices(
  self : BiotiteStericClashCheck,
) -> Array[Int32] {
  let clashes = self.get_clashes()
  let seen : Map[Int32, Bool] = Map([], capacity=clashes.length() * 2)
  let out : Array[Int32] = []
  for cl in clashes {
    if seen.get(cl.atom_i).is_empty() {
      seen[cl.atom_i] = true
      out.push(cl.atom_i)
    }
    if seen.get(cl.atom_j).is_empty() {
      seen[cl.atom_j] = true
      out.push(cl.atom_j)
    }
  }
  out.sort_by(fn(a, b) { Int32::compare(a, b) })
  out
}

///|
fn BiotiteStericClashCheck::is_bonded(
  self : BiotiteStericClashCheck,
  i : Int32,
  j : Int32,
) -> Bool {
  for (a, b) in self.bonded {
    if (a == i && b == j) || (a == j && b == i) {
      return true
    }
  }
  false
}

///|
/// A backbone geometry violation.
pub struct BiotiteBackboneViolation {
  kind : String
  atom_index : Int32
  detail : String
}

///|
/// Validate the protein backbone of an AtomArray.
///
/// Checks performed (per protein residue):
/// - each residue contains N, CA and C exactly once (`missing_backbone_atom`)
/// - adjacent residues are linked with 1.0-1.9 A C-N distance
///   (`invalid_peptide_bond`)
/// - the N-CA-C bond angle lies within 95-115 degrees
///   (`invalid_bond_angle`)
pub fn biotite_check_backbone_structure(
  arr : AtomArray,
) -> Array[BiotiteBackboneViolation] {
  let n = arr.len()
  let out : Array[BiotiteBackboneViolation] = []
  // Per-residue backbone atom collection (res_id + chain key).
  let res_keys : Map[String, Int32] = Map([], capacity=16)
  let res_order : Array[String] = []
  let res_n : Map[String, Int32] = Map([], capacity=16)
  let res_ca : Map[String, Int32] = Map([], capacity=16)
  let res_c : Map[String, Int32] = Map([], capacity=16)
  let is_protein : Array[Bool] = Array(capacity=n)
  let protein_res : Map[String, Bool] = Map([], capacity=16)
  for i = 0; i < n; i = i + 1 {
    let protein = biotite_is_protein_residue(arr.res_names[i])
    is_protein.push(protein)
    if !protein {
      continue
    }
    let key = arr.chain_ids[i].to_string() + ":" + arr.res_ids[i].to_string()
    if res_keys.get(key).is_empty() {
      res_keys[key] = res_order.length()
      res_order.push(key)
      res_n[key] = -1
      res_ca[key] = -1
      res_c[key] = -1
      protein_res[key] = true
    }
    let name = arr.names[i]
    if name == "N" {
      res_n[key] = i
    } else if name == "CA" {
      res_ca[key] = i
    } else if name == "C" {
      res_c[key] = i
    }
  }
  // 1) missing backbone atoms
  for key in res_order {
    if res_n.get(key).unwrap_or(-1) < 0 ||
      res_ca.get(key).unwrap_or(-1) < 0 ||
      res_c.get(key).unwrap_or(-1) < 0 {
      let idx = if res_n.get(key).unwrap_or(-1) >= 0 {
        res_n.get(key).unwrap_or(-1)
      } else if res_ca.get(key).unwrap_or(-1) >= 0 {
        res_ca.get(key).unwrap_or(-1)
      } else {
        res_c.get(key).unwrap_or(-1)
      }
      out.push(BiotiteBackboneViolation::{
        kind: "missing_backbone_atom",
        atom_index: idx,
        detail: key,
      })
    }
  }
  // 2) adjacent C-N distance + 3) N-CA-C bond angle
  for k = 1; k < res_order.length(); k = k + 1 {
    let prev = res_order[k - 1]
    let cur = res_order[k]
    let ci = res_c.get(prev).unwrap_or(-1)
    let ni = res_n.get(cur).unwrap_or(-1)
    if ci >= 0 && ni >= 0 {
      let d = arr.coords[ci].distance(arr.coords[ni])
      if d < 1.0 || d > 1.9 {
        out.push(BiotiteBackboneViolation::{
          kind: "invalid_peptide_bond",
          atom_index: ni,
          detail: "C-N distance " + d.to_string(),
        })
      }
    }
    let nai = res_n.get(cur).unwrap_or(-1)
    let cai = res_ca.get(cur).unwrap_or(-1)
    let cci = res_c.get(cur).unwrap_or(-1)
    if nai >= 0 && cai >= 0 && cci >= 0 {
      let angle = biotite_backbone_angle(
        arr.coords[nai],
        arr.coords[cai],
        arr.coords[cci],
      )
      let deg = angle * 180.0 / 3.14159265358979323846
      if deg < 95.0 || deg > 115.0 {
        out.push(BiotiteBackboneViolation::{
          kind: "invalid_bond_angle",
          atom_index: cai,
          detail: "N-CA-C angle " + deg.to_string(),
        })
      }
    }
  }
  out
}

///|
/// Angle (radians) at point b formed by a-b-c.
fn biotite_backbone_angle(a : Vector3, b : Vector3, c : Vector3) -> Double {
  let v1x = a.x - b.x
  let v1y = a.y - b.y
  let v1z = a.z - b.z
  let v2x = c.x - b.x
  let v2y = c.y - b.y
  let v2z = c.z - b.z
  let dot = v1x * v2x + v1y * v2y + v1z * v2z
  let n1 = @math.pow(v1x * v1x + v1y * v1y + v1z * v1z, 0.5)
  let n2 = @math.pow(v2x * v2x + v2y * v2y + v2z * v2z, 0.5)
  if n1 <= 0.0 || n2 <= 0.0 {
    return 0.0
  }
  let cos_val = dot / (n1 * n2)
  let clamped = if cos_val > 1.0 {
    1.0
  } else if cos_val < -1.0 {
    -1.0
  } else {
    cos_val
  }
  @math.acos(clamped)
}

///|
/// Whether a residue name belongs to the standard 20 protein residues.
fn biotite_is_protein_residue(res_name : String) -> Bool {
  match res_name {
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
pub fn _biotite_structure_check_dummy() -> Unit {
  let arr = AtomArray::new()
  let check = biotite_steric_clash_check(arr)
  let _ = check.get_clashes()
  let _ = check.get_clash_indices()
  let _ = biotite_check_backbone_structure(arr)
  let _ = biotite_vdw_radius("C")
  ignore(check.cell_size)
  ignore(check.respect_stereochemistry)
}
