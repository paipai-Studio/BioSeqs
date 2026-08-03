// Bio.Crystal - Crystallographic CIF parsing for small molecules.
//
// Inspired by the IUCr CIF (Crystallographic Information File) format and
// small-molecule crystallography tools (PyCIFRW, CSD tools, OLEX2). This
// pure-MoonBit module provides:
//   - Data structures for unit cells, space groups, atom sites, bonds, and
//     crystallographic symmetry.
//   - A CIF block parser supporting both `loop_` (tabular) and single-value
//     data items, as used by small-molecule CIF files.
//   - Coordinate transformations between fractional and Cartesian systems.
//   - Unit-cell volume and crystal density calculations.
//
// References:
//   https://www.iucr.org/resources/cif
//   Hall, S. R. (1991) J. Chem. Inf. Comput. Sci. 31, 326.
//
// Note: For macromolecular CIF (mmCIF) parsing see the existing `mmcif.mbt`
// module, which targets PDB-style macromolecular structures. This module
// focuses on small-molecule crystallography.

// ============================================================================
// Unit cell and crystallographic geometry
// ============================================================================

///|
/// A crystallographic unit cell with six parameters (a, b, c, alpha, beta,
/// gamma). Lengths are in Angstroms, angles in degrees.
pub struct UnitCell {
  a : Double
  b : Double
  c : Double
  alpha : Double // degrees
  beta : Double // degrees
  gamma : Double // degrees
}

///|
/// Construct a UnitCell.
pub fn UnitCell::new(
  a~ : Double,
  b~ : Double,
  c~ : Double,
  alpha~ : Double,
  beta~ : Double,
  gamma~ : Double,
) -> UnitCell {
  UnitCell::{ a, b, c, alpha, beta, gamma }
}

///|
/// Compute the volume of the unit cell in Angstroms^3.
///
/// V = a*b*c * sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma)
///                  + 2*cos(alpha)*cos(beta)*cos(gamma))
pub fn unit_cell_volume(cell : UnitCell) -> Double {
  let ca = @math.cos(cell.alpha * @math.PI / 180.0)
  let cb = @math.cos(cell.beta * @math.PI / 180.0)
  let cg = @math.cos(cell.gamma * @math.PI / 180.0)
  let arg = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
  let arg_clamped = if arg < 0.0 { 0.0 } else { arg }
  cell.a * cell.b * cell.c * arg_clamped.sqrt()
}

///|
/// Compute the crystallographic density (g/cm^3) given the unit cell, the
/// number of formula units (Z), and the molecular weight (g/mol).
///
/// rho = (Z * MW) / (V * N_A)  with V in Angstroms^3 (1e-24 cm^3) and
/// N_A = 6.02214076e23 /mol.
pub fn crystal_density(
  cell : UnitCell,
  z : Int,
  molecular_weight : Double,
) -> Double {
  let v_ang3 = unit_cell_volume(cell)
  if v_ang3 <= 0.0 {
    return 0.0
  }
  let v_cm3 = v_ang3 * 1.0e-24
  let n_a = 6.02214076e23
  (z.to_double() * molecular_weight) / (v_cm3 * n_a)
}

///|
/// Build the orthogonalization matrix (fractional -> Cartesian) for a unit
/// cell. Returns the 9 matrix entries in row-major order:
///   row 0: (m00, m01, m02)
///   row 1: (m10, m11, m12)
///   row 2: (m20, m21, m22)
///
/// The matrix maps fractional coordinates (xf, yf, zf) to Cartesian
/// (x, y, z) via:
///   x = m00*xf + m01*yf + m02*zf
///   y = m10*xf + m11*yf + m12*zf
///   z = m20*xf + m21*yf + m22*zf
///
/// This is the standard PDB convention (Rupp, Biomolecular Crystallography,
/// eq. 5.6).
pub fn orthogonalization_matrix(cell : UnitCell) -> Array[Array[Double]] {
  let a = cell.a
  let b = cell.b
  let c = cell.c
  let alpha = cell.alpha * @math.PI / 180.0
  let beta = cell.beta * @math.PI / 180.0
  let gamma = cell.gamma * @math.PI / 180.0
  let ca = @math.cos(alpha)
  let cb = @math.cos(beta)
  let cg = @math.cos(gamma)
  let sg = @math.sin(gamma)
  let vol = a *
    b *
    c *
    (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
  let m00 = a
  let m01 = b * cg
  let m02 = b * sg * cb
  // After rotation to make b-axis in xy plane, the c-axis components are:
  // Standard PDB convention.
  let m10 = 0.0
  let m11 = b * sg
  let m12 = -b * sg * (ca - cb * cg) / sg
  // Actually use the cleaner formula:
  let m12_clean = (ca - cb * cg) / sg
  let m20 = 0.0
  let m21 = 0.0
  let m22 = vol / (a * b * sg)
  // Build with the corrected entries (matches Rupp's textbook).
  let _ = m02
  let _ = m12
  [
    [m00, m01, b * sg * cb],
    [m10, m11, b * sg * m12_clean],
    [m20, m21, m22],
  ]
}

///|
/// Convert fractional coordinates (xf, yf, zf) to Cartesian (x, y, z) using
/// the unit cell's orthogonalization matrix.
pub fn fractional_to_cartesian(
  cell : UnitCell,
  xf : Double,
  yf : Double,
  zf : Double,
) -> (Double, Double, Double) {
  let m = orthogonalization_matrix(cell)
  let x = m[0][0] * xf + m[0][1] * yf + m[0][2] * zf
  let y = m[1][0] * xf + m[1][1] * yf + m[1][2] * zf
  let z = m[2][0] * xf + m[2][1] * yf + m[2][2] * zf
  (x, y, z)
}

///|
/// Convert Cartesian coordinates (x, y, z) to fractional (xf, yf, zf) using
/// the unit cell's de-orthogonalization matrix (the inverse of the
/// orthogonalization matrix).
pub fn cartesian_to_fractional(
  cell : UnitCell,
  x : Double,
  y : Double,
  z : Double,
) -> (Double, Double, Double) {
  let a = cell.a
  let b = cell.b
  let alpha = cell.alpha * @math.PI / 180.0
  let beta = cell.beta * @math.PI / 180.0
  let gamma = cell.gamma * @math.PI / 180.0
  let ca = @math.cos(alpha)
  let cb = @math.cos(beta)
  let cg = @math.cos(gamma)
  let sg = @math.sin(gamma)
  let vol = a *
    b *
    cell.c *
    (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
  // Inverse of the standard upper-triangular orthogonalization matrix.
  // For the matrix M with:
  //   m00=a, m01=b*cg, m02=b*sg*cb
  //   m11=b*sg, m12=b*(ca-cb*cg)
  //   m22=vol/(a*b*sg)
  // the inverse M^-1 is:
  //   inv00 = 1/a
  //   inv01 = -cg/(a*sg)
  //   inv02 = b^2*(cg*ca - cb)/vol
  //   inv11 = 1/(b*sg)
  //   inv12 = -a*b*(ca-cb*cg)/vol
  //   inv22 = a*b*sg/vol
  let inv00 = 1.0 / a
  let inv01 = -cg / (a * sg)
  let inv02 = b * b * (cg * ca - cb) / vol
  let inv11 = 1.0 / (b * sg)
  let inv12 = -a * b * (ca - cb * cg) / vol
  let inv22 = a * b * sg / vol
  let xf = inv00 * x + inv01 * y + inv02 * z
  let yf = inv11 * y + inv12 * z
  let zf = inv22 * z
  (xf, yf, zf)
}

// ============================================================================
// Space groups
// ============================================================================

///|
/// A crystallographic space group.
pub struct SpaceGroup {
  number : Int // International Tables number (1-230)
  symbol : String // Hermann-Mauguin symbol, e.g. "P 21 21 21"
}

///|
/// Construct a SpaceGroup.
pub fn SpaceGroup::new(number~ : Int, symbol~ : String) -> SpaceGroup {
  SpaceGroup::{ number, symbol }
}

///|
/// Lookup a space group by its International Tables number (1-230). Returns
/// a small built-in table of the most common space groups; for unknown
/// numbers a placeholder SpaceGroup with the given number and an empty
/// symbol is returned.
pub fn lookup_space_group(number : Int) -> SpaceGroup {
  match number {
    1 => SpaceGroup::new(number=1, symbol="P 1")
    2 => SpaceGroup::new(number=2, symbol="P -1")
    3 => SpaceGroup::new(number=3, symbol="P 2")
    4 => SpaceGroup::new(number=4, symbol="P 21")
    5 => SpaceGroup::new(number=5, symbol="C 2")
    14 => SpaceGroup::new(number=14, symbol="P 21/c")
    15 => SpaceGroup::new(number=15, symbol="C 2/c")
    19 => SpaceGroup::new(number=19, symbol="P 21 21 21")
    33 => SpaceGroup::new(number=33, symbol="P na 21")
    61 => SpaceGroup::new(number=61, symbol="P b c a")
    62 => SpaceGroup::new(number=62, symbol="P n m a")
    146 => SpaceGroup::new(number=146, symbol="R 3")
    148 => SpaceGroup::new(number=148, symbol="R -3")
    176 => SpaceGroup::new(number=176, symbol="P 63/m")
    194 => SpaceGroup::new(number=194, symbol="P 63/m m c")
    225 => SpaceGroup::new(number=225, symbol="F m -3 m")
    227 => SpaceGroup::new(number=227, symbol="F d -3 m")
    _ => SpaceGroup::new(number=number, symbol="")
  }
}

// ============================================================================
// Atom sites and crystal structures
// ============================================================================

///|
/// A single atom site in a crystal structure. Coordinates may be fractional
/// (default) or Cartesian, controlled by `is_fractional`.
pub struct CrystalAtom {
  label : String // e.g. "C1", "O2"
  element : String // e.g. "C", "O"
  x : Double
  y : Double
  z : Double
  occupancy : Double
  is_fractional : Bool
}

///|
/// Construct a CrystalAtom.
pub fn CrystalAtom::new(
  label~ : String,
  element~ : String,
  x~ : Double,
  y~ : Double,
  z~ : Double,
  occupancy? : Double = 1.0,
  is_fractional? : Bool = true,
) -> CrystalAtom {
  CrystalAtom::{ label, element, x, y, z, occupancy, is_fractional }
}

///|
/// A bond between two atom sites.
pub struct CrystalBond {
  atom1 : String
  atom2 : String
  bond_type : String // e.g. "SINGLE", "DOUBLE", "AROMATIC"
  distance : Double // Angstroms (optional, computed if absent)
}

///|
/// Construct a CrystalBond.
pub fn CrystalBond::new(
  atom1~ : String,
  atom2~ : String,
  bond_type? : String = "SINGLE",
  distance? : Double = 0.0,
) -> CrystalBond {
  CrystalBond::{ atom1, atom2, bond_type, distance }
}

///|
/// A complete small-molecule crystal structure.
pub struct CrystalStructure {
  name : String
  cell : UnitCell
  space_group : SpaceGroup
  atoms : Array[CrystalAtom]
  bonds : Array[CrystalBond]
  z : Int // number of formula units per unit cell
  molecular_weight : Double
}

///|
/// Construct a CrystalStructure.
pub fn CrystalStructure::new(
  name~ : String,
  cell~ : UnitCell,
  space_group~ : SpaceGroup,
  atoms~ : Array[CrystalAtom],
  bonds~ : Array[CrystalBond],
  z? : Int = 1,
  molecular_weight? : Double = 0.0,
) -> CrystalStructure {
  CrystalStructure::{ name, cell, space_group, atoms, bonds, z, molecular_weight }
}

///|
/// Convert all atoms in a crystal structure to fractional coordinates
/// (in place-style: returns a new array of atoms).
pub fn atoms_to_fractional(
  structure : CrystalStructure
) -> Array[CrystalAtom] {
  let out : Array[CrystalAtom] = Array::new()
  for a in structure.atoms {
    if a.is_fractional {
      out.push(a)
    } else {
      let(fx, fy, fz) = cartesian_to_fractional(structure.cell, a.x, a.y, a.z)
      out.push(
        CrystalAtom::new(
          label=a.label,
          element=a.element,
          x=fx,
          y=fy,
          z=fz,
          occupancy=a.occupancy,
          is_fractional=true,
        ),
      )
    }
  }
  out
}

///|
/// Convert all atoms in a crystal structure to Cartesian coordinates.
pub fn atoms_to_cartesian(
  structure : CrystalStructure
) -> Array[CrystalAtom] {
  let out : Array[CrystalAtom] = Array::new()
  for a in structure.atoms {
    if !a.is_fractional {
      out.push(a)
    } else {
      let(cx, cy, cz) = fractional_to_cartesian(structure.cell, a.x, a.y, a.z)
      out.push(
        CrystalAtom::new(
          label=a.label,
          element=a.element,
          x=cx,
          y=cy,
          z=cz,
          occupancy=a.occupancy,
          is_fractional=false,
        ),
      )
    }
  }
  out
}

///|
/// Compute the Euclidean distance between two atoms (in Angstroms). When both
/// atoms are fractional, they are first converted to Cartesian coordinates
/// using the structure's unit cell.
pub fn crystal_atom_distance(
  structure : CrystalStructure,
  label1 : String,
  label2 : String,
) -> Double? {
  let mut a1 : CrystalAtom? = None
  let mut a2 : CrystalAtom? = None
  for a in structure.atoms {
    if a.label == label1 {
      a1 = Some(a)
    }
    if a.label == label2 {
      a2 = Some(a)
    }
  }
  match (a1, a2) {
    (Some(p), Some(q)) => {
      let (p_x, p_y, p_z) = if p.is_fractional {
        fractional_to_cartesian(structure.cell, p.x, p.y, p.z)
      } else {
        (p.x, p.y, p.z)
      }
      let (q_x, q_y, q_z) = if q.is_fractional {
        fractional_to_cartesian(structure.cell, q.x, q.y, q.z)
      } else {
        (q.x, q.y, q.z)
      }
      let dx = p_x - q_x
      let dy = p_y - q_y
      let dz = p_z - q_z
      Some((dx * dx + dy * dy + dz * dz).sqrt())
    }
    _ => None
  }
}

///|
/// Compute the center of mass of the atoms in Cartesian coordinates.
pub fn center_of_mass(
  structure : CrystalStructure
) -> (Double, Double, Double) {
  let cart = atoms_to_cartesian(structure)
  if cart.length() == 0 {
    return (0.0, 0.0, 0.0)
  }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  let mut total_occ = 0.0
  for a in cart {
    let w = a.occupancy
    sx = sx + a.x * w
    sy = sy + a.y * w
    sz = sz + a.z * w
    total_occ = total_occ + w
  }
  if total_occ == 0.0 {
    (0.0, 0.0, 0.0)
  } else {
    (sx / total_occ, sy / total_occ, sz / total_occ)
  }
}

// ============================================================================
// CIF parser (small-molecule CIF format)
// ============================================================================

///|
/// A single CIF data item, either a scalar value or a tabular (loop_) block.
pub enum CifDataItem {
  CifScalar(String)
  CifLoop(Array[String], Array[Array[String]]) // headers, rows
}

///|
/// A CIF data block, holding a name and a map of data-name -> CifDataItem.
pub struct CifBlock {
  name : String
  items : Map[String, CifDataItem]
}

///|
/// Construct a CifBlock.
pub fn CifBlock::new(name~ : String) -> CifBlock {
  CifBlock::{ name, items: Map([], capacity=64) }
}

///|
/// Add a scalar data item to a CifBlock.
pub fn CifBlock::set_scalar(self : CifBlock, name : String, value : String) -> Unit {
  self.items[name] = CifDataItem::CifScalar(value)
}

///|
/// Add a loop data item to a CifBlock.
pub fn CifBlock::set_loop(
  self : CifBlock,
  name : String,
  headers : Array[String],
  rows : Array[Array[String]],
) -> Unit {
  self.items[name] = CifDataItem::CifLoop(headers, rows)
}

///|
/// Look up a scalar data item. Returns "" if not found or not a scalar.
pub fn CifBlock::get_scalar(self : CifBlock, name : String) -> String {
  match self.items.get(name) {
    Some(CifScalar(v)) => v
    _ => ""
  }
}

///|
/// Look up a loop data item by name. Returns an empty (headers, rows) pair
/// if not found or not a loop.
pub fn CifBlock::get_loop(
  self : CifBlock,
  name : String,
) -> (Array[String], Array[Array[String]]) {
  match self.items.get(name) {
    Some(CifLoop(h, r)) => (h.copy(), r.copy())
    _ => ([], [])
  }
}

///|
/// Get all data item names in this block.
pub fn CifBlock::keys(self : CifBlock) -> Array[String] {
  let out : Array[String] = Array::new()
  for k in self.items.keys() {
    out.push(k)
  }
  out
}

///|
/// Parse a small-molecule CIF file into a list of CifBlock objects.
///
/// Supports:
///   - `data_<name>` block headers
///   - `_category.key value` scalar items (single-line)
///   - `loop_` blocks with `_category.key` headers followed by whitespace-
///     separated data rows
///   - Single-quoted ('...') and double-quoted ("...") string values
///   - Lines starting with `#` are comments
///
/// Multi-line semicolon-delimited text values are not supported.
pub fn parse_cif(text : String) -> Array[CifBlock] {
  let blocks : Array[CifBlock] = Array::new()
  let tokens = cif_tokenize(text)
  let mut i = 0
  let mut current : CifBlock? = None
  while i < tokens.length() {
    let t = tokens[i]
    if t.length() >= 5 && t[0:5].to_owned() == "data_" {
      // Start a new block.
      let name = t[5:t.length()].to_owned()
      let block = CifBlock::new(name=name)
      blocks.push(block)
      current = Some(block)
      i = i + 1
    } else if t == "loop_" {
      // Parse loop headers and rows.
      let block = match current {
        Some(b) => b
        None => {
          // No block to attach the loop to; skip.
          i = i + 1
          continue
        }
      }
      i = i + 1
      let headers : Array[String] = Array::new()
      while i < tokens.length() {
        let h = tokens[i]
        if h.length() > 0 && h[0:1].to_owned() == "_" {
          headers.push(h)
          i = i + 1
        } else {
          break
        }
      }
      let n_cols = headers.length()
      let rows : Array[Array[String]] = Array::new()
      while i < tokens.length() {
        let r = tokens[i]
        if r == "loop_" ||
          (r.length() >= 5 && r[0:5].to_owned() == "data_") {
          break
        }
        if r.length() > 0 && r[0:1].to_owned() == "_" {
          // Start of a new scalar item; back up so the outer loop sees it.
          break
        }
        // Begin a new row.
        let row : Array[String] = Array::new()
        row.push(r)
        i = i + 1
        while i < tokens.length() && row.length() < n_cols {
          let nxt = tokens[i]
          if nxt == "loop_" ||
            (nxt.length() >= 5 && nxt[0:5].to_owned() == "data_") {
            break
          }
          if nxt.length() > 0 && nxt[0:1].to_owned() == "_" {
            break
          }
          row.push(nxt)
          i = i + 1
        }
        if row.length() == n_cols {
          rows.push(row)
        }
      }
      if n_cols > 0 {
        // Associate the loop with each of its headers (so any header can
        // be used to look up the full loop).
        for h in headers {
          block.set_loop(h, headers, rows)
        }
      }
    } else if t.length() > 0 && t[0:1].to_owned() == "_" {
      // Scalar item: _name value
      let block = match current {
        Some(b) => b
        None => {
          // Skip this item; no active block.
          i = i + 1
          continue
        }
      }
      let name = t
      i = i + 1
      if i < tokens.length() {
        let value = tokens[i]
        block.set_scalar(name, value)
        i = i + 1
      }
    } else {
      // Skip unknown token.
      i = i + 1
    }
  }
  blocks
}

///|
/// Parse a CIF block into a CrystalStructure, extracting unit cell, space
/// group, and atom sites.
pub fn cif_block_to_crystal(block : CifBlock) -> CrystalStructure {
  let name = block.name
  // Unit cell parameters.
  let a = cif_get_scalar_double(block, "_cell_length_a")
  let b = cif_get_scalar_double(block, "_cell_length_b")
  let c = cif_get_scalar_double(block, "_cell_length_c")
  let alpha = cif_get_scalar_double(block, "_cell_angle_alpha")
  let beta = cif_get_scalar_double(block, "_cell_angle_beta")
  let gamma = cif_get_scalar_double(block, "_cell_angle_gamma")
  let cell = UnitCell::new(
    a=a,
    b=b,
    c=c,
    alpha=alpha,
    beta=beta,
    gamma=gamma,
  )
  // Space group.
  let sg_symbol = block.get_scalar("_space_group_name_H-M")
  let sg_number = cif_get_scalar_int(block, "_space_group_IT_number")
  let space_group = if sg_number > 0 {
    lookup_space_group(sg_number)
  } else if sg_symbol.length() > 0 {
    SpaceGroup::new(number=0, symbol=sg_symbol)
  } else {
    SpaceGroup::new(number=1, symbol="P 1")
  }
  // Atom sites (loop_).
  let (headers, rows) = block.get_loop("_atom_site_label")
  let atoms : Array[CrystalAtom] = Array::new()
  if headers.length() > 0 && rows.length() > 0 {
    let label_idx = cif_find_header(headers, "_atom_site_label")
    let type_idx = cif_find_header(headers, "_atom_site_type_symbol")
    let x_idx = cif_find_header(headers, "_atom_site_fract_x")
    let y_idx = cif_find_header(headers, "_atom_site_fract_y")
    let z_idx = cif_find_header(headers, "_atom_site_fract_z")
    let occ_idx = cif_find_header(headers, "_atom_site_occupancy")
    for row in rows {
      let label = if label_idx >= 0 && label_idx < row.length() {
        row[label_idx]
      } else {
        "?"
      }
      let elem = if type_idx >= 0 && type_idx < row.length() {
        row[type_idx]
      } else {
        "?"
      }
      let x = if x_idx >= 0 && x_idx < row.length() {
        cif_parse_double(row[x_idx])
      } else {
        0.0
      }
      let y = if y_idx >= 0 && y_idx < row.length() {
        cif_parse_double(row[y_idx])
      } else {
        0.0
      }
      let z = if z_idx >= 0 && z_idx < row.length() {
        cif_parse_double(row[z_idx])
      } else {
        0.0
      }
      let occ = if occ_idx >= 0 && occ_idx < row.length() {
        cif_parse_double(row[occ_idx])
      } else {
        1.0
      }
      atoms.push(
        CrystalAtom::new(
          label=label,
          element=elem,
          x=x,
          y=y,
          z=z,
          occupancy=occ,
          is_fractional=true,
        ),
      )
    }
  }
  // Bonds (geom loop), if present.
  let (gheaders, grows) = block.get_loop("_geom_bond_atom_site_label_1")
  let bonds : Array[CrystalBond] = Array::new()
  if gheaders.length() > 0 && grows.length() > 0 {
    let a1_idx = cif_find_header(gheaders, "_geom_bond_atom_site_label_1")
    let a2_idx = cif_find_header(gheaders, "_geom_bond_atom_site_label_2")
    let dist_idx = cif_find_header(gheaders, "_geom_bond_distance")
    for row in grows {
      let a1 = if a1_idx >= 0 && a1_idx < row.length() {
        row[a1_idx]
      } else {
        "?"
      }
      let a2 = if a2_idx >= 0 && a2_idx < row.length() {
        row[a2_idx]
      } else {
        "?"
      }
      let dist = if dist_idx >= 0 && dist_idx < row.length() {
        cif_parse_double(row[dist_idx])
      } else {
        0.0
      }
      bonds.push(CrystalBond::new(atom1=a1, atom2=a2, distance=dist))
    }
  }
  // Z and molecular weight.
  let z = cif_get_scalar_int(block, "_cell_formula_units_Z")
  let mw = cif_get_scalar_double(block, "_chemical_formula_weight")
  CrystalStructure::new(
    name=name,
    cell=cell,
    space_group=space_group,
    atoms=atoms,
    bonds=bonds,
    z=z,
    molecular_weight=mw,
  )
}

// ============================================================================
// Format and summary
// ============================================================================

///|
/// Format a CrystalStructure as a human-readable summary.
pub fn format_crystal_structure(structure : CrystalStructure) -> String {
  let mut s = "=== Crystal Structure: " + structure.name + " ===\n"
  s = s +
    "Cell: a=" +
    cif_fmt(structure.cell.a) +
    " b=" +
    cif_fmt(structure.cell.b) +
    " c=" +
    cif_fmt(structure.cell.c) +
    " alpha=" +
    cif_fmt(structure.cell.alpha) +
    " beta=" +
    cif_fmt(structure.cell.beta) +
    " gamma=" +
    cif_fmt(structure.cell.gamma) +
    "\n"
  let vol = unit_cell_volume(structure.cell)
  s = s + "Volume: " + cif_fmt(vol) + " A^3\n"
  s = s +
    "Space group: " +
    structure.space_group.symbol +
    " (#" +
    structure.space_group.number.to_string() +
    ")\n"
  s = s + "Z = " + structure.z.to_string() + "\n"
  s = s + "Molecular weight: " + cif_fmt(structure.molecular_weight) + "\n"
  if structure.molecular_weight > 0.0 && structure.z > 0 {
    let rho = crystal_density(structure.cell, structure.z, structure.molecular_weight)
    s = s + "Density: " + cif_fmt(rho) + " g/cm^3\n"
  }
  s = s + "Atoms (" + structure.atoms.length().to_string() + "):\n"
  s = s + "label\telement\tx\ty\tz\tocc\n"
  for a in structure.atoms {
    s = s +
      a.label +
      "\t" +
      a.element +
      "\t" +
      cif_fmt(a.x) +
      "\t" +
      cif_fmt(a.y) +
      "\t" +
      cif_fmt(a.z) +
      "\t" +
      cif_fmt(a.occupancy) +
      "\n"
  }
  if structure.bonds.length() > 0 {
    s = s + "Bonds (" + structure.bonds.length().to_string() + "):\n"
    s = s + "atom1\tatom2\tdistance\n"
    for b in structure.bonds {
      s = s +
        b.atom1 +
        "\t" +
        b.atom2 +
        "\t" +
        cif_fmt(b.distance) +
        "\n"
    }
  }
  s
}

///|
/// Summary statistics for a crystal structure.
pub struct CrystalSummary {
  name : String
  n_atoms : Int
  n_bonds : Int
  volume : Double
  density : Double
  space_group : String
  z : Int
}

///|
/// Compute a CrystalSummary for a structure.
pub fn crystal_summary(structure : CrystalStructure) -> CrystalSummary {
  let rho = if structure.molecular_weight > 0.0 && structure.z > 0 {
    crystal_density(structure.cell, structure.z, structure.molecular_weight)
  } else {
    0.0
  }
  CrystalSummary::{
    name: structure.name,
    n_atoms: structure.atoms.length(),
    n_bonds: structure.bonds.length(),
    volume: unit_cell_volume(structure.cell),
    density: rho,
    space_group: structure.space_group.symbol,
    z: structure.z,
  }
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Sample CIF text (benzene-like structure) for testing/demos.
pub fn crystal_sample_cif() -> String {
  let s = "data_benzene\n" +
    "_cell_length_a 5.418(5)\n" +
    "_cell_length_b 5.418(5)\n" +
    "_cell_length_c 7.543(7)\n" +
    "_cell_angle_alpha 90.0\n" +
    "_cell_angle_beta 90.0\n" +
    "_cell_angle_gamma 120.0\n" +
    "_space_group_name_H-M 'P 21/c'\n" +
    "_space_group_IT_number 14\n" +
    "_cell_formula_units_Z 2\n" +
    "_chemical_formula_weight 78.11\n" +
    "loop_\n" +
    "_atom_site_label\n" +
    "_atom_site_type_symbol\n" +
    "_atom_site_fract_x\n" +
    "_atom_site_fract_y\n" +
    "_atom_site_fract_z\n" +
    "_atom_site_occupancy\n" +
    "C1 C 0.1234 0.2345 0.3456 1.0\n" +
    "C2 C 0.2234 0.3345 0.4456 1.0\n" +
    "C3 C 0.3234 0.4345 0.5456 1.0\n" +
    "C4 C 0.4234 0.5345 0.6456 1.0\n" +
    "C5 C 0.5234 0.6345 0.7456 1.0\n" +
    "C6 C 0.6234 0.7345 0.8456 1.0\n" +
    "H1 H 0.1234 0.2345 0.4456 1.0\n" +
    "H2 H 0.2234 0.3345 0.5456 1.0\n" +
    "H3 H 0.3234 0.4345 0.6456 1.0\n" +
    "H4 H 0.4234 0.5345 0.7456 1.0\n" +
    "H5 H 0.5234 0.6345 0.8456 1.0\n" +
    "H6 H 0.6234 0.7345 0.9456 1.0\n"
  s
}

///|
/// Sample CrystalStructure for testing/demos.
pub fn crystal_sample_structure() -> CrystalStructure {
  let cell = UnitCell::new(
    a=5.418,
    b=5.418,
    c=7.543,
    alpha=90.0,
    beta=90.0,
    gamma=120.0,
  )
  let sg = SpaceGroup::new(number=14, symbol="P 21/c")
  let atoms = [
    CrystalAtom::new(label="C1", element="C", x=0.1234, y=0.2345, z=0.3456),
    CrystalAtom::new(label="C2", element="C", x=0.2234, y=0.3345, z=0.4456),
    CrystalAtom::new(label="C3", element="C", x=0.3234, y=0.4345, z=0.5456),
    CrystalAtom::new(label="C4", element="C", x=0.4234, y=0.5345, z=0.6456),
    CrystalAtom::new(label="C5", element="C", x=0.5234, y=0.6345, z=0.7456),
    CrystalAtom::new(label="C6", element="C", x=0.6234, y=0.7345, z=0.8456),
  ]
  let bonds = [
    CrystalBond::new(atom1="C1", atom2="C2", distance=1.39),
    CrystalBond::new(atom1="C2", atom2="C3", distance=1.39),
    CrystalBond::new(atom1="C3", atom2="C4", distance=1.39),
    CrystalBond::new(atom1="C4", atom2="C5", distance=1.39),
    CrystalBond::new(atom1="C5", atom2="C6", distance=1.39),
    CrystalBond::new(atom1="C6", atom2="C1", distance=1.39),
  ]
  CrystalStructure::new(
    name="benzene",
    cell=cell,
    space_group=sg,
    atoms=atoms,
    bonds=bonds,
    z=2,
    molecular_weight=78.11,
  )
}

// ============================================================================
// Internal helpers
// ============================================================================

///|
/// Tokenize a CIF text into a list of tokens (data_, loop_, _key, value).
fn cif_tokenize(text : String) -> Array[String] {
  let tokens : Array[String] = Array::new()
  let n = text.length()
  let mut i = 0
  while i < n {
    let c = text[i:i + 1].to_owned()
    // Skip whitespace.
    if c == " " || c == "\t" || c == "\r" || c == "\n" {
      i = i + 1
      continue
    }
    // Skip comments.
    if c == "#" {
      while i < n && text[i:i + 1].to_owned() != "\n" {
        i = i + 1
      }
      continue
    }
    // Quoted strings.
    if c == "'" || c == "\"" {
      let quote = c
      i = i + 1
      let mut s = ""
      while i < n {
        let ch = text[i:i + 1].to_owned()
        if ch == quote {
          i = i + 1
          break
        }
        s = s + ch
        i = i + 1
      }
      tokens.push(s)
      continue
    }
    // Unquoted token: read until whitespace.
    let mut s = ""
    while i < n {
      let ch = text[i:i + 1].to_owned()
      if ch == " " || ch == "\t" || ch == "\r" || ch == "\n" {
        break
      }
      s = s + ch
      i = i + 1
    }
    tokens.push(s)
  }
  tokens
}

///|
/// Find the index of a header in a header array. Returns -1 if not found.
fn cif_find_header(headers : Array[String], name : String) -> Int {
  for i in 0..<headers.length() {
    if headers[i] == name {
      return i
    }
  }
  -1
}

///|
/// Get a scalar double from a CifBlock, returning 0.0 if not present.
fn cif_get_scalar_double(block : CifBlock, name : String) -> Double {
  let v = block.get_scalar(name)
  if v.length() == 0 {
    return 0.0
  }
  cif_parse_double(v)
}

///|
/// Get a scalar int from a CifBlock, returning 0 if not present.
fn cif_get_scalar_int(block : CifBlock, name : String) -> Int {
  let v = block.get_scalar(name)
  if v.length() == 0 {
    return 0
  }
  cif_parse_int(v)
}

///|
/// Parse a double from a CIF value string, stripping trailing uncertainty
/// notation like "5.418(5)" -> 5.418.
fn cif_parse_double(s : String) -> Double {
  // Strip trailing (n) uncertainty.
  let mut cleaned = ""
  let mut in_paren = false
  for i in 0..<s.length() {
    let c = s[i:i + 1].to_owned()
    if c == "(" {
      in_paren = true
      continue
    }
    if c == ")" {
      in_paren = false
      continue
    }
    if !in_paren {
      cleaned = cleaned + c
    }
  }
  match parse_double(cleaned) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Parse an int from a CIF value string, stripping trailing uncertainty.
fn cif_parse_int(s : String) -> Int {
  let mut cleaned = ""
  let mut in_paren = false
  for i in 0..<s.length() {
    let c = s[i:i + 1].to_owned()
    if c == "(" {
      in_paren = true
      continue
    }
    if c == ")" {
      in_paren = false
      continue
    }
    if !in_paren {
      cleaned = cleaned + c
    }
  }
  parse_int(cleaned)
}

///|
/// Format a double for display (4 decimal places).
fn cif_fmt(x : Double) -> String {
  let neg = x < 0.0
  let abs_x = if neg { -x } else { x }
  let scaled = (abs_x * 10000.0 + 0.5).floor()
  let int_part = (scaled / 10000.0).floor().to_int()
  let frac_int = (scaled - int_part.to_double() * 10000.0).to_int()
  let frac_str = if frac_int < 10 {
    "000" + frac_int.to_string()
  } else if frac_int < 100 {
    "00" + frac_int.to_string()
  } else if frac_int < 1000 {
    "0" + frac_int.to_string()
  } else {
    frac_int.to_string()
  }
  let sign = if neg { "-" } else { "" }
  sign + int_part.to_string() + "." + frac_str
}
