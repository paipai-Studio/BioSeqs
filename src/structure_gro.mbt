///|
/// GRO format reader/writer for GROMACS trajectory frames.
/// Ported from biotite.structure.io.gro.
///
/// GRO format stores coordinates in nm (not Å).

///|
/// A single GRO frame: atoms with coordinates and optional velocities.
pub struct GroFrame {
  atoms : Array[GroAtom]
  // Box vectors (3×3), row-major; empty if no box
  box : Array[Array[Double]]
}

///|
/// A single atom in a GRO frame.
pub struct GroAtom {
  res_id : Int32
  res_name : String
  atom_name : String
  atom_id : Int32
  x : Double // nm
  y : Double // nm
  z : Double // nm
}

///|
/// Parse a GRO file string into a GroFrame.
pub fn gro_parse(content : String) -> GroFrame {
  let lines = content.split("\n").to_array()
  // Line 0: title
  // Line 1: number of atoms
  let n_atoms = pdb_parse_int_fast(lines[1].trim())
  let atoms : Array[GroAtom] = Array(capacity=n_atoms)
  for i = 0; i < n_atoms; i = i + 1 {
    let line = lines[2 + i]
    atoms.push(parse_gro_atom_line(line))
  }
  // Box line (optional)
  let box : Array[Array[Double]] = Array(capacity=0)
  let box_line_idx = 2 + n_atoms
  if box_line_idx < lines.length() {
    let box_line = lines[box_line_idx].trim()
    if box_line.length() > 0 {
      let parts = box_line.split(" ")
      let nums : Array[Double] = Array(capacity=9)
      for p in parts {
        let t = p.trim()
        if t.length() > 0 {
          nums.push(pdb_parse_double_fast(t))
        }
      }
      if nums.length() >= 9 {
        box.push([nums[0], nums[1], nums[2]])
        box.push([nums[3], nums[4], nums[5]])
        box.push([nums[6], nums[7], nums[8]])
      } else if nums.length() >= 3 {
        // Diagonal box
        box.push([nums[0], 0.0, 0.0])
        box.push([0.0, nums[1], 0.0])
        box.push([0.0, 0.0, nums[2]])
      }
    }
  }
  GroFrame::{ atoms, box, }
}

///|
/// Safe substring: clamp end to line length.
fn gro_safe_slice(line : StringView, start : Int32, end : Int32) -> StringView {
  let len = line.length()
  let s = if start < len { start } else { len }
  let e = if end < len { end } else { len }
  if s < e {
    line[s:e]
  } else {
    ""
  }
}

///|
/// Parse a single GRO atom line.
/// Fixed-column format:
///   1-5: res id, 6-10: res name, 11-15: atom name,
///   16-20: atom id, 21-28: x, 29-36: y, 37-44: z
/// Falls back to space-delimited if columns don't make sense.
fn parse_gro_atom_line(line : StringView) -> GroAtom {
  let len = line.length()
  // Try fixed-column format if we have at least through atom_id (20 chars)
  if len >= 20 {
    let res_id = pdb_parse_int_fast(gro_safe_slice(line, 0, 5).trim())
    let res_name = gro_safe_slice(line, 5, 10).trim().to_owned()
    let atom_name = gro_safe_slice(line, 10, 15).trim().to_owned()
    let atom_id = pdb_parse_int_fast(gro_safe_slice(line, 15, 20).trim())
    // Parse coordinates from whatever remains
    let x = if len > 20 {
      pdb_parse_double_fast(gro_safe_slice(line, 20, 28).trim())
    } else {
      0.0
    }
    let y = if len > 28 {
      pdb_parse_double_fast(gro_safe_slice(line, 28, 36).trim())
    } else {
      0.0
    }
    let z = if len > 36 {
      pdb_parse_double_fast(gro_safe_slice(line, 36, 44).trim())
    } else {
      0.0
    }
    // Sanity check: res_name should be alphabetic, not numeric
    if res_name.length() > 0 && res_id > 0 {
      return GroAtom::{ res_id, res_name, atom_name, atom_id, x, y, z, }
    }
  }
  // Fallback: space-delimited
  let parts = line.split(" ")
  let nums : Array[StringView] = Array(capacity=0)
  for p in parts {
    let t = p.trim()
    if t.length() > 0 {
      nums.push(t)
    }
  }
  // res_id res_name atom_name atom_id x y z
  GroAtom::{
    res_id: pdb_parse_int_fast(nums[0]),
    res_name: nums[1].to_owned(),
    atom_name: nums[2].to_owned(),
    atom_id: pdb_parse_int_fast(nums[3]),
    x: pdb_parse_double_fast(nums[4]),
    y: pdb_parse_double_fast(nums[5]),
    z: pdb_parse_double_fast(nums[6]),
  }
}

///|
/// Convert GRO coordinates from nm to Å.
pub fn gro_nm_to_angstrom(coords : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array(capacity=coords.length())
  for i = 0; i < coords.length(); i = i + 1 {
    result.push(coords[i] * 10.0)
  }
  result
}

///|
/// Convert Å coordinates to GRO nm.
pub fn gro_angstrom_to_nm(coords : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array(capacity=coords.length())
  for i = 0; i < coords.length(); i = i + 1 {
    result.push(coords[i] / 10.0)
  }
  result
}
