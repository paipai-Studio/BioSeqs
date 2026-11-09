///|
/// LAMMPS custom dump text format reader.
///
/// Supports orthogonal (`ITEM: BOX BOUNDS`) and triclinic boxes
/// (`ITEM: BOX BOUNDS xy xz yz`, where each bound line carries the tilt
/// factor as its third value: `xlo xhi xy`, `ylo yhi xz`, `zlo zhi yz`).
/// Atom columns are parsed from the `ITEM: ATOMS ...` header; `x/y/z`
/// (and the scaled `xs/ys/zs` / unwrapped `xu/yu/zu` aliases) are treated
/// as cartesian coordinates, every other numeric column is collected
/// into `extra` in header order. `mol` defaults to 0 when absent, `type`
/// defaults to 0, and `id` defaults to the 1-based row index.

///|
/// Error type for LAMMPS dump parsing.
pub suberror LammpsDumpError {
  LammpsDumpError(String)
}

///|
/// One atom entry of a LAMMPS dump frame.
pub(all) struct LammpsAtom {
  id : Int32
  mol : Int32
  atom_type : Int32
  x : Double
  y : Double
  z : Double
  extra : Array[Double]
}

///|
/// One parsed LAMMPS dump frame. `box_lo`/`box_hi` hold 3 values each;
/// `tilt` holds xy/xz/yz for triclinic boxes and is empty for orthogonal
/// ones.
pub(all) struct LammpsFrame {
  timestep : Int32
  natoms : Int32
  box_lo : Array[Double]
  box_hi : Array[Double]
  tilt : Array[Double]
  atoms : Array[LammpsAtom]
}

///|
/// Split a line on spaces, dropping empty tokens.
fn lammps_tokens(line : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for tok in split_by_char(line, ' '.to_int()) {
    if tok.length() > 0 {
      out.push(tok)
    }
  }
  out
}

///|
/// Return the next non-empty (trimmed) line starting at `start`, together
/// with the index just past it; returns ("", lines.length()) at end of
/// input.
fn lammps_next_line(lines : Array[String], start : Int32) -> (String, Int32) {
  let n = lines.length()
  let mut k = start
  while k < n {
    let t = lines[k].trim().to_owned()
    if t.length() > 0 {
      return (t, k + 1)
    }
    k = k + 1
  }
  ("", n)
}

///|
/// Parse token `j` as a double, returning 0.0 when missing or malformed.
fn lammps_num(toks : Array[String], j : Int32) -> Double {
  if j >= 0 && j < toks.length() {
    parse_double(toks[j]).unwrap_or(0.0)
  } else {
    0.0
  }
}

///|
/// Parse a LAMMPS custom dump text (possibly holding many frames) into
/// frames. Raises LammpsDumpError on structural errors such as a missing
/// section header, a truncated box/atom block, or an atom block without
/// x/y/z columns.
pub fn lammps_dump_read(
  text : String,
) -> Array[LammpsFrame] raise LammpsDumpError {
  let lines = split_by_char(text, '\n'.to_int())
  let frames : Array[LammpsFrame] = Array(capacity=0)
  let n = lines.length()
  let mut i = 0
  while i < n {
    let t = lines[i].trim().to_owned()
    if t.has_prefix("ITEM: TIMESTEP") {
      let (ts_line, p1) = lammps_next_line(lines, i + 1)
      if ts_line.length() == 0 {
        raise LammpsDumpError("unexpected end of file after ITEM: TIMESTEP")
      }
      let timestep = lammps_num(lammps_tokens(ts_line), 0).to_int()
      let (na_line, p2) = lammps_next_line(lines, p1)
      if !na_line.has_prefix("ITEM: NUMBER OF ATOMS") {
        raise LammpsDumpError("expected 'ITEM: NUMBER OF ATOMS'")
      }
      let (nn_line, p3) = lammps_next_line(lines, p2)
      if nn_line.length() == 0 {
        raise LammpsDumpError(
          "unexpected end of file after ITEM: NUMBER OF ATOMS",
        )
      }
      let natoms = lammps_num(lammps_tokens(nn_line), 0).to_int()
      let (bb_line, p4) = lammps_next_line(lines, p3)
      if !bb_line.has_prefix("ITEM: BOX BOUNDS") {
        raise LammpsDumpError("expected 'ITEM: BOX BOUNDS'")
      }
      let triclinic = bb_line.contains("xy xz yz")
      let box_lo : Array[Double] = Array(capacity=3)
      let box_hi : Array[Double] = Array(capacity=3)
      let tilt : Array[Double] = Array(capacity=0)
      let mut p = p4
      for k = 0; k < 3; k = k + 1 {
        let (bl, pn) = lammps_next_line(lines, p)
        if bl.length() == 0 {
          raise LammpsDumpError("unexpected end of file in box bounds")
        }
        let toks = lammps_tokens(bl)
        if toks.length() < 2 {
          raise LammpsDumpError("invalid box bounds line '\{bl}'")
        }
        box_lo.push(lammps_num(toks, 0))
        box_hi.push(lammps_num(toks, 1))
        if triclinic {
          if toks.length() < 3 {
            raise LammpsDumpError(
              "triclinic box bounds line missing tilt factor",
            )
          }
          tilt.push(lammps_num(toks, 2))
        }
        p = pn
      }
      let (at_line, p5) = lammps_next_line(lines, p)
      if !at_line.has_prefix("ITEM: ATOMS") {
        raise LammpsDumpError("expected 'ITEM: ATOMS'")
      }
      let head_toks = lammps_tokens(at_line)
      let cols : Array[String] = Array(capacity=0)
      for j = 2; j < head_toks.length(); j = j + 1 {
        cols.push(head_toks[j])
      }
      let mut id_col = -1
      let mut mol_col = -1
      let mut type_col = -1
      let mut x_col = -1
      let mut y_col = -1
      let mut z_col = -1
      let extra_cols : Array[Int32] = Array(capacity=0)
      for j = 0; j < cols.length(); j = j + 1 {
        let c = cols[j]
        if c == "id" {
          id_col = j
        } else if c == "mol" {
          mol_col = j
        } else if c == "type" {
          type_col = j
        } else if c == "x" || c == "xs" || c == "xu" {
          x_col = j
        } else if c == "y" || c == "ys" || c == "yu" {
          y_col = j
        } else if c == "z" || c == "zs" || c == "zu" {
          z_col = j
        } else {
          extra_cols.push(j)
        }
      }
      if x_col < 0 || y_col < 0 || z_col < 0 {
        raise LammpsDumpError("missing x/y/z columns in ITEM: ATOMS")
      }
      p = p5
      let atoms : Array[LammpsAtom] = Array(capacity=natoms)
      for a = 0; a < natoms; a = a + 1 {
        let (al, pn) = lammps_next_line(lines, p)
        if al.length() == 0 {
          raise LammpsDumpError("unexpected end of file in atom block")
        }
        p = pn
        let toks = lammps_tokens(al)
        let id = if id_col >= 0 {
          lammps_num(toks, id_col).to_int()
        } else {
          a + 1
        }
        let mol = if mol_col >= 0 {
          lammps_num(toks, mol_col).to_int()
        } else {
          0
        }
        let atype = if type_col >= 0 {
          lammps_num(toks, type_col).to_int()
        } else {
          0
        }
        let extra : Array[Double] = Array(capacity=0)
        for ec in extra_cols {
          extra.push(lammps_num(toks, ec))
        }
        atoms.push(LammpsAtom::{
          id,
          mol,
          atom_type: atype,
          x: lammps_num(toks, x_col),
          y: lammps_num(toks, y_col),
          z: lammps_num(toks, z_col),
          extra,
        })
      }
      frames.push(LammpsFrame::{
        timestep,
        natoms,
        box_lo,
        box_hi,
        tilt,
        atoms,
      })
      i = p
    } else {
      i = i + 1
    }
  }
  frames
}

///|
/// Convert parsed LAMMPS frames into a Trajectory using the cartesian
/// atom coordinates.
pub fn lammps_to_trajectory(frames : Array[LammpsFrame]) -> Trajectory {
  let out : Array[Array[Array[Double]]] = Array(capacity=0)
  for f in frames {
    let frame : Array[Array[Double]] = Array(capacity=f.atoms.length())
    for atom in f.atoms {
      frame.push([atom.x, atom.y, atom.z])
    }
    out.push(frame)
  }
  Trajectory::new(out)
}
