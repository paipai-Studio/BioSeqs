/// Biotite `structure.integrity` — checks for errors in atom arrays.
///
/// Faithful port of the official implementation
/// (`src/biotite/structure/integrity.py` in biotite-dev/biotite):
/// - continuity checks flag positions where the difference between
///   successive IDs is neither 0 nor 1; the returned indices point to
///   the atom *after* a gap,
/// - `check_linear_continuity` flags atoms whose bond to the
///   *preceding* atom is outside `[min_len, max_len]`
///   (default 1.2-1.8 A), mirroring the official
///   `filter_linear_bond_continuity` + shift construction,
/// - `check_backbone_continuity` restricts the same distance check to
///   peptide backbone atoms ("N", "CA", "C" of amino acids) and
///   phosphate backbone atoms ("P", "O5'", "C5'", "C4'", "C3'", "O3'"
///   of nucleotides),
/// - `check_duplicate_atoms` flags atoms whose annotations equal the
///   annotations of an earlier atom (the first occurrence is kept;
///   coordinates may differ),
/// - `check_in_box` returns the indices of atoms inside the box
///   (fractional coordinates within [0, 1) on all axes) — note the
///   official docstring claims "outside" indices while its
///   implementation returns the inside ones; the implementation
///   semantics are ported here,
/// - a missing box raises upstream (`TypeError`); this port has no
///   exceptions, so calling `check_in_box` without a box yields an
///   empty result.

///|
/// Generic continuity check over an ID sequence: positions where the
/// difference between successive values is neither 0 nor 1
/// (official `_check_continuity`).
fn biotite_integrity_continuity(ids : Array[Int32]) -> Array[Int32] {
  let out : Array[Int32] = []
  for i = 0; i < ids.length() - 1; i = i + 1 {
    let diff = ids[i + 1] - ids[i]
    if diff != 0 && diff != 1 {
      out.push(i + 1)
    }
  }
  out
}

///|
/// Check if the atom IDs are incremented by more than 1 or
/// decremented from one atom to the next one. The local AtomArray has
/// no atom_id annotation, so the ID sequence is passed explicitly
/// (officially read from `array.atom_id`).
pub fn biotite_check_atom_id_continuity(
  atom_ids : Array[Int32],
) -> Array[Int32] {
  biotite_integrity_continuity(atom_ids)
}

///|
/// Check if the residue IDs are incremented by more than 1 or
/// decremented from one atom to the next one.
pub fn biotite_check_res_id_continuity(arr : AtomArray) -> Array[Int32] {
  biotite_integrity_continuity(arr.res_ids)
}

///|
/// Check linear (consecutive) bond continuity: flags atoms where the
/// distance to the *preceding* atom is outside `[min_len, max_len]`.
pub fn biotite_check_linear_continuity(
  arr : AtomArray,
  min_len~ : Double = 1.2,
  max_len~ : Double = 1.8,
) -> Array[Int32] {
  let out : Array[Int32] = []
  for i = 1; i < arr.len(); i = i + 1 {
    let d = arr.coords[i - 1].distance(arr.coords[i])
    if d < min_len || d > max_len {
      out.push(i)
    }
  }
  out
}

///|
/// Whether atom `i` is a peptide backbone atom ("N", "CA" or "C" of
/// an amino acid; official `filter_peptide_backbone`).
fn biotite_integrity_is_peptide_backbone(arr : AtomArray, i : Int32) -> Bool {
  let name = arr.names[i]
  if name != "N" && name != "CA" && name != "C" {
    return false
  }
  biotite_is_amino_acid(arr.res_names[i])
}

///|
/// Whether atom `i` is a phosphate backbone atom ("P", "O5'", "C5'",
/// "C4'", "C3'" or "O3'" of a nucleotide; official
/// `filter_phosphate_backbone`).
fn biotite_integrity_is_phosphate_backbone(arr : AtomArray, i : Int32) -> Bool {
  let name = arr.names[i]
  if name != "P" &&
    name != "O5'" &&
    name != "C5'" &&
    name != "C4'" &&
    name != "C3'" &&
    name != "O3'" {
    return false
  }
  biotite_is_nucleotide(arr.res_names[i])
}

///|
/// Check if the (peptide or phosphate) backbone atoms have
/// non-reasonable distance to the preceding backbone atom. Backbone
/// atoms form a subsequence of the array; consecutive atoms of that
/// subsequence are checked against `[min_len, max_len]` and
/// discontinuities are reported at their original indices.
pub fn biotite_check_backbone_continuity(
  arr : AtomArray,
  min_len~ : Double = 1.2,
  max_len~ : Double = 1.8,
) -> Array[Int32] {
  let out : Array[Int32] = []
  let mut prev = -1
  for i = 0; i < arr.len(); i = i + 1 {
    if !biotite_integrity_is_peptide_backbone(arr, i) &&
      !biotite_integrity_is_phosphate_backbone(arr, i) {
      continue
    }
    if prev >= 0 {
      let d = arr.coords[prev].distance(arr.coords[i])
      if d < min_len || d > max_len {
        out.push(i)
      }
    }
    prev = i
  }
  out
}

///|
/// Check if a structure contains duplicate atoms, i.e. two atoms
/// whose annotations are all equal (coordinates may differ). All
/// annotation fields of the local AtomArray participate: chain_id,
/// res_id, res_name, atom name, element, b_factor, occupancy and
/// altloc. The first occurrence of an atom is not reported.
pub fn biotite_check_duplicate_atoms(arr : AtomArray) -> Array[Int32] {
  let seen : Map[String, Bool] = Map([], capacity=16)
  let out : Array[Int32] = []
  for i = 0; i < arr.len(); i = i + 1 {
    let key = arr.chain_ids[i].to_string() +
      "|" +
      arr.res_ids[i].to_string() +
      "|" +
      arr.res_names[i] +
      "|" +
      arr.names[i] +
      "|" +
      arr.elements[i] +
      "|" +
      arr.b_factors[i].to_string() +
      "|" +
      arr.occupancies[i].to_string() +
      "|" +
      arr.altlocs[i].to_string()
    if seen.get(key).is_empty() {
      seen[key] = true
    } else {
      out.push(i)
    }
  }
  out
}

///|
/// Check which atoms lie inside the box: those whose fractional
/// coordinates are within [0, 1) on all axes (official
/// `check_in_box` implementation semantics; its docstring claims the
/// complement). A missing box raises upstream (`TypeError`); here the
/// result is empty.
pub fn biotite_check_in_box(
  arr : AtomArray,
  box~ : BiotiteBox? = None,
) -> Array[Int32] {
  let out : Array[Int32] = []
  match box {
    Some(b) => {
      let flat : Array[Double] = Array(capacity=arr.len() * 3)
      for p in arr.coords {
        flat.push(p.x)
        flat.push(p.y)
        flat.push(p.z)
      }
      let fractions = b.coord_to_fraction(flat)
      for i = 0; i < arr.len(); i = i + 1 {
        let fx = fractions[i * 3]
        let fy = fractions[i * 3 + 1]
        let fz = fractions[i * 3 + 2]
        if fx >= 0.0 &&
          fx < 1.0 &&
          fy >= 0.0 &&
          fy < 1.0 &&
          fz >= 0.0 &&
          fz < 1.0 {
          out.push(i)
        }
      }
    }
    None => ()
  }
  out
}

///|
pub fn _biotite_structure_integrity_dummy() -> Unit {
  let arr = AtomArray::new()
  let _ = biotite_check_atom_id_continuity([])
  let _ = biotite_check_res_id_continuity(arr)
  let _ = biotite_check_linear_continuity(arr)
  let _ = biotite_check_backbone_continuity(arr)
  let _ = biotite_check_duplicate_atoms(arr)
  let _ = biotite_check_in_box(arr)
}
