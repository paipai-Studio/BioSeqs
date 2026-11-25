///|
/// biotite_structure_io_general.mbt — format dispatch on top of the
/// repository's existing readers/writers.

///|
/// Load an AtomArray from text, selecting the reader by `format`.
///
/// Supported: pdb, xyz, cif/mmcif, gro, pdbqt.
/// Unknown formats return an empty AtomArray.
pub fn biotite_load_structure(text : String, format : String) -> AtomArray {
  match format {
    "pdb" => AtomArray::from_pdb_text(text)
    "xyz" => biotite_xyz_parse(text)
    "cif" | "mmcif" => AtomArray::from_mmcif(text)
    "gro" => biotite_gro_to_atom_array(gro_parse(text))
    "pdbqt" => biotite_pdbqt_to_atom_array(pdbqt_parse(text))
    _ => AtomArray::new()
  }
}

///|
/// Save an AtomArray to text in the requested format.
///
/// Supported: pdb, xyz, cif/mmcif. Other formats return an empty string.
pub fn biotite_save_structure(arr : AtomArray, format : String) -> String {
  match format {
    "pdb" => arr.to_pdb_string()
    "xyz" => biotite_xyz_write(arr)
    "cif" | "mmcif" => arr.to_mmcif()
    _ => ""
  }
}

///|
/// Guess a format from a file extension.
pub fn biotite_guess_format(path : String) -> String {
  if path.has_suffix(".pdb") {
    "pdb"
  } else if path.has_suffix(".xyz") {
    "xyz"
  } else if path.has_suffix(".cif") || path.has_suffix(".mmcif") {
    "mmcif"
  } else if path.has_suffix(".gro") {
    "gro"
  } else if path.has_suffix(".pdbqt") {
    "pdbqt"
  } else {
    ""
  }
}

///|
/// Convert a GRO frame (nm) to an AtomArray (Å).
fn biotite_gro_to_atom_array(frame : GroFrame) -> AtomArray {
  let atoms : Array[Atom] = Array(capacity=frame.atoms.length())
  for a in frame.atoms {
    atoms.push(
      Atom::new(
        name=a.atom_name,
        coord=Vector3::new(a.x * 10.0, a.y * 10.0, a.z * 10.0),
        resname=a.res_name,
        resseq=a.res_id,
        element="",
      ),
    )
  }
  AtomArray::from_atoms(atoms)
}

///|
/// Convert PDBQT atoms to an AtomArray.
fn biotite_pdbqt_to_atom_array(atoms_in : Array[PdbqtAtom]) -> AtomArray {
  let atoms : Array[Atom] = Array(capacity=atoms_in.length())
  for a in atoms_in {
    let mut chain = ' '
    for ch in a.chain {
      chain = ch
      break
    }
    atoms.push(
      Atom::new(
        name=a.atom_name,
        coord=Vector3::new(a.x, a.y, a.z),
        resname=a.res_name,
        chainid=chain,
        resseq=a.res_id,
        element=a.atom_type,
      ),
    )
  }
  AtomArray::from_atoms(atoms)
}
