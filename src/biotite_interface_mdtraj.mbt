///|
/// biotite_interface_mdtraj.mbt — lightweight MDTraj-style topology and
/// coordinate adapters without importing MDTraj.

///|
/// One line per atom: index, residue name, residue ID, atom name, element.
pub fn biotite_mdtraj_topology(arr : AtomArray) -> String {
  let sb = StringBuilder()
  for i in 0..<arr.len() {
    let _ = sb.write(
      i.to_string() +
      " " +
      arr.get_res_name(i) +
      " " +
      arr.res_ids[i].to_string() +
      " " +
      arr.get_name(i) +
      " " +
      arr.elements[i] +
      "\n",
    )
  }
  sb.to_string()
}

///|
/// Rebuild an AtomArray from a coordinate block using a template's
/// annotation (MDTraj topology + xyz coordinates).
pub fn biotite_mdtraj_from_coords(
  template : AtomArray,
  coords : Array[Array[Double]],
) -> AtomArray {
  let n = if template.len() < coords.length() {
    template.len()
  } else {
    coords.length()
  }
  let atoms : Array[Atom] = Array(capacity=n)
  for i in 0..<n {
    let c = coords[i]
    let x = if c.length() > 0 { c[0] } else { 0.0 }
    let y = if c.length() > 1 { c[1] } else { 0.0 }
    let z = if c.length() > 2 { c[2] } else { 0.0 }
    atoms.push(
      Atom::new(
        name=template.get_name(i),
        coord=Vector3::new(x, y, z),
        resname=template.get_res_name(i),
        chainid=template.get_chain_id(i),
        resseq=template.res_ids[i],
        occupancy=template.occupancies[i],
        bfactor=template.b_factors[i],
        element=template.elements[i],
      ),
    )
  }
  AtomArray::from_atoms(atoms)
}

///|
/// Extract coordinates as [atom][x, y, z] for MDTraj-style consumers.
pub fn biotite_mdtraj_coords(arr : AtomArray) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array(capacity=arr.len())
  for i in 0..<arr.len() {
    let c = arr.get_coord(i)
    result.push([c.x, c.y, c.z])
  }
  result
}
