///|
/// biotite_structure_mmcif_assembly.mbt — mmCIF biological assembly
/// operators and coordinate application.

///|
pub(all) struct BiotiteAssemblyOperator {
  id : String
  matrix : Array[Double]
  translation : Array[Double]
}

///|
/// Parse `_pdbx_struct_oper_list` operators from an mmCIF text.
pub fn biotite_parse_assembly_operators(
  text : String,
) -> Array[BiotiteAssemblyOperator] {
  let block = parse_mmcif(text)
  let category = "pdbx_struct_oper_list"
  let headers = block.get_category_headers(category)
  if headers.length() == 0 {
    return []
  }
  let rows = block.get_category_rows(category)
  let id_idx = biotite_find_header(headers, "id")
  let m11 = biotite_find_header(headers, "matrix[1][1]")
  let m12 = biotite_find_header(headers, "matrix[1][2]")
  let m13 = biotite_find_header(headers, "matrix[1][3]")
  let m21 = biotite_find_header(headers, "matrix[2][1]")
  let m22 = biotite_find_header(headers, "matrix[2][2]")
  let m23 = biotite_find_header(headers, "matrix[2][3]")
  let m31 = biotite_find_header(headers, "matrix[3][1]")
  let m32 = biotite_find_header(headers, "matrix[3][2]")
  let m33 = biotite_find_header(headers, "matrix[3][3]")
  let v1 = biotite_find_header(headers, "vector[1]")
  let v2 = biotite_find_header(headers, "vector[2]")
  let v3 = biotite_find_header(headers, "vector[3]")
  let result : Array[BiotiteAssemblyOperator] = []
  for row in rows {
    if id_idx < 0 || id_idx >= row.length() {
      continue
    }
    result.push({
      id: row[id_idx],
      matrix: [
        biotite_row_double(row, m11),
        biotite_row_double(row, m12),
        biotite_row_double(row, m13),
        biotite_row_double(row, m21),
        biotite_row_double(row, m22),
        biotite_row_double(row, m23),
        biotite_row_double(row, m31),
        biotite_row_double(row, m32),
        biotite_row_double(row, m33),
      ],
      translation: [
        biotite_row_double(row, v1),
        biotite_row_double(row, v2),
        biotite_row_double(row, v3),
      ],
    })
  }
  result
}

///|
fn biotite_find_header(headers : Array[String], suffix : String) -> Int32 {
  for i in 0..<headers.length() {
    if headers[i].has_suffix("." + suffix) {
      return i
    }
  }
  -1
}

///|
fn biotite_row_double(row : Array[String], idx : Int32) -> Double {
  if idx < 0 || idx >= row.length() {
    return 0.0
  }
  match parse_double(row[idx].trim().to_owned()) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Apply an assembly operator (fractional/cartesian rotation + translation)
/// to every coordinate.
pub fn biotite_apply_assembly_operator(
  arr : AtomArray,
  op : BiotiteAssemblyOperator,
) -> AtomArray {
  let m = op.matrix
  let t = op.translation
  let atoms = arr.to_atoms()
  let result : Array[Atom] = Array(capacity=atoms.length())
  for atom in atoms {
    let c = atom.coord
    let x = m[0] * c.x + m[1] * c.y + m[2] * c.z + t[0]
    let y = m[3] * c.x + m[4] * c.y + m[5] * c.z + t[1]
    let z = m[6] * c.x + m[7] * c.y + m[8] * c.z + t[2]
    result.push(
      Atom::new(
        name=atom.name,
        coord=Vector3::new(x, y, z),
        resname=atom.resname,
        chainid=atom.chainid,
        resseq=atom.resseq,
        occupancy=atom.occupancy,
        bfactor=atom.bfactor,
        element=atom.element,
      ),
    )
  }
  AtomArray::from_atoms(result)
}

///|
/// Read `_pdbx_struct_assembly.oligomeric_count` (0 if absent).
pub fn biotite_assembly_oligomeric_count(text : String) -> Int32 {
  let block = parse_mmcif(text)
  let category = "pdbx_struct_assembly"
  let headers = block.get_category_headers(category)
  let rows = block.get_category_rows(category)
  if headers.length() == 0 || rows.length() == 0 {
    return 0
  }
  let idx = biotite_find_header(headers, "oligomeric_count")
  if idx < 0 || idx >= rows[0].length() {
    return 0
  }
  parse_int(rows[0][idx].trim().to_owned())
}
