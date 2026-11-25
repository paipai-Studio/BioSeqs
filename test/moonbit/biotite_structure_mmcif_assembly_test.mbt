///|
/// Tests for biotite_structure_mmcif_assembly.mbt.
fn bma_mmcif() -> String {
  "data_test\n" +
  "loop_\n" +
  "_pdbx_struct_assembly.id\n" +
  "_pdbx_struct_assembly.oligomeric_count\n" +
  "1 2\n" +
  "#\n" +
  "loop_\n" +
  "_pdbx_struct_oper_list.id\n" +
  "_pdbx_struct_oper_list.matrix[1][1]\n" +
  "_pdbx_struct_oper_list.matrix[1][2]\n" +
  "_pdbx_struct_oper_list.matrix[1][3]\n" +
  "_pdbx_struct_oper_list.matrix[2][1]\n" +
  "_pdbx_struct_oper_list.matrix[2][2]\n" +
  "_pdbx_struct_oper_list.matrix[2][3]\n" +
  "_pdbx_struct_oper_list.matrix[3][1]\n" +
  "_pdbx_struct_oper_list.matrix[3][2]\n" +
  "_pdbx_struct_oper_list.matrix[3][3]\n" +
  "_pdbx_struct_oper_list.vector[1]\n" +
  "_pdbx_struct_oper_list.vector[2]\n" +
  "_pdbx_struct_oper_list.vector[3]\n" +
  "1 1 0 0 0 1 0 0 0 1 10 0 0\n" +
  "#\n"
}

///|
test "biotite_parse_assembly_operators" {
  let ops = @src.biotite_parse_assembly_operators(bma_mmcif())
  assert_eq(ops.length(), 1)
  assert_eq(ops[0].id, "1")
  assert_true(@moonbit.assert_num_close(ops[0].matrix[0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ops[0].translation[0], 10.0, 1.0e-12))
}

///|
test "biotite_apply_assembly_operator" {
  let ops = @src.biotite_parse_assembly_operators(bma_mmcif())
  let arr = @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(1.0, 2.0, 3.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
  ])
  let moved = @src.biotite_apply_assembly_operator(arr, ops[0])
  assert_eq(moved.len(), 1)
  assert_true(@moonbit.assert_num_close(moved.get_coord(0).x, 11.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(moved.get_coord(0).y, 2.0, 1.0e-12))
}

///|
test "biotite_assembly_oligomeric_count" {
  assert_eq(@src.biotite_assembly_oligomeric_count(bma_mmcif()), 2)
}

///|
test "biotite_parse_assembly_operators_absent" {
  assert_eq(@src.biotite_parse_assembly_operators("data_empty\n").length(), 0)
  assert_eq(@src.biotite_assembly_oligomeric_count("data_empty\n"), 0)
}
