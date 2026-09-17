///|
/// Tests for the Bio.PDB.MMCIF2Dict port. The sample document is the same one
/// cross-checked against Biopython 1.88.

///|
fn eq_strs(a : Array[String], b : Array[String]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return false
    }
  }
  true
}

///|
fn sample_cif() -> String {
  "data_test\n" +
  "# a comment line\n" +
  "_struct.entry_id 1ABC\n" +
  "_struct.title 'A \"quoted\" title'\n" +
  "loop_\n" +
  "_atom_site.group_PDB\n" +
  "_atom_site.id\n" +
  "_atom_site.Cartn_x\n" +
  "ATOM 1 1.5\n" +
  "ATOM 2 2.5\n" +
  "HETATM 3 3.5\n" +
  "_single.multi\n" +
  ";first line\n" +
  "second line\n" +
  ";\n" +
  "_last value2\n"
}

///|
test "mmcif2dict_parse" {
  let d = @src.mmcif2dict_parse(sample_cif()) catch {
    _ => abort("mmCIF sample should parse")
  }
  assert_true(eq_strs(d.get("data_").unwrap_or([]), ["test"]))
  assert_true(eq_strs(d.get("_struct.entry_id").unwrap_or([]), ["1ABC"]))
  assert_true(eq_strs(d.get("_struct.title").unwrap_or([]), ["A \"quoted\" title"]))
  assert_true(eq_strs(d.get("_atom_site.group_PDB").unwrap_or([]), ["ATOM", "ATOM", "HETATM"]))
  assert_true(eq_strs(d.get("_atom_site.id").unwrap_or([]), ["1", "2", "3"]))
  assert_true(eq_strs(d.get("_atom_site.Cartn_x").unwrap_or([]), ["1.5", "2.5", "3.5"]))
  assert_true(eq_strs(d.get("_single.multi").unwrap_or([]), ["first line\nsecond line"]))
  assert_true(eq_strs(d.get("_last").unwrap_or([]), ["value2"]))
}

///|
test "mmcif2dict_empty_file" {
  let r = try {
    ignore(@src.mmcif2dict_parse(""))
    false
  } catch {
    Mmcif2DictError(msg) => msg == "Empty file."
  }
  assert_true(r)
}

///|
test "mmcif2dict_missing_data_directive" {
  let r = try {
    ignore(@src.mmcif2dict_parse("_struct.entry_id 1ABC\n"))
    false
  } catch {
    Mmcif2DictError(msg) => msg == "The input mmCIF file must begin with a 'data_' directive."
  }
  assert_true(r)
}

///|
test "mmcif2dict_missing_closing_semicolon" {
  let r = try {
    ignore(@src.mmcif2dict_parse("data_x\n_key\n;never closed\n"))
    false
  } catch {
    Mmcif2DictError(msg) => msg == "Missing closing semicolon"
  }
  assert_true(r)
}