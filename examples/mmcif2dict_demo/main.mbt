///|
/// Bio.PDB.MMCIF2Dict Demo.
///
/// Demonstrates the MMCIF2Dict port: tokenizing a small mmCIF document
/// (quoted strings, a `loop_` block and a `;`-multi-line value) into a flat
/// token -> values dictionary.

///|
fn main raise {
  let cif =
    "data_1abc\n" +
    "# a comment line\n" +
    "_struct.entry_id 1ABC\n" +
    "_struct.title 'A small structure'\n" +
    "loop_\n" +
    "_atom_site.group_PDB\n" +
    "_atom_site.id\n" +
    "_atom_site.Cartn_x\n" +
    "ATOM 1 1.5\n" +
    "ATOM 2 2.5\n" +
    "HETATM 3 3.5\n" +
    "_note.text\n" +
    ";a multi-line\n" +
    "value here\n" +
    ";\n"
  let d = @src.mmcif2dict_parse(cif)
  println("data block = \{d.get("data_").unwrap_or([])[0]}")
  for key in ["_struct.entry_id", "_struct.title", "_atom_site.id", "_atom_site.Cartn_x", "_note.text"] {
    match d.get(key) {
      Some(vals) => println("\{key} -> \{vals.length()} value(s)")
      None => println("\{key} -> (absent)")
    }
  }
  println("_note.text = \{d.get("_note.text").unwrap_or([])[0]}")
}