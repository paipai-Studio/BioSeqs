///|
/// mmCIF format parsing example.
///
/// Run: moon run examples/mmcif_demo/main.mbt
fn main {
  let mmcif_content = "data_1ABC\n#\n_cell_length_a   50.0\n_cell_length_b   50.0\n_cell_length_c   50.0\n_atom_site.id 1\n_atom_site.type_symbol N\n_atom_site.label_atom_id N\n_atom_site.label_comp_id MET\n_atom_site.label_asym_id A\n_atom_site.label_seq_id 1\n_atom_site.Cartn_x 1.0\n_atom_site.Cartn_y 2.0\n_atom_site.Cartn_z 3.0\n"

  println("=== mmCIF Parser Demo ===")
  println("")

  let block = @bio.parse_mmcif(mmcif_content)
  println("Data block name: \{block.name}")

  let atom_sites = block.get_atom_sites()
  println("Atom sites found: \{atom_sites.length()}")

  let entity_data = block.get_entity()
  println("Entity data found: \{entity_data.length()}")

  let struct_data = block.get_struct()
  println("Structure data found: \{struct_data.length()}")

  let atoms = block.to_atoms()
  println("Atoms converted: \{atoms.length()}")

  println("")
  println("=== Done ===")
}
