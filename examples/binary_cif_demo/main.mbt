///|
fn main {
  println("=== Biopython Bio.PDB.binary_cif Demo ===")

  let bytes = @src.binary_cif_sample_bytes()
  let file = @src.binary_cif_parse(bytes) catch {
    _ => abort("failed to parse BinaryCIF sample")
  }

  println("\n1. MessagePack BinaryCIF document")
  println("  Encoded bytes: " + bytes.length().to_string())
  println("  Encoder: " + file.encoder)
  println("  " + file.summary())

  let block = file.get_block(0).unwrap()
  let atom_site = block.get_category("atom_site").unwrap()
  println("\n2. Category and column access")
  println("  Block header: " + block.header)
  println("  _atom_site rows: " + atom_site.row_count.to_string())
  let first_component = atom_site
    .get_column("label_comp_id")
    .unwrap()
    .string_at(0) catch {
    _ => abort("failed to read component column")
  }
  println(
    "  First component: " + first_component.unwrap(),
  )
  let first_x = atom_site
    .get_column("Cartn_x")
    .unwrap()
    .double_at(0) catch {
    _ => abort("failed to read coordinate column")
  }
  println(
    "  First x coordinate: " + first_x.unwrap().to_string(),
  )

  println("\n3. CIF missing-value masks")
  let altloc = atom_site.get_column("label_alt_id").unwrap()
  let token0 = altloc.cif_token_at(0) catch {
    _ => abort("failed to read mask row 0")
  }
  let token3 = altloc.cif_token_at(3) catch {
    _ => abort("failed to read mask row 3")
  }
  let token6 = altloc.cif_token_at(6) catch {
    _ => abort("failed to read mask row 6")
  }
  println("  Row 0 token: " + token0)
  println("  Row 3 token: " + token3)
  println("  Row 6 token: " + token6)

  println("\n4. PDB structure construction")
  let structure = @src.binary_cif_to_structure(file) catch {
    _ => abort("failed to construct PDB structure")
  }
  println("  Structure ID: " + structure.get_id())
  println("  Models: " + structure.get_num_models().to_string())
  println("  Chains: " + structure.get_num_chains().to_string())
  println("  Residues: " + structure.get_num_residues().to_string())
  println("  Atoms: " + structure.get_num_atoms().to_string())
  for model in structure.get_models() {
    println(
      "  Model " +
      model.get_id().to_string() +
      ": " +
      model.get_num_chains().to_string() +
      " chain(s), " +
      model.get_atoms().length().to_string() +
      " atom(s)",
    )
  }

  println("\n5. Atom annotations")
  let atoms = structure.get_atoms()
  let first = atoms[0].get_coord()
  println(
    "  First atom " +
    atoms[0].name +
    ": (" +
    first.x.to_string() +
    ", " +
    first.y.to_string() +
    ", " +
    first.z.to_string() +
    ")",
  )
  println(
    "  Alternate location on atom 4: " + atoms[3].altloc.to_string(),
  )
  println("  Water residue is hetero: " + structure.get_residues()[2].is_het().to_string())

  println("\nNotes")
  println("  Gzip-compressed .bcif input must be decompressed before parsing.")
  println("  PDB Structure chain IDs are one character; category data keeps full IDs.")
  println("\n=== Demo Complete ===")
}
