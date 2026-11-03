///|
/// SASA (Solvent Accessible Surface Area) demo.
///
/// Builds a small sample peptide, calculates per-atom SASA via the
/// Shrake-Rupley rolling-ball algorithm, aggregates per-residue results,
/// and prints an overall summary.
///
/// Run: moon run examples/sasa_demo/main.mbt

///|
fn main {
  let structure = @bio.sasa_sample_structure()

  println("=== Sample Structure ===")
  println("  ID: \{structure.get_id()}")
  println("  Models: \{structure.get_num_models()}")
  println("  Chains: \{structure.get_num_chains()}")
  println("  Residues: \{structure.get_num_residues()}")
  println("  Atoms: \{structure.get_num_atoms()}")

  // Calculate per-atom SASA (100 sphere points, 1.4 Å water probe).
  let atom_results = @bio.sasa_calc(structure, 100, 1.4)

  println("\n=== Per-Atom SASA ===")
  println("  idx  atom  res   sasa(A^2)")
  for r in atom_results {
    println(
      "  \{r.get_atom_index()}\t\{r.get_atom_name()}\t\{r.get_resname()}\t\{r.get_sasa()}",
    )
  }

  // Aggregate per-residue SASA.
  let residue_results = @bio.sasa_calc_residue(atom_results, structure)

  println("\n=== Per-Residue SASA ===")
  println("  res      total    backbone  sidechain  atoms")
  for r in residue_results {
    println(
      "  \{r.get_resname()}\t\{r.get_total_sasa()}\t\{r.get_backbone_sasa()}\t\{r.get_sidechain_sasa()}\t\{r.get_num_atoms()}",
    )
  }

  // Overall summary.
  let summary = @bio.sasa_calc_total(atom_results)

  println("\n=== Total SASA Summary ===")
  println("  Total SASA: \{summary.get_total_sasa()} A^2")
  println("  Backbone : \{summary.get_backbone_sasa()} A^2")
  println("  Sidechain: \{summary.get_sidechain_sasa()} A^2")
  println("  Atoms    : \{summary.get_num_atoms()}")
  println("  Residues : \{summary.get_num_residues()}")
}
