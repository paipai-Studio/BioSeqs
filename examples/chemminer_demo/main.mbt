///|
/// ChemmineR module demo - Cheminformatics toolkit for small-molecule analysis.
///
/// Run: moon run examples/chemminer_demo
fn main {
  println("=== ChemmineR Module Demo ===")
  println("Cheminformatics: SDF parsing, molecular descriptors, fingerprints, similarity")
  println("Reference: Backman et al. (2011) Bioinformatics 27:2400-1")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample molecules and compute descriptors
  // ------------------------------------------------------------------
  println("=== 1. Molecular Descriptors ===")
  let molecules = [
    @bio.chem_sample_water(),
    @bio.chem_sample_methane(),
    @bio.chem_sample_ethanol(),
    @bio.chem_sample_benzene(),
    @bio.chem_sample_aspirin(),
  ]
  for mol in molecules {
    let summary = @bio.chem_molecule_summary(mol)
    println(summary)
    println("")
  }

  // ------------------------------------------------------------------
  // 2. Atom-pair fingerprints and similarity
  // ------------------------------------------------------------------
  println("=== 2. Atom-Pair Fingerprints & Similarity ===")
  let methane = @bio.chem_sample_methane()
  let ethanol = @bio.chem_sample_ethanol()
  let benzene = @bio.chem_sample_benzene()
  let aspirin = @bio.chem_sample_aspirin()
  let fp_methane = @bio.chem_atom_pair_fingerprint(methane)
  let fp_ethanol = @bio.chem_atom_pair_fingerprint(ethanol)
  let fp_benzene = @bio.chem_atom_pair_fingerprint(benzene)
  let fp_aspirin = @bio.chem_atom_pair_fingerprint(aspirin)
  println("Methane fingerprint size: " + fp_methane.length().to_string())
  println("Ethanol fingerprint size: " + fp_ethanol.length().to_string())
  println("Benzene fingerprint size: " + fp_benzene.length().to_string())
  println("Aspirin fingerprint size: " + fp_aspirin.length().to_string())
  println("")
  let tanimoto_bb = @bio.chem_tanimoto(fp_benzene, fp_aspirin)
  let dice_bb = @bio.chem_dice(fp_benzene, fp_aspirin)
  println("Benzene vs Aspirin Tanimoto: " + tanimoto_bb.to_string())
  println("Benzene vs Aspirin Dice: " + dice_bb.to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Substructure search
  // ------------------------------------------------------------------
  println("=== 3. Substructure Search ===")
  let is_benzene_in_aspirin = @bio.chem_substructure_match(benzene, aspirin)
  println("Benzene is substructure of Aspirin: " + is_benzene_in_aspirin.to_string())
  let is_water_in_ethanol = @bio.chem_substructure_match(
    @bio.chem_sample_water(),
    ethanol,
  )
  println("Water is substructure of Ethanol: " + is_water_in_ethanol.to_string())
  println("")

  // ------------------------------------------------------------------
  // 4. SDF parsing and round-trip
  // ------------------------------------------------------------------
  println("=== 4. SDF Parsing & Round-Trip ===")
  let sdf_text = @bio.chem_sample_sdf_text()
  let parsed = @bio.chem_parse_sdf(sdf_text)
  println("Parsed " + parsed.length().to_string() + " molecule(s) from SDF")
  for mol in parsed {
    println("  " + mol.name() + ": " + mol.n_atoms().to_string() + " atoms, " + mol.n_bonds().to_string() + " bonds")
  }
  println("")
  let reserialized = @bio.chem_molecule_to_sdf(parsed[0])
  let reparsed = @bio.chem_parse_sdf(reserialized)
  println("Round-trip: parsed " + reparsed.length().to_string() + " molecule(s)")
  println("  Name: " + reparsed[0].name())
  println("  Atoms: " + reparsed[0].n_atoms().to_string())
  println("  Bonds: " + reparsed[0].n_bonds().to_string())
  println("")

  println("=== Demo complete ===")
}
