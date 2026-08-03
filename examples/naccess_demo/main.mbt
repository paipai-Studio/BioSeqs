///|
/// NACCESS (Bio.PDB.NACCESS) demo.
///
/// Loads sample NACCESS output (RSA residue-level and ASA atom-level data),
/// inspects per-residue and per-atom accessibility, classifies residues as
/// exposed vs buried, reports chain totals, and demonstrates round-trip
/// serialization back to RSA/ASA file content.

///|
fn main {
  println("=== Bio.PDB.NACCESS Demo ===")
  println("Inspired by Biopython Bio.PDB.NACCESS module")
  println("")

  // 1. Load sample NACCESS data (RSA + ASA).
  println("1. Load sample NACCESS data")
  let rsa_content = @bio.naccess_sample_rsa()
  let asa_content = @bio.naccess_sample_asa()
  let result = @bio.naccess_parse(rsa_content, asa_content)
  println("  Parsed RSA and ASA content into a NaccessResult")
  println("")

  // 2. Print result info.
  println("2. Result Info")
  println("  Residues: \{result.get_num_residues()}")
  println("  Atoms   : \{result.get_num_atoms()}")
  println("  Total ASA: \{result.total_abs_asa} A^2")
  println("  Chains  : \{result.chain_totals.length()}")
  println("")

  // 3. Show residue details.
  println("3. Residue Details")
  println("  res    num  chn  abs_all  rel_all  abs_side  rel_side")
  for r in result.get_residues() {
    println(
      "  \{r.res_name}\t\{r.res_num}\t\{r.chain_id}\t\{r.abs_all}\t\{r.rel_all}\t\{r.abs_side_chain}\t\{r.rel_side_chain}",
    )
  }
  println("")

  // 4. Show atom details.
  println("4. Atom Details (ASA)")
  println("  num  atom  res    chn  abs_asa  rel_asa")
  for a in result.get_atoms() {
    println(
      "  \{a.atom_num}\t\{a.atom_name}\t\{a.res_name}\t\{a.chain_id}\t\{a.abs_asa}\t\{a.rel_asa}",
    )
  }
  println("")

  // 5. Show exposed vs buried residues (threshold = 20% relative ASA).
  println("5. Exposed vs Buried Residues (threshold = 20%)")
  let exposed = @bio.naccess_exposed_residues(result, 20.0)
  let buried = @bio.naccess_buried_residues(result, 20.0)
  println("  Exposed (\{exposed.length()}):")
  for r in exposed {
    println("    \{r.res_name} \{r.res_num} (rel_all=\{r.rel_all}%)")
  }
  println("  Buried (\{buried.length()}):")
  for r in buried {
    println("    \{r.res_name} \{r.res_num} (rel_all=\{r.rel_all}%)")
  }
  println("  Count exposed: \{@bio.naccess_count_exposed(result, 20.0)}")
  println("  Mean rel ASA : \{@bio.naccess_mean_rel_asa(result)}%")
  println("")

  // 6. Show chain totals.
  println("6. Chain Totals")
  println("  Chain A: \{@bio.naccess_chain_total_asa(result, "A")} A^2")
  println("  Chain B: \{@bio.naccess_chain_total_asa(result, "B")} A^2")
  println("  Chain A residues: \{@bio.naccess_get_chain_residues(result, "A").length()}")
  println("  Chain B residues: \{@bio.naccess_get_chain_residues(result, "B").length()}")
  println("")

  // 7. Demonstrate serialization (write RSA/ASA).
  println("7. Serialization (round-trip)")
  let rsa_out = @bio.naccess_write_rsa(result)
  let asa_out = @bio.naccess_write_asa(result)
  // Show first few lines of serialized RSA.
  let rsa_lines = rsa_out.split("\n").to_array()
  println("  RSA output (\{rsa_lines.length()} lines, first 3):")
  let mut i = 0
  while i < 3 && i < rsa_lines.length() {
    println("    \{rsa_lines[i]}")
    i = i + 1
  }
  let asa_lines = asa_out.split("\n").to_array()
  println("  ASA output (\{asa_lines.length()} lines, first 3):")
  let mut j = 0
  while j < 3 && j < asa_lines.length() {
    println("    \{asa_lines[j]}")
    j = j + 1
  }
  // Verify round-trip: re-parse and check counts match.
  let reparsed = @bio.naccess_parse(rsa_out, asa_out)
  println("  Round-trip: \{reparsed.get_num_residues()} residues, \{reparsed.get_num_atoms()} atoms")
  println("")

  // 8. Standard ASA reference values.
  println("8. Standard ASA Reference Values")
  let std_res = ["ALA", "GLY", "MET", "CYS", "TRP"]
  for name in std_res {
    println("  \{name}: \{@bio.naccess_standard_asa(name)} A^2")
  }
  println("")

  // 9. Print summary.
  println("9. Summary")
  let summary = @bio.naccess_summary(result)
  println(summary)
}
