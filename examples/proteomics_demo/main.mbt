///|
/// Proteomics demo - protein digestion, mass calculation, and isotope pattern analysis.

///|
fn main {
  println("=== Proteomics Demo ===\n")

  let protein_seq = "PEPTIDEKTESTKRPGGLLGRLKDEQMFWYVN"
  println("Protein sequence: " + protein_seq)
  println("Length: " + protein_seq.length().to_string())

  println("\n=== Available Enzymes ===")
  let enzymes = @bio.proteomics_available_enzymes()
  println("Supported enzymes (" + enzymes.length().to_string() + "):")
  let mut i = 0
  while i < enzymes.length() {
    println("  - " + enzymes[i])
    i = i + 1
  }

  println("\n=== Trypsin Digestion ===")
  let trypsin_result = @bio.proteomics_digest(protein_seq, "trypsin")
  println("Number of peptides: " + trypsin_result.num_peptides.to_string())
  println("Protein mass (avg): " + trypsin_result.protein_mass.to_string())
  let mut p = 0
  while p < trypsin_result.peptides.length() {
    let pep = trypsin_result.peptides[p]
    println(
      "  Peptide " +
      (p + 1).to_string() +
      ": " +
      pep.peptide_sequence +
      " (mono: " +
      pep.monoisotopic_mass.to_string() +
      ", avg: " +
      pep.average_mass.to_string() +
      ")"
    )
    p = p + 1
  }

  println("\n=== Trypsin Digestion with Missed Cleavages (1) ===")
  let mc_result = @bio.proteomics_trypsin_digest(protein_seq, 1)
  println("Number of peptides with 1 missed cleavage: " + mc_result.num_peptides.to_string())

  println("\n=== Chymotrypsin Digestion ===")
  let chymo_result = @bio.proteomics_digest(protein_seq, "chymotrypsin")
  println("Number of peptides: " + chymo_result.num_peptides.to_string())
  p = 0
  while p < chymo_result.peptides.length() {
    let pep = chymo_result.peptides[p]
    println("  " + pep.peptide_sequence)
    p = p + 1
  }

  println("\n=== Cysteine Mass Calculation ===")
  let mono_mass = @bio.proteomics_calculate_mass("C", true)
  let avg_mass = @bio.proteomics_calculate_mass("C", false)
  println("Cysteine monoisotopic mass: " + mono_mass.to_string())
  println("Cysteine average mass: " + avg_mass.to_string())

  println("\n=== Isotope Pattern for 'GGR' (charge 2) ===")
  let pattern = @bio.proteomics_isotope_pattern("GGR", 2)
  println("Charge: " + pattern.charge.to_string())
  let mut i2 = 0
  while i2 < pattern.masses.length() {
    println(
      "  Peak " +
      i2.to_string() +
      ": m/z = " +
      pattern.masses[i2].to_string() +
      ", abundance = " +
      pattern.abundances[i2].to_string()
    )
    i2 = i2 + 1
  }

  println("\n=== Fragment Ions for 'PEPT' (charge 1) ===")
  let frags = @bio.proteomics_fragment_ions("PEPT", 1)
  println("Number of fragment ions: " + frags.length().to_string())
  let mut f = 0
  while f < frags.length() {
    let frag = frags[f]
    println(
      "  " +
      frag.peptide_sequence +
      ": mono = " +
      frag.monoisotopic_mass.to_string() +
      ", avg = " +
      frag.average_mass.to_string()
    )
    f = f + 1
  }

  println("\n=== Select Enzyme ===")
  let trypsin = @bio.proteomics_select_enzyme("trypsin")
  match trypsin {
    Some(agent) => {
      println("Enzyme: " + agent.name)
      println("Cleavage specificities:")
      let mut s = 0
      while s < agent.cleavage_sites.length() {
        let site = agent.cleavage_sites[s]
        println("  - " + site.residue)
        s = s + 1
      }
    }
    None => println("Enzyme not found")
  }

  println("\n=== Complete Proteomics Demo ===")
}