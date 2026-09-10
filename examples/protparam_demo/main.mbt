///|
/// ProtParam demo - Protein parameter analysis example.

///|
fn main raise {
  println("=== ProtParam Protein Analysis Demo ===\n")

  let protein_seq = "MKKLLLISVLLFLSSAYSRGVVVDQQCGGNIFRPEQLVSGSEIHARLGVLGSGGGFRLVAVQ"

  println("Protein sequence: " + protein_seq)
  println("Sequence length: " + protein_seq.length().to_string() + "\n")

  let protein = @bio.ProteinAnalysis::new(protein_seq)

  println("=== Basic Properties ===")
  println("Length: " + protein.length().to_string())
  println("Molecular Weight: " + protein.molecular_weight().to_string())
  println("GRAVY Score: " + protein.gravy().to_string())
  println("Aromaticity: " + protein.aromaticity().to_string())
  println("Instability Index: " + protein.instability_index().to_string())

  println("\n=== Isoelectric Point ===")
  let pi = protein.isoelectric_point()
  println("pI: " + pi.to_string())
  println("Charge at pH 7.0: " + protein.charge_at_ph(7.0).to_string())

  println("\n=== Secondary Structure Propensity ===")
  let (helix, sheet, coil) = protein.secondary_structure()
  println("Helix: " + helix.to_string())
  println("Sheet: " + sheet.to_string())
  println("Coil: " + coil.to_string())

  println("\n=== Amino Acid Percentages ===")
  let percent = protein.get_amino_acids_percent()
  let aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I']
  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    match percent.get(aa) {
      Some(p) => println("  " + aa.to_string() + ": " + p.to_string())
      None => ()
    }
    i = i + 1
  }

  println("\n=== Signal Peptide Prediction ===")
  let sp_result = protein.predict_signal_peptide()
  if sp_result.has_signal_peptide {
    println("Has signal peptide: Yes")
  } else {
    println("Has signal peptide: No")
  }
  println("Score: " + sp_result.score.to_string())
  match sp_result.cleavage_site {
    Some(pos) => println("Cleavage site at position: " + pos.to_string())
    None => println("Cleavage site: Not found")
  }

  println("\n=== Stability Prediction ===")
  let instability = protein.instability_index()
  if instability < 40.0 {
    println(
      "The protein is predicted to be STABLE (instability index: " +
      instability.to_string() +
      ")",
    )
  } else {
    println(
      "The protein is predicted to be UNSTABLE (instability index: " +
      instability.to_string() +
      ")",
    )
  }

  println("\n=== Secondary Structure Fractions ===")
  let (frac_helix, frac_turn, frac_sheet) = protein.secondary_structure_fraction()
  println(
    "Helix (EMALK): " + ((frac_helix * 1000.0).round() / 1000.0).to_string(),
  )
  println(
    "Turn (NPGSD): " + ((frac_turn * 1000.0).round() / 1000.0).to_string(),
  )
  println(
    "Sheet (VIYFWLT): " + ((frac_sheet * 1000.0).round() / 1000.0).to_string(),
  )

  println("\n=== Molar Extinction Coefficient (280 nm) ===")
  let (ext_reduced, ext_cystines) = protein.molar_extinction_coefficient()
  println("Reduced cysteines: " + ext_reduced.to_string())
  println("With cystine bonds: " + ext_cystines.to_string())

  println("\n=== Flexibility Profile (Vihinen 1994, window 9) ===")
  let flex = protein.flexibility()
  println(
    "Profile length: " + flex.length().to_string() + " (sequence length - 9)",
  )
  let mut first_scores = ""
  let mut fi = 0
  while fi < 5 {
    first_scores = first_scores +
      ((flex[fi] * 1000.0).round() / 1000.0).to_string()
    if fi < 4 {
      first_scores = first_scores + " "
    }
    fi = fi + 1
  }
  println("First 5 scores: " + first_scores)

  println("\n=== Hydrophobicity Profiles (protein_scale) ===")
  let kd = @bio.gravy_scale_table("KyteDoolitle")
  let profile9 = protein.protein_scale(kd, 9)
  println(
    "Kyte-Doolittle window 9: " +
    profile9.length().to_string() +
    " positions, mean-scaled first value " +
    ((profile9[0] * 1000.0).round() / 1000.0).to_string(),
  )
  let profile5 = protein.protein_scale(kd, 5, edge=0.4)
  println(
    "Kyte-Doolittle window 5, edge 0.4: " +
    profile5.length().to_string() +
    " positions (tapered window weights)",
  )
  let hw = @bio.hydrophilicity_table()
  let profile_hw = protein.protein_scale(hw, 19)
  println(
    "Hopp-Woods hydrophilicity window 19: " +
    profile_hw.length().to_string() +
    " positions",
  )

  println("\n=== GRAVY with Registered Scales ===")
  println("Registered scales: " + @bio.gravy_scale_names().length().to_string())
  let show_scales = [
    "KyteDoolitle", "Eisenberg", "Engelman", "Parker", "Zimmerman",
  ]
  let mut si = 0
  while si < show_scales.length() {
    let name = show_scales[si]
    let g = protein.gravy_scale(name)
    println("  " + name + ": " + ((g * 1000.0).round() / 1000.0).to_string())
    si = si + 1
  }

  println("\n=== Complete Analysis Done ===")
}
