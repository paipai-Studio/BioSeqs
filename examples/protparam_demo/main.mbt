///|
/// ProtParam demo - Protein parameter analysis example.

///|
fn main {
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
    println("The protein is predicted to be STABLE (instability index: " + instability.to_string() + ")")
  } else {
    println("The protein is predicted to be UNSTABLE (instability index: " + instability.to_string() + ")")
  }
  
  println("\n=== Complete Analysis Done ===")
}