///|
/// isoelectric_point demo - Protein isoelectric point (pI) calculation.
/// Inspired by Biopython Bio.SeqUtils.IsoelectricPoint.

fn main {
  println("=== isoelectric_point Demo: Protein pI Calculation ===")
  println("")

  // Step 1: Single amino acids
  println("--- Step 1: Single Amino Acids ---")
  let a_pi = @src.mt_isoelectric_point("A")
  println("Alanine (A): pI = " + a_pi.to_string())
  let d_pi = @src.mt_isoelectric_point("D")
  println("Aspartic acid (D): pI = " + d_pi.to_string())
  let k_pi = @src.mt_isoelectric_point("K")
  println("Lysine (K): pI = " + k_pi.to_string())
  println("")

  // Step 2: Charge at different pH
  println("--- Step 2: Charge vs pH ---")
  let seq = "DEKR"
  let ph_values = [2.0, 4.0, 6.0, 7.0, 8.0, 10.0, 12.0]
  let mut i = 0
  while i < ph_values.length() {
    let charge = @src.mt_charge_at_ph_public(seq, ph_values[i])
    println("  pH " + ph_values[i].to_string() + ": charge = " + charge.to_string())
    i = i + 1
  }
  println("")

  // Step 3: Acidic vs basic proteins
  println("--- Step 3: Acidic vs Basic Proteins ---")
  let acidic = "DDD"
  let basic = "KKK"
  let mixed = "DEKR"
  println("Acidic (DDD): pI = " + @src.mt_isoelectric_point(acidic).to_string())
  println("Basic (KKK): pI = " + @src.mt_isoelectric_point(basic).to_string())
  println("Mixed (DEKR): pI = " + @src.mt_isoelectric_point(mixed).to_string())
  println("")

  // Step 4: Example protein
  println("--- Step 4: Example Protein ---")
  let example = @src.mt_example_protein()
  println("Sequence: " + example)
  let pi = @src.mt_isoelectric_point(example)
  println("pI = " + pi.to_string())
  println("")

  // Step 5: Histidine special case
  println("--- Step 5: Histidine (pKa ~ 6.0) ---")
  let h_seq = "H"
  let h_charge_low = @src.mt_charge_at_ph_public(h_seq, 5.0)
  let h_charge_high = @src.mt_charge_at_ph_public(h_seq, 7.0)
  println("Charge at pH 5.0: " + h_charge_low.to_string() + " (protonated)")
  println("Charge at pH 7.0: " + h_charge_high.to_string() + " (deprotonated)")
  println("")

  println("=== isoelectric_point Demo Complete ===")
}