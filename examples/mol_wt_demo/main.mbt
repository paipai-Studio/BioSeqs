///|
/// Example: mol_wt_demo
/// Demonstrates the Bio.SeqUtils.MolWt module (molecular weight calculation).

fn main {
  println("=== Molecular Weight Calculator Demo ===")
  println("")
  
  // Example 1: DNA molecular weight
  println("=== DNA Molecular Weight ===")
  let dna_seq = "ACGTACGTACGT"
  let dna_len = dna_seq.length()
  let dna_mw = @bio.mol_weight_dna(dna_seq)
  println("Sequence: \{dna_seq}")
  println("Length: \{dna_len} bp")
  println("Molecular Weight: \{dna_mw} Da")
  println("")
  
  // Example 2: RNA molecular weight
  println("=== RNA Molecular Weight ===")
  let rna_seq = "ACGUACGUACGU"
  let rna_mw = @bio.mol_weight_rna(rna_seq)
  println("Sequence: \{rna_seq}")
  println("Molecular Weight: \{rna_mw} Da")
  println("")
  
  // Example 3: Protein molecular weight
  println("=== Protein Molecular Weight ===")
  let protein_seq = "MVLSPADKTNVKAAWGCK"
  let protein_len = protein_seq.length()
  let protein_mw = @bio.mol_weight_protein(protein_seq)
  println("Sequence: \{protein_seq}")
  println("Length: \{protein_len} aa")
  println("Molecular Weight: \{protein_mw} Da")
  println("")
  
  // Example 4: Extinction coefficient at 280nm
  println("=== Extinction Coefficient (280nm) ===")
  let ext = @bio.extinction_coefficient_280(protein_seq)
  println("Protein: \{protein_seq}")
  println("Extinction Coefficient: \{ext} M^-1 cm^-1")
  println("")
  
  // Example 5: Absorbance coefficient (1 mg/mL)
  println("=== Absorbance Coefficient ===")
  let abs = @bio.absorbance_280(protein_seq)
  println("Absorbance (1 mg/mL): \{abs}")
  println("")
  
  // Example 6: Isoelectric point
  println("=== Isoelectric Point ===")
  let pI = @bio.mol_wt_isoelectric_point(protein_seq)
  println("Protein: \{protein_seq}")
  println("pI: \{pI}")
  
  // Basic protein
  let basic = "KKKKKK"
  let pI_basic = @bio.mol_wt_isoelectric_point(basic)
  println("Basic (\{basic}): pI = \{pI_basic}")
  
  // Acidic protein
  let acidic = "DEDEDE"
  let pI_acidic = @bio.mol_wt_isoelectric_point(acidic)
  println("Acidic (\{acidic}): pI = \{pI_acidic}")
  println("")
  
  // Example 7: Molecular weight difference
  println("=== Molecular Weight Difference ===")
  let diff = @bio.mol_weight_diff("ACGT", "TGCA")
  println("Diff (ACGT vs TGCA): \{diff} Da")
  println("")
  
  // Example 8: Summary
  println("=== Full Summary ===")
  println(@bio.mol_weight_summary(protein_seq, is_protein=true))
  println("")
  println(@bio.mol_weight_summary(dna_seq))
}
