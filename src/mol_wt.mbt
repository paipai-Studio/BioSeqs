///|
/// mol_wt module - Molecular weight calculation for sequences.
/// Inspired by Biopython Bio.SeqUtils.MolWt.
///
/// References:
///   - https://biopython.org/wiki/SeqUtils
///   - Biopython Bio.SeqUtils.MolWt: Molecular weight calculation for
///     DNA, RNA, and protein sequences.
///
/// This module calculates molecular weights for biological sequences
/// using standard atomic weights for nucleotides and amino acids.

///|
/// Atomic weights for common elements (average isotopic masses).
/// These values are from standard biochemistry references.
pub struct AtomicWeights {
  /// Hydrogen
  h : Double
  /// Carbon
  c : Double
  /// Nitrogen
  n : Double
  /// Oxygen
  o : Double
  /// Phosphorus
  p : Double
  /// Sulfur
  s : Double
}

///|
/// Create default atomic weights.
pub fn AtomicWeights::new() -> AtomicWeights {
  AtomicWeights::{
    h: 1.00794,
    c: 12.0107,
    n: 14.0067,
    o: 15.9994,
    p: 30.97376,
    s: 32.065,
  }
}

///|
/// Default atomic weights instance.
pub fn default_atomic_weights() -> AtomicWeights {
  AtomicWeights::new()
}

///|
/// Molecular weight contribution of each DNA nucleotide.
/// Delegates to the unified table in data.mbt.
pub fn nucleotide_weight(nucleotide : String) -> Double {
  let table = dna_nucleotide_mass()
  table.get(nucleotide.to_upper()).unwrap_or(335.2444)
}

///|
/// Molecular weight contribution of each amino acid (average isotopic).
/// These are full amino acid molecule masses (complete mass, residue + water).
/// Delegates to the unified table in data.mbt, falls back to 118.8860 (average).
pub fn amino_acid_weight(aa : String) -> Double {
  let mass = amino_complete_mass(aa)
  if mass > 0.0 {
    mass
  } else {
    118.8860
  }
}

///|
/// Calculate molecular weight of a DNA sequence.
/// Formula: sum of nucleotide weights - (n-1) * water
/// where n is the length of the sequence.
/// Water weight: 18.01524
pub fn mol_weight_dna(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let mut total = 0.0
  let mut i = 0
  while i < n {
    let base = sequence[i:i + 1].to_owned()
    total = total + nucleotide_weight(base)
    i = i + 1
  }
  // Subtract water for each phosphodiester bond (n-1 bonds)
  total - (n - 1).to_double() * WATER_MASS_AVERAGE
}

///|
/// Calculate molecular weight of an RNA sequence.
/// Same as DNA but with uracil instead of thymine.
pub fn mol_weight_rna(sequence : String) -> Double {
  // RNA uses same weights as DNA but T->U conversion
  mol_weight_dna(sequence.to_upper().replace(old="T", new="U"))
}

///|
/// Calculate molecular weight of a protein sequence.
/// Formula: sum of amino acid weights - (n-1) * water
/// where n is the length of the sequence.
pub fn mol_weight_protein(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let mut total = 0.0
  let mut i = 0
  while i < n {
    let aa = sequence[i:i + 1].to_owned()
    total = total + amino_acid_weight(aa)
    i = i + 1
  }
  // Subtract water for each peptide bond (n-1 bonds)
  total - (n - 1).to_double() * WATER_MASS_AVERAGE
}

///|
/// Calculate the extinction coefficient of a protein at 280 nm.
/// Uses Pace et al. (1995) method:
/// E(280) = nTrp * 5500 + nTyr * 1490 + nCystine * 125
/// where nCystine = nCys / 2
pub fn extinction_coefficient_280(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let mut n_trp = 0
  let mut n_tyr = 0
  let mut n_cys = 0
  let mut i = 0
  while i < n {
    let aa = sequence[i:i + 1]
    match aa.to_upper() {
      "W" => n_trp = n_trp + 1
      "Y" => n_tyr = n_tyr + 1
      "C" => n_cys = n_cys + 1
      _ => ()
    }
    i = i + 1
  }
  let n_cystine = n_cys / 2
  n_trp.to_double() * 5500.0 +
  n_tyr.to_double() * 1490.0 +
  n_cystine.to_double() * 125.0
}

///|
/// Calculate the absorbance coefficient (1 mg/mL) at 280 nm.
/// A(280) = E(280) / molecular_weight * 10
pub fn absorbance_280(sequence : String) -> Double {
  let mw = mol_weight_protein(sequence)
  if mw < 1.0e-10 {
    return 0.0
  }
  let ext = extinction_coefficient_280(sequence)
  ext / mw * 10.0
}

///|
/// Calculate the isoelectric point (pI) using the pK values from
/// Lehninger. This is a simplified version.
/// pI = (pK1 + pK2) / 2 for simple cases
/// For proteins, uses pKa values for charged residues.
/// Named mol_wt_isoelectric_point to avoid conflict with emboss module.
pub fn mol_wt_isoelectric_point(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 7.0
  }
  // Count charged residues
  let mut n_asp = 0 // D
  let mut n_glu = 0 // E
  let mut n_lys = 0 // K
  let mut n_arg = 0 // R
  let mut n_cys = 0 // C
  let mut n_tyr = 0 // Y
  let mut i = 0
  while i < n {
    let aa = sequence[i:i + 1]
    match aa.to_upper() {
      "D" => n_asp = n_asp + 1
      "E" => n_glu = n_glu + 1
      "K" => n_lys = n_lys + 1
      "R" => n_arg = n_arg + 1
      "C" => n_cys = n_cys + 1
      "Y" => n_tyr = n_tyr + 1
      _ => ()
    }
    i = i + 1
  }
  // Simplified pI calculation using Henderson-Hasselbalch
  // pK values (simplified):
  let pK_nterm = 9.69 // N-terminus (alpha-amino)
  let pK_cterm = 2.34 // C-terminus (alpha-carboxyl)
  let pK_D = 3.65 // Asp
  let pK_E = 4.25 // Glu
  let pK_C = 8.18 // Cys
  let pK_Y = 10.07 // Tyr
  let pK_K = 10.54 // Lys
  let pK_R = 12.48 // Arg
  // Simple approximation
  let pos_charges = 1.0 + n_lys.to_double() + n_arg.to_double()
  let neg_charges = 1.0 +
    n_asp.to_double() +
    n_glu.to_double() +
    n_cys.to_double() +
    n_tyr.to_double()
  if pos_charges + neg_charges < 1.0e-10 {
    return 7.0
  }
  // Average pI approximation
  let pI_acidic = (pK_cterm + pK_D + pK_E + pK_C + pK_Y) / 5.0
  let pI_basic = (pK_nterm + pK_K + pK_R) / 3.0
  if pos_charges > neg_charges {
    pI_basic
  } else if neg_charges > pos_charges {
    pI_acidic
  } else {
    (pI_acidic + pI_basic) / 2.0
  }
}

///|
/// Compute the difference in molecular weight between two sequences.
pub fn mol_weight_diff(
  seq1 : String,
  seq2 : String,
  is_protein~ : Bool = false,
) -> Double {
  let mw1 = if is_protein {
    mol_weight_protein(seq1)
  } else {
    mol_weight_dna(seq1)
  }
  let mw2 = if is_protein {
    mol_weight_protein(seq2)
  } else {
    mol_weight_dna(seq2)
  }
  mw1 - mw2
}

///|
/// Generate a summary string for a sequence's molecular properties.
pub fn mol_weight_summary(
  sequence : String,
  is_protein~ : Bool = false,
) -> String {
  let length = sequence.length()
  let mw = if is_protein {
    mol_weight_protein(sequence)
  } else {
    mol_weight_dna(sequence)
  }
  let ext = if is_protein { extinction_coefficient_280(sequence) } else { 0.0 }
  let pI = if is_protein { mol_wt_isoelectric_point(sequence) } else { 0.0 }
  let type_str = if is_protein { "Protein" } else { "DNA" }
  "MolWt Summary:\n" +
  "  Type: " +
  type_str +
  "\n" +
  "  Length: " +
  length.to_string() +
  "\n" +
  "  Molecular Weight: " +
  mw.to_string() +
  " Da\n" +
  (if is_protein {
    "  Extinction Coefficient (280): " +
    ext.to_string() +
    " M^-1 cm^-1\n" +
    "  Isoelectric Point: " +
    pI.to_string()
  } else {
    ""
  })
}

///|
/// Generate sample DNA sequence data for testing.
pub fn mol_wt_sample_dna() -> String {
  "ACGTACGTACGTACGTACGT"
}

///|
/// Generate sample protein sequence data for testing.
pub fn mol_wt_sample_protein() -> String {
  "MVLSPADKTNVKAAWG"
}
