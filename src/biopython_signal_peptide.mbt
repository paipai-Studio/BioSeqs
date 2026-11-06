/// biopython_signal_peptide.mbt — Signal peptide prediction
/// Inspired by Biopython SignalP-like heuristics

///|
/// Signal peptide prediction result.
pub struct BiopythonSignalPeptide {
  is_signal : Bool
  cleavage_site : Int32 // position of cleavage site (0-based)
  n_region : Int32 // length of n-region (positively charged)
  h_region : Int32 // length of h-region (hydrophobic core)
  c_region : Int32 // length of c-region
  score : Double // confidence score 0-1
}

///|
/// Hydrophobicity values (Kyte-Doolittle scale).
fn biopython_sp_hydrophobicity(aa : String) -> Double {
  match aa {
    "A" => 1.8
    "R" => -4.5
    "N" => -3.5
    "D" => -3.5
    "C" => 2.5
    "Q" => -3.5
    "E" => -3.5
    "G" => -0.4
    "H" => -3.2
    "I" => 4.5
    "L" => 3.8
    "K" => -3.9
    "M" => 1.9
    "F" => 2.8
    "P" => -1.6
    "S" => -0.8
    "T" => -0.7
    "W" => -0.9
    "Y" => -1.3
    "V" => 4.2
    _ => 0.0
  }
}

///|
/// Predict signal peptide using a simplified heuristic.
/// Signal peptides: 1) positive n-region near N-terminus
/// 2) hydrophobic h-region 3) polar c-region with cleavage site
pub fn biopython_predict_signal_peptide(
  sequence : String,
  max_length~ : Int32 = 70, // max signal peptide length
) -> BiopythonSignalPeptide {
  let n = sequence.length()
  let max_check = if n < max_length { n } else { max_length }
  // Check for basic signal peptide criteria:
  // 1. Must start within first few residues
  // 2. Hydrophobic stretch of ~7-15 residues after n-region
  // 3. Cleavage site typically at position 18-30

  // Calculate hydrophobicity profile
  let _window = 7
  let mut best_score = 0.0
  let mut best_cleavage = 0
  let mut best_h_start = 0
  let mut best_h_end = 0

  // Scan for hydrophobic core (h-region)
  for h_start = 1; h_start < 15 && h_start < max_check; h_start = h_start + 1 {
    for h_len = 7; h_len <= 15; h_len = h_len + 1 {
      let h_end = h_start + h_len
      if h_end > max_check {
        break
      }
      let mut h_score = 0.0
      for i = h_start; i < h_end; i = i + 1 {
        h_score = h_score +
          biopython_sp_hydrophobicity(sequence[i:i + 1].to_owned())
      }
      h_score = h_score / h_len.to_double()
      // h-region should be hydrophobic (score > 1.0)
      if h_score > 1.0 {
        // Look for cleavage site after h-region
        for c_site = h_end
            c_site < h_end + 13 && c_site < max_check
            c_site = c_site + 1 {
          // Cleavage site preference: small residues at -1 and +1
          let pos_m1 = if c_site > 0 {
            sequence[c_site - 1:c_site].to_owned()
          } else {
            "X"
          }
          let pos_1 = if c_site < n - 1 {
            sequence[c_site:c_site + 1].to_owned()
          } else {
            "X"
          }
          let cleavage_score = biopython_sp_cleavage_preference(pos_m1, pos_1)
          let total_score = h_score * 0.5 + cleavage_score * 0.5
          if total_score > best_score {
            best_score = total_score
            best_cleavage = c_site
            best_h_start = h_start
            best_h_end = h_end
          }
        }
      }
    }
  }

  // Determine n-region (from start to h-region)
  let n_region = best_h_start
  let h_region = best_h_end - best_h_start
  let c_region = best_cleavage - best_h_end

  // Normalize score
  let normalized_score = if best_score > 1.0 { 1.0 } else { best_score }
  let is_signal = best_score > 0.5 && best_cleavage > 10 && best_cleavage < 40

  {
    is_signal,
    cleavage_site: best_cleavage,
    n_region,
    h_region,
    c_region,
    score: normalized_score,
  }
}

///|
/// Cleavage site preference: Ala, Gly, Ser, Cys common at -1 and +1.
fn biopython_sp_cleavage_preference(pos_m1 : String, pos_1 : String) -> Double {
  let score_m1 = match pos_m1 {
    "A" => 1.0
    "G" => 0.8
    "S" => 0.7
    "C" => 0.6
    "T" => 0.5
    "V" => 0.4
    _ => 0.2
  }
  let score_1 = match pos_1 {
    "A" => 1.0
    "G" => 0.8
    "S" => 0.7
    "E" => 0.6
    "K" => 0.5
    "D" => 0.4
    _ => 0.3
  }
  (score_m1 + score_1) / 2.0
}

///|
/// Extract the mature protein (after signal peptide cleavage).
pub fn biopython_mature_protein(
  sequence : String,
  cleavage_site : Int32,
) -> String {
  if cleavage_site >= sequence.length() {
    ""
  } else {
    sequence[cleavage_site:].to_owned()
  }
}
