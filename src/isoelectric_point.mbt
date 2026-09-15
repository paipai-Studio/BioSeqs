///|
/// isoelectric_point - Protein isoelectric point (pI) calculation.
/// Inspired by Biopython Bio.SeqUtils.IsoelectricPoint.
/// Uses pKa table and binary search to find pH where net charge = 0.

///|
/// Amino acid pKa values for the terminal groups and side chains.
/// Source: Biopython's IUPAC data (standard pKa table)
struct IonizableGroup {
  pKa : Double
  n_term : Bool
  c_term : Bool
  side_chain : Bool
}

///|
/// pKa values for ionizable groups in proteins.
/// N-terminal: pKa ~ 8.0, C-terminal: pKa ~ 3.1
/// Side chains: standard values from Biopython
fn mt_get_pka(aa : String) -> Array[IonizableGroup] {
  let groups : Array[IonizableGroup] = []

  // N-terminal amino group (always present)
  groups.push(IonizableGroup::{ pKa: 8.0, n_term: true, c_term: false, side_chain: false })

  // C-terminal carboxyl group (always present)
  groups.push(IonizableGroup::{ pKa: 3.1, n_term: false, c_term: true, side_chain: false })

  // Side chain pKa values
  match aa {
    "D" => {
      // Aspartic acid: side chain carboxyl
      groups.push(IonizableGroup::{ pKa: 3.9, n_term: false, c_term: false, side_chain: true })
    }
    "E" => {
      // Glutamic acid: side chain carboxyl
      groups.push(IonizableGroup::{ pKa: 4.3, n_term: false, c_term: false, side_chain: true })
    }
    "H" => {
      // Histidine: imidazole
      groups.push(IonizableGroup::{ pKa: 6.0, n_term: false, c_term: false, side_chain: true })
    }
    "C" => {
      // Cysteine: thiol
      groups.push(IonizableGroup::{ pKa: 8.3, n_term: false, c_term: false, side_chain: true })
    }
    "Y" => {
      // Tyrosine: phenol
      groups.push(IonizableGroup::{ pKa: 10.5, n_term: false, c_term: false, side_chain: true })
    }
    "K" => {
      // Lysine: amine
      groups.push(IonizableGroup::{ pKa: 10.5, n_term: false, c_term: false, side_chain: true })
    }
    "R" => {
      // Arginine: guanidinium
      groups.push(IonizableGroup::{ pKa: 12.0, n_term: false, c_term: false, side_chain: true })
    }
    _ => ()
  }

  groups
}

///|
/// Calculate net charge of a protein at a given pH.
/// Charge = sum of charges from all ionizable groups.
/// For each group: charge = +1 if pH < pKa (protonated, positive)
///                 charge = 0 if pH > pKa (deprotonated, neutral/negative)
/// For carboxyl groups: charge = 0 if pH < pKa (protonated, neutral)
///                      charge = -1 if pH > pKa (deprotonated, negative)
fn mt_charge_at_ph(seq : String, ph : Double) -> Double {
  let seq_upper = seq.to_upper()
  let mut charge = 0.0

  let mut i = 0
  while i < seq_upper.length() {
    let aa = seq_upper[i:i + 1].to_string()
    let groups = mt_get_pka(aa)

    let mut g = 0
    while g < groups.length() {
      let group = groups[g]
      if group.n_term {
        // N-terminal: +1 when protonated (pH < pKa)
        if ph < group.pKa {
          charge = charge + 1.0
        }
      } else if group.c_term {
        // C-terminal: -1 when deprotonated (pH > pKa)
        if ph > group.pKa {
          charge = charge - 1.0
        }
      } else if group.side_chain {
        // Side chains
        match aa {
          "D" | "E" => {
            // Acidic: -1 when deprotonated (pH > pKa)
            if ph > group.pKa {
              charge = charge - 1.0
            }
          }
          "H" => {
            // Histidine: +1 when protonated (pH < pKa)
            if ph < group.pKa {
              charge = charge + 1.0
            }
          }
          "C" => {
            // Cysteine: -1 when deprotonated (pH > pKa)
            if ph > group.pKa {
              charge = charge - 1.0
            }
          }
          "Y" => {
            // Tyrosine: -1 when deprotonated (pH > pKa)
            if ph > group.pKa {
              charge = charge - 1.0
            }
          }
          "K" => {
            // Lysine: +1 when protonated (pH < pKa)
            if ph < group.pKa {
              charge = charge + 1.0
            }
          }
          "R" => {
            // Arginine: +1 when protonated (pH < pKa)
            if ph < group.pKa {
              charge = charge + 1.0
            }
          }
          _ => ()
        }
      }
      g = g + 1
    }
    i = i + 1
  }

  charge
}

///|
/// Calculate isoelectric point (pI) using binary search.
/// pI = pH where net charge = 0.
/// Uses bisection method to find root of charge(pH) = 0.
pub fn mt_isoelectric_point(seq : String) -> Double {
  let seq_upper = seq.to_upper()
  if seq_upper.length() == 0 {
    return 0.0
  }

  // Search range: pH 0 to 14
  let mut low = 0.0
  let mut high = 14.0
  let tolerance = 0.001
  let max_iter = 100

  let mut iter = 0
  while iter < max_iter {
    let mid = (low + high) / 2.0
    let charge = mt_charge_at_ph(seq_upper, mid)

    if charge.abs() < tolerance {
      return mid
    }

    if charge > 0.0 {
      // Positive charge: need higher pH to deprotonate
      low = mid
    } else {
      // Negative charge: need lower pH to protonate
      high = mid
    }

    if high - low < tolerance {
      return mid
    }

    iter = iter + 1
  }

  (low + high) / 2.0
}

///|
/// Calculate pI with a custom pKa table.
pub fn mt_isoelectric_point_custom(
  seq : String,
  pka_table : Map[String, Array[IonizableGroup]]
) -> Double {
  let seq_upper = seq.to_upper()
  if seq_upper.length() == 0 {
    return 0.0
  }

  let mut low = 0.0
  let mut high = 14.0
  let tolerance = 0.001
  let max_iter = 100

  let mut iter = 0
  while iter < max_iter {
    let mid = (low + high) / 2.0
    let mut charge = 0.0

    let mut i = 0
    while i < seq_upper.length() {
      let aa = seq_upper[i:i + 1].to_string()
      match pka_table.get(aa) {
        Some(groups) => {
          let mut g = 0
          while g < groups.length() {
            let group = groups[g]
            if group.n_term {
              if mid < group.pKa { charge = charge + 1.0 }
            } else if group.c_term {
              if mid > group.pKa { charge = charge - 1.0 }
            } else if group.side_chain {
              match aa {
                "D" | "E" => { if mid > group.pKa { charge = charge - 1.0 } }
                "H" => { if mid < group.pKa { charge = charge + 1.0 } }
                "C" | "Y" => { if mid > group.pKa { charge = charge - 1.0 } }
                "K" | "R" => { if mid < group.pKa { charge = charge + 1.0 } }
                _ => ()
              }
            }
            g = g + 1
          }
        }
        None => ()
      }
      i = i + 1
    }

    if charge.abs() < tolerance {
      return mid
    }
    if charge > 0.0 {
      low = mid
    } else {
      high = mid
    }
    if high - low < tolerance {
      return mid
    }
    iter = iter + 1
  }

  (low + high) / 2.0
}

///|
/// Calculate charge at a specific pH.
pub fn mt_charge_at_ph_public(seq : String, ph : Double) -> Double {
  mt_charge_at_ph(seq, ph)
}

///|
/// Example protein sequence.
pub fn mt_example_protein() -> String {
  "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSS"
}

///|
/// Quick pI calculation for a protein.
pub fn mt_protein_pi(seq : String) -> Double {
  mt_isoelectric_point(seq)
}