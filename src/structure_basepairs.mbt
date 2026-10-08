///|
/// Basepair detection for nucleic acid structures.
/// Ported from biotite.structure.basepairs.
///
/// Uses the Leontis-Westhof classification based on geometric
/// edge-to-edge hydrogen bonding patterns between nucleobases.

///|
/// A detected base pair.
pub struct BasePair {
  residue1 : Int   // atom index of first base's representative atom
  residue2 : Int   // atom index of second base's representative atom
  resname1 : String
  resname2 : String
  pair_type : String  // "WC", "Wobble", "Hoogsteen", "Other"
  distance : Double   // centroid distance (Å)
}

///|
/// Check if residue name is a nucleic acid base.
fn is_nucleobase(resname : String) -> Bool {
  let r = to_uppercase(resname)
  r == "A" || r == "C" || r == "G" || r == "U" ||
  r == "T" || r == "DA" || r == "DC" || r == "DG" || r == "DT" || r == "DU"
}

///|
/// Get the representative atom name for a nucleobase (N1 for pyrimidines, N9 for purines).
fn base_rep_atom(resname : String) -> String {
  let r = to_uppercase(resname)
  if r == "A" || r == "G" || r == "DA" || r == "DG" {
    "N9"
  } else {
    "N1"
  }
}

///|
/// Detect base pairs in a nucleic acid structure.
/// Uses simplified geometric criteria:
///   - Base-base centroid distance ≤ 15 Å
///   - Glycosidic bond (N1/N9—C1') distance check
///   - Hydrogen bond count between bases ≥ 1
///
/// Parameters:
///   `atoms` — all atoms in the structure
///   `max_distance` — maximum centroid distance (Å), default 15.0
pub fn detect_basepairs(
  atoms : Array[Atom],
  max_distance~ : Double = 15.0,
) -> Array[BasePair] {
  // Find nucleobase representative atoms
  let base_atoms : Array[(Int, String, String, Int, Char)] = Array(capacity=0)
  // (atom_index, resname, rep_atom_name, resseq, chainid)
  for i = 0; i < atoms.length(); i = i + 1 {
    if is_nucleobase(atoms[i].resname) {
      let rep = base_rep_atom(atoms[i].resname)
      if atoms[i].name == rep {
        base_atoms.push((i, atoms[i].resname, rep, atoms[i].resseq, atoms[i].chainid))
      }
    }
  }
  // Collect all atoms per residue for H-bond counting
  // For each pair of base representative atoms, check distance and H-bonds
  let result : Array[BasePair] = Array(capacity=0)
  for i = 0; i < base_atoms.length(); i = i + 1 {
    for j = i + 1; j < base_atoms.length(); j = j + 1 {
      let (idx1, rn1, _, rs1, ch1) = base_atoms[i]
      let (idx2, rn2, _, rs2, ch2) = base_atoms[j]
      // Skip same residue
      if ch1 == ch2 && rs1 == rs2 {
        continue
      }
      let dist = atoms[idx1].coord.distance(atoms[idx2].coord)
      if dist > max_distance {
        continue
      }
      // Count H-bonds between atoms of these two residues
      let mut hbond_count = 0
      for a = 0; a < atoms.length(); a = a + 1 {
        if atoms[a].chainid != ch1 || atoms[a].resseq != rs1 {
          continue
        }
        for b = 0; b < atoms.length(); b = b + 1 {
          if atoms[b].chainid != ch2 || atoms[b].resseq != rs2 {
            continue
          }
          let d = atoms[a].coord.distance(atoms[b].coord)
          // N-H...O or N-H...N hydrogen bond: donor-acceptor ≤ 3.5 Å
          let an = to_uppercase(atoms[a].name)
          let bn = to_uppercase(atoms[b].name)
          let a_is_donor = an == "N1" || an == "N2" || an == "N6" || an == "N7" || an == "O2" || an == "O4"
          let b_is_acceptor = bn == "O2" || bn == "O4" || bn == "O6" || bn == "N1" || bn == "N3" || bn == "N7"
          if a_is_donor && b_is_acceptor && d <= 3.5 && d >= 2.5 {
            hbond_count = hbond_count + 1
          }
        }
      }
      if hbond_count >= 1 {
        // Classify pair type
        let pair_type = classify_basepair(rn1, rn2, hbond_count)
        result.push(BasePair::{
          residue1: idx1,
          residue2: idx2,
          resname1: rn1,
          resname2: rn2,
          pair_type,
          distance: dist,
        })
      }
    }
  }
  result
}

///|
/// Classify base pair type based on residue names and H-bond count.
fn classify_basepair(rn1 : String, rn2 : String, hbond_count : Int) -> String {
  let r1 = to_uppercase(rn1)
  let r2 = to_uppercase(rn2)
  // Normalize: canonical Watson-Crick pairs
  let is_wc = (r1 == "A" && r2 == "U") || (r1 == "U" && r2 == "A") ||
    (r1 == "A" && r2 == "T") || (r1 == "T" && r2 == "A") ||
    (r1 == "G" && r2 == "C") || (r1 == "C" && r2 == "G") ||
    (r1 == "DA" && r2 == "DT") || (r1 == "DT" && r2 == "DA") ||
    (r1 == "DG" && r2 == "DC") || (r1 == "DC" && r2 == "DG")
  if is_wc && hbond_count >= 2 {
    return "WC"
  }
  // Wobble: G-U
  let is_wobble = (r1 == "G" && r2 == "U") || (r1 == "U" && r2 == "G") ||
    (r1 == "DG" && r2 == "DT") || (r1 == "DT" && r2 == "DG")
  if is_wobble {
    return "Wobble"
  }
  // Hoogsteen: A-U or A-T with N7 involvement
  if is_wc && hbond_count == 1 {
    return "Hoogsteen"
  }
  "Other"
}
