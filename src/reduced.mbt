///|
/// reduced module - Reduced amino acid alphabets for simplified alignments.
/// Inspired by Biopython Bio.Align.Reduced module.
///
/// References:
///   - https://biopython.org/wiki/Align
///   - Biopython Bio.Align.Reduced: Reduced amino acid alphabets (RAD, Dayhoff,
///     CHARM, SDM12, etc.) that group similar amino acids together for
///     simplified sequence comparison and alignment.
///
/// Reduced alphabets are useful for detecting remote homology,
/// simplifying sequence comparison, and reducing computational complexity
/// in multiple sequence alignment.

///|
/// A reduced amino acid alphabet with its grouping scheme.
pub struct ReducedAlphabet {
  /// Name of the reduced alphabet
  name : String
  /// Number of groups
  n_groups : Int
  /// Mapping from single-letter amino acid code to group letter
  mapping : Map[String, String]
  /// Group definitions (which amino acids are in each group)
  groups : Array[String]
}

///|
/// Create a new reduced alphabet.
pub fn ReducedAlphabet::new(
  name : String,
  mapping : Map[String, String],
  groups : Array[String],
) -> ReducedAlphabet {
  ReducedAlphabet::{
    name,
    n_groups: groups.length(),
    mapping,
    groups,
  }
}

///|
/// Get the alphabet name.
pub fn ReducedAlphabet::get_name(self : ReducedAlphabet) -> String {
  self.name
}

///|
/// Get the number of groups.
pub fn ReducedAlphabet::get_n_groups(self : ReducedAlphabet) -> Int {
  self.n_groups
}

///|
/// Get the reduced representation for an amino acid.
pub fn ReducedAlphabet::reduce_aa(self : ReducedAlphabet, aa : String) -> String {
  let upper = aa.to_upper()
  self.mapping.get(upper).unwrap_or("X")
}

///|
/// Reduce a protein sequence to the reduced alphabet.
pub fn ReducedAlphabet::reduce_sequence(self : ReducedAlphabet, sequence : String) -> String {
  let mut result = "".to_string()
  let n = sequence.length()
  let mut i = 0
  while i < n {
    let aa = sequence.substring(start=i, end=i + 1)
    result = result + self.reduce_aa(aa)
    i = i + 1
  }
  result
}

///|
/// Get the group definition.
pub fn ReducedAlphabet::get_groups(self : ReducedAlphabet) -> Array[String] {
  self.groups.copy()
}

///|
/// Generate a summary string.
pub fn ReducedAlphabet::summary(self : ReducedAlphabet) -> String {
  let mut s = "ReducedAlphabet: " + self.name + "\n"
  s = s + "  Groups: " + self.n_groups.to_string() + "\n"
  let mut i = 0
  while i < self.groups.length() {
    s = s + "  Group " + i.to_string() + ": " + self.groups[i] + "\n"
    i = i + 1
  }
  s
}

///|
/// RAD (Reduced Alphabet Database) - 6 groups.
/// From Wang & Wang (1999) "How similar are protein sequences?"
/// Groups: (A,G), (C), (D,E,N,Q), (I,L,M,V), (F,Y,W), (H,K,R,S,T)
pub fn rad_alphabet() -> ReducedAlphabet {
  let mapping : Map[String, String] = Map([], capacity=20)
  // Group 0: (A,G)
  mapping.set("A", "0")
  mapping.set("G", "0")
  // Group 1: (C)
  mapping.set("C", "1")
  // Group 2: (D,E,N,Q)
  mapping.set("D", "2")
  mapping.set("E", "2")
  mapping.set("N", "2")
  mapping.set("Q", "2")
  // Group 3: (I,L,M,V)
  mapping.set("I", "3")
  mapping.set("L", "3")
  mapping.set("M", "3")
  mapping.set("V", "3")
  // Group 4: (F,Y,W)
  mapping.set("F", "4")
  mapping.set("Y", "4")
  mapping.set("W", "4")
  // Group 5: (H,K,R,S,T)
  mapping.set("H", "5")
  mapping.set("K", "5")
  mapping.set("R", "5")
  mapping.set("S", "5")
  mapping.set("T", "5")
  let groups = ["(A,G)", "(C)", "(D,E,N,Q)", "(I,L,M,V)", "(F,Y,W)", "(H,K,R,S,T)"]
  ReducedAlphabet::new("RAD", mapping, groups)
}

///|
/// Dayhoff - 6 groups based on chemical similarity.
/// From Dayhoff et al. (1978) "Atlas of Protein Sequence and Structure"
/// Groups: (A,G,P,S,T), (C), (D,E,N,Q), (I,L,M,V), (F,W,Y), (H,K,R)
pub fn dayhoff_alphabet() -> ReducedAlphabet {
  let mapping : Map[String, String] = Map([], capacity=20)
  // Group 0: (A,G,P,S,T) - small/hydroxyl
  mapping.set("A", "0")
  mapping.set("G", "0")
  mapping.set("P", "0")
  mapping.set("S", "0")
  mapping.set("T", "0")
  // Group 1: (C) - cysteine
  mapping.set("C", "1")
  // Group 2: (D,E,N,Q) - charged/polar
  mapping.set("D", "2")
  mapping.set("E", "2")
  mapping.set("N", "2")
  mapping.set("Q", "2")
  // Group 3: (I,L,M,V) - hydrophobic/aliphatic
  mapping.set("I", "3")
  mapping.set("L", "3")
  mapping.set("M", "3")
  mapping.set("V", "3")
  // Group 4: (F,W,Y) - aromatic
  mapping.set("F", "4")
  mapping.set("W", "4")
  mapping.set("Y", "4")
  // Group 5: (H,K,R) - positively charged
  mapping.set("H", "5")
  mapping.set("K", "5")
  mapping.set("R", "5")
  let groups = ["(A,G,P,S,T)", "(C)", "(D,E,N,Q)", "(I,L,M,V)", "(F,W,Y)", "(H,K,R)"]
  ReducedAlphabet::new("Dayhoff", mapping, groups)
}

///|
/// CHARM (CHARacter Alphabet for Redundancy Minimization) - 5 groups.
/// From Seibert et al. (2002) "Structure-based alignment of remote homologs"
/// Groups: (A,C,L,M,F,I,V,W,Y,H,R,K), (S,T,N,Q,G,P,D,E), (S,T,N,Q,G,P), (D,E), (K,R,H)
/// Note: This is a simplified version.
pub fn charm_alphabet() -> ReducedAlphabet {
  let mapping : Map[String, String] = Map([], capacity=20)
  // Group 0: hydrophobic + aromatic + basic
  mapping.set("A", "0")
  mapping.set("C", "0")
  mapping.set("I", "0")
  mapping.set("L", "0")
  mapping.set("M", "0")
  mapping.set("F", "0")
  mapping.set("V", "0")
  mapping.set("W", "0")
  mapping.set("Y", "0")
  // Group 1: polar + small
  mapping.set("S", "1")
  mapping.set("T", "1")
  mapping.set("N", "1")
  mapping.set("Q", "1")
  mapping.set("G", "1")
  mapping.set("P", "1")
  // Group 2: charged
  mapping.set("D", "2")
  mapping.set("E", "2")
  // Group 3: basic
  mapping.set("H", "3")
  mapping.set("K", "3")
  mapping.set("R", "3")
  // Note: This is simplified (CHARM has overlapping assignments)
  let groups = ["(A,C,I,L,M,F,V,W,Y)", "(S,T,N,Q,G,P)", "(D,E)", "(H,K,R)"]
  ReducedAlphabet::new("CHARM", mapping, groups)
}

///|
/// SDM12 - 12 groups from Stanfel et al. (1996).
/// Groups amino acids into 12 classes based on physico-chemical properties.
pub fn sdm12_alphabet() -> ReducedAlphabet {
  let mapping : Map[String, String] = Map([], capacity=20)
  // 12 groups
  mapping.set("A", "0") // Ala - hydrophobic
  mapping.set("R", "1") // Arg - positive
  mapping.set("N", "2") // Asn - polar
  mapping.set("D", "3") // Asp - negative
  mapping.set("C", "4") // Cys - special
  mapping.set("Q", "5") // Gln - polar
  mapping.set("E", "6") // Glu - negative
  mapping.set("G", "7") // Gly - special
  mapping.set("H", "8") // His - positive
  mapping.set("I", "9") // Ile - hydrophobic
  mapping.set("L", "10") // Leu - hydrophobic
  mapping.set("K", "11") // Lys - positive
  mapping.set("M", "0") // Met - hydrophobic (same as Ala)
  mapping.set("F", "10") // Phe - hydrophobic (same as Leu)
  mapping.set("P", "7") // Pro - special (same as Gly)
  mapping.set("S", "2") // Ser - polar (same as Asn)
  mapping.set("T", "2") // Thr - polar (same as Asn)
  mapping.set("W", "10") // Trp - hydrophobic (same as Leu)
  mapping.set("Y", "5") // Tyr - polar (same as Gln)
  mapping.set("V", "9") // Val - hydrophobic (same as Ile)
  let groups = [
    "(A,M)", "(R)", "(N,S,T)", "(D)", "(C)",
    "(Q,Y)", "(E)", "(G,P)", "(H)", "(I,V)", "(L,F,W)", "(K)",
  ]
  ReducedAlphabet::new("SDM12", mapping, groups)
}

///|
/// Generate all built-in reduced alphabets.
pub fn all_reduced_alphabets() -> Array[ReducedAlphabet] {
  [rad_alphabet(), dayhoff_alphabet(), charm_alphabet(), sdm12_alphabet()]
}

///|
/// Get a reduced alphabet by name.
pub fn get_reduced_alphabet(name : String) -> ReducedAlphabet? {
  let lower = name.to_lower()
  if lower == "rad" {
    Some(rad_alphabet())
  } else if lower == "dayhoff" {
    Some(dayhoff_alphabet())
  } else if lower == "charm" {
    Some(charm_alphabet())
  } else if lower == "sdm12" {
    Some(sdm12_alphabet())
  } else {
    None
  }
}

///|
/// Compare two sequences using a reduced alphabet.
/// Returns the number of matching groups divided by total length.
pub fn reduced_identity(
  seq1 : String,
  seq2 : String,
  alphabet : ReducedAlphabet,
) -> Double {
  let n = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }
  if n == 0 {
    return 0.0
  }
  let mut matches = 0
  let mut i = 0
  while i < n {
    let aa1 = seq1.substring(start=i, end=i + 1)
    let aa2 = seq2.substring(start=i, end=i + 1)
    let g1 = alphabet.reduce_aa(aa1)
    let g2 = alphabet.reduce_aa(aa2)
    if g1 == g2 {
      matches = matches + 1
    }
    i = i + 1
  }
  matches.to_double() / n.to_double()
}

///|
/// Generate sample sequence data for testing.
pub fn reduced_sample_sequence() -> String {
  "MVLSPADKTNVKAAWGCK"
}
