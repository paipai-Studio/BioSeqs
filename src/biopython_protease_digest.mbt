/// biopython_protease_digest.mbt — Proteolytic digestion simulation
/// Inspired by Biopython Bio.Restriction and protease cleavage rules

///|
/// Protease cleavage rule: cleaves after specific residues
/// under certain conditions.
pub struct BiopythonProtease {
  name : String
  cleave_after : Array[String] // residues to cleave after
  exclude_next : Array[String] // don't cleave if next residue is in this list
}

///|
/// Trypsin: cleaves after K and R, but not before P.
pub fn biopython_trypsin() -> BiopythonProtease {
  { name: "trypsin", cleave_after: ["K", "R"], exclude_next: ["P"], }
}

///|
/// Chymotrypsin: cleaves after F, Y, W (aromatic residues).
pub fn biopython_chymotrypsin() -> BiopythonProtease {
  { name: "chymotrypsin", cleave_after: ["F", "Y", "W"], exclude_next: ["P"], }
}

///|
/// Pepsin (pH 2): cleaves after F, Y, W, L, E.
pub fn biopython_pepsin() -> BiopythonProtease {
  { name: "pepsin", cleave_after: ["F", "Y", "W", "L", "E"], exclude_next: [], }
}

///|
/// LysC: cleaves after K.
pub fn biopython_lys_c() -> BiopythonProtease {
  { name: "LysC", cleave_after: ["K"], exclude_next: ["P"], }
}

///|
/// AspN: cleaves before D (special case).
pub fn biopython_asp_n() -> BiopythonProtease {
  {
    name: "AspN",
    cleave_after: ["D"], // Note: actually cleaves BEFORE D
    exclude_next: [],
  }
}

///|
/// Thermolysin: cleaves before A, F, I, L, M, V, Y.
pub fn biopython_thermolysin() -> BiopythonProtease {
  {
    name: "thermolysin",
    cleave_after: ["A", "F", "I", "L", "M", "V", "Y"],
    exclude_next: [],
  }
}

///|
/// Digest a protein sequence with a protease.
/// Returns array of peptide fragments.
pub fn biopython_digest(
  sequence : String,
  protease : BiopythonProtease,
  missed_cleavages~ : Int32 = 0,
  min_length~ : Int32 = 1,
) -> Array[String] {
  let n = sequence.length()
  if n == 0 {
    return []
  }
  // Find cleavage sites
  let sites : Array[Int32] = Array(capacity=0)
  for i = 0; i < n - 1; i = i + 1 {
    let aa = sequence[i:i + 1].to_owned()
    let next_aa = sequence[i + 1:i + 2].to_owned()
    let mut should_cleave = false
    for j = 0; j < protease.cleave_after.length(); j = j + 1 {
      if aa == protease.cleave_after[j] {
        should_cleave = true
        break
      }
    }
    if should_cleave {
      // Check exclusion
      let mut excluded = false
      for j = 0; j < protease.exclude_next.length(); j = j + 1 {
        if next_aa == protease.exclude_next[j] {
          excluded = true
          break
        }
      }
      if !excluded {
        sites.push(i + 1) // cleavage site is after position i
      }
    }
  }
  // Generate peptides
  let peptides : Array[String] = Array(capacity=0)
  // Add start and end as boundaries
  let bounds : Array[Int32] = [0]
  for i = 0; i < sites.length(); i = i + 1 {
    bounds.push(sites[i])
  }
  bounds.push(n)
  // Generate peptides with missed cleavages
  for mc = 0; mc <= missed_cleavages; mc = mc + 1 {
    for i = 0; i + 1 + mc < bounds.length(); i = i + 1 {
      let start = bounds[i]
      let end = bounds[i + 1 + mc]
      let peptide = sequence[start:end].to_owned()
      if peptide.length() >= min_length {
        // Check for duplicates
        let mut found = false
        for j = 0; j < peptides.length(); j = j + 1 {
          if peptides[j] == peptide {
            found = true
            break
          }
        }
        if !found {
          peptides.push(peptide)
        }
      }
    }
  }
  peptides
}

///|
/// Calculate peptide mass (monoisotopic).
pub fn biopython_peptide_mass(peptide : String) -> Double {
  // Monoisotopic masses of amino acid residues
  let mut mass = 18.0105646863 // H2O
  for i = 0; i < peptide.length(); i = i + 1 {
    let aa = peptide[i:i + 1].to_owned()
    mass = mass +
      (match aa {
        "G" => 57.021463735
        "A" => 71.037113805
        "S" => 87.032028435
        "P" => 97.052763875
        "V" => 99.068413945
        "T" => 101.047678505
        "C" => 103.009184775
        "I" | "L" => 113.084064015
        "N" => 114.042927470
        "D" => 115.026943065
        "Q" => 128.058577540
        "K" => 128.094963050
        "E" => 129.042592610
        "M" => 131.040484605
        "H" => 137.058911880
        "F" => 147.068413945
        "R" => 156.101111050
        "Y" => 163.053336640
        "W" => 186.079312950
        _ => 0.0
      })
  }
  mass
}

///|
/// Get peptide properties (mass, length, pI estimate).
pub fn biopython_peptide_properties(
  peptide : String,
) -> (Double, Int32, Double) {
  let mass = biopython_peptide_mass(peptide)
  let length = peptide.length()
  // Simple pI estimate based on charged residues
  let mut n_pos = 0.0 // K, R, H
  let mut n_neg = 0.0 // D, E
  for i = 0; i < length; i = i + 1 {
    let aa = peptide[i:i + 1].to_owned()
    match aa {
      "K" | "R" => n_pos = n_pos + 1.0
      "H" => n_pos = n_pos + 0.5
      "D" | "E" => n_neg = n_neg + 1.0
      _ => ()
    }
  }
  let pi = if n_pos == n_neg {
    7.0
  } else if n_pos > n_neg {
    7.0 + (n_pos - n_neg) * 1.5
  } else {
    7.0 - (n_neg - n_pos) * 1.5
  }
  (mass, length, pi)
}
