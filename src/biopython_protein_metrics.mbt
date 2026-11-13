/// biopython_protein_metrics.mbt — Protein sequence analysis metrics
/// Inspired by Biopython Bio.SeqUtils.ProteinAnalysis

///|
/// Kyte-Doolittle hydropathy values for each amino acid.
fn biopython_kd_value(aa : String) -> Double {
  match aa {
    "I" => 4.5
    "V" => 4.2
    "L" => 3.8
    "F" => 2.8
    "C" => 2.5
    "M" => 1.9
    "A" => 1.8
    "G" => -0.4
    "T" => -0.7
    "S" => -0.8
    "W" => -0.9
    "Y" => -1.3
    "P" => -1.6
    "H" => -3.2
    "E" => -3.5
    "Q" => -3.5
    "D" => -3.5
    "N" => -3.5
    "K" => -3.9
    "R" => -4.5
    _ => 0.0
  }
}

///|
/// GRAVY (Grand Average of Hydropathy).
pub fn biopython_gravy(protein : String) -> Double {
  let n = protein.length()
  if n == 0 {
    return 0.0
  }
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    total = total + biopython_kd_value(protein[i:i + 1].to_owned())
  }
  total / n.to_double()
}

///|
/// Aromaticity: fraction of aromatic amino acids (F, W, Y).
pub fn biopython_aromaticity(protein : String) -> Double {
  let n = protein.length()
  if n == 0 {
    return 0.0
  }
  let mut count = 0
  for i = 0; i < n; i = i + 1 {
    let aa = protein[i:i + 1].to_owned()
    if aa == "F" || aa == "W" || aa == "Y" {
      count = count + 1
    }
  }
  count.to_double() / n.to_double()
}

///|
/// Instability index (Guruprasad et al.).
/// DIWV (dipeptide instability weight values) for 400 dipeptide combinations.
/// Simplified: use a subset of the weight table.
pub fn biopython_instability_index(protein : String) -> Double {
  let n = protein.length()
  if n < 2 {
    return 0.0
  }
  // Simplified instability weights for common dipeptides
  let weight = fn(dp : String) -> Double {
    match dp {
      "QP" | "PQ" => 0.0
      "RP" | "PR" => 0.0
      "HP" | "PH" => 0.0
      "WP" | "PW" => 0.0
      "YP" | "PY" => 0.0
      "SP" | "PS" => 0.0
      "TP" | "PT" => 0.0
      "GP" | "PG" => 0.0
      "AP" | "PA" => 0.0
      "MP" | "PM" => 0.0
      "CP" | "PC" => 0.0
      "LP" | "PL" => 0.0
      "IP" | "PI" => 0.0
      "VP" | "PV" => 0.0
      "FP" | "PF" => 0.0
      "KP" | "PK" => 0.0
      "EP" | "PE" => 0.0
      "DP" | "PD" => 0.0
      "NP" | "PN" => 0.0
      "QR" | "RQ" => 1.0
      "QH" | "HQ" => 1.0
      "QW" | "WQ" => 1.0
      "QY" | "YQ" => 1.0
      "QS" | "SQ" => 1.0
      "QT" | "TQ" => 1.0
      "QG" | "GQ" => 1.0
      "QA" | "AQ" => 1.0
      "QM" | "MQ" => 1.0
      "QC" | "CQ" => 1.0
      "QL" | "LQ" => 1.0
      "QI" | "IQ" => 1.0
      "QV" | "VQ" => 1.0
      "QF" | "FQ" => 1.0
      "QK" | "KQ" => 1.0
      "QE" | "EQ" => 1.0
      "QD" | "DQ" => 1.0
      "QN" | "NQ" => 1.0
      _ => 10.0 // Default weight for other dipeptides
    }
  }
  let mut total = 0.0
  for i = 0; i < n - 1; i = i + 1 {
    let dp = protein[i:i + 2].to_owned()
    total = total + weight(dp)
  }
  // Scale factor: 10 / N
  total * 10.0 / n.to_double()
}

///|
/// Amino acid composition: returns array of (aa, count) sorted by count.
pub fn biopython_amino_acid_composition(
  protein : String,
) -> Array[(String, Int32)] {
  let aas = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S",
    "T", "V", "W", "Y",
  ]
  let counts : Array[Int32] = Array::make(aas.length(), 0)
  for i = 0; i < protein.length(); i = i + 1 {
    let aa = protein[i:i + 1].to_owned()
    for j = 0; j < aas.length(); j = j + 1 {
      if aa == aas[j] {
        counts[j] = counts[j] + 1
        break
      }
    }
  }
  let result : Array[(String, Int32)] = Array(capacity=0)
  for j = 0; j < aas.length(); j = j + 1 {
    result.push((aas[j], counts[j]))
  }
  result.sort_by(fn(a, b) {
    if b.1 < a.1 {
      -1
    } else if b.1 > a.1 {
      1
    } else {
      0
    }
  })
  result
}

///|
/// Molecular weight of a protein (average weights).
/// Delegates to the unified calculator in data.mbt.
pub fn biopython_protein_molecular_weight(protein : String) -> Double {
  protein_molecular_weight(protein)
}

///|
/// Count disulfide bond candidates (Cys pairs).
pub fn biopython_count_cysteines(protein : String) -> Int32 {
  let mut count = 0
  for i = 0; i < protein.length(); i = i + 1 {
    if protein[i:i + 1].to_owned() == "C" {
      count = count + 1
    }
  }
  count
}

///|
/// Flexibility scale (Vihinen, 1987).
/// Returns an array of flexibility values for each position.
pub fn biopython_flexibility(protein : String) -> Array[Double] {
  let n = protein.length()
  if n == 0 {
    return []
  }
  // Flexibility parameters for each amino acid
  let flex = fn(aa : String) -> Double {
    match aa {
      "A" => 0.357
      "C" => 0.346
      "D" => 0.511
      "E" => 0.495
      "F" => 0.274
      "G" => 0.544
      "H" => 0.323
      "I" => 0.242
      "K" => 0.439
      "L" => 0.247
      "M" => 0.297
      "N" => 0.469
      "P" => 0.509
      "Q" => 0.493
      "R" => 0.529
      "S" => 0.462
      "T" => 0.450
      "V" => 0.243
      "W" => 0.252
      "Y" => 0.284
      _ => 0.35
    }
  }
  let result : Array[Double] = Array::make(n, 0.0)
  let window = 9
  for i = 0; i < n; i = i + 1 {
    let start = if i - window / 2 > 0 { i - window / 2 } else { 0 }
    let end = if i + window / 2 < n { i + window / 2 + 1 } else { n }
    let mut sum = 0.0
    let count = end - start
    for k = start; k < end; k = k + 1 {
      sum = sum + flex(protein[k:k + 1].to_owned())
    }
    result[i] = if count > 0 { sum / count.to_double() } else { 0.0 }
  }
  result
}
