///|
/// Protein sequence analysis demo.
///
/// Demonstrates the protein analysis module for:
///   - Hydropathy plotting (Kyte-Doolittle scale)
///   - Secondary structure prediction (GOR method)
///   - Antigenicity prediction (Hopp-Woods scale)
///   - Transmembrane segment prediction
///   - Amino acid composition analysis
///   - Sequence conservation scoring
///
/// Inspired by Biopython's Bio.SeqUtils.ProtParam and related modules.
/// Run: moon run examples/protein_analysis_demo
fn main {
  println("=== Protein Sequence Analysis Demo ===")
  println("Advanced protein sequence analysis tools")
  println("Reference: Biopython Bio.SeqUtils.ProtParam")
  println("")

  let globin = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLR"

  // ------------------------------------------------------------------
  // 1. Hydropathy analysis
  // ------------------------------------------------------------------
  println("=== 1. Hydropathy Analysis (Kyte-Doolittle) ===")
  println("Sequence: Hemoglobin beta subunit (first 50 aa)")
  let hydro = @bio.protein_hydropathy(globin, window=7)
  println("Hydropathy values (window=7):")
  let mut h = 0
  while h < hydro.length() && h < 20 {
    println("  Position " + (h + 1).to_string() + ": " + hydro[h].to_string())
    h = h + 1
  }
  println("  ... (showing first 20 of " + hydro.length().to_string() + " values)")
  println("")

  // Find hydrophobic regions (potential transmembrane domains)
  println("Hydrophobic regions (score > 1.6):")
  let mut ht = 0
  while ht < hydro.length() {
    if hydro[ht] > 1.6 {
      println("  Position " + (ht + 1).to_string() + "-" + (ht + 7).to_string() + ": " + hydro[ht].to_string())
    }
    ht = ht + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Secondary structure prediction
  // ------------------------------------------------------------------
  println("=== 2. Secondary Structure Prediction (GOR Method) ===")
  let ss = @bio.protein_secondary_structure(globin)
  println("Secondary structure assignments:")
  let mut s = 0
  while s < ss.length() && s < 50 {
    let residue = globin[s:s + 1].to_owned()
    println("  Position " + (s + 1).to_string() + ": " + residue + " -> " + ss[s])
    s = s + 1
  }
  println("  ... (showing first 50 of " + ss.length().to_string() + " residues)")

  // Count structure types
  let mut helix = 0
  let mut sheet = 0
  let mut turn = 0
  let mut coil = 0
  let mut c = 0
  while c < ss.length() {
    if ss[c] == "H" { helix = helix + 1 }
    else if ss[c] == "E" { sheet = sheet + 1 }
    else if ss[c] == "T" { turn = turn + 1 }
    else { coil = coil + 1 }
    c = c + 1
  }
  println("")
  println("Structure summary:")
  println("  Helix (H): " + helix.to_string() + " residues")
  println("  Sheet (E): " + sheet.to_string() + " residues")
  println("  Turn (T): " + turn.to_string() + " residues")
  println("  Coil (C): " + coil.to_string() + " residues")
  println("")

  // ------------------------------------------------------------------
  // 3. Antigenicity prediction
  // ------------------------------------------------------------------
  println("=== 3. Antigenicity Prediction (Hopp-Woods) ===")
  let antigen = @bio.protein_antigenicity(globin)
  println("Antigenicity profile:")
  let mut a = 0
  let mut positive_count = 0
  while a < antigen.length() {
    if antigen[a] > 0.0 { positive_count = positive_count + 1 }
    a = a + 1
  }
  println("  Total residues: " + antigen.length().to_string())
  println("  Positive antigenic sites: " + positive_count.to_string())
  println("  Percentage positive: " + ((positive_count.to_double() / antigen.length().to_double() * 100.0).to_string()) + "%")
  println("")

  // Show top antigenic regions
  println("Top antigenic regions (score > 0.2):")
  let mut ar = 0
  while ar < antigen.length() {
    if antigen[ar] > 0.2 {
      let start = ar
      while ar < antigen.length() && antigen[ar] > 0.2 {
        ar = ar + 1
      }
      let region = globin[start:ar].to_owned()
      println("  Positions " + (start + 1).to_string() + "-" + ar.to_string() + ": " + region)
    }
    ar = ar + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Transmembrane prediction
  // ------------------------------------------------------------------
  println("=== 4. Transmembrane Segment Prediction ===")
  let tm_seq = "MFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLR"
  let tm = @bio.protein_transmembrane(tm_seq)
  println("Transmembrane propensity profile:")
  let mut tm_print = 0
  while tm_print < tm.length() && tm_print < 20 {
    println("  Position " + (tm_print + 1).to_string() + ": " + tm[tm_print].to_string())
    tm_print = tm_print + 1
  }
  println("  ... (showing first 20 of " + tm.length().to_string() + " values)")
  println("")

  // ------------------------------------------------------------------
  // 5. Amino acid composition
  // ------------------------------------------------------------------
  println("=== 5. Amino Acid Composition ===")
  let comp = @bio.protein_aa_composition(globin)
  println("Amino acid frequencies:")
  let aa_list = comp.keys().collect()
  // Sort for consistent output
  aa_list.sort_by(fn(a : Char, b : Char) -> Int { a.to_int() - b.to_int() })
  let mut k = 0
  while k < aa_list.length() {
    let aa = aa_list[k]
    let freq = comp[aa]
    println("  " + aa.to_string() + ": " + (freq * 100.0).to_string() + "%")
    k = k + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Sequence conservation
  // ------------------------------------------------------------------
  println("=== 6. Sequence Conservation ===")
  println("Testing conservation across an alignment:")
  let alignment = ["MVLSPADKTNVKAAWGKVGAHAGEYGAEALER", "MVLSAADKTNVKAAWSKVGGHAGEYGAEALER", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALER"]
  let cons = @bio.protein_conservation(alignment)
  println("Conservation scores (0-1, higher = more conserved):")
  let mut c2 = 0
  while c2 < cons.length() {
    let bar = if cons[c2] > 0.8 { "*****" }
      else if cons[c2] > 0.6 { "****" }
      else if cons[c2] > 0.4 { "***" }
      else if cons[c2] > 0.2 { "**" }
      else { "*" }
    println("  Position " + (c2 + 1).to_string() + ": " + cons[c2].to_string() + " " + bar)
    c2 = c2 + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Summary report
  // ------------------------------------------------------------------
  println("=== 7. Summary Report ===")
  let summary = @bio.protein_summary(globin)
  println(summary)

  println("=== Protein Analysis Demo Complete ===")
}
