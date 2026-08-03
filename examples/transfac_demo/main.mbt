///|
/// transfac module demo - TRANSFAC motif format parsing.
///
/// Demonstrates parsing of TRANSFAC database files containing transcription
/// factor binding site matrices with two-letter field codes, position
/// frequency matrices (PFM), and rich metadata including references.
///
/// Reference: Biopython Bio.motifs.transfac
/// Run: moon run examples/transfac_demo
fn main {
  println("=== transfac Module Demo ===")
  println("TRANSFAC transcription factor binding profile parsing")
  println("Reference: Biopython Bio.motifs.transfac")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse sample TRANSFAC text
  // ------------------------------------------------------------------
  println("=== 1. Parse Sample TRANSFAC Text ===")
  let text = @bio.transfac_sample_text_full()
  let record = @bio.transfac_parse_full(text)
  println("Version: " + record.version)
  println("Number of motifs: " + record.num_motifs().to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. Motif metadata
  // ------------------------------------------------------------------
  println("=== 2. Motif Metadata ===")
  let mut i = 0
  while i < record.motifs.length() {
    let m = record.motifs[i]
    println("Motif " + (i + 1).to_string() + ":")
    println("  Accession (AC): " + m.accession)
    println("  Identifier (ID): " + m.identifier)
    println("  Name (NA): " + m.name)
    println("  Description (DE): " + m.description)
    println("  Species (OS): " + m.species)
    if m.taxonomic_class.length() > 0 {
      println("  Taxonomic class (OC): " + m.taxonomic_class)
    }
    println("  Binding factors (BF): " + m.binding_factors.length().to_string())
    let mut b = 0
    while b < m.binding_factors.length() {
      println("    - " + m.binding_factors[b])
      b = b + 1
    }
    if m.comments.length() > 0 {
      println("  Comments (CC): " + m.comments)
    }
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Position Frequency Matrix (PFM)
  // ------------------------------------------------------------------
  println("=== 3. Position Frequency Matrix (PFM) ===")
  let mut mi = 0
  while mi < record.motifs.length() {
    let m = record.motifs[mi]
    println("Motif " + m.accession + " (" + m.identifier + "):")
    println("  Width: " + m.width.to_string())
    println("  Alphabet: " + alphabet_to_string_impl(m.alphabet))
    // Print PFM header
    let header = StringBuilder::new()
    header.write_string("  pos  ")
    for a in m.alphabet {
      header.write_string(" " + a + " ")
    }
    println(header.to_string())
    // Print PFM rows (transposed: positions are columns)
    let mut col = 0
    while col < m.width {
      let row = StringBuilder::new()
      let row_num = if col + 1 < 10 { "0" + (col + 1).to_string() } else { (col + 1).to_string() }
      row.write_string("  " + row_num + "  ")
      let mut r = 0
      while r < m.pfm.length() {
        if col < m.pfm[r].length() {
          let v = m.pfm[r][col]
          let s = if v < 10 { " " + v.to_string() } else { v.to_string() }
          row.write_string(" " + s + " ")
        }
        r = r + 1
      }
      println(row.to_string())
      col = col + 1
    }
    mi = mi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Consensus sequence
  // ------------------------------------------------------------------
  println("=== 4. Consensus Sequence ===")
  for m in record.motifs {
    println("  " + m.accession + " (" + m.identifier + "): " + m.consensus())
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Total counts and frequencies
  // ------------------------------------------------------------------
  println("=== 5. Total Counts and Frequencies ===")
  if record.motifs.length() > 0 {
    let m = record.motifs[0]
    let mut pos = 0
    while pos < m.width {
      let total = m.total_at(pos)
      let cons_base = m.consensus()[pos:pos + 1].to_owned()
      let cons_freq = m.frequency(cons_base, pos)
      println("  pos " + pos.to_string() + ": total=" + total.to_string() + ", consensus=" + cons_base + " (" + cons_freq.to_string() + ")")
      pos = pos + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Find motif by accession
  // ------------------------------------------------------------------
  println("=== 6. Find Motif by Accession ===")
  match record.find("T00001") {
    Some(m) => println("  Found T00001: " + m.identifier + " (" + m.species + ")")
    None => println("  T00001 not found")
  }
  match record.find("nonexistent") {
    Some(_) => println("  Unexpected: found nonexistent")
    None => println("  Correctly did not find 'nonexistent'")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Alphabet index lookup
  // ------------------------------------------------------------------
  println("=== 7. Alphabet Index Lookup ===")
  if record.motifs.length() > 0 {
    let m = record.motifs[0]
    println("  A index: " + m.alphabet_index("A").to_string())
    println("  C index: " + m.alphabet_index("C").to_string())
    println("  G index: " + m.alphabet_index("G").to_string())
    println("  T index: " + m.alphabet_index("T").to_string())
    println("  N index: " + m.alphabet_index("N").to_string() + " (not in alphabet)")
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Serialize motif back to TRANSFAC format
  // ------------------------------------------------------------------
  println("=== 8. Serialize Motif to TRANSFAC Format ===")
  if record.motifs.length() > 0 {
    let m = record.motifs[0]
    let s = @bio.transfac_motif_full_to_string(m)
    let lines = s.split("\n").to_array()
    println("  Serialized motif " + m.accession + " (" + lines.length().to_string() + " lines):")
    let mut l = 0
    while l < lines.length() && l < 10 {
      println("    " + lines[l].to_string())
      l = l + 1
    }
    if lines.length() > 10 {
      println("    ...")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 9. References
  // ------------------------------------------------------------------
  println("=== 9. References ===")
  let mut ri = 0
  while ri < record.motifs.length() {
    let m = record.motifs[ri]
    println("  Motif " + m.accession + " references: " + m.references.length().to_string())
    let mut r = 0
    while r < m.references.length() {
      let ref_ = m.references[r]
      println("    RN: " + ref_.number)
      if ref_.authors.length() > 0 {
        println("    RA: " + ref_.authors)
      }
      if ref_.title.length() > 0 {
        println("    RT: " + ref_.title)
      }
      if ref_.journal.length() > 0 {
        println("    RL: " + ref_.journal)
      }
      r = r + 1
    }
    ri = ri + 1
  }
  println("")

  println("=== Demo complete ===")
}

///|
/// Helper to convert alphabet array to comma-separated string.
fn alphabet_to_string_impl(alphabet : Array[String]) -> String {
  let sb = StringBuilder::new()
  let mut first = true
  for a in alphabet {
    if !first {
      sb.write_string(", ")
    }
    sb.write_string(a)
    first = false
  }
  sb.to_string()
}
