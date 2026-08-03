///|
/// meme module demo - MEME motif format parsing.
///
/// Demonstrates parsing of MEME Minimal Motif Format files containing
/// position-specific probability matrices (PSPMs) for discovered motifs,
/// including consensus generation and information content calculation.
///
/// Reference: Biopython Bio.motifs.meme
/// Run: moon run examples/meme_demo
fn main {
  println("=== meme Module Demo ===")
  println("MEME motif format parsing")
  println("Reference: Biopython Bio.motifs.meme")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse sample MEME text
  // ------------------------------------------------------------------
  println("=== 1. Parse Sample MEME Text ===")
  let text = @bio.meme_sample_text()
  let record = @bio.meme_parse(text)
  println("MEME version: " + record.version)
  println("Alphabet: " + record.alphabet)
  println("Strands: " + record.strands)
  println("Number of motifs: " + record.num_motifs().to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. Background frequencies
  // ------------------------------------------------------------------
  println("=== 2. Background Letter Frequencies ===")
  let mut bi = 0
  while bi < record.background_frequencies.length() {
    let pair = record.background_frequencies[bi]
    println("  " + pair.0 + ": " + pair.1.to_string())
    bi = bi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Motif details
  // ------------------------------------------------------------------
  println("=== 3. Motif Details ===")
  let mut i = 0
  while i < record.motifs.length() {
    let m = record.motifs[i]
    println("Motif " + (i + 1).to_string() + ":")
    println("  Name: " + m.name)
    println("  Alt name: " + m.alt_name)
    println("  Width: " + m.width.to_string())
    println("  Nsites: " + m.nsites.to_string())
    println("  E-value: " + m.evalue.to_string())
    println("  Strands: " + m.strands)
    if m.url.length() > 0 {
      println("  URL: " + m.url)
    }
    println("  Consensus: " + m.consensus())
    println("  Information content: " + m.information_content().to_string() + " bits")
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Position-Specific Probability Matrix
  // ------------------------------------------------------------------
  println("=== 4. PSPM Matrix (first motif) ===")
  if record.motifs.length() > 0 {
    let m = record.motifs[0]
    let alphabet = m.alphabet_to_array()
    // Header
    let header = StringBuilder::new()
    header.write_string("  pos  ")
    for a in alphabet {
      header.write_string("   " + a + "   ")
    }
    println(header.to_string())
    // Rows
    let mut pos = 0
    while pos < m.pspm.length() {
      let row = m.pspm[pos]
      let line = StringBuilder::new()
      line.write_string("  " + pos.to_string() + "   ")
      let mut j = 0
      while j < row.length() {
        let val = row[j]
        let raw = val.to_string()
        let s = if raw.length() > 4 { raw[0:4].to_string() } else { raw }
        // Pad to 4 chars
        let padded = if s.length() >= 4 { s } else { s + "    ".substring(start=0, end=4 - s.length()) }
        line.write_string(" " + padded + " ")
        j = j + 1
      }
      println(line.to_string())
      pos = pos + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Probability lookup
  // ------------------------------------------------------------------
  println("=== 5. Probability Lookup ===")
  if record.motifs.length() > 0 {
    let m = record.motifs[0]
    println("P(A at pos 0) = " + m.probability("A", 0).to_string())
    println("P(C at pos 0) = " + m.probability("C", 0).to_string())
    println("P(G at pos 0) = " + m.probability("G", 0).to_string())
    println("P(T at pos 0) = " + m.probability("T", 0).to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Find motif by name
  // ------------------------------------------------------------------
  println("=== 6. Find Motif by Name ===")
  if record.motifs.length() > 0 {
    let name = record.motifs[0].name
    match record.find(name) {
      Some(m) => println("Found motif '" + name + "': width=" + m.width.to_string())
      None => println("Motif '" + name + "' not found")
    }
    match record.find("nonexistent") {
      Some(_) => println("Unexpected: found nonexistent motif")
      None => println("Correctly did not find 'nonexistent'")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Sample motif and information content
  // ------------------------------------------------------------------
  println("=== 7. Sample Motif ===")
  let sample = @bio.meme_sample_motif()
  println("Created motif: " + sample.name)
  println("  Width: " + sample.width.to_string())
  println("  Consensus: " + sample.consensus())
  println("  Information content: " + sample.information_content().to_string() + " bits")
  println("")

  // ------------------------------------------------------------------
  // 8. Serialize record to string
  // ------------------------------------------------------------------
  println("=== 8. Serialize Record ===")
  let serialized = @bio.meme_record_to_string(record)
  let lines = serialized.split("\n").to_array()
  println("Serialized record (" + lines.length().to_string() + " lines):")
  let mut l = 0
  while l < 5 && l < lines.length() {
    println("  " + lines[l].to_string())
    l = l + 1
  }
  println("  ...")
  println("")

  println("=== Demo complete ===")
}
