///|
/// mast_motif module demo - MAST motif search output parsing.
///
/// Demonstrates parsing of MAST (Motif Alignment and Search Tool) text
/// output files containing motif search results with p-values, E-values,
/// and positional information for matches in sequences.
///
/// Reference: Biopython Bio.motifs.mast
/// Run: moon run examples/mast_motif_demo
fn main {
  println("=== mast_motif Module Demo ===")
  println("MAST motif search output parsing")
  println("Reference: Biopython Bio.motifs.mast")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse sample MAST output
  // ------------------------------------------------------------------
  println("=== 1. Parse Sample MAST Output ===")
  let text = @bio.mast_motif_sample_text()
  let record = @bio.mast_motif_parse(text)
  println("MAST version: " + record.version)
  println("Database: " + record.database)
  println("Motif count: " + record.motif_count.to_string())
  println("Sequence hit count: " + record.hit_count.to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. Motif info
  // ------------------------------------------------------------------
  println("=== 2. Motifs Used in Search ===")
  let mut i = 0
  while i < record.motifs.length() {
    let m = record.motifs[i]
    println("  Motif " + (i + 1).to_string() + ": id=" + m.id + ", width=" + m.width.to_string())
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Sequence matches
  // ------------------------------------------------------------------
  println("=== 3. Sequence Matches ===")
  let mut s = 0
  while s < record.sequences.length() {
    let seq = record.sequences[s]
    println("  Sequence: " + seq.name)
    if seq.description.length() > 0 {
      println("    Description: " + seq.description)
    }
    println("    Length: " + seq.length.to_string())
    println("    Combined p-value: " + seq.combined_pvalue.to_string())
    println("    Combined E-value: " + seq.combined_evalue.to_string())
    println("    Number of matches: " + seq.num_matches().to_string())
    s = s + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Best sequence by E-value
  // ------------------------------------------------------------------
  println("=== 4. Best Sequence by E-value ===")
  match record.best_sequence() {
    Some(seq) => {
      println("  Best sequence: " + seq.name)
      println("  E-value: " + seq.combined_evalue.to_string())
      println("  p-value: " + seq.combined_pvalue.to_string())
    }
    None => println("  No sequences found")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Best match per sequence
  // ------------------------------------------------------------------
  println("=== 5. Best Match Per Sequence ===")
  let mut si = 0
  while si < record.sequences.length() {
    let seq = record.sequences[si]
    match seq.best_match() {
      Some(m) => println("  " + seq.name + ": motif=" + m.motif_id + ", pos=" + m.position.to_string() + ", p=" + m.pvalue.to_string())
      None => println("  " + seq.name + ": no matches")
    }
    si = si + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Find sequence by name
  // ------------------------------------------------------------------
  println("=== 6. Find Sequence by Name ===")
  if record.sequences.length() > 0 {
    let name = record.sequences[0].name
    match record.find_sequence(name) {
      Some(seq) => println("  Found '" + name + "': E=" + seq.combined_evalue.to_string())
      None => println("  Not found: " + name)
    }
    match record.find_sequence("nonexistent") {
      Some(_) => println("  Unexpected: found nonexistent")
      None => println("  Correctly did not find 'nonexistent'")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Motif matches across sequences
  // ------------------------------------------------------------------
  println("=== 7. Motif Matches Across Sequences ===")
  if record.motifs.length() > 0 {
    let motif_id = record.motifs[0].id
    let matches = record.motif_matches(motif_id)
    println("  Motif '" + motif_id + "' matches: " + matches.length().to_string())
    let mut m = 0
    while m < matches.length() {
      let mt = matches[m]
      println("    pos=" + mt.position.to_string() + ", p=" + mt.pvalue.to_string() + ", strand=" + mt.strand)
      m = m + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Number of hits and statistics
  // ------------------------------------------------------------------
  println("=== 8. Statistics ===")
  println("  Total sequences with hits: " + record.num_hits().to_string())
  println("  Total motifs in search: " + record.motifs.length().to_string())
  println("  Total sequences: " + record.sequences.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 9. Serialize record
  // ------------------------------------------------------------------
  println("=== 9. Serialize Record ===")
  let serialized = @bio.mast_motif_record_to_string(record)
  let lines = serialized.split("\n").to_array()
  println("  Serialized (" + lines.length().to_string() + " lines), first 5:")
  let mut l = 0
  while l < 5 && l < lines.length() {
    println("    " + lines[l].to_string())
    l = l + 1
  }
  println("")

  println("=== Demo complete ===")
}
