///|
/// Demo for xms module: XMS weighted site motif + microRNA motif formats.

fn main {
  println("=== xms Demo: XMS & microRNA Motif Formats ===")
  println("")

  // Example 1: Parse XMS format
  println("--- Step 1: Parse XMS Format ---")
  let xms_text = ">TF1\nACGTACGT\nTGCAACGT\nACGTACGT 2.0\n>TF2\nAAAAACCC\nGGGGGTTT\n"
  let motifs = @src.xms_parse(xms_text)
  println("Parsed motifs: " + motifs.length().to_string())
  for m in motifs {
    println("  " + m.name + ": " + m.sites.length().to_string() + " sites")
  }
  println("")

  // Example 2: Count matrix
  println("--- Step 2: Count Matrix ---")
  let counts = @src.xms_counts(motifs[0])
  println("Count matrix: " + counts.length().to_string() + " x " + counts[0].length().to_string())
  println("  A: " + counts[0].to_string())
  println("  C: " + counts[1].to_string())
  println("")

  // Example 3: PFM (frequency matrix)
  println("--- Step 3: PFM (Frequency Matrix) ---")
  let pfm = @src.xms_pfm(counts)
  println("PFM: " + pfm.length().to_string() + " positions")
  for i in 0..<pfm.length() {
    println("  pos " + (i+1).to_string() + ": A=" + pfm[i][0].to_string() + " C=" + pfm[i][1].to_string() + " G=" + pfm[i][2].to_string() + " T=" + pfm[i][3].to_string())
  }
  println("")

  // Example 4: PWM (log-odds)
  println("--- Step 4: PWM (Log-Odds) ---")
  let pwm = @src.xms_pwm(pfm)
  println("PWM: " + pwm.length().to_string() + " positions")
  for i in 0..<pwm.length() {
    println("  pos " + (i+1).to_string() + ": A=" + pwm[i][0].to_string() + " C=" + pwm[i][1].to_string() + " G=" + pwm[i][2].to_string() + " T=" + pwm[i][3].to_string())
  }
  println("")

  // Example 5: Consensus & information content
  println("--- Step 5: Consensus & Information Content ---")
  let cons = @src.xms_consensus(pfm)
  let (per_pos, total) = @src.xms_info_content(pfm)
  println("Consensus: " + cons)
  println("Total information content: " + total.to_string() + " bits")
  println("Per-position IC: " + per_pos.to_string())
  println("")

  // Example 6: Round-trip write/parse
  println("--- Step 6: Round-Trip Write/Parse ---")
  let written = @src.xms_write(motifs[0])
  let reparsed = @src.xms_parse(written)
  println("Round-trip OK: " + reparsed.length().to_string() + " motif(s)")
  println("")

  // Example 7: microRNA seed classification
  println("--- Step 7: microRNA Seed Classification ---")
  let mature = "UGAGGUAGUAGGUUGUAUAGUU"
  let seeds = @src.mirs_classify_seeds(mature)
  println("Mature miRNA: " + mature)
  println("Seed types: " + seeds.length().to_string())
  for entry in seeds {
    let (label, seq) = entry
    println("  " + label + ": " + seq)
  }
  println("")

  // Example 8: Reverse complement
  println("--- Step 8: Reverse Complement ---")
  let rc = @src.mirs_reverse_complement("ACGTACGT")
  println("ACGTACGT -> " + rc)
  println("")

  // Example 9: microRNA PFM from sequences
  println("--- Step 9: microRNA PFM ---")
  let mir_seqs = ["UGA", "UGA", "UGC"]
  let mir_pfm = @src.mirs_pfm(mir_seqs)
  println("microRNA PFM: " + mir_pfm.length().to_string() + " positions")
  for i in 0..<mir_pfm.length() {
    println("  pos " + (i+1).to_string() + ": A=" + mir_pfm[i][0].to_string() + " C=" + mir_pfm[i][1].to_string() + " G=" + mir_pfm[i][2].to_string() + " T=" + mir_pfm[i][3].to_string())
  }
  println("")

  println("=== xms Demo Complete ===")
}
