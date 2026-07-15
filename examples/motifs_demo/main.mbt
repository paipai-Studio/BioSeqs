///|
/// Motifs analysis demo.

fn main {
  println("=== Motifs Analysis Demo ===\n")
  
  let sites = ["ACGT", "ACAT", "ACCT", "AGGT"]
  println("1. Create motif from sites:")
  let motif = @src.create_motif("TFBS", sites)
  println("   Motif name: \{motif.name}")
  println("   Motif width: \{motif.pwm.length()}")
  
  println("\n2. Motif consensus:")
  let consensus = @src.motif_consensus(motif)
  println("   Consensus: \{consensus}")
  
  println("\n3. Motif information content:")
  let ic = @src.motif_ic(motif)
  println("   IC: \{ic}")
  
  println("\n4. Search motif in sequence:")
  let sequence = "GACGTACGTACGTACGT"
  let hits = @src.search_motif(motif, sequence, 0.0)
  println("   Sequence: \{sequence}")
  println("   Hits found: \{hits.length()}")
  
  println("\n5. MEME format export:")
  let meme = @src.motif_to_meme(motif)
  println("   MEME output (first 5 lines):")
  let lines = meme.split("\n").to_array()
  for i = 0; i < 5 && i < lines.length(); i = i + 1 {
    println("     \{lines[i]}")
  }
  
  println("\n=== Demo Complete ===")
}