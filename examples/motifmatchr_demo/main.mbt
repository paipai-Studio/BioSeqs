///|
fn main {
  println("=== motifmatchr Demo ===")

  let motif = @src.motifmatchr_sample_motif()
  let seqs = @src.motifmatchr_sample_sequences()
  println("\n1. Motif: " + motif.name + " (width=" + motif.len().to_string() + ")")
  println("   Sequences: " + seqs.length().to_string())

  println("\n2. Matching motifs (cutoff = -5.0)")
  let results = @src.match_motifs([motif], seqs, cutoff=(-5.0))
  let r = results[0]
  for i in 0..<r.seq_names.length() {
    println(
      "   " + r.seq_names[i] + ": hits=" + r.counts[i].to_string() +
      ", max_score=" + r.max_scores[i].to_string(),
    )
  }

  println("\n3. Hit positions")
  let positions = r.hit_positions()
  for p in positions {
    println("   seq " + p.0.to_string() + " at position " + p.1.to_string())
  }

  println("\n4. Summary")
  println("   Total hits: " + r.total_hits().to_string())
  println("   Best score: " + r.best_score().to_string())

  println("\n5. Reverse complement")
  println("   revcomp(ACTG) = " + @src.motifmatchr_revcomp("ACTG"))
  println("   revcomp(ACGT) = " + @src.motifmatchr_revcomp("ACGT"))
}