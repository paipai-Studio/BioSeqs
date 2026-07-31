///|
/// Markov module example - Higher-order Markov chain modeling of biological sequences.
///
/// Run: moon run examples/markov_demo
fn main {
  println("=== Bio.Markov Module Demo ===")
  println("Higher-order Markov Chain Modeling for Biological Sequences")
  println("")

  // ---- Section 1: Build first-order DNA Markov model ----
  println("=== 1. Build 1st-order Markov Model (DNA) ===")
  let dna_seqs = [
    "ACGTACGTACGTACGT",
    "CGTACGTACGTACGTA",
    "GTACGTACGTACGTAC",
    "TACGTACGTACGTACG",
    "AAACCCAAACCCAAAC",
    "CCCGGGCCCGGGCCCG",
  ]
  let model1 = @bio.markov_build_model(dna_seqs, order=1)
  println("Training sequences: " + dna_seqs.length().to_string())
  println("States (alphabet): " + model1.states.length().to_string())
  println("Transition contexts: " + model1.transition_probs.length().to_string())
  println("Pseudo-count smoothing: " + model1.pseudo_count.to_string())
  println("")

  // ---- Section 2: Build 2nd order model ----
  println("=== 2. Build 2nd-order Markov Model ===")
  let model2 = @bio.markov_build_model(dna_seqs, order=2)
  println("Transition contexts (2-mers): " + model2.transition_probs.length().to_string())
  let test_dna = "ACGTACGT"
  let score1 = @bio.markov_score_sequence(model1, test_dna)
  let score2 = @bio.markov_score_sequence(model2, test_dna)
  println("Score(seq=\"" + test_dna + "\"):")
  println("  Order-1: " + score1.to_string())
  println("  Order-2: " + score2.to_string())
  println("")

  // ---- Section 3: Per-position scores ----
  println("=== 3. Per-Position Log Probabilities ===")
  let per_pos = @bio.markov_score_per_base(model1, test_dna)
  println("Sequence: " + test_dna)
  let mut ppi = 0
  while ppi < per_pos.length() {
    let position = model1.order + ppi
    let prev = test_dna[position - 1:position].to_string()
    let next = test_dna[position:position + 1].to_string()
    println("  pos " + position.to_string() + ": P(" + next + " | " + prev + ")_log = " + per_pos[ppi].to_string())
    ppi = ppi + 1
  }
  println("")

  // ---- Section 4: CpG island detection with log-odds ratio ----
  println("=== 4. CpG Island Detection (Log-Odds Ratio) ===")
  // Build two models: one from CpG-rich sequences, one from background DNA
  let cpg_seqs = [
    "CGCGCGCGCGCGCGCGC",
    "GCGCGCGCGCGCGCGCG",
    "ACGTACGTACGTACGTACGT",
    "CGATCGATCGATCGATCGAT",
    "CGCCGGCGCCGGCGCCGGCG",
  ]
  let bg_seqs = [
    "ATATATATATATATATA",
    "TATATATATATATATAT",
    "AAATTTAAATTTAAA",
    "GGCCTTAAGGCCGGCC",
    "TTCCCAAAACCCTTTT",
  ]
  let cpg_model = @bio.markov_build_model(cpg_seqs, order=1, pseudo_count=0.5)
  let bg_model  = @bio.markov_build_model(bg_seqs, order=1, pseudo_count=0.5)
  let test_cpg = "CGCGCGCGCGCGC"
  let test_at  = "ATATATATATATA"
  let test_mix = "ACGTACGTACGTACGT"
  println("Sequence: " + test_cpg + " (CpG-rich)")
  println("  Log-odds (CpG / BG) = " + @bio.markov_log_odds(cpg_model, bg_model, test_cpg).to_string() + " (expected > 0)")
  println("Sequence: " + test_at + " (AT-rich)")
  println("  Log-odds (CpG / BG) = " + @bio.markov_log_odds(cpg_model, bg_model, test_at).to_string() + " (expected < 0)")
  println("Sequence: " + test_mix + " (mixed)")
  println("  Log-odds (CpG / BG) = " + @bio.markov_log_odds(cpg_model, bg_model, test_mix).to_string())
  println("")

  // ---- Section 5: Stationary distribution ----
  println("=== 5. Stationary Distribution of the DNA Markov Chain ===")
  let stat = @bio.markov_stationary_distribution(model1)
  let keys = stat.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    let k = keys[ki]
    let p = stat.get(k).unwrap_or(0.0)
    println("  P(\u{03C0} = " + k + ") = " + p.to_string())
    ki = ki + 1
  }
  println("")

  // ---- Section 6: Generate synthetic sequences ----
  println("=== 6. Synthetic Sequence Generation ===")
  let dna_model = @bio.markov_build_model(dna_seqs, order=1, pseudo_count=1.0)
  let len30 = @bio.markov_generate_sequence(dna_model, 30, seed=42)
  let len50 = @bio.markov_generate_sequence(dna_model, 50, seed=99)
  println("Generated (len=30): " + len30)
  println("Generated (len=50): " + len50)
  println("")

  // ---- Section 7: Protein Markov model ----
  println("=== 7. Protein Markov Model (order=1) ===")
  let protein_seqs = [
    "MKVLFSAVLILFFVSTQAD",
    "MAGWQVLTALIGAMLCGGP",
    "MALMFAVLISVLLPQGEG",
    "MKILAVLISLLLSPWAGS",
  ]
  let prot_model = @bio.MarkovModel::default_protein()
    .set_pseudo_count(1.0)
  let trained_prot = @bio.markov_build_model(
    protein_seqs,
    order=1,
    states=prot_model.states,
    pseudo_count=1.0,
  )
  let unknown_prot = "XXXXXXXX"
  let test_prot = "MKVLFSAVLILF"
  println("States loaded: " + trained_prot.states.length().to_string())
  println("Score on unseen prot: " + test_prot + " -> " +
    @bio.markov_score_sequence(trained_prot, test_prot).to_string())
  println("Score on non-AA chars: " + unknown_prot + " -> " +
    @bio.markov_score_sequence(trained_prot, unknown_prot).to_string())
  println("")

  // ---- Section 8: Order comparison with ChainType helpers ----
  println("=== 8. Model Order Comparison ===")
  let m1 = @bio.markov_build_model(dna_seqs, order=1)
  let m2 = @bio.markov_build_model(dna_seqs, order=2)
  let m3 = @bio.markov_build_model(dna_seqs, order=3)
  let seq = dna_seqs[0]
  println("Reference sequence: " + seq)
  println("  Order-1 score: " + @bio.markov_score_sequence(m1, seq).to_string())
  println("  Order-2 score: " + @bio.markov_score_sequence(m2, seq).to_string())
  println("  Order-3 score: " + @bio.markov_score_sequence(m3, seq).to_string())
  println("")

  println("=== Demo complete ===")
}
