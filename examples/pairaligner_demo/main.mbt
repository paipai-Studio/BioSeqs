///|
/// PairwiseAligner module example - Unified pairwise sequence alignment.
///
/// Run: moon run examples/pairaligner_demo
fn main {
  println("=== Bio.Align.PairwiseAligner Module Demo ===")
  println("Unified Pairwise Sequence Alignment (Needleman-Wunsch / Smith-Waterman)")
  println("")

  // ---- Section 1: Global DNA alignment (Needleman-Wunsch) ----
  println("=== 1. Global DNA Alignment (Needleman-Wunsch) ===")
  let dna_cfg = @bio.PairwiseAlignerConfig::default_dna()
  let dna_target = "ACGTACGTACGT"
  let dna_query  = "ACGTAGGTACGT"
  let aln_global = @bio.pairaligner_align(dna_target, dna_query, config=dna_cfg)
  println("Target:  " + dna_target)
  println("Query:   " + dna_query)
  println("Score:   " + aln_global.score().to_string())
  println("Aligned:")
  println("  T:     " + aln_global.aligned1())
  println("  Match: " + aln_global.match_line)
  println("  Q:     " + aln_global.aligned2())
  println("Identities: " + aln_global.identities().to_string() + " / " +
    aln_global.alignment_length().to_string() + " (" + aln_global.identity_pct().to_string() + "%)")
  println("")

  // ---- Section 2: Local DNA alignment (Smith-Waterman) ----
  println("=== 2. Local DNA Alignment (Smith-Waterman) ===")
  let local_cfg = @bio.PairwiseAlignerConfig::default_dna()
    .set_mode(@bio.pairaligner_local())
    .set_match_score(3.0)
    .set_mismatch_score(-2.0)
  let long_target = "XXXXXXXXXXACGTACGTXXXXXXXXXX"
  let short_query = "ACGTACGT"
  let aln_local = @bio.pairaligner_align(long_target, short_query, config=local_cfg)
  println("Target context: ..." + long_target[10:18].to_string() + "...")
  println("Query:          " + short_query)
  println("Score: " + aln_local.score().to_string())
  println("Local alignment:")
  println("  T[" + aln_local.target_start.to_string() + ":" + aln_local.target_end.to_string() + "] " + aln_local.aligned1())
  println("  Match: " + aln_local.match_line)
  println("  Q[" + aln_local.query_start.to_string() + ":" + aln_local.query_end.to_string() + "] " + aln_local.aligned2())
  println("")

  // ---- Section 3: Protein alignment with BLOSUM62 ----
  println("=== 3. Protein Alignment (BLOSUM62 Substitution Matrix) ===")
  let (prot_target, prot_query) = @bio.pairaligner_sample_data()
  let prot_cfg = @bio.PairwiseAlignerConfig::default_protein()
  let aln_prot = @bio.pairaligner_align(prot_target, prot_query, config=prot_cfg)
  println("Protein 1: " + prot_target)
  println("Protein 2: " + prot_query)
  println("Score:     " + aln_prot.score().to_string())
  println("Aligned:")
  println("  " + aln_prot.aligned1())
  println("  " + aln_prot.match_line)
  println("  " + aln_prot.aligned2())
  println("Identity: " + aln_prot.identity_pct().to_string() + "%")
  println("Gaps count: " + aln_prot.gaps_count().to_string())
  println("")

  // ---- Section 4: Affine gaps vs linear gaps comparison ----
  println("=== 4. Affine vs Linear Gap Penalties ===")
  let t = "AAAAAAAAAA"
  let q = "AAAA"
  let affine = @bio.PairwiseAlignerConfig::default_dna()
    .set_mode(@bio.pairaligner_global())
    .set_gap_open(-5.0)
    .set_gap_extend(-1.0)
  let linear = @bio.PairwiseAlignerConfig::default_dna()
    .set_mode(@bio.pairaligner_global())
    .set_gap_open(-1.0)
    .set_gap_extend(-1.0)
  let a_aln = @bio.pairaligner_align(t, q, config=affine)
  let l_aln = @bio.pairaligner_align(t, q, config=linear)
  println("Target len=" + t.length().to_string() + ", Query len=" + q.length().to_string())
  println("Affine (open=-5, extend=-1): score=" + a_aln.score().to_string() +
    ", alignment_len=" + a_aln.alignment_length().to_string())
  println("  " + a_aln.aligned1())
  println("  " + a_aln.aligned2())
  println("Linear (open=-1, extend=-1): score=" + l_aln.score().to_string() +
    ", alignment_len=" + l_aln.alignment_length().to_string())
  println("  " + l_aln.aligned1())
  println("  " + l_aln.aligned2())
  println("")

  println("=== Demo complete ===")
}
