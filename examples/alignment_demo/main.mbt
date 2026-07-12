///|
/// Sequence alignment example: global & local pairwise alignment.
///
/// Demonstrates Needleman-Wunsch (global) and Smith-Waterman (local)
/// alignment algorithms for both DNA and protein sequences.
///
/// Run: moon run examples/alignment_demo/main.mbt
fn main {
  // === Global pairwise alignment (Needleman-Wunsch) ===
  let dna1 = @bio.DNA::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
  let dna2 = @bio.DNA::new("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAC")

  let (msa_global, score_g, start_g, end_g) = @bio.global_pairwise_align_nucleotide(
    dna1, dna2,
  )
  println("=== Global Pairwise Alignment (Needleman-Wunsch) ===")
  println("  Score: \{score_g}")
  println("  Start: \{start_g}")
  println("  End:   \{end_g}")
  let recs = msa_global.records
  for i = 0; i < recs.length(); i = i + 1 {
    println("  Seq \{i}: \{recs[i].seq.to_string()}")
  }

  // === Local pairwise alignment (Smith-Waterman) ===
  let (msa_local, score_l, start_l, end_l) = @bio.local_pairwise_align_nucleotide(
    dna1, dna2,
  )
  println("\n=== Local Pairwise Alignment (Smith-Waterman) ===")
  println("  Score: \{score_l}")
  println("  Start: \{start_l}")
  println("  End:   \{end_l}")
  let local_recs = msa_local.records
  for i = 0; i < local_recs.length(); i = i + 1 {
    println("  Seq \{i}: \{local_recs[i].seq.to_string()}")
  }

  // === Protein alignment ===
  let prot1 = @bio.Protein::new("MKKLLLISVLLFLSSAYSR")
  let prot2 = @bio.Protein::new("MKKLLLISVLLFLSSAYSQ")

  let (prot_msa, prot_score, _, _) = @bio.global_pairwise_align_protein(
    prot1, prot2,
  )
  println("\n=== Protein Global Alignment ===")
  println("  Score: \{prot_score}")
  let p_recs = prot_msa.records
  for i = 0; i < p_recs.length(); i = i + 1 {
    println("  Seq \{i}: \{p_recs[i].seq.to_string()}")
  }

  // === Parse Clustal alignment ===
  let clustal = "CLUSTAL X (1.81) multiple sequence alignment\n\n\nseq1                          ATGCCGATCGATCGATCGAT\nseq2                          ATGCCGATCGATCGATCGAT\nseq3                          ATGCCGATCGATCGATCGAT\n"

  try {
    let align = @bio.parse_clustal(clustal)
    println("\n=== Clustal Alignment ===")
    println("  Records: \{align.num_records()}")
    println("  Length:   \{align.get_alignment_length()}")

    // === Write Clustal ===
    let clustal_out = @bio.write_clustal(align)
    println("\n=== Clustal Output ===")
    println(clustal_out)
  } catch {
    _ => println("Clustal error occurred")
  }

  // === Parse Phylip ===
  let phylip = " 3 20\nseq1       ATGCCGATCG ATCGATCGAT \nseq2       ATGCCGATCG ATCGATCGAT \nseq3       ATGCCGATCG ATCGATCGAT \n"
  try {
    let phy_align = @bio.parse_phylip(phylip)
    println("=== Phylip Alignment ===")
    println("  Records: \{phy_align.num_records()}")
  } catch {
    _ => println("Phylip error occurred")
  }

  // === TabularMSA consensus & conservation ===
  let dna_arr = [
    @bio.DNA::new("ATGC"),
    @bio.DNA::new("ATGC"),
    @bio.DNA::new("ATGT"),
  ]
  let tab_msa = @bio.TabularMSA::from_dna(dna_arr)
  let consensus = tab_msa.consensus()
  let conservation = tab_msa.conservation()
  println("\n=== TabularMSA ===")
  println("  Consensus: \{consensus.to_string()}")
  println("  Conservation: \{conservation}")
}
