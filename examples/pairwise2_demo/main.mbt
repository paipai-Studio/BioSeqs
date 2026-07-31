///|
/// Bio.Pairwise2 module example - Flexible pairwise sequence alignment.
///
/// Run: moon run examples/pairwise2_demo/main.mbt
fn main {
  println("=== Bio.Pairwise2 Module Demo ===")
  println("")

  // Global alignment (identical sequences)
  println("=== Global Alignment (Identical Sequences) ===")
  let r1 = @bio.pairwise_globalxx("ACGT", "ACGT", 2.0, -1.0, -5.0, -0.5)
  println("Alignment: \{r1.get_aligned_seq1()} vs \{r1.get_aligned_seq2()}")
  println("Score: \{r1.get_score()}")
  println("Mode: \{r1.get_mode()}")
  println("")

  // Global alignment (with gaps)
  println("=== Global Alignment (With Gaps) ===")
  let r2 = @bio.pairwise_globalxx("ACGT", "ACG", 2.0, -1.0, -5.0, -0.5)
  println("Alignment: \{r2.get_aligned_seq1()} vs \{r2.get_aligned_seq2()}")
  println("Score: \{r2.get_score()}")
  println("")

  // Local alignment (identical sequences)
  println("=== Local Alignment (Identical Sequences) ===")
  let r3 = @bio.pairwise_localxx("ACGTACGT", "ACGTACGT", 2.0, -1.0, -5.0, -0.5)
  println("Alignment: \{r3.get_aligned_seq1()} vs \{r3.get_aligned_seq2()}")
  println("Score: \{r3.get_score()}")
  println("")

  // Local alignment (substring)
  println("=== Local Alignment (Substring) ===")
  let r4 = @bio.pairwise_localxx("XXXACGTXXX", "ACGT", 2.0, -1.0, -5.0, -0.5)
  println("Alignment: \{r4.get_aligned_seq1()} vs \{r4.get_aligned_seq2()}")
  println("Score: \{r4.get_score()}")
  println("")

  // Local alignment (no match)
  println("=== Local Alignment (No Match) ===")
  let r5 = @bio.pairwise_localxx("AAAA", "TTTT", 2.0, -1.0, -5.0, -0.5)
  println("Score: \{r5.get_score()}")
  println("Found match: \{r5.get_score() > 0.0}")
  println("")

  // Convenience functions
  println("=== Convenience: Global Alignment ===")
  let r6 = @bio.pairwise_global("ACGT", "ACGT")
  println("Alignment: \{r6.get_aligned_seq1()} vs \{r6.get_aligned_seq2()}")
  println("Score: \{r6.get_score()}")
  println("")

  println("=== Convenience: Local Alignment ===")
  let r7 = @bio.pairwise_local("ACGTACGT", "ACGT")
  println("Alignment: \{r7.get_aligned_seq1()} vs \{r7.get_aligned_seq2()}")
  println("Score: \{r7.get_score()}")
  println("")

  // Scoring functions
  println("=== Simple Scoring Function ===")
  let scorer = @bio.simple_score(2.0, -1.0)
  println("Match(A,A): \{scorer("A", "A")}")
  println("Mismatch(A,T): \{scorer("A", "T")}")
  println("Case insensitive(A,a): \{scorer("A", "a")}")
  println("")

  println("=== Identity Scoring Function ===")
  let identity = @bio.identity_score()
  println("Match(A,A): \{identity("A", "A")}")
  println("Mismatch(A,T): \{identity("A", "T")}")
  println("")

  // DNA matrix
  println("=== DNA Scoring Matrix ===")
  let dna_mat = @bio.dna_matrix(2.0, -1.0)
  println("Matrix size: \{dna_mat.size()}")
  println("AA score: \{dna_mat.get("AA").unwrap()}")
  println("AT score: \{dna_mat.get("AT").unwrap()}")
  println("")

  // Matrix-based scoring
  println("=== Matrix-Based Scoring ===")
  let mat_scorer = @bio.matrix_score(dna_mat, 0.0)
  println("Match(A,A): \{mat_scorer("A", "A")}")
  println("Mismatch(A,T): \{mat_scorer("A", "T")}")
  println("Unknown(X,X): \{mat_scorer("X", "X")}")
  println("")

  // GlobalMS alignment
  println("=== GlobalMS Alignment ===")
  let r8 = @bio.pairwise_globalms("ACGT", "ACGT", dna_mat, -5.0, -0.5)
  println("Alignment: \{r8.get_aligned_seq1()} vs \{r8.get_aligned_seq2()}")
  println("Score: \{r8.get_score()}")
  println("Mode: \{r8.get_mode()}")
  println("")

  // LocalMS alignment
  println("=== LocalMS Alignment ===")
  let r9 = @bio.pairwise_localms("XXXACGTXXX", "ACGT", dna_mat, -5.0, -0.5)
  println("Alignment: \{r9.get_aligned_seq1()} vs \{r9.get_aligned_seq2()}")
  println("Score: \{r9.get_score()}")
  println("Mode: \{r9.get_mode()}")
  println("")

  // Alignment summary
  println("=== Alignment Summary ===")
  let summary = @bio.alignment_summary(r1)
  println(summary)

  // Formatted alignment
  println("=== Formatted Alignment ===")
  let formatted = r1.format_alignment(60)
  println(formatted)

  println("=== Done ===")
}