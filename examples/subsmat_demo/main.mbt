///|
/// Substitution matrix demo - demonstrates using BLOSUM and PAM matrices
/// for protein sequence alignment.
///
/// Run with: moon run examples/subsmat_demo/main.mbt
fn main {
  println("=== Substitution Matrix Demo ===\n")

  println("1. BLOSUM62 Matrix Lookup:")
  let blosum62 = @bio.blosum62()
  println("   Matrix name: " + blosum62.name())
  println("   A-A score: " + blosum62.get_score("A", "A").to_string())
  println("   R-R score: " + blosum62.get_score("R", "R").to_string())
  println("   A-R score: " + blosum62.get_score("A", "R").to_string())
  println("   W-W score: " + blosum62.get_score("W", "W").to_string())

  println("\n2. PAM250 Matrix Lookup:")
  let pam250 = @bio.pam250()
  println("   Matrix name: " + pam250.name())
  println("   A-A score: " + pam250.get_score("A", "A").to_string())
  println("   W-W score: " + pam250.get_score("W", "W").to_string())

  println("\n3. Case Insensitivity Test:")
  println(
    "   blosum62.get_score('a', 'A') = " +
    blosum62.get_score_case_insensitive("a", "A").to_string(),
  )

  println("\n4. Protein Alignment with BLOSUM62:")
  let seq1 = "MVKLSQVAPFASG"
  let seq2 = "MVKLSEVAPFSSG"
  let sw_result = @bio.smith_waterman_blosum62(seq1, seq2, -1.0)
  println("   Sequence 1: " + seq1)
  println("   Sequence 2: " + seq2)
  println("   Alignment score: " + sw_result.score().to_string())
  println(
    "   Percent identity: " + sw_result.percent_identity().to_string() + "%",
  )
  println("   Alignment:")
  println("     " + sw_result.aligned1())
  println("     " + sw_result.match_line())
  println("     " + sw_result.aligned2())

  println("\n5. Global Alignment with PAM250:")
  let nw_result = @bio.needleman_wunsch_pam250(seq1, seq2, -2.0)
  println("   Alignment score: " + nw_result.score().to_string())
  println("   Alignment:")
  println("     " + nw_result.aligned1())
  println("     " + nw_result.match_line())
  println("     " + nw_result.aligned2())

  println("\n=== Demo Complete ===")
}
