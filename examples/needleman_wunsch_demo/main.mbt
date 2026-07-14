///|
/// Needleman-Wunsch global sequence alignment demo.
///
/// Demonstrates dynamic programming-based global alignment for DNA and
/// protein sequences with custom scoring schemes, score matrix inspection,
/// and formatted output.
///
/// Run: moon run examples/needleman_wunsch_demo/main.mbt
fn main {
  println("=== Needleman-Wunsch Global Sequence Alignment ===")
  println("")

  // --- 1. DNA Alignment (identical) ---
  println("1. DNA Alignment (identical sequences):")
  let dna1 = "ACGTACGT"
  let dna2 = "ACGTACGT"
  let result1 = @bio.needleman_wunsch_dna(dna1, dna2)
  println("   Seq1: \{dna1}")
  println("   Seq2: \{dna2}")
  println("   Score: \{result1.score()}")
  println("   Identity: \{result1.percent_identity()}%")
  println("   Alignment:")
  let lines1 = result1.format().split("\n").to_array()
  println("   \{lines1[0]}")
  println("   \{lines1[1]}")
  println("   \{lines1[2]}")
  println("")

  // --- 2. DNA Alignment with gaps ---
  println("2. DNA Alignment with gaps:")
  let result2 = @bio.needleman_wunsch_dna("ACGTACGT", "ACGTT")
  println("   Seq1: ACGTACGT")
  println("   Seq2: ACGTT")
  println("   Score: \{result2.score()}")
  println("   Matches: \{result2.num_matches()}")
  println("   Mismatches: \{result2.num_mismatches()}")
  println("   Gaps in seq1: \{result2.num_gaps_seq1()}")
  println("   Gaps in seq2: \{result2.num_gaps_seq2()}")
  println("   Identity: \{result2.percent_identity()}%")
  let lines2 = result2.format().split("\n").to_array()
  println("   \{lines2[0]}")
  println("   \{lines2[1]}")
  println("   \{lines2[2]}")
  println("")

  // --- 3. DNA Alignment (different lengths) ---
  println("3. DNA Alignment (different lengths):")
  let result3 = @bio.needleman_wunsch_dna("GATTACA", "GCATGCU")
  println("   Seq1: GATTACA")
  println("   Seq2: GCATGCU")
  println("   Score: \{result3.score()}")
  let lines3 = result3.format().split("\n").to_array()
  println("   \{lines3[0]}")
  println("   \{lines3[1]}")
  println("   \{lines3[2]}")
  println("")

  // --- 4. Custom scoring ---
  println("4. Custom Scoring (match=3, mismatch=-2, gap=-2):")
  let scoring = @bio.NWScoring::new(3.0, -2.0, -2.0)
  let result4 = @bio.needleman_wunsch("ACGTACGT", "ACGT", scoring)
  println("   Score: \{result4.score()}")
  let lines4 = result4.format().split("\n").to_array()
  println("   \{lines4[0]}")
  println("   \{lines4[1]}")
  println("   \{lines4[2]}")
  println("")

  // --- 5. Protein alignment ---
  println("5. Protein Alignment:")
  let prot1 = "HEAGAWGHEE"
  let prot2 = "PAWHEAE"
  let result5 = @bio.needleman_wunsch_protein(prot1, prot2)
  println("   Seq1: \{prot1}")
  println("   Seq2: \{prot2}")
  println("   Score: \{result5.score()}")
  println("   Alignment length: \{result5.alignment_length()}")
  println("   Matches: \{result5.num_matches()}")
  println("   Identity: \{result5.percent_identity()}%")
  let lines5 = result5.format().split("\n").to_array()
  for i = 0; i < lines5.length(); i = i + 1 {
    println("   \{lines5[i]}")
  }
  println("")

  // --- 6. Completely different sequences ---
  println("6. Completely Different (AAAA vs TTTT):")
  let result6 = @bio.needleman_wunsch_dna("AAAA", "TTTT")
  println("   Score: \{result6.score()}")
  println("   Note: NW forces global alignment, so score is negative")
  let lines6 = result6.format().split("\n").to_array()
  println("   \{lines6[0]}")
  println("   \{lines6[1]}")
  println("   \{lines6[2]}")
  println("")

  // --- 7. Long sequence with wrapping ---
  println("7. Long Sequence Alignment (wrapped at 40 chars):")
  let long1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
  let long2 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
  let result7 = @bio.needleman_wunsch_dna(long1, long2)
  println("   Score: \{result7.score()}")
  println("   Matches: \{result7.num_matches()}")
  println("   Identity: \{result7.percent_identity()}%")
  println("   Alignment (wrapped):")
  println(result7.format_wrapped(40))
  println("")

  // --- 8. Score matrix inspection ---
  println("8. Score Matrix (5x5):")
  let result8 = @bio.needleman_wunsch_dna("ACGT", "ACGT")
  let matrix = result8.score_matrix()
  let rows = if matrix.length() > 5 { 5 } else { matrix.length() }
  for i = 0; i < rows; i = i + 1 {
    let row = matrix[i]
    let cols = if row.length() > 5 { 5 } else { row.length() }
    let buf = StringBuilder::new()
    buf.write_string("   [")
    for j = 0; j < cols; j = j + 1 {
      if j > 0 {
        buf.write_string(", ")
      }
      buf.write_string(row[j].to_string())
    }
    buf.write_string("]")
    println(buf.to_string())
  }
  println("")

  // --- 9. Empty sequence ---
  println("9. Empty Sequence (ACGT vs empty):")
  let result9 = @bio.needleman_wunsch_dna("ACGT", "")
  println("   Score: \{result9.score()}")
  let lines9 = result9.format().split("\n").to_array()
  println("   \{lines9[0]}")
  println("   \{lines9[1]}")
  println("   \{lines9[2]}")
  println("")

  println("=== Needleman-Wunsch Demo Complete ===")
}
