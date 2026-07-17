///|
/// Smith-Waterman local sequence alignment demo.
///
/// Demonstrates dynamic programming-based local alignment for DNA and protein
/// sequences with custom scoring schemes, score matrix inspection, and
/// formatted output.
///
/// Run: moon run examples/smith_waterman_demo/main.mbt
fn main {
  println("=== Smith-Waterman Local Sequence Alignment ===")
  println("")

  // --- 1. DNA Alignment (exact match) ---
  println("1. DNA Alignment (exact match):")
  let dna1 = "ACGTACGT"
  let dna2 = "ACGTACGT"
  let result1 = @bio.smith_waterman_dna(dna1, dna2)
  println("   Score: \{result1.score()}")
  println("   Identity: \{result1.percent_identity()}%")
  println("   Alignment:")
  println("   \{result1.format().split("\n").to_array()[0]}")
  println("   \{result1.format().split("\n").to_array()[1]}")
  println("   \{result1.format().split("\n").to_array()[2]}")
  println("")

  // --- 2. DNA Alignment (substring) ---
  println("2. DNA Alignment (find substring):")
  let long_seq = "GGTACGTAACTGGTACGATC"
  let short_seq = "ACGT"
  let result2 = @bio.smith_waterman_dna(long_seq, short_seq)
  println("   Seq1: \{long_seq}")
  println("   Seq2: \{short_seq}")
  println("   Score: \{result2.score()}")
  let (s1, e1) = result2.pos1()
  let (s2, e2) = result2.pos2()
  println("   Match in seq1: [\{s1}, \{e1})")
  println("   Match in seq2: [\{s2}, \{e2})")
  println("   Alignment:")
  let lines2 = result2.format().split("\n").to_array()
  println("   \{lines2[0]}")
  println("   \{lines2[1]}")
  println("   \{lines2[2]}")
  println("")

  // --- 3. DNA Alignment with gaps ---
  println("3. DNA Alignment with gaps:")
  let result3 = @bio.smith_waterman_dna("ACGTACGT", "ACGTT")
  println("   Score: \{result3.score()}")
  let lines3 = result3.format().split("\n").to_array()
  println("   \{lines3[0]}")
  println("   \{lines3[1]}")
  println("   \{lines3[2]}")
  println("   Matches: \{result3.num_matches()}")
  println("   Identity: \{result3.percent_identity()}%")
  println("")

  // --- 4. Custom scoring ---
  println("4. Custom Scoring (match=3, mismatch=-2, gap=-2):")
  let scoring = @bio.Scoring::new(3.0, -2.0, -2.0)
  let result4 = @bio.smith_waterman("ACGTACGT", "ACGT", scoring)
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
  let result5 = @bio.smith_waterman_protein(prot1, prot2)
  println("   Seq1: \{prot1}")
  println("   Seq2: \{prot2}")
  println("   Score: \{result5.score()}")
  println("   Alignment length: \{result5.alignment_length()}")
  println("   Alignment:")
  let lines5 = result5.format().split("\n").to_array()
  for i = 0; i < lines5.length(); i = i + 1 {
    println("   \{lines5[i]}")
  }
  println("")

  // --- 6. Protein alignment with BLOSUM62 matrix ---
  println("6. Protein Alignment with BLOSUM62 Matrix:")
  let prot3 = "MVKLSQVAPFASG"
  let prot4 = "MVKLSEVAPFSSG"
  let result6 = @bio.smith_waterman_blosum62(prot3, prot4, -1.0)
  println("   Seq1: \{prot3}")
  println("   Seq2: \{prot4}")
  println("   Score: \{result6.score()}")
  println("   Identity: \{result6.percent_identity()}%")
  println("   Alignment:")
  let lines6 = result6.format().split("\n").to_array()
  for i = 0; i < lines6.length(); i = i + 1 {
    println("   \{lines6[i]}")
  }
  println("")

  // --- 7. Protein alignment with PAM250 matrix ---
  println("7. Protein Alignment with PAM250 Matrix:")
  let prot5 = "ACDEFGHIKLMNPQRSTVWY"
  let prot6 = "ACDEFGHIKLMNPQRSTVWY"
  let result7 = @bio.smith_waterman_pam250(prot5, prot6, -2.0)
  println("   Seq1: \{prot5}")
  println("   Seq2: \{prot6}")
  println("   Score: \{result7.score()}")
  println("   Identity: \{result7.percent_identity()}%")
  println("")

  // --- 8. Long sequence with wrapping ---
  println("8. Long Sequence Alignment (wrapped at 40 chars):")
  let long1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
  let long2 = "GGGGATCGATCGATCGGGGGATCGATCGATCGGGGGATCGATCGATCG"
  let result6 = @bio.smith_waterman_dna(long1, long2)
  println("   Score: \{result6.score()}")
  println("   Matches: \{result6.num_matches()}")
  println("   Identity: \{result6.percent_identity()}%")
  println("   Alignment (wrapped):")
  println(result6.format_wrapped(40))
  println("")

  // --- 9. Score matrix inspection ---
  println("9. Score Matrix (first 5x5):")
  let result7 = @bio.smith_waterman_dna("ACGT", "ACGT")
  let matrix = result7.score_matrix()
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

  // --- 10. No match case ---
  println("10. No Match (AAAA vs TTTT):")
  let result8 = @bio.smith_waterman_dna("AAAA", "TTTT")
  println("   Score: \{result8.score()}")
  println("   Alignment length: \{result8.alignment_length()}")
  println("")

  println("=== Smith-Waterman Demo Complete ===")
}
