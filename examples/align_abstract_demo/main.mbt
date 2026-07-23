///|
/// AlignAbstract Demo - Demonstrates abstract alignment types and validation.
/// Inspired by Biopython's Bio.Align.AlignAbstract module.
///
/// Run with: moon run examples/align_abstract_demo/main.mbt
fn main {
  println("=== AlignAbstract Demo ===\n")

  // 1. Create an alignment
  println("1. Creating nucleotide alignment:")
  let sequences = ["ACGTACGT", "ACGTACGA", "ACGTACGT"]
  let identifiers = ["seq1", "seq2", "seq3"]
  let aln = @bio.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@bio.align_abstract_type_from_str("nucleotide"),
  )
  println("   Sequences: \{aln.abstract_n_seqs()}")
  println("   Length: \{aln.abstract_alignment_length()}")
  println("   Type: \{@bio.alignment_type_str(aln.abstract_type())}")

  // 2. Validate alignment
  println("\n2. Validation:")
  let (valid, msg) = aln.validate()
  println("   Valid: \{valid}, Message: \{msg}")

  // 3. Consensus sequence
  println("\n3. Consensus sequence:")
  let consensus = @bio.abstract_consensus_sequence(aln)
  println("   Consensus: \{consensus}")
  println("   Consensus length: \{consensus.length()}")

  // 4. Column statistics
  println("\n4. Column statistics:")
  let col0_stats = @bio.abstract_column_stats(aln, 0)
  println("   Column 0 - conservation: \{col0_stats.conservation}, diversity: \{col0_stats.diversity}")
  println("   Column 0 - gap_fraction: \{col0_stats.gap_fraction}, unique_chars: \{col0_stats.n_unique_chars}")

  // 5. All column stats
  println("\n5. All column statistics:")
  let all_stats = @bio.abstract_all_column_stats(aln)
  let mut i = 0
  while i < all_stats.length() {
    let s = all_stats[i]
    println("   Col \{i}: conservation=\{s.conservation}, unique=\{s.n_unique_chars}")
    i = i + 1
  }

  // 6. Identity matrix
  println("\n6. Identity matrix:")
  let id_matrix = @bio.abstract_identity_matrix(aln)
  let mut r = 0
  while r < id_matrix.length() {
    let row = id_matrix[r]
    let mut c = 0
    let mut row_str = "   "
    while c < row.length() {
      row_str = row_str + row[c].to_string() + " "
      c = c + 1
    }
    println(row_str)
    r = r + 1
  }

  // 7. Overall identity
  println("\n7. Overall identity:")
  let overall_id = @bio.abstract_overall_identity(aln)
  println("   Overall identity: \{overall_id}")

  // 8. Coverage
  println("\n8. Coverage:")
  let coverage = @bio.abstract_coverage(aln)
  println("   Alignment coverage: \{coverage}")

  // 9. Variable and parsimony sites
  println("\n9. Informative sites:")
  let variable = @bio.abstract_variable_sites(aln)
  println("   Variable sites: \{variable}")
  let parsimony = @bio.abstract_parsimony_sites(aln)
  println("   Parsimony-informative sites: \{parsimony}")
  let singleton = @bio.abstract_singleton_sites(aln)
  println("   Singleton sites: \{singleton}")

  // 10. Distance matrix
  println("\n10. Distance matrix (p-distance):")
  let dist_matrix = @bio.abstract_distance_matrix(aln)
  let mut dr = 0
  while dr < dist_matrix.length() {
    let row = dist_matrix[dr]
    let mut dc = 0
    let mut row_str = "   "
    while dc < row.length() {
      row_str = row_str + row[dc].to_string() + " "
      dc = dc + 1
    }
    println(row_str)
    dr = dr + 1
  }

  // 11. Filter gaps
  println("\n11. Filter columns by gap fraction:")
  let gapped_seqs = ["A-CGT", "A-CGA", "A-CGT"]
  let gapped_ids = ["seq1", "seq2", "seq3"]
  let gapped_aln = @bio.AbstractAlignment::new(
    sequences=gapped_seqs,
    identifiers=gapped_ids,
    alignment_type=@bio.align_abstract_type_from_str("nucleotide"),
  )
  let filtered = @bio.abstract_filter_gaps(gapped_aln, 0.5)
  println("   Original length: \{gapped_aln.abstract_alignment_length()}")
  println("   Filtered length: \{filtered.abstract_alignment_length()}")

  // 12. Trim alignment
  println("\n12. Trim alignment:")
  let trimmed = @bio.abstract_trim(aln, 2, 5)
  println("   Original length: \{aln.abstract_alignment_length()}")
  println("   Trimmed length: \{trimmed.abstract_alignment_length()}")

  // 13. Summary
  println("\n13. Alignment summary:")
  let summary = @bio.abstract_summary(aln)
  println(summary)

  // 14. Protein alignment
  println("14. Protein alignment validation:")
  let prot_seqs = ["ALWMRLLPLL", "ALWMRLLPLL", "ALWLRLLPLL"]
  let prot_ids = ["prot1", "prot2", "prot3"]
  let prot_aln = @bio.AbstractAlignment::new(
    sequences=prot_seqs,
    identifiers=prot_ids,
    alignment_type=@bio.align_abstract_type_from_str("protein"),
  )
  let (prot_valid, prot_msg) = prot_aln.validate()
  println("   Valid: \{prot_valid}, Message: \{prot_msg}")
  println("   Consensus: \{@bio.abstract_consensus_sequence(prot_aln)}")

  println("\n=== Demo Complete ===")
}