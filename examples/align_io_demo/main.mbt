///|
/// Alignment format parsing example inspired by Biopython Bio.AlignIO.
///
/// Run: moon run examples/align_io_demo/main.mbt
fn main {
  println("=== Bio.AlignIO Demo ===")

  let clustalw_content = "CLUSTAL W (1.83) multiple sequence alignment\n\nSeq1        ACGTACGT\nSeq2        ACATACAT\nSeq3        AGGTAGGT\n"

  println("\n1. Parse ClustalW alignment:")
  let align = @bio.parse_clustalw(clustalw_content)
  println("   Number of sequences: \{@bio.num_sequences(align)}")
  println("   Alignment length: \{@bio.alignment_length(align)}")

  println("\n2. Get sequence by name:")
  match @bio.get_sequence(align, "Seq1") {
    Some(seq) => println("   Seq1: \{seq}")
    None => println("   Seq1 not found")
  }

  println("\n3. Consensus sequence:")
  let consensus = @bio.consensus_sequence(align)
  println("   Consensus: \{consensus}")

  println("\n4. Column consensus (position 2):")
  let col_cons = @bio.column_consensus(align, 2)
  println("   Position 2: \{col_cons}")

  println("\n5. Write ClustalW format:")
  let clustalw_out = @bio.write_clustalw(align)
  println(clustalw_out)

  let fasta_content = ">Seq1\nACGTACGT\n>Seq2\nACATACAT\n"

  println("\n6. Parse FASTA alignment:")
  let fasta_align = @bio.parse_fasta_alignment(fasta_content)
  println("   Number of sequences: \{@bio.num_sequences(fasta_align)}")

  println("\n7. Write FASTA format:")
  let fasta_out = @bio.write_fasta_alignment(fasta_align)
  println(fasta_out)
}
