///|
fn main {
  println("=== Biopython Alignment.map/mapall Demo ===")

  println("\n1. Compose chromosome -> transcript -> read")
  let chromosome_to_transcript = @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA",
    "transcript",
    "CCCCCCCGGGGGG",
    [8, 15, 26, 32],
    [0, 7, 7, 13],
  ) catch {
    _ => abort("failed to build chromosome-to-transcript alignment")
  }
  let transcript_to_read = @src.coordinate_pairwise_alignment(
    "transcript",
    "CCCCCCCGGGGGG",
    "read",
    "CCCCGGGG",
    [3, 11],
    [0, 8],
  ) catch {
    _ => abort("failed to build transcript-to-read alignment")
  }
  let chromosome_to_read = chromosome_to_transcript.map(transcript_to_read) catch {
    _ => abort("failed to compose coordinate alignments")
  }
  println("  " + chromosome_to_read.summary())
  println(
    "  target coordinates: " + chromosome_to_read.target_coordinates.to_string(),
  )
  println(
    "  query coordinates:  " + chromosome_to_read.query_coordinates.to_string(),
  )

  println("\n2. Render the exon blocks and intron gap")
  match chromosome_to_read.aligned_rows() {
    Some((target, query)) => {
      println("  chromosome: " + target)
      println("  read:       " + query)
    }
    None => abort("concrete alignment rows should render")
  }
  let counts = chromosome_to_read.counts()
  println(
    "  blocks=" +
    counts.blocks.to_string() +
    ", aligned=" +
    counts.aligned.to_string() +
    ", query-gap bases=" +
    counts.query_gap_bases.to_string(),
  )

  println("\n3. Query coordinates and serialize PSL")
  println(
    "  chromosome position 14 -> read: " +
    chromosome_to_read.target_to_query(14).to_string(),
  )
  println(
    "  read position 4 -> chromosome: " +
    chromosome_to_read.query_to_target(4).to_string(),
  )
  let psl = chromosome_to_read.to_psl() catch {
    _ => abort("failed to serialize PSL")
  }
  println("  PSL: " + psl)

  println("4. Preserve strand through composition")
  let chromosome_to_reverse_transcript = @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA",
    "reverse_transcript",
    "TCCCCCCGGGGGCCCCCCC",
    [12, 31],
    [19, 0],
  ) catch {
    _ => abort("failed to build reverse transcript alignment")
  }
  let reverse_transcript_to_read = @src.coordinate_pairwise_alignment(
    "reverse_transcript",
    "TCCCCCCGGGGGCCCCCCC",
    "reverse_read",
    "CCCGGGGGCC",
    [4, 14],
    [0, 10],
  ) catch {
    _ => abort("failed to build reverse read alignment")
  }
  let reverse_projection = chromosome_to_reverse_transcript.map(
    reverse_transcript_to_read,
  ) catch {
    _ => abort("failed to compose reverse-strand alignment")
  }
  println("  " + reverse_projection.summary())
  println(
    "  projected query coordinates: " +
    reverse_projection.query_coordinates.to_string(),
  )

  println("\n5. Project a protein MSA to a codon-aware nucleotide MSA")
  let protein_msa = @src.coordinate_multiple_alignment(
    ["protein_1", "protein_2"],
    ["M-K", "MTK"],
  ) catch {
    _ => abort("failed to build protein MSA")
  }
  let protein_1_to_dna = @src.coordinate_pairwise_alignment(
    "protein_1",
    "MK",
    "dna_1",
    "ATGAAA",
    [0, 2],
    [0, 6],
  ) catch {
    _ => abort("failed to build first protein-to-DNA mapping")
  }
  let protein_2_to_dna = @src.coordinate_pairwise_alignment(
    "protein_2",
    "MTK",
    "dna_2",
    "ATGACCAAA",
    [0, 3],
    [0, 9],
  ) catch {
    _ => abort("failed to build second protein-to-DNA mapping")
  }
  let nucleotide_msa = protein_msa.mapall([protein_1_to_dna, protein_2_to_dna]) catch {
    _ => abort("failed to project protein MSA")
  }
  println("  " + nucleotide_msa.summary())
  println(nucleotide_msa.to_fasta())

  println("=== Demo Complete ===")
}
