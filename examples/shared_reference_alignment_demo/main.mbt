///|
fn main {
  println("=== Biopython Shared-reference Alignment Demo ===")

  println("\n1. Adapt an existing pairwise alignment")
  let pairwise = @src.pairaligner_align("ACGT", "ACT")
  let pairwise_input = @src.shared_reference_input_from_pairwise(
    pairwise,
    reference_name="reference",
    reference_description="shared DNA template",
    query_name="pairwise_read",
    query_description="one deletion",
  ) catch {
    _ => abort("failed to adapt pairwise alignment")
  }
  println("  reference: " + pairwise_input.aligned_reference)
  println("  query:     " + pairwise_input.queries[0].aligned_sequence)

  println("\n2. Build a multi-query alignment with an insertion")
  let inserted = @src.shared_reference_sequence(
    "inserted_read",
    "ACGGT",
    description="one insertion",
  ) catch {
    _ => abort("failed to build inserted query")
  }
  let deleted = @src.shared_reference_sequence(
    "deleted_read",
    "A---T",
    description="three deletions",
  ) catch {
    _ => abort("failed to build deleted query")
  }
  let multiple_input = @src.shared_reference_input(
    "ACGT",
    "ACG-T",
    [inserted, deleted],
    reference_name="reference",
    reference_description="shared DNA template",
  ) catch {
    _ => abort("failed to build multi-query input")
  }

  println("\n3. Synchronize reference-boundary insertion slots")
  let alignment = @src.alignments_with_same_reference([
    pairwise_input, multiple_input,
  ]) catch {
    _ => abort("failed to merge shared-reference alignments")
  }
  println("  " + alignment.summary())
  println("  insertion widths: " + alignment.insertion_widths.to_string())
  for index in 0..<alignment.num_sequences() {
    let name = alignment.row_name(index).unwrap()
    let row = alignment.row(index).unwrap()
    println("  " + name + ": " + row)
  }

  println("\n4. Map query and reference coordinates")
  println(
    "  inserted_read query position 3 -> reference: " +
    alignment.query_to_reference(1, 3).to_string(),
  )
  println(
    "  inserted_read query position 4 -> reference: " +
    alignment.query_to_reference(1, 4).to_string(),
  )
  println(
    "  reference position 3 -> inserted_read query: " +
    alignment.reference_to_query(1, 3).to_string(),
  )

  println("\n5. Inspect alignment statistics")
  for index in 0..<alignment.num_queries() {
    let query = alignment.queries[index]
    let counts = alignment.counts(index).unwrap()
    println(
      "  " +
      query.name +
      ": identity=" +
      counts.identity_percent().to_string() +
      "%, mismatches=" +
      counts.mismatches.to_string() +
      ", insertions=" +
      counts.insertions.to_string() +
      ", deletions=" +
      counts.deletions.to_string(),
    )
  }

  println("\n6. Convert to MultipleSeqAlignment and aligned FASTA")
  let msa = alignment.to_msa() catch {
    _ => abort("failed to convert merged alignment")
  }
  println(
    "  MSA dimensions: " +
    msa.num_records().to_string() +
    " x " +
    msa.get_alignment_length().to_string(),
  )
  println(alignment.to_fasta())

  println("=== Demo Complete ===")
}
