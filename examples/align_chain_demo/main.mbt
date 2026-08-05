///|
fn chain_coordinate_path_text(
  coordinates : Array[@src.AlignChainCoordinate],
) -> String {
  let output = StringBuilder::new()
  for index = 0; index < coordinates.length(); index = index + 1 {
    if index > 0 {
      output.write_string(" -> ")
    }
    output.write_string(
      "(" +
      coordinates[index].target.to_string() +
      ", " +
      coordinates[index].query.to_string() +
      ")",
    )
  }
  output.to_string()
}

///|
fn main {
  println("=== Biopython Bio.Align.chain Demo ===")

  let alignments = @src.align_chain_parse_all(@src.align_chain_example_text()) catch {
    AlignChainError(message) => abort("Chain parsing failed: " + message)
  }
  let alignment = alignments[0]

  println("\n1. Parse consecutive UCSC Chain records")
  println("  records: " + alignments.length().to_string())
  for record in alignments {
    println("  " + record.summary())
  }

  println("\n2. Inspect absolute coordinates and canonical blocks")
  println("  operations: " + alignment.operation_path())
  println("  coordinates: " + chain_coordinate_path_text(alignment.coordinates))
  let blocks = alignment.blocks() catch {
    AlignChainError(message) =>
      abort("Chain block reconstruction failed: " + message)
  }
  for index = 0; index < blocks.length(); index = index + 1 {
    let block = blocks[index]
    println(
      "  block " +
      (index + 1).to_string() +
      ": size=" +
      block.size.to_string() +
      ", dt=" +
      block.target_gap.to_string() +
      ", dq=" +
      block.query_gap.to_string(),
    )
  }

  println("\n3. Calculate coordinate-only counts")
  let counts = alignment.counts()
  println(
    "  aligned=" +
    counts.aligned.to_string() +
    ", target-gap bases=" +
    counts.target_gap_bases.to_string() +
    ", query-gap bases=" +
    counts.query_gap_bases.to_string(),
  )
  println(
    "  alignment length=" +
    alignment.alignment_length().to_string() +
    ", aligned blocks=" +
    counts.aligned_blocks.to_string(),
  )

  println("\n4. Map residues and intervals in both directions")
  println(
    "  target 42530895 -> query " +
    alignment.map_target_position(42530895).unwrap_or(-1).to_string(),
  )
  println(
    "  query 180 -> target " +
    alignment.map_query_position(180).unwrap_or(-1).to_string(),
  )
  let pieces = alignment.map_target_range(42530950, 42532030) catch {
    AlignChainError(message) => abort("Chain range mapping failed: " + message)
  }
  for piece in pieces {
    println(
      "  target[" +
      piece.target_start.to_string() +
      ", " +
      piece.target_end.to_string() +
      ") -> query[" +
      piece.query_start.to_string() +
      ", " +
      piece.query_end.to_string() +
      ")",
    )
  }

  println("\n5. Invert and round-trip canonical Chain output")
  let inverted = alignment.invert() catch {
    AlignChainError(message) => abort("Chain inversion failed: " + message)
  }
  println("  inverted: " + inverted.summary())
  let canonical = @src.align_chain_write_all(alignments) catch {
    AlignChainError(message) => abort("Chain serialization failed: " + message)
  }
  let round_trip = @src.align_chain_parse_all(canonical) catch {
    AlignChainError(message) =>
      abort("Chain round-trip parsing failed: " + message)
  }
  println("  output bytes: " + canonical.length().to_string())
  println("  round-trip preserved: " + (round_trip == alignments).to_string())
}
