///|
fn main {
  println("=== Biopython Bio.Align.a2m Demo ===")

  let alignment = @src.a2m_parse(@src.a2m_example_text()) catch {
    _ => abort("failed to parse A2M sample")
  }

  println("\n1. State-aware A2M parsing")
  println("  " + alignment.summary())
  println("  Column states: " + alignment.states)
  println("  Reference sequence: " + alignment.sequences[0].sequence)

  println("\n2. Insertion slots and match consensus")
  for insertion in alignment.insertion_runs() {
    println(
      "  Slot " +
      insertion.slot.to_string() +
      ": columns [" +
      insertion.start_column.to_string() +
      ", " +
      insertion.end_column.to_string() +
      "), width " +
      insertion.width.to_string(),
    )
  }
  let consensus = alignment.consensus(include_insertions=false) catch {
    _ => abort("failed to calculate A2M consensus")
  }
  println("  Match-state consensus: " + consensus)

  println("\n3. Coordinate mapping")
  match alignment.map_position(0, 1, 5) {
    Some(position) =>
      println(
        "  Reference residue 5 maps to query_one residue " +
        position.to_string() +
        " (zero-based)",
      )
    None => println("  Reference residue 5 maps to a query_one gap")
  }

  println("\n4. Pairwise alignment statistics")
  let counts = alignment.pair_counts(0, 1) catch {
    _ => abort("failed to calculate A2M pair counts")
  }
  println(
    "  Aligned/identity/mismatch: " +
    counts.aligned.to_string() +
    "/" +
    counts.identities.to_string() +
    "/" +
    counts.mismatches.to_string(),
  )
  println("  Identity fraction: " + counts.identity().to_string())

  println("\n5. Match-only projection")
  let projected = alignment.match_projection() catch {
    _ => abort("failed to project A2M match columns")
  }
  println("  " + projected.summary())
  println("  Projected states: " + projected.states)

  println("\n6. Wrapped canonical serialization")
  let serialized = @src.a2m_write(alignment, line_width=6) catch {
    _ => abort("failed to write A2M alignment")
  }
  let reparsed = @src.a2m_parse(serialized) catch {
    _ => abort("failed to reparse A2M alignment")
  }
  println("  Serialized bytes: " + serialized.length().to_string())
  println(
    "  Round trip preserved alignment: " + (alignment == reparsed).to_string(),
  )

  println("\n=== Demo Complete ===")
}
