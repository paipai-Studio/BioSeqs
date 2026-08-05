///|
fn stockholm_integer_array_text(values : Array[Int]) -> String {
  let output = StringBuilder::new()
  output.write_char('[')
  for index = 0; index < values.length(); index = index + 1 {
    if index > 0 {
      output.write_string(", ")
    }
    output.write_string(values[index].to_string())
  }
  output.write_char(']')
  output.to_string()
}

///|
fn main {
  println("=== Biopython Bio.Align.stockholm Demo ===")

  let alignment = @src.align_stockholm_parse(
    @src.align_stockholm_example_text(),
  ) catch {
    AlignStockholmError(message) =>
      abort("Stockholm parsing failed: " + message)
  }

  println("\n1. Parse the official Pfam HAT record")
  println("  " + alignment.summary())
  println(
    "  family: " +
    alignment.annotation("identifier").unwrap_or("unknown") +
    " (" +
    alignment.annotation("accession").unwrap_or("unknown") +
    ")",
  )
  println(
    "  references=" +
    alignment.references.length().to_string() +
    ", database references=" +
    alignment.database_references.length().to_string(),
  )
  for sequence in alignment.sequences {
    println(
      "  " +
      sequence.id +
      ": " +
      sequence.aligned_sequence +
      " -> " +
      sequence.sequence,
    )
  }

  println("\n2. Inspect insertion operations and coordinates")
  println("  operations: " + alignment.operations)
  let path = alignment.coordinate_path()
  for row = 0; row < alignment.num_sequences(); row = row + 1 {
    println(
      "  " +
      alignment.sequences[row].id +
      ": " +
      stockholm_integer_array_text(path[row]),
    )
  }
  println(
    "  row 3 residue 17 -> row 1 residue " +
    alignment.map_position(2, 0, 17).unwrap_or(-1).to_string() +
    " (-1 denotes a gap)",
  )

  println("\n3. Calculate alignment statistics")
  let counts = alignment.counts()
  println(
    "  pairs=" +
    counts.pairs.to_string() +
    ", aligned=" +
    counts.aligned.to_string() +
    ", identities=" +
    counts.identities.to_string() +
    ", gap columns=" +
    counts.gap_columns.to_string(),
  )
  let consensus = alignment.consensus() catch {
    AlignStockholmError(message) => abort("consensus failed: " + message)
  }
  println("  consensus: " + consensus)
  println(
    "  insertion-column occupancy: " + alignment.occupancy()[17].to_string(),
  )

  println("\n4. Write canonical Stockholm and parse it again")
  let canonical = @src.align_stockholm_write(alignment) catch {
    AlignStockholmError(message) =>
      abort("Stockholm serialization failed: " + message)
  }
  let round_trip = @src.align_stockholm_parse(canonical) catch {
    AlignStockholmError(message) =>
      abort("Stockholm round-trip parsing failed: " + message)
  }
  println("  output bytes: " + canonical.length().to_string())
  println(
    "  rows preserved: " +
    (round_trip.sequences == alignment.sequences).to_string(),
  )
  println(
    "  operations preserved: " +
    (round_trip.operations == alignment.operations).to_string(),
  )

  println("\n5. Remove all-gap source columns during construction")
  let compact = @src.align_stockholm_from_aligned(["alpha", "beta"], [
    "A.-C--G", "AT-C--G",
  ]) catch {
    AlignStockholmError(message) =>
      abort("Stockholm construction failed: " + message)
  }
  println("  " + compact.summary())
  println("  operations: " + compact.operations)
  println("  alpha: " + compact.sequences[0].aligned_sequence)
  println("  beta:  " + compact.sequences[1].aligned_sequence)
}
