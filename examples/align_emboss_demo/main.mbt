///|
fn coordinate_text(values : Array[Int]) -> String {
  let output = StringBuilder()
  output.write_char('[')
  for index in 0..<values.length() {
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
  println("=== Biopython Bio.Align.emboss Demo ===")

  let document = @src.align_emboss_parse(@src.align_emboss_example_text()) catch {
    AlignEmbossError(message) => abort("EMBOSS parsing failed: " + message)
  }
  let alignment = document.alignments[0]

  println("\n1. Parse report and alignment metadata")
  println("  " + document.summary())
  println("  " + alignment.summary())
  println(
    "  matrix=" +
    alignment.annotations.matrix +
    ", score=" +
    alignment.annotations.score.unwrap_or(0.0).to_string(),
  )

  println("\n2. Preserve local coordinates and consensus")
  for sequence in alignment.sequences {
    println(
      "  " +
      sequence.id +
      ": boundaries [" +
      sequence.start.to_string() +
      ", " +
      sequence.end.to_string() +
      "), residues=" +
      sequence.sequence.length().to_string(),
    )
  }
  println("  consensus: " + alignment.consensus)

  println("\n3. Map absolute residue coordinates")
  match alignment.map_position(0, 1, 87) {
    Some(position) =>
      println(
        "  reference residue 87 maps to query residue " +
        position.to_string() +
        " (zero-based absolute coordinates)",
      )
    None => println("  reference residue 87 maps to a query gap")
  }
  let path = alignment.coordinate_path()
  println("  reference path: " + coordinate_text(path[0]))
  println("  query path:     " + coordinate_text(path[1]))

  println("\n4. Calculate pairwise EMBOSS statistics")
  let counts = alignment.pair_counts(0, 1) catch {
    AlignEmbossError(message) => abort("pair counting failed: " + message)
  }
  println(
    "  identities/mismatches/gaps: " +
    counts.identities.to_string() +
    "/" +
    counts.mismatches.to_string() +
    "/" +
    (counts.insertions + counts.deletions).to_string(),
  )
  println(
    "  positives=" +
    counts.positives.to_string() +
    ", identity=" +
    counts.identity().to_string(),
  )

  println("\n5. Handle reverse-strand boundaries")
  let annotations = @src.AlignEmbossAnnotations::create(4) catch {
    AlignEmbossError(message) =>
      abort("annotation construction failed: " + message)
  }
  let forward = @src.AlignEmbossSequence::create("target", "ACGT", 2, 6) catch {
    AlignEmbossError(message) =>
      abort("forward row construction failed: " + message)
  }
  let reverse = @src.AlignEmbossSequence::create("reverse", "ACGT", 20, 16) catch {
    AlignEmbossError(message) =>
      abort("reverse row construction failed: " + message)
  }
  let reverse_alignment = @src.AlignEmbossAlignment::create(
    [forward, reverse],
    annotations,
    consensus="||||",
  ) catch {
    AlignEmbossError(message) =>
      abort("reverse alignment construction failed: " + message)
  }
  println(
    "  reverse columns 0..3 map to residues " +
    reverse_alignment
    .column_to_sequence_position(1, 0)
    .unwrap_or(-1)
    .to_string() +
    ".." +
    reverse_alignment
    .column_to_sequence_position(1, 3)
    .unwrap_or(-1)
    .to_string(),
  )

  println("\n6. Wrapped canonical serialization")
  let serialized = @src.align_emboss_write(document, line_width=7) catch {
    AlignEmbossError(message) => abort("EMBOSS writing failed: " + message)
  }
  let reparsed = @src.align_emboss_parse(serialized) catch {
    AlignEmbossError(message) =>
      abort("serialized EMBOSS parsing failed: " + message)
  }
  println("  serialized bytes: " + serialized.length().to_string())
  println(
    "  round trip preserved report: " + (document == reparsed).to_string(),
  )

  println("\n=== Demo Complete ===")
}
