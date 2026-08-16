///|
fn phylip_demo_parse(text : String) -> @bio.AlignPhylipAlignment {
  @bio.align_phylip_parse(text) catch {
    AlignPhylipError(message) => abort("PHYLIP parse failed: " + message)
  }
}

///|
fn phylip_demo_write(
  alignment : @bio.AlignPhylipAlignment,
  layout : @bio.AlignPhylipLayout,
) -> String {
  @bio.align_phylip_write(alignment, layout~, block_width=4, group_width=0) catch {
    AlignPhylipError(message) => abort("PHYLIP write failed: " + message)
  }
}

///|
fn phylip_demo_path(path : Array[Int]) -> String {
  let output = StringBuilder::new()
  output.write_char('[')
  for index = 0; index < path.length(); index = index + 1 {
    if index > 0 {
      output.write_string(", ")
    }
    output.write_string(path[index].to_string())
  }
  output.write_char(']')
  output.to_string()
}

///|
fn main {
  println("=== Biopython Bio.Align.phylip Offline Demo ===")
  let alignment = phylip_demo_parse(@bio.align_phylip_example_text())

  println("\n1. Automatic layout detection and fixed-width rows")
  println("  " + alignment.summary())
  for sequence in alignment.sequences {
    println(
      "  " +
      sequence.id +
      ": aligned=" +
      sequence.aligned_sequence +
      ", ungapped=" +
      sequence.sequence,
    )
  }

  println("\n2. Biopython-compatible identifier normalization")
  let normalized = @bio.align_phylip_normalize_id("Taxon[1]:extra") catch {
    AlignPhylipError(message) =>
      abort("PHYLIP name normalization failed: " + message)
  }
  println("  Taxon[1]:extra -> " + normalized)

  println("\n3. Coordinate path and cross-row projection")
  let paths = alignment.coordinate_path()
  for row = 0; row < paths.length(); row = row + 1 {
    println(
      "  " + alignment.sequences[row].id + ": " + phylip_demo_path(paths[row]),
    )
  }
  println(
    "  reference residue 3 -> query_one residue " +
    alignment.map_position(0, 1, 3).unwrap_or(-1).to_string(),
  )

  println("\n4. Pairwise statistics and majority consensus")
  let counts = alignment.pair_counts(0, 1) catch {
    AlignPhylipError(message) => abort("PHYLIP counting failed: " + message)
  }
  let consensus = alignment.majority_consensus(minimum_fraction=0.5) catch {
    AlignPhylipError(message) => abort("PHYLIP consensus failed: " + message)
  }
  println(
    "  aligned=" +
    counts.aligned.to_string() +
    ", identities=" +
    counts.identities.to_string() +
    ", gap opens=" +
    counts.gap_opens.to_string() +
    ", identity=" +
    counts.identity().to_string(),
  )
  println("  majority=" + consensus)

  println("\n5. Canonical sequential and interleaved writing")
  let sequential = phylip_demo_write(alignment, @bio.PhylipSequential)
  let interleaved = phylip_demo_write(alignment, @bio.PhylipInterleaved)
  let reparsed = phylip_demo_parse(interleaved)
  let sequential_lines = sequential.split("\n").to_array()
  println("  sequential first row=" + sequential_lines[1].to_owned())
  println(
    "  stable interleaved round-trip=" +
    (phylip_demo_write(reparsed, @bio.PhylipInterleaved) == interleaved).to_string(),
  )
  println(
    "\nThe demo parses and writes existing artifacts; it does not launch an aligner.",
  )
}
