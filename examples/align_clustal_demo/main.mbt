///|
fn clustal_demo_parse(text : String) -> @bio.AlignClustalAlignment {
  @bio.align_clustal_parse(text) catch {
    AlignClustalError(message) => abort("CLUSTAL parse failed: " + message)
  }
}

///|
fn clustal_demo_write(alignment : @bio.AlignClustalAlignment) -> String {
  @bio.align_clustal_write(
    alignment,
    block_width=5,
    name_width=16,
    include_counts=true,
    include_consensus=true,
  ) catch {
    AlignClustalError(message) => abort("CLUSTAL write failed: " + message)
  }
}

///|
fn clustal_demo_path(path : Array[Int]) -> String {
  let output = StringBuilder()
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
  println("=== Biopython Bio.Align.clustal Offline Demo ===")
  let alignment = clustal_demo_parse(@bio.align_clustal_example_text())

  println("\n1. Generator metadata and interleaved rows")
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
  println(
    "  stored consensus=" + alignment.consensus.unwrap_or("<not provided>"),
  )
  println("  calculated consensus=" + alignment.calculated_consensus())

  println("\n2. Compressed coordinate path")
  let paths = alignment.coordinate_path()
  for row = 0; row < paths.length(); row = row + 1 {
    println(
      "  " + alignment.sequences[row].id + ": " + clustal_demo_path(paths[row]),
    )
  }
  println(
    "  reference residue 3 -> query_one residue " +
    alignment.map_position(0, 1, 3).unwrap_or(-1).to_string(),
  )

  println("\n3. Pairwise statistics and majority residues")
  let counts = alignment.pair_counts(0, 1) catch {
    AlignClustalError(message) => abort("CLUSTAL counting failed: " + message)
  }
  let majority = alignment.majority_consensus(minimum_fraction=0.5) catch {
    AlignClustalError(message) => abort("CLUSTAL consensus failed: " + message)
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
  println("  majority=" + majority)

  println("\n4. Canonical writing with cumulative residue counts")
  let canonical = clustal_demo_write(alignment)
  let reparsed = clustal_demo_parse(canonical)
  println(
    "  canonical generator=" +
    reparsed.metadata.program +
    " " +
    reparsed.metadata.version,
  )
  println(
    "  stable round-trip=" +
    (clustal_demo_write(reparsed) == canonical).to_string(),
  )
  println(
    "\nThe demo parses and writes existing artifacts; it does not launch an aligner.",
  )
}
