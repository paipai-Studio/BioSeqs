///|
fn align_mauve_demo_document() -> @bio.AlignMauveDocument {
  @bio.align_mauve_parse(@bio.align_mauve_example_text()) catch {
    AlignMauveError(message) =>
      abort("failed to parse XMFA example: " + message)
  }
}

///|
fn align_mauve_demo_path(values : Array[Int]) -> String {
  "[" + values.map(fn(value) { value.to_string() }).join(", ") + "]"
}

///|
fn main {
  println("=== Biopython Bio.Align.mauve XMFA Demo ===")
  let document = align_mauve_demo_document()
  let summary = @bio.align_mauve_summary(document)
  println("\n1. Strict document and locally collinear blocks")
  println(
    "  format=" +
    document.format_version() +
    ", identifiers=" +
    document.identifiers().join(", "),
  )
  println(
    "  sources=" +
    summary.sources().to_string() +
    ", blocks=" +
    summary.blocks().to_string() +
    ", rows=" +
    summary.rows().to_string(),
  )

  let first = document.blocks()[0]
  let compact = first.compact_coordinates()
  println("\n2. Forward and reverse coordinate paths")
  println(
    "  first block: " +
    first.row_count().to_string() +
    " rows x " +
    first.width().to_string() +
    " columns",
  )
  println("  reverse row compact path: " + align_mauve_demo_path(compact[0]))
  println("  forward row compact path: " + align_mauve_demo_path(compact[2]))

  let counts = first.counts()
  println("\n3. Pairwise alignment statistics")
  println(
    "  aligned=" +
    counts.aligned().to_string() +
    ", identities=" +
    counts.identities().to_string() +
    ", mismatches=" +
    counts.mismatches().to_string() +
    ", gaps=" +
    counts.gaps().to_string() +
    ", gap opens=" +
    counts.gap_opens().to_string(),
  )

  let index = @bio.align_mauve_index(document)
  let overlaps = index.query("0", 1, 10) catch {
    AlignMauveError(message) => abort("interval query failed: " + message)
  }
  let positions = @bio.align_mauve_map_position(document, "0", "2", 48) catch {
    AlignMauveError(message) => abort("position projection failed: " + message)
  }
  println("\n4. Interval index and cross-sequence projection")
  println("  source 0 [1,10) overlaps: " + overlaps.length().to_string())
  if positions.length() > 0 {
    println(
      "  source 0 position 48 -> source 2 position " +
      positions[0].target_position().to_string() +
      " in block " +
      positions[0].block_index().to_string(),
    )
  }

  let reconstructed = @bio.align_mauve_reconstruct_sequence(document, "0") catch {
    AlignMauveError(message) =>
      abort("sequence reconstruction failed: " + message)
  }
  println("\n5. Forward source reconstruction")
  println(
    "  source 0 length=" +
    reconstructed.length().to_string() +
    ", sequence=" +
    reconstructed,
  )

  let written = @bio.align_mauve_write(document, line_width=20) catch {
    AlignMauveError(message) => abort("XMFA write failed: " + message)
  }
  let reparsed = @bio.align_mauve_parse(written) catch {
    AlignMauveError(message) => abort("XMFA round-trip failed: " + message)
  }
  println("\n6. Canonical wrapped writer")
  println(
    "  round-trip blocks=" +
    reparsed.block_count().to_string() +
    ", equivalent=" +
    (reparsed == document).to_string(),
  )
}
