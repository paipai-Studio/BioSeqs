///|
fn maf_integer_array_text(values : Array[Int]) -> String {
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
  println("=== Biopython Bio.Align.maf Demo ===")

  let document = @src.align_maf_parse(@src.align_maf_example_text()) catch {
    AlignMafError(message) => abort("MAF parsing failed: " + message)
  }
  let track = match document.track {
    Some(value) => value
    None => abort("example MAF track metadata is missing")
  }

  println("\n1. Parse MAF header, track metadata, and a/s/i/e/q records")
  println(
    "  version=" +
    document.version +
    ", scoring=" +
    document.scoring.unwrap_or("none") +
    ", program=" +
    document.program.unwrap_or("none"),
  )
  println(
    "  track=" +
    track.name.unwrap_or("unnamed") +
    ", visibility=" +
    track.visibility.unwrap_or("default") +
    ", species=" +
    track.species_order.length().to_string(),
  )
  println(
    "  blocks=" +
    document.blocks.length().to_string() +
    ", comments=" +
    document.comments.length().to_string(),
  )
  let first = document.blocks[0]
  let reference = first.components[0]
  let insertion = match reference.insertion {
    Some(value) =>
      value.left_status +
      ":" +
      value.left_count.to_string() +
      "/" +
      value.right_status +
      ":" +
      value.right_count.to_string()
    None => "none"
  }
  println(
    "  reference=" +
    reference.source +
    ", sequence=" +
    reference.sequence() +
    ", insertion=" +
    insertion +
    ", quality=" +
    reference.quality.unwrap_or("none"),
  )
  println(
    "  empty component=" +
    first.empty_components[0].source +
    ", status=" +
    first.empty_components[0].status,
  )

  println("\n2. Build absolute coordinate paths across plus/minus strands")
  let path = @src.align_maf_coordinate_path(first)
  for point in path {
    println(
      "  column " +
      point.column.to_string() +
      " -> " +
      maf_integer_array_text(point.positions),
    )
  }
  let to_mouse = @src.align_maf_map_position(
    first, "hg38.chr7", 100, "mm39.chr5",
  ) catch {
    AlignMafError(message) => abort("MAF position mapping failed: " + message)
  }
  let to_human = @src.align_maf_map_position(
    first, "mm39.chr5", 449, "hg38.chr7",
  ) catch {
    AlignMafError(message) => abort("MAF reverse mapping failed: " + message)
  }
  println(
    "  hg38:100 -> mm39:" +
    to_mouse.unwrap_or(-1).to_string() +
    "; mm39:449 -> hg38:" +
    to_human.unwrap_or(-1).to_string(),
  )

  println("\n3. Summarize the multi-block document")
  let summary = document.summary()
  println(
    "  components=" +
    summary.component_count.to_string() +
    ", empty=" +
    summary.empty_component_count.to_string() +
    ", aligned columns=" +
    summary.aligned_columns.to_string(),
  )
  println(
    "  reference bases=" +
    summary.reference_bases.to_string() +
    ", distinct sources=" +
    summary.source_count.to_string(),
  )

  println("\n4. Query a MafIndex-style half-open reference interval")
  let index = @src.AlignMafIndex::create(document, "hg38.chr7") catch {
    AlignMafError(message) => abort("MAF index creation failed: " + message)
  }
  let hits = index.search(103, 111) catch {
    AlignMafError(message) => abort("MAF interval search failed: " + message)
  }
  println(
    "  entries=" +
    index.entry_count().to_string() +
    ", overlaps [103, 111)=" +
    hits.length().to_string(),
  )

  println("\n5. Splice multiple exons with missing-data filling")
  let spliced = index.get_spliced([102, 110], [106, 113]) catch {
    AlignMafError(message) => abort("MAF splicing failed: " + message)
  }
  println("  hg38.chr7: " + spliced.sequence("hg38.chr7").unwrap_or("missing"))
  println("  mm39.chr5: " + spliced.sequence("mm39.chr5").unwrap_or("missing"))
  println("  output columns: " + spliced.columns.to_string())

  println("\n6. Write canonical MAF and parse it again")
  let canonical = @src.align_maf_write(document) catch {
    AlignMafError(message) => abort("MAF serialization failed: " + message)
  }
  let round_trip = @src.align_maf_parse(canonical) catch {
    AlignMafError(message) => abort("MAF round-trip parsing failed: " + message)
  }
  println("  output bytes: " + canonical.length().to_string())
  println("  round-trip preserved: " + (round_trip == document).to_string())
}
