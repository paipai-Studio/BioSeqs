///|
fn parse_sam(text : String) -> @src.AlignSamDocument {
  @src.align_sam_parse(text) catch {
    AlignSamError(message) => abort("SAM parsing failed: " + message)
  }
}

///|
fn demo_reference() -> String {
  "N".repeat(4) + "ACGTTCGGAA" + "N".repeat(26)
}

///|
fn optional_int(value : Int?) -> String {
  match value {
    Some(item) => item.to_string()
    None => "unknown"
  }
}

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
  println("=== Biopython Bio.Align.sam Demo ===")

  println("\n1. Parse headers, references, and records")
  let document = parse_sam(@src.align_sam_example_data())
  println(
    "  headers=" +
    document.headers.length().to_string() +
    ", references=" +
    document.references.length().to_string() +
    ", alignments=" +
    document.alignments.length().to_string(),
  )
  match document.reference("chr1") {
    Some(reference) =>
      println(
        "  reference " +
        reference.name +
        ", length=" +
        reference.length.to_string(),
      )
    None => abort("chr1 was not indexed")
  }

  println("\n2. Inspect a clipped and gapped forward alignment")
  let forward = document.alignments[0]
  let stats = forward.stats()
  println("  " + forward.summary())
  println(
    "  CIGAR=" +
    @src.align_sam_cigar_string(forward.cigar) +
    ", query bases=" +
    stats.query_consumed.to_string() +
    ", reference bases=" +
    stats.reference_consumed.to_string(),
  )
  match forward.target_to_query(8) {
    Some(position) =>
      println("  target base 8 maps to query base " + position.to_string())
    None => println("  target base 8 lies in a gap")
  }
  let md = match forward.string_tag("MD") {
    Some(value) => value
    None => "unknown"
  }
  println("  NM=" + optional_int(forward.computed_nm()) + ", MD=" + md)

  println("\n3. Reconstruct aligned rows and calculate MD")
  let rows = forward.aligned_rows(reference_sequence=demo_reference()) catch {
    AlignSamError(message) => abort("row reconstruction failed: " + message)
  }
  let calculated_md = @src.align_sam_calculate_md(forward, demo_reference()) catch {
    AlignSamError(message) => abort("MD calculation failed: " + message)
  }
  println("  target  " + rows.0)
  println("  query   " + rows.1)
  println("  calculated MD=" + calculated_md)

  println("\n4. Inspect reverse-strand and spliced coordinates")
  let reverse = document.alignments[1]
  println(
    "  reverse=" +
    reverse.is_reverse().to_string() +
    ", biological query=" +
    reverse.query_sequence,
  )
  println(
    "  query path " +
    coordinate_text(reverse.query_coordinates) +
    ", skipped reference bases=" +
    reverse.stats().skipped.to_string(),
  )
  match reverse.target_to_query(27) {
    Some(position) =>
      println(
        "  target base 27 maps to reverse query base " + position.to_string(),
      )
    None => println("  target base 27 lies in the skipped interval")
  }

  println("\n5. Construct and round-trip a typed SAM record")
  let tags = [
    @src.align_sam_integer_tag("NM", 0),
    @src.align_sam_string_tag("RG", "demo"),
    @src.align_sam_float_array_tag("BF", [0.25, 0.75]),
  ] catch {
    AlignSamError(message) => abort("tag construction failed: " + message)
  }
  let created = @src.align_sam_create(
    "created_read",
    "chr1",
    40,
    15,
    "AACGTT",
    "1S4M1S",
    reverse=true,
    mapq=50,
    qualities=[30, 31, 32, 33, 34, 35],
    tags~,
  ) catch {
    AlignSamError(message) => abort("alignment construction failed: " + message)
  }
  let record = @src.align_sam_format(created) catch {
    AlignSamError(message) => abort("SAM formatting failed: " + message)
  }
  let reparsed = parse_sam("@SQ\tSN:chr1\tLN:40\n" + record).alignments[0]
  println("  " + record)
  println(
    "  round-trip query=" +
    reparsed.query_sequence +
    ", mapq=" +
    optional_int(reparsed.mapq) +
    ", tags=" +
    reparsed.tags.length().to_string(),
  )

  println("=== Demo Complete ===")
}
