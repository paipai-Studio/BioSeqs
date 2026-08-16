///|
fn bed_coordinate_path_text(
  coordinates : Array[@src.AlignBedCoordinate],
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
  println("=== Biopython Bio.Align.bed Demo ===")

  let document = @src.align_bed_parse(@src.align_bed_example_text()) catch {
    AlignBedError(message) => abort("BED parsing failed: " + message)
  }
  let plus = document.alignments[0]
  let minus = document.alignments[1]

  println("\n1. Parse mixed BED3 and BED12 records")
  println(
    "  records=" +
    document.alignments.length().to_string() +
    ", targets=" +
    document.targets().length().to_string() +
    ", queries=" +
    document.queries().length().to_string(),
  )
  println(
    "  numeric score=" +
    plus.score.unwrap().number().unwrap_or(-1.0).to_string(),
  )
  let text_score_document = @src.align_bed_parse(
    "chr1\t10\t20\tcurated-hit\thigh-confidence\t+",
  ) catch {
    AlignBedError(message) => abort("text BED score parsing failed: " + message)
  }
  println(
    "  text score=" +
    text_score_document.alignments[0].score.unwrap().text_value().unwrap_or(
      "missing",
    ),
  )

  println("\n2. Inspect plus/minus target-query coordinate paths")
  println("  plus:  " + bed_coordinate_path_text(plus.coordinates))
  println("  minus: " + bed_coordinate_path_text(minus.coordinates))
  println(
    "  plus strand=" +
    plus.strand() +
    ", minus strand=" +
    minus.strand(),
  )

  println("\n3. Reconstruct exon blocks and transcript sizes")
  for index = 0; index < plus.blocks().length(); index = index + 1 {
    let block = plus.blocks()[index]
    println(
      "  exon " +
      (index + 1).to_string() +
      ": target[" +
      block.target_start.to_string() +
      ", " +
      block.target_end.to_string() +
      ") query[" +
      block.query_start.to_string() +
      ", " +
      block.query_end.to_string() +
      ")",
    )
  }
  println(
    "  plus query size=" +
    plus.query_size().to_string() +
    ", minus query size=" +
    minus.query_size().to_string(),
  )

  println("\n4. Map aligned residues in both directions")
  println(
    "  plus target 1000 -> query " +
    plus.map_target_position(1000).unwrap_or(-1).to_string(),
  )
  println(
    "  plus query 1054 -> target " +
    plus.map_query_position(1054).unwrap_or(-1).to_string(),
  )
  println(
    "  minus target 2000 -> query " +
    minus.map_target_position(2000).unwrap_or(-1).to_string(),
  )
  println(
    "  intron target 3000 -> query " +
    plus.map_target_position(3000).unwrap_or(-1).to_string(),
  )

  println("\n5. Count path operations and query half-open intervals")
  let counts = plus.counts()
  let summary = document.summary()
  let hits = document.search("chr22", 4900, 5100) catch {
    AlignBedError(message) => abort("BED interval search failed: " + message)
  }
  println(
    "  aligned=" +
    counts.aligned.to_string() +
    ", intron bases=" +
    counts.target_skip_bases.to_string() +
    ", blocks=" +
    counts.blocks.to_string(),
  )
  println(
    "  document aligned=" +
    summary.aligned_bases.to_string() +
    ", plus/minus=" +
    summary.plus_strand_count.to_string() +
    "/" +
    summary.minus_strand_count.to_string(),
  )
  println("  chr22 overlaps [4900, 5100): " + hits.length().to_string())

  println("\n6. Project one alignment to BED3 through BED12")
  for columns = 3; columns <= 12; columns = columns + 1 {
    let line = @src.align_bed_format(plus, bed_columns=columns) catch {
      AlignBedError(message) => abort("BED formatting failed: " + message)
    }
    println(
      "  BED" +
      columns.to_string() +
      ": " +
      line[0:line.length() - 1].to_owned(),
    )
  }

  println("\n7. Write canonical BED12 and parse it again")
  let single = @src.AlignBedDocument::create([plus]) catch {
    AlignBedError(message) => abort("BED document creation failed: " + message)
  }
  let canonical = @src.align_bed_write(single) catch {
    AlignBedError(message) => abort("BED serialization failed: " + message)
  }
  let round_trip = @src.align_bed_parse(canonical) catch {
    AlignBedError(message) => abort("BED round-trip parsing failed: " + message)
  }
  println("  " + canonical[0:canonical.length() - 1].to_owned())
  println("  round-trip preserved: " + (round_trip == single).to_string())
}
