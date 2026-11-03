///|
fn exonerate_text_demo_parse(text : String) -> @bio.ExonerateTextDocument {
  @bio.exonerate_text_parse(text) catch {
    ExonerateTextError(message) => {
      println("Exonerate text parse failed: " + message)
      abort(message)
    }
  }
}

///|
fn exonerate_text_demo_range(range : @bio.ExonerateTextRange) -> String {
  "[" + range.start.to_string() + ", " + range.end.to_string() + ")"
}

///|
fn main {
  println("=== Biopython Bio.SearchIO.ExonerateIO Offline Demo ===")
  let document = exonerate_text_demo_parse(@bio.exonerate_text_example())
  println("\n1. Report metadata and SearchIO hierarchy")
  println("  command=" + document.metadata.command_line)
  println("  host=" + document.metadata.hostname)
  println("  " + document.summary())

  let query = document.queries[0]
  let hit = query.hits[0]
  let hsp = hit.hsps[0]
  println(
    "  query=" +
    query.id +
    ", target=" +
    hit.id +
    ", model=" +
    hsp.model +
    ", score=" +
    hsp.score.to_string(),
  )

  println("\n2. Spliced fragments and normalized coordinates")
  let query_ranges = hsp.query_ranges()
  let hit_ranges = hsp.hit_ranges()
  for index = 0; index < hsp.fragments.length(); index = index + 1 {
    let fragment = hsp.fragments[index]
    println(
      "  fragment " +
      (index + 1).to_string() +
      ": query=" +
      exonerate_text_demo_range(query_ranges[index]) +
      ", target=" +
      exonerate_text_demo_range(hit_ranges[index]) +
      ", phase=" +
      fragment.phase.to_string(),
    )
    println(
      "    " +
      fragment.query_sequence +
      "\n    " +
      fragment.similarity +
      "\n    " +
      fragment.hit_sequence,
    )
  }

  println("\n3. Inter-fragment intron ranges")
  println(
    "  query=" +
    exonerate_text_demo_range(hsp.query_inter_ranges()[0]) +
    ", target=" +
    exonerate_text_demo_range(hsp.hit_inter_ranges()[0]),
  )

  let counts = hsp.counts()
  println("\n4. Alignment statistics and position projection")
  println(
    "  columns=" +
    counts.alignment_columns.to_string() +
    ", identities=" +
    counts.identities.to_string() +
    ", gaps=" +
    counts.gap_columns.to_string() +
    ", intervening gaps=" +
    counts.introns_or_ner_gaps.to_string(),
  )
  println(
    "  first aligned query base -> " +
    hsp.fragments[0].query_position(0).unwrap_or(-1).to_string(),
  )
  println(
    "\nThe demo parses an existing report artifact; it does not launch Exonerate.",
  )
}
