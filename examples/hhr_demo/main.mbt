///|
fn main {
  println("=== Biopython Bio.Align.hhr Demo ===")

  let record = @src.hhr_parse(@src.hhr_sample_text()) catch {
    _ => abort("failed to parse HHR sample")
  }

  println("\n1. HH-suite run metadata")
  println("  " + record.summary())
  println("  Effective sequences: " + record.metadata.neff.to_string())
  println("  Searched HMMs: " + record.metadata.searched_hmms.to_string())

  println("\n2. Ranked profile hits")
  for alignment in record.alignments {
    println("  " + alignment.summary())
    println(
      "    query/target coverage: " +
      alignment.query_coverage().to_string() +
      " / " +
      alignment.target_coverage().to_string(),
    )
  }

  println("\n3. Filtering and target lookup")
  let confident = record.filter_by_probability(90.0)
  println("  Hits with probability >= 90%: " + confident.length().to_string())
  match record.find_target("target_B") {
    Some(alignment) =>
      println(
        "  target_B aligned query: " +
        alignment.query_sequence +
        "\n  target_B aligned target: " +
        alignment.target_sequence,
      )
    None => abort("missing target_B")
  }

  println("\n4. Coordinate mapping")
  match record.best_alignment() {
    Some(alignment) => {
      println("  Best hit: " + alignment.target_name)
      match alignment.query_to_target(3) {
        Some(target_position) =>
          println(
            "  Query position 3 maps to target position " +
            target_position.to_string() +
            " (zero-based)",
          )
        None => println("  Query position 3 maps to a target gap")
      }
      println("  Preserved confidence columns: " + alignment.confidence)
    }
    None => abort("missing best alignment")
  }

  println("\n5. Canonical serialization")
  let serialized = record.to_string()
  let reparsed = @src.hhr_parse(serialized) catch {
    _ => abort("failed to reparse canonical HHR")
  }
  println("  Serialized bytes: " + serialized.length().to_string())
  println("  Round trip preserved record: " + (record == reparsed).to_string())

  println("\n=== Demo Complete ===")
}
