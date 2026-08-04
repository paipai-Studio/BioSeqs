///|
fn main {
  println("=== Bio.UniGene Demo ===")

  let record = @src.unigene_read(@src.unigene_sample_text()) catch {
    _ => abort("failed to parse UniGene sample")
  }

  println("\n1. Cluster")
  println("  " + record.summary())
  println("  Species prefix: " + record.species)
  println("  Cytoband: " + record.cytoband)
  println("  Expression: " + record.expression.join(", "))
  println("  Restricted expression: " + record.restricted_expression.join(", "))

  println("\n2. Protein similarities")
  for similarity in record.protein_similarities {
    println(
      "  " +
      similarity.protein_id +
      " (taxon " +
      similarity.organism +
      "): " +
      similarity.percent +
      "% over " +
      similarity.alignment_length +
      " aa",
    )
  }

  println("\n3. Sequence entries")
  println("  mRNA: " + record.sequences_of_type("mRNA").length().to_string())
  println("  EST: " + record.sequences_of_type("EST").length().to_string())
  for sequence in record.image_sequences() {
    println(
      "  IMAGE clone " +
      sequence.image_id +
      ": " +
      sequence.accession +
      " (" +
      sequence.read_end +
      " end)",
    )
  }

  println("\n4. Mapping information")
  for transcript_map in record.transcript_maps {
    println(
      "  Marker " +
      transcript_map.marker +
      " on panel " +
      transcript_map.radiation_hybrid_panel,
    )
  }
  for site in record.sts {
    println("  STS " + site.accession + ", UniSTS " + site.unists)
  }

  println("\n5. Canonical serialization")
  let serialized = record.to_string()
  let reparsed = @src.unigene_read(serialized) catch {
    _ => abort("failed to reparse serialized UniGene record")
  }
  println("  Serialized bytes: " + serialized.length().to_string())
  println("  Round trip preserved record: " + (record == reparsed).to_string())

  println("\n=== Demo Complete ===")
}
