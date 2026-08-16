///|
/// Cellosaurus example: parse, query, and serialize a cell-line record.
///
/// Run: moon run examples/cellosaurus_demo/main.mbt
fn main {
  println("=== Cellosaurus Record Demo ===")
  try {
    let record = match @bio.cellosaurus_read(@bio.cellosaurus_sample_text()) {
      Some(record) => record
      None => abort("expected one Cellosaurus record")
    }

    println("Record: \{record.repr()}")
    println("Accession: \{record.accession}")
    println("Human origin: \{record.has_species("Homo sapiens")}")

    println("Synonyms:")
    for synonym in record.synonym_list() {
      println("  - \{synonym}")
    }

    println("Database cross-references:")
    for reference in record.cross_references {
      println("  - \{reference.database}: \{reference.accession}")
    }
    let ecacc = record.cross_references_for("ECACC")
    println("ECACC accessions: \{ecacc.length()}")

    let serialized = record.to_string()
    let reparsed = @bio.cellosaurus_parse(serialized)
    println("Round-trip records: \{reparsed.length()}")
  } catch {
    _ => println("Failed to parse the Cellosaurus sample")
  }
}
