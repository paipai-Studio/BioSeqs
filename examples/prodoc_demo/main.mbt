///|
/// Bio.ExPASy.Prodoc Demo.
///
/// Demonstrates the port of Biopython's Bio.ExPASy.Prodoc module: parsing a
/// (possibly multi-record) prosite.doc document, reading a single record,
/// and the strict "exactly one record" behaviour of read().

///|
fn main raise {
  println("=== 1. Parse a two-record prosite.doc ===")
  let records = @src.prodoc_parse(@src.prodoc_sample_content())
  println("records: \{records.length()}")
  for i = 0; i < records.length(); i = i + 1 {
    let r = records[i]
    println("--- record \{i}: \{r.accession} ---")
    println("  text: \{r.text}")
    for j = 0; j < r.references.length(); j = j + 1 {
      let ref = r.references[j]
      println("  [\{ref.number}] \{ref.authors} \{ref.citation}")
    }
  }

  println("\n=== 2. Read a single record ===")
  let single =
    "{PDOC00001}\n{PS00001; ASN_GLYCOSYLATION}\n{BEGIN}\nSome text.\n{END}\n"
  match @src.prodoc_read(single) {
    Some(r) => println("accession = \{r.accession}, text = \{r.text}")
    None => println("(empty document)")
  }

  println("\n=== 3. read() rejects a multi-record document ===")
  let err = try {
    ignore(@src.prodoc_read(@src.prodoc_sample_content()))
    ""
  } catch {
    ProdocError(msg) => msg
  }
  println("error: \{err}")
}