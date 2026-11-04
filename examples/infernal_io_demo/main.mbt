///|
fn main {
  println("=== Biopython Bio.SearchIO.InfernalIO Demo ===")

  println("\n1. Auto-detect and parse Infernal tabular output")
  let tabular = @src.infernal_tabular_sample()
  let format = @src.infernal_tabular_format(tabular) catch {
    _ => abort("failed to detect Infernal tabular format")
  }
  let format_name = match format {
    Format1 => "1"
    Format2 => "2"
    Format3 => "3"
  }
  println("   detected format: " + format_name)
  let queries = @src.infernal_parse_tabular(tabular) catch {
    _ => abort("failed to parse Infernal tabular sample")
  }
  for query in queries {
    println("   " + query.summary())
    for hit in query.hits {
      println(
        "     hit=" +
        hit.id +
        ", description=" +
        hit.description +
        ", HSPs=" +
        hit.hsps.length().to_string(),
      )
    }
  }

  println("\n2. Reverse-strand normalized coordinates")
  let reverse = queries[0].hits[0].hsps[1].fragments[0]
  println(
    "   Infernal 480..412 -> zero-based half-open [" +
    reverse.hit_start.to_string() +
    ", " +
    reverse.hit_end.to_string() +
    "), strand=" +
    (if reverse.hit_strand == @src.strand_minus() { "-" } else { "+" }),
  )

  println("\n3. Parse plain text and split local-end alignment")
  let text_queries = @src.infernal_parse_text(@src.infernal_text_sample()) catch {
    _ => abort("failed to parse Infernal plain-text sample")
  }
  let text_query = text_queries[0]
  let local_hsp = text_query.hits[0].hsps[0]
  println("   " + text_query.summary())
  println("   local-end fragments: " + local_hsp.fragments.length().to_string())
  for index in 0..<local_hsp.fragments.length() {
    let fragment = local_hsp.fragments[index]
    println(
      "     fragment " +
      (index + 1).to_string() +
      ": model [" +
      fragment.query_start.to_string() +
      ", " +
      fragment.query_end.to_string() +
      "), sequence [" +
      fragment.hit_start.to_string() +
      ", " +
      fragment.hit_end.to_string() +
      "), omitted model=" +
      fragment.model_omission_before.to_string(),
    )
  }

  println("\n4. Filter and select the best HSP")
  let filtered = queries[0].filter(maximum_evalue=1.0e-10, included_only=true) catch {
    _ => abort("failed to filter Infernal results")
  }
  println("   retained HSPs: " + filtered.count_hsps().to_string())
  match queries[0].hits[0].best_hsp() {
    Some(best) =>
      println(
        "   best HSP: score=" +
        best.bitscore.to_string() +
        ", E-value=" +
        best.evalue.to_string(),
      )
    None => abort("expected a best Infernal HSP")
  }

  println("\n5. Convert to the generic SearchIO hierarchy")
  let generic = text_query.to_searchio()
  println(
    "   QueryResult id=" +
    generic.id +
    ", hits=" +
    generic.hits.length().to_string() +
    ", fragments=" +
    generic.hits[0].hsps[0].fragments.length().to_string(),
  )

  println(
    "\nScope: Infernal 1.0+ non-verbose cmscan/cmsearch tabular formats 1/2/3, plain text, --noali, CM and HMM-only output.",
  )
  println("\n=== Demo Complete ===")
}
