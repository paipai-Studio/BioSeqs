///|
/// Bio.KEGG.REST Demo.
///
/// Demonstrates the REST-style KEGG API client (info / list / find / get /
/// conv / link) with an injected transport, so the URL construction and the
/// three-queries-per-second rate limiting can be shown without any network
/// access.

///|
fn main raise {
  let log : Array[String] = []
  let client = @src.KeggRestClient::new(
    fetch=fn(url) {
      log.push(url)
      "(simulated response for \{url})"
    },
    now=fn() { 0.0 },
  )

  println("=== URL construction ===")
  println(@src.kegg_rest_url("list", "pathway", arg2="hsa"))
  println(@src.kegg_rest_url("get", "hsa:10458", arg2="aaseq"))

  println("\n=== 1. info: database statistics ===")
  println(@src.kegg_rest_info(client, "pathway"))

  println("\n=== 2. list: entry lists ===")
  println(@src.kegg_rest_list(client, "pathway", org="hsa"))
  println(@src.kegg_rest_list_entries(client, ["hsa:10458", "hsa:10459"]))

  println("\n=== 3. find: keyword and chemical searches ===")
  println(@src.kegg_rest_find(client, "genes", "shiga toxin"))
  println(@src.kegg_rest_find(client, "compound", "C00001", option="exact_mass"))

  println("\n=== 4. get: entry retrieval ===")
  println(@src.kegg_rest_get(client, "hsa:10458", option="aaseq"))
  println(@src.kegg_rest_get_multi(client, ["hsa:10458", "hsa:10459"], option="ntseq"))

  println("\n=== 5. conv: cross-database identifier conversion ===")
  println(@src.kegg_rest_conv(client, "uniprot", "hsa"))
  println(@src.kegg_rest_conv_multi(client, "compound", ["pubchem"]))

  println("\n=== 6. link: related entries via cross-references ===")
  println(@src.kegg_rest_link(client, "pathway", "hsa:10458"))
  println(@src.kegg_rest_link_multi(client, "genes", ["dru:D00109"]))

  println("\n=== Issued requests (rate limited to 3/second) ===")
  for url in log {
    println(url)
  }

  println("\nDone.")
}
