///|
/// Bio.ExPASy.ScanProsite Demo.
///
/// Demonstrates the port of Biopython's Bio.ExPASy.ScanProsite module:
/// query URL construction (urllib.parse.urlencode quote_plus semantics),
/// parsing of a `<scanprosite_response>` XML document, and the injected-fetch
/// scan path. Everything runs offline against a canned response.

///|
fn main raise {
  println("=== 1. Query URL construction (quote_plus semantics) ===")
  println(@src.scanprosite_url(seq="MKWVTF"))
  println(@src.scanprosite_url(seq="A C+G"))
  println(@src.scanprosite_url(seq="MKWVTF", params=[("db", "sp"), ("sig", "P-x(2)-G")]))

  println("\n=== 2. Parsing a canned ScanProsite XML response ===")
  let record = @src.scanprosite_read(@src.scanprosite_sample_xml())
  println("n_match = \{record.n_match}, n_seq = \{record.n_seq}")
  for i = 0; i < record.matches.length(); i = i + 1 {
    let m = record.matches[i]
    println(
      "  match \{i}: \{m.signature_ac} on \{m.sequence_id} (\{m.sequence_ac}) \{m.start}-\{m.stop}",
    )
  }

  println("\n=== 3. scan() through an injected (offline) fetch ===")
  let calls : Array[String] = []
  let xml = @src.scanprosite_sample_xml()
  let raw = @src.scanprosite_scan(
    fetch=fn(url : String, _post : String?) -> String {
      calls.push(url)
      xml
    },
    seq="MKWVTF",
    params=[("db", "sp")],
  )
  println("requested URL : \{calls[0]}")
  println("response bytes: \{raw.length()}")

  println("\n=== 4. Non-xml output is rejected ===")
  let err = try {
    ignore(@src.scanprosite_url(seq="MKWVTF", output="txt"))
    ""
  } catch {
    ScanPrositeError(msg) => msg
  }
  println("error: \{err}")
}