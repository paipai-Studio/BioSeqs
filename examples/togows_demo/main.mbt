///|
/// Bio.TogoWS Demo.
///
/// Demonstrates the port of Biopython's Bio.TogoWS module: percent quoting,
/// URL construction, entry/search/search_count/search_iter/convert calls,
/// and the built-in rate limiter. The HTTP transfer is injected, so this
/// demo runs against a small canned response table instead of the network.

///|
fn main raise {
  println("=== 1. Percent quoting (urllib.parse.quote semantics) ===")
  println("togows_quote(\"hello world\") = \{@src.togows_quote(text="hello world")}")
  println("togows_quote(\"1,2,3\")        = \{@src.togows_quote(text="1,2,3")}")
  println("togows_quote(\"a/b\", safe=\"\") = \{@src.togows_quote(text="a/b", safe="")}")

  println("\n=== 2. URL builders ===")
  println(
    @src.togows_entry_url(db="pubmed", id="123456", format="xml"),
  )
  println(
    @src.togows_entry_url(db="pubmed", id="123456", field="title"),
  )
  println(
    @src.togows_entry_multi_url(db="uniprot", ids=["P53_HUMAN", "KRAS_HUMAN"]),
  )
  println(
    @src.togows_search_count_url(db="pubmed", query="diabetes human"),
  )
  println(
    @src.togows_search_url(db="pubmed", query="diabetes human", offset=1, limit=100),
  )
  println(
    @src.togows_convert_url(in_format="genbank", out_format="fasta"),
  )

  println("\n=== 3. Client with a canned (offline) fetch ===")
  let calls : Array[String] = []
  let responses : Map[String, String] = Map([], capacity=16)
  responses.set("http://togows.dbcls.jp/entry", "pubmed nuccore uniprot compound")
  responses.set(
    "http://togows.dbcls.jp/entry/pubmed?fields",
    "abstract authors title",
  )
  responses.set("http://togows.dbcls.jp/entry/pubmed?formats", "xml json medline")
  responses.set(
    "http://togows.dbcls.jp/entry/pubmed/123456.xml",
    "<article>CRISPR in 2020</article>",
  )
  responses.set(
    "http://togows.dbcls.jp/entry/pubmed/123456/title",
    "CRISPR in 2020",
  )
  let client = @src.TogoWSClient::new(
    fetch=fn(url : String, _post : Option[String]) -> String {
      calls.push(url)
      responses.get(url).unwrap_or("")
    },
  )
  let record = client.entry(db="pubmed", id="123456", format="xml")
  println("entry(...)  -> \{record}")
  let title = client.entry(db="pubmed", id="123456", field="title")
  println("field ...   -> \{title}")
  // like the original module, the pubmed alias "ti" maps to "title"
  let aliased = client.entry(db="pubmed", id="123456", field="ti")
  println("alias ti    -> \{aliased}")

  println("\n=== 4. search_count and search ===")
  responses.set(
    "http://togows.dbcls.jp/search/pubmed/crispr/count",
    "4321",
  )
  let count = client.search_count(db="pubmed", query="crispr")
  println("search_count(pubmed, crispr) = \{count}")
  responses.set(
    "http://togows.dbcls.jp/search/pubmed/crispr/1,2",
    "36421234\n32109552\n",
  )
  let page = client.search(db="pubmed", query="crispr", offset=1, limit=2)
  println("first page: \{page}")

  println("\n=== 5. search_iter (batched iteration) ===")
  responses.set("http://togows.dbcls.jp/search/pubmed/crispr/3,2", "31909519\n30195406\n")
  let ids = client.search_iter(db="pubmed", query="crispr", batch=2)
  println("all ids: \{ids}")

  println("\n=== 6. convert (POST body) ===")
  responses.set("http://togows.dbcls.jp/convert/", "genbank.fasta fasta.genbank")
  responses.set("http://togows.dbcls.jp/convert/genbank.fasta", ">converted\nATGC")
  let fasta = client.convert(
    data="LOCUS       demo\n//",
    in_format="genbank",
    out_format="fasta",
  )
  println("convert(genbank -> fasta):\n\{fasta}")

  println("\n=== 7. Rate limiter (three queries per second) ===")
  println("calls made so far: \{calls.length()}")
  println(
    "scheduled next-call time (virtual clock): \{client.last_call_time()} s",
  )
  println("TOGOWS_DELAY = \{@src.TOGOWS_DELAY} s between requests")

  println("\nDone.")
}
