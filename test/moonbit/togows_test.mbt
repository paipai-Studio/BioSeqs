///|
/// Tests for the Bio.TogoWS port. The HTTP transfer and the clock are
/// injected, so all behaviour is verified without network access. Expected
/// URL shapes and error messages mirror the original Biopython module.

///|
/// Build a client backed by a canned response table. Every request is
/// recorded (URL and POST body) in the given arrays.
fn tgw_make_client(
  responses : Map[String, String],
  calls : Array[String],
  posts : Array[String?],
) -> @src.TogoWSClient raise {
  @src.TogoWSClient::new(fetch=fn(url : String, post : String?) -> String {
    calls.push(url)
    posts.push(post)
    // every URL used in the tests has a canned response
    responses.get(url).unwrap_or("")
  })
}

///|
/// Standard metadata for the fake TogoWS service.
fn tgw_default_responses() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=16)
  m.set("http://togows.dbcls.jp/entry", "pubmed nuccore uniprot compound")
  m.set("http://togows.dbcls.jp/entry/pubmed?fields", "abstract authors title")
  m.set("http://togows.dbcls.jp/entry/pubmed?formats", "xml json medline")
  m.set("http://togows.dbcls.jp/convert/", "genbank.fasta fasta.genbank")
  m
}

///|
test "togows_quote_basic" {
  // unreserved characters are kept
  assert_eq(@src.togows_quote(text="abcXYZ0189._-~"), "abcXYZ0189._-~")
  // space is percent encoded
  assert_eq(@src.togows_quote(text="hello world"), "hello%20world")
  // default safe is "/" (like the original module's usage)
  assert_eq(@src.togows_quote(text="a/b/c"), "a/b/c")
  // commas are encoded (identifiers joined with commas)
  assert_eq(@src.togows_quote(text="1,2"), "1%2C2")
  // custom safe list
  assert_eq(@src.togows_quote(text="a,b", safe=","), "a,b")
  // other ASCII characters
  assert_eq(@src.togows_quote(text="?=&#"), "%3F%3D%26%23")
  // non-ASCII is UTF-8 percent encoded
  assert_eq(@src.togows_quote(text="é"), "%C3%A9")
  assert_eq(@src.togows_quote(text="β"), "%CE%B2")
  assert_eq(@src.togows_quote(text="日"), "%E6%97%A5")
}

///|
test "togows_url_builders" {
  assert_eq(
    @src.togows_entry_url(db="pubmed", id="123"),
    "http://togows.dbcls.jp/entry/pubmed/123",
  )
  assert_eq(
    @src.togows_entry_url(db="pubmed", id="123", format="xml", field="title"),
    "http://togows.dbcls.jp/entry/pubmed/123/title.xml",
  )
  assert_eq(
    @src.togows_entry_url(db="nuccore", id="NM_000001"),
    "http://togows.dbcls.jp/entry/nuccore/NM_000001",
  )
  // multiple ids are joined with commas and the comma is quoted
  assert_eq(
    @src.togows_entry_multi_url(db="pubmed", ids=["1", "2", "3"]),
    "http://togows.dbcls.jp/entry/pubmed/1%2C2%2C3",
  )
  assert_eq(
    @src.togows_search_count_url(db="pubmed", query="diabetes human"),
    "http://togows.dbcls.jp/search/pubmed/diabetes%20human/count",
  )
  assert_eq(
    @src.togows_search_url(
      db="pubmed",
      query="diabetes human",
      offset=1,
      limit=100,
    ),
    "http://togows.dbcls.jp/search/pubmed/diabetes%20human/1,100",
  )
  assert_eq(
    @src.togows_search_url(db="pubmed", query="x", format="json"),
    "http://togows.dbcls.jp/search/pubmed/x.json",
  )
  assert_eq(
    @src.togows_search_url(db="pubmed", query="x"),
    "http://togows.dbcls.jp/search/pubmed/x",
  )
  assert_eq(
    @src.togows_convert_url(in_format="genbank", out_format="fasta"),
    "http://togows.dbcls.jp/convert/genbank.fasta",
  )
}

///|
test "togows_search_url_validation" {
  // offset without limit (or vice versa) is rejected
  let r1 = try {
    let _ = @src.togows_search_url(db="pubmed", query="x", offset=1)
    false
  } catch {
    TogoError(_) => true
  }
  assert_true(r1)
  let r2 = try {
    let _ = @src.togows_search_url(db="pubmed", query="x", limit=1)
    false
  } catch {
    TogoError(_) => true
  }
  assert_true(r2)
  // non-positive values are rejected, offset first
  let r3 = try {
    let _ = @src.togows_search_url(db="pubmed", query="x", offset=0, limit=1)
    false
  } catch {
    TogoError(msg) => msg == "Offset should be at least one, not 0"
  }
  assert_true(r3)
  let r4 = try {
    let _ = @src.togows_search_url(db="pubmed", query="x", offset=1, limit=0)
    false
  } catch {
    TogoError(msg) => msg == "Count should be at least one, not 0"
  }
  assert_true(r4)
}

///|
test "togows_entry_client" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses = tgw_default_responses()
  responses.set(
    "http://togows.dbcls.jp/entry/pubmed/123.xml", "<xml>record</xml>",
  )
  responses.set("http://togows.dbcls.jp/entry/pubmed/123/title", "Some title")
  let client = tgw_make_client(responses, calls, posts)
  // plain fetch with format; the db list is requested first
  let body = client.entry(db="pubmed", id="123", format="xml")
  assert_eq(body, "<xml>record</xml>")
  assert_eq(calls[0], "http://togows.dbcls.jp/entry")
  assert_eq(calls[1], "http://togows.dbcls.jp/entry/pubmed?formats")
  assert_eq(calls[2], "http://togows.dbcls.jp/entry/pubmed/123.xml")
  // field extraction; the fields list is fetched and the URL gains /title
  let title = client.entry(db="pubmed", id="123", field="title")
  assert_eq(title, "Some title")
  assert_eq(calls[3], "http://togows.dbcls.jp/entry/pubmed?fields")
  assert_eq(calls[4], "http://togows.dbcls.jp/entry/pubmed/123/title")
  // the pubmed alias "ti" is mapped to "title"
  let fixed_title = client.entry(db="pubmed", id="123", field="ti")
  assert_eq(fixed_title, "Some title")
  assert_eq(calls[5], "http://togows.dbcls.jp/entry/pubmed/123/title")
  // multiple ids are joined with commas before quoting
  let _ = client.entry_multi(db="pubmed", ids=["1", "2"], format="xml")
  assert_eq(calls[6], "http://togows.dbcls.jp/entry/pubmed/1%2C2.xml")
  // the db list and metadata were fetched once and then cached
  let mut db_list_calls = 0
  for c in calls {
    if c == "http://togows.dbcls.jp/entry" {
      db_list_calls = db_list_calls + 1
    }
  }
  assert_eq(db_list_calls, 1)
}

///|
test "togows_entry_validation" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let client = tgw_make_client(tgw_default_responses(), calls, posts)
  // unsupported database
  let r1 = try {
    let _ = client.entry(db="nosuchdb", id="123")
    false
  } catch {
    TogoError(msg) =>
      msg ==
      "TogoWS entry fetch does not officially support database 'nosuchdb'."
  }
  assert_true(r1)
  // unsupported field: the sorted field list is part of the message
  let r2 = try {
    let _ = client.entry(db="pubmed", id="123", field="genome")
    false
  } catch {
    TogoError(msg) =>
      msg ==
      "TogoWS entry fetch does not explicitly support field 'genome' for database 'pubmed'. Only: abstract, authors, title"
  }
  assert_true(r2)
  // unsupported format
  let r3 = try {
    let _ = client.entry(db="pubmed", id="123", format="gff")
    false
  } catch {
    TogoError(msg) =>
      msg ==
      "TogoWS entry fetch does not explicitly support format 'gff' for database 'pubmed'. Only: json, medline, xml"
  }
  assert_true(r3)
}

///|
test "togows_search_count" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses : Map[String, String] = Map([], capacity=8)
  responses.set(
    "http://togows.dbcls.jp/search/pubmed/diabetes%20human/count", "1234\n",
  )
  let client = tgw_make_client(responses, calls, posts)
  assert_eq(client.search_count(db="pubmed", query="diabetes human"), 1234)
  assert_eq(
    calls[0],
    "http://togows.dbcls.jp/search/pubmed/diabetes%20human/count",
  )
  // empty response
  let responses2 : Map[String, String] = Map([], capacity=8)
  responses2.set("http://togows.dbcls.jp/search/pubmed/x/count", "")
  let client2 = tgw_make_client(responses2, [], [])
  let r1 = try {
    let _ = client2.search_count(db="pubmed", query="x")
    false
  } catch {
    TogoError(msg) =>
      msg ==
      "TogoWS returned no data from URL http://togows.dbcls.jp/search/pubmed/x/count"
  }
  assert_true(r1)
  // non-integer response
  let responses3 : Map[String, String] = Map([], capacity=8)
  responses3.set("http://togows.dbcls.jp/search/pubmed/x/count", "abc\n")
  let client3 = tgw_make_client(responses3, [], [])
  let r2 = try {
    let _ = client3.search_count(db="pubmed", query="x")
    false
  } catch {
    TogoError(msg) =>
      msg ==
      "Expected an integer from URL http://togows.dbcls.jp/search/pubmed/x/count, got: 'abc\n'"
  }
  assert_true(r2)
}

///|
test "togows_search_client" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses : Map[String, String] = Map([], capacity=8)
  responses.set("http://togows.dbcls.jp/search/pubmed/x/1,2", "111\n222\n")
  let client = tgw_make_client(responses, calls, posts)
  let body = client.search(db="pubmed", query="x", offset=1, limit=2)
  assert_eq(body, "111\n222\n")
  assert_eq(calls[0], "http://togows.dbcls.jp/search/pubmed/x/1,2")
}

///|
test "togows_search_iter_batches" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses : Map[String, String] = Map([], capacity=8)
  responses.set("http://togows.dbcls.jp/search/pubmed/q/count", "5")
  responses.set("http://togows.dbcls.jp/search/pubmed/q/1,2", "id1\nid2\n")
  responses.set("http://togows.dbcls.jp/search/pubmed/q/3,2", "id3\nid4\n")
  responses.set("http://togows.dbcls.jp/search/pubmed/q/5,1", "id5\n")
  let client = tgw_make_client(responses, calls, posts)
  let ids = client.search_iter(db="pubmed", query="q", batch=2)
  assert_eq(ids, ["id1", "id2", "id3", "id4", "id5"])
  assert_eq(calls.length(), 4)
  assert_eq(calls[0], "http://togows.dbcls.jp/search/pubmed/q/count")
  assert_eq(calls[1], "http://togows.dbcls.jp/search/pubmed/q/1,2")
  assert_eq(calls[2], "http://togows.dbcls.jp/search/pubmed/q/3,2")
  assert_eq(calls[3], "http://togows.dbcls.jp/search/pubmed/q/5,1")
}

///|
test "togows_search_iter_limit" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses : Map[String, String] = Map([], capacity=8)
  responses.set("http://togows.dbcls.jp/search/pubmed/q/count", "100")
  responses.set("http://togows.dbcls.jp/search/pubmed/q/1,2", "id1\nid2\n")
  responses.set("http://togows.dbcls.jp/search/pubmed/q/3,1", "id3\n")
  let client = tgw_make_client(responses, calls, posts)
  // the overall limit trims the last batch to what remains
  let ids = client.search_iter(db="pubmed", query="q", limit=3, batch=2)
  assert_eq(ids, ["id1", "id2", "id3"])
  assert_eq(calls[2], "http://togows.dbcls.jp/search/pubmed/q/3,1")
}

///|
test "togows_search_iter_stalled_results" {
  // a batch identical to the previous one raises
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses : Map[String, String] = Map([], capacity=8)
  responses.set("http://togows.dbcls.jp/search/pubmed/q/count", "4")
  responses.set("http://togows.dbcls.jp/search/pubmed/q/1,2", "a\nb\n")
  responses.set("http://togows.dbcls.jp/search/pubmed/q/3,2", "a\nb\n")
  let client = tgw_make_client(responses, calls, posts)
  let r1 = try {
    let _ = client.search_iter(db="pubmed", query="q", batch=2)
    false
  } catch {
    TogoError(msg) => msg == "Same search results for previous offset"
  }
  assert_true(r1)
  // an identifier repeated from the previous batch raises
  let responses2 : Map[String, String] = Map([], capacity=8)
  responses2.set("http://togows.dbcls.jp/search/pubmed/q/count", "4")
  responses2.set("http://togows.dbcls.jp/search/pubmed/q/1,2", "a\nb\n")
  responses2.set("http://togows.dbcls.jp/search/pubmed/q/3,2", "b\nc\n")
  let client2 = tgw_make_client(responses2, [], [])
  let r2 = try {
    let _ = client2.search_iter(db="pubmed", query="q", batch=2)
    false
  } catch {
    TogoError(msg) => msg == "Result b was in previous batch"
  }
  assert_true(r2)
}

///|
test "togows_convert" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses = tgw_default_responses()
  responses.set("http://togows.dbcls.jp/convert/genbank.fasta", ">sp|p|k\nMKV")
  let client = tgw_make_client(responses, calls, posts)
  let out = client.convert(
    data="LOCUS       X\n//",
    in_format="genbank",
    out_format="fasta",
  )
  assert_eq(out, ">sp|p|k\nMKV")
  assert_eq(calls[0], "http://togows.dbcls.jp/convert/")
  assert_eq(calls[1], "http://togows.dbcls.jp/convert/genbank.fasta")
  // the input data is sent as the POST body
  assert_eq(posts[1], Some("LOCUS       X\n//"))
  // unsupported conversion lists the known pairs
  let r1 = try {
    let _ = client.convert(data="x", in_format="fasta", out_format="genbank")
    true
  } catch {
    TogoError(msg) =>
      msg ==
      "Unsupported conversion. Choose from:\ngenbank -> fasta\nfasta -> genbank"
  }
  assert_true(r1)
}

///|
test "togows_rate_limit_virtual_clock" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses : Map[String, String] = Map([], capacity=8)
  responses.set("http://togows.dbcls.jp/search/pubmed/x/1,1", "1\n")
  let client = tgw_make_client(responses, calls, posts)
  // default clock never advances, so each call is scheduled one delay later
  let _ = client.search(db="pubmed", query="x", offset=1, limit=1)
  assert_true((@src.TOGOWS_DELAY - client.last_call_time()).abs() < 1.0e-12)
  let _ = client.search(db="pubmed", query="x", offset=1, limit=1)
  assert_true(
    (2.0 * @src.TOGOWS_DELAY - client.last_call_time()).abs() < 1.0e-12,
  )
  let _ = client.search(db="pubmed", query="x", offset=1, limit=1)
  assert_true(
    (3.0 * @src.TOGOWS_DELAY - client.last_call_time()).abs() < 1.0e-12,
  )
}

///|
test "togows_rate_limit_injected_clock" {
  let calls : Array[String] = []
  let posts : Array[String?] = []
  let responses : Map[String, String] = Map([], capacity=8)
  responses.set("http://togows.dbcls.jp/search/pubmed/x/1,1", "1\n")
  let clock : Ref[Double] = Ref(100.0)
  let client = @src.TogoWSClient::new(
    fetch=fn(url : String, post : String?) -> String {
      calls.push(url)
      posts.push(post)
      responses.get(url).unwrap_or("")
    },
    now=fn() -> Double { clock.val },
  )
  // first call: nothing scheduled yet, previous = now
  let _ = client.search(db="pubmed", query="x", offset=1, limit=1)
  assert_true((client.last_call_time() - 100.0).abs() < 1.0e-12)
  // 0.4s later: enough time has passed, previous = now
  clock.val = 100.4
  let _ = client.search(db="pubmed", query="x", offset=1, limit=1)
  assert_true((client.last_call_time() - 100.4).abs() < 1.0e-12)
  // 0.1s later: the call must wait 0.233333333s, previous = now + wait
  clock.val = 100.5
  let _ = client.search(db="pubmed", query="x", offset=1, limit=1)
  assert_true((client.last_call_time() - 100.733333333).abs() < 1.0e-9)
}
