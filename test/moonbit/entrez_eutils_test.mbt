///|
/// Tests for the Bio.Entrez E-utilities request layer (einfo, esearch,
/// efetch, epost, esummary, espell, elink, egquery, ecitmatch). All HTTP
/// traffic goes through an injected transport, so every behaviour is
/// verified offline against the URL/parameter semantics and XML fixtures
/// (20060628/20190110 DTDs) of Biopython's Bio.Entrez 1.88.

///|
/// A captured request: the URL (GET) and the body (POST, if any).
struct CapturedRequest {
  url : String
  body : String
}

///|
/// Record requests on a virtual clock (time does not advance), echoing a
/// canned response.
fn capture_client(
  response : String,
  log : Array[CapturedRequest],
) -> @src.EntrezClient {
  @src.EntrezClient::new(
    fetch=fn(url, body) {
      log.push({
        url,
        body: match body {
          Some(b) => b
          None => ""
        },
      })
      response
    },
    now=fn() { 0.0 },
  )
}

///|
test "entrez_quote_plus_encoding" {
  assert_eq(@src.entrez_quote_plus("a b"), "a+b")
  // '/' is not in the safe set (urlencode uses quote_plus with safe='')
  assert_eq(@src.entrez_quote_plus("abc123_./-~"), "abc123_.%2F-~")
  assert_eq(@src.entrez_quote_plus("a&b=c"), "a%26b%3Dc")
  assert_eq(@src.entrez_quote_plus("x|y"), "x%7Cy")
  assert_eq(@src.entrez_quote_plus("(gene)"), "%28gene%29")
  // non-ASCII is percent-encoded over its UTF-8 bytes
  assert_eq(@src.entrez_quote_plus("\u{00e9}") == "%C3%A9", true)
}

///|
test "epost_posts_ids_and_parses_session" {
  let log : Array[CapturedRequest] = []
  let client = capture_client(
    "<?xml version=\"1.0\"?>\n ePostResult\n\n  <QueryKey>1</QueryKey>\n\n  <WebEnv>0A1B2C3D4E5F6</WebEnv>\n\n </ePostResult>",
    log,
  )
  let result = @src.epost(client, "nucleotide", ["568815597", "568815598"])
  assert_eq(result.query_key(), "1")
  assert_eq(result.web_env(), "0A1B2C3D4E5F6")
  // epost always uses HTTP POST
  assert_eq(log.length(), 1)
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi",
  )
  assert_eq(
    log[0].body,
    "db=nucleotide&id=568815597%2C568815598&tool=biopython",
  )
}

///|
test "epost_missing_fields_raise" {
  let r = try {
    let log : Array[CapturedRequest] = []
    let _ = @src.epost(
      capture_client("<ePostResult><QueryKey>1</QueryKey></ePostResult>", log),
      "pubmed",
      ["1"],
    )
    false
  } catch {
    EntrezUtilError(msg) => msg == "epost response is missing WebEnv"
  }
  assert_true(r)
}

///|
test "standard_params_order_and_overrides" {
  let log : Array[CapturedRequest] = []
  let client = @src.EntrezClient::new(
    fetch=fn(url, body) {
      log.push({
        url,
        body: match body {
          Some(b) => b
          None => ""
        },
      })
      "<QueryKey>1</QueryKey><WebEnv>W</WebEnv>"
    },
    now=fn() { 0.0 },
    email="user@example.com",
    api_key="KEY123",
  )
  let _ = @src.epost(client, "pubmed", ["1"])
  // user parameters first, then tool, email, api_key
  assert_eq(
    log[0].body,
    "db=pubmed&id=1&tool=biopython&email=user%40example.com&api_key=KEY123",
  )
  // user-supplied values are not overwritten
  let log2 : Array[CapturedRequest] = []
  let client2 = @src.EntrezClient::new(
    fetch=fn(url, body) {
      log2.push({
        url,
        body: match body {
          Some(b) => b
          None => ""
        },
      })
      "<QueryKey>1</QueryKey><WebEnv>W</WebEnv>"
    },
    now=fn() { 0.0 },
    email="user@example.com",
  )
  let _ = @src.epost(client2, "pubmed", ["1"])
  assert_eq(
    log2[0].body,
    "db=pubmed&id=1&tool=biopython&email=user%40example.com",
  )
}

///|
test "get_vs_post_selection" {
  // epost always uses HTTP POST (like the original module)
  let log : Array[CapturedRequest] = []
  let client = capture_client("<QueryKey>1</QueryKey><WebEnv>W</WebEnv>", log)
  let _ = @src.epost(client, "pubmed", ["1"])
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi",
  )
  assert_eq(log[0].body, "db=pubmed&id=1&tool=biopython")
  // a short espell request goes over GET
  let log2 : Array[CapturedRequest] = []
  let xml = "<eSpellResult><Database>pubmed</Database><Query>stem</Query></eSpellResult>"
  let client2 = capture_client(xml, log2)
  let _ = @src.espell(client2, db="pubmed", term="stem")
  assert_eq(
    log2[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi?db=pubmed&term=stem&tool=biopython",
  )
  assert_eq(log2[0].body, "")
  // more than 200 ids force a POST
  let many : Array[String] = []
  for i in 0..<200 {
    many.push(i.to_string())
  }
  let log3 : Array[CapturedRequest] = []
  let client3 = capture_client(xml, log3)
  let _ = @src.esummary(client3, "pubmed", many)
  // 200 ids force a POST over the plain cgi url
  assert_eq(
    log3[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
  )
  assert_eq(log3[0].body.contains("id=0%2C1%2C2%2C3"), true)
}

///|
test "esummary_parses_docsums" {
  let log : Array[CapturedRequest] = []
  let xml = "<eSummaryResult><DocSum><Id>19303033</Id><Item Name=\"PubDate\" Type=\"Date\">2009</Item><Item Name=\"Title\" Type=\"String\">Biopython</Item></DocSum><DocSum><Id>19303034</Id><Item Name=\"Title\" Type=\"String\">Next</Item></DocSum></eSummaryResult>"
  let client = capture_client(xml, log)
  let sums = @src.esummary(client, "pubmed", ["19303033", "19303034"])
  // short requests go over GET
  assert_eq(log.length(), 1)
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=19303033%2C19303034&tool=biopython",
  )
  assert_eq(sums.length(), 2)
  assert_eq(sums[0].id(), "19303033")
  assert_eq(sums[0].get("PubDate"), Some("2009"))
  assert_eq(sums[0].get("Title"), Some("Biopython"))
  assert_eq(sums[0].get("Missing"), None)
  let items = sums[0].items()
  assert_eq(items.length(), 2)
  assert_eq(items[0].name(), "PubDate")
  assert_eq(items[0].type_(), "Date")
  assert_eq(items[1].value(), "Biopython")
  assert_eq(sums[1].id(), "19303034")
}

///|
test "espell_parses_suggestions" {
  let xml = "<eSpellResult><Database>pubmed</Database><Query>biopythpn</Query><CorrectedQuery>biopython</CorrectedQuery><SpelledQuery><Replaced>biopython</Replaced></SpelledQuery></eSpellResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let result = @src.espell(client, db="pubmed", term="biopythpn")
  assert_eq(result.query(), "biopythpn")
  assert_eq(result.corrected_query(), "biopython")
  assert_eq(result.replaced(), ["biopython"])
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi?db=pubmed&term=biopythpn&tool=biopython",
  )
}

///|
test "espell_missing_correction_is_empty" {
  let xml = "<eSpellResult><Database>pubmed</Database><Query>biopython</Query></eSpellResult>"
  let client = capture_client(xml, [])
  let result = @src.espell(client, term="biopython")
  assert_eq(result.corrected_query(), "")
  assert_eq(result.replaced().length(), 0)
}

///|
test "elink_repeated_id_params_and_parse" {
  let xml = "<eLinkResult><LinkSet><DbFrom>pubmed</DbFrom><LinkSetDb><DbTo>pmc</DbTo><LinkName>pubmed_pmc</LinkName><Link><Id>27854912</Id></Link><Link><Id>27854913</Id></Link></LinkSetDb></LinkSet></eLinkResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let sets = @src.elink(client, dbfrom="pubmed", db="pmc", ids=[
    "27854912", "27854913",
  ])
  assert_eq(log.length(), 1)
  // multiple ids are transmitted as repeated "id" parameters
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pmc&id=27854912&id=27854913&tool=biopython",
  )
  assert_eq(sets.length(), 1)
  assert_eq(sets[0].db_from(), "pubmed")
  let lsdb = sets[0].link_set_dbs()
  assert_eq(lsdb.length(), 1)
  assert_eq(lsdb[0].db_to(), "pmc")
  assert_eq(lsdb[0].link_name(), "pubmed_pmc")
  assert_eq(lsdb[0].ids(), ["27854912", "27854913"])
}

///|
test "elink_cmd_and_linkname_params" {
  let xml = "<eLinkResult><LinkSet><DbFrom>pubmed</DbFrom></LinkSet></eLinkResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let _ = @src.elink(
    client,
    dbfrom="pubmed",
    ids=["1"],
    cmd="prlinks",
    linkname="pubmed_pubmed",
  )
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=1&cmd=prlinks&linkname=pubmed_pubmed&tool=biopython",
  )
}

///|
test "ecitmatch_sends_bdata_with_pipes" {
  let response = "proc natl acad sci u s a|1990|87|8227|Lifes,P|10fips|9025659\ncell|1991|64|671|Smith,J|11fips|20048972\n"
  let log : Array[CapturedRequest] = []
  let client = capture_client(response, log)
  let citations = [
    @src.EntrezCitation::new(
      journal_title="proc natl acad sci u s a",
      year="1990",
      volume="87",
      first_page="8227",
      author_name="Lifes P",
      key="10fips",
    ),
    @src.EntrezCitation::new(
      journal_title="cell",
      year="1991",
      volume="64",
      first_page="671",
      author_name="Smith J",
      key="11fips",
    ),
  ]
  let matches = @src.ecitmatch(client, "pubmed", citations~)
  assert_eq(matches.length(), 2)
  assert_eq(matches[0].0, "10fips")
  assert_eq(matches[0].1, "9025659")
  assert_eq(matches[1].0, "11fips")
  // citations are joined by \r (sent urlencoded as %0D); pipe characters
  // are not URL-encoded (ecitmatch expects literal "|"); short requests
  // go over GET, so the bdata is in the url
  assert_eq(
    log[0].url.contains(
      "bdata=proc+natl+acad+sci+u+s+a|1990|87|8227|Lifes+P|10fips%0Dcell|1991|64|671|Smith+J|11fips",
    ),
    true,
  )
}

///|
test "parse_ecitmatch_txt_and_xml_forms" {
  // plain-text pipe form: key and pmid are the last two fields
  let txt = "j|1990|87|8227|L P|key1|9025659\nj2|1991|64|671|S J|key2|20048972\n"
  assert_eq(@src.parse_ecitmatch(txt), [
    ("key1", "9025659"),
    ("key2", "20048972"),
  ])
  // XML form with ByteId elements
  let xml = "<eCitMatchResult><ByteId key=\"key1\" id=\"9025659\"/><ByteId key=\"key2\" id=\"20048972\"/></eCitMatchResult>"
  assert_eq(@src.parse_ecitmatch(xml), [
    ("key1", "9025659"),
    ("key2", "20048972"),
  ])
}

///|
test "rate_limiting_allows_bursts_on_virtual_clock" {
  // on the virtual clock the computed pause never blocks the caller;
  // five rapid requests must all go through
  let log : Array[CapturedRequest] = []
  let client = capture_client("<QueryKey>1</QueryKey><WebEnv>W</WebEnv>", log)
  for _ in 0..<5 {
    let _ = @src.epost(client, "pubmed", ["1"])
  }
  assert_eq(log.length(), 5)
}

///|
test "esearch_parses_ids_window_and_translation" {
  // mirrors the real esearch 20060628 DTD (see Biopython fixture
  // Tests/Entrez/esearch1.xml)
  let xml = "<eSearchResult><Count>63</Count><RetMax>20</RetMax><RetStart>0</RetStart><IdList>\n<Id>41282813</Id>\n<Id>41148224</Id>\n<Id>41011574</Id>\n</IdList><TranslationSet/><QueryTranslation>\"biopython\"[All Fields]</QueryTranslation></eSearchResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let result = client.esearch("pubmed", "biopython")
  assert_eq(result.count(), "63")
  assert_eq(result.retmax(), "20")
  assert_eq(result.retstart(), "0")
  assert_eq(result.ids(), ["41282813", "41148224", "41011574"])
  assert_eq(result.translation(), Some("\"biopython\"[All Fields]"))
  // no history without usehistory="y"
  assert_eq(result.query_key(), None)
  assert_eq(result.web_env(), None)
  // short request over GET, db/term first and tool last
  assert_eq(log.length(), 1)
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=biopython&tool=biopython",
  )
  assert_eq(log[0].body, "")
}

///|
test "esearch_optional_params_and_history" {
  let xml = "<eSearchResult><Count>71</Count><RetMax>3</RetMax><RetStart>10</RetStart><QueryKey>1</QueryKey><WebEnv>MCID_6aa204c5140f7f7b5404aca2</WebEnv><IdList><Id>42638005</Id></IdList><TranslationSet/><QueryTranslation>\"biopython\"[All Fields]</QueryTranslation></eSearchResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let result = client.esearch(
    "pubmed",
    "biopython",
    retstart="10",
    retmax="3",
    idtype="acc",
    usehistory="y",
  )
  assert_eq(result.count(), "71")
  assert_eq(result.query_key(), Some("1"))
  assert_eq(result.web_env(), Some("MCID_6aa204c5140f7f7b5404aca2"))
  assert_eq(result.ids(), ["42638005"])
  // parameters are appended in Biopython's order
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=biopython&retstart=10&retmax=3&idtype=acc&usehistory=y&tool=biopython",
  )
}

///|
test "esearch_missing_count_raises" {
  let r = try {
    let _ = capture_client("<eSearchResult><IdList/></eSearchResult>", []).esearch(
      "pubmed", "biopython",
    )
    false
  } catch {
    EntrezUtilError(msg) => msg == "esearch response is missing Count"
  }
  assert_true(r)
}

///|
test "efetch_returns_raw_text_with_rettype_retmode" {
  let record = "LOCUS       NM_001\nDEFINITION  demo record\n//\n"
  let log : Array[CapturedRequest] = []
  let client = capture_client(record, log)
  let text = client.efetch(
    "nuccore",
    ["568815597", "568815598"],
    rettype="gb",
    retmode="text",
  )
  assert_eq(text, record)
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=568815597%2C568815598&rettype=gb&retmode=text&tool=biopython",
  )
}

///|
test "efetch_supports_history_query_key_webenv" {
  let log : Array[CapturedRequest] = []
  let client = capture_client("RAW", log)
  let _ = client.efetch("protein", [], query_key="1", web_env="MCID_1")
  // Biopython spells the environment parameter WebEnv
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=&query_key=1&WebEnv=MCID_1&tool=biopython",
  )
}

///|
test "efetch_many_ids_force_post" {
  let many : Array[String] = []
  for i in 0..<200 {
    many.push(i.to_string())
  }
  let log : Array[CapturedRequest] = []
  let client = capture_client("RAW", log)
  let _ = client.efetch("nuccore", many, retmode="text")
  // 200 ids trigger POST over the plain cgi url (NCBI recommendation)
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
  )
  assert_eq(log[0].body.contains("id=0%2C1%2C2%2C3"), true)
  assert_eq(log[0].body.contains("retmode=text"), true)
}

///|
test "einfo_lists_databases_without_db" {
  // mirrors the real einfo 20190110 DTD (Biopython fixture einfo1.xml)
  let xml = "<eInfoResult>\n<DbList>\n<DbName>pubmed</DbName>\n<DbName>protein</DbName>\n<DbName>nuccore</DbName>\n</DbList>\n</eInfoResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let result = client.einfo()
  assert_eq(result.db_list(), ["pubmed", "protein", "nuccore"])
  // the database-specific fields stay empty in the list form
  assert_eq(result.dbname(), "")
  assert_eq(result.field_names().length(), 0)
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?tool=biopython",
  )
}

///|
test "einfo_db_stats_fields_and_links" {
  // <MenuName> is the 20190110 DTD spelling used by the Biopython 1.88
  // fixtures; the inner <Link> blocks use a different <Menu> tag
  let xml = "<eInfoResult><DbInfo><DbName>pubmed</DbName><MenuName>PubMed</MenuName><Description>PubMed bibliographic record</Description><Count>41127422</Count><LastUpdate>2026/09/09 11:56</LastUpdate><FieldList><Field><Name>ALL</Name><FullName>All Fields</FullName></Field><Field><Name>UID</Name><FullName>UID</FullName></Field></FieldList><LinkList><Link><Name>pubmed_pubmed</Name><Menu>Related Links</Menu></Link><Link><Name>pubmed_assembly</Name><Menu>Assembly</Menu></Link></LinkList></DbInfo></eInfoResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let result = client.einfo(db="pubmed")
  assert_eq(result.db_list().length(), 0)
  assert_eq(result.dbname(), "pubmed")
  assert_eq(result.menu(), "PubMed")
  assert_eq(result.description(), "PubMed bibliographic record")
  assert_eq(result.count(), "41127422")
  assert_eq(result.last_update(), "2026/09/09 11:56")
  assert_eq(result.field_names(), ["ALL", "UID"])
  assert_eq(result.link_names(), ["pubmed_pubmed", "pubmed_assembly"])
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=pubmed&tool=biopython",
  )
}

///|
test "egquery_parses_rows_and_lookup" {
  let xml = "<eGQueryResult><Term>biopython</Term><Result><DbName>pubmed</DbName><MenuName>PubMed</MenuName><Count>69</Count><Status>Ok</Status></Result><Result><DbName>pmc</DbName><MenuName>PubMed Central</MenuName><Count>22</Count><Status>Ok</Status></Result><Result><DbName>books</DbName><MenuName>Books</MenuName><Status>Term or database not supported</Status></Result></eGQueryResult>"
  let log : Array[CapturedRequest] = []
  let client = capture_client(xml, log)
  let result = client.egquery("biopython")
  assert_eq(result.term(), "biopython")
  let rows = result.rows()
  assert_eq(rows.length(), 3)
  assert_eq(rows[0].dbname(), "pubmed")
  assert_eq(rows[0].count(), 69)
  assert_eq(rows[1].dbname(), "pmc")
  assert_eq(rows[1].count(), 22)
  // a row without a Count parses as zero
  assert_eq(rows[2].dbname(), "books")
  assert_eq(rows[2].count(), 0)
  assert_eq(result.get("pmc"), Some(22))
  assert_eq(result.get("missing"), None)
  assert_eq(
    log[0].url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi?term=biopython&tool=biopython",
  )
}
