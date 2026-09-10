///|
/// Bio.Entrez E-utilities Demo.
///
/// Demonstrates the complete NCBI E-utilities request layer (einfo,
/// esearch, efetch, epost, esummary, espell, elink, egquery, ecitmatch)
/// with an injected transport, so the URL construction, parameter encoding
/// and XML parsing can all be shown without any network access.

///|
/// A canned-response transport that remembers every request URL/body.
let responses : Array[String] = [
  // einfo: list databases
  "<eInfoResult><DbList>\n<DbName>pubmed</DbName>\n<DbName>protein</DbName>\n<DbName>nuccore</DbName>\n<DbName>genome</DbName>\n</DbList></eInfoResult>",
  // esearch: ids plus a server-side history
   "<eSearchResult><Count>2</Count><RetMax>2</RetMax><RetStart>0</RetStart><QueryKey>1</QueryKey><WebEnv>MCID_0123456789abcdef</WebEnv><IdList><Id>19303033</Id><Id>19303034</Id></IdList><TranslationSet/><QueryTranslation>\"biopython\"[All Fields]</QueryTranslation></eSearchResult>",
  // efetch: raw record text, returned unchanged
   "PMID- 19303033\nTI  - Biopython: freely available Python tools\n//\n",
  // epost
   "<ePostResult><QueryKey>2</QueryKey><WebEnv>MCID_0123456789abcdef</WebEnv></ePostResult>",
  // esummary
   "<eSummaryResult><DocSum><Id>19303033</Id><Item Name=\"PubDate\" Type=\"Date\">2009</Item><Item Name=\"Title\" Type=\"String\">Biopython: freely available tools for computational molecular biology</Item></DocSum></eSummaryResult>",
  // espell
   "<eSpellResult><Database>pubmed</Database><Query>biopythpn</Query><CorrectedQuery>biopython</CorrectedQuery><SpelledQuery><Replaced>biopython</Replaced></SpelledQuery></eSpellResult>",
  // elink
   "<eLinkResult><LinkSet><DbFrom>pubmed</DbFrom><LinkSetDb><DbTo>pmc</DbTo><LinkName>pubmed_pmc</LinkName><Link><Id>27854912</Id></Link><Link><Id>27854913</Id></Link></LinkSetDb></LinkSet></eLinkResult>",
  // egquery
   "<eGQueryResult><Term>biopython</Term><Result><DbName>pubmed</DbName><MenuName>PubMed</MenuName><Count>63</Count><Status>Ok</Status></Result><Result><DbName>pmc</DbName><MenuName>PubMed Central</MenuName><Count>22</Count><Status>Ok</Status></Result></eGQueryResult>",
  // ecitmatch
   "proc natl acad sci u s a|1991|88|3248|mann bj|citation_1|2014248\n",
]

///|
fn main raise {
  let log : Array[String] = []
  let mut served = 0
  let client = @src.EntrezClient::new(
    fetch=fn(url, body) {
      let recorded = match body {
        Some(b) => "POST \{url}  body=\{b}"
        None => "GET  \{url}"
      }
      log.push(recorded)
      let response = responses[served]
      served = served + 1
      response
    },
    now=fn() { 0.0 },
    email="Your.Name.Here@example.org",
  )

  println("=== 1. EInfo: list the available databases ===")
  let info = client.einfo()
  println("databases: \{info.db_list()}")

  println("\n=== 2. ESearch: search PubMed and keep a history ===")
  let search = client.esearch("pubmed", "biopython", retmax="2", usehistory="y")
  println("count:     \{search.count()}")
  println("ids:       \{search.ids()}")
  println("translated: \{search.translation()}")
  println("QueryKey:  \{search.query_key()}")
  println("WebEnv:    \{search.web_env()}")

  println("\n=== 3. EFetch: retrieve the records (raw text) ===")
  let record = client.efetch(
    "pubmed",
    search.ids(),
    rettype="medline",
    retmode="text",
  )
  println(record)

  println("=== 4. EPost: store identifiers in the user environment ===")
  let posted = @src.epost(client, "nucleotide", ["568815597", "568815598"])
  println("QueryKey: \{posted.query_key()}")
  println("WebEnv:   \{posted.web_env()}")

  println("\n=== 5. ESummary: document summaries for primary ids ===")
  let sums = @src.esummary(client, "pubmed", ["19303033"])
  for doc in sums {
    let pub_date = doc.get("PubDate")
    let title = doc.get("Title")
    let pd = match pub_date {
      Some(v) => v
      None => ""
    }
    let t = match title {
      Some(v) => v
      None => ""
    }
    println("id:      \{doc.id()}")
    println("PubDate: \{pd}")
    println("Title:   \{t}")
  }

  println("\n=== 6. ESpell: spelling suggestions ===")
  let spell = @src.espell(client, db="pubmed", term="biopythpn")
  println("query:     \{spell.query()}")
  println("corrected: \{spell.corrected_query()}")
  println("replaced:  \{spell.replaced()}")

  println("\n=== 7. ELink: linked records across databases ===")
  let sets = @src.elink(client, dbfrom="pubmed", db="pmc", ids=[
    "27854912", "27854913",
  ])
  for set in sets {
    println("DbFrom: \{set.db_from()}")
    for lsdb in set.link_set_dbs() {
      println("DbTo: \{lsdb.db_to()}  LinkName: \{lsdb.link_name()}")
      println("linked ids: \{lsdb.ids()}")
    }
  }

  println("\n=== 8. EGQuery: global cross-database counts ===")
  let global = client.egquery("biopython")
  for row in global.rows() {
    println("\{row.dbname()}: \{row.count()}")
  }
  println("pubmed count via get: \{global.get("pubmed")}")

  println("\n=== 9. ECitMatch: PMIDs from citation strings ===")
  let citation = @src.EntrezCitation::new(
    journal_title="proc natl acad sci u s a",
    year="1991",
    volume="88",
    first_page="3248",
    author_name="mann bj",
    key="citation_1",
  )
  let matches = @src.ecitmatch(client, "pubmed", citations=[citation])
  for m in matches {
    println("key \{m.0} -> PMID \{m.1}")
  }

  println("\n=== 10. Recorded requests ===")
  for request in log {
    println(request)
  }

  println("\nDone.")
}
