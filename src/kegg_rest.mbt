///|
/// Bio.KEGG.REST - access to the REST-style KEGG online API.
///
/// Port of Biopython's Bio.KEGG.REST module, which provides simple access to
/// a range of KEGG databases via simple URLs (constructed for you here),
/// with any errors indicated via HTTP error levels:
/// https://www.kegg.jp/kegg/rest/keggapi.html
///
/// Provided operations:
///   - kegg_info: current statistics of a given database
///   - kegg_list: entry list for a database (or specified database entries)
///   - kegg_find: data search by keyword (or formula/exact_mass/mol_weight)
///   - kegg_get:  data retrieval (plain / aaseq / ntseq / mol / kcf /
///                image / kgml / json)
///   - kegg_conv: convert KEGG identifiers to/from outside identifiers
///   - kegg_link: find related entries via database cross-references
///
/// To avoid overloading the service, Biopython only allows three calls per
/// second; the client tracks the same 1/3 second delay via an injected
/// clock (no real sleeping - the caller's transport is responsible).

///|
/// Error raised by KEGG REST operations for invalid arguments or transport
/// failures.
pub suberror KeggRestError {
  KeggRestError(String)
}

///|
fn[T] kegg_rest_fail(msg : String) -> T raise KeggRestError {
  raise KeggRestError(msg)
}

///|
/// REST-style KEGG API client with an injectable transport.
pub struct KeggRestClient {
  priv fetch : (String) -> String raise KeggRestError
  priv now : () -> Double
  priv mut previous : Double
}

///|
/// Create a KEGG REST client.
///
/// - fetch: callback performing the HTTP transfer: url -> text.
/// - now: optional wall clock in seconds; defaults to a virtual clock.
pub fn KeggRestClient::new(
  fetch~ : (String) -> String raise KeggRestError,
  now? : () -> Double,
) -> KeggRestClient {
  { fetch, now: now.unwrap_or(fn() { 0.0 }), previous: 0.0, }
}

///|
/// Base URL of the KEGG REST API.
pub let kegg_rest_base_url : String = "https://rest.kegg.jp/"

///|
/// Enforce the "up to three queries per second" rule and issue the request.
fn KeggRestClient::call(
  self : KeggRestClient,
  path : String,
) -> String raise KeggRestError {
  // Biopython uses a delay of one third of a second between calls.
  let delay = 0.333333333
  let current = (self.now)()
  let wait = self.previous + delay - current
  if wait > 0 {
    self.previous = current + wait
  } else {
    self.previous = current
  }
  (self.fetch)(kegg_rest_base_url + path)
}

///|
/// Construct the canonical URL for an operation, e.g.
/// kegg_rest_url("list", "pathway", "hsa") ->
/// "https://rest.kegg.jp/list/pathway/hsa"
pub fn kegg_rest_url(
  op : String,
  arg1 : String,
  arg2? : String,
  arg3? : String,
) -> String {
  let url = kegg_rest_base_url + op + "/" + arg1
  let url = match arg2 {
    Some(a2) => url + "/" + a2
    None => url
  }
  let url = match arg3 {
    Some(a3) => url + "/" + a3
    None => url
  }
  url
}

///|
/// KEGG info - Displays the current statistics of a given database.
///
/// The argument db can be a KEGG database name (e.g. 'pathway' or its
/// official abbreviation, 'path'), or a KEGG organism code or T number
/// (e.g. 'hsa' or 'T01001' for human).
pub fn kegg_rest_info(
  client : KeggRestClient,
  database : String,
) -> String raise KeggRestError {
  client.call("info/" + database)
}

///|
/// KEGG list - Entry list for a database, or specified database entries.
///
/// For the pathway and module databases the optional organism can be used
/// to restrict the results. Passing another database together with an
/// organism raises an error (mirrors the original module).
pub fn kegg_rest_list(
  client : KeggRestClient,
  database : String,
  org? : String,
) -> String raise KeggRestError {
  if (database == "pathway" || database == "module") && org is Some(o) {
    return client.call("list/" + database + "/" + o)
  }
  if database.length() > 0 && org is Some(_) {
    kegg_rest_fail("Invalid database arg for kegg list request.")
  }
  client.call("list/" + database)
}

///|
/// KEGG list over a list of database entries (joined with "+").
/// The maximum number of databases is 100 for a kegg list query.
pub fn kegg_rest_list_entries(
  client : KeggRestClient,
  databases : Array[String],
) -> String raise KeggRestError {
  if databases.length() > 100 {
    kegg_rest_fail("Maximum number of databases is 100 for kegg list query")
  }
  client.call("list/" + join_plus(databases))
}

///|
/// KEGG find - Data search.
///
/// Finds entries with matching query keywords or other query data in a
/// given database. For the compound and drug database, set option to the
/// string 'formula', 'exact_mass' or 'mol_weight' to search on that field
/// only. Any other option combination raises an error.
pub fn kegg_rest_find(
  client : KeggRestClient,
  database : String,
  query : String,
  option? : String,
) -> String raise KeggRestError {
  let valid_chem = database == "compound" || database == "drug"
  match option {
    Some(opt) =>
      if valid_chem &&
        (opt == "formula" || opt == "exact_mass" || opt == "mol_weight") {
        return client.call("find/" + database + "/" + query + "/" + opt)
      } else {
        kegg_rest_fail("Invalid option arg for kegg find request.")
      }
    None => client.call("find/" + database + "/" + query)
  }
}

///|
/// KEGG find over a list of query terms (joined with "+").
pub fn kegg_rest_find_multi(
  client : KeggRestClient,
  database : String,
  queries : Array[String],
) -> String raise KeggRestError {
  kegg_rest_find(client, database, join_plus(queries))
}

///|
/// KEGG get - Data retrieval.
///
/// dbentries: KEGG database entries (single string). The option is one of
/// "aaseq", "ntseq", "mol", "kcf", "image", "kgml" or "json". The input is
/// limited up to 10 entries when using the array variant; the input is
/// limited to one pathway entry with the image or kgml option, and to one
/// compound/glycan/drug entry with the image option (not enforced beyond
/// the entry-count limit, like the original module).
pub fn kegg_rest_get(
  client : KeggRestClient,
  dbentries : String,
  option? : String,
) -> String raise KeggRestError {
  match option {
    Some(opt) =>
      if opt == "aaseq" ||
        opt == "ntseq" ||
        opt == "mol" ||
        opt == "kcf" ||
        opt == "image" ||
        opt == "kgml" ||
        opt == "json" {
        client.call("get/" + dbentries + "/" + opt)
      } else {
        kegg_rest_fail("Invalid option arg for kegg get request.")
      }
    None => client.call("get/" + dbentries)
  }
}

///|
/// KEGG get over a list of dbentries (joined with "+").
/// The maximum number of dbentries is 10 for a kegg get query.
pub fn kegg_rest_get_multi(
  client : KeggRestClient,
  dbentries : Array[String],
  option? : String,
) -> String raise KeggRestError {
  if dbentries.length() > 10 {
    kegg_rest_fail("Maximum number of dbentries is 10 for kegg get query")
  }
  let joined = join_plus(dbentries)
  match option {
    Some(opt) => kegg_rest_get(client, joined, option=opt)
    None => kegg_rest_get(client, joined)
  }
}

///|
/// KEGG conv - convert KEGG identifiers to/from outside identifiers.
///
/// For gene identifiers the pair (target, source) must involve one of
/// ncbi-gi / ncbi-geneid / uniprot and a KEGG organism; for chemical
/// substance identifiers one of drug / compound / glycan and pubchem /
/// chebi. option can be "turtle" or "n-triple".
pub fn kegg_rest_conv(
  client : KeggRestClient,
  target_db : String,
  source_db : String,
  option? : String,
) -> String raise KeggRestError {
  match option {
    Some(opt) =>
      if opt != "turtle" && opt != "n-triple" {
        kegg_rest_fail("Invalid option arg for kegg conv request.")
      }
    None => ()
  }
  if !kegg_rest_conv_valid(target_db, source_db) {
    kegg_rest_fail("Bad argument target_db or source_db for kegg conv request.")
  }
  match option {
    Some(opt) => client.call("conv/" + target_db + "/" + source_db + "/" + opt)
    None => client.call("conv/" + target_db + "/" + source_db)
  }
}

///|
/// KEGG conv over a list of source databases / dbentries (joined with "+").
pub fn kegg_rest_conv_multi(
  client : KeggRestClient,
  target_db : String,
  source_dbs : Array[String],
  option? : String,
) -> String raise KeggRestError {
  let joined = join_plus(source_dbs)
  match option {
    Some(opt) => kegg_rest_conv(client, target_db, joined, option=opt)
    None => kegg_rest_conv(client, target_db, joined)
  }
}

///|
/// Validate the (target_db, source_db) pair for kegg_conv.
fn kegg_rest_conv_valid(target_db : String, source_db : String) -> Bool {
  let gene_dbs = ["ncbi-gi", "ncbi-geneid", "uniprot"]
  let chem_kegg = ["drug", "compound", "glycan"]
  let chem_outside = ["pubchem", "chebi"]
  let is_gene = gene_dbs.contains(target_db) || gene_dbs.contains(source_db)
  if is_gene {
    return true
  }
  (chem_kegg.contains(target_db) && chem_outside.contains(source_db)) ||
  (chem_outside.contains(target_db) && chem_kegg.contains(source_db))
}

///|
/// KEGG link - find related entries by using database cross-references.
///
/// option can be "turtle" or "n-triple". Note: an invalid option raises
/// the "kegg conv request" message - faithful to the original module,
/// which reuses that error text in kegg_link.
pub fn kegg_rest_link(
  client : KeggRestClient,
  target_db : String,
  source_db : String,
  option? : String,
) -> String raise KeggRestError {
  match option {
    Some(opt) =>
      if opt != "turtle" && opt != "n-triple" {
        kegg_rest_fail("Invalid option arg for kegg conv request.")
      }
    None => ()
  }
  match option {
    Some(opt) => client.call("link/" + target_db + "/" + source_db + "/" + opt)
    None => client.call("link/" + target_db + "/" + source_db)
  }
}

///|
/// KEGG link over a list of source databases / dbentries (joined with "+").
pub fn kegg_rest_link_multi(
  client : KeggRestClient,
  target_db : String,
  source_dbs : Array[String],
  option? : String,
) -> String raise KeggRestError {
  let joined = join_plus(source_dbs)
  match option {
    Some(opt) => kegg_rest_link(client, target_db, joined, option=opt)
    None => kegg_rest_link(client, target_db, joined)
  }
}

///|
/// Join identifier strings with "+" (the KEGG REST multi-entry separator).
fn join_plus(items : Array[String]) -> String {
  let mut out = ""
  let mut first = true
  for item in items {
    if !first {
      out = out + "+"
    }
    first = false
    out = out + item
  }
  out
}
