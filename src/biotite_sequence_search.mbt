///|
/// biotite_sequence_search.mbt — BLAST command/URL construction and
/// tabular output parsing (no process execution).

///|
pub(all) struct BiotiteBlastQuery {
  program : String
  database : String
}

///|
pub fn biotite_blast_query(
  program : String,
  database : String,
) -> BiotiteBlastQuery {
  { program, database, }
}

///|
/// Build a local BLAST command line.
pub fn biotite_blast_local_command(
  query : BiotiteBlastQuery,
  query_file : String,
  out_file : String,
  outfmt~ : Int32 = 6,
) -> String {
  query.program +
  " -query " +
  query_file +
  " -db " +
  query.database +
  " -out " +
  out_file +
  " -outfmt " +
  outfmt.to_string()
}

///|
/// Build a remote NCBI BLAST URL for a query sequence.
pub fn biotite_blast_remote_url(
  query : BiotiteBlastQuery,
  sequence : String,
) -> String {
  "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=" +
  query.program +
  "&DATABASE=" +
  query.database +
  "&QUERY=" +
  sequence
}

///|
/// Parse BLAST tabular output (outfmt 6/7).
/// Returns (query_id, subject_id, percent_identity, e_value).
pub fn biotite_blast_parse_tabular(
  text : String,
) -> Array[(String, String, Double, Double)] {
  let lines = text.split("\n").to_array()
  let result : Array[(String, String, Double, Double)] = []
  for line in lines {
    if line.length() == 0 || line.has_prefix("#") {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let pident = match parse_double(fields[2].trim().to_owned()) {
      Some(v) => v
      None => 0.0
    }
    let evalue = match parse_double(fields[10].trim().to_owned()) {
      Some(v) => v
      None => 0.0
    }
    result.push((fields[0].to_owned(), fields[1].to_owned(), pident, evalue))
  }
  result
}

///|
/// Parse BLAST tabular output with subject coordinates.
/// Returns (query, subject, identity, alignment_length, evalue).
pub fn biotite_blast_parse_tabular_extended(
  text : String,
) -> Array[(String, String, Double, Int32, Double)] {
  let lines = text.split("\n").to_array()
  let result : Array[(String, String, Double, Int32, Double)] = []
  for line in lines {
    if line.length() == 0 || line.has_prefix("#") {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let pident = match parse_double(fields[2].trim().to_owned()) {
      Some(v) => v
      None => 0.0
    }
    let aln_len = parse_int(fields[3].trim().to_owned())
    let evalue = match parse_double(fields[10].trim().to_owned()) {
      Some(v) => v
      None => 0.0
    }
    result.push(
      (fields[0].to_owned(), fields[1].to_owned(), pident, aln_len, evalue),
    )
  }
  result
}
