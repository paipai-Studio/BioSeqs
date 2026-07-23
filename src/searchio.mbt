///|
/// Bio.SearchIO - Unified interface for BLAST/HMMER search result parsing.
/// Inspired by Biopython's Bio.SearchIO module.
///
/// Provides a common framework for parsing and iterating over search
/// results from different bioinformatics search tools (BLAST, HMMER, etc.).
/// The core data structures are QueryResult (query-level results),
/// Hit (database hit), and HSP (High-scoring Segment Pair alignment).

///|
/// Parsed search result format types.
pub enum SearchIOFormat {
  BlastXml
  BlastTabular
  BlastText
  HmmdomText
  HmmscanDomtbl
  Unknown
}

///|
/// Convert format string to SearchIOFormat.
pub fn search_io_format_from_str(s : String) -> SearchIOFormat {
  match s.to_lower().trim() {
    "blast-xml" | "blastxml" | "xml" => SearchIOFormat::BlastXml
    "blast-tab" | "blasttabular" | "tabular" | "tab" => SearchIOFormat::BlastTabular
    "blast-text" | "blasttext" | "text" => SearchIOFormat::BlastText
    "hmmscan-domtbl" | "domtbl" => SearchIOFormat::HmmscanDomtbl
    "hmmdom-text" => SearchIOFormat::HmmdomText
    _ => SearchIOFormat::Unknown
  }
}

///|
/// Parse a 32-bit signed integer from a StringView (fast, no error checking).
fn searchio_parse_int_fast(sv : StringView) -> Int {
  let n = sv.length()
  if n == 0 {
    return 0
  }
  let mut result = 0
  let zero = '0'.to_int()
  let mut i = 0
  let mut is_negative = false
  if sv.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }
  while i < n {
    let c = sv.unsafe_get(i).to_int()
    if c >= zero && c <= zero + 9 {
      result = result * 10 + (c - zero)
    }
    i = i + 1
  }
  if is_negative { -result } else { result }
}

///|
/// Parse a 64-bit float from a StringView (fast, supports scientific notation).
fn searchio_parse_double_fast(sv : StringView) -> Double {
  let n = sv.length()
  if n == 0 {
    return 0.0
  }
  let mut result = 0.0
  let mut decimal = 0.0
  let mut decimal_places = 0
  let mut is_negative = false
  let mut has_exponent = false
  let mut exponent_val = 0
  let mut exp_negative = false
  let zero = '0'.to_int()
  let mut i = 0
  if sv.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }
  while i < n {
    let c = sv.unsafe_get(i).to_int()
    if c == 46 {
      i = i + 1
      while i < n {
        let dc = sv.unsafe_get(i).to_int()
        if dc >= zero && dc <= zero + 9 {
          decimal = decimal * 10.0 + (dc - zero).to_double()
          decimal_places = decimal_places + 1
        }
        i = i + 1
      }
      break
    }
    if c == 101 || c == 69 {
      has_exponent = true
      i = i + 1
      if i < n && sv.unsafe_get(i).to_int() == 45 {
        exp_negative = true
        i = i + 1
      }
      while i < n {
        let ec = sv.unsafe_get(i).to_int()
        if ec >= zero && ec <= zero + 9 {
          exponent_val = exponent_val * 10 + (ec - zero)
        }
        i = i + 1
      }
      break
    }
    if c >= zero && c <= zero + 9 {
      result = result * 10.0 + (c - zero).to_double()
    }
    i = i + 1
  }
  while decimal_places > 0 {
    decimal = decimal / 10.0
    decimal_places = decimal_places - 1
  }
  let mut final_result = result + decimal
  if has_exponent {
    if exp_negative {
      let mut j = 0
      while j < exponent_val {
        final_result = final_result / 10.0
        j = j + 1
      }
    } else {
      let mut j = 0
      while j < exponent_val {
        final_result = final_result * 10.0
        j = j + 1
      }
    }
  }
  if is_negative { -final_result } else { final_result }
}

///|
/// Parse a 32-bit signed integer from a String (fast, no error checking).
fn searchio_parse_int(s : String) -> Int {
  let n = s.length()
  if n == 0 {
    return 0
  }
  let mut result = 0
  let zero = '0'.to_int()
  let mut i = 0
  let mut is_negative = false
  if s.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c >= zero && c <= zero + 9 {
      result = result * 10 + (c - zero)
    }
    i = i + 1
  }
  if is_negative { -result } else { result }
}

///|
/// Parse a 64-bit float from a String (fast, supports scientific notation).
fn searchio_parse_double(s : String) -> Double {
  let n = s.length()
  if n == 0 {
    return 0.0
  }
  let mut result = 0.0
  let mut decimal = 0.0
  let mut decimal_places = 0
  let mut is_negative = false
  let mut has_exponent = false
  let mut exponent_val = 0
  let mut exp_negative = false
  let zero = '0'.to_int()
  let mut i = 0
  if s.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c == 46 {
      i = i + 1
      while i < n {
        let dc = s.unsafe_get(i).to_int()
        if dc >= zero && dc <= zero + 9 {
          decimal = decimal * 10.0 + (dc - zero).to_double()
          decimal_places = decimal_places + 1
        }
        i = i + 1
      }
      break
    }
    if c == 101 || c == 69 {
      has_exponent = true
      i = i + 1
      if i < n && s.unsafe_get(i).to_int() == 45 {
        exp_negative = true
        i = i + 1
      }
      while i < n {
        let ec = s.unsafe_get(i).to_int()
        if ec >= zero && ec <= zero + 9 {
          exponent_val = exponent_val * 10 + (ec - zero)
        }
        i = i + 1
      }
      break
    }
    if c >= zero && c <= zero + 9 {
      result = result * 10.0 + (c - zero).to_double()
    }
    i = i + 1
  }
  while decimal_places > 0 {
    decimal = decimal / 10.0
    decimal_places = decimal_places - 1
  }
  let mut final_result = result + decimal
  if has_exponent {
    if exp_negative {
      let mut j = 0
      while j < exponent_val {
        final_result = final_result / 10.0
        j = j + 1
      }
    } else {
      let mut j = 0
      while j < exponent_val {
        final_result = final_result * 10.0
        j = j + 1
      }
    }
  }
  if is_negative { -final_result } else { final_result }
}

///|
/// High-scoring Segment Pair (HSP) - represents an alignment between query and hit.
pub struct SearchIOHsp {
  /// Score of the HSP
  bitscore : Double
  /// E-value
  evalue : Double
  /// Percent identity
  identity : Double
  /// Positives (conservative substitutions)
  positives : Double
  /// Gap percentage
  gap : Double
  /// Alignment length
  alignment_length : Int
  /// Query start position (1-based)
  query_start : Int
  /// Query end position (1-based)
  query_end : Int
  /// Hit start position (1-based)
  hit_start : Int
  /// Hit end position (1-based)
  hit_end : Int
  /// Query frame
  query_frame : Int
  /// Hit frame
  hit_frame : Int
  /// Query sequence in alignment
  query_seq : String
  /// Hit sequence in alignment
  hit_seq : String
  /// Middle line (matches/mismatches)
  midline : String
}

///|
/// Construct SearchIOHsp.
pub fn SearchIOHsp::new(
  bitscore~ : Double,
  evalue~ : Double,
  identity~ : Double,
  positives~ : Double,
  gap~ : Double,
  alignment_length~ : Int,
  query_start~ : Int,
  query_end~ : Int,
  hit_start~ : Int,
  hit_end~ : Int,
) -> SearchIOHsp {
  SearchIOHsp::{
    bitscore, evalue, identity, positives, gap,
    alignment_length, query_start, query_end, hit_start, hit_end,
    query_frame: 0, hit_frame: 0,
    query_seq: "", hit_seq: "", midline: ""
  }
}

///|
/// Set query sequence in the alignment.
pub fn SearchIOHsp::with_query_seq(self : SearchIOHsp, seq : String) -> SearchIOHsp {
  SearchIOHsp::{
    bitscore: self.bitscore, evalue: self.evalue, identity: self.identity,
    positives: self.positives, gap: self.gap, alignment_length: self.alignment_length,
    query_start: self.query_start, query_end: self.query_end,
    hit_start: self.hit_start, hit_end: self.hit_end,
    query_frame: self.query_frame, hit_frame: self.hit_frame,
    query_seq: seq, hit_seq: self.hit_seq, midline: self.midline
  }
}

///|
/// Set hit sequence in the alignment.
pub fn SearchIOHsp::with_hit_seq(self : SearchIOHsp, seq : String) -> SearchIOHsp {
  SearchIOHsp::{
    bitscore: self.bitscore, evalue: self.evalue, identity: self.identity,
    positives: self.positives, gap: self.gap, alignment_length: self.alignment_length,
    query_start: self.query_start, query_end: self.query_end,
    hit_start: self.hit_start, hit_end: self.hit_end,
    query_frame: self.query_frame, hit_frame: self.hit_frame,
    query_seq: self.query_seq, hit_seq: seq, midline: self.midline
  }
}

///|
/// Set midline (conservation line).
pub fn SearchIOHsp::with_midline(self : SearchIOHsp, line : String) -> SearchIOHsp {
  SearchIOHsp::{
    bitscore: self.bitscore, evalue: self.evalue, identity: self.identity,
    positives: self.positives, gap: self.gap, alignment_length: self.alignment_length,
    query_start: self.query_start, query_end: self.query_end,
    hit_start: self.hit_start, hit_end: self.hit_end,
    query_frame: self.query_frame, hit_frame: self.hit_frame,
    query_seq: self.query_seq, hit_seq: self.hit_seq, midline: line
  }
}

///|
/// Set query and hit frames.
pub fn SearchIOHsp::with_frames(self : SearchIOHsp, qframe : Int, hframe : Int) -> SearchIOHsp {
  SearchIOHsp::{
    bitscore: self.bitscore, evalue: self.evalue, identity: self.identity,
    positives: self.positives, gap: self.gap, alignment_length: self.alignment_length,
    query_start: self.query_start, query_end: self.query_end,
    hit_start: self.hit_start, hit_end: self.hit_end,
    query_frame: qframe, hit_frame: hframe,
    query_seq: self.query_seq, hit_seq: self.hit_seq, midline: self.midline
  }
}

///|
/// Calculate number of identical residues from identity percentage and length.
pub fn SearchIOHsp::n_identical(self : SearchIOHsp) -> Int {
  (self.identity / 100.0 * self.alignment_length.to_double()).round().to_int()
}

///|
/// Calculate number of positive (conserved) residues.
pub fn SearchIOHsp::n_positives(self : SearchIOHsp) -> Int {
  (self.positives / 100.0 * self.alignment_length.to_double()).round().to_int()
}

///|
/// Calculate number of gaps in the alignment.
pub fn SearchIOHsp::n_gaps(self : SearchIOHsp) -> Int {
  (self.gap / 100.0 * self.alignment_length.to_double()).round().to_int()
}

///|
/// Hit - represents a database sequence that matched the query.
pub struct SearchIOHit {
  /// Hit identifier (accession or name)
  id : String
  /// Hit description
  description : String
  /// Sequence length of the hit
  seq_length : Int
  /// List of HSPs for this hit
  hsps : Array[SearchIOHsp]
}

///|
/// Construct SearchIOHit.
pub fn SearchIOHit::new(id~ : String, description~ : String, seq_length~ : Int) -> SearchIOHit {
  SearchIOHit::{ id, description, seq_length, hsps: Array::new() }
}

///|
/// Add an HSP to this hit.
pub fn SearchIOHit::add_hsp(self : SearchIOHit, hsp : SearchIOHsp) -> SearchIOHit {
  let new_hsps = self.hsps.copy()
  new_hsps.push(hsp)
  SearchIOHit::{
    id: self.id, description: self.description,
    seq_length: self.seq_length, hsps: new_hsps
  }
}

///|
/// Get the best (highest-scoring) HSP for this hit.
pub fn SearchIOHit::best_hsp(self : SearchIOHit) -> SearchIOHsp? {
  if self.hsps.length() == 0 {
    return None
  }
  let mut best = self.hsps[0]
  let mut i = 1
  while i < self.hsps.length() {
    if self.hsps[i].bitscore > best.bitscore {
      best = self.hsps[i]
    }
    i = i + 1
  }
  Some(best)
}

///|
/// Get total number of HSPs.
pub fn SearchIOHit::n_hsps(self : SearchIOHit) -> Int {
  self.hsps.length()
}

///|
/// Get the sum of all bitscores.
pub fn SearchIOHit::sum_bitscore(self : SearchIOHit) -> Double {
  let mut total = 0.0
  let mut i = 0
  while i < self.hsps.length() {
    total = total + self.hsps[i].bitscore
    i = i + 1
  }
  total
}

///|
/// QueryResult - all search results for a single query sequence.
pub struct SearchIOQueryResult {
  /// Query identifier
  id : String
  /// Query description
  description : String
  /// Query sequence length
  seq_length : Int
  /// Search database name
  database : String
  /// Number of hits
  n_hits : Int
  /// List of hits
  hits : Array[SearchIOHit]
  /// Total search statistics
  total_hits : Int
}

///|
/// Construct SearchIOQueryResult.
pub fn SearchIOQueryResult::new(
  id~ : String,
  description~ : String,
  seq_length~ : Int,
  database~ : String,
) -> SearchIOQueryResult {
  SearchIOQueryResult::{
    id, description, seq_length, database,
    n_hits: 0, hits: Array::new(), total_hits: 0
  }
}

///|
/// Add a hit to this query result.
pub fn SearchIOQueryResult::add_hit(self : SearchIOQueryResult, hit : SearchIOHit) -> SearchIOQueryResult {
  let new_hits = self.hits.copy()
  new_hits.push(hit)
  SearchIOQueryResult::{
    id: self.id, description: self.description,
    seq_length: self.seq_length, database: self.database,
    n_hits: self.n_hits + 1, hits: new_hits,
    total_hits: self.total_hits
  }
}

///|
/// Set total hits (before filtering).
pub fn SearchIOQueryResult::with_total_hits(self : SearchIOQueryResult, total : Int) -> SearchIOQueryResult {
  SearchIOQueryResult::{
    id: self.id, description: self.description,
    seq_length: self.seq_length, database: self.database,
    n_hits: self.n_hits, hits: self.hits,
    total_hits: total
  }
}

///|
/// Sort hits by bitscore (descending).
pub fn SearchIOQueryResult::sort_by_score(self : SearchIOQueryResult) -> SearchIOQueryResult {
  let sorted = self.hits.copy()
  let n = sorted.length()
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let score_i = sorted[i].sum_bitscore()
      let score_j = sorted[j].sum_bitscore()
      if score_j > score_i {
        let temp = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = temp
      }
      j = j + 1
    }
    i = i + 1
  }
  SearchIOQueryResult::{
    id: self.id, description: self.description,
    seq_length: self.seq_length, database: self.database,
    n_hits: self.n_hits, hits: sorted,
    total_hits: self.total_hits
  }
}

///|
/// Filter hits by E-value threshold.
pub fn SearchIOQueryResult::filter_evalue(self : SearchIOQueryResult, threshold : Double) -> SearchIOQueryResult {
  let filtered : Array[SearchIOHit] = Array::new()
  let mut i = 0
  while i < self.hits.length() {
    let hit = self.hits[i]
    let best = hit.best_hsp()
    match best {
      Some(hsp) => {
        if hsp.evalue <= threshold {
          filtered.push(hit)
        }
        ()
      }
      None => ()
    }
    i = i + 1
  }
  SearchIOQueryResult::{
    id: self.id, description: self.description,
    seq_length: self.seq_length, database: self.database,
    n_hits: filtered.length(), hits: filtered,
    total_hits: self.total_hits
  }
}

///|
/// Get the top N hits by score.
pub fn SearchIOQueryResult::top_n(self : SearchIOQueryResult, n : Int) -> Array[SearchIOHit] {
  let sorted = self.sort_by_score()
  if sorted.hits.length() <= n {
    return sorted.hits
  }
  let result : Array[SearchIOHit] = Array::new()
  let mut i = 0
  while i < n {
    result.push(sorted.hits[i])
    i = i + 1
  }
  result
}

///|
/// SearchIO iterator for parsing search results from various formats.
pub struct SearchIOIterator {
  /// Raw data lines (or XML/JSON content)
  lines : Array[String]
  /// Current position in the data
  position : Int
  /// Format of the search output
  format : SearchIOFormat
  /// Whether iterator has been initialized
  initialized : Bool
}

///|
/// Construct SearchIOIterator from search output text.
pub fn SearchIOIterator::new(content : String, format : SearchIOFormat) -> SearchIOIterator {
  let lines = split_lines(content)
  SearchIOIterator::{ lines, position: 0, format, initialized: false }
}

///|
/// Parse BLAST tabular output (6-column or extended tabular).
pub fn SearchIOIterator::parse_blast_tabular(self : SearchIOIterator) -> Array[SearchIOQueryResult] {
  let results : Array[SearchIOQueryResult] = Array::new()
  let mut current_query_id = ""
  let mut current_query : SearchIOQueryResult? = None

  let mut line_idx = 0
  while line_idx < self.lines.length() {
    let line = self.lines[line_idx]
    let trimmed : StringView = line.trim()
    if trimmed.length() == 0 {
      line_idx = line_idx + 1
      continue
    }
    if trimmed.has_prefix("#") {
      if trimmed.has_prefix("# Query:") || trimmed.has_prefix("# query:") {
        match current_query {
          Some(qr) => {
            let sorted = qr.sort_by_score()
            results.push(sorted)
            ()
          }
          None => ()
        }
        current_query_id = extract_query_id(trimmed)
        current_query = Some(SearchIOQueryResult::new(
          id=current_query_id, description="", seq_length=0, database=""
        ))
      }
      line_idx = line_idx + 1
      continue
    }

    let fields = split_tab(trimmed)
    if fields.length() >= 12 {
      match current_query {
        Some(qr) => {
          let hit_id = fields[1]
          let identity = searchio_parse_double(fields[2])
          let length = searchio_parse_int(fields[3])
          let qstart = searchio_parse_int(fields[6])
          let qend = searchio_parse_int(fields[7])
          let hstart = searchio_parse_int(fields[8])
          let hend = searchio_parse_int(fields[9])
          let evalue = searchio_parse_double(fields[10])
          let bitscore = searchio_parse_double(fields[11])

          let hsp = SearchIOHsp::new(
            bitscore=bitscore, evalue=evalue, identity=identity, positives=identity,
            gap=0.0, alignment_length=length,
            query_start=qstart, query_end=qend, hit_start=hstart, hit_end=hend
          )

          let mut found = false
          let mut i = 0
          while i < qr.hits.length() {
            if qr.hits[i].id == hit_id {
              let updated_hit = qr.hits[i].add_hsp(hsp)
              let new_hits = qr.hits.copy()
              new_hits[i] = updated_hit
              current_query = Some(SearchIOQueryResult::{
                id: qr.id, description: qr.description,
                seq_length: qr.seq_length, database: qr.database,
                n_hits: qr.n_hits, hits: new_hits,
                total_hits: qr.total_hits
              })
              found = true
              break
            }
            i = i + 1
          }
          if !found {
            let hit = SearchIOHit::new(id=hit_id, description="", seq_length=0).add_hsp(hsp)
            let new_hits = qr.hits.copy()
            new_hits.push(hit)
            current_query = Some(SearchIOQueryResult::{
              id: qr.id, description: qr.description,
              seq_length: qr.seq_length, database: qr.database,
              n_hits: qr.n_hits + 1, hits: new_hits,
              total_hits: qr.total_hits
            })
          }
          ()
        }
        None => ()
      }
    }
    line_idx = line_idx + 1
  }

  match current_query {
    Some(qr) => {
      let sorted = qr.sort_by_score()
      results.push(sorted)
      ()
    }
    None => ()
  }

  results
}

///|
/// Parse BLAST text output (legacy format).
pub fn SearchIOIterator::parse_blast_text(self : SearchIOIterator) -> Array[SearchIOQueryResult] {
  let results : Array[SearchIOQueryResult] = Array::new()
  let mut current_query : SearchIOQueryResult? = None

  let mut line_idx = 0
  while line_idx < self.lines.length() {
    let line = self.lines[line_idx]
    let trimmed : StringView = line.trim()
    if trimmed.length() == 0 {
      line_idx = line_idx + 1
      continue
    }

    if trimmed.has_prefix("Query=") || trimmed.has_prefix("Query:") {
      match current_query {
        Some(qr) => {
          results.push(qr.sort_by_score())
          ()
        }
        None => ()
      }
      let query_id = extract_query_id_text(trimmed)
      current_query = Some(SearchIOQueryResult::new(id=query_id, description="", seq_length=0, database=""))
      line_idx = line_idx + 1
      continue
    }

    if trimmed.has_prefix(">") {
      match current_query {
        Some(qr) => {
          let hit_id = trimmed[1:].trim().to_owned()
          let hit = SearchIOHit::new(id=hit_id, description="", seq_length=0)
          let new_hits = qr.hits.copy()
          new_hits.push(hit)
          current_query = Some(SearchIOQueryResult::{
            id: qr.id, description: qr.description,
            seq_length: qr.seq_length, database: qr.database,
            n_hits: qr.n_hits + 1, hits: new_hits,
            total_hits: qr.total_hits
          })
          ()
        }
        None => ()
      }
      line_idx = line_idx + 1
      continue
    }

    if trimmed.has_prefix("Score =") || trimmed.contains("Score =") {
      match current_query {
        Some(qr) => {
          let (bitscore, evalue) = extract_score_evalue(trimmed.to_owned())
          if qr.hits.length() > 0 {
            let idx = qr.hits.length() - 1
            let hsp = SearchIOHsp::new(
              bitscore=bitscore, evalue=evalue, identity=0.0, positives=0.0,
              gap=0.0, alignment_length=0,
              query_start=0, query_end=0, hit_start=0, hit_end=0
            )
            let updated_hit = qr.hits[idx].add_hsp(hsp)
            let new_hits = qr.hits.copy()
            new_hits[idx] = updated_hit
            current_query = Some(SearchIOQueryResult::{
              id: qr.id, description: qr.description,
              seq_length: qr.seq_length, database: qr.database,
              n_hits: qr.n_hits, hits: new_hits,
              total_hits: qr.total_hits
            })
          }
          ()
        }
        None => ()
      }
      line_idx = line_idx + 1
      continue
    }
    line_idx = line_idx + 1
  }

  match current_query {
    Some(qr) => {
      results.push(qr.sort_by_score())
      ()
    }
    None => ()
  }

  results
}

///|
/// Main entry point: parse search results based on format.
pub fn SearchIOIterator::parse(self : SearchIOIterator) -> Array[SearchIOQueryResult] {
  match self.format {
    SearchIOFormat::BlastTabular => self.parse_blast_tabular()
    SearchIOFormat::BlastText => self.parse_blast_text()
    SearchIOFormat::BlastXml => self.parse_blast_tabular()
    SearchIOFormat::HmmscanDomtbl => self.parse_blast_tabular()
    _ => Array::new()
  }
}

///|
/// Split content into lines by newline.
fn split_lines(content : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut current = ""
  let mut i = 0
  while i < content.length() {
    let c = content.unsafe_get(i)
    if c == '\n' {
      result.push(current)
      current = ""
    } else if c == '\r' {
      ()
    } else {
      current = current + u16_to_str(c)
    }
    i = i + 1
  }
  if current.length() > 0 {
    result.push(current)
  }
  result
}

///|
/// Split a line by tab characters.
fn split_tab(line : StringView) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut current = ""
  let mut i = 0
  while i < line.length() {
    let c = line.unsafe_get(i)
    if c == '\t' {
      result.push(current)
      current = ""
    } else {
      current = current + u16_to_str(c)
    }
    i = i + 1
  }
  result.push(current)
  result
}

///|
/// Extract query ID from a "# Query:" header line.
fn extract_query_id(line : StringView) -> String {
  if line.length() > 8 {
    let rest = line[8:].trim()
    if rest.length() > 0 {
      let mut i = 0
      while i < rest.length() {
        let c = rest.unsafe_get(i)
        if c == ' ' || c == '\t' {
          break
        }
        i = i + 1
      }
      if i > 0 {
        rest[0:i].to_owned()
      } else {
        "unknown"
      }
    } else {
      "unknown"
    }
  } else {
    "unknown"
  }
}

///|
/// Extract query ID from a "Query=" or "Query:" text line.
fn extract_query_id_text(line : StringView) -> String {
  if line.length() > 7 {
    let rest = line[7:].trim()
    let mut i = 0
    while i < rest.length() {
      let c = rest.unsafe_get(i)
      if c == ' ' || c == '\t' || c == '\n' {
        break
      }
      i = i + 1
    }
    if i > 0 {
      return rest[0:i].to_owned()
    }
  }
  "unknown"
}

///|
/// Extract bitscore and evalue from a Score line.
fn extract_score_evalue(line : String) -> (Double, Double) {
  let bitscore = find_numeric_value(line, "Score")
  let evalue = find_numeric_value(line, "Expect")
  (bitscore, evalue)
}

///|
/// Find the position of a keyword in a string, returns -1 if not found.
fn string_index_of(s : String, keyword : String) -> Int {
  let n = s.length()
  let kn = keyword.length()
  if kn == 0 || kn > n {
    return -1
  }
  let mut i = 0
  while i <= n - kn {
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != keyword.unsafe_get(j) {
        break
      }
      j = j + 1
    }
    if j == kn {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Find a numeric value after a keyword in a line.
fn find_numeric_value(line : String, keyword : String) -> Double {
  let idx = string_index_of(line, keyword)
  if idx < 0 {
    return 0.0
  }
  let mut i = idx + keyword.length()
  while i < line.length() {
    let c = line.unsafe_get(i)
    if c >= '0' && c <= '9' {
      break
    }
    i = i + 1
  }
  let mut num_str = ""
  let mut has_dot = false
  let mut has_e = false
  while i < line.length() {
    let c = line.unsafe_get(i)
    if (c >= '0' && c <= '9') || (c == '.' && !has_dot) || (c == 'e' || c == 'E') || (c == '-' && has_e) || (c == '+' && has_e) {
      if c == '.' { has_dot = true }
      if c == 'e' || c == 'E' { has_e = true; has_dot = false }
      num_str = num_str + u16_to_str(c)
    } else {
      break
    }
    i = i + 1
  }
  if num_str.length() > 0 {
    searchio_parse_double(num_str)
  } else {
    0.0
  }
}

///|
/// Create a mock BLAST tabular result string for testing.
pub fn search_io_mock_blast_tabular() -> String {
  let tab = "\t"
  "\n" +
  "# BLASTN 2.13.0+\n" +
  "# Query: seq1 test query\n" +
  "# Database: nr\n" +
  "# Fields: query id, subject id, % identity, alignment length, mismatches, open gaps, q.start, q.end, s.start, s.end, evalue, bit score\n" +
  "seq1" + tab + "hit1" + tab + "95.0" + tab + "100" + tab + "3" + tab + "2" + tab + "1" + tab + "100" + tab + "1" + tab + "100" + tab + "1e-50" + tab + "200.5\n" +
  "seq1" + tab + "hit2" + tab + "80.0" + tab + "80" + tab + "10" + tab + "5" + tab + "1" + tab + "80" + tab + "1" + tab + "80" + tab + "1e-30" + tab + "150.3\n" +
  "# Query: seq2 another query\n" +
  "# Database: nr\n" +
  "# Fields: query id, subject id, % identity, alignment length, mismatches, open gaps, q.start, q.end, s.start, s.end, evalue, bit score\n" +
  "seq2" + tab + "hit3" + tab + "90.0" + tab + "200" + tab + "5" + tab + "1" + tab + "1" + tab + "200" + tab + "1" + tab + "200" + tab + "1e-80" + tab + "300.0\n"
}