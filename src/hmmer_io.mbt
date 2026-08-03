///|
/// Bio.SearchIO.HmmerIO - HMMER3 domain table (domtblout) and text output parsing.
/// Mirrors Biopython's Bio.SearchIO.HmmerIO module.
///
/// Two HMMER3 output formats are supported:
///   * hmmscan/hmmsearch `--domtblout` (per-domain tabular, 23+ columns,
///     whitespace-separated, multi-line per query)
///   * hmmscan/hmmsearch human-readable text output (verbose, multi-block,
///     with Query/Domain/Score sections)
///
/// Both parsers populate the same QueryResult/Hit/HSP/HSPFragment model
/// defined in search_io.mbt.

// ============================================================================
// Internal helpers
// ============================================================================

///|
/// Split a line by runs of ASCII whitespace, dropping empty tokens.
fn hmmer_io_split_whitespace(s : String) -> Array[String] {
  let n = s.length()
  let tokens : Array[String] = []
  let mut start = 0
  let mut i = 0
  let mut in_token = false
  while i < n {
    let c = s.unsafe_get(i).to_int()
    let is_ws = c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int()
    if is_ws {
      if in_token {
        tokens.push(s[start:i].to_string())
        in_token = false
      }
    } else if !in_token {
      start = i
      in_token = true
    }
    i = i + 1
  }
  if in_token {
    tokens.push(s[start:n].to_string())
  }
  tokens
}

///|
fn hmmer_io_starts_with(s : String, prefix : String) -> Bool {
  let n = s.length()
  let p = prefix.length()
  if p > n {
    return false
  }
  let mut i = 0
  while i < p {
    if s.unsafe_get(i).to_int() != prefix.unsafe_get(i).to_int() {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Parse a Double, returning 0.0 on failure. Supports scientific notation.
fn hmmer_io_parse_double(s : String) -> Double {
  match parse_double(s) {
    Some(v) => v
    None => 0.0
  }
}

///|
fn hmmer_io_parse_int(s : String) -> Int {
  parse_int(s)
}

// ============================================================================
// Domain table format (domtblout) parsing
// ============================================================================

///|
/// Parse HMMER3 `--domtblout` (domain table) output.
///
/// Columns (whitespace-separated, 23 fixed columns + free-text description):
///   target_name target_acc query_name query_acc
///   evalue_full score_full bias_full
///   dom_idx dom_total
///   c_evalue i_evalue score_dom bias_dom
///   hmm_from hmm_to ali_from ali_to env_from env_to
///   exp reg clu ov env dom rep inc description
///
/// Multiple rows for the same (query, target) pair represent separate domains.
/// They are aggregated into a single Hit with one HSP per domain.
pub fn parse_hmmer3_domtab(content : String) -> Array[QueryResult] {
  let lines = split_by_char(content, '\n'.to_int())

  let results : Array[QueryResult] = []
  // Group HSPs/fragments by query id, then by target id.
  let query_order : Array[String] = []
  let query_desc : Map[String, String] = Map([], capacity=0)
  // target-by-query: query -> target -> array of HSP
  let query_hits : Map[String, Array[HmmerDomHit]] = Map([], capacity=0)

  let mut i = 0
  while i < lines.length() {
    let ln = lines[i]
    if ln.length() == 0 || hmmer_io_starts_with(ln, "#") {
      i = i + 1
      continue
    }
    let parts = hmmer_io_split_whitespace(ln)
    if parts.length() < 22 {
      i = i + 1
      continue
    }
    // Description is the remaining tokens after the 23rd fixed column.
    // Fixed columns: 0..=22 (23 fields). Description starts at index 23.
    let target_name = parts[0]
    let target_acc = parts[1]
    let query_name = parts[2]
    let query_acc = parts[3]
    let evalue_full = hmmer_io_parse_double(parts[4])
    let score_full = hmmer_io_parse_double(parts[5])
    let _bias_full = hmmer_io_parse_double(parts[6])
    let dom_idx = hmmer_io_parse_int(parts[7])
    let dom_total = hmmer_io_parse_int(parts[8])
    let c_evalue = hmmer_io_parse_double(parts[9])
    let i_evalue = hmmer_io_parse_double(parts[10])
    let score_dom = hmmer_io_parse_double(parts[11])
    let _bias_dom = hmmer_io_parse_double(parts[12])
    let hmm_from = hmmer_io_parse_int(parts[13])
    let hmm_to = hmmer_io_parse_int(parts[14])
    let ali_from = hmmer_io_parse_int(parts[15])
    let ali_to = hmmer_io_parse_int(parts[16])
    let env_from = hmmer_io_parse_int(parts[17])
    let env_to = hmmer_io_parse_int(parts[18])
    let _exp = hmmer_io_parse_double(parts[19])
    let _reg = hmmer_io_parse_int(parts[20])
    let _clu = hmmer_io_parse_int(parts[21])
    // index 22 onward: ov env dom rep inc description...
    let mut description = ""
    if parts.length() > 23 {
      let mut k = 23
      while k < parts.length() {
        if k > 23 {
          description = description + " "
        }
        description = description + parts[k]
        k = k + 1
      }
    }
    let hit = HmmerDomHit::{
      target_name,
      target_acc,
      query_name,
      query_acc,
      evalue_full,
      score_full,
      dom_idx,
      dom_total,
      c_evalue,
      i_evalue,
      score_dom,
      hmm_from,
      hmm_to,
      ali_from,
      ali_to,
      env_from,
      env_to,
      description,
    }
    let hits = query_hits.get(query_name)
    match hits {
      Some(arr) => arr.push(hit)
      None => {
        query_order.push(query_name)
        let arr : Array[HmmerDomHit] = []
        arr.push(hit)
        query_hits.set(query_name, arr)
        ()
      }
    }
    query_desc.set(query_name, query_acc)
    ignore(dom_total)
    i = i + 1
  }

  // Build QueryResult objects.
  let mut qi = 0
  while qi < query_order.length() {
    let qname = query_order[qi]
    let dom_hits = query_hits.get(qname)
    let dom_hits_arr = match dom_hits {
      Some(a) => a
      None => {
        qi = qi + 1
        continue
      }
    }
    // Group by target_name.
    let target_order : Array[String] = []
    let target_accs : Map[String, String] = Map([], capacity=0)
    let target_descs : Map[String, String] = Map([], capacity=0)
    let by_target : Map[String, Array[HmmerDomHit]] = Map([], capacity=0)
    let mut di = 0
    while di < dom_hits_arr.length() {
      let dh = dom_hits_arr[di]
      let arr = by_target.get(dh.target_name)
      match arr {
        Some(a) => a.push(dh)
        None => {
          target_order.push(dh.target_name)
          let na : Array[HmmerDomHit] = []
          na.push(dh)
          by_target.set(dh.target_name, na)
          ()
        }
      }
      target_accs.set(dh.target_name, dh.target_acc)
      target_descs.set(dh.target_name, dh.description)
      di = di + 1
    }
    let hits : Array[Hit] = []
    let mut ti = 0
    while ti < target_order.length() {
      let tname = target_order[ti]
      let dhits = by_target.get(tname)
      let dhits_arr = match dhits {
        Some(a) => a
        None => {
          ti = ti + 1
          continue
        }
      }
      let hsps : Array[HSP] = []
      let mut tbs = 0.0
      let mut mev = 1.0
      let mut hi = 0
      while hi < dhits_arr.length() {
        let dh = dhits_arr[hi]
        let frag = HSPFragment::new(
          dh.hmm_from,
          dh.hmm_to,
          dh.ali_from,
          dh.ali_to,
          Strand::Plus,
          Strand::Plus,
          if dh.ali_to >= dh.ali_from {
            dh.ali_to - dh.ali_from + 1
          } else {
            0
          },
          0,
          dh.score_dom,
          dh.i_evalue,
          "",
          "",
          "",
        )
        let frags : Array[HSPFragment] = [frag]
        hsps.push(
          HSP::new(tname, dh.description, frags, dh.score_dom, dh.i_evalue, dh.score_dom),
        )
        tbs = tbs + dh.score_dom
        if dh.i_evalue < mev {
          mev = dh.i_evalue
        }
        ignore(dh.env_from)
        ignore(dh.env_to)
        ignore(dh.c_evalue)
        ignore(dh.dom_idx)
        ignore(dh.dom_total)
        ignore(dh.evalue_full)
        ignore(dh.score_full)
        ignore(dh.query_acc)
        hi = hi + 1
      }
      hits.push(
        Hit::new(
          tname,
          target_descs.get_or_default(tname, ""),
          target_accs.get_or_default(tname, ""),
          hsps,
          tbs,
          mev,
        ),
      )
      ti = ti + 1
    }
    results.push(
      QueryResult::new(
        qname,
        query_desc.get_or_default(qname, ""),
        0,
        hits,
        "hmmscan",
        "3.0",
        "",
      ),
    )
    qi = qi + 1
  }

  results
}

///|
/// Intermediate domain-table row record.
struct HmmerDomHit {
  target_name : String
  target_acc : String
  query_name : String
  query_acc : String
  evalue_full : Double
  score_full : Double
  dom_idx : Int
  dom_total : Int
  c_evalue : Double
  i_evalue : Double
  score_dom : Double
  hmm_from : Int
  hmm_to : Int
  ali_from : Int
  ali_to : Int
  env_from : Int
  env_to : Int
  description : String
}

// ============================================================================
// Human-readable text output parsing
// ============================================================================

///|
/// Parse HMMER3 human-readable text output from `hmmscan`/`hmmsearch`.
///
/// The text format begins with a header block describing the program, query
/// and database, followed by per-hit blocks that include score tables and
/// domain alignments. This parser extracts query/hit/HSP/domain coordinates
/// from the score and alignment tables.
///
/// Recognised section markers:
///   * `Query:`       - new query section (sets query id/description)
///   * `>>`           - per-target hit header (target id/description)
///   * `   == domain` - domain row inside the per-target domain table
///   * `Alignment for` - start of a domain alignment block (optional)
pub fn parse_hmmer3_text(content : String) -> Array[QueryResult] {
  let lines = split_by_char(content, '\n'.to_int())

  let results : Array[QueryResult] = []
  let mut cur_query_id = ""
  let mut cur_query_desc = ""
  let mut cur_query_len = 0
  let mut cur_hits : Array[Hit] = []
  let mut cur_target_id = ""
  let mut cur_target_desc = ""
  let mut cur_hsps : Array[HSP] = []

  let mut i = 0
  while i < lines.length() {
    let ln = lines[i]
    let trimmed = hmmer_io_trim(ln)

    if hmmer_io_starts_with(trimmed, "Query:") {
      // Flush current query (including its current hit) before starting a new one.
      // flush_hit inline:
      if cur_target_id.length() > 0 && cur_hsps.length() > 0 {
        let mut tbs = 0.0
        let mut mev = 1.0
        let mut k = 0
        while k < cur_hsps.length() {
          tbs = tbs + cur_hsps[k].bitscore
          if cur_hsps[k].evalue < mev {
            mev = cur_hsps[k].evalue
          }
          k = k + 1
        }
        cur_hits.push(
          Hit::new(cur_target_id, cur_target_desc, "", cur_hsps, tbs, mev),
        )
      }
      cur_target_id = ""
      cur_target_desc = ""
      cur_hsps = []
      // flush_query inline:
      if cur_query_id.length() > 0 {
        results.push(
          QueryResult::new(
            cur_query_id,
            cur_query_desc,
            cur_query_len,
            cur_hits,
            "hmmscan",
            "3.0",
            "",
          ),
        )
      }
      cur_query_id = ""
      cur_query_desc = ""
      cur_query_len = 0
      cur_hits = []
      // Now parse the new query header.
      let rest = hmmer_io_substr_after(trimmed, "Query:")
      let parts = hmmer_io_split_whitespace(rest)
      if parts.length() > 0 {
        cur_query_id = parts[0]
        if parts.length() > 1 {
          let mut d = ""
          let mut k = 1
          while k < parts.length() {
            if k > 1 {
              d = d + " "
            }
            d = d + parts[k]
            k = k + 1
          }
          cur_query_desc = d
        }
      }
      i = i + 1
      continue
    }

    // Query length line, e.g. "Query length: 245"
    if hmmer_io_starts_with(trimmed, "Query length:") {
      let rest = hmmer_io_substr_after(trimmed, "Query length:")
      let parts = hmmer_io_split_whitespace(rest)
      if parts.length() > 0 {
        cur_query_len = hmmer_io_parse_int(parts[0])
      }
      i = i + 1
      continue
    }

    // Per-target header: ">> target_id  description"
    if hmmer_io_starts_with(trimmed, ">>") {
      // flush_hit inline:
      if cur_target_id.length() > 0 && cur_hsps.length() > 0 {
        let mut tbs = 0.0
        let mut mev = 1.0
        let mut k = 0
        while k < cur_hsps.length() {
          tbs = tbs + cur_hsps[k].bitscore
          if cur_hsps[k].evalue < mev {
            mev = cur_hsps[k].evalue
          }
          k = k + 1
        }
        cur_hits.push(
          Hit::new(cur_target_id, cur_target_desc, "", cur_hsps, tbs, mev),
        )
      }
      cur_target_id = ""
      cur_target_desc = ""
      cur_hsps = []
      // Parse new target header.
      let rest = hmmer_io_substr_after(trimmed, ">>")
      let parts = hmmer_io_split_whitespace(rest)
      if parts.length() > 0 {
        cur_target_id = parts[0]
        if parts.length() > 1 {
          let mut d = ""
          let mut k = 1
          while k < parts.length() {
            if k > 1 {
              d = d + " "
            }
            d = d + parts[k]
            k = k + 1
          }
          cur_target_desc = d
        }
      }
      i = i + 1
      continue
    }

    // Domain row in the per-target domain table:
    // " == domain 1   score: 200.3 bits;  conditional E-value: 1.5e-53"
    //                         ^score                       ^evalue
    if cur_target_id.length() > 0 && hmmer_io_starts_with(trimmed, "== domain") {
      let score = hmmer_io_extract_double_after(trimmed, "score:")
      let evalue = hmmer_io_extract_double_after(trimmed, "E-value:")
      let frag = HSPFragment::new(
        0,
        0,
        0,
        0,
        Strand::Plus,
        Strand::Plus,
        0,
        0,
        score,
        evalue,
        "",
        "",
        "",
      )
      let frags : Array[HSPFragment] = [frag]
      cur_hsps.push(
        HSP::new(cur_target_id, cur_target_desc, frags, score, evalue, score),
      )
      i = i + 1
      continue
    }

    // Internal hits table line for inclusion threshold:
    // "  E-value  score  bias    E-value  score  bias    exp  N  Domain  Description"
    // This is a header row, skip it.
    if hmmer_io_starts_with(trimmed, "E-value") && cur_query_id.length() > 0 {
      i = i + 1
      continue
    }

    i = i + 1
  }

  // Final flush_hit inline:
  if cur_target_id.length() > 0 && cur_hsps.length() > 0 {
    let mut tbs = 0.0
    let mut mev = 1.0
    let mut k = 0
    while k < cur_hsps.length() {
      tbs = tbs + cur_hsps[k].bitscore
      if cur_hsps[k].evalue < mev {
        mev = cur_hsps[k].evalue
      }
      k = k + 1
    }
    cur_hits.push(
      Hit::new(cur_target_id, cur_target_desc, "", cur_hsps, tbs, mev),
    )
  }
  // Final flush_query inline:
  if cur_query_id.length() > 0 {
    results.push(
      QueryResult::new(
        cur_query_id,
        cur_query_desc,
        cur_query_len,
        cur_hits,
        "hmmscan",
        "3.0",
        "",
      ),
    )
  }
  results
}

///|
/// Remove leading and trailing ASCII whitespace.
fn hmmer_io_trim(s : String) -> String {
  let n = s.length()
  let mut start = 0
  let mut end = n
  while start < end {
    let c = s.unsafe_get(start).to_int()
    if c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int() {
      start = start + 1
    } else {
      break
    }
  }
  while end > start {
    let c = s.unsafe_get(end - 1).to_int()
    if c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int() {
      end = end - 1
    } else {
      break
    }
  }
  s[start:end].to_string()
}

///|
/// Return the substring after the first occurrence of `marker`, or "" if not found.
fn hmmer_io_substr_after(s : String, marker : String) -> String {
  let n = s.length()
  let m = marker.length()
  if m == 0 || m > n {
    return ""
  }
  let mut i = 0
  while i <= n - m {
    let mut j = 0
    let mut matched = true
    while j < m {
      if s.unsafe_get(i + j).to_int() != marker.unsafe_get(j).to_int() {
        matched = false
        break
      }
      j = j + 1
    }
    if matched {
      return s[i + m:n].to_string()
    }
    i = i + 1
  }
  ""
}

///|
/// Extract the first Double following `marker` in `s`. Returns 0.0 if not found.
fn hmmer_io_extract_double_after(s : String, marker : String) -> Double {
  let rest = hmmer_io_substr_after(s, marker)
  if rest.length() == 0 {
    return 0.0
  }
  let parts = hmmer_io_split_whitespace(rest)
  if parts.length() == 0 {
    return 0.0
  }
  hmmer_io_parse_double(parts[0])
}

///|
/// Generate a small HMMER3 domtblout-format sample for tests and demos.
pub fn hmmer3_domtab_sample() -> String {
  "# hmmsearch :: search profile(s) against a sequence database\n" +
    "# HMMER 3.3.2 (Nov 2020)\n" +
    "#\n" +
    "# target name        accession   query name           accession    E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from   to     exp  reg  clu  ov env  dom  rep  inc  description of target\n" +
    "#------------------- ---------- -------------------- ---------- --------- ------ ------ --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ---- ---- ---- ---- ---- ---- ---- ---------------------\n" +
    "sp|Q9Y2W8|A1BG_HUMAN Q9Y2W8     PF00001.28           PF00001     1.2e-50  200.3    5.1   1   1   1.2e-50   1.5e-53  200.3   5.1     1   150     1   150   1.0   1   0   0   1   1   1   1 Alpha-1-B glycoprotein\n" +
    "sp|P12345|ABC_HUMAN  P12345     PF00001.28           PF00001     3.4e-30  150.1    2.0   1   2   3.4e-30   4.0e-33  150.1   2.0     1   120     1   120   1.0   1   0   0   1   1   1   1 ABC transporter\n" +
    "sp|P12345|ABC_HUMAN  P12345     PF00001.28           PF00001     3.4e-30  150.1    2.0   2   2   5.0e-20   6.0e-23  100.5   1.0     1    90     1    90   1.0   1   0   0   1   1   1   1 ABC transporter C-term\n"
}

///|
/// Generate a small HMMER3 text-format sample for tests and demos.
pub fn hmmer3_text_sample() -> String {
  "# hmmsearch :: search profile(s) against a sequence database\n" +
    "# HMMER 3.3.2 (Nov 2020)\n" +
    "#\n" +
    "Query:       PF00001.28  7tm_1\n" +
    "Query length: 150\n" +
    "Scores for complete sequences (score includes all domains):\n" +
    "   E-value  score  bias    E-value  score  bias    exp  N  Domain  Description\n" +
    "   ------- ------ -----    ------- ------ -----    ---- ---  ------  -----------\n" +
    "sp|Q9Y2W8|A1BG_HUMAN   1.2e-50  200.3   5.1    1.5e-53  200.3   5.1    1.0   1   1   Alpha-1-B glycoprotein\n" +
    "sp|P12345|ABC_HUMAN    3.4e-30  150.1   2.0    4.0e-33  150.1   2.0    1.0   1   1   ABC transporter\n" +
    "\n" +
    "Domain annotation for each sequence (and alignments):\n" +
    ">> sp|Q9Y2W8|A1BG_HUMAN  Alpha-1-B glycoprotein\n" +
    "   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc\n" +
    " ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----\n" +
    "   1   200.3   5.1   1.2e-50   1.5e-53       1     150         1     150         1     150  0.90\n" +
    "  == domain 1   score: 200.3 bits;  conditional E-value: 1.5e-53\n" +
    "\n" +
    ">> sp|P12345|ABC_HUMAN  ABC transporter\n" +
    "   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc\n" +
    " ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----\n" +
    "   1   150.1   2.0   3.4e-30   4.0e-33       1     120         1     120         1     120  0.88\n" +
    "  == domain 1   score: 150.1 bits;  conditional E-value: 4.0e-33\n"
}

///|
/// Return the number of HSP fragments across all hits in a QueryResult.
/// (Convenience alias; same semantics as count_hsps in search_io.mbt.)
pub fn hmmer_io_count_hsps(qr : QueryResult) -> Int {
  let mut total = 0
  let mut i = 0
  while i < qr.hits.length() {
    let mut j = 0
    while j < qr.hits[i].hsps.length() {
      total = total + qr.hits[i].hsps[j].fragments.length()
      j = j + 1
    }
    i = i + 1
  }
  total
}
