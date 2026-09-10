///|
/// Bio.SearchIO.BlatIO - parse and write BLAT PSL / PSLX output, building the
/// standard SearchIO QueryResult/Hit/HSP model.
///
/// Ported 1:1 from Biopython Bio/SearchIO/BlatIO.py:
/// - 21-column PSL and 23-column PSLX (with per-block sequences) rows,
/// - optional `psLayout version 3` header is skipped while parsing,
/// - minus-strand block starts are reoriented to plus-strand coordinates,
/// - protein versus nucleotide queries are auto-detected (hit blocks then
///   span 3 nucleotides per residue),
/// - UCSC millibad identity and score formulas,
/// - rows are grouped into QueryResult / Hit objects exactly like Biopython
///   (consecutive hit rows; HSPs of hits seen earlier are merged back in).
///
/// Coordinates follow the Biopython convention: 0-based, half-open.

///|
/// One aligned block (fragment) of a BLAT HSP.
pub struct BlatBlock {
  query_start : Int
  query_end : Int
  hit_start : Int
  hit_end : Int
  /// -1 / 0 / 1 for minus / unknown(protein) / plus
  query_strand : Int
  hit_strand : Int
  query_seq : String
  hit_seq : String
}

///|
/// One PSL row as a SearchIO-style HSP with block-level fragments.
pub struct BlatHsp {
  blocks : Array[BlatBlock]
  match_num : Int
  mismatch_num : Int
  match_rep_num : Int
  n_num : Int
  query_gapopen_num : Int
  query_gap_num : Int
  hit_gapopen_num : Int
  hit_gap_num : Int
  ident_num : Int
  gapopen_num : Int
  gap_num : Int
  query_is_protein : Bool
  ident_pct : Double
  score : Int
  /// Whether the source row carried a hit strand character
  has_hit_strand : Bool
}

///|
/// One target sequence hit, with one or more HSP rows.
pub struct BlatHit {
  id : String
  seq_len : Int
  hsps : Array[BlatHsp]
}

///|
/// All BLAT hits reported for one query sequence.
pub struct BlatQueryResult {
  id : String
  seq_len : Int
  program : String
  hits : Array[BlatHit]
}

///|
/// Writer counters (queries, hits, HSPs, fragments written).
pub struct BlatWriteStats {
  qresults : Int
  hits : Int
  hsps : Int
  fragments : Int
}

///|
/// Raw values of one PSL/PSLX row, exactly as written in the file.
priv struct BlatRow {
  qname : String
  qsize : Int
  tname : String
  tsize : Int
  matches : Int
  mismatches : Int
  repmatches : Int
  ncount : Int
  qnuminsert : Int
  qbaseinsert : Int
  tnuminsert : Int
  tbaseinsert : Int
  strand : String
  qstart : Int
  qend : Int
  tstart : Int
  tend : Int
  blocksizes : Array[Int]
  qstarts : Array[Int]
  tstarts : Array[Int]
  qseqs : Array[String]
  tseqs : Array[String]
  has_seqs : Bool
}

///|
fn blat_log_nat(x : Double) -> Double {
  if x <= 0.0 {
    @double.not_a_number
  } else if x == 1.0 {
    0.0
  } else {
    let y = (x - 1.0) / (x + 1.0)
    let y2 = y * y
    let mut sum = y
    let mut term = y
    let mut i = 1
    while i <= 20 {
      term = term * y2
      sum = sum + term / (2 * i + 1).to_double()
      i = i + 1
    }
    2.0 * sum
  }
}

///|
fn blat_is_int(s : String) -> Bool {
  if s.length() == 0 {
    return false
  }
  let zero = '0'.to_int()
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 {
      return false
    }
    i = i + 1
  }
  true
}

///|
fn blat_trim(s : String) -> String {
  s.trim().to_owned()
}

///|
/// Split a tab-separated row into non-empty fields (Biopython filters empties).
fn blat_fields(line : String) -> Array[String] {
  let fields : Array[String] = []
  for field in line.split("\t") {
    if field.length() > 0 {
      fields.push(field.to_owned())
    }
  }
  fields
}

///|
/// Integer cast of a comma-separated value list ("10,20,").
fn blat_csv_ints(csv : String) -> Array[Int] {
  let out : Array[Int] = []
  for item in csv.split(",") {
    if item.length() > 0 {
      out.push(parse_int(item.trim().to_owned()))
    }
  }
  out
}

///|
/// String cast of a comma-separated value list.
fn blat_csv_strs(csv : String) -> Array[String] {
  let out : Array[String] = []
  for item in csv.split(",") {
    if item.length() > 0 {
      out.push(item.to_owned())
    }
  }
  out
}

///|
/// Row check: first four tab fields must be integers (matches/mis/rep/N).
fn blat_is_data_row(line : String) -> Bool {
  let trimmed = blat_trim(line)
  if trimmed.length() == 0 {
    return false
  }
  let cols = blat_fields(trimmed)
  if cols.length() < 4 {
    return false
  }
  blat_is_int(cols[0]) &&
  blat_is_int(cols[1]) &&
  blat_is_int(cols[2]) &&
  blat_is_int(cols[3])
}

///|
/// Reorient minus-strand block starts into plus-strand coordinates.
/// plus_coord = length - minus_coord - block_size
fn blat_reorient_starts(
  starts : Array[Int],
  block_sizes : Array[Int],
  seq_len : Int,
  strand : Int,
) -> Array[Int] raise SeqError {
  if starts.length() != block_sizes.length() {
    raise SeqError(
      "Unequal start coordinates and block sizes list (\{starts.length()} vs \{block_sizes.length()})",
    )
  }
  let out : Array[Int] = []
  if strand >= 0 {
    for start in starts {
      out.push(start)
    }
  } else {
    let mut i = 0
    while i < starts.length() {
      out.push(seq_len - starts[i] - block_sizes[i])
      i = i + 1
    }
  }
  out
}

///|
/// Detect whether a row is a protein query, using UCSC's last-block test.
fn blat_is_protein(psl : BlatRow) -> Bool {
  if psl.strand.length() == 2 {
    let hit_strand = psl.strand.unsafe_get(1)
    let last_tstart = psl.tstarts[psl.tstarts.length() - 1]
    let last_span = 3 * psl.blocksizes[psl.blocksizes.length() - 1]
    if hit_strand == '+' {
      return psl.tend == last_tstart + last_span
    } else if hit_strand == '-' {
      return psl.tstart == psl.tsize - (last_tstart + last_span)
    } else {
      return false
    }
  }
  false
}

///|
/// UCSC millibad (thousandths of mismatching bases), protein-aware.
fn blat_calc_millibad(psl : BlatRow, is_protein : Bool) -> Double {
  let size_mul : Int = if is_protein { 3 } else { 1 }
  let qali_size = size_mul * (psl.qend - psl.qstart)
  let tali_size = psl.tend - psl.tstart
  let ali_size = if qali_size < tali_size { qali_size } else { tali_size }
  if ali_size <= 0 {
    return 0.0
  }
  let mut size_dif = qali_size - tali_size
  if size_dif < 0 {
    size_dif = 0
  }
  let total = size_mul * (psl.matches + psl.repmatches + psl.mismatches)
  if total == 0 {
    return 0.0
  }
  // Python round() is half-to-even; the log argument never hits a .5 case in
  // practice, floor(x+0.5) is equivalent here.
  let log_term = @math.floor(
    3.0 * blat_log_nat((1 + size_dif).to_double()) + 0.5,
  ).to_int()
  let numerator = psl.mismatches * size_mul + psl.qnuminsert + log_term
  1000.0 * numerator.to_double() / total.to_double()
}

///|
/// UCSC BLAT score: matches + repMatches/2 - mismatches - inserts.
fn blat_calc_score(psl : BlatRow, is_protein : Bool) -> Int {
  let size_mul : Int = if is_protein { 3 } else { 1 }
  size_mul * (psl.matches + (psl.repmatches >> 1)) -
  size_mul * psl.mismatches -
  psl.qnuminsert -
  psl.tnuminsert
}

///|
/// Parse one PSL/PSLX line into a BlatRow.
fn blat_parse_row(line : String, pslx : Bool) -> BlatRow raise SeqError {
  let raw = blat_trim(line)
  let cols = blat_fields(raw)
  let expected = if pslx { 23 } else { 21 }
  let fmt = if pslx { "PSLX" } else { "PSL" }
  if cols.length() != expected {
    raise SeqError(
      "Invalid \{fmt} line: \{raw}. Expected \{expected} tab-separated columns, found \{cols.length()}",
    )
  }
  let blocksizes = blat_csv_ints(cols[18])
  let qstarts = blat_csv_ints(cols[19])
  let tstarts = blat_csv_ints(cols[20])
  let (qseqs, tseqs, has_seqs) = if pslx {
    (blat_csv_strs(cols[21]), blat_csv_strs(cols[22]), true)
  } else {
    ([], [], false)
  }
  let blockcount = parse_int(cols[17])
  if blockcount != blocksizes.length() ||
    blockcount != qstarts.length() ||
    blockcount != tstarts.length() {
    raise SeqError(
      "Invalid \{fmt} line: \{raw}. blockCount \{blockcount} does not match block coordinate lists",
    )
  }
  {
    qname: cols[9],
    qsize: parse_int(cols[10]),
    tname: cols[13],
    tsize: parse_int(cols[14]),
    matches: parse_int(cols[0]),
    mismatches: parse_int(cols[1]),
    repmatches: parse_int(cols[2]),
    ncount: parse_int(cols[3]),
    qnuminsert: parse_int(cols[4]),
    qbaseinsert: parse_int(cols[5]),
    tnuminsert: parse_int(cols[6]),
    tbaseinsert: parse_int(cols[7]),
    strand: cols[8],
    qstart: parse_int(cols[11]),
    qend: parse_int(cols[12]),
    tstart: parse_int(cols[15]),
    tend: parse_int(cols[16]),
    blocksizes,
    qstarts,
    tstarts,
    qseqs,
    tseqs,
    has_seqs,
  }
}

///|
/// Build the SearchIO HSP (with per-block fragments) from one row.
fn blat_create_hsp(psl : BlatRow) -> BlatHsp raise SeqError {
  let is_protein = blat_is_protein(psl)
  // protein queries have no query strand (0)
  let qstrand : Int = if is_protein {
    0
  } else if psl.strand.unsafe_get(0) == '+' {
    1
  } else {
    -1
  }
  // hit strand defaults to plus if absent
  let has_hit_strand = psl.strand.length() == 2
  let hstrand : Int = if has_hit_strand && psl.strand.unsafe_get(1) == '-' {
    -1
  } else {
    1
  }
  let block_mul : Int = if is_protein { 3 } else { 1 }
  let qstarts = blat_reorient_starts(
    psl.qstarts,
    psl.blocksizes,
    psl.qsize,
    qstrand,
  )
  let hstarts = if has_hit_strand {
    let hit_sizes : Array[Int] = []
    for size in psl.blocksizes {
      hit_sizes.push(block_mul * size)
    }
    blat_reorient_starts(psl.tstarts, hit_sizes, psl.tsize, hstrand)
  } else {
    psl.tstarts.copy()
  }
  if qstarts.length() != hstarts.length() ||
    qstarts.length() != psl.blocksizes.length() {
    raise SeqError("PSL block coordinate lists have unequal lengths")
  }
  let blocks : Array[BlatBlock] = []
  let mut i = 0
  while i < qstarts.length() {
    let qseq = if psl.has_seqs && i < psl.qseqs.length() {
      psl.qseqs[i]
    } else {
      ""
    }
    let hseq = if psl.has_seqs && i < psl.tseqs.length() {
      psl.tseqs[i]
    } else {
      ""
    }
    blocks.push({
      query_start: qstarts[i],
      query_end: qstarts[i] + psl.blocksizes[i],
      hit_start: hstarts[i],
      hit_end: hstarts[i] + block_mul * psl.blocksizes[i],
      query_strand: qstrand,
      hit_strand: hstrand,
      query_seq: qseq,
      hit_seq: hseq,
    })
    i = i + 1
  }
  let query_start = blat_range_start(blocks, true)
  let query_end = blat_range_end(blocks, true)
  let hit_start = blat_range_start(blocks, false)
  let hit_end = blat_range_end(blocks, false)
  if query_start != psl.qstart || query_end != psl.qend {
    raise SeqError(
      "PSL query coordinates inconsistent with block starts: \{query_start}-\{query_end} vs \{psl.qstart}-\{psl.qend}",
    )
  }
  if hit_start != psl.tstart || hit_end != psl.tend {
    raise SeqError(
      "PSL hit coordinates inconsistent with block starts: \{hit_start}-\{hit_end} vs \{psl.tstart}-\{psl.tend}",
    )
  }
  {
    blocks,
    match_num: psl.matches,
    mismatch_num: psl.mismatches,
    match_rep_num: psl.repmatches,
    n_num: psl.ncount,
    query_gapopen_num: psl.qnuminsert,
    query_gap_num: psl.qbaseinsert,
    hit_gapopen_num: psl.tnuminsert,
    hit_gap_num: psl.tbaseinsert,
    ident_num: psl.matches + psl.repmatches,
    gapopen_num: psl.qnuminsert + psl.tnuminsert,
    gap_num: psl.qbaseinsert + psl.tbaseinsert,
    query_is_protein: is_protein,
    ident_pct: 100.0 - blat_calc_millibad(psl, is_protein) * 0.1,
    score: blat_calc_score(psl, is_protein),
    has_hit_strand,
  }
}

///|
fn blat_range_start(blocks : Array[BlatBlock], query : Bool) -> Int {
  let mut best = blocks[0]
  for block in blocks {
    let value = if query { block.query_start } else { block.hit_start }
    let best_value = if query { best.query_start } else { best.hit_start }
    if value < best_value {
      best = block
    }
  }
  if query {
    best.query_start
  } else {
    best.hit_start
  }
}

///|
fn blat_range_end(blocks : Array[BlatBlock], query : Bool) -> Int {
  let mut best = blocks[0]
  for block in blocks {
    let value = if query { block.query_end } else { block.hit_end }
    let best_value = if query { best.query_end } else { best.hit_end }
    if value > best_value {
      best = block
    }
  }
  if query {
    best.query_end
  } else {
    best.hit_end
  }
}

///|
fn blat_merge_hit(hits : Array[BlatHit], hit : BlatHit) -> Unit {
  let mut i = 0
  while i < hits.length() {
    if hits[i].id == hit.id {
      for hsp in hit.hsps {
        hits[i].hsps.push(hsp)
      }
      return
    }
    i = i + 1
  }
  hits.push(hit)
}

///|
/// Parse BLAT PSL (`pslx=false`, 21 columns) or PSLX (`pslx=true`,
/// 23 columns) text into one BlatQueryResult per query, mirroring
/// `Bio.SearchIO.parse(handle, "blat-psl")`. Header lines are skipped.
pub fn blat_psl_parse(
  content : String,
  pslx? : Bool = false,
) -> Array[BlatQueryResult] raise SeqError {
  let lines : Array[String] = []
  for line in content.split("\n") {
    lines.push(line.to_owned())
  }
  // locate the first data row
  let mut first = 0
  let mut found = false
  while first < lines.length() {
    if blat_is_data_row(lines[first]) {
      found = true
      break
    }
    first = first + 1
  }
  if !found {
    return []
  }
  // a trailing empty line (file ends with "\n") is EOF, like Python readline
  let mut end = lines.length()
  if end > 0 && lines[end - 1].length() == 0 {
    end = end - 1
  }
  // grouping state machine, same states as BlatIO._parse_qresult
  let qresults : Array[BlatQueryResult] = []
  let mut hit_list : Array[BlatHit] = []
  let mut hsp_list : Array[BlatHsp] = []
  let mut prev_row : BlatRow? = None
  let mut idx = first
  while idx <= end {
    // at EOF, emit with None state; otherwise parse the current row
    let cur_row : BlatRow? = if idx < end {
      Some(blat_parse_row(lines[idx], pslx))
    } else {
      None
    }
    let new_qresult = match (prev_row, cur_row) {
      (Some(p), Some(c)) => p.qname != c.qname
      (_, _) => true
    }
    let new_hit = match (prev_row, cur_row) {
      (Some(p), Some(c)) => p.tname != c.tname || new_qresult
      (_, _) => true
    }
    match prev_row {
      Some(prev) => {
        // copy identifying fields before the row is consumed by create_hsp
        let prev_hid = "\{prev.tname}"
        let prev_qid = "\{prev.qname}"
        let hit_size = prev.tsize
        let qsize = prev.qsize
        let hsp = blat_create_hsp(prev)
        hsp_list.push(hsp)
        if new_hit {
          let hit : BlatHit = {
            id: prev_hid,
            seq_len: hit_size,
            hsps: hsp_list,
          }
          blat_merge_hit(hit_list, hit)
          hsp_list = []
        }
        if new_qresult || idx == end {
          let qr : BlatQueryResult = {
            id: prev_qid,
            seq_len: qsize,
            program: "blat",
            hits: hit_list,
          }
          qresults.push(qr)
          if idx == end {
            return qresults
          }
          hit_list = []
        }
      }
      None => ()
    }
    prev_row = cur_row
    idx = idx + 1
  }
  qresults
}

///|
/// Build the standard psLayout version 3 header emitted by the writer.
pub fn blat_psl_header() -> String {
  let dashes = StringBuilder()
  let mut i = 0
  while i < 159 {
    dashes.write_char('-')
    i = i + 1
  }
  "psLayout version 3\n" +
  "\nmatch\tmis- \trep. \tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ        \tQ   \tQ    \tQ  \tT        \tT   \tT    \tT  \tblock\tblockSizes \tqStarts\t tStarts" +
  "\n     \tmatch\tmatch\t   \tcount\tbases\tcount\tbases\t      \tname     \tsize\tstart\tend\tname     \tsize\tstart\tend\tcount\n" +
  dashes.to_string() +
  "\n"
}

///|
fn blat_csv_join_ints(values : Array[Int]) -> String {
  let sb = StringBuilder()
  let mut i = 0
  while i < values.length() {
    if i > 0 {
      sb.write_char(',')
    }
    sb.write_string(values[i].to_string())
    i = i + 1
  }
  sb.write_char(',')
  sb.to_string()
}

///|
fn blat_csv_join_strs(values : Array[String]) -> String {
  let sb = StringBuilder()
  let mut i = 0
  while i < values.length() {
    if i > 0 {
      sb.write_char(',')
    }
    sb.write_string(values[i])
    i = i + 1
  }
  sb.write_char(',')
  sb.to_string()
}

///|
fn blat_hsp_query_span_all(hsp : BlatHsp) -> Array[Int] {
  let out : Array[Int] = []
  for block in hsp.blocks {
    out.push(block.query_end - block.query_start)
  }
  out
}

///|
fn blat_hsp_hit_span_all(hsp : BlatHsp) -> Array[Int] {
  let out : Array[Int] = []
  for block in hsp.blocks {
    out.push(block.hit_end - block.hit_start)
  }
  out
}

///|
fn blat_build_row(
  qresult : BlatQueryResult,
  hit : BlatHit,
  hsp : BlatHsp,
  pslx : Bool,
) -> String raise SeqError {
  let block_mul : Int = if hsp.query_is_protein { 3 } else { 1 }
  let query_spans = blat_hsp_query_span_all(hsp)
  let hit_spans = blat_hsp_hit_span_all(hsp)
  // span check: hit block spans must be blocksize_multiplier * query spans
  let mut i = 0
  while i < query_spans.length() {
    if hit_spans[i] != block_mul * query_spans[i] {
      raise SeqError("HSP hit span and query span values do not match.")
    }
    i = i + 1
  }
  let qstrand0 = hsp.blocks[0].query_strand
  let mut strand = if qstrand0 >= 0 { "+" } else { "-" }
  let qstart_list : Array[Int] = []
  for block in hsp.blocks {
    qstart_list.push(block.query_start)
  }
  let qstarts = blat_reorient_starts(
    qstart_list,
    query_spans,
    qresult.seq_len,
    qstrand0,
  )
  let hstrand0 = hsp.blocks[0].hit_strand
  if hstrand0 == 1 {
    if hsp.has_hit_strand {
      strand = strand + "+"
    }
  } else {
    strand = strand + "-"
  }
  let hstart_list : Array[Int] = []
  for block in hsp.blocks {
    hstart_list.push(block.hit_start)
  }
  let hstarts = blat_reorient_starts(
    hstart_list,
    hit_spans,
    hit.seq_len,
    hstrand0,
  )
  let q_start = blat_range_start(hsp.blocks, true)
  let q_end = blat_range_end(hsp.blocks, true)
  let h_start = blat_range_start(hsp.blocks, false)
  let h_end = blat_range_end(hsp.blocks, false)
  let line : Array[String] = []
  line.push(hsp.match_num.to_string())
  line.push(hsp.mismatch_num.to_string())
  line.push(hsp.match_rep_num.to_string())
  line.push(hsp.n_num.to_string())
  line.push(hsp.query_gapopen_num.to_string())
  line.push(hsp.query_gap_num.to_string())
  line.push(hsp.hit_gapopen_num.to_string())
  line.push(hsp.hit_gap_num.to_string())
  line.push(strand)
  line.push(qresult.id)
  line.push(qresult.seq_len.to_string())
  line.push(q_start.to_string())
  line.push(q_end.to_string())
  line.push(hit.id)
  line.push(hit.seq_len.to_string())
  line.push(h_start.to_string())
  line.push(h_end.to_string())
  line.push(hsp.blocks.length().to_string())
  line.push(blat_csv_join_ints(query_spans))
  line.push(blat_csv_join_ints(qstarts))
  line.push(blat_csv_join_ints(hstarts))
  if pslx {
    let qseqs : Array[String] = []
    let hseqs : Array[String] = []
    for block in hsp.blocks {
      qseqs.push(block.query_seq)
      hseqs.push(block.hit_seq)
    }
    line.push(blat_csv_join_strs(qseqs))
    line.push(blat_csv_join_strs(hseqs))
  }
  let sb = StringBuilder()
  let mut j = 0
  while j < line.length() {
    if j > 0 {
      sb.write_char('\t')
    }
    sb.write_string(line[j])
    j = j + 1
  }
  sb.to_string()
}

///|
/// Write QueryResults back to BLAT PSL/PSLX text, mirroring
/// `Bio.SearchIO.write(..., "blat-psl")`. A trailing newline is always
/// emitted; pass `header=true` for the psLayout version 3 preamble.
pub fn blat_psl_write(
  qresults : Array[BlatQueryResult],
  header? : Bool = false,
  pslx? : Bool = false,
) -> String raise SeqError {
  let sb = StringBuilder()
  if header {
    sb.write_string(blat_psl_header())
  }
  for qresult in qresults {
    if qresult.hits.length() > 0 {
      for hit in qresult.hits {
        for hsp in hit.hsps {
          sb.write_string(blat_build_row(qresult, hit, hsp, pslx))
          sb.write_char('\n')
        }
      }
    }
  }
  sb.to_string()
}

///|
/// Return (queries, hits, HSPs, fragments) that would be written.
pub fn blat_psl_write_stats(
  qresults : Array[BlatQueryResult],
) -> BlatWriteStats {
  let mut q_count = 0
  let mut hit_count = 0
  let mut hsp_count = 0
  let mut frag_count = 0
  for qresult in qresults {
    if qresult.hits.length() > 0 {
      q_count = q_count + 1
      for hit in qresult.hits {
        hit_count = hit_count + 1
        for hsp in hit.hsps {
          hsp_count = hsp_count + 1
          frag_count = frag_count + hsp.blocks.length()
        }
      }
    }
  }
  {
    qresults: q_count,
    hits: hit_count,
    hsps: hsp_count,
    fragments: frag_count,
  }
}

///|
/// Convert a BLAT query result into the generic SearchIO model
/// (`SearchIOQueryResult`). UCSC score maps to bitscore, ident_pct to
/// identity; coordinates are converted to 1-based inclusive convention used
/// by the generic model; BLAT provides no E-value, so evalue is 0.
pub fn BlatQueryResult::to_search_io(
  self : BlatQueryResult,
) -> SearchIOQueryResult {
  let mut qr = SearchIOQueryResult::new(
    id=self.id,
    description="",
    seq_length=self.seq_len,
    database="",
  )
  for hit in self.hits {
    let mut new_hit = SearchIOHit::new(
      id=hit.id,
      description="",
      seq_length=hit.seq_len,
    )
    for hsp in hit.hsps {
      let aligned = hsp.match_num +
        hsp.match_rep_num +
        hsp.mismatch_num +
        hsp.n_num +
        hsp.gap_num
      let gap_pct = if aligned == 0 {
        0.0
      } else {
        100.0 * hsp.gap_num.to_double() / aligned.to_double()
      }
      let qs = blat_range_start(hsp.blocks, true) + 1
      let qe = blat_range_end(hsp.blocks, true)
      let hs = blat_range_start(hsp.blocks, false) + 1
      let he = blat_range_end(hsp.blocks, false)
      let mut generic = SearchIOHsp::new(
        bitscore=hsp.score.to_double(),
        evalue=0.0,
        identity=hsp.ident_pct,
        positives=hsp.ident_pct,
        gap=gap_pct,
        alignment_length=aligned,
        query_start=qs,
        query_end=qe,
        hit_start=hs,
        hit_end=he,
      )
      let qseq = StringBuilder()
      let hseq = StringBuilder()
      let mut has_seq = false
      for block in hsp.blocks {
        if block.query_seq.length() > 0 {
          qseq.write_string(block.query_seq)
          has_seq = true
        }
        if block.hit_seq.length() > 0 {
          hseq.write_string(block.hit_seq)
          has_seq = true
        }
      }
      if has_seq {
        generic = generic.with_query_seq(qseq.to_string())
        generic = generic.with_hit_seq(hseq.to_string())
      }
      new_hit = new_hit.add_hsp(generic)
    }
    qr = qr.add_hit(new_hit)
  }
  qr
}
