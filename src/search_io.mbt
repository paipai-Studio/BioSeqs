///|
/// Bio.SearchIO module, inspired by Biopython Bio.SearchIO.

///|
pub struct HSPFragment {
  query_start : Int
  query_end : Int
  hit_start : Int
  hit_end : Int
  query_strand : Strand
  hit_strand : Strand
  aln_span : Int
  gap_num : Int
  bitscore : Double
  evalue : Double
  query_seq : String
  hit_seq : String
  match_seq : String
}

///|
pub fn HSPFragment::new(
  query_start : Int,
  query_end : Int,
  hit_start : Int,
  hit_end : Int,
  query_strand : Strand,
  hit_strand : Strand,
  aln_span : Int,
  gap_num : Int,
  bitscore : Double,
  evalue : Double,
  query_seq : String,
  hit_seq : String,
  match_seq : String,
) -> HSPFragment {
  HSPFragment::{
    query_start,
    query_end,
    hit_start,
    hit_end,
    query_strand,
    hit_strand,
    aln_span,
    gap_num,
    bitscore,
    evalue,
    query_seq,
    hit_seq,
    match_seq,
  }
}

///|
pub struct HSP {
  hit_id : String
  hit_def : String
  fragments : Array[HSPFragment]
  bitscore : Double
  evalue : Double
  score : Double
}

///|
pub fn HSP::new(
  hit_id : String,
  hit_def : String,
  fragments : Array[HSPFragment],
  bitscore : Double,
  evalue : Double,
  score : Double,
) -> HSP {
  HSP::{ hit_id, hit_def, fragments, bitscore, evalue, score }
}

///|
pub struct Hit {
  id : String
  description : String
  accession : String
  hsps : Array[HSP]
  bitscore : Double
  evalue : Double
}

///|
pub fn Hit::new(
  id : String,
  description : String,
  accession : String,
  hsps : Array[HSP],
  bitscore : Double,
  evalue : Double,
) -> Hit {
  Hit::{ id, description, accession, hsps, bitscore, evalue }
}

///|
pub struct QueryResult {
  id : String
  description : String
  seq_len : Int
  hits : Array[Hit]
  program : String
  version : String
  target : String
}

///|
pub fn QueryResult::new(
  id : String,
  description : String,
  seq_len : Int,
  hits : Array[Hit],
  program : String,
  version : String,
  target : String,
) -> QueryResult {
  QueryResult::{ id, description, seq_len, hits, program, version, target }
}

///|
pub fn parse_hmmer3_tab(content : String) -> Array[QueryResult] {
  let lines = split_by_char(content, '\n'.to_int())

  let results : Array[QueryResult] = Array::new()
  let mut cq = ""
  let mut hc = false
  let mut ch : Array[Hit] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let ln = lines[i]
    if ln.length() == 0 || starts_with(ln, 0, "#") {
      i = i + 1
      continue
    }

    let parts = split_by_whitespace(ln)
    if parts.length() < 6 {
      i = i + 1
      continue
    }

    let qn = parts[0]
    let tn = parts[1]
    let ta = parts[2]
    let ev = parse_double(parts[3]).unwrap_or(0.0)
    let sc = parse_double(parts[4]).unwrap_or(0.0)
    let bs = sc
    let mut desc = ""
    let mut j = 5
    while j < parts.length() {
      desc = desc + parts[j] + " "
      j = j + 1
    }
    desc = desc

    if qn != cq {
      if hc {
        results.push(QueryResult::new(cq, "", 0, ch, "hmmscan", "3.0", ta))
      }
      cq = qn
      ch = Array::new()
      hc = true
    }

    let frags : Array[HSPFragment] = Array::new()
    frags.push(
      HSPFragment::new(
        0,
        0,
        0,
        0,
        Strand::Plus,
        Strand::Plus,
        0,
        0,
        bs,
        ev,
        "",
        "",
        "",
      ),
    )
    let hsps_arr : Array[HSP] = Array::new()
    hsps_arr.push(HSP::new(tn, desc, frags, bs, ev, sc))
    ch.push(Hit::new(tn, desc, ta, hsps_arr, bs, ev))

    i = i + 1
  }

  if hc {
    results.push(QueryResult::new(cq, "", 0, ch, "hmmscan", "3.0", ""))
  }

  results
}

///|
pub fn parse_blat_psl(content : String) -> Array[QueryResult] {
  let lines = split_by_char(content, '\n'.to_int())

  let results : Array[QueryResult] = Array::new()
  let mut cq = ""
  let mut hc = false
  let mut ch : Array[Hit] = Array::new()
  let mut cqs = 0

  let mut i = 0
  while i < lines.length() {
    let ln = lines[i]
    if ln.length() == 0 || starts_with(ln, 0, "#") {
      i = i + 1
      continue
    }

    let parts = split_by_char(ln, '\t'.to_int())
    if parts.length() < 20 {
      i = i + 1
      continue
    }

    let strand = parts[8]
    let qs = if starts_with(strand, 0, "-") {
      Strand::Minus
    } else {
      Strand::Plus
    }
    let ts = if starts_with(strand, strand.length() - 1, "-") {
      Strand::Minus
    } else {
      Strand::Plus
    }

    let qn = parts[9]
    let qsz = parse_int(parts[10])
    let tn = parts[13]

    let bc = parse_int(parts[17])

    let bsizes = split_by_char(parts[18], ','.to_int())
    let qstarts = split_by_char(parts[19], ','.to_int())
    let tstarts = split_by_char(parts[20], ','.to_int())

    let matches = parse_int(parts[0])
    let score = matches.to_double()

    if qn != cq {
      if hc {
        results.push(QueryResult::new(cq, "", cqs, ch, "blat", "35", tn))
      }
      cq = qn
      cqs = qsz
      ch = Array::new()
      hc = true
    }

    let frags : Array[HSPFragment] = Array::new()

    let mut j = 0
    while j < bc {
      let bsz = if j < bsizes.length() { parse_int(bsizes[j]) } else { 0 }
      let qst = if j < qstarts.length() { parse_int(qstarts[j]) } else { 0 }
      let tst = if j < tstarts.length() { parse_int(tstarts[j]) } else { 0 }

      frags.push(
        HSPFragment::new(
          qst,
          qst + bsz,
          tst,
          tst + bsz,
          qs,
          ts,
          bsz,
          0,
          score,
          0.0,
          "",
          "",
          "",
        ),
      )
      j = j + 1
    }

    let hsps_arr : Array[HSP] = Array::new()
    hsps_arr.push(HSP::new(tn, "", frags, score, 0.0, score))
    ch.push(Hit::new(tn, "", "", hsps_arr, score, 0.0))

    i = i + 1
  }

  if hc {
    results.push(QueryResult::new(cq, "", cqs, ch, "blat", "35", ""))
  }

  results
}

///|
pub fn blast_to_searchio(
  blast_records : Array[BlastRecord],
) -> Array[QueryResult] {
  let results : Array[QueryResult] = Array::new()

  let mut i = 0
  while i < blast_records.length() {
    let record = blast_records[i]
    let hits : Array[Hit] = Array::new()

    let mut j = 0
    while j < record.hits.length() {
      let bh = record.hits[j]
      let hsps_arr : Array[HSP] = Array::new()

      let mut k = 0
      while k < bh.hsps.length() {
        let bhsp = bh.hsps[k]

        let frags : Array[HSPFragment] = Array::new()
        frags.push(
          HSPFragment::new(
            bhsp.query_start(),
            bhsp.query_end(),
            bhsp.subject_start(),
            bhsp.subject_end(),
            bhsp.query_strand(),
            bhsp.subject_strand(),
            bhsp.align_length(),
            bhsp.gaps(),
            bhsp.bitscore(),
            bhsp.evalue(),
            bhsp.query_seq(),
            bhsp.subject_seq(),
            bhsp.match_seq(),
          ),
        )

        hsps_arr.push(
          HSP::new(
            bh.subject_id(),
            bh.subject_def(),
            frags,
            bhsp.bitscore(),
            bhsp.evalue(),
            bhsp.bitscore(),
          ),
        )
        k = k + 1
      }

      let mut tbs = 0.0
      let mut mev = 1.0
      let mut l = 0
      while l < hsps_arr.length() {
        tbs = tbs + hsps_arr[l].bitscore
        if hsps_arr[l].evalue < mev {
          mev = hsps_arr[l].evalue
        }
        l = l + 1
      }

      hits.push(
        Hit::new(bh.subject_id(), bh.subject_def(), "", hsps_arr, tbs, mev),
      )
      j = j + 1
    }

    results.push(
      QueryResult::new(
        record.query_id(),
        record.query_def(),
        record.query_length(),
        hits,
        "blast",
        "2.14.0",
        record.database(),
      ),
    )
    i = i + 1
  }

  results
}

///|
pub fn top_hits(qr : QueryResult, n : Int, threshold : Double) -> Array[Hit] {
  let filtered : Array[Hit] = Array::new()

  let mut i = 0
  while i < qr.hits.length() {
    if qr.hits[i].evalue < threshold {
      filtered.push(qr.hits[i])
    }
    i = i + 1
  }

  let cnt = filtered.length()
  let sorted : Array[Hit] = Array::new()
  let mut j = 0
  while j < cnt {
    sorted.push(filtered[j])
    j = j + 1
  }

  let mut k = 0
  while k < cnt - 1 {
    let mut l = k + 1
    while l < cnt {
      if sorted[k].evalue > sorted[l].evalue {
        let temp = sorted[k]
        sorted[k] = sorted[l]
        sorted[l] = temp
      }
      l = l + 1
    }
    k = k + 1
  }

  let take = if n < cnt { n } else { cnt }
  let result : Array[Hit] = Array::new()
  let mut m = 0
  while m < take {
    result.push(sorted[m])
    m = m + 1
  }
  result
}

///|
pub fn count_hsps(qr : QueryResult) -> Int {
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
