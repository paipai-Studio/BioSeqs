// G-quadruplex (PQS) detection - port of Bioconductor `pqsfinder`.
//
// Identifies DNA sequence patterns likely to fold into an intramolecular
// G-quadruplex. The algorithm scans for four consecutive G-runs separated by
// loops, allowing bulges and mismatches inside the runs, and assigns each
// candidate a quantitative stability score. Reported PQS (potential
// quadruplex-forming sequences) are non-overlapping and sorted by start.
//
// Reference: Hon J, Martinek T, Lexa M (2020), Bioinformatics 36:5494-5500.
// This port implements the exhaustive (deep) search semantics of
// pqsfinder 2.x with the default scoring system.

///|
/// Scoring and search options, mirroring pqsfinder defaults.
pub struct PqsOptions {
  max_len : Int
  min_score : Int
  run_min_len : Int
  run_max_len : Int
  loop_min_len : Int
  loop_max_len : Int
  max_bulges : Int
  max_mismatches : Int
  max_defects : Int
  tetrad_bonus : Int
  mismatch_penalty : Int
  bulge_penalty : Int
  bulge_len_factor : Double
  bulge_len_exponent : Double
  loop_mean_factor : Double
  loop_mean_exponent : Double
}

///|
/// Default options identical to `pqsfinder()` defaults.
pub fn PqsOptions::default() -> PqsOptions {
  {
    max_len: 50,
    min_score: 47,
    run_min_len: 2,
    run_max_len: 11,
    loop_min_len: 0,
    loop_max_len: 30,
    max_bulges: 3,
    max_mismatches: 3,
    max_defects: 3,
    tetrad_bonus: 40,
    mismatch_penalty: 28,
    bulge_penalty: 20,
    bulge_len_factor: 0.2,
    bulge_len_exponent: 1.0,
    loop_mean_factor: 6.6,
    loop_mean_exponent: 0.8,
  }
}

///|
/// Construct options with explicit overrides; unset parameters keep the
/// pqsfinder defaults.
pub fn PqsOptions::of(
  max_len? : Int = 50,
  min_score? : Int = 47,
  run_min_len? : Int = 2,
  run_max_len? : Int = 11,
  loop_min_len? : Int = 0,
  loop_max_len? : Int = 30,
  max_bulges? : Int = 3,
  max_mismatches? : Int = 3,
  max_defects? : Int = 3,
  tetrad_bonus? : Int = 40,
  mismatch_penalty? : Int = 28,
  bulge_penalty? : Int = 20,
  bulge_len_factor? : Double = 0.2,
  bulge_len_exponent? : Double = 1.0,
  loop_mean_factor? : Double = 6.6,
  loop_mean_exponent? : Double = 0.8,
) -> PqsOptions {
  {
    max_len,
    min_score,
    run_min_len,
    run_max_len,
    loop_min_len,
    loop_max_len,
    max_bulges,
    max_mismatches,
    max_defects,
    tetrad_bonus,
    mismatch_penalty,
    bulge_penalty,
    bulge_len_factor,
    bulge_len_exponent,
    loop_mean_factor,
    loop_mean_exponent,
  }
}

///|
/// One detected potential quadruplex-forming sequence.
///
/// `start` is 0-based inclusive and `end` is 0-based exclusive.
/// Features: `nt` tetrads, `nb` bulges, `nm` mismatches, `rl1..rl3` run
/// lengths and `ll1..ll3` loop lengths.
pub struct Pqs {
  start : Int
  end : Int
  strand : String
  score : Int
  nt : Int
  nb : Int
  nm : Int
  rl1 : Int
  rl2 : Int
  rl3 : Int
  ll1 : Int
  ll2 : Int
  ll3 : Int
}

///|
priv struct PqsFeatures {
  mut nt : Int
  mut nb : Int
  mut nm : Int
  mut rl1 : Int
  mut rl2 : Int
  mut rl3 : Int
  mut ll1 : Int
  mut ll2 : Int
  mut ll3 : Int
}

///|
fn PqsFeatures::zero() -> PqsFeatures {
  { nt: 0, nb: 0, nm: 0, rl1: 0, rl2: 0, rl3: 0, ll1: 0, ll2: 0, ll3: 0 }
}

///|
/// Round half away from zero (C `round` semantics used by pqsfinder tables).
fn pqs_round(x : Double) -> Int {
  if x >= 0.0 {
    (x + 0.5).floor().to_int()
  } else {
    -((-x + 0.5).floor().to_int())
  }
}

///|
fn pqs_run_penalties(sc : PqsOptions) -> Array[Int] {
  let out : Array[Int] = []
  let mut i = 0
  while i < sc.run_max_len - 1 {
    out.push(
      pqs_round(
        sc.bulge_len_factor * @math.pow(i.to_double(), sc.bulge_len_exponent),
      ),
    )
    i = i + 1
  }
  out
}

///|
fn pqs_loop_penalties(sc : PqsOptions) -> Array[Int] {
  let out : Array[Int] = []
  let len_sum_max = sc.loop_max_len * 3
  let mut i = 0
  while i <= len_sum_max {
    let mean = i.to_double() / 3.0
    out.push(
      pqs_round(sc.loop_mean_factor * @math.pow(mean, sc.loop_mean_exponent)),
    )
    i = i + 1
  }
  out
}

///|
/// Count leading + trailing G's of a run window (pqsfinder `count_g_num`).
fn pqs_count_g(seq : String, first : Int, second : Int) -> Int {
  let mut cnt = 0
  let mut s = first
  let mut e = second - 1
  while s < second && seq.unsafe_get(s) == 'G' {
    s = s + 1
    cnt = cnt + 1
  }
  while e > s && seq.unsafe_get(e) == 'G' {
    e = e - 1
    cnt = cnt + 1
  }
  cnt
}

///|
/// Find a G-run inside [s, limit): the window starts at the first G and ends
/// at the last G within `run_max_len`, with length >= `run_min_len`.
/// Port of pqsfinder's non-regex `find_run`.
fn pqs_find_run(
  seq : String,
  start : Int,
  limit : Int,
  run_min_len : Int,
  run_max_len : Int,
) -> (Int, Int)? {
  let mut s = start
  while s < limit {
    while s < limit && seq.unsafe_get(s) != 'G' {
      s = s + 1
    }
    if s >= limit {
      break
    }
    let cap = if s + run_max_len < limit { s + run_max_len } else { limit }
    let mut e = cap - 1
    while e > s && seq.unsafe_get(e) != 'G' {
      e = e - 1
    }
    if e - s + 1 >= run_min_len {
      return Some((s, e + 1))
    }
    s = s + 1
  }
  None
}

///|
priv struct PqsCandidate {
  score : Int
  start : Int
  len : Int
  features : PqsFeatures
}

///|
/// Score four runs forming a quadruplex (port of `score_pqs` +
/// `score_run_defects`).
fn pqs_score_pqs(
  seq : String,
  m : Array[(Int, Int)],
  sc : PqsOptions,
  loop_penalties : Array[Int],
  bulge_penalties : Array[Int],
  features : PqsFeatures,
) -> Int {
  let run_cnt = 4
  let l : Array[Int] = []
  let w : Array[Int] = []
  let g : Array[Int] = []
  let mut i = 0
  while i < run_cnt {
    w.push(m[i].1 - m[i].0)
    g.push(pqs_count_g(seq, m[i].0, m[i].1))
    if i < run_cnt - 1 {
      l.push(m[i + 1].0 - m[i].1)
    }
    i = i + 1
  }
  // at most one zero-length loop is scoreable
  if sc.loop_min_len == 0 &&
    ((l[0] == 0 && l[1] == 0) || (l[0] == 0 && l[2] == 0) || (l[1] == 0 && l[2] == 0)) {
    return 0
  }
  // reference run: minimal width among pure-G runs
  let mut min_pw = 2147483647
  let mut pi = -1
  let mut i2 = 0
  while i2 < run_cnt {
    if w[i2] < min_pw && g[i2] == w[i2] {
      min_pw = w[i2]
      pi = i2
    }
    i2 = i2 + 1
  }
  if pi < 0 {
    return 0
  }
  // score run defects relative to the reference run
  let mut mismatches = 0
  let mut bulges = 0
  let mut score = 0
  let mut i3 = 0
  while i3 < run_cnt {
    if w[i3] == w[pi] && g[i3] == g[pi] {
      // perfect run
    } else if w[i3] == w[pi] && g[i3] == g[pi] - 1 {
      mismatches = mismatches + 1
    } else if w[i3] > w[pi] && g[i3] >= g[pi] {
      bulges = bulges + 1
      score = score - bulge_penalties[w[i3] - w[pi]]
    } else {
      return 0
    }
    i3 = i3 + 1
  }
  if mismatches > sc.max_mismatches ||
    bulges > sc.max_bulges ||
    mismatches + bulges > sc.max_defects {
    return 0
  }
  score = score +
    (w[pi] - 1) * sc.tetrad_bonus -
    mismatches * sc.mismatch_penalty -
    bulges * sc.bulge_penalty
  if score <= 0 {
    return 0
  }
  features.nt = w[pi]
  features.nb = bulges
  features.nm = mismatches
  features.rl1 = w[0]
  features.rl2 = w[1]
  features.rl3 = w[2]
  features.ll1 = l[0]
  features.ll2 = l[1]
  features.ll3 = l[2]
  let final = score - loop_penalties[l[0] + l[1] + l[2]]
  if final > 0 {
    final
  } else {
    0
  }
}

///|
/// Recursive search for four consecutive runs (port of `find_all_runs`,
/// exhaustive mode without score pruning).
fn pqs_find_all_runs(
  seq : String,
  sc : PqsOptions,
  loop_penalties : Array[Int],
  bulge_penalties : Array[Int],
  runs : Array[(Int, Int)],
  i : Int,
  start : Int,
  limit : Int,
  zero_loop : Bool,
  candidates : Array[PqsCandidate],
) -> Unit {
  let mut start = start
  if i > 0 {
    let prev_loop_len = start - runs[i - 1].1
    if prev_loop_len < sc.loop_min_len {
      let moved = runs[i - 1].1 + sc.loop_min_len
      start = if moved < limit { moved } else { limit }
    } else if prev_loop_len == 0 && zero_loop {
      let moved = runs[i - 1].1 + 1
      start = if moved < limit { moved } else { limit }
    }
  }
  let mut s = start
  while s < limit {
    let min_e = s + sc.run_min_len
    let mut found_any = false
    let mut e = limit
    let mut continue_inner = true
    while continue_inner && e >= min_e {
      match pqs_find_run(seq, s, e, sc.run_min_len, sc.run_max_len) {
        None => {
          continue_inner = false
        }
        Some(r) => {
          runs[i] = r
          found_any = true
          s = r.0
          e = r.1
          if i > 0 {
            let loop_len = s - runs[i - 1].1
            if loop_len > sc.loop_max_len {
              return
            }
            let next_zero = if loop_len == 0 { true } else { zero_loop }
            if i < 3 {
              pqs_find_all_runs(
                seq,
                sc,
                loop_penalties,
                bulge_penalties,
                runs,
                i + 1,
                r.1,
                limit,
                next_zero,
                candidates,
              )
            } else {
              let features = PqsFeatures::zero()
              let score = pqs_score_pqs(
                seq,
                runs,
                sc,
                loop_penalties,
                bulge_penalties,
                features,
              )
              if score >= sc.min_score {
                candidates.push(
                  {
                    score,
                    start: runs[0].0,
                    len: runs[3].1 - runs[0].0,
                    features,
                  },
                )
              }
            }
          } else {
            // enforce G4 total length limit relative to the first run start
            let bound = s + sc.max_len
            pqs_find_all_runs(
              seq,
              sc,
              loop_penalties,
              bulge_penalties,
              runs,
              1,
              r.1,
              if bound < limit { bound } else { limit },
              false,
              candidates,
            )
          }
          e = e - 1
        }
      }
    }
    if !found_any {
      break
    }
    s = s + 1
  }
}

///|
/// Resolve overlapping candidates preferring higher score, then longer hit.
/// Port of `revised_non_overlapping_storage::export_pqs`.
fn pqs_resolve_non_overlapping(
  candidates : Array[PqsCandidate],
) -> Array[PqsCandidate] {
  if candidates.is_empty() {
    return []
  }
  // collect distinct candidate scores
  let scores : Array[Int] = []
  for c in candidates {
    let mut found = false
    for s in scores {
      if s == c.score {
        found = true
        break
      }
    }
    if !found {
      scores.push(c.score)
    }
  }
  scores.sort_by(fn(a : Int, b : Int) { b.compare(a) })
  let result : Array[PqsCandidate] = []
  for score in scores {
    let level : Array[PqsCandidate] = []
    for c in candidates {
      if c.score == score {
        level.push(c)
      }
    }
    // resolve overlaps within the same score level, keeping the longer hit
    let kept : Array[PqsCandidate] = []
    for c in level {
      match kept.last() {
        None => kept.push(c)
        Some(prev) => {
          if c.start < prev.start + prev.len {
            if c.len < prev.len {
              // previous is longer or equal: drop current
            } else {
              kept.pop() |> ignore
              kept.push(c)
            }
          } else {
            kept.push(c)
          }
        }
      }
    }
    // save hits that do not overlap an already saved, higher scoring hit
    for c in kept {
      let mut overlaps = false
      for r in result {
        if c.start < r.start + r.len && r.start < c.start + c.len {
          overlaps = true
          break
        }
      }
      if !overlaps {
        result.push(c)
      }
    }
  }
  result
}

///|
fn pqs_revcomp(seq : String) -> String {
  let n = seq.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let c = seq.unsafe_get(i)
    let comp : Char = match c {
      'A' => 'T'
      'T' => 'A'
      'U' => 'A'
      'G' => 'C'
      'C' => 'G'
      _ => 'N'
    }
    buf[n - 1 - i] = comp.to_int().to_uint16()
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Detect potential quadruplex-forming sequences.
///
/// `strand` is "+", "-" or "*" (default, both strands). Antisense hits are
/// mapped back to the original coordinates with the run/loop features
/// reversed, matching pqsfinder output.
pub fn pqs_finder(
  seq : String,
  strand? : String = "*",
  options? : PqsOptions = PqsOptions::default(),
) -> Array[Pqs] raise SeqError {
  if strand != "+" && strand != "-" && strand != "*" {
    raise SeqError("strand must be \"+\", \"-\" or \"*\"")
  }
  if options.run_min_len < 1 {
    raise SeqError("run_min_len must be at least 1")
  }
  if options.run_max_len < options.run_min_len {
    raise SeqError("run_max_len must be >= run_min_len")
  }
  let normalized = seq.to_upper()
  let n = normalized.length()
  if n == 0 {
    return []
  }
  let loop_penalties = pqs_loop_penalties(options)
  let bulge_penalties = pqs_run_penalties(options)
  let search : (String, String) -> Array[Pqs] = fn(subject, strand_label) {
    let candidates : Array[PqsCandidate] = []
    let runs : Array[(Int, Int)] = [(0, 0), (0, 0), (0, 0), (0, 0)]
    pqs_find_all_runs(
      subject,
      options,
      loop_penalties,
      bulge_penalties,
      runs,
      0,
      0,
      subject.length(),
      false,
      candidates,
    )
    let resolved = pqs_resolve_non_overlapping(candidates)
    let items : Array[Pqs] = []
    for c in resolved {
      items.push(
        {
          start: c.start,
          end: c.start + c.len,
          strand: strand_label,
          score: c.score,
          nt: c.features.nt,
          nb: c.features.nb,
          nm: c.features.nm,
          rl1: c.features.rl1,
          rl2: c.features.rl2,
          rl3: c.features.rl3,
          ll1: c.features.ll1,
          ll2: c.features.ll2,
          ll3: c.features.ll3,
        },
      )
    }
    items
  }
  let merged : Array[Pqs] = []
  if strand == "+" || strand == "*" {
    for p in search(normalized, "+") {
      merged.push(p)
    }
  }
  if strand == "-" || strand == "*" {
    let rc = pqs_revcomp(normalized)
    for p in search(rc, "-") {
      // map coordinates back to the original sequence and reverse features
      let width = p.end - p.start
      let start = n - p.end
      merged.push(
        {
          start,
          end: start + width,
          strand: "-",
          score: p.score,
          nt: p.nt,
          nb: p.nb,
          nm: p.nm,
          rl1: width - (p.rl1 + p.rl2 + p.rl3 + p.ll1 + p.ll2 + p.ll3),
          rl2: p.rl3,
          rl3: p.rl2,
          ll1: p.ll3,
          ll2: p.ll2,
          ll3: p.ll1,
        },
      )
    }
  }
  merged.sort_by(fn(a : Pqs, b : Pqs) { a.start.compare(b.start) })
  merged
}
