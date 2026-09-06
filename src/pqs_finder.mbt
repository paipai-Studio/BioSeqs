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
/// exhaustive/deep mode without score pruning). Every scoreable conformation
/// (score >= min_score) is also tallied into the per-position `density`
/// counter and `max_scores` vector, matching pqsfinder's deep-search
/// `density`/`maxScores` slots.
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
  density : Array[Int],
  max_scores : Array[Int],
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
                density,
                max_scores,
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
                let cstart = runs[0].0
                let clen = runs[3].1 - runs[0].0
                candidates.push(
                  {
                    score,
                    start: cstart,
                    len: clen,
                    features,
                  },
                )
                // tally every scoreable conformation on every position it
                // covers (deep-search density / maxScores semantics)
                for k = 0; k < clen; k = k + 1 {
                  let pos = cstart + k
                  density[pos] = density[pos] + 1
                  if score > max_scores[pos] {
                    max_scores[pos] = score
                  }
                }
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
              density,
              max_scores,
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
/// Export all scoreable conformations (port of `overlapping_storage`).
/// Candidates are enumerated in non-decreasing start order; for each start
/// group, conformations ending at the same position are deduplicated keeping
/// the highest score.
fn pqs_resolve_overlapping(
  candidates : Array[PqsCandidate],
) -> Array[PqsCandidate] {
  let result : Array[PqsCandidate] = []
  let mut group_start = -1
  let group_ends : Array[Int] = []
  let group_idx : Array[Int] = []
  for c in candidates {
    if c.start != group_start {
      group_start = c.start
      group_ends.clear()
      group_idx.clear()
    }
    let end = c.start + c.len
    let mut found = -1
    for gi = 0; gi < group_ends.length(); gi = gi + 1 {
      if group_ends[gi] == end {
        found = gi
        break
      }
    }
    if found < 0 {
      group_ends.push(end)
      group_idx.push(result.length())
      result.push(c)
    } else if c.score > result[group_idx[found]].score {
      result[group_idx[found]] = c
    }
  }
  // C++ exports each start group through a std::map keyed by end position,
  // i.e. the full output is ordered by (start, end).
  result.sort_by(fn(a : PqsCandidate, b : PqsCandidate) {
    if a.start == b.start {
      (a.start + a.len).compare(b.start + b.len)
    } else {
      a.start.compare(b.start)
    }
  })
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
/// Result of a PQS search: the reported (non-overlapping or all
/// overlapping) hits plus the deep-search per-position vectors.
///
/// `density[i]` is the number of scoreable PQS conformations covering
/// sequence position `i` (summed over both searched strands);
/// `max_scores[i]` is the highest score among all conformations covering
/// position `i`. These mirror pqsfinder's `density()` and `maxScores()`.
pub struct PqsResult {
  hits : Array[Pqs]
  density : Array[Int]
  max_scores : Array[Int]
}

///|
/// Extract the substring spanned by a PQS from the (original-strand)
/// sequence. Mirrors `as(pqs, "DNAStringSet")`: minus-strand hits return
/// the corresponding C-rich substring of the input.
pub fn Pqs::sequence(self : Pqs, seq : String) -> String {
  seq[self.start:self.end].to_string()
}

///|
/// Detect potential quadruplex-forming sequences and compute the
/// deep-search density/max-score vectors.
///
/// `strand` is "+", "-" or "*" (default, both strands). When
/// `overlapping` is true, every scoreable conformation is reported
/// (deduplicated per start/end position, best score kept), matching
/// pqsfinder's `overlapping = TRUE`. Antisense hits are mapped back to
/// the original coordinates with the run/loop features reversed.
pub fn pqs_search(
  seq : String,
  strand? : String = "*",
  overlapping? : Bool = false,
  options? : PqsOptions = PqsOptions::default(),
) -> PqsResult raise SeqError {
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
    return { hits: [], density: [], max_scores: [] }
  }
  let loop_penalties = pqs_loop_penalties(options)
  let bulge_penalties = pqs_run_penalties(options)
  let merged : Array[Pqs] = []
  let density : Array[Int] = Array::make(n, 0)
  let max_scores : Array[Int] = Array::make(n, 0)
  // Run the exhaustive search on one strand; `flip` indicates the
  // reverse-complement mapping back to original coordinates.
  let run_strand : (String, String, Bool) -> Unit = fn(subject, label, flip) {
    let candidates : Array[PqsCandidate] = []
    let sd : Array[Int] = Array::make(subject.length(), 0)
    let sm : Array[Int] = Array::make(subject.length(), 0)
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
      sd,
      sm,
    )
    let resolved = if overlapping {
      pqs_resolve_overlapping(candidates)
    } else {
      pqs_resolve_non_overlapping(candidates)
    }
    for c in resolved {
      if flip {
        let width = c.len
        let start = n - c.start - c.len
        let f = c.features
        merged.push(
          {
            start,
            end: start + width,
            strand: label,
            score: c.score,
            nt: f.nt,
            nb: f.nb,
            nm: f.nm,
            rl1: width - (f.rl1 + f.rl2 + f.rl3 + f.ll1 + f.ll2 + f.ll3),
            rl2: f.rl3,
            rl3: f.rl2,
            ll1: f.ll3,
            ll2: f.ll2,
            ll3: f.ll1,
          },
        )
      } else {
        merged.push(
          {
            start: c.start,
            end: c.start + c.len,
            strand: label,
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
    }
    // merge the per-position vectors; antisense positions are mirrored
    for j = 0; j < subject.length(); j = j + 1 {
      let pos = if flip { n - 1 - j } else { j }
      density[pos] = density[pos] + sd[j]
      if sm[j] > max_scores[pos] {
        max_scores[pos] = sm[j]
      }
    }
  }
  if strand == "+" || strand == "*" {
    run_strand(normalized, "+", false)
  }
  if strand == "-" || strand == "*" {
    run_strand(pqs_revcomp(normalized), "-", true)
  }
  merged.sort_by(fn(a : Pqs, b : Pqs) { a.start.compare(b.start) })
  { hits: merged, density, max_scores }
}

///|
/// Detect potential quadruplex-forming sequences.
///
/// `strand` is "+", "-" or "*" (default, both strands). When
/// `overlapping` is true, all scoreable conformations are reported
/// instead of the non-overlapping selection. Antisense hits are mapped
/// back to the original coordinates with the run/loop features reversed,
/// matching pqsfinder output. Use `pqs_search` to also obtain the
/// `density` and `max_scores` vectors.
pub fn pqs_finder(
  seq : String,
  strand? : String = "*",
  overlapping? : Bool = false,
  options? : PqsOptions = PqsOptions::default(),
) -> Array[Pqs] raise SeqError {
  pqs_search(seq, strand=strand, overlapping=overlapping, options=options).hits
}
