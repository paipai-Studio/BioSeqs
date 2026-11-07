/// seq_filter_catalog.mbt — explainable sequence filter catalog
/// (FilterCatalog analogue, item B2)
///
/// A curated set of named rules (low complexity, poly-A/T tails, N runs,
/// adapter/vector contamination, restriction sites, stop codons) applied
/// in one pass. Every hit carries the rule id, category, severity,
/// genomic interval and a human-readable reason, so filtering decisions
/// are auditable.

///|
/// A named filter rule in the catalog.
pub struct SeqFilterRule {
  id : String
  category : String // "low_complexity" | "contamination" | "artifact" | "feature"
  severity : String // "warn" | "fail"
}

///|
/// A rule hit with its interval and explanation.
pub struct SeqFilterHit {
  rule_id : String
  category : String
  severity : String
  start : Int32
  end : Int32 // exclusive
  reason : String
}

///|
/// The default catalog: low complexity (DUST-style), homopolymer runs,
/// common Illumina/Nextera adapter prefixes, N runs, frequent restriction
/// sites and premature stop codon scan.
pub fn seq_filter_catalog_default() -> Array[SeqFilterRule] {
  [
    { id: "low_complexity_dust", category: "low_complexity", severity: "warn", },
    { id: "poly_a_tail", category: "artifact", severity: "warn", },
    { id: "poly_t_head", category: "artifact", severity: "warn", },
    { id: "n_run", category: "artifact", severity: "fail", },
    { id: "illumina_adapter", category: "contamination", severity: "fail", },
    { id: "nextera_adapter", category: "contamination", severity: "fail", },
    { id: "ecori_site", category: "feature", severity: "warn", },
    { id: "bamhi_site", category: "feature", severity: "warn", },
    { id: "pstI_site", category: "feature", severity: "warn", },
    { id: "stop_codon_scan", category: "feature", severity: "warn", },
  ]
}

///|
/// Count occurrences of a literal substring starting at or after `from`.
fn seq_filter_find_from(s : String, motif : String, from : Int32) -> Int32 {
  let n = s.length()
  let m = motif.length()
  if m == 0 {
    return -1
  }
  let mut i = from
  while i <= n - m {
    let mut ok = true
    for j = 0; j < m; j = j + 1 {
      if s.unsafe_get(i + j) != motif.unsafe_get(j) {
        ok = false
      }
    }
    if ok {
      return i
    }
    i = i + 1
  }
  -1
}

///|
fn seq_filter_homopolymer_start(
  s : String,
  base : Char,
  min_run : Int32,
) -> Int32 {
  let n = s.length()
  let mut run = 0
  let mut i = 0
  while i < n {
    if s.unsafe_get(i).unsafe_to_char() == base {
      run = run + 1
      if run >= min_run {
        return i - min_run + 1
      }
    } else {
      run = 0
    }
    i = i + 1
  }
  -1
}

///|
/// Apply the rules to one sequence and return all hits (sorted by start
/// position, stable).
pub fn seq_apply_filter_catalog(
  seq : Seq,
  rules~ : Array[SeqFilterRule] = seq_filter_catalog_default(),
) -> Array[SeqFilterHit] {
  let hits : Array[SeqFilterHit] = Array(capacity=0)
  let s = seq.to_string()
  let n = s.length()
  for rule in rules {
    match rule.id {
      "low_complexity_dust" => {
        // sliding 16bp windows, 8bp step, dust score threshold 0.6
        let w = 16
        let step = 8
        let mut i = 0
        while i + w <= n {
          let window = s[i:i + w].to_owned()
          let d = dust_score(window)
          if d > 0.6 {
            hits.push({
              rule_id: rule.id,
              category: rule.category,
              severity: rule.severity,
              start: i,
              end: i + w,
              reason: "dust score " + d.to_string() + " above 0.6",
            })
          }
          i = i + step
        }
      }
      "poly_a_tail" => {
        let st = seq_filter_homopolymer_start(s, 'A', 8)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 8,
            reason: "at least 8 consecutive A",
          })
        }
      }
      "poly_t_head" => {
        let st = seq_filter_homopolymer_start(s, 'T', 8)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 8,
            reason: "at least 8 consecutive T",
          })
        }
      }
      "n_run" => {
        let st = seq_filter_homopolymer_start(s, 'N', 5)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 5,
            reason: "at least 5 consecutive N",
          })
        }
      }
      "illumina_adapter" => {
        let st = seq_filter_find_from(s, "AGATCGGAAGAGC", 0)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 13,
            reason: "Illumina universal adapter prefix found",
          })
        }
      }
      "nextera_adapter" => {
        let st = seq_filter_find_from(s, "CTGTCTCTTATA", 0)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 12,
            reason: "Nextera transposase adapter prefix found",
          })
        }
      }
      "ecori_site" => {
        let st = seq_filter_find_from(s, "GAATTC", 0)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 6,
            reason: "EcoRI restriction site GAATTC",
          })
        }
      }
      "bamhi_site" => {
        let st = seq_filter_find_from(s, "GGATCC", 0)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 6,
            reason: "BamHI restriction site GGATCC",
          })
        }
      }
      "pstI_site" => {
        let st = seq_filter_find_from(s, "CTGCAG", 0)
        if st >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: st,
            end: st + 6,
            reason: "PstI restriction site CTGCAG",
          })
        }
      }
      "stop_codon_scan" => {
        // report the first in-frame-independent stop codon occurrence
        let mut best = -1
        for stop in ["TAA", "TAG", "TGA"] {
          let st = seq_filter_find_from(s, stop, 0)
          if st >= 0 && (best < 0 || st < best) {
            best = st
          }
        }
        if best >= 0 {
          hits.push({
            rule_id: rule.id,
            category: rule.category,
            severity: rule.severity,
            start: best,
            end: best + 3,
            reason: "stop codon (TAA/TAG/TGA) present",
          })
        }
      }
      _ => ()
    }
  }
  // stable sort by start
  for i = 1; i < hits.length(); i = i + 1 {
    let key = hits[i]
    let mut j = i - 1
    while j >= 0 && hits[j].start > key.start {
      hits[j + 1] = hits[j]
      j = j - 1
    }
    hits[j + 1] = key
  }
  hits
}

///|
/// Convenience predicate: does the sequence trigger any "fail"-severity
/// rule?
pub fn seq_filter_has_fatal(
  seq : Seq,
  rules~ : Array[SeqFilterRule] = seq_filter_catalog_default(),
) -> Bool {
  for hit in seq_apply_filter_catalog(seq, rules~) {
    if hit.severity == "fail" {
      return true
    }
  }
  false
}
