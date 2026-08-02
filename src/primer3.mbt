// Bio.Emboss.Primer3 - Primer design and analysis utilities.
//
// Inspired by Biopython's Bio.Emboss.Primer3 module, which wraps the
// EMBOSS primer3 tool for PCR primer design. This pure-MoonBit module
// provides the core primer-design primitives: melting temperature (Tm)
// via Wallace rule and a salt-adjusted formula, GC content, self- and
// cross-complementarity scoring, primer scoring, and a simple primer
// design pipeline that scans a template for primers satisfying
// user-supplied constraints. No external command-line tools are used.

// ============================================================================
// Primer3Config - design parameters
// ============================================================================

///|
/// Configuration parameters for primer design.
pub struct Primer3Config {
  primer_length_min : Int
  primer_length_max : Int
  tm_min : Double
  tm_max : Double
  gc_min : Double
  gc_max : Double
  max_self_complement : Int
  max_cross_complement : Int
  salt_conc : Double // mol/L (Na+ equivalent)
  flank : Int // minimal distance from template edges
}

///|
/// Build a `Primer3Config` with the given parameters.
pub fn Primer3Config::new(
  primer_length_min : Int,
  primer_length_max : Int,
  tm_min : Double,
  tm_max : Double,
  gc_min : Double,
  gc_max : Double,
  max_self_complement : Int,
  max_cross_complement : Int,
  salt_conc : Double,
  flank : Int
) -> Primer3Config {
  Primer3Config::{
    primer_length_min,
    primer_length_max,
    tm_min,
    tm_max,
    gc_min,
    gc_max,
    max_self_complement,
    max_cross_complement,
    salt_conc,
    flank
  }
}

///|
/// A sensible default configuration: 18-22 nt primers, Tm 55-65C,
/// GC 40-60%, salt 50 mM, 50 bp flanks.
pub fn primer3_default_config() -> Primer3Config {
  Primer3Config::new(18, 22, 55.0, 65.0, 0.4, 0.6, 6, 6, 0.05, 50)
}

///|
/// Minimum primer length.
pub fn Primer3Config::primer_length_min(self : Primer3Config) -> Int {
  self.primer_length_min
}

///|
/// Maximum primer length.
pub fn Primer3Config::primer_length_max(self : Primer3Config) -> Int {
  self.primer_length_max
}

///|
/// Minimum acceptable Tm in degrees Celsius.
pub fn Primer3Config::tm_min(self : Primer3Config) -> Double {
  self.tm_min
}

///|
/// Maximum acceptable Tm in degrees Celsius.
pub fn Primer3Config::tm_max(self : Primer3Config) -> Double {
  self.tm_max
}

///|
/// Minimum acceptable GC fraction (0..1).
pub fn Primer3Config::gc_min(self : Primer3Config) -> Double {
  self.gc_min
}

///|
/// Maximum acceptable GC fraction (0..1).
pub fn Primer3Config::gc_max(self : Primer3Config) -> Double {
  self.gc_max
}

///|
/// Maximum allowed self-complementarity score.
pub fn Primer3Config::max_self_complement(
  self : Primer3Config
) -> Int {
  self.max_self_complement
}

///|
/// Maximum allowed cross-complementarity score.
pub fn Primer3Config::max_cross_complement(
  self : Primer3Config
) -> Int {
  self.max_cross_complement
}

///|
/// Salt concentration (mol/L Na+ equivalent).
pub fn Primer3Config::salt_conc(self : Primer3Config) -> Double {
  self.salt_conc
}

///|
/// Minimum distance from either end of the template.
pub fn Primer3Config::flank(self : Primer3Config) -> Int {
  self.flank
}

// ============================================================================
// Primer3Result - a designed primer
// ============================================================================

///|
/// A primer candidate produced by `primer3_design_primers`.
pub struct Primer3Result {
  sequence : String
  start : Int // 0-based start in template
  length : Int
  tm : Double
  gc : Double
  self_complement : Int
  is_forward : Bool
  score : Double
}

///|
/// Construct a `Primer3Result`.
pub fn Primer3Result::new(
  sequence : String,
  start : Int,
  length : Int,
  tm : Double,
  gc : Double,
  self_complement : Int,
  is_forward : Bool,
  score : Double
) -> Primer3Result {
  Primer3Result::{
    sequence,
    start,
    length,
    tm,
    gc,
    self_complement,
    is_forward,
    score
  }
}

///|
/// Primer sequence (5' -> 3').
pub fn Primer3Result::sequence(self : Primer3Result) -> String {
  self.sequence
}

///|
/// 0-based start position in the template.
pub fn Primer3Result::start(self : Primer3Result) -> Int {
  self.start
}

///|
/// Primer length in nucleotides.
pub fn Primer3Result::length(self : Primer3Result) -> Int {
  self.length
}

///|
/// Melting temperature in degrees Celsius.
pub fn Primer3Result::tm(self : Primer3Result) -> Double {
  self.tm
}

///|
/// GC fraction (0..1).
pub fn Primer3Result::gc(self : Primer3Result) -> Double {
  self.gc
}

///|
/// Self-complementarity score (lower is better).
pub fn Primer3Result::self_complement(
  self : Primer3Result
) -> Int {
  self.self_complement
}

///|
/// True if forward primer, false if reverse.
pub fn Primer3Result::is_forward(self : Primer3Result) -> Bool {
  self.is_forward
}

///|
/// Overall primer quality score (higher is better).
pub fn Primer3Result::score(self : Primer3Result) -> Double {
  self.score
}

///|
/// Render a one-line summary of the primer.
pub fn Primer3Result::to_string(self : Primer3Result) -> String {
  let strand = if self.is_forward { "F" } else { "R" }
  let gc_pct = self.gc * 100.0
  "Primer[" +
  strand +
  "] len=" +
  self.length.to_string() +
  " start=" +
  self.start.to_string() +
  " Tm=" +
  self.tm.to_string() +
  " GC=" +
  gc_pct.to_string() +
  "% self=" +
  self.self_complement.to_string() +
  " score=" +
  self.score.to_string() +
  " seq=" +
  self.sequence
}

// ============================================================================
// Core calculations: GC content, Tm
// ============================================================================

///|
/// Compute the GC fraction (0..1) of a DNA sequence.
/// Non-ACGT characters are ignored.
pub fn primer3_gc_content(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  let up = seq.to_upper()
  let mut gc = 0
  let mut i = 0
  while i < n {
    let c = substring(up, i, 1)
    if c == "G" || c == "C" {
      gc = gc + 1
    }
    i = i + 1
  }
  gc.to_double() / n.to_double()
}

///|
/// Basic melting temperature using the Wallace rule:
///   Tm = 2 * (A+T) + 4 * (G+C)
/// Suitable for short oligos (< 14 nt).
pub fn primer3_tm_wallace(seq : String) -> Double {
  let up = seq.to_upper()
  let mut at = 0
  let mut gc = 0
  let mut i = 0
  while i < up.length() {
    let c = substring(up, i, 1)
    if c == "G" || c == "C" {
      gc = gc + 1
    } else if c == "A" || c == "T" {
      at = at + 1
    }
    i = i + 1
  }
  2.0 * at.to_double() + 4.0 * gc.to_double()
}

///|
/// Salt-adjusted melting temperature:
///   Tm = 81.5 + 16.6 * log10([Na+]) + 0.41 * %GC - 600 / length
/// Suitable for longer primers. `salt_conc` is in mol/L.
pub fn primer3_tm_salt(
  seq : String,
  salt_conc : Double
) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  let gc_frac = primer3_gc_content(seq)
  let pct_gc = gc_frac * 100.0
  let safe_salt = if salt_conc > 0.0 { salt_conc } else { 0.05 }
  let log10_salt = @math.ln(safe_salt) / @math.ln(10.0)
  81.5 + 16.6 * log10_salt + 0.41 * pct_gc - 600.0 / n.to_double()
}

///|
/// Choose an appropriate Tm formula based on length.
/// Short oligos (< 14 nt) use Wallace; longer use the salt-adjusted formula.
pub fn primer3_tm(
  seq : String,
  salt_conc : Double
) -> Double {
  if seq.length() < 14 {
    primer3_tm_wallace(seq)
  } else {
    primer3_tm_salt(seq, salt_conc)
  }
}

// ============================================================================
// Complementarity: hairpin and dimer scoring
// ============================================================================

///|
/// Complementarity base pairing: A-T, C-G.
fn primer3_pair(a : String, b : String) -> Bool {
  (a == "A" && b == "T") ||
  (a == "T" && b == "A") ||
  (a == "C" && b == "G") ||
  (a == "G" && b == "C")
}

///|
/// Score self-complementarity of a primer by sliding it against its own
/// reverse complement. Returns the maximum number of paired bases found
/// across all alignments. Lower scores are better.
pub fn primer3_self_complement(seq : String) -> Int {
  let n = seq.length()
  if n == 0 {
    return 0
  }
  let up = seq.to_upper()
  let rc = reverse_complement(up)
  let mut best = 0
  let mut offset = -(n - 1)
  while offset < n {
    let mut pairs = 0
    let mut i = 0
    while i < n {
      let j = i + offset
      if j >= 0 && j < n {
        let a = substring(up, i, 1)
        let b = substring(rc, j, 1)
        if primer3_pair(a, b) {
          pairs = pairs + 1
        }
      }
      i = i + 1
    }
    if pairs > best {
      best = pairs
    }
    offset = offset + 1
  }
  best
}

///|
/// Score cross-dimer complementarity between two primers by sliding
/// the reverse complement of seq2 against seq1 and counting paired bases.
/// Returns the maximum paired count across all offsets.
pub fn primer3_cross_dimer(
  seq1 : String,
  seq2 : String
) -> Int {
  let n1 = seq1.length()
  let n2 = seq2.length()
  if n1 == 0 || n2 == 0 {
    return 0
  }
  let up1 = seq1.to_upper()
  let rc2 = reverse_complement(seq2.to_upper())
  let m = rc2.length()
  let mut best = 0
  let mut offset = -(m - 1)
  while offset < n1 {
    let mut pairs = 0
    let mut i = 0
    while i < n1 {
      let j = i - offset
      if j >= 0 && j < m {
        let a = substring(up1, i, 1)
        let b = substring(rc2, j, 1)
        if primer3_pair(a, b) {
          pairs = pairs + 1
        }
      }
      i = i + 1
    }
    if pairs > best {
      best = pairs
    }
    offset = offset + 1
  }
  best
}

///|
/// Detect hairpin potential by checking self-complementarity at the 3' end.
/// Returns the self-complementarity score restricted to the last
/// `end_len` bases of the primer.
pub fn primer3_hairpin_3prime(
  seq : String,
  end_len : Int
) -> Int {
  let n = seq.length()
  if n == 0 || end_len <= 0 {
    return 0
  }
  let take = if end_len > n { n } else { end_len }
  let start = n - take
  let sub = substring(seq, start, take)
  primer3_self_complement(sub)
}

// ============================================================================
// Primer scoring
// ============================================================================

///|
/// Compute a quality score for a primer given a config.
/// Higher is better. The score rewards proximity to the Tm mid-point,
/// GC mid-point, and penalises self-complementarity.
pub fn primer3_score(
  seq : String,
  tm : Double,
  gc : Double,
  self_comp : Int,
  config : Primer3Config
) -> Double {
  let tm_mid = (config.tm_min + config.tm_max) / 2.0
  let tm_span = (config.tm_max - config.tm_min) / 2.0
  let tm_penalty = if tm_span > 0.0 {
    (tm - tm_mid).abs() / tm_span
  } else {
    0.0
  }
  let gc_mid = (config.gc_min + config.gc_max) / 2.0
  let gc_span = (config.gc_max - config.gc_min) / 2.0
  let gc_penalty = if gc_span > 0.0 {
    (gc - gc_mid).abs() / gc_span
  } else {
    0.0
  }
  let self_penalty = self_comp.to_double() / config.primer_length_min.to_double()
  // Score in [0, 100]; 100 is ideal.
  100.0 - 40.0 * tm_penalty - 30.0 * gc_penalty - 30.0 * self_penalty
}

// ============================================================================
// Primer design pipeline
// ============================================================================

///|
/// Container for primer-design results.
pub struct Primer3Design {
  forward : Array[Primer3Result]
  reverse : Array[Primer3Result]
}

///|
/// Construct a `Primer3Design`.
pub fn Primer3Design::new(
  forward : Array[Primer3Result],
  reverse : Array[Primer3Result]
) -> Primer3Design {
  Primer3Design::{ forward, reverse }
}

///|
/// Forward primer candidates sorted by score (best first).
pub fn Primer3Design::forward(self : Primer3Design) -> Array[Primer3Result] {
  self.forward
}

///|
/// Reverse primer candidates sorted by score (best first).
pub fn Primer3Design::reverse(self : Primer3Design) -> Array[Primer3Result] {
  self.reverse
}

///|
/// Number of forward candidates.
pub fn Primer3Design::n_forward(self : Primer3Design) -> Int {
  self.forward.length()
}

///|
/// Number of reverse candidates.
pub fn Primer3Design::n_reverse(self : Primer3Design) -> Int {
  self.reverse.length()
}

///|
/// True if at least one forward and one reverse candidate exist.
pub fn Primer3Design::has_pair(self : Primer3Design) -> Bool {
  self.forward.length() > 0 && self.reverse.length() > 0
}

///|
/// Best forward primer (highest score), or None if no candidate.
pub fn Primer3Design::best_forward(
  self : Primer3Design
) -> Primer3Result? {
  if self.forward.length() == 0 {
    None
  } else {
    Some(self.forward[0])
  }
}

///|
/// Best reverse primer (highest score), or None if no candidate.
pub fn Primer3Design::best_reverse(
  self : Primer3Design
) -> Primer3Result? {
  if self.reverse.length() == 0 {
    None
  } else {
    Some(self.reverse[0])
  }
}

///|
/// Best primer pair (forward, reverse), or None if either is missing.
pub fn Primer3Design::best_pair(
  self : Primer3Design
) -> (Primer3Result, Primer3Result)? {
  if self.has_pair() {
    Some((self.forward[0], self.reverse[0]))
  } else {
    None
  }
}

///|
/// Design forward and reverse primers for a template sequence satisfying
/// the constraints in `config`. Returns the best forward and reverse
/// primer candidates sorted by score (best first).
pub fn primer3_design_primers(
  template : String,
  config : Primer3Config
) -> Primer3Design {
  let n = template.length()
  let forward = Array::new()
  let reverse = Array::new()
  if n < config.primer_length_max + config.flank * 2 {
    return Primer3Design::new(forward, reverse)
  }
  // Forward primers: scan from left.
  let fwd_end_max = n - config.primer_length_max - config.flank
  let mut start = config.flank
  while start <= fwd_end_max {
    let mut len = config.primer_length_min
    while len <= config.primer_length_max {
      let seq = substring(template, start, len)
      let tm = primer3_tm(seq, config.salt_conc)
      let gc = primer3_gc_content(seq)
      let sc = primer3_self_complement(seq)
      if tm >= config.tm_min &&
        tm <= config.tm_max &&
        gc >= config.gc_min &&
        gc <= config.gc_max &&
        sc <= config.max_self_complement {
        let score = primer3_score(seq, tm, gc, sc, config)
        forward.push(
          Primer3Result::new(seq, start, len, tm, gc, sc, true, score)
        )
      }
      len = len + 1
    }
    start = start + 1
  }
  // Reverse primers: scan from right. Reverse primers are reverse
  // complements of the template's 3' region.
  let rev_start_min = config.primer_length_max + config.flank - 1
  let mut rstart = n - config.flank - 1
  while rstart >= rev_start_min {
    let mut len = config.primer_length_min
    while len <= config.primer_length_max {
      let ts = rstart - len + 1
      if ts < 0 {
        len = len + 1
        continue
      }
      let sub = substring(template, ts, len)
      let seq = reverse_complement(sub)
      let tm = primer3_tm(seq, config.salt_conc)
      let gc = primer3_gc_content(seq)
      let sc = primer3_self_complement(seq)
      if tm >= config.tm_min &&
        tm <= config.tm_max &&
        gc >= config.gc_min &&
        gc <= config.gc_max &&
        sc <= config.max_self_complement {
        let score = primer3_score(seq, tm, gc, sc, config)
        reverse.push(
          Primer3Result::new(seq, ts, len, tm, gc, sc, false, score)
        )
      }
      len = len + 1
    }
    rstart = rstart - 1
  }
  // Sort by score descending (highest/best first).
  forward.sort_by(fn(a, b) { if a.score() > b.score() { -1 } else if a.score() < b.score() { 1 } else { 0 } })
  reverse.sort_by(fn(a, b) { if a.score() > b.score() { -1 } else if a.score() < b.score() { 1 } else { 0 } })
  Primer3Design::new(forward, reverse)
}

// ============================================================================
// Amplicon utilities
// ============================================================================

///|
/// Compute the amplicon length given a forward primer start/length and
/// a reverse primer template start/length.
pub fn primer3_amplicon_length(
  fwd_start : Int,
  fwd_len : Int,
  rev_start : Int,
  rev_len : Int
) -> Int {
  let amp_end = rev_start + rev_len // exclusive end
  amp_end - fwd_start
}

///|
/// Predicted amplicon length for a primer pair.
pub fn primer3_pair_amplicon_length(
  fwd : Primer3Result,
  rev : Primer3Result
) -> Int {
  primer3_amplicon_length(fwd.start(), fwd.length(), rev.start(), rev.length())
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Build a sample DNA template (~600 bp) for primer design demos and
/// tests. The template is a mix of AT-rich and GC-rich regions so that
/// primer candidates with a range of Tm/GC values exist.
pub fn primer3_sample_template() -> String {
  // A 12-bp motif with ~50% GC, repeated to build a 600 bp template.
  // We use a mix of motifs to create variability.
  let motifs = [
    "ATGCGTACGTAG",
    "CTAGCATCGATC",
    "GATCGATCGATC",
    "ATCGATCGATCG",
    "CGTACGTACGTA",
    "TACGTACGTACG",
    "AGCTAGCTAGCT",
    "GCGCGATATATA",
    "ATATGCGCGCGC",
    "TCTAGATCTAGA",
    "GAGCTCGAGCTC",
    "AATTGGCCAATT",
    "TTAACCGGTTAA",
    "GGCCAATTGGCC",
    "AATTAATTGGCC",
    "CGGCCGAATTCC",
    "AAGCTTAAGCTT",
    "TTCGAATTCGAA",
    "GGATCCGGATCC",
    "AAGATCTAGATC",
    "CTAGAGATCTAG",
    "GATCTCTAGATC",
    "ATCGATCGATCG",
    "GCTAGCTAGCTA",
    "TACGTACGTACG",
    "CATGCATGCATG",
    "GTACGTACGTAC",
    "AGCTAGCTAGCT",
    "CGATCGATCGAT",
    "TAGCTAGCTAGC",
    "ATGcatgcATGC",
    "GCATgcATGCAT",
    "TTGGCCAATTGG",
    "CCAATTGGCCAA",
    "AATTGGCCAATT",
    "GGCCAATTGGCC",
    "AATTAATTGGCC",
    "CGGCCGAATTCC",
    "AAGCTTAAGCTT",
    "TTCGAATTCGAA",
    "GGATCCGGATCC",
    "AAGATCTAGATC",
    "CTAGAGATCTAG",
    "GATCTCTAGATC",
    "ATCGATCGATCG",
    "GCTAGCTAGCTA",
    "TACGTACGTACG",
    "CATGCATGCATG",
    "GTACGTACGTAC"
  ]
  let mut out = ""
  for m in motifs {
    out = out + m
  }
  out
}

///|
/// A short sample template (~200 bp) for quick tests.
pub fn primer3_short_template() -> String {
  let motifs = [
    "ATGCGTACGTAG",
    "CTAGCATCGATC",
    "GATCGATCGATC",
    "ATCGATCGATCG",
    "CGTACGTACGTA",
    "TACGTACGTACG",
    "AGCTAGCTAGCT",
    "GCGCGATATATA",
    "ATATGCGCGCGC",
    "TCTAGATCTAGA",
    "GAGCTCGAGCTC",
    "AATTGGCCAATT",
    "TTAACCGGTTAA",
    "GGCCAATTGGCC",
    "AATTAATTGGCC",
    "CGGCCGAATTCC"
  ]
  let mut out = ""
  for m in motifs {
    out = out + m
  }
  out
}

///|
/// Render a summary string for a `Primer3Design` showing the count of
/// forward/reverse candidates and the best primer in each direction.
pub fn primer3_design_summary(design : Primer3Design) -> String {
  let mut s = "Primer3Design Summary\n"
  s = s + "====================\n"
  s = s + "Forward candidates: \{design.n_forward()}\n"
  s = s + "Reverse candidates: \{design.n_reverse()}\n"
  s = s + "Has pair: \{design.has_pair()}\n"
  match design.best_forward() {
    Some(f) => s = s + "Best forward: \{f.to_string()}\n"
    None => s = s + "Best forward: (none)\n"
  }
  match design.best_reverse() {
    Some(r) => s = s + "Best reverse: \{r.to_string()}\n"
    None => s = s + "Best reverse: (none)\n"
  }
  match design.best_pair() {
    Some((f, r)) =>
      s = s +
        "Amplicon length: \{primer3_pair_amplicon_length(f, r)} bp\n"
    None => s = s + "Amplicon length: (no pair)\n"
  }
  s
}
