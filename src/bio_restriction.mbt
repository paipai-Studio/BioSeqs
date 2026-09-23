///|
/// Bio.Restriction — Biopython 1.88 compatible restriction enzyme analysis.
///
/// 1088 REBASE enzymes (see `restriction_db.mbt`), IUPAC ambiguous sites,
/// linear and circular search, one-cut / two-cut / uncharacterised enzymes,
/// catalyse fragments, elucidate annotations, compatible ends and
/// isoschizomer/neoschizomer relationships, RestrictionBatch search and
/// RestrictionAnalysis filters.
///
/// Cut positions follow the Biopython convention: 1-based index of the first
/// base of the 3' fragment, e.g. EcoRI.search("GAATTC") == [2].

///|
/// A restriction enzyme record, faithful to Biopython's Bio.Restriction.
///
/// Enzymes are value records, so equality is structural: two records are equal
/// exactly when all fields (including the name) are equal.
pub struct BioRestrictionEnzyme {
  name : String
  site : String
  size : Int
  fst5 : Int?
  fst3 : Int?
  scd5 : Int?
  scd3 : Int?
  /// Signed overhang length: negative = 5' overhang, positive = 3' overhang,
  /// 0 = blunt, None = unknown.
  ovhg : Int?
  ovhgseq : String
  /// Theoretical frequency: one cut per x bases.
  freq : Double
  supplier_codes : String
  palindromic : Bool
  /// 0 = defined, 1 = ambiguous overhang, 2 = not defined (uncharacterised).
  spec : Int
  elucidation : String
} derive(Eq)
pub extend BioRestrictionEnzyme with Eq::{equal, not_equal}

///|
/// All 1088 enzymes as records.
pub let bio_restriction_enzyme_list : Array[BioRestrictionEnzyme] = restriction_db.map(fn(
    r,
  ) {
    BioRestrictionEnzyme::{
      name: r.0,
      site: r.1,
      size: r.2,
      fst5: r.3,
      fst3: r.4,
      scd5: r.5,
      scd3: r.6,
      ovhg: r.7,
      ovhgseq: r.8,
      freq: r.9,
      supplier_codes: r.10,
      palindromic: r.11,
      spec: r.12,
      elucidation: r.13,
    }
  },
)

///|
/// Name -> enzyme record.
pub let bio_restriction_enzyme_map : Map[String, BioRestrictionEnzyme] = {
  let m : Map[String, BioRestrictionEnzyme] = Map(
    [],
    capacity=restriction_db_size,
  )
  for e in bio_restriction_enzyme_list {
    m[e.name] = e
  }
  m
}

///|
/// Number of characterised / ambiguous / uncharacterised enzymes.
pub fn bio_restriction_enzyme_count() -> Int {
  restriction_db_size
}

///|
/// Look up an enzyme by name (e.g. "EcoRI").
pub fn get_bio_restriction_enzyme(name : String) -> BioRestrictionEnzyme? {
  bio_restriction_enzyme_map.get(name)
}

///|
/// Look up an enzyme by name or raise `SeqError`.
pub fn require_bio_restriction_enzyme(
  name : String,
) -> BioRestrictionEnzyme raise SeqError {
  match get_bio_restriction_enzyme(name) {
    Some(e) => e
    None => raise SeqError("not a restriction enzyme: " + name)
  }
}

///|
/// All enzymes carried by a given supplier code (e.g. "N" for NEB, "B").
pub fn bio_restriction_enzymes_by_supplier(
  code : String,
) -> Array[BioRestrictionEnzyme] {
  let r : Array[BioRestrictionEnzyme] = []
  if code.length() != 1 {
    return r
  }
  let c = code.unsafe_get(0)
  for e in bio_restriction_enzyme_list {
    let mut found = false
    for i in 0..<e.supplier_codes.length() {
      if e.supplier_codes.unsafe_get(i) == c {
        found = true
        break
      }
    }
    if found {
      r.push(e)
    }
  }
  r.sort_by(fn(a, b) { a.name.compare(b.name) })
  r
}

// ---------------------------------------------------------------------------
// Sequence formatting and IUPAC matching
// ---------------------------------------------------------------------------

///|
/// Format a sequence like Biopython's FormattedSeq: delete whitespace and
/// digits, uppercase lowercase letters, reject any other character.
fn restriction_format_seq(seq : String) -> String raise SeqError {
  let sb = StringBuilder()
  for i in 0..<seq.length() {
    let c = seq.unsafe_get(i)
    let code = c.to_int()
    if code >= 97 && code <= 122 {
      sb.write_char((code - 32).unsafe_to_char())
    } else if code >= 65 && code <= 90 {
      sb.write_char(code.unsafe_to_char())
    } else if code == 32 ||
      code == 9 ||
      code == 10 ||
      code == 13 ||
      code == 11 ||
      code == 12 {
      // whitespace removed
    } else if code >= 48 && code <= 57 {
      // digits removed
    } else {
      raise SeqError("Invalid character found in " + seq)
    }
  }
  sb.to_string()
}

///|
/// Does a degenerate IUPAC pattern character match the concrete sequence
/// character? Pattern N matches any character, mirroring the regex "." of
/// Biopython compsite.
fn restriction_iupac_match(p : Int, s : Int) -> Bool {
  match p {
    65 => s == 65 // A
    67 => s == 67 // C
    71 => s == 71 // G
    84 => s == 84 // T
    // R = A|G, Y = C|T, S = G|C, W = A|T, K = G|T, M = A|C
    82 => s == 65 || s == 71
    89 => s == 67 || s == 84
    83 => s == 71 || s == 67
    87 => s == 65 || s == 84
    75 => s == 71 || s == 84
    77 => s == 65 || s == 67
    // B = C|G|T, D = A|G|T, H = A|C|T, V = A|C|G
    66 => s == 67 || s == 71 || s == 84
    68 => s == 65 || s == 71 || s == 84
    72 => s == 65 || s == 67 || s == 84
    86 => s == 65 || s == 67 || s == 71
    78 => true // N matches any character (regex ".")
    _ => s == p
  }
}

///|
fn restriction_complement_code(c : Int) -> Int {
  match c {
    65 => 84 // A -> T
    84 => 65 // T -> A
    67 => 71 // C -> G
    71 => 67 // G -> C
    82 => 89 // R -> Y
    89 => 82 // Y -> R
    83 => 83 // S
    87 => 87 // W
    75 => 77 // K -> M
    77 => 75 // M -> K
    66 => 86 // B -> V
    86 => 66 // V -> B
    68 => 72 // D -> H
    72 => 68 // H -> D
    78 => 78 // N
    _ => c
  }
}

///|
/// Reverse complement of an IUPAC recognition site alternative.
fn restriction_rc(site : StringView) -> String {
  let sb = StringBuilder()
  let mut i = site.length() - 1
  while i >= 0 {
    let rc = restriction_complement_code(site.unsafe_get(i).to_int())
    sb.write_char(rc.unsafe_to_char())
    i = i - 1
  }
  sb.to_string()
}

///|
fn restriction_window_match(
  scan : String,
  pos : Int,
  pattern : StringView,
) -> Bool {
  for j in 0..<pattern.length() {
    let p = pattern.unsafe_get(j).to_int()
    let s = scan.unsafe_get(pos + j).to_int()
    if !restriction_iupac_match(p, s) {
      return false
    }
  }
  true
}

///|
fn restriction_any_alt(scan : String, pos : Int, alts : Array[String]) -> Bool {
  for a in alts {
    if restriction_window_match(scan, pos, a) {
      return true
    }
  }
  false
}

///|
/// Number of double strand cuts: 0 (uncharacterised), 1 or 2.
fn restriction_ncuts(e : BioRestrictionEnzyme) -> Int {
  match e.fst5 {
    None => 0
    Some(_) =>
      match e.scd5 {
        None => 1
        Some(_) => 2
      }
  }
}

///|
/// Overhang kind: 0 unknown, 1 blunt, 2 = 5' overhang, 3 = 3' overhang.
fn restriction_end_kind(e : BioRestrictionEnzyme) -> Int {
  match e.ovhg {
    None => 0
    Some(0) => 1
    Some(v) => if v < 0 { 2 } else { 3 }
  }
}

///|
/// Core search on an already formatted sequence.
fn restriction_search_formatted(
  e : BioRestrictionEnzyme,
  f : String,
  linear : Bool,
) -> Array[Int] {
  let n = f.length()
  let size = e.size
  let alts : Array[String] = []
  for a in e.site.split("|") {
    alts.push(a.to_owned())
  }
  let rc_alts : Array[String] = []
  for a in alts {
    rc_alts.push(restriction_rc(a))
  }
  let scan = if linear || size <= 1 {
    f
  } else {
    let extra = if n >= size - 1 { size - 1 } else { n }
    f + substring(f, 0, extra)
  }
  let ncuts = restriction_ncuts(e)
  let forward : Array[Int] = []
  let reverse : Array[Int] = []
  let last = scan.length() - size
  let mut i = 0
  while i <= last {
    let fwd = restriction_any_alt(scan, i, alts)
    let rev = if e.palindromic {
      false
    } else {
      restriction_any_alt(scan, i, rc_alts)
    }
    // The compsite regex alternation makes the two branches mutually
    // exclusive at a given position.
    if fwd {
      let loc = i + 1 // 1-based start (leading-space convention)
      if ncuts == 0 {
        forward.push(loc)
      } else {
        match e.fst5 {
          Some(v) => forward.push(loc + v)
          None => ()
        }
        if ncuts == 2 {
          match e.scd5 {
            Some(v) => forward.push(loc + v)
            None => ()
          }
        }
      }
    } else if rev {
      let loc = i + 1
      if ncuts == 0 {
        reverse.push(loc)
      } else {
        match e.fst3 {
          Some(v) => reverse.push(loc - v)
          None => ()
        }
        if ncuts == 2 {
          match e.scd3 {
            Some(v) => reverse.push(loc - v)
            None => ()
          }
        }
      }
    }
    i = i + 1
  }
  let results : Array[Int] = []
  for v in forward {
    results.push(v)
  }
  for v in reverse {
    results.push(v)
  }
  // Palindromic results keep generation order; non-palindromic results are
  // sorted before applying the boundary rules.
  if !e.palindromic {
    results.sort()
  }
  // _drop: remove / wrap cuts outside the molecule.
  if linear {
    if ncuts != 0 {
      match e.ovhg {
        Some(ov) => {
          let kept : Array[Int] = []
          for c in results {
            if 1 < c && c <= n && 1 < c - ov && c - ov <= n {
              kept.push(c)
            }
          }
          return kept
        }
        None => ()
      }
    }
    return results
  }
  // circular: wrap the small prefix (<1) and large suffix (>n), no resort.
  let mut idx = 0
  while idx < results.length() && results[idx] < 1 {
    results[idx] = results[idx] + n
    idx = idx + 1
  }
  let mut tail = 0
  while tail < results.length() && results[results.length() - 1 - tail] > n {
    results[results.length() - 1 - tail] = results[results.length() - 1 - tail] -
      n
    tail = tail + 1
  }
  results
}

// ---------------------------------------------------------------------------
// Enzyme methods
// ---------------------------------------------------------------------------

///|
/// Search the cutting sites of the enzyme. Returns 1-based positions of the
/// first base of each 3' fragment. `linear=false` searches across the edges
/// of circular molecules.
pub fn BioRestrictionEnzyme::search(
  self : BioRestrictionEnzyme,
  seq : String,
  linear? : Bool = true,
) -> Array[Int] raise SeqError {
  let f = restriction_format_seq(seq)
  restriction_search_formatted(self, f, linear)
}

///|
/// Restrict `seq` into fragments. Uncharacterised enzymes raise `SeqError`.
pub fn BioRestrictionEnzyme::catalyse(
  self : BioRestrictionEnzyme,
  seq : String,
  linear? : Bool = true,
) -> Array[String] raise SeqError {
  let f = restriction_format_seq(seq)
  if restriction_ncuts(self) == 0 {
    raise SeqError(self.name + " restriction is unknown.")
  }
  let cuts = restriction_search_formatted(self, f, linear)
  let n = f.length()
  if cuts.length() == 0 {
    return [f]
  }
  let fragments : Array[String] = []
  if linear {
    let mut prev = 0
    for c in cuts {
      let start = c - 1
      fragments.push(substring(f, prev, start - prev))
      prev = start
    }
    fragments.push(substring(f, prev, n - prev))
  } else {
    // bridge the last site through the origin to the first site
    let first_start = cuts[0] - 1
    let last_start = cuts[cuts.length() - 1] - 1
    fragments.push(
      substring(f, last_start, n - last_start) + substring(f, 0, first_start),
    )
    let mut k = 1
    while k < cuts.length() {
      let a = cuts[k - 1] - 1
      let b = cuts[k] - 1
      fragments.push(substring(f, a, b - a))
      k = k + 1
    }
  }
  fragments
}

///|
/// Representation of the site with '^' for the (+)-strand cut and '_' for the
/// (-)-strand cut, e.g. EcoRI) -> "G^AATT_C".
pub fn BioRestrictionEnzyme::elucidate(self : BioRestrictionEnzyme) -> String {
  self.elucidation
}

///|
/// Theoretical cutting frequency (one cut per x bases).
pub fn BioRestrictionEnzyme::frequency(self : BioRestrictionEnzyme) -> Double {
  self.freq
}

///|
/// Characteristic tuple `(fst5, fst3, scd5, scd3, site)`.
pub fn BioRestrictionEnzyme::characteristic(
  self : BioRestrictionEnzyme,
) -> (Int?, Int?, Int?, Int?, String) {
  (self.fst5, self.fst3, self.scd5, self.scd3, self.site)
}

///|
pub fn BioRestrictionEnzyme::cut_once(self : BioRestrictionEnzyme) -> Bool {
  restriction_ncuts(self) == 1
}

///|
pub fn BioRestrictionEnzyme::cut_twice(self : BioRestrictionEnzyme) -> Bool {
  restriction_ncuts(self) == 2
}

///|
pub fn BioRestrictionEnzyme::is_palindromic(
  self : BioRestrictionEnzyme,
) -> Bool {
  self.palindromic
}

///|
pub fn BioRestrictionEnzyme::is_blunt(self : BioRestrictionEnzyme) -> Bool {
  restriction_end_kind(self) == 1
}

///|
pub fn BioRestrictionEnzyme::is_5overhang(self : BioRestrictionEnzyme) -> Bool {
  restriction_end_kind(self) == 2
}

///|
pub fn BioRestrictionEnzyme::is_3overhang(self : BioRestrictionEnzyme) -> Bool {
  restriction_end_kind(self) == 3
}

///|
pub fn BioRestrictionEnzyme::is_unknown_end(
  self : BioRestrictionEnzyme,
) -> Bool {
  restriction_end_kind(self) == 0
}

///|
pub fn BioRestrictionEnzyme::is_defined(self : BioRestrictionEnzyme) -> Bool {
  self.spec == 0
}

///|
pub fn BioRestrictionEnzyme::is_ambiguous(self : BioRestrictionEnzyme) -> Bool {
  self.spec == 1
}

///|
pub fn BioRestrictionEnzyme::is_not_defined(
  self : BioRestrictionEnzyme,
) -> Bool {
  self.spec == 2
}

///|
pub fn BioRestrictionEnzyme::is_commercial(self : BioRestrictionEnzyme) -> Bool {
  self.supplier_codes.length() > 0
}

///|
/// Overhang kind as a string: "blunt", "5' overhang", "3' overhang",
/// "unknown".
pub fn BioRestrictionEnzyme::overhang_type(
  self : BioRestrictionEnzyme,
) -> String {
  match restriction_end_kind(self) {
    1 => "blunt"
    2 => "5' overhang"
    3 => "3' overhang"
    _ => "unknown"
  }
}

///|
/// Supplier codes as a list of one-letter strings.
pub fn BioRestrictionEnzyme::supplier_list(
  self : BioRestrictionEnzyme,
) -> Array[String] {
  let r : Array[String] = []
  for i in 0..<self.supplier_codes.length() {
    r.push(self.supplier_codes[i:i + 1].to_owned())
  }
  r
}

///|
fn restriction_same_charac(
  a : BioRestrictionEnzyme,
  b : BioRestrictionEnzyme,
) -> Bool {
  a.fst5 == b.fst5 &&
  a.fst3 == b.fst3 &&
  a.scd5 == b.scd5 &&
  a.scd3 == b.scd3 &&
  a.site == b.site
}

///|
/// True if the two enzymes are equischizomers (same site, same cuts).
pub fn BioRestrictionEnzyme::is_equischizomer(
  self : BioRestrictionEnzyme,
  other : BioRestrictionEnzyme,
) -> Bool {
  restriction_same_charac(self, other)
}

///|
/// True if the two enzymes are neoschizomers (same site, different cuts).
pub fn BioRestrictionEnzyme::is_neoschizomer(
  self : BioRestrictionEnzyme,
  other : BioRestrictionEnzyme,
) -> Bool {
  self.site == other.site && !restriction_same_charac(self, other)
}

///|
/// True if the two enzymes are isoschizomers (same recognition site).
pub fn BioRestrictionEnzyme::is_isoschizomer(
  self : BioRestrictionEnzyme,
  other : BioRestrictionEnzyme,
) -> Bool {
  restriction_same_charac(self, other) ||
  (self.site == other.site && !restriction_same_charac(self, other))
}

///|
pub fn BioRestrictionEnzyme::equischizomers(
  self : BioRestrictionEnzyme,
) -> Array[BioRestrictionEnzyme] {
  let r : Array[BioRestrictionEnzyme] = []
  for e in bio_restriction_enzyme_list {
    if e.name != self.name && restriction_same_charac(self, e) {
      r.push(e)
    }
  }
  r.sort_by(fn(a, b) { a.name.compare(b.name) })
  r
}

///|
pub fn BioRestrictionEnzyme::neoschizomers(
  self : BioRestrictionEnzyme,
) -> Array[BioRestrictionEnzyme] {
  let r : Array[BioRestrictionEnzyme] = []
  for e in bio_restriction_enzyme_list {
    if e.name != self.name &&
      self.site == e.site &&
      !restriction_same_charac(self, e) {
      r.push(e)
    }
  }
  r.sort_by(fn(a, b) { a.name.compare(b.name) })
  r
}

///|
pub fn BioRestrictionEnzyme::isoschizomers(
  self : BioRestrictionEnzyme,
) -> Array[BioRestrictionEnzyme] {
  let r : Array[BioRestrictionEnzyme] = []
  for e in bio_restriction_enzyme_list {
    if e.name != self.name &&
      (
        restriction_same_charac(self, e) ||
        (self.site == e.site && !restriction_same_charac(self, e))
      ) {
      r.push(e)
    }
  }
  r.sort_by(fn(a, b) { a.name.compare(b.name) })
  r
}

// ---------------------------------------------------------------------------
// Compatible overhangs (Biopython % operator)
// ---------------------------------------------------------------------------

///|
/// Overhang matching classes used by Biopython's Ambiguous._mod2.
fn restriction_overhang_class(c : Int) -> String? {
  match c {
    82 => Some("ABDGHKMNSRWV") // R
    89 => Some("CBDHKMNSTWVY") // Y
    87 => Some("ABDHKMNRTWVY") // W
    83 => Some("CBDGHKMNSRVY") // S
    77 => Some("ACBDHMNSRWVY") // M
    75 => Some("BDGHKNSRTWVY") // K
    72 => Some("ACBDHKMNSRTWVY") // H
    66 => Some("CBDGHKMNSRTWVY") // B
    86 => Some("ACBDGHKMNSRWVY") // V
    68 => Some("ABDGHKMNSRTWVY") // D
    _ => None
  }
}

///|
/// Interpret the (possibly rewritten) pattern against the target overhang,
/// equivalent to re.match(pattern, target) with equal-length strings.
fn restriction_regex_like_match(pattern : String, target : String) -> Bool {
  let mut pi = 0
  let mut ti = 0
  while pi < pattern.length() {
    let c = pattern.unsafe_get(pi).to_int()
    if c == 91 {
      // '[' opens a character class
      let mut end = pi + 1
      let mut closed = false
      while end < pattern.length() {
        if pattern.unsafe_get(end).to_int() == 93 {
          closed = true
          break
        }
        end = end + 1
      }
      if !closed {
        return false
      }
      if ti >= target.length() {
        return false
      }
      let tc = target.unsafe_get(ti).to_int()
      let mut ok = false
      for k in (pi + 1)..<end {
        if pattern.unsafe_get(k).to_int() == tc {
          ok = true
          break
        }
      }
      if !ok {
        return false
      }
      pi = end + 1
    } else if c == 46 {
      // '.' matches one character
      if ti >= target.length() {
        return false
      }
      pi = pi + 1
    } else {
      if ti >= target.length() || target.unsafe_get(ti).to_int() != c {
        return false
      }
      pi = pi + 1
    }
    ti = ti + 1
  }
  ti == target.length()
}

///|
/// Literal transcription of Biopython Ambiguous._mod2's string rewriting,
/// including its quirk of splitting previously inserted character classes.
fn restriction_amb_overhang_match(pattern : String, target : String) -> Bool {
  if pattern.length() != target.length() {
    return false
  }
  let mut se = pattern
  for i in 0..<pattern.length() {
    let base = pattern.unsafe_get(i).to_int()
    if base == 65 || base == 84 || base == 67 || base == 71 {
      // concrete base, nothing to rewrite
      ()
    } else if base == 78 {
      se = se.replace(old="N", new=".")
    } else {
      match restriction_overhang_class(base) {
        Some(cls) => {
          let expanded = "[" + cls + "]"
          let one = pattern[i:i + 1].to_owned()
          se = se.replace(old=one, new=expanded)
        }
        None => ()
      }
    }
  }
  restriction_regex_like_match(se, target)
}

///|
/// True if the two enzymes produce compatible overhangs (Biopython `%`).
pub fn BioRestrictionEnzyme::compatible_overhang(
  self : BioRestrictionEnzyme,
  other : BioRestrictionEnzyme,
) -> Bool {
  let ke = restriction_end_kind(self)
  let ko = restriction_end_kind(other)
  if ke == 0 || ko == 0 {
    return false
  }
  if ke == 1 && ko == 1 {
    return true
  }
  if ke != ko {
    return false
  }
  if self.spec == 1 {
    return restriction_amb_overhang_match(self.ovhgseq, other.ovhgseq)
  }
  if self.ovhgseq == other.ovhgseq {
    return true
  }
  if other.spec == 1 {
    return restriction_amb_overhang_match(other.ovhgseq, self.ovhgseq)
  }
  false
}

///|
/// All enzymes in the database producing ends compatible with this enzyme.
pub fn BioRestrictionEnzyme::compatible_ends(
  self : BioRestrictionEnzyme,
) -> Array[BioRestrictionEnzyme] {
  let r : Array[BioRestrictionEnzyme] = []
  let ke = restriction_end_kind(self)
  if ke == 0 {
    return r
  }
  for e in bio_restriction_enzyme_list {
    let should_include = if ke == 1 {
      restriction_end_kind(e) == 1
    } else {
      // Biopython iterates x and evaluates x % self.
      restriction_end_kind(e) == ke && e.compatible_overhang(self)
    }
    if should_include {
      r.push(e)
    }
  }
  r.sort_by(fn(a, b) { a.name.compare(b.name) })
  r
}

///|
/// Free-function form of `enzyme.search`.
pub fn bio_restriction_search(
  enzyme : BioRestrictionEnzyme,
  seq : String,
  linear? : Bool = true,
) -> Array[Int] raise SeqError {
  enzyme.search(seq, linear~)
}

///|
/// Free-function form of `enzyme.catalyse`.
pub fn bio_restriction_catalyse(
  enzyme : BioRestrictionEnzyme,
  seq : String,
  linear? : Bool = true,
) -> Array[String] raise SeqError {
  enzyme.catalyse(seq, linear~)
}

// ---------------------------------------------------------------------------
// RestrictionBatch
// ---------------------------------------------------------------------------

///|
/// A batch of restriction enzymes (Biopython RestrictionBatch).
pub struct RestrictionBatch {
  names : Array[String]
}

///|
fn restriction_batch_dedup(names : Array[String]) -> Array[String] {
  let seen : Map[String, Unit] = Map([], capacity=names.length() + 1)
  let r : Array[String] = []
  for n in names {
    if !seen.contains(n) {
      seen[n] = ()
      r.push(n)
    }
  }
  r
}

///|
/// Build a batch from enzyme names. Raises `SeqError` on unknown names.
pub fn RestrictionBatch::new(
  names~ : Array[String],
) -> RestrictionBatch raise SeqError {
  for n in names {
    if get_bio_restriction_enzyme(n) is None {
      raise SeqError("not a restriction enzyme: " + n)
    }
  }
  RestrictionBatch::{ names: restriction_batch_dedup(names), }
}

///|
/// Batch containing every enzyme in the database.
pub fn all_enzymes() -> RestrictionBatch {
  let names : Array[String] = []
  for e in bio_restriction_enzyme_list {
    names.push(e.name)
  }
  RestrictionBatch::{ names, }
}

///|
/// Batch of commercially available enzymes.
pub fn commercial_enzymes() -> RestrictionBatch {
  let names : Array[String] = []
  for e in bio_restriction_enzyme_list {
    if e.is_commercial() {
      names.push(e.name)
    }
  }
  RestrictionBatch::{ names, }
}

///|
/// Batch of enzymes without a listed supplier.
pub fn noncommercial_enzymes() -> RestrictionBatch {
  let names : Array[String] = []
  for e in bio_restriction_enzyme_list {
    if !e.is_commercial() {
      names.push(e.name)
    }
  }
  RestrictionBatch::{ names, }
}

///|
pub fn RestrictionBatch::length(self : RestrictionBatch) -> Int {
  self.names.length()
}

///|
pub fn RestrictionBatch::contains(
  self : RestrictionBatch,
  name : String,
) -> Bool {
  self.names.contains(name)
}

///|
pub fn RestrictionBatch::add(
  self : RestrictionBatch,
  name : String,
) -> Unit raise SeqError {
  if get_bio_restriction_enzyme(name) is None {
    raise SeqError("not a restriction enzyme: " + name)
  }
  if !self.names.contains(name) {
    self.names.push(name)
  }
}

///|
pub fn RestrictionBatch::add_supplier(
  self : RestrictionBatch,
  code : String,
) -> Unit {
  for e in bio_restriction_enzymes_by_supplier(code) {
    if !self.names.contains(e.name) {
      self.names.push(e.name)
    }
  }
}

///|
pub fn RestrictionBatch::remove(self : RestrictionBatch, name : String) -> Unit {
  let mut idx = -1
  let mut i = 0
  while i < self.names.length() {
    if self.names[i] == name {
      idx = i
      i = self.names.length()
    }
    i = i + 1
  }
  if idx >= 0 {
    let last = self.names.length() - 1
    for j in idx..<last {
      self.names[j] = self.names[j + 1]
    }
    let _ = self.names.pop()
  }
}

///|
/// Lexicographic (code point) comparison, matching Python's ``sorted`` for
/// enzyme names. MoonBit's native ``String`` ordering is length-first, so the
/// default ``sort`` cannot be used here.
fn restriction_name_cmp(a : String, b : String) -> Int {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut cmp = 0
  let mut i = 0
  while i < n && cmp == 0 {
    let ca = a.unsafe_get(i).to_int()
    let cb = b.unsafe_get(i).to_int()
    if ca < cb {
      cmp = -1
    } else if ca > cb {
      cmp = 1
    }
    i = i + 1
  }
  if cmp != 0 {
    cmp
  } else if a.length() < b.length() {
    -1
  } else if a.length() > b.length() {
    1
  } else {
    0
  }
}

///|
/// Enzyme names, sorted alphabetically.
pub fn RestrictionBatch::elements(self : RestrictionBatch) -> Array[String] {
  let r = self.names.copy()
  r.sort_by(restriction_name_cmp)
  r
}

///|
/// Enzyme records of the batch.
pub fn RestrictionBatch::enzymes(
  self : RestrictionBatch,
) -> Array[BioRestrictionEnzyme] {
  let r : Array[BioRestrictionEnzyme] = []
  for n in self.names {
    match get_bio_restriction_enzyme(n) {
      Some(e) => r.push(e)
      None => ()
    }
  }
  r
}

///|
/// Search `seq` once for every enzyme of the batch, returning a map
/// enzyme name -> 1-based cutting positions.
pub fn RestrictionBatch::search(
  self : RestrictionBatch,
  seq : String,
  linear? : Bool = true,
) -> Map[String, Array[Int]] raise SeqError {
  let f = restriction_format_seq(seq)
  let m : Map[String, Array[Int]] = Map([], capacity=self.names.length() + 1)
  for n in self.names {
    let e = require_bio_restriction_enzyme(n)
    m[n] = restriction_search_formatted(e, f, linear)
  }
  m
}

// ---------------------------------------------------------------------------
// RestrictionAnalysis
// ---------------------------------------------------------------------------

///|
/// Enhanced analysis of a sequence with a RestrictionBatch, mirroring
/// Biopython's Analysis class.
pub struct RestrictionAnalysis {
  batch : RestrictionBatch
  sequence_length : Int
  results : Map[String, Array[Int]]
}

///|
pub fn RestrictionAnalysis::new(
  batch~ : RestrictionBatch,
  seq~ : String,
  linear? : Bool = true,
) -> RestrictionAnalysis raise SeqError {
  let f = restriction_format_seq(seq)
  let results : Map[String, Array[Int]] = Map([], capacity=batch.length() + 1)
  for n in batch.names {
    let e = require_bio_restriction_enzyme(n)
    results[n] = restriction_search_formatted(e, f, linear)
  }
  RestrictionAnalysis::{ batch, sequence_length: f.length(), results, }
}

///|
fn restriction_copy_map(m : Map[String, Array[Int]]) -> Map[String, Array[Int]] {
  let r : Map[String, Array[Int]] = Map([], capacity=m.length() + 1)
  for k, v in m {
    r[k] = v
  }
  r
}

///|
/// Full enzyme -> cut sites map.
pub fn RestrictionAnalysis::full(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  restriction_copy_map(self.results)
}

///|
pub fn RestrictionAnalysis::with_sites(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  let r : Map[String, Array[Int]] = Map([], capacity=self.results.length() + 1)
  for k, v in self.results {
    if v.length() > 0 {
      r[k] = v
    }
  }
  r
}

///|
pub fn RestrictionAnalysis::without_site(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  let r : Map[String, Array[Int]] = Map([], capacity=self.results.length() + 1)
  for k, v in self.results {
    if v.length() == 0 {
      r[k] = v
    }
  }
  r
}

///|
pub fn RestrictionAnalysis::with_n_sites(
  self : RestrictionAnalysis,
  n : Int,
) -> Map[String, Array[Int]] {
  let r : Map[String, Array[Int]] = Map([], capacity=self.results.length() + 1)
  for k, v in self.results {
    if v.length() == n {
      r[k] = v
    }
  }
  r
}

///|
fn restriction_analysis_filter(
  analysis : RestrictionAnalysis,
  pred : (BioRestrictionEnzyme) -> Bool,
) -> Map[String, Array[Int]] {
  let r : Map[String, Array[Int]] = Map(
    [],
    capacity=analysis.results.length() + 1,
  )
  for k, v in analysis.results {
    match get_bio_restriction_enzyme(k) {
      Some(e) => if pred(e) { r[k] = v }
      None => ()
    }
  }
  r
}

///|
pub fn RestrictionAnalysis::blunt(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  restriction_analysis_filter(self, fn(e) { e.is_blunt() })
}

///|
pub fn RestrictionAnalysis::overhang5(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  restriction_analysis_filter(self, fn(e) { e.is_5overhang() })
}

///|
pub fn RestrictionAnalysis::overhang3(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  restriction_analysis_filter(self, fn(e) { e.is_3overhang() })
}

///|
pub fn RestrictionAnalysis::defined(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  restriction_analysis_filter(self, fn(e) { e.is_defined() })
}

///|
pub fn RestrictionAnalysis::ambiguous(
  self : RestrictionAnalysis,
) -> Map[String, Array[Int]] {
  restriction_analysis_filter(self, fn(e) { e.is_ambiguous() })
}

///|
/// Normalise 1-based boundaries (negative indices count from the end,
/// swap if reversed).
fn restriction_boundaries(
  analysis : RestrictionAnalysis,
  start : Int,
  end : Int,
) -> (Int, Int)? {
  let mut s = start
  let mut e = end
  if s < 1 {
    s = s + analysis.sequence_length
  }
  if e < 1 {
    e = e + analysis.sequence_length
  }
  if s >= e {
    let t = s
    s = e
    e = t
  }
  if s < e {
    Some((s, e))
  } else {
    None
  }
}

///|
/// Enzymes whose cuts are all within `[start, end)` (at least one cut).
pub fn RestrictionAnalysis::only_between(
  self : RestrictionAnalysis,
  start : Int,
  end : Int,
) -> Map[String, Array[Int]] {
  match restriction_boundaries(self, start, end) {
    Some((s, e)) => {
      let r : Map[String, Array[Int]] = Map(
        [],
        capacity=self.results.length() + 1,
      )
      for k, v in self.results {
        if v.length() == 0 {
          continue
        }
        let mut all_inside = true
        for site in v {
          if !(s <= site && site < e) {
            all_inside = false
            break
          }
        }
        if all_inside {
          r[k] = v
        }
      }
      r
    }
    None => Map([])
  }
}

///|
/// Enzymes that cut at least once within `[start, end)` (all cuts returned).
pub fn RestrictionAnalysis::between(
  self : RestrictionAnalysis,
  start : Int,
  end : Int,
) -> Map[String, Array[Int]] {
  match restriction_boundaries(self, start, end) {
    Some((s, e)) => {
      let r : Map[String, Array[Int]] = Map(
        [],
        capacity=self.results.length() + 1,
      )
      for k, v in self.results {
        let mut inside = false
        for site in v {
          if s <= site && site < e {
            inside = true
            break
          }
        }
        if inside {
          r[k] = v
        }
      }
      r
    }
    None => Map([])
  }
}

///|
/// Only the cut positions inside `[start, end]` for enzymes cutting there.
/// Reversed boundaries (`start > end`) select the circular complement
/// region: `site >= start || site <= end`.
pub fn RestrictionAnalysis::show_only_between(
  self : RestrictionAnalysis,
  start : Int,
  end : Int,
) -> Map[String, Array[Int]] {
  let candidate = self.between(start, end)
  let r : Map[String, Array[Int]] = Map([], capacity=candidate.length() + 1)
  for k, v in candidate {
    let kept : Array[Int] = []
    for site in v {
      let ok = if start <= end {
        start <= site && site <= end
      } else {
        site >= start || site <= end
      }
      if ok {
        kept.push(site)
      }
    }
    r[k] = kept
  }
  r
}

///|
/// Enzymes with at least one cut, none of which lie within `[start, end)`.
pub fn RestrictionAnalysis::only_outside(
  self : RestrictionAnalysis,
  start : Int,
  end : Int,
) -> Map[String, Array[Int]] {
  match restriction_boundaries(self, start, end) {
    Some((s, e)) => {
      let r : Map[String, Array[Int]] = Map(
        [],
        capacity=self.results.length() + 1,
      )
      for k, v in self.results {
        if v.length() == 0 {
          continue
        }
        let mut any_inside = false
        for site in v {
          if s <= site && site < e {
            any_inside = true
            break
          }
        }
        if !any_inside {
          r[k] = v
        }
      }
      r
    }
    None => Map([])
  }
}

///|
/// Enzymes with at least one cut outside `[start, end)`.
pub fn RestrictionAnalysis::outside(
  self : RestrictionAnalysis,
  start : Int,
  end : Int,
) -> Map[String, Array[Int]] {
  match restriction_boundaries(self, start, end) {
    Some((s, e)) => {
      let r : Map[String, Array[Int]] = Map(
        [],
        capacity=self.results.length() + 1,
      )
      for k, v in self.results {
        let mut any_outside = false
        for site in v {
          if !(s <= site && site < e) {
            any_outside = true
            break
          }
        }
        if any_outside {
          r[k] = v
        }
      }
      r
    }
    None => Map([])
  }
}

///|
/// Enzymes that do not cut between `start` and `end`: enzymes without any
/// site plus enzymes cutting only outside the region.
pub fn RestrictionAnalysis::do_not_cut(
  self : RestrictionAnalysis,
  start : Int,
  end : Int,
) -> Map[String, Array[Int]] {
  let r = self.without_site()
  let outside_map = self.only_outside(start, end)
  for k, v in outside_map {
    r[k] = v
  }
  r
}
