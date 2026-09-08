///|
/// universalmotif module, inspired by the Bioconductor `universalmotif` package.
///
/// Port of the unified motif object and its core operations:
///   - Matrix type conversions (PCM <-> PPM <-> PWM <-> ICM)
///   - Information content scoring (Schneider default + relative-entropy/KL mode)
///   - IUPAC consensus string (probability-threshold based, matching `get_consensusC`)
///   - Motif comparison with 12 metrics (PCC/EUCL/SW/KL/WEUCL/ALLR/BHAT/HELL/WPCC/
///     SEUCL/MAN/ALLR_LL), sliding alignment, score strategies, reverse-complement
///     try, score normalization
///   - Motif merging (weighted arithmetic mean of PPMs)
///   - Motif trimming (information-content based edge trimming)
///   - MEME and JASPAR read/write
///
/// Matrix convention: `motif` is stored row-major with **row = motif position,
/// column = alphabet letter** (matching the rest of this library). The alphabet
/// order for DNA is `[A, C, G, T]`, RNA `[A, C, G, U]`, amino acids the 20
/// standard residues, custom alphabets one letter per character.

///|
/// Type of the matrix held in a `UniversalMotif`.
pub(all) enum UMotifType {
  PCM
  PPM
  PWM
  ICM
} derive(Eq)

///|
/// Comparison metrics supported by `compare_motifs`.
///
/// Similarity metrics (higher is better): `PCC`, `SW`, `BHAT`, `WPCC`.
/// Distance metrics (lower is better): `EUCL`, `KL`, `WEUCL`, `ALLR`,
/// `ALLR_LL`, `HELL`, `SEUCL`, `MAN`.
pub(all) enum UMetric {
  PCC
  EUCL
  SW
  KL
  WEUCL
  ALLR
  BHAT
  HELL
  WPCC
  SEUCL
  MAN
  ALLR_LL
}

///|
/// Strategy for combining per-position comparison scores.
pub(all) enum UScoreStrat {
  USum
  UAMean
  UGMean
  UMedian
  UFZT
  UWAMean
  UWGMean
}

///|
/// Which edges to trim in `trim_motif`.
pub(all) enum UTrimFrom {
  UBoth
  ULeft
  URight
} derive(Eq)

///|
/// The unified motif object.
///
/// Fields mirror the slots of the Bioconductor S4 `universalmotif` class. The
/// `motif` matrix is the single source of truth, disambiguated by `motif_type`.
pub struct UniversalMotif {
  pub name : String
  pub altname : String
  pub family : String
  pub organism : String
  pub motif : Array[Array[Double]]
  pub alphabet : String
  pub motif_type : UMotifType
  pub icscore : Double
  pub nsites : Double
  pub pseudocount : Double
  pub bkg : Array[Double]
  pub consensus : String
  pub strand : String
  pub bkgsites : Double
  pub extrainfo : String
}

///|
/// Construct a `UniversalMotif`. `consensus` and `icscore` are computed
/// automatically from the matrix (the IC is computed on the PPM conversion of
/// the matrix according to `motif_type`).
pub fn UniversalMotif::new(
  name : String,
  alphabet : String,
  motif : Array[Array[Double]],
  motif_type : UMotifType,
  nsites? : Double = 100.0,
  pseudocount? : Double = 0.8,
  bkg? : Array[Double] = [],
  strand? : String = "+-",
  altname? : String = "",
  family? : String = "",
  organism? : String = "",
  bkgsites? : Double = 0.0,
  extrainfo? : String = "",
) -> UniversalMotif {
  let letters = um_alphabet_letters(alphabet)
  let n = letters.length()
  let bg = if bkg.length() == n { bkg } else { um_default_bkg(n) }
  let ppm = um_to_ppm(motif, motif_type, pseudocount, nsites, bg)
  let icscore = universal_icscore(ppm, bg, false)
  let consensus = consensus_string(ppm, alphabet)
  UniversalMotif::{
    name,
    altname,
    family,
    organism,
    motif,
    alphabet,
    motif_type,
    icscore,
    nsites,
    pseudocount,
    bkg: bg,
    consensus,
    strand,
    bkgsites,
    extrainfo,
  }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

///|
/// Return the ordered single-letter alphabet for the given alphabet string.
///
/// DNA -> `[A,C,G,T]`, RNA -> `[A,C,G,U]`, AA -> the 20 standard residues,
/// otherwise each character of the string.
pub fn um_alphabet_letters(alphabet : String) -> Array[String] {
  if alphabet == "DNA" {
    return ["A", "C", "G", "T"]
  }
  if alphabet == "RNA" {
    return ["A", "C", "G", "U"]
  }
  if alphabet == "AA" {
    return [
      "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F",
      "P", "S", "T", "W", "Y", "V",
    ]
  }
  let result = Array::new()
  let mut i = 0
  while i < alphabet.length() {
    result.push(alphabet[i:i + 1].to_owned())
    i = i + 1
  }
  result
}

///|
/// Number of letters in the alphabet.
pub fn um_alphabet_size(alphabet : String) -> Int {
  um_alphabet_letters(alphabet).length()
}

///|
/// Default uniform background (1/n per letter).
pub fn um_default_bkg(n : Int) -> Array[Double] {
  if n <= 0 {
    return []
  }
  let p = 1.0 / n.to_double()
  Array::make(n, p)
}

///|
/// `ln` that returns 0 for non-positive input (used where `0 * ln(0)` should
/// vanish in entropy-like sums).
pub fn um_ln(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  @math.ln(x)
}

///|
/// Square root (thin wrapper around `Double::sqrt`).
pub fn um_sqrt(x : Double) -> Double {
  x.sqrt()
}

///|
/// Whether a metric is a similarity (higher is better) rather than a distance.
pub fn um_is_similarity(metric : UMetric) -> Bool {
  match metric {
    PCC => true
    SW => true
    BHAT => true
    WPCC => true
    _ => false
  }
}

///|
/// Deep-copy a matrix.
fn um_copy_matrix(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let out = Array::new()
  for i = 0; i < m.length(); i = i + 1 {
    out.push(m[i].copy())
  }
  out
}

///|
/// Allocate an `rows x cols` Double matrix with independent rows.
fn um_make_matrix(
  rows : Int,
  cols : Int,
  default : Double,
) -> Array[Array[Double]] {
  let m = Array::new()
  let mut i = 0
  while i < rows {
    m.push(Array::make(cols, default))
    i = i + 1
  }
  m
}

///|
/// Allocate an `rows x cols` Int matrix with independent rows.
fn um_make_matrix_int(
  rows : Int,
  cols : Int,
  default : Int,
) -> Array[Array[Int]] {
  let m = Array::new()
  let mut i = 0
  while i < rows {
    m.push(Array::make(cols, default))
    i = i + 1
  }
  m
}

// ---------------------------------------------------------------------------
// Matrix type conversions
// ---------------------------------------------------------------------------

///|
/// PCM -> PPM. Pseudocount is split evenly across letters: PPM[j] = (C[j] +
/// p/n) / (sum + p). When `pseudocount == 0` no smoothing is applied.
pub fn pcm_to_ppm(
  pcm : Array[Array[Double]],
  pseudocount : Double,
) -> Array[Array[Double]] {
  if pcm.length() == 0 {
    return []
  }
  let n = pcm[0].length()
  let ppm = um_make_matrix(pcm.length(), n, 0.0)
  for pos = 0; pos < pcm.length(); pos = pos + 1 {
    let row = pcm[pos]
    let mut total = 0.0
    for j = 0; j < n; j = j + 1 {
      total = total + row[j]
    }
    let denom = if pseudocount > 0.0 {
      total + pseudocount
    } else {
      total
    }
    if denom <= 0.0 {
      let p = 1.0 / n.to_double()
      for j = 0; j < n; j = j + 1 {
        ppm[pos][j] = p
      }
    } else {
      for j = 0; j < n; j = j + 1 {
        let num = if pseudocount > 0.0 {
          row[j] + pseudocount / n.to_double()
        } else {
          row[j]
        }
        ppm[pos][j] = num / denom
      }
    }
  }
  ppm
}

///|
/// PPM -> PCM: `PCM[j] = round(PPM[j] * nsites)`; the column is then adjusted
/// so it sums to exactly `nsites` by adding the residual to the largest entry.
pub fn ppm_to_pcm(
  ppm : Array[Array[Double]],
  nsites : Double,
) -> Array[Array[Int]] {
  let n_sites = if nsites <= 1.0 { 100.0 } else { nsites }
  if ppm.length() == 0 {
    return []
  }
  let n = ppm[0].length()
  let pcm = um_make_matrix_int(ppm.length(), n, 0)
  for pos = 0; pos < ppm.length(); pos = pos + 1 {
    let row = ppm[pos]
    let col = Array::make(n, 0)
    let mut sum = 0
    let mut max_idx = 0
    let mut max_val = -1.0
    for j = 0; j < n; j = j + 1 {
      let v = (row[j] * n_sites + 0.5).floor().to_int()
      col[j] = v
      sum = sum + v
      if row[j] > max_val {
        max_val = row[j]
        max_idx = j
      }
    }
    let target = n_sites.to_int()
    col[max_idx] = col[max_idx] + (target - sum)
    for j = 0; j < n; j = j + 1 {
      pcm[pos][j] = col[j]
    }
  }
  pcm
}

///|
/// PPM -> PWM (log-odds). Faithful to `ppm_to_pwmC`: when `pseudocount > 0`
/// the PPM is round-tripped through PCM (with `nsites`) and back to PPM (with
/// `pseudocount`) before taking `log2(PPM / bkg)`. Non-finite results (from a
/// zero probability) are clamped to a large negative number.
pub fn ppm_to_pwm(
  ppm : Array[Array[Double]],
  bkg : Array[Double],
  pseudocount : Double,
  nsites : Double,
) -> Array[Array[Double]] {
  if ppm.length() == 0 {
    return []
  }
  let n = ppm[0].length()
  let bg = if bkg.length() == n { bkg } else { um_default_bkg(n) }
  let pos = if pseudocount > 0.0 {
    let pcm = ppm_to_pcm(ppm, nsites)
    let pcm_d = um_make_matrix(pcm.length(), n, 0.0)
    for i = 0; i < pcm.length(); i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        pcm_d[i][j] = pcm[i][j].to_double()
      }
    }
    pcm_to_ppm(pcm_d, pseudocount)
  } else {
    um_copy_matrix(ppm)
  }
  let pwm = um_make_matrix(ppm.length(), n, 0.0)
  for i = 0; i < ppm.length(); i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if pos[i][j] > 0.0 && bg[j] > 0.0 {
        pwm[i][j] = @math.log2(pos[i][j] / bg[j])
      } else {
        pwm[i][j] = -1.0e10
      }
    }
  }
  pwm
}

///|
/// PWM -> PPM (`pwm_to_ppmC`): `2^PWM`, multiply by background, then renormalize
/// the column to sum to 1.
pub fn pwm_to_ppm(
  pwm : Array[Array[Double]],
  bkg : Array[Double],
) -> Array[Array[Double]] {
  if pwm.length() == 0 {
    return []
  }
  let n = pwm[0].length()
  let bg = if bkg.length() == n { bkg } else { um_default_bkg(n) }
  let ppm = um_make_matrix(pwm.length(), n, 0.0)
  for pos = 0; pos < pwm.length(); pos = pos + 1 {
    let mut total = 0.0
    for j = 0; j < n; j = j + 1 {
      let v = @math.exp(pwm[pos][j] * 0.6931471805599453) // 2^x = e^(x ln2)
      ppm[pos][j] = v
      total = total + v
    }
    if total > 0.0 && (total < 0.99 || total > 1.01) {
      let mut total2 = 0.0
      for j = 0; j < n; j = j + 1 {
        ppm[pos][j] = ppm[pos][j] * bg[j]
        total2 = total2 + ppm[pos][j]
      }
      if total2 > 0.0 {
        for j = 0; j < n; j = j + 1 {
          ppm[pos][j] = ppm[pos][j] / total2
        }
      }
    }
  }
  ppm
}

///|
/// PPM -> ICM (information content matrix). Default Schneider mode:
/// `ICM[j] = PPM[j] * (log2(n) - H)` where `H = -sum(PPM * log2(PPM))`.
/// When `relative_entropy = true`, the KL/relative-entropy form is used:
/// `ICM[j] = PPM[j] * log2(PPM[j] / bkg[j])` (clamped to >= 0).
pub fn ppm_to_icm(
  ppm : Array[Array[Double]],
  bkg : Array[Double],
  relative_entropy : Bool,
) -> Array[Array[Double]] {
  if ppm.length() == 0 {
    return []
  }
  let n = ppm[0].length()
  let bg = if bkg.length() == n { bkg } else { um_default_bkg(n) }
  let icm = um_make_matrix(ppm.length(), n, 0.0)
  for pos = 0; pos < ppm.length(); pos = pos + 1 {
    let row = ppm[pos]
    if relative_entropy {
      for j = 0; j < n; j = j + 1 {
        let val = if row[j] > 0.0 && bg[j] > 0.0 {
          row[j] * @math.log2(row[j] / bg[j])
        } else {
          0.0
        }
        icm[pos][j] = if val < 0.0 { 0.0 } else { val }
      }
    } else {
      let mut ent = 0.0 // sum of p*log2(p) (<= 0)
      for j = 0; j < n; j = j + 1 {
        ent = ent + (if row[j] > 0.0 { row[j] * @math.log2(row[j]) } else { 0.0 })
      }
      let h = -ent // H = entropy
      let total_ic = @math.log2(n.to_double()) - h
      for j = 0; j < n; j = j + 1 {
        icm[pos][j] = row[j] * total_ic
      }
    }
  }
  icm
}

///|
/// ICM -> PPM: renormalize each column to sum to 1.
pub fn icm_to_ppm(icm : Array[Array[Double]]) -> Array[Array[Double]] {
  if icm.length() == 0 {
    return []
  }
  let n = icm[0].length()
  let ppm = um_make_matrix(icm.length(), n, 0.0)
  for pos = 0; pos < icm.length(); pos = pos + 1 {
    let mut total = 0.0
    for j = 0; j < n; j = j + 1 {
      total = total + icm[pos][j]
    }
    if total <= 0.0 {
      let p = 1.0 / n.to_double()
      for j = 0; j < n; j = j + 1 {
        ppm[pos][j] = p
      }
    } else {
      for j = 0; j < n; j = j + 1 {
        ppm[pos][j] = icm[pos][j] / total
      }
    }
  }
  ppm
}

///|
/// Per-position information content of a single PPM column.
/// Default: `log2(n) - H`. Relative-entropy mode: `sum(PPM * log2(PPM/bkg))`.
pub fn position_icscore_ppm(
  row : Array[Double],
  bkg : Array[Double],
  relative_entropy : Bool,
) -> Double {
  let n = row.length()
  if n == 0 {
    return 0.0
  }
  let bg = if bkg.length() == n { bkg } else { um_default_bkg(n) }
  if relative_entropy {
    let mut s = 0.0
    for j = 0; j < n; j = j + 1 {
      let v = if row[j] > 0.0 && bg[j] > 0.0 {
        row[j] * @math.log2(row[j] / bg[j])
      } else {
        0.0
      }
      s = s + (if v < 0.0 { 0.0 } else { v })
    }
    s
  } else {
    let mut ent = 0.0
    for j = 0; j < n; j = j + 1 {
      ent = ent + (if row[j] > 0.0 { row[j] * @math.log2(row[j]) } else { 0.0 })
    }
    let h = -ent
    @math.log2(n.to_double()) - h
  }
}

///|
/// Total information content of a PPM (sum of per-position IC).
pub fn universal_icscore(
  ppm : Array[Array[Double]],
  bkg : Array[Double],
  relative_entropy : Bool,
) -> Double {
  let mut total = 0.0
  for pos = 0; pos < ppm.length(); pos = pos + 1 {
    total = total + position_icscore_ppm(ppm[pos], bkg, relative_entropy)
  }
  total
}

///|
/// Convert a motif matrix to PPM according to its declared type.
fn um_to_ppm(
  motif : Array[Array[Double]],
  motif_type : UMotifType,
  pseudocount : Double,
  nsites : Double,
  bkg : Array[Double],
) -> Array[Array[Double]] {
  match motif_type {
    PCM => pcm_to_ppm(motif, pseudocount)
    PPM => um_copy_matrix(motif)
    PWM => pwm_to_ppm(motif, bkg)
    ICM => icm_to_ppm(motif)
  }
}

// ---------------------------------------------------------------------------
// Consensus / IUPAC
// ---------------------------------------------------------------------------

///|
/// Consensus letter for a single PPM position (probability-threshold IUPAC,
/// matching `get_consensusC`). Only DNA/RNA use the full IUPAC table; for other
/// alphabets the single highest-probability letter is returned.
pub fn get_consensus(row : Array[Double], alphabet : String) -> String {
  let letters = um_alphabet_letters(alphabet)
  let n = if row.length() < letters.length() {
    row.length()
  } else {
    letters.length()
  }
  if n == 0 {
    return ""
  }
  if alphabet != "DNA" && alphabet != "RNA" {
    let mut max_idx = 0
    let mut max_val = -1.0
    for j = 0; j < n; j = j + 1 {
      if row[j] > max_val {
        max_val = row[j]
        max_idx = j
      }
    }
    return letters[max_idx]
  }
  // pos = [A, C, G, T]
  let a = if n > 0 { row[0] } else { 0.0 }
  let c = if n > 1 { row[1] } else { 0.0 }
  let g = if n > 2 { row[2] } else { 0.0 }
  let t = if n > 3 { row[3] } else { 0.0 }
  let u = if alphabet == "RNA" { "U" } else { "T" }
  // single-letter strict dominant
  if a > 0.5 && c <= 0.25 && g <= 0.25 && t <= 0.25 {
    return "A"
  }
  if c > 0.5 && a <= 0.25 && g <= 0.25 && t <= 0.25 {
    return "C"
  }
  if g > 0.5 && a <= 0.25 && c <= 0.25 && t <= 0.25 {
    return "G"
  }
  if t > 0.5 && a <= 0.25 && c <= 0.25 && g <= 0.25 {
    return u
  }
  // two-letter (sum > 0.75)
  if a + c > 0.75 {
    return "M"
  }
  if a + g > 0.75 {
    return "R"
  }
  if a + t > 0.75 {
    return "W"
  }
  if c + g > 0.75 {
    return "S"
  }
  if c + t > 0.75 {
    return "Y"
  }
  if g + t > 0.75 {
    return "K"
  }
  // three-letter (each > 0.25)
  if a > 0.25 && c > 0.25 && t > 0.25 {
    return "H"
  }
  if c > 0.25 && g > 0.25 && t > 0.25 {
    return "B"
  }
  if a > 0.25 && c > 0.25 && g > 0.25 {
    return "V"
  }
  if a > 0.25 && g > 0.25 && t > 0.25 {
    return "D"
  }
  "N"
}

///|
/// Inverse of `get_consensus` for a single IUPAC letter, returning the fixed
/// PPM column (alphabet order `[A, C, G, T]` for DNA/RNA). `A/C/G/T/U` map to a
/// near-pure column, two-letter codes share `0.499` each, three-letter codes
/// share `0.333` each, `N` maps to uniform `0.25`.
pub fn consensus_letter_to_ppm(letter : String) -> Array[Double] {
  let l = if letter.length() == 1 {
    letter
  } else {
    letter[0:1].to_owned()
  }
  let pure = 0.001
  match l {
    "A" => [0.997, pure, pure, pure]
    "C" => [pure, 0.997, pure, pure]
    "G" => [pure, pure, 0.997, pure]
    "T" => [pure, pure, pure, 0.997]
    "U" => [pure, pure, pure, 0.997]
    "M" => [0.499, 0.499, pure, pure]
    "R" => [0.499, pure, 0.499, pure]
    "W" => [0.499, pure, pure, 0.499]
    "S" => [pure, 0.499, 0.499, pure]
    "Y" => [pure, 0.499, pure, 0.499]
    "K" => [pure, pure, 0.499, 0.499]
    "H" => [0.333, 0.333, pure, 0.333]
    "B" => [pure, 0.333, 0.333, 0.333]
    "V" => [0.333, 0.333, 0.333, pure]
    "D" => [0.333, pure, 0.333, 0.333]
    _ => [0.25, 0.25, 0.25, 0.25]
  }
}

///|
/// Build the consensus string for an entire PPM motif.
pub fn consensus_string(ppm : Array[Array[Double]], alphabet : String) -> String {
  let mut s = ""
  for pos = 0; pos < ppm.length(); pos = pos + 1 {
    s = s + get_consensus(ppm[pos], alphabet)
  }
  s
}

// ---------------------------------------------------------------------------
// Motif-level type conversion (return new UniversalMotif)
// ---------------------------------------------------------------------------

///|
/// Convert a motif to PPM form (recomputes IC and consensus).
pub fn UniversalMotif::to_ppm(self : UniversalMotif) -> UniversalMotif {
  let ppm = um_to_ppm(
    self.motif, self.motif_type, self.pseudocount, self.nsites, self.bkg,
  )
  UniversalMotif::new(
    self.name,
    self.alphabet,
    ppm,
    PPM,
    nsites=self.nsites,
    pseudocount=self.pseudocount,
    bkg=self.bkg,
    strand=self.strand,
    altname=self.altname,
    family=self.family,
    organism=self.organism,
    bkgsites=self.bkgsites,
    extrainfo=self.extrainfo,
  )
}

///|
/// Convert a motif to PWM form.
pub fn UniversalMotif::to_pwm(self : UniversalMotif) -> UniversalMotif {
  let ppm = um_to_ppm(
    self.motif, self.motif_type, self.pseudocount, self.nsites, self.bkg,
  )
  let pwm = ppm_to_pwm(ppm, self.bkg, self.pseudocount, self.nsites)
  UniversalMotif::new(
    self.name,
    self.alphabet,
    pwm,
    PWM,
    nsites=self.nsites,
    pseudocount=self.pseudocount,
    bkg=self.bkg,
    strand=self.strand,
    altname=self.altname,
    family=self.family,
    organism=self.organism,
    bkgsites=self.bkgsites,
    extrainfo=self.extrainfo,
  )
}

///|
/// Convert a motif to ICM form.
pub fn UniversalMotif::to_icm(
  self : UniversalMotif,
  relative_entropy? : Bool = false,
) -> UniversalMotif {
  let ppm = um_to_ppm(
    self.motif, self.motif_type, self.pseudocount, self.nsites, self.bkg,
  )
  let icm = ppm_to_icm(ppm, self.bkg, relative_entropy)
  UniversalMotif::new(
    self.name,
    self.alphabet,
    icm,
    ICM,
    nsites=self.nsites,
    pseudocount=self.pseudocount,
    bkg=self.bkg,
    strand=self.strand,
    altname=self.altname,
    family=self.family,
    organism=self.organism,
    bkgsites=self.bkgsites,
    extrainfo=self.extrainfo,
  )
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

///|
/// Compare two single PPM columns (positions) with the given metric.
/// `bkg1`/`bkg2` default to uniform; `nsites1`/`nsites2` default to 100.
/// Zero-handling for KL/ALLR/ALLR_LL adds 0.01 to every cell.
pub fn compare_positions(
  metric : UMetric,
  row1 : Array[Double],
  row2 : Array[Double],
  bkg1 : Array[Double],
  bkg2 : Array[Double],
  nsites1 : Double,
  nsites2 : Double,
) -> Double {
  let n = if row1.length() < row2.length() {
    row1.length()
  } else {
    row2.length()
  }
  if n == 0 {
    return 0.0
  }
  let bg1 = if bkg1.length() == n { bkg1 } else { um_default_bkg(n) }
  let bg2 = if bkg2.length() == n { bkg2 } else { um_default_bkg(n) }
  let ns1 = if nsites1 <= 0.0 { 100.0 } else { nsites1 }
  let ns2 = if nsites2 <= 0.0 { 100.0 } else { nsites2 }
  let needs_zero_fix = match metric {
    KL | ALLR | ALLR_LL => true
    _ => false
  }
  let (x, y, b1, b2) = if needs_zero_fix {
    let x2 = Array::make(n, 0.0)
    let y2 = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      x2[j] = row1[j] + 0.01
      y2[j] = row2[j] + 0.01
    }
    let mut bkg_zero = false
    for j = 0; j < n; j = j + 1 {
      if bg1[j] == 0.0 || bg2[j] == 0.0 {
        bkg_zero = true
      }
    }
    let add = 0.01 / n.to_double()
    let bb1 = Array::make(n, 0.0)
    let bb2 = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      bb1[j] = bg1[j] + (if bkg_zero { add } else { 0.0 })
      bb2[j] = bg2[j] + (if bkg_zero { add } else { 0.0 })
    }
    (x2, y2, bb1, bb2)
  } else {
    (row1, row2, bg1, bg2)
  }
  match metric {
    PCC => {
      let mut sx = 0.0
      let mut sy = 0.0
      let mut sxy = 0.0
      let mut sxx = 0.0
      let mut syy = 0.0
      for j = 0; j < n; j = j + 1 {
        sx = sx + x[j]
        sy = sy + y[j]
        sxy = sxy + x[j] * y[j]
        sxx = sxx + x[j] * x[j]
        syy = syy + y[j] * y[j]
      }
      let nd = n.to_double()
      let denom = (nd * sxx - sx * sx) * (nd * syy - sy * sy)
      if denom <= 0.0 {
        return 1.0
      }
      (nd * sxy - sx * sy) / um_sqrt(denom)
    }
    WPCC => {
      let w = Array::make(n, 0.0)
      let mut sw = 0.0
      for j = 0; j < n; j = j + 1 {
        w[j] = b1[j] + b2[j]
        sw = sw + w[j]
      }
      if sw <= 0.0 {
        return 1.0
      }
      let mut mx = 0.0
      let mut my = 0.0
      for j = 0; j < n; j = j + 1 {
        mx = mx + w[j] * x[j]
        my = my + w[j] * y[j]
      }
      mx = mx / sw
      my = my / sw
      let mut cov = 0.0
      let mut vx = 0.0
      let mut vy = 0.0
      for j = 0; j < n; j = j + 1 {
        let dx = x[j] - mx
        let dy = y[j] - my
        cov = cov + w[j] * dx * dy
        vx = vx + w[j] * dx * dx
        vy = vy + w[j] * dy * dy
      }
      if vx <= 0.0 || vy <= 0.0 {
        return 1.0
      }
      cov / um_sqrt(vx * vy)
    }
    EUCL => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        let d = x[j] - y[j]
        s = s + d * d
      }
      um_sqrt(s)
    }
    SEUCL => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        let d = x[j] - y[j]
        s = s + d * d
      }
      s
    }
    MAN => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        s = s + (x[j] - y[j]).abs()
      }
      s
    }
    KL => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        s = s + 0.5 * (x[j] * um_ln(x[j] / y[j]) + y[j] * um_ln(y[j] / x[j]))
      }
      s
    }
    HELL => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        let d = um_sqrt(x[j]) - um_sqrt(y[j])
        s = s + d * d
      }
      um_sqrt(s) / um_sqrt(2.0)
    }
    BHAT => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        s = s + um_sqrt(x[j] * y[j])
      }
      s
    }
    SW => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        let d = x[j] - y[j]
        s = s + d * d
      }
      2.0 - s
    }
    ALLR => {
      let mut t1 = 0.0
      let mut t2 = 0.0
      for j = 0; j < n; j = j + 1 {
        t1 = t1 + y[j] * um_ln(x[j] / b1[j])
        t2 = t2 + x[j] * um_ln(y[j] / b2[j])
      }
      (ns2 * t1 + ns1 * t2) / (ns1 + ns2)
    }
    ALLR_LL => {
      let mut t1 = 0.0
      let mut t2 = 0.0
      for j = 0; j < n; j = j + 1 {
        t1 = t1 + y[j] * um_ln(x[j] / b1[j])
        t2 = t2 + x[j] * um_ln(y[j] / b2[j])
      }
      let ans = (ns2 * t1 + ns1 * t2) / (ns1 + ns2)
      if ans < -2.0 {
        -2.0
      } else {
        ans
      }
    }
    WEUCL => {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        let d = x[j] - y[j]
        s = s + d * d * (b1[j] + b2[j])
      }
      um_sqrt(s)
    }
  }
}

///|
/// Combine per-position scores according to a score strategy.
fn combine_scores(scores : Array[Double], strat : UScoreStrat) -> Double {
  if scores.length() == 0 {
    return 0.0
  }
  match strat {
    USum => {
      let mut s = 0.0
      for i = 0; i < scores.length(); i = i + 1 {
        s = s + scores[i]
      }
      s
    }
    UAMean | UWAMean => {
      let mut s = 0.0
      for i = 0; i < scores.length(); i = i + 1 {
        s = s + scores[i]
      }
      s / scores.length().to_double()
    }
    UGMean | UWGMean => {
      let mut sum_log = 0.0
      let mut cnt = 0
      for i = 0; i < scores.length(); i = i + 1 {
        if scores[i] > 0.0 {
          sum_log = sum_log + @math.ln(scores[i])
          cnt = cnt + 1
        }
      }
      if cnt == 0 {
        return 0.0
      }
      @math.exp(sum_log / cnt.to_double())
    }
    UMedian => {
      let sorted = scores.copy()
      sorted.sort_by(fn(a : Double, b : Double) -> Int {
        if a < b {
          -1
        } else if a > b {
          1
        } else {
          0
        }
      })
      let m = sorted.length()
      if m % 2 == 1 {
        sorted[m / 2]
      } else {
        (sorted[m / 2 - 1] + sorted[m / 2]) / 2.0
      }
    }
    UFZT => {
      let mut s = 0.0
      for i = 0; i < scores.length(); i = i + 1 {
        let v = scores[i]
        let clamped = if v >= 0.9999 {
          0.9999
        } else if v <= -0.9999 {
          -0.9999
        } else {
          v
        }
        s = s + 0.5 * @math.ln((1.0 + clamped) / (1.0 - clamped))
      }
      let mean = s / scores.length().to_double()
      (@math.exp(2.0 * mean) - 1.0) / (@math.exp(2.0 * mean) + 1.0)
    }
  }
}

///|
/// Reverse-complement of a DNA/RNA motif (rows are reversed and the letter
/// columns are complemented: A<->T/U, C<->G). For other alphabets the motif
/// is returned unchanged.
pub fn reverse_complement_motif(motif : UniversalMotif) -> UniversalMotif {
  let alphabet = motif.alphabet
  if alphabet != "DNA" && alphabet != "RNA" {
    return motif
  }
  let ppm = um_to_ppm(
    motif.motif, motif.motif_type, motif.pseudocount, motif.nsites, motif.bkg,
  )
  let len = ppm.length()
  let n = 4
  let rc = um_make_matrix(len, n, 0.0)
  for pos = 0; pos < len; pos = pos + 1 {
    let src = ppm[len - 1 - pos]
    rc[pos][0] = src[3]
    rc[pos][1] = src[2]
    rc[pos][2] = src[1]
    rc[pos][3] = src[0]
  }
  UniversalMotif::new(
    motif.name + "_rc",
    alphabet,
    rc,
    motif.motif_type,
    nsites=motif.nsites,
    pseudocount=motif.pseudocount,
    bkg=motif.bkg,
    strand=motif.strand,
    altname=motif.altname,
    family=motif.family,
    organism=motif.organism,
    bkgsites=motif.bkgsites,
    extrainfo=motif.extrainfo,
  )
}

///|
/// Compare two motifs by sliding one across the other and scoring every
/// alignment offset. The best offset is selected (max for similarity metrics,
/// min for distance metrics). Optionally tries the reverse complement of the
/// second motif. Scores are normalized by alignment length when `normalise` is
/// true (distances scaled by `tlen/alignlen`, similarities by `alignlen/tlen`).
/// Offsets whose overlapping mean information content (either side) falls below
/// `min_mean_ic` are rejected.
pub fn compare_motifs(
  m1 : UniversalMotif,
  m2 : UniversalMotif,
  metric : UMetric,
  min_overlap? : Int = 6,
  min_mean_ic? : Double = 0.25,
  score_strat? : UScoreStrat = UScoreStrat::USum,
  normalise? : Bool = true,
  try_rc? : Bool = true,
  relative_entropy? : Bool = false,
) -> Double {
  let ppm1 = um_to_ppm(
    m1.motif, m1.motif_type, m1.pseudocount, m1.nsites, m1.bkg,
  )
  let ppm2 = um_to_ppm(
    m2.motif, m2.motif_type, m2.pseudocount, m2.nsites, m2.bkg,
  )
  let best = um_score_alignment(
    ppm1, ppm2, m1.bkg, m2.bkg, m1.nsites, m2.nsites, metric, min_overlap,
    min_mean_ic, score_strat, normalise, relative_entropy,
  )
  if try_rc && (m1.alphabet == "DNA" || m1.alphabet == "RNA") {
    let m2rc = reverse_complement_motif(m2)
    let ppm2rc = um_to_ppm(
      m2rc.motif, m2rc.motif_type, m2rc.pseudocount, m2rc.nsites, m2rc.bkg,
    )
    let best_rc = um_score_alignment(
      ppm1, ppm2rc, m1.bkg, m2rc.bkg, m1.nsites, m2rc.nsites, metric,
      min_overlap, min_mean_ic, score_strat, normalise, relative_entropy,
    )
    return um_pick_best(metric, best, best_rc)
  }
  best
}

/// Score every alignment offset between two PPMs and return the best.
fn um_score_alignment(
  ppm1 : Array[Array[Double]],
  ppm2 : Array[Array[Double]],
  bkg1 : Array[Double],
  bkg2 : Array[Double],
  nsites1 : Double,
  nsites2 : Double,
  metric : UMetric,
  min_overlap : Int,
  min_mean_ic : Double,
  score_strat : UScoreStrat,
  normalise : Bool,
  relative_entropy : Bool,
) -> Double {
  let l1 = ppm1.length()
  let l2 = ppm2.length()
  if l1 == 0 || l2 == 0 {
    return if um_is_similarity(metric) { 0.0 } else { 1.0e18 }
  }
  let tlen = if l1 > l2 { l1 } else { l2 }
  let is_sim = um_is_similarity(metric)
  let mut best = if is_sim { -1.0e18 } else { 1.0e18 }
  let mut found = false
  let mut offset = -(l2 - 1)
  while offset <= l1 - 1 {
    let start1 = if offset > 0 { offset } else { 0 }
    let end1 = if offset + l2 < l1 { offset + l2 } else { l1 }
    let overlap_len = end1 - start1
    if overlap_len >= min_overlap && overlap_len > 0 {
      let pos_scores = Array::new()
      let mut ic1_sum = 0.0
      let mut ic2_sum = 0.0
      let mut i = start1
      while i < end1 {
        let j = i - offset
        pos_scores.push(
          compare_positions(metric, ppm1[i], ppm2[j], bkg1, bkg2, nsites1, nsites2),
        )
        ic1_sum = ic1_sum + position_icscore_ppm(ppm1[i], bkg1, relative_entropy)
        ic2_sum = ic2_sum + position_icscore_ppm(ppm2[j], bkg2, relative_entropy)
        i = i + 1
      }
      let mean_ic1 = ic1_sum / overlap_len.to_double()
      let mean_ic2 = ic2_sum / overlap_len.to_double()
      if mean_ic1 >= min_mean_ic && mean_ic2 >= min_mean_ic {
        let combined = combine_scores(pos_scores, score_strat)
        let alignlen = overlap_len.to_double()
        let tlend = tlen.to_double()
        let normalized = if normalise {
          if is_sim {
            combined * alignlen / tlend
          } else {
            combined * tlend / alignlen
          }
        } else {
          combined
        }
        if !found {
          best = normalized
          found = true
        } else if is_sim && normalized > best {
          best = normalized
        } else if !is_sim && normalized < best {
          best = normalized
        }
      }
    }
    offset = offset + 1
  }
  if !found {
    if is_sim {
      0.0
    } else {
      1.0e18
    }
  } else {
    best
  }
}

fn um_pick_best(metric : UMetric, a : Double, b : Double) -> Double {
  if um_is_similarity(metric) {
    if a > b {
      a
    } else {
      b
    }
  } else if a < b {
    a
  } else {
    b
  }
}

// ---------------------------------------------------------------------------
// Merge
// ---------------------------------------------------------------------------

///|
/// Merge a list of motifs into a single motif using a weighted arithmetic
/// mean of PPM columns (matching `merge_motifs`). All inputs are converted to
/// PPM; the best pairwise alignment (by `metric`) is found and the overlapping
/// columns are averaged with the running accumulator weighted by how many
/// motifs it has already absorbed; non-overlapping columns are taken verbatim.
/// `nsites`/`bkgsites` become the max across inputs, `pseudocount` the mean;
/// `name` is inputs joined with `/`.
pub fn merge_motifs(
  motifs : Array[UniversalMotif],
  metric? : UMetric = UMetric::ALLR,
  min_overlap? : Int = 6,
  min_mean_ic? : Double = 0.25,
  try_rc? : Bool = true,
) -> UniversalMotif {
  if motifs.length() == 0 {
    return UniversalMotif::new("", "DNA", [], PCM)
  }
  if motifs.length() == 1 {
    return motifs[0].to_ppm()
  }
  let mut acc_ppm = um_to_ppm(
    motifs[0].motif,
    motifs[0].motif_type,
    motifs[0].pseudocount,
    motifs[0].nsites,
    motifs[0].bkg,
  )
  let acc_bkg = motifs[0].bkg.copy()
  let mut acc_name = motifs[0].name
  let mut acc_nsites = motifs[0].nsites
  let mut acc_bkgsites = motifs[0].bkgsites
  let mut acc_pseudocount = motifs[0].pseudocount
  let acc_alphabet = motifs[0].alphabet
  let acc_strand = motifs[0].strand
  let mut weight = 1
  for k = 1; k < motifs.length(); k = k + 1 {
    let next = motifs[k]
    let next_ppm = um_to_ppm(
      next.motif, next.motif_type, next.pseudocount, next.nsites, next.bkg,
    )
    let merged = um_merge_pair(
      acc_ppm,
      next_ppm,
      acc_bkg,
      next.bkg,
      weight,
      metric,
      min_overlap,
      min_mean_ic,
      try_rc,
      acc_alphabet,
    )
    acc_ppm = merged
    acc_name = acc_name + "/" + next.name
    if next.nsites > acc_nsites {
      acc_nsites = next.nsites
    }
    if next.bkgsites > acc_bkgsites {
      acc_bkgsites = next.bkgsites
    }
    acc_pseudocount = (acc_pseudocount + next.pseudocount) / 2.0
    weight = weight + 1
  }
  // Update the merged background (um_merge_pair returns merged PPM only).
  UniversalMotif::new(
    acc_name,
    acc_alphabet,
    acc_ppm,
    PPM,
    nsites=acc_nsites,
    pseudocount=acc_pseudocount,
    bkg=acc_bkg,
    strand=acc_strand,
  )
}

/// Merge a pair of PPMs with a weighted arithmetic mean at the best offset.
/// Returns the merged PPM (background is merged in-place by the caller via
/// `acc_bkg`, so this helper keeps the caller's `acc_bkg`).
fn um_merge_pair(
  ppm1 : Array[Array[Double]],
  ppm2 : Array[Array[Double]],
  bkg1 : Array[Double],
  bkg2 : Array[Double],
  weight : Int,
  metric : UMetric,
  min_overlap : Int,
  min_mean_ic : Double,
  try_rc : Bool,
  alphabet : String,
) -> Array[Array[Double]] {
  let n = if ppm1.length() == 0 {
    if ppm2.length() == 0 { 0 } else { ppm2[0].length() }
  } else {
    ppm1[0].length()
  }
  let best_offset = um_find_best_offset(
    ppm1, ppm2, bkg1, bkg2, 100.0, 100.0, metric, min_overlap, min_mean_ic,
    try_rc, alphabet,
  )
  let l1 = ppm1.length()
  let ppm2_use = if best_offset.rc { um_rc_ppm(ppm2) } else { ppm2 }
  let l2 = ppm2_use.length()
  let left = if best_offset.offset > 0 { best_offset.offset } else { 0 }
  let right1 = best_offset.offset + l2
  let right = if right1 > l1 { right1 } else { l1 }
  let merged_width = right - left
  if merged_width <= 0 || n == 0 {
    return ppm1
  }
  let merged = um_make_matrix(merged_width, n, 0.0)
  for pos = 0; pos < merged_width; pos = pos + 1 {
    let i1 = best_offset.offset + pos // index into ppm1
    let i2 = pos // index into ppm2 (col 0 aligns to offset)
    let in1 = i1 >= 0 && i1 < l1
    let in2 = i2 >= 0 && i2 < l2
    for j = 0; j < n; j = j + 1 {
      if in1 && in2 {
        merged[pos][j] = (ppm1[i1][j] * weight.to_double() +
          ppm2_use[i2][j]) /
          (weight.to_double() + 1.0)
      } else if in1 {
        merged[pos][j] = ppm1[i1][j]
      } else if in2 {
        merged[pos][j] = ppm2_use[i2][j]
      } else {
        merged[pos][j] = 0.0
      }
    }
  }
  // Merge background in-place is the caller's job; here we leave bkg alone.
  // (The caller keeps `acc_bkg` unchanged; background merging is done below.)
  let nb = if bkg1.length() < bkg2.length() {
    bkg1.length()
  } else {
    bkg2.length()
  }
  let mut bi = 0
  while bi < nb {
    acc_bkg_merge(bkg1, bkg2, bi, weight)
    bi = bi + 1
  }
  merged
}

// In-place weighted background merge helper.
fn acc_bkg_merge(
  bkg1 : Array[Double],
  bkg2 : Array[Double],
  idx : Int,
  weight : Int,
) -> Unit {
  bkg1[idx] = (bkg1[idx] * weight.to_double() + bkg2[idx]) /
    (weight.to_double() + 1.0)
}

/// Reverse-complement of a PPM (DNA/RNA, row=position col=letter).
fn um_rc_ppm(ppm : Array[Array[Double]]) -> Array[Array[Double]] {
  let len = ppm.length()
  if len == 0 {
    return []
  }
  let n = ppm[0].length()
  let rc = um_make_matrix(len, n, 0.0)
  for pos = 0; pos < len; pos = pos + 1 {
    let src = ppm[len - 1 - pos]
    if n >= 4 {
      rc[pos][0] = src[3]
      rc[pos][1] = src[2]
      rc[pos][2] = src[1]
      rc[pos][3] = src[0]
    } else {
      for j = 0; j < n; j = j + 1 {
        rc[pos][j] = src[n - 1 - j]
      }
    }
  }
  rc
}

struct UmOffset {
  offset : Int
  rc : Bool
  score : Double
}

/// Find the alignment offset that yields the best comparison score.
fn um_find_best_offset(
  ppm1 : Array[Array[Double]],
  ppm2 : Array[Array[Double]],
  bkg1 : Array[Double],
  bkg2 : Array[Double],
  nsites1 : Double,
  nsites2 : Double,
  metric : UMetric,
  min_overlap : Int,
  min_mean_ic : Double,
  try_rc : Bool,
  alphabet : String,
) -> UmOffset {
  let is_sim = um_is_similarity(metric)
  let mut best = UmOffset::{ offset: 0, rc: false, score: 0.0 }
  let mut found = false
  let candidates : Array[Array[Array[Double]]] = [ppm2]
  let rc_flags : Array[Bool] = [false]
  if try_rc && (alphabet == "DNA" || alphabet == "RNA") {
    candidates.push(um_rc_ppm(ppm2))
    rc_flags.push(true)
  }
  for c = 0; c < candidates.length(); c = c + 1 {
    let ppm2c = candidates[c]
    let rc = rc_flags[c]
    let l1 = ppm1.length()
    let l2 = ppm2c.length()
    let tlen = if l1 > l2 { l1 } else { l2 }
    let mut offset = -(l2 - 1)
    while offset <= l1 - 1 {
      let start1 = if offset > 0 { offset } else { 0 }
      let end1 = if offset + l2 < l1 { offset + l2 } else { l1 }
      let overlap_len = end1 - start1
      if overlap_len >= min_overlap && overlap_len > 0 {
        let pos_scores = Array::new()
        let mut ic1_sum = 0.0
        let mut ic2_sum = 0.0
        let mut i = start1
        while i < end1 {
          let j = i - offset
          pos_scores.push(
            compare_positions(
              metric, ppm1[i], ppm2c[j], bkg1, bkg2, nsites1, nsites2,
            ),
          )
          ic1_sum = ic1_sum + position_icscore_ppm(ppm1[i], bkg1, false)
          ic2_sum = ic2_sum + position_icscore_ppm(ppm2c[j], bkg2, false)
          i = i + 1
        }
        let mean_ic1 = ic1_sum / overlap_len.to_double()
        let mean_ic2 = ic2_sum / overlap_len.to_double()
        if mean_ic1 >= min_mean_ic && mean_ic2 >= min_mean_ic {
          let combined = combine_scores(pos_scores, UScoreStrat::USum)
          let alignlen = overlap_len.to_double()
          let tlend = tlen.to_double()
          let normalized = if is_sim {
            combined * alignlen / tlend
          } else {
            combined * tlend / alignlen
          }
          if !found {
            best = UmOffset::{ offset, rc, score: normalized }
            found = true
          } else if is_sim && normalized > best.score {
            best = UmOffset::{ offset, rc, score: normalized }
          } else if !is_sim && normalized < best.score {
            best = UmOffset::{ offset, rc, score: normalized }
          }
        }
      }
      offset = offset + 1
    }
  }
  best
}

// ---------------------------------------------------------------------------
// Trim
// ---------------------------------------------------------------------------

///|
/// Trim low-information-content columns from the flanks of a motif (matching
/// `trim_motifs`). Per-column IC is computed from the motif's own type/bkg;
/// columns with IC below `min_ic` are removed from the edges only (`trim_from`
/// controls left/right/both). `icscore` and `consensus` are recomputed.
pub fn trim_motif(
  motif : UniversalMotif,
  min_ic? : Double = 0.25,
  trim_from? : UTrimFrom = UTrimFrom::UBoth,
) -> UniversalMotif {
  let ppm = um_to_ppm(
    motif.motif, motif.motif_type, motif.pseudocount, motif.nsites, motif.bkg,
  )
  if ppm.length() == 0 {
    return motif
  }
  let ic_scores = Array::make(ppm.length(), 0.0)
  for pos = 0; pos < ppm.length(); pos = pos + 1 {
    ic_scores[pos] = position_icscore_ppm(ppm[pos], motif.bkg, false)
  }
  let mut cut_left = 0
  let mut cut_right = 0
  let len = ppm.length()
  let trim_left = match trim_from {
    UBoth | ULeft => true
    _ => false
  }
  let trim_right = match trim_from {
    UBoth | URight => true
    _ => false
  }
  if trim_left {
    while cut_left < len && ic_scores[cut_left] < min_ic {
      cut_left = cut_left + 1
    }
  }
  if trim_right {
    while cut_right < len && ic_scores[len - 1 - cut_right] < min_ic {
      cut_right = cut_right + 1
    }
  }
  if cut_left + cut_right >= len {
    return UniversalMotif::new(
      motif.name,
      motif.alphabet,
      [],
      PPM,
      nsites=motif.nsites,
      pseudocount=motif.pseudocount,
      bkg=motif.bkg,
      strand=motif.strand,
    )
  }
  let new_ppm = Array::new()
  for pos = cut_left; pos < len - cut_right; pos = pos + 1 {
    new_ppm.push(ppm[pos].copy())
  }
  UniversalMotif::new(
    motif.name,
    motif.alphabet,
    new_ppm,
    PPM,
    nsites=motif.nsites,
    pseudocount=motif.pseudocount,
    bkg=motif.bkg,
    strand=motif.strand,
    altname=motif.altname,
    family=motif.family,
    organism=motif.organism,
    bkgsites=motif.bkgsites,
    extrainfo=motif.extrainfo,
  )
}

// ---------------------------------------------------------------------------
// MEME / JASPAR I/O
// ---------------------------------------------------------------------------

/// Format a probability to a fixed number of decimal places (round half-up).
fn um_format_fixed(d : Double, decimals : Int) -> String {
  let neg = d < 0.0
  let abs_d = if neg { -d } else { d }
  let scale = @math.pow(10.0, decimals.to_double())
  let scaled = (abs_d * scale + 0.5).floor().to_int()
  let mut pow10 = 1
  let mut k = 0
  while k < decimals {
    pow10 = pow10 * 10
    k = k + 1
  }
  let whole = scaled / pow10
  let frac = scaled % pow10
  let mut frac_str = frac.to_string()
  while frac_str.length() < decimals {
    frac_str = "0" + frac_str
  }
  (if neg { "-" } else { "" }) + whole.to_string() + "." + frac_str
}

///|
/// Render motifs as a MEME-format string. Motifs are written as PPMs (the
/// universalmotif default). The alphabet line, strands and background letter
/// frequencies are emitted.
pub fn write_meme(motifs : Array[UniversalMotif]) -> String {
  let letters = if motifs.length() > 0 {
    um_alphabet_letters(motifs[0].alphabet)
  } else {
    um_alphabet_letters("DNA")
  }
  let alpha_str = letters.join("")
  let mut s = "MEME version 5\n\n"
  s = s + "ALPHABET= " + alpha_str + "\n\n"
  s = s + "strands: + -\n\n"
  s = s + "Background letter frequencies\n"
  for i = 0; i < letters.length(); i = i + 1 {
    if i > 0 {
      s = s + " "
    }
    let p = if motifs.length() > 0 {
      let b = motifs[0].bkg
      if i < b.length() { b[i] } else { 1.0 / letters.length().to_double() }
    } else {
      1.0 / letters.length().to_double()
    }
    s = s + letters[i] + " " + um_format_fixed(p, 6)
  }
  s = s + "\n"
  for k = 0; k < motifs.length(); k = k + 1 {
    let m = motifs[k].to_ppm()
    let altname = m.altname
    s = s + "\nMOTIF " + m.name
    if altname.length() > 0 {
      s = s + " " + altname
    }
    s = s + "\n"
    let w = m.motif.length()
    let al = letters.length()
    let ns = if m.nsites > 0.0 { m.nsites } else { 100.0 }
    s = s + "letter-probability matrix: alength= " + al.to_string() + " w= " + w.to_string() + " nsites= " + um_format_fixed(ns, 3) + " E= 0\n"
    for pos = 0; pos < w; pos = pos + 1 {
      let row = m.motif[pos]
      for j = 0; j < al; j = j + 1 {
        if j > 0 {
          s = s + "  "
        }
        s = s + um_format_fixed(row[j], 6)
      }
      s = s + "\n"
    }
  }
  s
}

///|
/// Parse a MEME-format string into `UniversalMotif` objects (PPM type). The
/// parser recognizes the `MEME version`, `ALPHABET=`, `strands:`,
/// `Background letter frequencies` and `MOTIF` / `letter-probability matrix:`
/// blocks.
pub fn read_meme(text : String) -> Array[UniversalMotif] {
  let motifs = Array::new()
  let lines = text.split("\n").to_array()
  let n = lines.length()
  let mut alphabet = "DNA"
  let mut i = 0
  while i < n {
    let line = lines[i].trim().to_owned()
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if line.has_prefix("ALPHABET=") {
      let alpha = line.substring(start=9).trim().to_owned()
      alphabet = um_infer_alphabet(alpha)
      i = i + 1
      continue
    }
    if line.has_prefix("MOTIF ") {
      let header = line.substring(start=6).trim().to_owned()
      let parts = header.split(" ").to_array()
      let name = if parts.length() >= 1 {
        parts[0].trim().to_owned()
      } else {
        ""
      }
      let altname = if parts.length() >= 2 {
        parts[1].trim().to_owned()
      } else {
        ""
      }
      i = i + 1
      let mut alen = 4
      let mut width = 0
      let mut nsites = 100.0
      while i < n {
        let ml = lines[i].trim().to_owned()
        if ml.has_prefix("letter-probability matrix:") {
          let fields = um_parse_lpm_header(ml)
          alen = fields.alength
          width = fields.width
          nsites = fields.nsites
          i = i + 1
          break
        }
        if ml.length() == 0 || ml.has_prefix("MOTIF") {
          break
        }
        i = i + 1
      }
      if width == 0 {
        continue
      }
      let matrix : Array[Array[Double]] = Array::new()
      let mut read = 0
      while i < n && read < width {
        let ml = lines[i].trim().to_owned()
        if ml.length() == 0 {
          i = i + 1
          continue
        }
        if ml.has_prefix("MOTIF") {
          break
        }
        let parts = ml.split(" ").to_array()
        let row = Array::new()
        for p = 0; p < parts.length(); p = p + 1 {
          let t = parts[p].trim().to_owned()
          if t.length() > 0 {
            match parse_double(t) {
              Some(v) => row.push(v)
              None => ()
            }
          }
        }
        if row.length() >= alen {
          matrix.push(row)
          read = read + 1
        }
        i = i + 1
      }
      if matrix.length() > 0 {
        let m = UniversalMotif::new(
          name,
          alphabet,
          matrix,
          PPM,
          nsites=nsites,
          pseudocount=0.8,
          altname=altname,
        )
        motifs.push(m)
      }
      continue
    }
    i = i + 1
  }
  motifs
}

struct UmLpmHeader {
  alength : Int
  width : Int
  nsites : Double
}

fn um_parse_lpm_header(line : String) -> UmLpmHeader {
  let mut alen = 4
  let mut width = 0
  let mut nsites = 100.0
  let parts = line.split(" ").to_array()
  // Collect non-empty tokens to handle variable spacing
  let tokens = Array::new()
  for p = 0; p < parts.length(); p = p + 1 {
    let t = parts[p].trim().to_owned()
    if t.length() > 0 {
      tokens.push(t)
    }
  }
  for p = 0; p < tokens.length(); p = p + 1 {
    let t = tokens[p]
    if t.has_prefix("alength=") {
      let rest = t.substring(start=8).trim().to_owned()
      if rest.length() > 0 {
        alen = parse_int(rest)
      } else if p + 1 < tokens.length() {
        alen = parse_int(tokens[p + 1])
      }
    } else if t.has_prefix("w=") {
      let rest = t.substring(start=2).trim().to_owned()
      if rest.length() > 0 {
        width = parse_int(rest)
      } else if p + 1 < tokens.length() {
        width = parse_int(tokens[p + 1])
      }
    } else if t.has_prefix("nsites=") {
      let rest = t.substring(start=7).trim().to_owned()
      if rest.length() > 0 {
        match parse_double(rest) {
          Some(x) => nsites = x
          None => ()
        }
      } else if p + 1 < tokens.length() {
        match parse_double(tokens[p + 1]) {
          Some(x) => nsites = x
          None => ()
        }
      }
    }
  }
  UmLpmHeader::{ alength: alen, width, nsites }
}

/// Infer DNA/RNA/AA alphabet from an alphabet string.
fn um_infer_alphabet(alpha : String) -> String {
  let a = alpha.trim().to_owned()
  if a == "ACGT" {
    return "DNA"
  }
  if a == "ACGU" {
    return "RNA"
  }
  let aa = "ACDEFGHIKLMNPQRSTVWY"
  if a.length() >= 20 {
    let prefix = a.substring(start=0, end=20)
    if aa.contains(prefix) {
      return "AA"
    }
  }
  "DNA"
}

///|
/// Render a single motif as a JASPAR-format string (PCM counts, one row per
/// letter). The motif is converted to PCM with `nsites` columns.
pub fn write_jaspar(motif : UniversalMotif) -> String {
  let ppm = um_to_ppm(
    motif.motif, motif.motif_type, motif.pseudocount, motif.nsites, motif.bkg,
  )
  let letters = um_alphabet_letters(motif.alphabet)
  let n = letters.length()
  let width = ppm.length()
  let ns = if motif.nsites > 0.0 { motif.nsites } else { 100.0 }
  let pcm = ppm_to_pcm(ppm, ns)
  let mut s = ">" + motif.name
  if motif.altname.length() > 0 {
    s = s + " " + motif.altname
  }
  s = s + "\n"
  for j = 0; j < n; j = j + 1 {
    s = s + letters[j] + " ["
    for pos = 0; pos < width; pos = pos + 1 {
      if pos > 0 {
        s = s + " "
      }
      s = s + pcm[pos][j].to_string()
    }
    s = s + " ]\n"
  }
  s
}

///|
/// Parse a JASPAR-format string into `UniversalMotif` objects (PCM type). The
/// header line `>name altname` is followed by one row per letter:
/// `LETTER [ c1 c2 ... cn ]`.
pub fn read_jaspar(text : String) -> Array[UniversalMotif] {
  let motifs = Array::new()
  let lines = text.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  while i < n {
    let line = lines[i].trim().to_owned()
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if line.has_prefix(">") {
      let header = line.substring(start=1).trim().to_owned()
      let parts = header.split(" ").to_array()
      let name = if parts.length() >= 1 {
        parts[0].trim().to_owned()
      } else {
        ""
      }
      let altname = if parts.length() >= 2 {
        parts[1].trim().to_owned()
      } else {
        ""
      }
      let row_labels : Array[String] = Array::new()
      let row_counts : Array[Array[Double]] = Array::new()
      i = i + 1
      while i < n {
        let rl = lines[i].trim().to_owned()
        if rl.length() == 0 {
          i = i + 1
          break
        }
        if rl.has_prefix(">") {
          break
        }
        if rl.contains("[") {
          let after = rl.split("[").to_array()
          let tail = if after.length() >= 2 {
            after[1].to_owned()
          } else {
            ""
          }
          let label_parts = rl.split(" ").to_array()
          row_labels.push(label_parts[0].trim().to_owned())
          let cleaned_parts = tail.split("]").to_array()
          let cleaned = if cleaned_parts.length() >= 1 {
            cleaned_parts[0].to_owned()
          } else {
            ""
          }
          let nums = cleaned.split(" ").to_array()
          let counts = Array::new()
          for p = 0; p < nums.length(); p = p + 1 {
            let t = nums[p].trim().to_owned()
            if t.length() > 0 {
              match parse_double(t) {
                Some(v) => counts.push(v)
                None => ()
              }
            }
          }
          row_counts.push(counts)
        }
        i = i + 1
      }
      if row_counts.length() == 0 {
        continue
      }
      let n_letters = row_counts.length()
      let width = row_counts[0].length()
      let pcm = um_make_matrix(width, n_letters, 0.0)
      for letter = 0; letter < n_letters; letter = letter + 1 {
        for pos = 0; pos < width; pos = pos + 1 {
          if pos < row_counts[letter].length() {
            pcm[pos][letter] = row_counts[letter][pos]
          }
        }
      }
      let alphabet = um_infer_alphabet_from_labels(row_labels)
      let mut total_sites = 0.0
      if width > 0 {
        for letter = 0; letter < n_letters; letter = letter + 1 {
          total_sites = total_sites + pcm[0][letter]
        }
      }
      if total_sites <= 0.0 {
        total_sites = 100.0
      }
      let m = UniversalMotif::new(
        name,
        alphabet,
        pcm,
        PCM,
        nsites=total_sites,
        pseudocount=0.8,
        altname=altname,
      )
      motifs.push(m)
    } else {
      i = i + 1
    }
  }
  motifs
}

/// Infer alphabet from JASPAR row labels (A/C/G/T, A/C/G/U, or AA).
fn um_infer_alphabet_from_labels(labels : Array[String]) -> String {
  let mut joined = ""
  for i = 0; i < labels.length(); i = i + 1 {
    joined = joined + labels[i]
  }
  if joined == "ACGT" {
    return "DNA"
  }
  if joined == "ACGU" {
    return "RNA"
  }
  let aa = "ACDEFGHIKLMNPQRSTVWY"
  if joined.length() >= 20 {
    let prefix = joined.substring(start=0, end=20)
    if aa.contains(prefix) {
      return "AA"
    }
  }
  "DNA"
}

// ---------------------------------------------------------------------------
// Example
// ---------------------------------------------------------------------------

///|
/// A small TATA-box example motif (DNA, PPM).
pub fn create_example_universal_motif() -> UniversalMotif {
  let ppm = [
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
  ]
  UniversalMotif::new("TATA-box", "DNA", ppm, PPM, nsites=20.0)
}
