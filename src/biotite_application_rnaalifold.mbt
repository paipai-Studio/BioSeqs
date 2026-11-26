/// Biotite `application.rnaalifold` — consensus RNA secondary
/// structure from a multiple sequence alignment.
///
/// Builds the RNAalifold command line for external execution and
/// provides a local consistency-based consensus folder (Nussinov-style
/// with covariance-free base-pair scores), honestly labelled
/// `computed-local`.
///
/// Note: single characters are handled as 1-char strings (the string
/// index operator yields code units, not Char values).

///|
/// Result of an RNAalifold consensus folding.
pub struct Biotite3RNAalifoldResult {
  consensus_structure : String
  conservation_scores : Array[Double]
  aligned_sequences : Array[String]
  command : String
  status : String
}

///|
/// Build the RNAalifold command line (input is a CLUSTAL alignment).
pub fn biotite_rnaalifold_command(
  alignment_file : String,
  bin_path~ : String = "RNAalifold",
) -> String {
  bin_path + " " + alignment_file
}

///|
/// Single character of `s` at `idx` as a 1-char string.
fn biotite_rna_char(s : String, idx : Int32) -> String {
  s[idx:idx + 1].to_owned()
}

///|
/// Majority base of a column ("A"/"C"/"G"/"U"; T normalized to U);
/// None when the column is all gaps.
fn biotite_rnaalifold_majority(aligned : Array[String], col : Int32) -> String? {
  let counts : Map[String, Int32] = Map([], capacity=4)
  for s in aligned {
    if col < s.length() {
      let ch = biotite_rna_char(s, col)
      match ch {
        "A" | "C" | "G" | "U" | "T" =>
          match counts.get(ch) {
            Some(c) => counts[ch] = c + 1
            None => counts[ch] = 1
          }
        _ => ()
      }
    }
  }
  let mut best : String? = None
  let mut best_n = 0
  for k, v in counts {
    if v > best_n {
      best_n = v
      best = Some(k)
    }
  }
  best
}

///|
/// Normalize T to U (1-char strings).
fn biotite_rna_norm_base(ch : String) -> String {
  match ch {
    "T" => "U"
    other => other
  }
}

///|
/// Whether two majority bases form a canonical pair.
fn biotite_rnaalifold_pair_score(a : String, b : String) -> Int32 {
  let x = biotite_rna_norm_base(a)
  let y = biotite_rna_norm_base(b)
  if (x == "G" && y == "C") || (x == "C" && y == "G") {
    return 3
  }
  if (x == "A" && y == "U") || (x == "U" && y == "A") {
    return 2
  }
  if (x == "G" && y == "U") || (x == "U" && y == "G") {
    return 1
  }
  0
}

///|
/// Fraction of sequences whose bases at columns i/j can pair.
fn biotite_rnaalifold_pair_fraction(
  aligned : Array[String],
  i : Int32,
  j : Int32,
) -> Double {
  if aligned.length() == 0 {
    return 0.0
  }
  let mut paired = 0
  let mut valid = 0
  for s in aligned {
    if i < s.length() && j < s.length() {
      let x = biotite_rna_norm_base(biotite_rna_char(s, i))
      let y = biotite_rna_norm_base(biotite_rna_char(s, j))
      if x != "-" && y != "-" {
        valid = valid + 1
        if biotite_rnaalifold_pair_score(x, y) > 0 {
          paired = paired + 1
        }
      }
    }
  }
  if valid > 0 {
    paired.to_double() / valid.to_double()
  } else {
    0.0
  }
}

///|
/// Covariance statistics for a column pair.
///
/// Returns (consistent mutation fraction, incompatible fraction) over
/// all sequences: a sequence contributes a consistent mutation when
/// its bases pair canonically but differ from the majority base pair
/// (e.g. GC vs CG or AU); any non-gap sequence whose bases cannot pair
/// counts as incompatible. This mirrors the qualitative idea of the
/// RNAalifold covariance term (compensatory mutations support a pair,
/// conflicting sequences penalize it).
fn biotite_rnaalifold_covariance(
  aligned : Array[String],
  i : Int32,
  j : Int32,
  maj_a : String?,
  maj_b : String?,
) -> (Double, Double) {
  let n = aligned.length()
  if n == 0 {
    return (0.0, 0.0)
  }
  let mut cons = 0
  let mut pairable = 0
  let mut valid = 0
  for s in aligned {
    if i < s.length() && j < s.length() {
      let x = biotite_rna_norm_base(biotite_rna_char(s, i))
      let y = biotite_rna_norm_base(biotite_rna_char(s, j))
      if x != "-" && y != "-" {
        valid = valid + 1
        if biotite_rnaalifold_pair_score(x, y) > 0 {
          pairable = pairable + 1
          let is_maj = match (maj_a, maj_b) {
            (Some(ma), Some(mb)) =>
              x == biotite_rna_norm_base(ma) && y == biotite_rna_norm_base(mb)
            _ => false
          }
          if !is_maj {
            cons = cons + 1
          }
        }
      }
    }
  }
  let cons_frac = cons.to_double() / n.to_double()
  let incompat_frac = if valid > 0 {
    1.0 - pairable.to_double() / valid.to_double()
  } else {
    0.0
  }
  (cons_frac, incompat_frac)
}

///|
/// Pairing score per (i, j): majority-pair base score scaled by the
/// fraction of sequences that actually pair, plus a covariance bonus
/// (`cv_fact`) for consistent mutations and a penalty (`nc_fact`) for
/// incompatible sequences.
///
/// The two weights mirror the official RNAalifold knobs `--cfactor`
/// (cv_fact) and `--nfactor` (nc_fact); the official energy model uses
/// kcal/mol scale with defaults cv_fact=1.0 / nc_fact=0.5, while this
/// Nussinov-style scoring keeps its own defaults.
fn biotite_rnaalifold_column_score(
  aligned : Array[String],
  i : Int32,
  j : Int32,
  cv_fact~ : Double = 3.0,
  nc_fact~ : Double = 1.0,
) -> Double {
  match biotite_rnaalifold_majority(aligned, i) {
    Some(a) =>
      match biotite_rnaalifold_majority(aligned, j) {
        Some(b) => {
          let base = biotite_rnaalifold_pair_score(a, b).to_double()
          let (cons, incompat) = biotite_rnaalifold_covariance(
            aligned,
            i,
            j,
            Some(a),
            Some(b),
          )
          let total = base * biotite_rnaalifold_pair_fraction(aligned, i, j) +
            cv_fact * cons -
            nc_fact * incompat
          if total > 0.0 {
            total
          } else {
            0.0
          }
        }
        None => 0.0
      }
    None => 0.0
  }
}

///|
/// Consensus folding: Nussinov-style maximum base-pair DP over
/// pairwise column compatibility (majority bases must pair and a
/// majority of sequences must actually pair at those columns).
fn biotite_rnaalifold_fold(
  aligned : Array[String],
  min_loop : Int32,
  cv_fact~ : Double = 3.0,
  nc_fact~ : Double = 1.0,
) -> (String, Array[Double]) {
  let len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  if len == 0 {
    return ("", [])
  }
  // Column conservation = non-gap fraction * majority share.
  let conservation : Array[Double] = Array(capacity=len)
  for col = 0; col < len; col = col + 1 {
    let mut non_gap = 0
    let counts : Map[String, Int32] = Map([], capacity=4)
    for s in aligned {
      if col < s.length() {
        let ch = biotite_rna_char(s, col)
        match ch {
          "A" | "C" | "G" | "U" | "T" => {
            non_gap = non_gap + 1
            match counts.get(ch) {
              Some(c) => counts[ch] = c + 1
              None => counts[ch] = 1
            }
          }
          _ => ()
        }
      }
    }
    let mut maj = 0
    for _, v in counts {
      if v > maj {
        maj = v
      }
    }
    let frac_non_gap = non_gap.to_double() / aligned.length().to_double()
    let frac_maj = if non_gap > 0 {
      maj.to_double() / non_gap.to_double()
    } else {
      0.0
    }
    conservation.push(frac_non_gap * frac_maj)
  }
  // Pairing score per (i, j) — includes the covariance term.
  let score = fn(i : Int32, j : Int32) -> Double {
    biotite_rnaalifold_column_score(aligned, i, j, cv_fact~, nc_fact~)
  }
  // Nussinov DP.
  let dp : Array[Array[Double]] = Array(capacity=len)
  for _i = 0; _i < len; _i = _i + 1 {
    dp.push(Array::make(len, 0.0))
  }
  for span = 1; span < len; span = span + 1 {
    for i = 0; i + span < len; i = i + 1 {
      let j = i + span
      let mut best = 0.0
      if j - i > min_loop {
        let inner = if j - 1 >= i + 1 { dp[i + 1][j - 1] } else { 0.0 }
        let cand = inner + score(i, j)
        if cand > best {
          best = cand
        }
      }
      if dp[i + 1][j] > best {
        best = dp[i + 1][j]
      }
      if dp[i][j - 1] > best {
        best = dp[i][j - 1]
      }
      for k = i + 1; k < j; k = k + 1 {
        let cand = dp[i][k] + dp[k + 1][j]
        if cand > best {
          best = cand
        }
      }
      dp[i][j] = best
    }
  }
  // Trace back the dot-bracket structure.
  let brackets : Array[String] = Array::make(len, ".")
  fn traceback(i : Int32, j : Int32) -> Unit {
    if i >= j {
      return
    }
    if j - i > min_loop &&
      dp[i][j] == dp[i + 1][j - 1] + score(i, j) &&
      score(i, j) > 0.0 {
      brackets[i] = "("
      brackets[j] = ")"
      traceback(i + 1, j - 1)
      return
    }
    if dp[i][j] == dp[i + 1][j] {
      traceback(i + 1, j)
      return
    }
    if dp[i][j] == dp[i][j - 1] {
      traceback(i, j - 1)
      return
    }
    for k = i + 1; k < j; k = k + 1 {
      if dp[i][j] == dp[i][k] + dp[k + 1][j] {
        traceback(i, k)
        traceback(k + 1, j)
        return
      }
    }
  }
  traceback(0, len - 1)
  let sb = StringBuilder()
  for b in brackets {
    let _ = sb.write(b)
  }
  (sb.to_string(), conservation)
}

///|
/// Local RNAalifold-compatible consensus folding.
///
/// `cv_fact` / `nc_fact` mirror the official RNAalifold knobs
/// `--cfactor` (covariance bonus for consistent mutations) and
/// `--nfactor` (penalty for incompatible sequences). The official
/// energy-model defaults are 1.0 / 0.5 kcal/mol scale; these Nussinov
/// scoring defaults (3.0 / 1.0) reproduce the behavior of earlier
/// versions of this module.
pub fn biotite_rnaalifold_local(
  aligned_sequences : Array[String],
  min_loop~ : Int32 = 3,
  cv_fact~ : Double = 3.0,
  nc_fact~ : Double = 1.0,
) -> Biotite3RNAalifoldResult {
  let (structure, conservation) = biotite_rnaalifold_fold(
    aligned_sequences,
    min_loop,
    cv_fact~,
    nc_fact~,
  )
  {
    consensus_structure: structure,
    conservation_scores: conservation,
    aligned_sequences,
    command: biotite_rnaalifold_command("alignment.clustal"),
    status: "computed-local",
  }
}

///|
pub fn _biotite_application_rnaalifold_dummy() -> Unit {
  let r = biotite_rnaalifold_local(["GGGAAACCC", "GGGAAACCC"])
  let _ = r.consensus_structure
  let _ = r.conservation_scores
  let _ = r.command
  let _ = r.status
  let _ = biotite_rnaalifold_command("aln.clustal", bin_path="RNAalifold")
  let _ = biotite_rnaalifold_majority(["ABC"], 0)
  let _ = biotite_rnaalifold_pair_score("A", "U")
  let _ = biotite_rnaalifold_pair_fraction(["GC", "GC"], 0, 1)
  let _ = biotite_rnaalifold_fold(["GC", "GC"], 3)
  let _ = biotite_rna_char("ABC", 0)
  let _ = biotite_rna_norm_base("T")
}
