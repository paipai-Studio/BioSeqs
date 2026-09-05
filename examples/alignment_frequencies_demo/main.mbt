///|
/// alignment_frequencies module example - MSA substitution and column
/// statistics ported from Biopython `Bio.Align`:
/// MultipleSeqAlignment.substitutions (symmetric weighted pair-count
/// matrix) and Alignment.frequencies (per-column letter counts with a
/// '-' gap row), plus observed/expected relative pair frequencies.
///
/// Run: moon run examples/alignment_frequencies_demo
fn main {
  println("=== Biopython Bio.Align MSA Statistics Demo ===")
  println("Substitution counts and column frequencies of multiple alignments")
  println("")

  // ---- Section 1: substitution count matrix (Biopython doc example) ----
  let rows : Array[String] = ["ACGT", "A--A", "ACGT", "TTTC"]
  let sub = @bio.alignment_substitution_matrix(rows) catch {
    _ => abort("substitution matrix should build")
  }
  println("1) Substitution count matrix (MultipleSeqAlignment.substitutions)")
  println("   MSA: ACGT / A--A / ACGT / TTTC")
  println(sub.to_table())
  println("   total pair mass: \{sub.total()}")
  println("")

  // ---- Section 2: per-column counts and consensus --------------------
  let counts = @bio.alignment_letter_counts(rows) catch {
    _ => abort("letter counts should build")
  }
  let freq = counts.column_frequencies()
  let gap_k = match counts.letter_index("-") {
    Some(i) => i
    None => -1
  }
  println("2) Per-column letter counts (Alignment.frequencies; '-' = gaps)")
  println(counts.to_table())
  let mut consensus = ""
  let mut conserved = 0
  let mut j = 0
  while j < counts.ncols {
    let mut best_k = -1
    let mut best_f = 0.0
    let mut k = 0
    while k < counts.letters.length() {
      if k != gap_k && freq[k][j] > best_f {
        best_f = freq[k][j]
        best_k = k
      }
      k = k + 1
    }
    let gap_f = if gap_k >= 0 { freq[gap_k][j] } else { 0.0 }
    if best_k >= 0 {
      consensus += counts.letters[best_k]
    } else {
      consensus += "-"
    }
    if best_f >= 0.75 && gap_f < 0.25 {
      conserved = conserved + 1
    }
    j = j + 1
  }
  println("   consensus sequence : \{consensus}")
  println("   conserved columns  : \{conserved} / \{counts.ncols} (>=75% identity)")
  println("")

  // ---- Section 3: sequence weights ------------------------------------
  let wrows : Array[String] = ["GTATC", "AT--C", "CTGTC"]
  let wweights : Array[Double] = [0.5, 0.8, 1.0]
  let wsub = @bio.alignment_substitution_matrix(wrows, weights=wweights) catch {
    _ => abort("weighted substitution matrix should build")
  }
  println("3) Weighted alignment (sequence weights 0.5 / 0.8 / 1.0)")
  println("   MSA: GTATC / AT--C / CTGTC")
  println(wsub.to_table())
  println("   total weighted pair mass: \{lod_demo_fmt(wsub.total()).trim()} (unweighted would be 9.000)")
  println("")

  // ---- Section 4: observed vs expected pair frequencies --------------
  let bg = sub.background_frequencies()
  let rel = sub.relative_frequencies() catch {
    _ => abort("relative frequencies should build")
  }
  println("4) Log-odds scores log2(q_ij / e_ij) with MLE background")
  let mut header = "        "
  let mut li = 0
  while li < sub.letters.length() {
    header += demo_pad_left(sub.letters[li], 7)
    li = li + 1
  }
  println(header)
  let mut a = 0
  while a < sub.letters.length() {
    let mut line = demo_pad_left(sub.letters[a], 7) + " "
    let mut b = 0
    while b < sub.letters.length() {
      if sub.matrix[a][b] == 0.0 {
        line += "      ."
      } else {
        let lod = @math.log2(rel[a][b])
        line += lod_demo_fmt(lod)
      }
      b = b + 1
    }
    println(line)
    a = a + 1
  }
  let mut bg_line = "   background: "
  li = 0
  while li < sub.letters.length() {
    let pv = match bg.get(sub.letters[li]) {
      Some(v) => v
      None => abort("background should cover all letters")
    }
    bg_line += sub.letters[li] + "=" + lod_demo_fmt(pv) + " "
    li = li + 1
  }
  println(bg_line)
}

///|
/// Left-pad a string with spaces to the requested width.
fn demo_pad_left(text : String, width : Int) -> String {
  let mut out = text
  while out.length() < width {
    out = " " + out
  }
  out
}

///|
/// Fixed-point formatting helper (3 decimals, 7-wide) for the demo.
fn lod_demo_fmt(x : Double) -> String {
  let neg = x < 0.0
  let v = x.abs() + 0.0005
  let whole = v.floor().to_int()
  let frac = ((v - v.floor()) * 1000.0).floor().to_int()
  let frac_s = if frac >= 100 {
    frac.to_string()
  } else if frac >= 10 {
    "0" + frac.to_string()
  } else {
    "00" + frac.to_string()
  }
  let mut s = (if neg { "-" } else { " " }) + whole.to_string() + "." + frac_s
  while s.length() < 7 {
    s = " " + s
  }
  s
}
