/// Back-translation: convert a protein sequence back to a nucleotide sequence.
///
/// Biopython equivalent: Bio.SeqUtils.back_translate / Bio.Seq.back_transcribe
///
/// Given a protein sequence, produces a DNA sequence by selecting codons for each
/// amino acid. Three modes:
///   - "first":   pick the first available codon for each AA (deterministic)
///   - "random":  pick a random codon for each AA (uniform)
///   - "biased":  pick codon based on a codon usage table (most frequent codon wins)
///
/// This is the inverse operation of translation. It is inherently degenerate
/// (most amino acids map to multiple codons), so the result is one of many
/// possible back-translations.

///|
/// Back-translation result
pub struct BackTranslateResult {
  protein : String
  nucleotide : String
  table_id : Int32
  mode : String
}

///|
pub fn BackTranslateResult::protein(self : BackTranslateResult) -> String {
  self.protein
}

///|
pub fn BackTranslateResult::nucleotide(self : BackTranslateResult) -> String {
  self.nucleotide
}

///|
pub fn BackTranslateResult::table_id(self : BackTranslateResult) -> Int32 {
  self.table_id
}

///|
pub fn BackTranslateResult::mode(self : BackTranslateResult) -> String {
  self.mode
}

///|
/// Simple LCG random for deterministic "random" mode
priv struct BtRng {
  mut state : Int32
}

///|
fn BtRng::next(self : BtRng) -> Double {
  self.state = (self.state * 1103515245 + 12345) % 2000000000
  self.state.to_double() / 2000000000.0
}

///|
/// Back-translate a protein sequence to a nucleotide sequence.
///
/// mode: "first" (default), "random", or "biased"
/// For "biased" mode, pass a codon usage table (codon -> count).
/// Returns the DNA sequence as a string.
pub fn back_translate(
  protein : String,
  mode~ : String = "first",
  table_id~ : Int32 = 1,
  usage_table~ : Map[String, Double] = Map([]),
) -> String {
  let reverse_map = codon_table_reverse(table_id~)
  let seq_upper = protein.to_upper()
  let n = seq_upper.length()
  let buf = StringBuilder()
  let rng = BtRng::{ state: 42, }
  let mut i = 0
  while i < n {
    let aa = seq_upper[i:i + 1].to_owned()
    if aa == "*" {
      // Stop codon: pick first available stop codon
      let stops = codon_table_stop_codons(table_id)
      if stops.length() > 0 {
        buf.write_string(stops[0])
      }
    } else {
      match reverse_map.get(aa) {
        Some(codons) =>
          if codons.length() == 0 {
            // Unknown amino acid, skip
            ()
          } else if mode == "first" {
            buf.write_string(codons[0])
          } else if mode == "random" {
            let r = rng.next()
            let idx = (r * codons.length().to_double()).floor().to_int()
            let safe_idx = if idx >= codons.length() {
              codons.length() - 1
            } else if idx < 0 {
              0
            } else {
              idx
            }
            buf.write_string(codons[safe_idx])
          } else if mode == "biased" {
            // Pick the codon with highest usage frequency
            let mut best_codon = codons[0]
            let mut best_freq = -1.0
            let mut j = 0
            while j < codons.length() {
              let freq = match usage_table.get(codons[j]) {
                Some(f) => f
                None => 0.0
              }
              if freq > best_freq {
                best_freq = freq
                best_codon = codons[j]
              }
              j = j + 1
            }
            buf.write_string(best_codon)
          } else {
            // Default: first
            buf.write_string(codons[0])
          }
        None => ()
        // Unknown amino acid, skip
      }
    }
    i = i + 1
  }
  buf.to_string()
}

///|
/// Back-translate with full result object
pub fn back_translate_detailed(
  protein : String,
  mode~ : String = "first",
  table_id~ : Int32 = 1,
  usage_table~ : Map[String, Double] = Map([]),
) -> BackTranslateResult {
  let nt = back_translate(protein, mode~, table_id~, usage_table~)
  { protein, nucleotide: nt, table_id, mode, }
}

///|
/// Back-translate using the most frequent codon for each amino acid
/// (convenience wrapper for "biased" mode with a usage table)
pub fn back_translate_optimal(
  protein : String,
  usage_table : Map[String, Double],
  table_id~ : Int32 = 1,
) -> String {
  back_translate(protein, mode="biased", table_id~, usage_table~)
}

///|
/// Generate all possible back-translations for a short protein sequence.
/// Warning: the number of possibilities grows exponentially (up to 6^n).
/// Only practical for sequences shorter than ~15 amino acids.
pub fn back_translate_all(
  protein : String,
  table_id~ : Int32 = 1,
  max_results~ : Int32 = 10000,
) -> Array[String] {
  let reverse_map = codon_table_reverse(table_id~)
  let seq_upper = protein.to_upper()
  let n = seq_upper.length()
  let mut results : Array[String] = [""]
  let mut i = 0
  while i < n {
    let aa = seq_upper[i:i + 1].to_owned()
    let codons : Array[String] = if aa == "*" {
      codon_table_stop_codons(table_id)
    } else {
      match reverse_map.get(aa) {
        Some(cs) => cs
        None => []
      }
    }
    if codons.length() == 0 {
      // Unknown AA, keep results as-is
      ()
    } else {
      let new_results : Array[String] = []
      let mut j = 0
      while j < results.length() {
        let mut k = 0
        while k < codons.length() {
          new_results.push(results[j] + codons[k])
          if new_results.length() >= max_results {
            return new_results
          }
          k = k + 1
        }
        j = j + 1
      }
      results = new_results
    }
    i = i + 1
  }
  results
}

///|
/// Count the number of possible back-translations (without generating them).
/// Product of degeneracy for each amino acid.
pub fn back_translate_count(protein : String, table_id~ : Int32 = 1) -> Int64 {
  let reverse_map = codon_table_reverse(table_id~)
  let seq_upper = protein.to_upper()
  let n = seq_upper.length()
  let mut count : Int64 = 1L
  let mut i = 0
  while i < n {
    let aa = seq_upper[i:i + 1].to_owned()
    if aa == "*" {
      let stops = codon_table_stop_codons(table_id)
      count = count * stops.length().to_int64()
    } else {
      match reverse_map.get(aa) {
        Some(cs) => count = count * cs.length().to_int64()
        None => ()
      }
    }
    i = i + 1
  }
  count
}
