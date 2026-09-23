///|
/// primersearch - EMBOSS primersearch style primer matching and PCR product prediction.
/// Inspired by Biopython Bio.Emboss.PrimerSearch module.
/// Provides: primer pair file parsing, mismatch-tolerant search,
/// PCR product size inference, and primersearch output formatting.

///|
/// A primer pair: name + forward + reverse sequences.
pub struct PrimerPair {
  name : String
  forward : String
  reverse : String
}

///|
/// A single primer hit on a target sequence.
pub struct PrimerHit {
  strand : String // "+" forward match, "-" reverse-complement match
  position : Int
  mismatches : Int
}

///|
/// PCR product inferred from a primer pair on a sequence.
pub struct PcrProduct {
  primer_name : String
  sequence_id : String
  fwd_hit : PrimerHit
  rev_hit : PrimerHit
  product_size : Int
}

///|
/// Results of running primersearch: all products per sequence.
pub struct PrimerSearchResult {
  products : Array[PcrProduct]
  n_primers : Int
  n_sequences : Int
}

///|
/// Construct a primer pair.
pub fn PrimerPair::new(name : String, forward : String, reverse : String) -> PrimerPair {
  PrimerPair::{ name, forward, reverse }
}

///|
/// Parse an EMBOSS primers file: one primer pair per line "name FWD REV".
/// Lines starting with '#' are comments.
pub fn primersearch_read_primers(content : String) -> Array[PrimerPair] {
  let pairs : Array[PrimerPair] = []
  let lines = content.split("\n").to_array()
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim(chars=" \t\r").to_string()
    if line.length() > 0 && line[0] != 35 {
      // '#' comment
      let fields : Array[String] = []
      let mut start = 0
      let mut j = 0
      while j <= line.length() {
        let is_sep = j == line.length() || line[j] == 32 || line[j] == 9
        if is_sep {
          if j > start {
            fields.push(line[start:j].to_string())
          }
          start = j + 1
        }
        j = j + 1
      }
      if fields.length() >= 3 {
        pairs.push(PrimerPair::new(fields[0], fields[1], fields[2]))
      }
    }
    i = i + 1
  }
  pairs
}

///|
/// Reverse complement a DNA sequence.
pub fn primersearch_revcomp(seq : String) -> String {
  let sb = StringBuilder()
  let mut i = seq.length() - 1
  while i >= 0 {
    let c = seq[i]
    let comp : UInt16 = match c {
      65 => 84 // A -> T
      84 => 65 // T -> A
      67 => 71 // C -> G
      71 => 67 // G -> C
      97 => 116 // a -> t
      116 => 97 // t -> a
      99 => 103 // c -> g
      103 => 99 // g -> c
      _ => c
    }
    sb.write_char(comp.unsafe_to_char())
    i = i - 1
  }
  sb.to_string()
}

///|
/// Count mismatches between two equal-length strings.
fn primersearch_mismatches(a : String, b : String) -> Int {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut mm = 0
  let mut i = 0
  while i < n {
    if a[i] != b[i] {
      mm = mm + 1
    }
    i = i + 1
  }
  mm
}

///|
/// Search for a primer (both strands) in a sequence allowing up to max_mismatches.
///   Returns all hits (position, strand, mismatches).
pub fn primersearch_find_hits(
  seq : String,
  primer : String,
  max_mismatches~ : Int = 0
) -> Array[PrimerHit] {
  let hits : Array[PrimerHit] = []
  let plen = primer.length()
  if plen == 0 || seq.length() < plen {
    return hits
  }
  let rc = primersearch_revcomp(primer)
  let mut i = 0
  while i <= seq.length() - plen {
    let window = seq[i:i + plen].to_string()
    let mm_fwd = primersearch_mismatches(window, primer)
    if mm_fwd <= max_mismatches {
      hits.push(PrimerHit::{ strand: "+", position: i, mismatches: mm_fwd })
    } else {
      let mm_rev = primersearch_mismatches(window, rc)
      if mm_rev <= max_mismatches {
        hits.push(PrimerHit::{ strand: "-", position: i, mismatches: mm_rev })
      }
    }
    i = i + 1
  }
  hits
}

///|
/// Infer PCR products: for each primer pair, find forward and reverse hits
/// on opposite strands and compute product size between primer 3' ends.
///   EMBOSS convention: product spans from fwd primer start to rev primer end.
pub fn primersearch_run(
  primers : Array[PrimerPair],
  sequences : Array[(String, String)],
  max_mismatches~ : Int = 0
) -> PrimerSearchResult {
  let products : Array[PcrProduct] = []
  let mut p = 0
  while p < primers.length() {
    let pair = primers[p]
    let rev = primersearch_revcomp(pair.reverse)
    let mut s = 0
    while s < sequences.length() {
      let (seq_id, seq) = sequences[s]
      let fwd_hits = primersearch_find_hits(seq, pair.forward, max_mismatches~)
      let rev_hits = primersearch_find_hits(seq, rev, max_mismatches~)
      // Pair opposite-strand hits: fwd "+", rev "-" (or fwd "-", rev "+")
      let mut fi = 0
      while fi < fwd_hits.length() {
        let fh = fwd_hits[fi]
        let mut ri = 0
        while ri < rev_hits.length() {
          let rh = rev_hits[ri]
          // Need opposite strands and fwd before rev (or rev before fwd)
          let fwd_start = if fh.strand == "+" { fh.position } else { fh.position + pair.forward.length() - 1 }
          let rev_end = if rh.strand == "-" { rh.position + pair.reverse.length() - 1 } else { rh.position }
          let rev_start = rh.position
          if fh.strand == "+" && rh.strand == "-" && rh.position >= fh.position {
            let size = rev_end - fwd_start + 1
            products.push(PcrProduct::{ primer_name: pair.name, sequence_id: seq_id, fwd_hit: fh, rev_hit: rh, product_size: size })
          } else if fh.strand == "-" && rh.strand == "+" && fh.position >= rh.position {
            // Reverse orientation: product from rev primer to fwd primer
            let fwd_end = fh.position + pair.forward.length() - 1
            let size = fwd_end - rev_start + 1
            products.push(PcrProduct::{ primer_name: pair.name, sequence_id: seq_id, fwd_hit: fh, rev_hit: rh, product_size: size })
          }
          ri = ri + 1
        }
        fi = fi + 1
      }
      s = s + 1
    }
    p = p + 1
  }
  PrimerSearchResult::{ products, n_primers: primers.length(), n_sequences: sequences.length() }
}

///|
/// Format results in EMBOSS primersearch output style.
pub fn primersearch_to_string(result : PrimerSearchResult) -> String {
  let sb = StringBuilder()
  sb.write_string("Primer searching\n\n")
  let mut i = 0
  while i < result.products.length() {
    let prod = result.products[i]
    sb.write_string(prod.primer_name)
    sb.write_string(" sequence ")
    sb.write_string(prod.sequence_id)
    sb.write_string("; amplicon length ")
    sb.write_string(prod.product_size.to_string())
    sb.write_string(" bp\n  forward hit: strand=")
    sb.write_string(prod.fwd_hit.strand)
    sb.write_string(" pos=")
    sb.write_string(prod.fwd_hit.position.to_string())
    sb.write_string(" mm=")
    sb.write_string(prod.fwd_hit.mismatches.to_string())
    sb.write_string("\n  reverse hit: strand=")
    sb.write_string(prod.rev_hit.strand)
    sb.write_string(" pos=")
    sb.write_string(prod.rev_hit.position.to_string())
    sb.write_string(" mm=")
    sb.write_string(prod.rev_hit.mismatches.to_string())
    sb.write_string("\n")
    i = i + 1
  }
  sb.to_string()
}

///|
/// Filter products by size range [min_size, max_size].
pub fn primersearch_filter_by_size(
  result : PrimerSearchResult,
  min_size : Int,
  max_size : Int
) -> Array[PcrProduct] {
  let filtered : Array[PcrProduct] = []
  let mut i = 0
  while i < result.products.length() {
    let prod = result.products[i]
    if prod.product_size >= min_size && prod.product_size <= max_size {
      filtered.push(prod)
    }
    i = i + 1
  }
  filtered
}

///|
/// Example primers file content.
pub fn primersearch_example_primers() -> Array[PrimerPair] {
  [
    PrimerPair::new("GeneA_F", "ATGAAACGTT", "TTACAGCTAA"),
  ]
}

///|
/// Example target sequences.
pub fn primersearch_example_sequences() -> Array[(String, String)] {
  [
    ("chr1", "GGGATGAAACGTTCCCTTTACAGCTAAAGGG"),
  ]
}