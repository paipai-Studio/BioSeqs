/// Biotite `sequence.rules` — rule-based prokaryotic gene finding.
///
/// `biotite_find_genes` scans both strands for gene candidates that
/// follow the bacterial transcription/translation pattern:
/// promoter (e.g. Pribnow box) -> ribosome binding site -> start codon
/// -> in-frame stop codon -> terminator downstream of the stop codon.

///|
/// A gene candidate in the coordinate system of the original (forward)
/// strand. The gene spans the half-open interval [start, stop).
pub struct BiotiteGeneProduct {
  start : Int32
  stop : Int32
  codon_start : Int32
  direction : Int32
  stop_far_before_start : Bool
}

///|
/// Find the first occurrence of `pattern` in `seq` with the search
/// confined to [from, to). Returns None when absent.
fn biotite_rules_find(
  seq : String,
  pattern : String,
  from : Int32,
  to : Int32,
) -> Int32? {
  let n = seq.length()
  let m = pattern.length()
  if m == 0 || from < 0 || to > n || from + m > to {
    return None
  }
  let first = pattern[0]
  for i = from; i + m <= to; i = i + 1 {
    if seq[i] == first && seq[i:i + m].to_owned() == pattern {
      return Some(i)
    }
  }
  None
}

///|
/// Whether `pattern` occurs anywhere in [from, to).
fn biotite_rules_contains(
  seq : String,
  pattern : String,
  from : Int32,
  to : Int32,
) -> Bool {
  biotite_rules_find(seq, pattern, from, to) is Some(_)
}

///|
/// Scan one strand (already oriented) for genes.
fn biotite_find_genes_strand(
  seq : String,
  promoter : String,
  terminator : String,
  rbs : String,
  start_codons : Array[String],
  stop_codons : Array[String],
  min_protein_len : Int32,
  max_protein_len : Int32,
  rbs_min_gap : Int32,
  rbs_max_gap : Int32,
  rbs_search_dist : Int32,
  direction : Int32,
  strand_len : Int32,
) -> Array[BiotiteGeneProduct] {
  let out : Array[BiotiteGeneProduct] = []
  let n = seq.length()
  for p = 0; p < n; p = p + 1 {
    let prom_end = p + promoter.length()
    if prom_end > n {
      break
    }
    if seq[p:prom_end].to_owned() != promoter {
      continue
    }
    // RBS within rbs_search_dist downstream of the promoter.
    let rbs_hi = if p + rbs_search_dist + rbs.length() < n {
      p + rbs_search_dist + rbs.length()
    } else {
      n
    }
    let r_pos = biotite_rules_find(seq, rbs, p, rbs_hi)
    let r = match r_pos {
      Some(r) => r
      None => continue
    }
    // Start codons in [r + rbs_min_gap, r + rbs_max_gap].
    let s_lo = r + rbs.length() + rbs_min_gap
    let s_hi = r + rbs.length() + rbs_max_gap
    for codon in start_codons {
      let mut s = s_lo
      while s <= s_hi {
        let found = biotite_rules_find(seq, codon, s, s_hi + codon.length())
        let start_pos = match found {
          Some(sp) => sp
          None => break
        }
        // Scan in-frame codons for the first stop codon.
        let codon_start = start_pos + codon.length()
        let mut stop_pos = -1
        let mut pos = codon_start
        while pos + 3 <= n {
          let codon3 = seq[pos:pos + 3].to_owned()
          let mut is_stop = false
          for sc in stop_codons {
            if codon3 == sc {
              is_stop = true
              break
            }
          }
          if is_stop {
            stop_pos = pos
            break
          }
          pos = pos + 3
        }
        if stop_pos < 0 {
          // No in-frame stop codon: incomplete gene, runs to the end.
          let gene_end = n
          let protein_len = (gene_end - codon_start) / 3
          if protein_len >= min_protein_len {
            out.push(BiotiteGeneProduct::{
              start: biotite_rules_map_coord(
                if direction > 0 {
                  start_pos
                } else {
                  gene_end
                },
                strand_len,
                direction,
              ),
              stop: biotite_rules_map_coord(
                if direction > 0 {
                  gene_end
                } else {
                  start_pos
                },
                strand_len,
                direction,
              ),
              codon_start: biotite_rules_map_coord(
                codon_start, strand_len, direction,
              ),
              direction,
              stop_far_before_start: true,
            })
          }
          s = start_pos + 1
          continue
        }
        let gene_end = stop_pos + 3
        let protein_len = (stop_pos - codon_start) / 3
        if protein_len < min_protein_len ||
          (max_protein_len > 0 && protein_len > max_protein_len) {
          s = start_pos + 1
          continue
        }
        // Terminator downstream of the stop codon.
        let term_hi = if gene_end + 50 < n { gene_end + 50 } else { n }
        if biotite_rules_contains(seq, terminator, gene_end, term_hi) {
          out.push(BiotiteGeneProduct::{
            start: biotite_rules_map_coord(
              if direction > 0 {
                start_pos
              } else {
                gene_end
              },
              strand_len,
              direction,
            ),
            stop: biotite_rules_map_coord(
              if direction > 0 {
                gene_end
              } else {
                start_pos
              },
              strand_len,
              direction,
            ),
            codon_start: biotite_rules_map_coord(
              codon_start, strand_len, direction,
            ),
            direction,
            stop_far_before_start: false,
          })
        }
        s = start_pos + 1
      }
    }
  }
  out
}

///|
/// Map a forward-strand-local coordinate to the original sequence.
fn biotite_rules_map_coord(
  pos_local : Int32,
  strand_len : Int32,
  direction : Int32,
) -> Int32 {
  if direction > 0 {
    pos_local
  } else {
    strand_len - pos_local
  }
}

///|
/// Rule-based prokaryotic gene finding (both strands).
///
/// Parameters:
/// - `promoter`: transcription initiation motif (e.g. "TATAAT")
/// - `terminator`: transcription termination motif downstream of the
///   stop codon (e.g. poly-A run)
/// - `rbs`: ribosome binding site (Shine-Dalgarno, e.g. "AGGAGGT")
/// - `start_codons` / `stop_codons`: codon lists
/// - `min_protein_len` / `max_protein_len`: amino-acid length bounds
///   (`max_protein_len = None` means unbounded)
/// - `rbs_min_gap` / `rbs_max_gap`: allowed gap between RBS end and
///   the start codon
/// - `rbs_search_dist`: max distance from the promoter to the RBS
pub fn biotite_find_genes(
  sequence : String,
  promoter : String,
  terminator : String,
  rbs : String,
  start_codons : Array[String],
  stop_codons : Array[String],
  min_protein_len~ : Int32 = 30,
  max_protein_len~ : Int32? = None,
  rbs_min_gap~ : Int32 = 5,
  rbs_max_gap~ : Int32 = 13,
  rbs_search_dist~ : Int32 = 100,
) -> Array[BiotiteGeneProduct] {
  let n = sequence.length()
  let max_len = match max_protein_len {
    Some(v) => v
    None => 0
  }
  // Forward strand.
  let forward = biotite_find_genes_strand(
    sequence, promoter, terminator, rbs, start_codons, stop_codons, min_protein_len,
    max_len, rbs_min_gap, rbs_max_gap, rbs_search_dist, 1, n,
  )
  // Reverse strand: search the reverse complement and map back.
  let rc = biostrings_reverse_complement(sequence)
  let reverse = biotite_find_genes_strand(
    rc, promoter, terminator, rbs, start_codons, stop_codons, min_protein_len, max_len,
    rbs_min_gap, rbs_max_gap, rbs_search_dist, -1, n,
  )
  let out : Array[BiotiteGeneProduct] = []
  for g in forward {
    out.push(g)
  }
  for g in reverse {
    out.push(g)
  }
  out.sort_by(fn(a, b) { Int32::compare(a.start, b.start) })
  out
}

///|
pub fn _biotite_sequence_rules_dummy() -> Unit {
  let genes = biotite_find_genes(
    "ATGAAATAA",
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA"],
  )
  let _ = genes
  let _ = biotite_rules_find("ATGAAATAA", "ATG", 0, 9)
  let _ = biotite_rules_contains("ATGAAATAA", "TAA", 0, 9)
  let _ = biotite_rules_map_coord(3, 9, -1)
}
