// ============================================================================
// Prodigal-lite: Bacterial gene prediction
// ============================================================================
// A lite reimplementation of Prodigal (PROkaryotic DYnamic programming
// Gene-finding ALgorithm) for bacterial/archaeal gene prediction.
//
// Pipeline:
//   1. Find all ORFs in 6 reading frames (3 forward + 3 reverse)
//   2. Score each ORF: codon usage bias, GC content, length, RBS motif
//   3. Dynamic programming to select non-overlapping gene set
//   4. Output as GFF3 / FASTA protein
//
// References:
//   Hyatt et al. "Prodigal: prokaryotic gene recognition and translation
//   initiation site identification." BMC Bioinformatics 11, 119 (2010).

// ============================================================================
// Data structures
// ============================================================================

///|
/// A predicted gene.
pub struct ProdigalGene {
  pub seqid : String
  pub start : Int         // 1-based start position (inclusive)
  pub end_ : Int          // 1-based end position (inclusive)
  pub strand : String     // "+" or "-"
  pub frame : Int         // 0, 1, or 2
  pub start_codon : String
  pub stop_codon : String
  pub length : Int        // gene length in nucleotides (including stop)
  pub score : Double      // overall gene score
  pub partial : Bool      // partial gene at sequence boundary
  pub rbs_score : Double  // RBS motif score
  pub gc_content : Double // GC content of the gene
  pub protein : String    // translated protein sequence
  pub gene_id : String    // assigned gene ID
}

///|
/// Gene prediction result.
pub struct ProdigalResult {
  pub seqid : String
  pub sequence_length : Int
  pub genes : Array[ProdigalGene]
  pub method : String
  pub num_genes : Int
  pub mean_gc : Double
}

///|
/// Parameters for gene prediction.
pub struct ProdigalParams {
  pub min_gene_len : Int      // minimum gene length (nt)
  pub max_overlap : Int       // maximum allowed overlap (nt)
  pub start_codons : Array[String]
  pub stop_codons : Array[String]
  pub use_rbs : Bool          // use RBS scoring
  pub allow_partial : Bool    // allow partial genes at boundaries
  pub closed_ends : Bool      // closed ends (no partial genes)
}

///|
/// Default parameters for bacterial gene prediction.
pub fn prodigal_default_params() -> ProdigalParams {
  {
    min_gene_len: 60,
    max_overlap: 60,
    start_codons: ["ATG"],
    stop_codons: ["TAA", "TAG", "TGA"],
    use_rbs: true,
    allow_partial: true,
    closed_ends: false,
  }
}

///|
/// Create custom parameters.
pub fn prodigal_params(
  min_gene_len : Int,
  max_overlap : Int,
) -> ProdigalParams {
  {
    min_gene_len,
    max_overlap,
    start_codons: ["ATG"],
    stop_codons: ["TAA", "TAG", "TGA"],
    use_rbs: true,
    allow_partial: true,
    closed_ends: false,
  }
}

///|
/// Create ProdigalParams with all fields specified.
pub fn prodigal_params_full(
  min_gene_len : Int,
  max_overlap : Int,
  start_codons : Array[String],
  stop_codons : Array[String],
  use_rbs : Bool,
  allow_partial : Bool,
  closed_ends : Bool,
) -> ProdigalParams {
  {
    min_gene_len,
    max_overlap,
    start_codons,
    stop_codons,
    use_rbs,
    allow_partial,
    closed_ends,
  }
}

// ============================================================================
// 1. ORF finding in 6 reading frames
// ============================================================================

///|
/// An ORF candidate found during scanning.
struct OrfCandidate {
  start_pos : Int    // 0-based start position in the original sequence
  end_pos : Int      // 0-based end position (exclusive) in the original sequence
  strand : String    // "+" or "-"
  frame : Int        // 0, 1, or 2
  start_codon : String
  stop_codon : String
  partial_start : Bool
  partial_end : Bool
}

///|
/// Check if a codon is in a list.
fn is_codon(codon : String, codons : Array[String]) -> Bool {
  for c in codons {
    if codon == c {
      return true
    }
  }
  false
}

///|
/// Extract a substring of length n starting at position i.
fn extract_codon(seq : String, i : Int, n : Int) -> String {
  if i < 0 || i + n > seq.length() {
    return ""
  }
  let buf : StringBuilder = StringBuilder()
  for k = 0; k < n; k = k + 1 {
    buf.write_char(seq[i + k].unsafe_to_char())
  }
  buf.to_string()
}

///|
/// Find all ORFs in one forward frame.
fn find_orfs_forward_frame(
  seq : String,
  frame : Int,
  params : ProdigalParams,
) -> Array[OrfCandidate] {
  let n = seq.length()
  let result : Array[OrfCandidate] = Array(capacity=0)
  let mut i = frame
  while i < n - 2 {
    let codon = extract_codon(seq, i, 3)
    if is_codon(codon, params.start_codons) {
      // Search for stop codon
      let mut j = i + 3
      let mut found_stop = false
      while j < n - 2 {
        let stop_codon = extract_codon(seq, j, 3)
        if is_codon(stop_codon, params.stop_codons) {
          let gene_len = j + 3 - i
          if gene_len >= params.min_gene_len {
            result.push({
              start_pos: i,
              end_pos: j + 3,
              strand: "+",
              frame,
              start_codon: codon,
              stop_codon: stop_codon,
              partial_start: false,
              partial_end: false,
            })
          }
          found_stop = true
          // Move to after this stop codon
          i = j + 3
          break
        }
        j += 3
      }
      if not(found_stop) {
        // Partial ORF (no stop codon found before end)
        if params.allow_partial && n - i >= params.min_gene_len {
          result.push({
            start_pos: i,
            end_pos: n,
            strand: "+",
            frame,
            start_codon: codon,
            stop_codon: "",
            partial_start: false,
            partial_end: true,
          })
        }
        i = n
        break
      }
    } else {
      i += 3
    }
  }
  result
}

///|
/// Find all ORFs in one reverse frame.
fn find_orfs_reverse_frame(
  seq : String,
  rc_seq : String,
  frame : Int,
  params : ProdigalParams,
) -> Array[OrfCandidate] {
  let n = rc_seq.length()
  let result : Array[OrfCandidate] = Array(capacity=0)
  let mut i = frame
  while i < n - 2 {
    let codon = extract_codon(rc_seq, i, 3)
    if is_codon(codon, params.start_codons) {
      let mut j = i + 3
      let mut found_stop = false
      while j < n - 2 {
        let stop_codon = extract_codon(rc_seq, j, 3)
        if is_codon(stop_codon, params.stop_codons) {
          let gene_len = j + 3 - i
          if gene_len >= params.min_gene_len {
            // Map back to forward strand coordinates
            // rc position i corresponds to forward position n - 1 - (i + len - 1)
            let fwd_start = n - (j + 3) // 0-based start on forward strand
            let fwd_end = n - i         // 0-based end (exclusive) on forward strand
            result.push({
              start_pos: fwd_start,
              end_pos: fwd_end,
              strand: "-",
              frame,
              start_codon: codon,
              stop_codon: stop_codon,
              partial_start: false,
              partial_end: false,
            })
          }
          found_stop = true
          i = j + 3
          break
        }
        j += 3
      }
      if not(found_stop) {
        if params.allow_partial && n - i >= params.min_gene_len {
          let fwd_start = 0
          let fwd_end = n - i
          result.push({
            start_pos: fwd_start,
            end_pos: fwd_end,
            strand: "-",
            frame,
            start_codon: codon,
            stop_codon: "",
            partial_start: false,
            partial_end: true,
          })
        }
        i = n
        break
      }
    } else {
      i += 3
    }
  }
  result
}

///|
/// Find all ORFs in all 6 reading frames.
fn find_all_orfs(
  seq : String,
  params : ProdigalParams,
) -> Array[OrfCandidate] {
  let rc_seq = reverse_complement(seq)
  let result : Array[OrfCandidate] = Array(capacity=0)
  // Forward frames 0, 1, 2
  for frame in [0, 1, 2] {
    let orfs = find_orfs_forward_frame(seq, frame, params)
    for orf in orfs {
      result.push(orf)
    }
  }
  // Reverse frames 0, 1, 2
  for frame in [0, 1, 2] {
    let orfs = find_orfs_reverse_frame(seq, rc_seq, frame, params)
    for orf in orfs {
      result.push(orf)
    }
  }
  result
}

// ============================================================================
// 2. ORF scoring
// ============================================================================

///|
/// Standard genetic code: codon -> amino acid.
fn genetic_code() -> Map[String, String] {
  let code : Map[String, String] = Map([], capacity=0)
  let pairs = [
    ("TTT", "F"), ("TTC", "F"), ("TTA", "L"), ("TTG", "L"),
    ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"),
    ("ATT", "I"), ("ATC", "I"), ("ATA", "I"), ("ATG", "M"),
    ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
    ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"),
    ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
    ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
    ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
    ("TAT", "Y"), ("TAC", "Y"), ("TAA", "*"), ("TAG", "*"),
    ("CAT", "H"), ("CAC", "H"), ("CAA", "Q"), ("CAG", "Q"),
    ("AAT", "N"), ("AAC", "N"), ("AAA", "K"), ("AAG", "K"),
    ("GAT", "D"), ("GAC", "D"), ("GAA", "E"), ("GAG", "E"),
    ("TGT", "C"), ("TGC", "C"), ("TGA", "*"), ("TGG", "W"),
    ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"),
    ("AGT", "S"), ("AGC", "S"), ("AGA", "R"), ("AGG", "R"),
    ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
  ]
  for (codon, aa) in pairs {
    code[codon] = aa
  }
  code
}

///|
/// Translate a DNA sequence to protein.
fn translate_seq(seq : String) -> String {
  let code = genetic_code()
  let buf : StringBuilder = StringBuilder()
  let n = seq.length()
  let mut i = 0
  while i + 2 < n {
    let codon = extract_codon(seq, i, 3)
    let aa = match code.get(codon) {
      Some(a) => a
      None => "X"
    }
    if aa == "*" {
      buf.write_char('*')
    } else {
      buf.write_string(aa)
    }
    i += 3
  }
  buf.to_string()
}

///|
/// Extract a subsequence from a string.
fn subseq(seq : String, start : Int, end_ : Int) -> String {
  let buf : StringBuilder = StringBuilder()
  for i = start; i < end_; i = i + 1 {
    if i >= 0 && i < seq.length() {
      buf.write_char(seq[i].unsafe_to_char())
    }
  }
  buf.to_string()
}

///|
/// Compute GC content of a sequence.
fn prodigal_gc_content(seq : String) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  let mut gc = 0
  for i = 0; i < n; i = i + 1 {
    let c = seq[i].unsafe_to_char()
    if c == 'G' || c == 'C' {
      gc += 1
    }
  }
  gc.to_double() / n.to_double()
}

///|
/// Compute codon usage bias score.
/// Uses a simple model: count codons and compare to expected uniform frequency.
fn codon_bias_score(seq : String) -> Double {
  let n = seq.length()
  if n < 6 {
    return 0.0
  }
  // Count codons
  let counts : Map[String, Int] = Map([], capacity=0)
  let mut total_codons = 0
  let mut i = 0
  while i + 2 < n {
    let codon = extract_codon(seq, i, 3)
    if codon.length() == 3 {
      match counts.get(codon) {
        Some(c) => counts[codon] = c + 1
        None => counts[codon] = 1
      }
      total_codons += 1
    }
    i += 3
  }
  if total_codons == 0 {
    return 0.0
  }
  // Score: reward ATG starts and common codons
  // Simple heuristic: prefer genes with more "common" bacterial codons
  let common_codons = [
    "ATG", "GTG", "TTG",  // start codons
    "CTG", "ATC", "GTC",  // common Leu, Ile, Val
    "GCT", "GCC", "ACC",  // common Ala, Thr
    "TTC", "TAC", "TGC",  // common Phe, Tyr, Cys
    "CAG", "AAC", "GAC",  // common Gln, Asn, Asp
    "GTG", "GCG", "CCG",  // common Val, Ala, Pro
  ]
  let mut score = 0.0
  for (codon, count) in counts {
    let freq = count.to_double() / total_codons.to_double()
    // Base score: information content
    if freq > 0.0 {
      score += freq * @math.ln(freq * 64.0)
    }
    // Bonus for common codons
    for cc in common_codons {
      if codon == cc {
        score += freq * 0.5
        break
      }
    }
  }
  score
}

///|
/// Score RBS (Shine-Dalgarno) motif upstream of start codon.
/// Looks for AGGAGG-like motifs in the -20 to -4 region.
fn rbs_score(
  seq : String,
  start_pos : Int,
  strand : String,
) -> Double {
  // Extract upstream region (20 nt before start)
  let upstream_start = if start_pos >= 20 { start_pos - 20 } else { 0 }
  let upstream_end = start_pos
  let upstream = subseq(seq, upstream_start, upstream_end)
  let rc_upstream = if strand == "-" {
    reverse_complement(upstream)
  } else {
    upstream
  }
  // Search for Shine-Dalgarno motifs
  let rbs_motifs = [
    "AGGAGG", "AGGAG", "GGAGG", "GAGGT", "AGGUG",
    "GGAG", "AGGA", "GAGG", "TAAGG", "AGGU",
  ]
  let mut best_score = 0.0
  for motif in rbs_motifs {
    let mlen = motif.length()
    let ulen = rc_upstream.length()
    if ulen < mlen {
      continue
    }
    for i = 0; i <= ulen - mlen; i = i + 1 {
      let mut matches = 0
      for k = 0; k < mlen; k = k + 1 {
        if rc_upstream[i + k].unsafe_to_char() == motif[k].unsafe_to_char() {
          matches += 1
        }
      }
      // Score: fraction of matching positions, weighted by motif strength
      let frac = matches.to_double() / mlen.to_double()
      // Position bonus: closer to -6 is better
      let pos_in_upstream = i
      let dist_from_start = ulen - pos_in_upstream - mlen
      let pos_bonus = if dist_from_start >= 3 && dist_from_start <= 12 {
        1.0
      } else if dist_from_start >= 2 && dist_from_start <= 15 {
        0.5
      } else {
        0.2
      }
      let s = frac * pos_bonus
      if s > best_score {
        best_score = s
      }
    }
  }
  best_score
}

///|
/// Score an ORF candidate and produce a ProdigalGene.
fn score_orf(
  orf : OrfCandidate,
  seq : String,
  rc_seq : String,
  params : ProdigalParams,
  gene_idx : Int,
  seqid : String,
) -> (ProdigalGene, Double) {
  // Extract the gene sequence
  let gene_seq = if orf.strand == "+" {
    subseq(seq, orf.start_pos, orf.end_pos)
  } else {
    subseq(rc_seq, seq.length() - orf.end_pos, seq.length() - orf.start_pos)
  }
  // Translate to protein
  let protein = translate_seq(gene_seq)
  // Compute scores
  let gc = prodigal_gc_content(gene_seq)
  let codon_score = codon_bias_score(gene_seq)
  let rbs = if params.use_rbs {
    rbs_score(seq, orf.start_pos, orf.strand)
  } else {
    0.0
  }
  // Length score: longer genes are more likely real
  let len_score = @math.ln((orf.end_pos - orf.start_pos).to_double()) / 10.0
  // Combined score
  let total_score = codon_score + rbs * 2.0 + len_score + gc * 0.5
  let gene_id = "gene_" + seqid + "_" + gene_idx.to_string()
  let gene : ProdigalGene = {
    seqid,
    start: orf.start_pos + 1,  // 1-based
    end_: orf.end_pos,          // 1-based (end_pos is exclusive 0-based = inclusive 1-based)
    strand: orf.strand,
    frame: orf.frame,
    start_codon: orf.start_codon,
    stop_codon: orf.stop_codon,
    length: orf.end_pos - orf.start_pos,
    score: total_score,
    partial: orf.partial_start || orf.partial_end,
    rbs_score: rbs,
    gc_content: gc,
    protein,
    gene_id,
  }
  (gene, total_score)
}

// ============================================================================
// 3. Gene selection: non-overlapping set via greedy + DP
// ============================================================================

///|
/// Check if two ORF candidates overlap.
fn orfs_overlap(a : OrfCandidate, b : OrfCandidate) -> Bool {
  // Different strands can overlap in Prodigal (but we avoid it for simplicity)
  if a.strand != b.strand {
    // On different strands, overlap if they share any position
    return a.start_pos < b.end_pos && b.start_pos < a.end_pos
  }
  a.start_pos < b.end_pos && b.start_pos < a.end_pos
}

///|
/// Compute overlap length between two ORF candidates.
fn overlap_len(a : OrfCandidate, b : OrfCandidate) -> Int {
  let overlap_start = if a.start_pos > b.start_pos { a.start_pos } else { b.start_pos }
  let overlap_end = if a.end_pos < b.end_pos { a.end_pos } else { b.end_pos }
  if overlap_end > overlap_start {
    overlap_end - overlap_start
  } else {
    0
  }
}

///|
/// Select non-overlapping gene set using greedy approach.
/// Sort by score descending, greedily pick genes that don't overlap
/// (or overlap by less than max_overlap).
fn select_genes(
  candidates : Array[(ProdigalGene, Double)],
  max_overlap : Int,
) -> Array[ProdigalGene] {
  // Copy and sort by score descending
  let sorted : Array[(ProdigalGene, Double)] = Array(capacity=candidates.length())
  for i = 0; i < candidates.length(); i = i + 1 {
    sorted.push(candidates[i])
  }
  // Simple sort by score descending
  for i = 0; i < sorted.length(); i = i + 1 {
    for j = i + 1; j < sorted.length(); j = j + 1 {
      let (gi, si) = sorted[i]
      let (gj, sj) = sorted[j]
      if sj > si {
        sorted[i] = (gj, sj)
        sorted[j] = (gi, si)
      }
    }
  }
  // Greedily select non-overlapping genes
  let selected : Array[ProdigalGene] = Array(capacity=0)
  for i = 0; i < sorted.length(); i = i + 1 {
    let (gene, _) = sorted[i]
    let mut ok = true
    for s in selected {
      let a_start = gene.start
      let a_end = gene.end_
      let b_start = s.start
      let b_end = s.end_
      let ov_start = if a_start > b_start { a_start } else { b_start }
      let ov_end = if a_end < b_end { a_end } else { b_end }
      let ov = if ov_end > ov_start { ov_end - ov_start } else { 0 }
      if ov > max_overlap {
        ok = false
        break
      }
    }
    if ok {
      selected.push(gene)
    }
  }
  // Sort selected genes by position
  for i = 0; i < selected.length(); i = i + 1 {
    for j = i + 1; j < selected.length(); j = j + 1 {
      if selected[j].start < selected[i].start {
        let tmp = selected[i]
        selected[i] = selected[j]
        selected[j] = tmp
      }
    }
  }
  selected
}

// ============================================================================
// 4. Public API
// ============================================================================

///|
/// Predict genes in a bacterial sequence.
pub fn prodigal_predict(
  sequence : String,
  seqid : String,
) -> ProdigalResult {
  prodigal_predict_with_params(sequence, seqid, prodigal_default_params())
}

///|
/// Predict genes with custom parameters.
pub fn prodigal_predict_with_params(
  sequence : String,
  seqid : String,
  params : ProdigalParams,
) -> ProdigalResult {
  let seq = sequence.to_upper()
  let n = seq.length()
  if n < 6 {
    return {
      seqid,
      sequence_length: n,
      genes: [],
      method: "prodigal",
      num_genes: 0,
      mean_gc: 0.0,
    }
  }
  let rc_seq = reverse_complement(seq)
  // Step 1: Find all ORFs
  let orfs = find_all_orfs(seq, params)
  // Step 2: Score each ORF
  let candidates : Array[(ProdigalGene, Double)] = Array(capacity=0)
  for i = 0; i < orfs.length(); i = i + 1 {
    let (gene, score) = score_orf(orfs[i], seq, rc_seq, params, i + 1, seqid)
    candidates.push((gene, score))
  }
  // Step 3: Select non-overlapping genes
  let selected = select_genes(candidates, params.max_overlap)
  // Compute mean GC
  let mut total_gc = 0.0
  for gene in selected {
    total_gc += gene.gc_content
  }
  let mean_gc = if selected.length() > 0 {
    total_gc / selected.length().to_double()
  } else {
    prodigal_gc_content(seq)
  }
  {
    seqid,
    sequence_length: n,
    genes: selected,
    method: "prodigal",
    num_genes: selected.length(),
    mean_gc,
  }
}

// ============================================================================
// 5. Output formatting
// ============================================================================

///|
/// Convert prediction result to GFF3 format.
pub fn prodigal_to_gff(result : ProdigalResult) -> String {
  let buf : StringBuilder = StringBuilder()
  buf.write_string("##gff-version 3\\n")
  buf.write_string("##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.2/Features.md\\n")
  buf.write_string("##sequence-region ")
  buf.write_string(result.seqid)
  buf.write_string(" 1 ")
  buf.write_string(result.sequence_length.to_string())
  buf.write_string("\\n")
  for gene in result.genes {
    // Gene feature
    buf.write_string(result.seqid)
    buf.write_string("\\tProdigal\\tgene\\t")
    buf.write_string(gene.start.to_string())
    buf.write_string("\\t")
    buf.write_string(gene.end_.to_string())
    buf.write_string("\\t")
    buf.write_string(gene.score.to_string())
    buf.write_string("\\t")
    buf.write_string(gene.strand)
    buf.write_string("\\t.\\tID=")
    buf.write_string(gene.gene_id)
    buf.write_string(";partial=")
    buf.write_string(if gene.partial { "true" } else { "false" })
    buf.write_string("\\n")
    // CDS feature
    buf.write_string(result.seqid)
    buf.write_string("\\tProdigal\\tCDS\\t")
    buf.write_string(gene.start.to_string())
    buf.write_string("\\t")
    buf.write_string(gene.end_.to_string())
    buf.write_string("\\t")
    buf.write_string(gene.score.to_string())
    buf.write_string("\\t")
    buf.write_string(gene.strand)
    buf.write_string("\\t0\\tID=")
    buf.write_string(gene.gene_id)
    buf.write_string("_cds;Parent=")
    buf.write_string(gene.gene_id)
    buf.write_string(";gc_content=")
    buf.write_string(gene.gc_content.to_string())
    buf.write_string(";rbs_score=")
    buf.write_string(gene.rbs_score.to_string())
    buf.write_string(";start_codon=")
    buf.write_string(gene.start_codon)
    buf.write_string(";stop_codon=")
    buf.write_string(gene.stop_codon)
    buf.write_string("\\n")
  }
  buf.to_string()
}

///|
/// Convert predicted genes to FASTA protein format.
pub fn prodigal_to_fasta(result : ProdigalResult) -> String {
  let buf : StringBuilder = StringBuilder()
  for gene in result.genes {
    buf.write_string(">")
    buf.write_string(gene.gene_id)
    buf.write_string(" | ")
    buf.write_string(result.seqid)
    buf.write_string(":")
    buf.write_string(gene.start.to_string())
    buf.write_string("-")
    buf.write_string(gene.end_.to_string())
    buf.write_string("(")
    buf.write_string(gene.strand)
    buf.write_string(") | score=")
    buf.write_string(gene.score.to_string())
    buf.write_string(" | length=")
    buf.write_string(gene.length.to_string())
    buf.write_string("\\n")
    buf.write_string(gene.protein)
    buf.write_string("\\n")
  }
  buf.to_string()
}

///|
/// Convert predicted genes to nucleotide FASTA format.
pub fn prodigal_to_fasta_nt(
  result : ProdigalResult,
  sequence : String,
) -> String {
  let seq = sequence.to_upper()
  let rc_seq = reverse_complement(seq)
  let buf : StringBuilder = StringBuilder()
  for gene in result.genes {
    buf.write_string(">")
    buf.write_string(gene.gene_id)
    buf.write_string(" | ")
    buf.write_string(result.seqid)
    buf.write_string(":")
    buf.write_string(gene.start.to_string())
    buf.write_string("-")
    buf.write_string(gene.end_.to_string())
    buf.write_string("(")
    buf.write_string(gene.strand)
    buf.write_string(")\\n")
    let nt_seq = if gene.strand == "+" {
      subseq(seq, gene.start - 1, gene.end_)
    } else {
      subseq(rc_seq, seq.length() - gene.end_, seq.length() - gene.start + 1)
    }
    buf.write_string(nt_seq)
    buf.write_string("\\n")
  }
  buf.to_string()
}

///|
/// Summary statistics as text.
pub fn prodigal_summary(result : ProdigalResult) -> String {
  let buf : StringBuilder = StringBuilder()
  buf.write_string("Prodigal gene prediction summary\\n")
  buf.write_string("================================\\n")
  buf.write_string("Sequence: ")
  buf.write_string(result.seqid)
  buf.write_string(" (")
  buf.write_string(result.sequence_length.to_string())
  buf.write_string(" nt)\\n")
  buf.write_string("Genes predicted: ")
  buf.write_string(result.num_genes.to_string())
  buf.write_string("\\n")
  buf.write_string("Mean GC content: ")
  buf.write_string((result.mean_gc * 100.0).to_string())
  buf.write_string("%\\n")
  let mut fwd = 0
  let mut rev = 0
  let mut partial = 0
  let mut total_len = 0
  for gene in result.genes {
    if gene.strand == "+" {
      fwd += 1
    } else {
      rev += 1
    }
    if gene.partial {
      partial += 1
    }
    total_len += gene.length
  }
  buf.write_string("Forward strand: ")
  buf.write_string(fwd.to_string())
  buf.write_string("\\n")
  buf.write_string("Reverse strand: ")
  buf.write_string(rev.to_string())
  buf.write_string("\\n")
  buf.write_string("Partial genes: ")
  buf.write_string(partial.to_string())
  buf.write_string("\\n")
  if result.num_genes > 0 {
    buf.write_string("Mean gene length: ")
    buf.write_string((total_len / result.num_genes).to_string())
    buf.write_string(" nt\\n")
    buf.write_string("Coding density: ")
    buf.write_string((total_len.to_double() / result.sequence_length.to_double() * 100.0).to_string())
    buf.write_string("%\\n")
  }
  buf.to_string()
}
