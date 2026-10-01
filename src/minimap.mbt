///|
/// minimap2-lite: minimizer-based approximate aligner.
///
/// Implements the core pipeline of minimap2:
///   1. Minimizer extraction (windowed minimum of k-mer hashes)
///   2. Minimizer index (hash table from reference)
///   3. Seeding (lookup query minimizers in index)
///   4. Chaining (DP with gap/skip constraints)
///   5. Base-level alignment (banded Smith-Waterman)
///
/// This is a pure-MoonBit reimplementation without SIMD optimization.
/// Algorithm follows Heng Li, "Minimap2: pairwise alignment for nucleotide
/// sequences" (Bioinformatics 2018).

// ============================================================================
// 1. Minimizer extraction
// ============================================================================

///|
/// A minimizer: (hash, position, strand).
/// strand = 0 for forward, 1 for reverse complement.
pub struct Minimizer {
  hash : Int
  pos : Int
  strand : Int
}

///|
/// DNA complement table.
fn dna_complement(c : Char) -> Char {
  match c {
    'A' => 'T'
    'a' => 't'
    'T' => 'A'
    't' => 'a'
    'C' => 'G'
    'c' => 'g'
    'G' => 'C'
    'g' => 'c'
    _ => c
  }
}

///|
/// Reverse complement is defined in seq_utils.mbt

///|
/// Hash a k-mer using a simple rolling hash (2-bit encoding).
/// Returns an integer hash. Non-ACGT characters get hash = -1.
fn kmer_hash(seq : String, start : Int, k : Int) -> Int {
  let mut h = 0
  for i = 0; i < k; i = i + 1 {
    let c = seq[start + i].unsafe_to_char()
    let v = match c {
      'A' | 'a' => 0
      'C' | 'c' => 1
      'G' | 'g' => 2
      'T' | 't' => 3
      _ => return -1
    }
    h = (h << 2) | v
  }
  h
}

///|
/// Hash the reverse complement of a k-mer.
fn kmer_hash_rc(seq : String, start : Int, k : Int) -> Int {
  let mut h = 0
  for i = k - 1; i >= 0; i = i - 1 {
    let c = seq[start + i].unsafe_to_char()
    let v = match c {
      'A' | 'a' => 3 // complement of A is T=3
      'T' | 't' => 0 // complement of T is A=0
      'C' | 'c' => 2 // complement of C is G=2
      'G' | 'g' => 1 // complement of G is C=1
      _ => return -1
    }
    h = (h << 2) | v
  }
  h
}

///|
/// Extract minimizers from a sequence.
/// `k` = k-mer size, `w` = window size.
/// For each window of w consecutive k-mers, select the one with the
/// smallest hash (canonical: min of forward and RC hash).
/// Returns a list of minimizers, sorted by position.
pub fn minimap_minimizers(seq : String, k : Int, w : Int) -> Array[Minimizer] {
  let n = seq.length()
  if n < k {
    return []
  }
  let n_kmers = n - k + 1
  if n_kmers < 1 {
    return []
  }
  // Precompute all k-mer hashes (canonical = min of fwd and rc)
  let hashes = Array::make(n_kmers, (-1, 0)) // (hash, strand)
  for i = 0; i < n_kmers; i = i + 1 {
    let h_fwd = kmer_hash(seq, i, k)
    let h_rc = kmer_hash_rc(seq, i, k)
    if h_fwd < 0 || h_rc < 0 {
      hashes[i] = (-1, 0)
    } else {
      if h_fwd <= h_rc {
        hashes[i] = (h_fwd, 0)
      } else {
        hashes[i] = (h_rc, 1)
      }
    }
  }
  // Sliding window minimum
  let result : Array[Minimizer] = Array(capacity=0)
  let mut last_min_hash = -2
  let mut last_min_pos = -1
  for win_start = 0; win_start <= n_kmers - w; win_start = win_start + 1 {
    let win_end = win_start + w
    let mut best_hash = 0x7FFFFFFF
    let mut best_pos = -1
    let mut best_strand = 0
    for i = win_start; i < win_end && i < n_kmers; i = i + 1 {
      let (h, s) = hashes[i]
      if h >= 0 && h < best_hash {
        best_hash = h
        best_pos = i
        best_strand = s
      }
    }
    if best_pos >= 0 && best_pos != last_min_pos {
      result.push({ hash: best_hash, pos: best_pos, strand: best_strand })
      last_min_hash = best_hash
      last_min_pos = best_pos
    }
  }
  // Handle case where window > n_kmers: take the single best
  if result.length() == 0 && n_kmers >= 1 {
    let mut best_hash = 0x7FFFFFFF
    let mut best_pos = -1
    let mut best_strand = 0
    for i = 0; i < n_kmers; i = i + 1 {
      let (h, s) = hashes[i]
      if h >= 0 && h < best_hash {
        best_hash = h
        best_pos = i
        best_strand = s
      }
    }
    if best_pos >= 0 {
      result.push({ hash: best_hash, pos: best_pos, strand: best_strand })
    }
  }
  result
}

// ============================================================================
// 2. Minimizer index
// ============================================================================

///|
/// A reference sequence in the index.
pub struct MinimapRefSeq {
  pub name : String
  pub seq : String
}

///|
/// Constructor for MinimapRefSeq.
pub fn MinimapRefSeq::new(name : String, seq : String) -> MinimapRefSeq {
  { name, seq }
}

///|
/// A seed match: (query_pos, ref_seq_id, ref_pos, strand, k).
pub struct MinimapSeed {
  q_pos : Int
  r_seq_id : Int
  r_pos : Int
  strand : Int
  k : Int
}

///|
/// Minimizer index: maps minimizer hashes to reference positions.
pub struct MinimapIndex {
  // hash → list of (seq_id, pos, strand)
  table : Map[Int, Array[(Int, Int, Int)]]
  refs : Array[MinimapRefSeq]
  k : Int
  w : Int
}

///|
/// Build a minimizer index from reference sequences.
pub fn minimap_index_build(
  refs : Array[MinimapRefSeq],
  k : Int,
  w : Int,
) -> MinimapIndex {
  let table : Map[Int, Array[(Int, Int, Int)]] = Map([], capacity=0)
  for seq_id = 0; seq_id < refs.length(); seq_id = seq_id + 1 {
    let minimizers = minimap_minimizers(refs[seq_id].seq, k, w)
    for m in minimizers {
      let key = m.hash
      let entry = match table.get(key) {
        Some(arr) => arr
        None => {
          let arr : Array[(Int, Int, Int)] = Array(capacity=0)
          table[key] = arr
          arr
        }
      }
      entry.push((seq_id, m.pos, m.strand))
    }
  }
  { table, refs, k, w }
}

// ============================================================================
// 3. Seeding
// ============================================================================

///|
/// Find seed matches between a query sequence and the index.
pub fn minimap_seed(
  index : MinimapIndex,
  query : String,
) -> Array[MinimapSeed] {
  let q_mins = minimap_minimizers(query, index.k, index.w)
  let seeds : Array[MinimapSeed] = Array(capacity=0)
  for qm in q_mins {
    match index.table.get(qm.hash) {
      Some(arr) => {
        for (seq_id, r_pos, r_strand) in arr {
          // Only keep seeds where strands are compatible
          let strand = if qm.strand == r_strand { 0 } else { 1 }
          seeds.push({
            q_pos: qm.pos,
            r_seq_id: seq_id,
            r_pos,
            strand,
            k: index.k,
          })
        }
      }
      None => continue
    }
  }
  // Sort seeds by reference position then query position
  seeds.sort_by(fn(a, b) {
    if a.r_seq_id != b.r_seq_id {
      a.r_seq_id - b.r_seq_id
    } else if a.r_pos != b.r_pos {
      a.r_pos - b.r_pos
    } else {
      a.q_pos - b.q_pos
    }
  })
  seeds
}

// ============================================================================
// 4. Chaining
// ============================================================================

///|
/// A chain: a list of seeds that form a colinear alignment.
pub struct MinimapChain {
  seeds : Array[MinimapSeed]
  score : Double
  q_start : Int
  q_end : Int
  r_start : Int
  r_end : Int
  r_seq_id : Int
  strand : Int
}

///|
/// Chain seeds using a DP approach.
/// Two seeds (q1,r1) and (q2,r2) can be chained if:
///   - q2 > q1 and r2 > r1 (colinear)
///   - gaps are within max_gap
///   - q2 - q1 and r2 - r1 are approximately equal (consistent indel)
pub fn minimap_chain(
  seeds : Array[MinimapSeed],
  max_gap : Int,
  min_seeds : Int,
) -> Array[MinimapChain] {
  let n = seeds.length()
  if n == 0 {
    return []
  }
  // DP: f[i] = best chain score ending at seed i
  let f = Array::make(n, 0.0)
  let prev = Array::make(n, -1)
  for i = 0; i < n; i = i + 1 {
    f[i] = seeds[i].k.to_double()
    prev[i] = -1
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let si = seeds[i]
      let sj = seeds[j]
      // Must be on same reference sequence and strand
      if si.r_seq_id != sj.r_seq_id || si.strand != sj.strand {
        continue
      }
      // Colinear: q and r both increase
      let dq = sj.q_pos - si.q_pos
      let dr = sj.r_pos - si.r_pos
      if dq <= 0 || dr <= 0 {
        continue
      }
      // Gap constraint
      let gap = (dq - dr).abs()
      if dq > max_gap || dr > max_gap {
        continue
      }
      // Chaining score: reward match, penalize gap
      let chain_score = f[i] + sj.k.to_double() - 0.01 * gap.to_double()
      if chain_score > f[j] {
        f[j] = chain_score
        prev[j] = i
      }
    }
  }
  // Backtrack to find chains
  let used = Array::make(n, false)
  let chains : Array[MinimapChain] = Array(capacity=0)
  // Sort indices by score descending
  let order = Array::makei(n, fn(i) { i })
  order.sort_by(fn(a, b) {
    if f[a] > f[b] { -1 } else if f[a] < f[b] { 1 } else { 0 }
  })
  for idx in order {
    if used[idx] {
      continue
    }
    // Backtrack chain
    let chain_seeds : Array[MinimapSeed] = Array(capacity=0)
    let mut cur = idx
    while cur >= 0 && !used[cur] {
      used[cur] = true
      chain_seeds.push(seeds[cur])
      cur = prev[cur]
    }
    if chain_seeds.length() < min_seeds {
      continue
    }
    // Reverse to get seeds in order
    {
    let n = chain_seeds.length()
    for i = 0; i < n / 2; i = i + 1 {
      let tmp = chain_seeds[i]
      chain_seeds[i] = chain_seeds[n - 1 - i]
      chain_seeds[n - 1 - i] = tmp
    }
  }
    // Compute chain bounds
    let mut q_min = chain_seeds[0].q_pos
    let mut q_max = chain_seeds[0].q_pos + chain_seeds[0].k
    let mut r_min = chain_seeds[0].r_pos
    let mut r_max = chain_seeds[0].r_pos + chain_seeds[0].k
    for s in chain_seeds {
      if s.q_pos < q_min {
        q_min = s.q_pos
      }
      if s.q_pos + s.k > q_max {
        q_max = s.q_pos + s.k
      }
      if s.r_pos < r_min {
        r_min = s.r_pos
      }
      if s.r_pos + s.k > r_max {
        r_max = s.r_pos + s.k
      }
    }
    chains.push({
      seeds: chain_seeds,
      score: f[idx],
      q_start: q_min,
      q_end: q_max,
      r_start: r_min,
      r_end: r_max,
      r_seq_id: chain_seeds[0].r_seq_id,
      strand: chain_seeds[0].strand,
    })
  }
  chains
}

// ============================================================================
// 5. Base-level alignment (banded Smith-Waterman)
// ============================================================================

///|
/// Alignment result: (score, ref_start, query_start, cigar_string).
pub struct MinimapAlignment {
  score : Int
  r_start : Int
  r_end : Int
  q_start : Int
  q_end : Int
  cigar : String  // e.g. "10M2D5M"
  nm : Int        // number of mismatches
}

///|
/// Banded Smith-Waterman alignment.
/// Aligns query[q_start:q_end] to ref[r_start:r_end] within a band.
/// `band_width` = half-width of the band.
fn banded_sw(
  ref_seq : String,
  query_seq : String,
  band_width : Int,
  match_score : Int,
  mismatch_penalty : Int,
  gap_penalty : Int,
) -> MinimapAlignment {
  let m = query_seq.length()
  let n = ref_seq.length()
  if m == 0 || n == 0 {
    return { score: 0, r_start: 0, r_end: 0, q_start: 0, q_end: 0, cigar: "", nm: 0 }
  }
  // DP matrix (banded)
  let w = band_width
  let mut best_score = 0
  let mut best_i = 0
  let mut best_j = 0
  // Use full matrix for simplicity (sequences are short within a chain)
  let dp = Array::make((m + 1) * (n + 1), 0)
  let get = fn(i : Int, j : Int) -> Int {
    dp[i * (n + 1) + j]
  }
  let set_dp = fn(i : Int, j : Int, v : Int) -> Unit {
    dp[i * (n + 1) + j] = v
  }
  for i = 1; i <= m; i = i + 1 {
    let j_start = if i - w < 1 { 1 } else { i - w }
    let j_end = if i + w > n { n } else { i + w }
    for j = j_start; j <= j_end; j = j + 1 {
      let diag = get(i - 1, j - 1) + (
        if query_seq[i - 1] == ref_seq[j - 1] { match_score } else { -mismatch_penalty }
      )
      let up = if i - 1 >= 0 && j >= j_start && j <= j_end {
        get(i - 1, j) - gap_penalty
      } else {
        -999999
      }
      let left = if j - 1 >= j_start && j - 1 <= j_end {
        get(i, j - 1) - gap_penalty
      } else {
        -999999
      }
      let val = if diag >= 0 { diag } else { 0 }
      let val = if up > val { up } else { val }
      let val = if left > val { left } else { val }
      set_dp(i, j, val)
      if val > best_score {
        best_score = val
        best_i = i
        best_j = j
      }
    }
  }
  // Backtrack for CIGAR
  let cigar_ops : Array[Char] = Array(capacity=0)
  let mut i = best_i
  let mut j = best_j
  let mut nm = 0
  while i > 0 && j > 0 && get(i, j) > 0 {
    let cur = get(i, j)
    let diag = get(i - 1, j - 1) + (
      if query_seq[i - 1] == ref_seq[j - 1] { match_score } else { -mismatch_penalty }
    )
    if cur == diag {
      if query_seq[i - 1] != ref_seq[j - 1] {
        nm = nm + 1
      }
      cigar_ops.push('M')
      i = i - 1
      j = j - 1
    } else if cur == get(i - 1, j) - gap_penalty {
      cigar_ops.push('I')
      i = i - 1
    } else {
      cigar_ops.push('D')
      j = j - 1
    }
  }
  {
    let n = cigar_ops.length()
    for i = 0; i < n / 2; i = i + 1 {
      let tmp = cigar_ops[i]
      cigar_ops[i] = cigar_ops[n - 1 - i]
      cigar_ops[n - 1 - i] = tmp
    }
  }
  // Compress CIGAR
  let cigar = compress_cigar(cigar_ops)
  {
    score: best_score,
    r_start: j,
    r_end: best_j,
    q_start: i,
    q_end: best_i,
    cigar,
    nm,
  }
}

///|
/// Compress CIGAR operations into run-length encoded string.
/// e.g. ['M','M','M','D','M','M'] → "3M1D2M"
fn compress_cigar(ops : Array[Char]) -> String {
  if ops.length() == 0 {
    return ""
  }
  let buf : StringBuilder = StringBuilder()
  let mut count = 1
  let mut last = ops[0]
  for i = 1; i < ops.length(); i = i + 1 {
    if ops[i] == last {
      count = count + 1
    } else {
      buf.write_string(count.to_string())
      buf.write_char(last)
      count = 1
      last = ops[i]
    }
  }
  buf.write_string(count.to_string())
  buf.write_char(last)
  buf.to_string()
}

// ============================================================================
// 6. End-to-end alignment
// ============================================================================

///|
/// Full alignment result for one query.
pub struct MinimapResult {
  query_name : String
  ref_name : String
  ref_start : Int     // 0-based
  ref_end : Int       // exclusive
  query_start : Int   // 0-based
  query_end : Int     // exclusive
  strand : Int        // 0=forward, 1=reverse
  score : Int
  cigar : String
  nm : Int            // edit distance
  mapq : Int          // mapping quality
}

///|
/// Align a query sequence against the index.
/// This is the main entry point: index → seed → chain → align.
pub fn minimap_align(
  index : MinimapIndex,
  query : String,
  query_name : String,
) -> Array[MinimapResult] {
  let seeds = minimap_seed(index, query)
  let chains = minimap_chain(seeds, index.k * 20, 1)
  let results : Array[MinimapResult] = Array(capacity=0)
  for chain in chains {
    let ref_seq = index.refs[chain.r_seq_id].seq
    let ref_name = index.refs[chain.r_seq_id].name
    // Extract regions for base-level alignment
    let r_pad = index.k * 2
    let r_start = if chain.r_start - r_pad < 0 { 0 } else { chain.r_start - r_pad }
    let r_end = if chain.r_end + r_pad > ref_seq.length() { ref_seq.length() } else { chain.r_end + r_pad }
    let q_pad = index.k * 2
    let q_start = if chain.q_start - q_pad < 0 { 0 } else { chain.q_start - q_pad }
    let q_end = if chain.q_end + q_pad > query.length() { query.length() } else { chain.q_end + q_pad }
    let ref_region = ref_seq[r_start:r_end].to_string()
    let query_region = if chain.strand == 0 {
      query[q_start:q_end].to_string()
    } else {
      reverse_complement(query[q_start:q_end].to_string())
    }
    let aln = banded_sw(
      ref_region, query_region, index.k * 5,
      2, 1, 1,
    )
    if aln.score > 0 {
      results.push({
        query_name,
        ref_name,
        ref_start: r_start + aln.r_start,
        ref_end: r_start + aln.r_end,
        query_start: q_start + aln.q_start,
        query_end: q_start + aln.q_end,
        strand: chain.strand,
        score: aln.score,
        cigar: aln.cigar,
        nm: aln.nm,
        mapq: compute_mapq(aln.score, chains.length()),
      })
    }
  }
  // Sort by score descending
  results.sort_by(fn(a, b) {
    if a.score > b.score { -1 } else if a.score < b.score { 1 } else { 0 }
  })
  results
}

///|
/// Compute mapping quality (simplified).
fn compute_mapq(score : Int, n_chains : Int) -> Int {
  let base = score * 5
  let penalty = if n_chains > 1 { (n_chains - 1) * 10 } else { 0 }
  let mq = base - penalty
  if mq < 0 { 0 } else if mq > 60 { 60 } else { mq }
}

///|
/// Convert alignment result to SAM line.
pub fn minimap_to_sam(result : MinimapResult, query_seq : String) -> String {
  let flag = if result.strand == 1 { 16 } else { 0 }
  let rname = result.ref_name
  let pos = result.ref_start + 1 // SAM is 1-based
  let mapq = result.mapq
  let cigar = if result.cigar.length() == 0 { "*" } else { result.cigar }
  let seq = if result.strand == 1 {
    reverse_complement(query_seq)
  } else {
    query_seq
  }
  let qual = "*"
  let nm_tag = "NM:i:" + result.nm.to_string()
  let qname = result.query_name
  let line = qname + "\t" + flag.to_string() + "\t" + rname + "\t" + pos.to_string() + "\t" + mapq.to_string() + "\t" + cigar + "\t*\t0\t0\t" + seq + "\t" + qual + "\t" + nm_tag
  line
}

///|
/// Align multiple queries against the index and produce SAM output.
pub fn minimap_align_batch(
  index : MinimapIndex,
  queries : Array[(String, String)], // (name, seq)
) -> String {
  let buf : StringBuilder = StringBuilder()
  // SAM header
  for ref in index.refs {
    buf.write_string("@SQ\tSN:" + ref.name + "\tLN:" + ref.seq.length().to_string() + "\n")
  }
  for (name, seq) in queries {
    let results = minimap_align(index, seq, name)
    if results.length() == 0 {
      // Unmapped
      buf.write_string(name + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + seq + "\t*\n")
    } else {
      for r in results {
        buf.write_string(minimap_to_sam(r, seq) + "\n")
      }
    }
  }
  buf.to_string()
}
