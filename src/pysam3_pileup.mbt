// Pysam third-round: Pileup engine (items 1-5)
// pileup, PileupRead, PileupColumn, mpileup, per-position base counts

///|
pub struct Pysam3PileupRead {
  query_pos : Int32 // position in read (0-based)
  ref_pos : Int32 // reference position (0-based)
  indel : Int32 // indel length: 0=no indel, >0=insertion, <0=deletion
  level : Int32 // stacking level for visualization
  is_del : Bool // whether this read has a deletion at this position
  is_refskip : Bool // whether this is a reference skip (N in CIGAR)
  base : Int32 // base as ASCII code
  qual : Int32 // base quality
  mapq : Int32 // mapping quality
  strand : Int32 // 0=forward, 1=reverse
}

///|
pub struct Pysam3PileupColumn {
  ref_id : Int32 // reference sequence ID
  pos : Int32 // reference position (0-based)
  depth : Int32 // number of reads covering this position
  reads : Array[Pysam3PileupRead]
}

///|
pub struct Pysam3BaseCounts {
  a : Int32
  c : Int32
  g : Int32
  t : Int32
  ins : Int32 // insertions
  del : Int32 // deletions
  n : Int32 // N bases
  ref_skip : Int32 // reference skips
}

///|
pub fn pysam3_base_counts_empty() -> Pysam3BaseCounts {
  { a: 0, c: 0, g: 0, t: 0, ins: 0, del: 0, n: 0, ref_skip: 0, }
}

///|
pub fn pysam3_base_counts_add(
  counts : Pysam3BaseCounts,
  base : Int32,
) -> Pysam3BaseCounts {
  match base {
    65 | 97 => { ..counts, a: counts.a + 1, }
    67 | 99 => { ..counts, c: counts.c + 1, }
    71 | 103 => { ..counts, g: counts.g + 1, }
    84 | 116 => { ..counts, t: counts.t + 1, }
    78 | 110 => { ..counts, n: counts.n + 1, }
    _ => counts
  }
}

///|
pub fn pysam3_base_counts_total(counts : Pysam3BaseCounts) -> Int32 {
  counts.a + counts.c + counts.g + counts.t + counts.n
}

///|
pub fn pysam3_base_counts_majority(counts : Pysam3BaseCounts) -> Int32 {
  let pairs = [(65, counts.a), (67, counts.c), (71, counts.g), (84, counts.t)]
  let mut best_base = 78
  let mut best_count = counts.n
  for i in 0..<4 {
    let (b, c) = pairs[i]
    if c > best_count {
      best_base = b
      best_count = c
    }
  }
  best_base
}

///|
pub fn pysam3_base_counts_to_string(counts : Pysam3PileupColumn) -> String {
  let sb = StringBuilder()
  let _ = sb.write("pos=")
  let _ = sb.write(counts.pos.to_string())
  let _ = sb.write(" depth=")
  let _ = sb.write(counts.depth.to_string())
  let _ = sb.write(" reads=")
  let _ = sb.write(counts.reads.length().to_string())
  sb.to_string()
}

// --- Pileup engine ---

///|
pub struct Pysam3PileupAlign {
  ref_pos : Int32 // reference start position (0-based)
  cigar_ops : Array[(Int32, Int32)] // (op_code, length): 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
  seq : String // read sequence
  qual : String // base qualities (Phred+33)
  mapq : Int32 // mapping quality
  is_reverse : Bool // is reverse strand
}

///|
pub fn pysam3_pileup_column(
  alignments : Array[Pysam3PileupAlign],
  target_pos : Int32,
) -> Pysam3PileupColumn {
  let reads : Array[Pysam3PileupRead] = []
  let mut level = 0
  for i in 0..<alignments.length() {
    let aln = alignments[i]
    let mut ref_idx = aln.ref_pos
    let mut seq_idx = 0
    let mut found = false
    for j in 0..<aln.cigar_ops.length() {
      let (op, length) = aln.cigar_ops[j]
      if op == 0 || op == 7 || op == 8 {
        // M, =, X: match/mismatch
        for k in 0..<length {
          if ref_idx + k == target_pos {
            let base = if seq_idx + k < aln.seq.length() {
              aln.seq[seq_idx + k].to_int()
            } else {
              78
            }
            let q = if seq_idx + k < aln.qual.length() {
              aln.qual[seq_idx + k].to_int() - 33
            } else {
              0
            }
            reads.push({
              query_pos: seq_idx + k,
              ref_pos: target_pos,
              indel: 0,
              level,
              is_del: false,
              is_refskip: false,
              base,
              qual: q,
              mapq: aln.mapq,
              strand: if aln.is_reverse {
                1
              } else {
                0
              },
            })
            found = true
            level = level + 1
          }
        }
        ref_idx = ref_idx + length
        seq_idx = seq_idx + length
      } else if op == 1 {
        // I: insertion
        if ref_idx == target_pos && found {
          // mark insertion on the last read
          if reads.length() > 0 {
            let last = reads[reads.length() - 1]
            reads[reads.length() - 1] = { ..last, indel: length, }
          }
        }
        seq_idx = seq_idx + length
      } else if op == 2 {
        // D: deletion
        for k in 0..<length {
          if ref_idx + k == target_pos {
            reads.push({
              query_pos: seq_idx,
              ref_pos: target_pos,
              indel: -length,
              level,
              is_del: true,
              is_refskip: false,
              base: 0,
              qual: 0,
              mapq: aln.mapq,
              strand: if aln.is_reverse {
                1
              } else {
                0
              },
            })
            level = level + 1
          }
        }
        ref_idx = ref_idx + length
      } else if op == 3 {
        // N: reference skip
        for k in 0..<length {
          if ref_idx + k == target_pos {
            reads.push({
              query_pos: seq_idx,
              ref_pos: target_pos,
              indel: 0,
              level,
              is_del: false,
              is_refskip: true,
              base: 0,
              qual: 0,
              mapq: aln.mapq,
              strand: if aln.is_reverse {
                1
              } else {
                0
              },
            })
            level = level + 1
          }
        }
        ref_idx = ref_idx + length
      } else if op == 4 {
        // S: soft clip
        seq_idx = seq_idx + length
      }
    }
  }
  { ref_id: 0, pos: target_pos, depth: reads.length(), reads, }
}

///|
pub fn pysam3_pileup(
  alignments : Array[Pysam3PileupAlign],
  start : Int32,
  end_ : Int32,
) -> Array[Pysam3PileupColumn] {
  let columns : Array[Pysam3PileupColumn] = []
  for pos in start..<end_ {
    columns.push(pysam3_pileup_column(alignments, pos))
  }
  columns
}

///|
pub fn pysam3_pileup_base_counts(
  column : Pysam3PileupColumn,
) -> Pysam3BaseCounts {
  let mut counts = pysam3_base_counts_empty()
  for i in 0..<column.reads.length() {
    let r = column.reads[i]
    if r.is_del {
      counts = { ..counts, del: counts.del + 1, }
    } else if r.is_refskip {
      counts = { ..counts, ref_skip: counts.ref_skip + 1, }
    } else if r.indel > 0 {
      counts = { ..counts, ins: counts.ins + 1, }
    } else if r.indel < 0 {
      counts = { ..counts, del: counts.del + 1, }
    } else {
      counts = pysam3_base_counts_add(counts, r.base)
    }
  }
  counts
}

// --- Multi-sample pileup (mpileup) ---

///|
pub struct Pysam3MPileupResult {
  pos : Int32
  samples : Array[(String, Pysam3BaseCounts)] // (sample_name, counts)
}

///|
pub fn pysam3_mpileup(
  samples : Array[(String, Array[Pysam3PileupAlign])],
  start : Int32,
  end_ : Int32,
) -> Array[Pysam3MPileupResult] {
  let results : Array[Pysam3MPileupResult] = []
  for pos in start..<end_ {
    let sample_counts : Array[(String, Pysam3BaseCounts)] = []
    for s in 0..<samples.length() {
      let (name, alns) = samples[s]
      let col = pysam3_pileup_column(alns, pos)
      let counts = pysam3_pileup_base_counts(col)
      sample_counts.push((name, counts))
    }
    results.push({ pos, samples: sample_counts, })
  }
  results
}

///|
pub fn pysam3_pileup_depth_array(
  columns : Array[Pysam3PileupColumn],
) -> Array[Int32] {
  let depths : Array[Int32] = []
  for i in 0..<columns.length() {
    depths.push(columns[i].depth)
  }
  depths
}

///|
pub fn pysam3_pileup_consensus(
  columns : Array[Pysam3PileupColumn],
  min_depth : Int32,
) -> String {
  let sb = StringBuilder()
  for i in 0..<columns.length() {
    if columns[i].depth < min_depth {
      let _ = sb.write_char('N')
    } else {
      let counts = pysam3_pileup_base_counts(columns[i])
      let base = pysam3_base_counts_majority(counts)
      let _ = sb.write_char(Int32::unsafe_to_char(base))
    }
  }
  sb.to_string()
}

///|
pub fn pysam3_pileup_align_dummy() -> Pysam3PileupAlign {
  { ref_pos: 0, cigar_ops: [], seq: "", qual: "", mapq: 0, is_reverse: false, }
}
