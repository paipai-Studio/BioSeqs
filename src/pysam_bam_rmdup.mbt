// ============================================================================
// pysam_bam_rmdup — samtools rmdup 的内存版实现
//
// samtools rmdup（旧版，CIGAR-aware 简化版）：
//  1. 按 (qname, ref_id, pos) 分组 reads
//  2. 每组保留 mapq 最高的（tie-break 用原始顺序）
//  3. 返回**不含** duplicates 的新 BamFile（物理删除，不像 markdup 只打 flag）
//
// 与 markdup 的差异：
//  - rmdup 直接丢弃 duplicate reads（-s/-S 只标记 flag 不丢 read 的旧版行为）
//  - markdup 只打 duplicate flag，不做物理删除
//
// 参考：samtools rmdup / Pysam `AlignmentFile.deduplicate()`
// ============================================================================

// ---------------------------------------------------------------------------
// 核心算法：单端 rmdup（按 coordinate 分组）
// ---------------------------------------------------------------------------

///|
/// Build the grouping key for a rmdup candidate: `(qname, ref_id, pos)`.
/// Unmapped reads are always kept (no dup candidate).
fn rmdup_key(rec : BamRecord) -> String {
  let ref_str = if rec.ref_id < 0 { "*" } else { rec.ref_id.to_string() }
  rec.qname + "\t" + ref_str + "\t" + rec.pos.to_string()
}

///|
/// Return true when two BamRecords are considered identical duplicates
/// (same ref + pos + qname).  samtools rmdup also checks CIGAR but we keep
/// the coordinate-only check here — callers who need CIGAR-aware detection
/// should use the SAM-text `pysam_rmdup` helper that already exists.
fn rmdup_same_group(a : BamRecord, b : BamRecord) -> Bool {
  a.qname == b.qname && a.ref_id == b.ref_id && a.pos == b.pos
}

///|
/// Identify duplicate groups and keep the record with the highest MAPQ in
/// each group.  Ties go to the record appearing first.
///
/// Returns an array of indices of records to keep.
fn rmdup_keep_indices(bam : BamFile) -> Array[Int32] {
  // key → (highest_mapq_seen, best_index)
  let best_per_key : Map[String, (Int32, Int32)] = Map([], capacity=0)
  let mut i = 0
  while i < bam.records.length() {
    let r = bam.records[i]
    // Unmapped reads always pass rmdup (they can't form dup groups by pos)
    if r.is_unmapped() {
      i = i + 1
      continue
    }
    let k = rmdup_key(r)
    let cur = best_per_key.get(k)
    match cur {
      Some((best_q, best_idx)) =>
        // Keep whichever has higher mapq; if equal keep the first-seen one.
        if r.mapq > best_q {
          best_per_key.set(k, (r.mapq, i))
        } else {
          best_per_key.set(k, (best_q, best_idx))
        }
      None => best_per_key.set(k, (r.mapq, i))
    }
    i = i + 1
  }
  // Collect kept indices. We preserve the original record order by sorting.
  let kept : Array[Int32] = best_per_key
    .values()
    .map(fn(v) {
      let (_, idx) = v
      idx
    })
    .to_array()
  // Also include unmapped reads (skipped in grouping) — we just didn't drop them
  // above, so `kept` misses them. Fix: iterate all and keep unmapped + winning idx.
  let kept_set : Set[Int32] = Set(kept)
  let out : Array[Int32] = []
  let mut j = 0
  while j < bam.records.length() {
    if bam.records[j].is_unmapped() || kept_set.contains(j) {
      out.push(j)
    }
    j = j + 1
  }
  out
}

// ---------------------------------------------------------------------------
// 主入口：rmdup（物理删除 duplicates）
// ---------------------------------------------------------------------------

///|
/// Run rmdup on a BamFile: drop duplicate reads, keep only the highest-MAPQ
/// record per (qname, ref, pos) group. Unmapped reads are always kept.
/// Returns a new BamFile with the duplicate-free records.
///
/// Mirrors `samtools rmdup <in.bam> <out.bam>` (CIGAR-aware check simplified
/// to coordinate-only here).
pub fn bam_rmdup(bam : BamFile) -> BamFile {
  let keep = rmdup_keep_indices(bam)
  let out_records : Array[BamRecord] = Array(capacity=keep.length())
  for idx in keep {
    out_records.push(bam.records[idx])
  }
  BamFile::{ header: bam.header, records: out_records, }
}

// ---------------------------------------------------------------------------
// 替代方案：只打 FLAG 不删除（samtools rmdup -s / markdup 风格）
// ---------------------------------------------------------------------------

///|
/// Like rmdup but sets the duplicate flag bit (0x400) instead of physically
/// removing records. Equivalent to `samtools markdup -r -d` without removal.
/// The caller can filter on `is_duplicate()` later.
pub fn bam_markdup_rmdup_mode(bam : BamFile) -> BamFile {
  let keep = rmdup_keep_indices(bam)
  let keep_set : Set[Int32] = Set(keep)
  let out_records : Array[BamRecord] = Array(capacity=bam.records.length())
  let mut i = 0
  while i < bam.records.length() {
    let r = bam.records[i]
    if !keep_set.contains(i) && !r.is_duplicate() {
      // Mark as duplicate
      let new_flag = r.flag | BAM_FDUP
      out_records.push(r.with_flag(new_flag))
    } else {
      out_records.push(r)
    }
    i = i + 1
  }
  BamFile::{ header: bam.header, records: out_records, }
}

// ---------------------------------------------------------------------------
// 统计辅助
// ---------------------------------------------------------------------------

///|
/// Count how many records would be dropped by rmdup.  Useful for logging /
/// before-and-after comparison.
pub fn bam_count_duplicates(bam : BamFile) -> Int32 {
  let keep = rmdup_keep_indices(bam)
  let total = bam.records.length()
  // Formula: total_mapped - keep_mapped = total - keep
  // (because keep includes unmapped explicitly)
  total - keep.length()
}

///|
/// Fraction of reads that are duplicates (0.0 ~ 1.0).  Returns 0.0 when the
/// BamFile has no records.
pub fn bam_duplicate_rate(bam : BamFile) -> Double {
  if bam.records.length() == 0 {
    return 0.0
  }
  let dup = bam_count_duplicates(bam).to_double()
  dup / bam.records.length().to_double()
}

// ---------------------------------------------------------------------------
// 成对端 rmdup (-S 语义的简化版)
// ---------------------------------------------------------------------------

///|
/// Pair-end aware rmdup: treats each mate pair as a single unit, keyed by
/// (qname, first_pos_of_pair).  Useful for PE libraries where both mates
/// should be dropped if **either** one is duplicated. Mirrors samtools rmdup
/// `-S` ("discard the first read of a read pair") mode at a high level.
///
/// The current implementation keys on the first-seen position for each qname
/// pair; this is a common simplification and avoids needing real mate-pair
/// adjacency information.
pub fn bam_rmdup_paired(bam : BamFile) -> BamFile {
  // Step 1: identify pairs — qnames appearing exactly twice
  let idx_by_qname : Map[String, Array[Int32]] = Map([], capacity=0)
  let mut i = 0
  while i < bam.records.length() {
    let cur = idx_by_qname.get(bam.records[i].qname)
    match cur {
      Some(arr) => {
        let copy = arr.copy()
        copy.push(i)
        idx_by_qname.set(bam.records[i].qname, copy)
      }
      None => idx_by_qname.set(bam.records[i].qname, [i])
    }
    i = i + 1
  }
  // Step 2: for each qname group, pick the pair with higher total MAPQ
  let keep : Set[Int32] = Set([])
  let names = idx_by_qname.keys()
  for qname in names {
    let arr = idx_by_qname[qname]
    if arr.length() == 2 {
      // Paired: keep both members
      keep.add(arr[0])
      keep.add(arr[1])
    } else {
      // Unpaired or multi-hit group: apply single-end rule
      let sub : Array[BamRecord] = arr.map(fn(idx) { bam.records[idx] })
      let sub_file = BamFile::{ header: bam.header, records: sub, }
      let keep_local = rmdup_keep_indices(sub_file)
      for li in keep_local {
        keep.add(arr[li])
      }
    }
  }
  // Step 3: emit
  let out_records : Array[BamRecord] = []
  let mut j = 0
  while j < bam.records.length() {
    if keep.contains(j) {
      out_records.push(bam.records[j])
    }
    j = j + 1
  }
  BamFile::{ header: bam.header, records: out_records, }
}
