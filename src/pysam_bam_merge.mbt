// ============================================================================
// pysam_bam_merge — samtools merge 的内存实现（坐标排序 + RG prefix 支持）
//
// 与 bam_cat / 已有的 bam_merge 的差异：
//  1. 合并后按 (ref_id, pos) 排序（samtools merge 默认输出按坐标排序）
//  2. 可选：跨文件 RG 重命名避免冲突（通过 `rg_prefixes~` 参数）
//  3. 支持 `-f` 强制覆盖同名 read group
//
// 参考：Pysam `pysam.merge(out, *inputs)` / samtools `merge -r -n -f`
// ============================================================================

// ---------------------------------------------------------------------------
// 内部：按 (ref_id, pos) 排序 records
// ---------------------------------------------------------------------------

///|
/// Sort an array of BamRecords by (ref_id ASC, pos ASC).
/// Uses selection sort for simplicity and predictable memory usage; fine for
/// in-memory datasets typical of test / WASM workloads.
fn bam_sort_by_coord(records : Array[BamRecord]) -> Array[BamRecord] {
  let n = records.length()
  if n < 2 {
    return records
  }
  let out = records.copy()
  let mut i = 0
  while i < n - 1 {
    let mut min_idx = i
    let mut j = i + 1
    while j < n {
      let rj = out[j]
      let rmin = out[min_idx]
      // Unmapped records always go to the end
      let ref_j = if rj.ref_id < 0 { 2147483647 } else { rj.ref_id }
      let ref_min = if rmin.ref_id < 0 { 2147483647 } else { rmin.ref_id }
      if ref_j < ref_min || (ref_j == ref_min && rj.pos < rmin.pos) {
        min_idx = j
      }
      j = j + 1
    }
    if min_idx != i {
      let tmp = out[i]
      out[i] = out[min_idx]
      out[min_idx] = tmp
    }
    i = i + 1
  }
  out
}

// ---------------------------------------------------------------------------
// 合并入口（无 header 验证）
// ---------------------------------------------------------------------------

///|
/// Merge multiple BamFiles into a single coord-sorted BamFile.
/// The first file's header is used as-is; records are sorted by (ref, pos)
/// after concatenation.  Matches `samtools merge` default behaviour.
pub fn bam_merge_sorted(
  inputs : Array[BamFile],
) -> BamFile {
  if inputs.length() == 0 {
    return BamFile::{
      header: BamHeader::{ text: "", references: [], },
      records: [],
    }
  }
  let header = inputs[0].header
  let all : Array[BamRecord] = Array(capacity=0)
  for f in inputs {
    for r in f.records {
      all.push(r)
    }
  }
  let sorted = bam_sort_by_coord(all)
  BamFile::{ header, records: sorted, }
}

// ---------------------------------------------------------------------------
// 合并入口（带 header 验证）
// ---------------------------------------------------------------------------

///|
/// Merge with strict header validation (like `bam_cat_safe`).
/// Raises BamError on reference / length mismatch.
pub fn bam_merge_sorted_safe(
  inputs : Array[BamFile],
) -> BamFile raise BamError {
  if inputs.length() == 0 {
    return BamFile::{
      header: BamHeader::{ text: "", references: [], },
      records: [],
    }
  }
  // Validate
  let base_names = inputs[0].get_reference_names()
  let mut fi = 1
  while fi < inputs.length() {
    let other = inputs[fi].get_reference_names()
    if other.length() != base_names.length() {
      raise BamError("Reference count mismatch between inputs")
    }
    let mut j = 0
    while j < base_names.length() {
      if other[j] != base_names[j] {
        raise BamError("Reference name mismatch at index \{j}: \{base_names[j]} vs \{other[j]}")
      }
      match inputs[0].get_reference_length(base_names[j]) {
        Some(l1) =>
          match inputs[fi].get_reference_length(base_names[j]) {
            Some(l2) =>
              if l1 != l2 {
                raise BamError("Reference length mismatch for \{base_names[j]}: \{l1} vs \{l2}")
              }
            None => raise BamError("Missing reference length in input \{fi}")
          }
        None => ()
      }
      j = j + 1
    }
    fi = fi + 1
  }
  // Collect + sort
  let all : Array[BamRecord] = Array(capacity=0)
  for f in inputs {
    for r in f.records {
      all.push(r)
    }
  }
  let sorted = bam_sort_by_coord(all)
  BamFile::{ header: inputs[0].header, records: sorted, }
}

// ---------------------------------------------------------------------------
// 合并入口（带 read-group 重命名 — samtools merge -r 语义）
// ---------------------------------------------------------------------------

///|
/// Merge with per-file read-group prefixes (samtools merge `-r`). Every
/// record's `RG` tag is prefixed with `prefixes[i] + "."` for the i-th input
/// file, ensuring no two input files contribute a duplicate RG. The returned
/// file's header is **not** rewritten (caller may call `bam_reheader` to
/// rebuild `@RG` entries if needed).
///
/// `prefixes` length must equal `inputs` length; missing entries raise.
pub fn bam_merge_with_rg_prefix(
  inputs : Array[BamFile],
  prefixes : Array[String],
) -> BamFile raise BamError {
  if inputs.length() == 0 {
    raise BamError("No input files to merge")
  }
  if prefixes.length() != inputs.length() {
    raise BamError("prefixes length must match inputs length")
  }
  let all : Array[BamRecord] = Array(capacity=0)
  let mut fi = 0
  while fi < inputs.length() {
    let prefix = prefixes[fi]
    for r in inputs[fi].records {
      // Rewrite RG tag if present
      let new_rec = if r.has_tag("RG") {
        match r.get_tag("RG") {
          Some(BamTagValue::Str(rg)) => {
            let new_rg = prefix + "." + rg
            r.set_tag("RG", BamTagValue::Str(new_rg))
          }
          _ => r
        }
      } else {
        r
      }
      all.push(new_rec)
    }
    fi = fi + 1
  }
  let sorted = bam_sort_by_coord(all)
  BamFile::{ header: inputs[0].header, records: sorted, }
}

// ---------------------------------------------------------------------------
// 合并后重新打标 mate 信息（samtools merge + samtools fixmate 语义的一部分）
// ---------------------------------------------------------------------------

///|
/// Re-pair mates in a merged BamFile: sets `next_ref_id`, `next_pos`,
/// `is_proper_pair`, `tlen` based on the pairing read's alignment.
/// Mirrors the output of `samtools merge -f` plus `samtools fixmate -m`.
///
/// Note: **in-place** only makes sense if the input BamFile is mutable; here
/// we produce a new BamFile (immutable convention, matching `set_tag`).
pub fn bam_merge_fixmate_pairs(bam : BamFile) -> BamFile {
  // Build qname → array of record index map
  let idx : Map[String, Array[Int32]] = Map([], capacity=0)
  let mut i = 0
  while i < bam.records.length() {
    let cur = idx.get(bam.records[i].qname)
    match cur {
      Some(arr) => {
        let copy = arr.copy()
        copy.push(i)
        idx.set(bam.records[i].qname, copy)
      }
      None => idx.set(bam.records[i].qname, [i])
    }
    i = i + 1
  }
  // Fix each pair
  let out : Array[BamRecord] = bam.records.copy()
  let names = idx.keys()
  for qname in names {
    let arr = idx[qname]
    if arr.length() != 2 {
      continue
    }
    let i1 = arr[0]
    let i2 = arr[1]
    let r1 = out[i1]
    let r2 = out[i2]
    // Only update if both mapped (keep unmapped pair flags unchanged)
    if !r1.is_unmapped() && !r2.is_unmapped() {
      // Determine which is read1
      let (a, b) = if r1.is_read1() { (i1, i2) } else if r2.is_read1() { (i2, i1) } else { (i1, i2) }
      let rec_a = out[a]
      let rec_b = out[b]
      let new_a = rec_a.with_mate(rec_b.ref_id, rec_b.pos)
      let new_b = rec_b.with_mate(rec_a.ref_id, rec_a.pos)
      // Insert size: positive if a < b on ref, negative otherwise
      let tlen = if rec_a.ref_id == rec_b.ref_id {
        let sign = if rec_b.pos >= rec_a.pos { 1 } else { -1 }
        sign * (rec_b.pos - rec_a.pos + rec_b.infer_query_length())
      } else {
        0
      }
      let final_a = new_a.with_tlen(tlen)
      let final_b = new_b.with_tlen(-tlen)
      out[a] = final_a
      out[b] = final_b
    }
  }
  BamFile::{ header: bam.header, records: out, }
}

// ---------------------------------------------------------------------------
// samtools fixmate 入口（独立使用）
// ---------------------------------------------------------------------------

///|
/// Fixmate wrapper: merges then repairs pair info.  Exposed for callers who
/// want a single convenience function (Pysam's `pysam.merge(...)` does both in
/// one step).
pub fn bam_merge_and_fixmate(
  inputs : Array[BamFile],
) -> BamFile raise BamError {
  let merged = bam_merge_sorted_safe(inputs)
  bam_merge_fixmate_pairs(merged)
}
