// ============================================================================
// pysam_vcf_merge — bcftools merge 的内存实现（VCF 纵向合并）
//
// 核心操作：
//  1. 拼接多份 VCF 的 records（可选按 (chrom, pos) 重新排序）
//  2. 合并 meta headers（contigs / INFO / FILTER / FORMAT 去重）
//  3. 合并 SAMPLE 列（同名样本加后缀或重命名）
//
// 参考：bcftools merge -m none / -O z / --missing-to-ref
// ============================================================================

///|
/// Public helper — merge arrays of VcfMetaHeader deduplicating by
/// `(key, value)` pair. Later duplicates (across inputs) are dropped;
/// the first occurrence wins. Also merges contigs so only unique ones
/// survive.
pub fn merge_vcf_meta(inputs : Array[VcfFile]) -> Array[VcfMetaHeader] {
  let out : Array[VcfMetaHeader] = Array(capacity=0)
  for vcf in inputs {
    for mh in vcf.meta_headers {
      let mut dup = false
      for existing in out {
        if existing.key == mh.key && existing.value == mh.value {
          dup = true
          break
        }
      }
      if !dup { out.push(mh) }
    }
  }
  out
}

///|
/// Merge the #CHROM header columns, joining SAMPLE names from every input.
/// Returns the new columns array.
///
/// If `sample_suffixes` is non-empty, its entries are appended to SAMPLE
/// names from the corresponding input to disambiguate duplicates.
pub fn merge_vcf_columns(
  inputs : Array[VcfFile],
  sample_suffixes~ : Array[String] = [],
) -> Array[String] {
  // Base columns (always from first input, or fall back to a synthetic
  // default when inputs have no fixed header)
  let base : Array[String] = if inputs.length() > 0 {
    match inputs[0].fixed_header {
      Some(fh) => fh.columns.copy()
      None => ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    }
  } else {
    ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
  }
  // If no input has FORMAT/SAMPLE, stop early — nothing to merge
  let has_samples = {
    let mut any = false
    for vcf in inputs {
      match vcf.fixed_header {
        Some(fh) => {
          let mut i = 0
          while i < fh.columns.length() {
            if fh.columns[i] == "FORMAT" { any = true; break }
            i = i + 1
          }
        }
        None => ()
      }
      if any { break }
    }
    any
  }
  if !has_samples { return base }

  // Collect SAMPLE names with optional suffixes
  let out : Array[String] = base.copy()
  let suffix_defined = sample_suffixes.length() > 0
  let seen_samples : Array[String] = Array(capacity=0)
  let mut vi = 0
  while vi < inputs.length() {
    let suffix = if suffix_defined && vi < sample_suffixes.length() {
      sample_suffixes[vi]
    } else {
      ""
    }
    match inputs[vi].fixed_header {
      Some(fh) => {
        // 找 FORMAT 的索引
        let mut fmt_idx : Int32 = -1
        let mut ci = 0
        while ci < fh.columns.length() {
          if fh.columns[ci] == "FORMAT" { fmt_idx = ci; break }
          ci = ci + 1
        }
        if fmt_idx >= 0 {
          let mut si = fmt_idx + 1
          while si < fh.columns.length() {
            let orig = fh.columns[si]
            let name = if suffix.length() > 0 {
              orig + "_" + suffix
            } else {
              orig
            }
            let mut final_name = name
            let mut dup = false
            for s in seen_samples {
              if s == final_name { dup = true; break }
            }
            if dup {
              final_name = name + "_" + vi.to_string()
            }
            seen_samples.push(final_name)
            si = si + 1
          }
        }
      }
      None => ()
    }
    vi = vi + 1
  }
  for s in seen_samples { out.push(s) }
  out
}

// ---------------------------------------------------------------------------
// Records 合并
// ---------------------------------------------------------------------------

///|
/// Concatenate records from every input VCF, preserving input order.
/// Does **not** sort — callers that need coordinate-sorted output must
/// sort afterwards (see `vcf_merge_sorted`).
pub fn vcf_merge_unsorted(inputs : Array[VcfFile]) -> VcfFile {
  let all_records : Array[VcfRecord] = Array(capacity=0)
  let mut vi = 0
  while vi < inputs.length() {
    let rs = inputs[vi].records
    let mut ri = 0
    while ri < rs.length() {
      all_records.push(rs[ri])
      ri = ri + 1
    }
    vi = vi + 1
  }
  VcfFile::{
    meta_headers: merge_vcf_meta(inputs),
    fixed_header: Some(VcfFixedHeader::{ columns: merge_vcf_columns(inputs) }),
    records: all_records,
  }
}

///|
/// Stable merge + sort by (chrom, pos). This mirrors `bcftools merge`
/// default behaviour.
pub fn vcf_merge_sorted(inputs : Array[VcfFile]) -> VcfFile {
  let merged = vcf_merge_unsorted(inputs)
  let records = merged.records
  if records.length() <= 1 {
    return merged
  }
  // 按 chrom 分桶
  let grouped : Map[String, Array[VcfRecord]] = Map([], capacity=0)
  let mut ri = 0
  while ri < records.length() {
    let r = records[ri]
    match grouped.get(r.chrom) {
      Some(list) => {
        let copy = list.copy()
        copy.push(r)
        grouped.set(r.chrom, copy)
      }
      None => grouped.set(r.chrom, [r])
    }
    ri = ri + 1
  }
  // Chromosome 顺序：取 inputs 里第一次出现的顺序
  let chrom_order : Array[String] = Array(capacity=0)
  for vcf in inputs {
    let vrs = vcf.records
    let mut vri = 0
    while vri < vrs.length() {
      let c = vrs[vri].chrom
      let mut seen = false
      for sc in chrom_order { if sc == c { seen = true; break } }
      if !seen { chrom_order.push(c) }
      vri = vri + 1
    }
  }
  let sorted_records : Array[VcfRecord] = Array(capacity=0)
  let mut ci = 0
  while ci < chrom_order.length() {
    let chrom = chrom_order[ci]
    match grouped.get(chrom) {
      Some(list) => {
        // 按 pos 冒泡排
        let a = list.copy()
        let nn = a.length()
        let mut p = 0
        while p < nn {
          let mut q = p + 1
          while q < nn {
            if a[p].pos > a[q].pos {
              let tmp = a[p]
              a[p] = a[q]
              a[q] = tmp
            }
            q = q + 1
          }
          p = p + 1
        }
        let mut ai = 0
        while ai < a.length() {
          sorted_records.push(a[ai])
          ai = ai + 1
        }
      }
      None => ()
    }
    ci = ci + 1
  }
  VcfFile::{
    meta_headers: merged.meta_headers,
    fixed_header: merged.fixed_header,
    records: sorted_records,
  }
}

// ---------------------------------------------------------------------------
// 辅助变体
// ---------------------------------------------------------------------------

///|
/// Merge with explicit SAMPLE suffixes — useful when two inputs carry
/// identically-named samples that must be disambiguated in the output.
pub fn vcf_merge_with_sample_suffix(
  inputs : Array[VcfFile],
  suffixes : Array[String],
  sort~ : Bool = true,
) -> VcfFile {
  let merged_columns = merge_vcf_columns(inputs, sample_suffixes=suffixes)
  let raw = vcf_merge_unsorted(inputs)
  let base = VcfFile::{
    meta_headers: raw.meta_headers,
    fixed_header: Some(VcfFixedHeader::{ columns: merged_columns }),
    records: raw.records,
  }
  if sort {
    vcf_merge_sort_records_only(base)
  } else {
    base
  }
}

///|
/// Sort a VcfFile's records **only** (header untouched). Reused by the
/// suffix-variant above and exposed for callers who already built the
/// header themselves.
pub fn vcf_merge_sort_records_only(vcf : VcfFile) -> VcfFile {
  let records = vcf.records
  if records.length() <= 1 {
    return vcf
  }
  let grouped : Map[String, Array[VcfRecord]] = Map([], capacity=0)
  let mut ri = 0
  while ri < records.length() {
    let r = records[ri]
    match grouped.get(r.chrom) {
      Some(list) => {
        let copy = list.copy()
        copy.push(r)
        grouped.set(r.chrom, copy)
      }
      None => grouped.set(r.chrom, [r])
    }
    ri = ri + 1
  }
  let chrom_order : Array[String] = Array(capacity=0)
  ri = 0
  while ri < records.length() {
    let c = records[ri].chrom
    let mut seen = false
    for sc in chrom_order { if sc == c { seen = true; break } }
    if !seen { chrom_order.push(c) }
    ri = ri + 1
  }
  let sorted_records : Array[VcfRecord] = Array(capacity=0)
  let mut ci = 0
  while ci < chrom_order.length() {
    let chrom = chrom_order[ci]
    match grouped.get(chrom) {
      Some(list) => {
        let a = list.copy()
        let nn = a.length()
        let mut p = 0
        while p < nn {
          let mut q = p + 1
          while q < nn {
            if a[p].pos > a[q].pos {
              let tmp = a[p]
              a[p] = a[q]
              a[q] = tmp
            }
            q = q + 1
          }
          p = p + 1
        }
        let mut ai = 0
        while ai < a.length() {
          sorted_records.push(a[ai])
          ai = ai + 1
        }
      }
      None => ()
    }
    ci = ci + 1
  }
  VcfFile::{
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: sorted_records,
  }
}

// ---------------------------------------------------------------------------
// 统计辅助
// ---------------------------------------------------------------------------

///|
/// Return `(total_input_records, merged_records, unique_contigs_after)`.
pub fn vcf_merge_stats(merged : VcfFile, inputs : Array[VcfFile]) -> (Int32, Int32, Int32) {
  let mut total_input_records : Int32 = 0
  for vcf in inputs { total_input_records = total_input_records + vcf.records.length() }
  // unique contigs in merged
  let contigs : Set[String] = Set(Array(capacity=merged.records.length()))
  let mut ri = 0
  while ri < merged.records.length() {
    let _ = contigs.add(merged.records[ri].chrom)
    ri = ri + 1
  }
  (total_input_records, merged.records.length(), contigs.length())
}
