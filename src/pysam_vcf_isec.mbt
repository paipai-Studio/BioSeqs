// ============================================================================
// pysam_vcf_isec — bcftools isec 的内存版实现
//
// 功能：
//   多 VCF 间按位置精确匹配或重叠匹配，计算集合运算：
//     - 交集 (intersection): 所有输入都有的 records
//     - 只在第一个 (unique to first)
//     - 只在第二个 (unique to second)
//     - 并集 (union)
//     - 差集 (first minus second, second minus first)
//
// 参考：bcftools isec -n/-C/-w1/-w2 模式
// ============================================================================

///|
/// How records are matched between VCFs.
pub enum IsecMatchMode {
  /// `(chrom, pos)` 完全相等
  Exact
  /// `(chrom, pos, ref_seq)` 完全相等 — 适合精确 variant 匹配
  ExactWithRef
  /// 同一 contig 上 pos 距离 ≤ `tolerance`（默认 0）
  Nearby
}

pub fn isec_mode_exact() -> IsecMatchMode { IsecMatchMode::Exact }
pub fn isec_mode_exact_with_ref() -> IsecMatchMode { IsecMatchMode::ExactWithRef }
pub fn isec_mode_nearby() -> IsecMatchMode { IsecMatchMode::Nearby }

// ---------------------------------------------------------------------------
// 内部：构建 (chrom, pos) 键 Map 用于加速查找
// ---------------------------------------------------------------------------

///|
/// Build a Map from "chrom\tpos" → 该位置所有 record 的 Array 索引。
/// 支持同一位置多条 records。
fn build_pos_index(vcf : VcfFile) -> Map[String, Array[Int32]] {
  let m : Map[String, Array[Int32]] = Map([], capacity=vcf.records.length())
  let mut i = 0
  while i < vcf.records.length() {
    let r = vcf.records[i]
    let key = r.chrom + "\t" + r.pos.to_string()
    match m.get(key) {
      Some(arr) => {
        let copy = arr.copy()
        copy.push(i)
        m.set(key, copy)
      }
      None => {
        m.set(key, [i])
      }
    }
    i = i + 1
  }
  m
}

///|
/// 判断两条 record 是否在指定模式下匹配。
fn isec_match(
  a : VcfRecord,
  b : VcfRecord,
  mode : IsecMatchMode,
  tolerance~ : Int32 = 0,
) -> Bool {
  if a.chrom != b.chrom { return false }
  match mode {
    Exact => a.pos == b.pos
    ExactWithRef => a.pos == b.pos && a.ref_seq == b.ref_seq
    Nearby => {
      let d = a.pos - b.pos
      if d < 0 { -d } else { d }
    } <= tolerance
  }
}

///|
/// 从 `vcf` 中找出所有与 `target_rec` 匹配的 records。
fn find_matches(
  vcf : VcfFile,
  target_rec : VcfRecord,
  mode : IsecMatchMode,
  tolerance~ : Int32 = 0,
) -> Array[VcfRecord] {
  let out : Array[VcfRecord] = Array(capacity=0)
  for r in vcf.records {
    if isec_match(r, target_rec, mode, tolerance=tolerance) {
      out.push(r)
    }
  }
  out
}

// ---------------------------------------------------------------------------
// 主入口 1：两个 VCF 的集合运算
// ---------------------------------------------------------------------------

///|
/// 返回只在 `a` 中的 records（不在 `b` 里匹配到任何一条）。
/// 等价 `bcftools isec -C a.vcf b.vcf` 中第一个输出文件。
pub fn vcf_isec_unique_first(
  a : VcfFile,
  b : VcfFile,
  mode~ : IsecMatchMode = Exact,
  tolerance~ : Int32 = 0,
) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=0)
  for ra in a.records {
    let mut found = false
    for rb in b.records {
      if isec_match(ra, rb, mode, tolerance=tolerance) {
        found = true
        break
      }
    }
    if !found { out.push(ra) }
  }
  VcfFile::{
    meta_headers: a.meta_headers,
    fixed_header: a.fixed_header,
    records: out,
  }
}

///|
/// 返回只在 `b` 中的 records。
pub fn vcf_isec_unique_second(
  a : VcfFile,
  b : VcfFile,
  mode~ : IsecMatchMode = Exact,
  tolerance~ : Int32 = 0,
) -> VcfFile {
  vcf_isec_unique_first(b, a, mode=mode, tolerance=tolerance)
}

///|
/// 返回 `a` 与 `b` 的交集 —— `a` 中每条在 `b` 里能找到匹配的 record。
/// 注意：输出 records 来自 `a`（保留 a 的样本信息）。
/// 等价 `bcftools isec -n=2 a.vcf b.vcf`。
pub fn vcf_isec_intersection(
  a : VcfFile,
  b : VcfFile,
  mode~ : IsecMatchMode = Exact,
  tolerance~ : Int32 = 0,
) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=0)
  for ra in a.records {
    let mut found = false
    for rb in b.records {
      if isec_match(ra, rb, mode, tolerance=tolerance) {
        found = true
        break
      }
    }
    if found { out.push(ra) }
  }
  VcfFile::{
    meta_headers: a.meta_headers,
    fixed_header: a.fixed_header,
    records: out,
  }
}

///|
/// 返回 `a` 与 `b` 的并集（`a` 的所有 records + `b` 里没在 `a` 里匹配到的 records）。
/// 等价 `bcftools isec -n~=0 a.vcf b.vcf` 的合并输出。
pub fn vcf_isec_union(
  a : VcfFile,
  b : VcfFile,
  mode~ : IsecMatchMode = Exact,
  tolerance~ : Int32 = 0,
) -> VcfFile {
  // a 的全部 records
  let out_a = a.records.copy()
  // 加上 b 里不在 a 中的 records
  let only_b = vcf_isec_unique_second(a, b, mode=mode, tolerance=tolerance)
  for rb in only_b.records {
    out_a.push(rb)
  }
  // 合并 meta headers（去重）
  let meta = merge_meta_headers(a.meta_headers, b.meta_headers)
  let fixed = match a.fixed_header {
    Some(fh) => Some(fh)
    None => b.fixed_header
  }
  VcfFile::{
    meta_headers: meta,
    fixed_header: fixed,
    records: out_a,
  }
}

// ---------------------------------------------------------------------------
// 主入口 2：多 VCF 交集（N 路）
// ---------------------------------------------------------------------------

///|
/// 计算所有输入 VCFs 的 N 路交集：每条 record 必须在 **每一个** VCF 中
/// 都能找到匹配才算命中。返回的 records 来自 `inputs[0]`。
///
/// 等价 `bcftools isec -n=<N>` 多输入情况。
pub fn vcf_isec_intersect_all(
  inputs : Array[VcfFile],
  mode~ : IsecMatchMode = Exact,
  tolerance~ : Int32 = 0,
) -> VcfFile {
  if inputs.length() == 0 {
    return VcfFile::{ meta_headers: [], fixed_header: None, records: [] }
  }
  if inputs.length() == 1 {
    return inputs[0]
  }
  let first = inputs[0]
  let rest_sorted : Array[VcfFile] = Array(capacity=inputs.length() - 1)
  {
    let mut fi = 1
    while fi < inputs.length() {
      rest_sorted.push(inputs[fi])
      fi = fi + 1
    }
  }
  let out : Array[VcfRecord] = Array(capacity=0)
  for ra in first.records {
    let mut all_found = true
    let mut fi = 0
    while fi < rest_sorted.length() {
      let b = rest_sorted[fi]
      let mut found = false
      for rb in b.records {
        if isec_match(ra, rb, mode, tolerance=tolerance) {
          found = true
          break
        }
      }
      if !found { all_found = false; break }
      fi = fi + 1
    }
    if all_found { out.push(ra) }
  }
  VcfFile::{
    meta_headers: first.meta_headers,
    fixed_header: first.fixed_header,
    records: out,
  }
}

///|
/// 计算所有输入的差集 —— 只在 `inputs[0]` 中出现、在其余所有 VCFs 中
/// 都找不到匹配的 records。
pub fn vcf_isec_unique_to_first(
  inputs : Array[VcfFile],
  mode~ : IsecMatchMode = Exact,
  tolerance~ : Int32 = 0,
) -> VcfFile {
  if inputs.length() == 0 {
    return VcfFile::{ meta_headers: [], fixed_header: None, records: [] }
  }
  if inputs.length() == 1 {
    return inputs[0]
  }
  let first = inputs[0]
  let rest : Array[VcfFile] = Array(capacity=inputs.length() - 1)
  {
    let mut fi = 1
    while fi < inputs.length() {
      rest.push(inputs[fi])
      fi = fi + 1
    }
  }
  let out : Array[VcfRecord] = Array(capacity=0)
  for ra in first.records {
    let mut all_missing = true
    let mut fi = 0
    while fi < rest.length() {
      let b = rest[fi]
      let mut hit_b = false
      for rb in b.records {
        if isec_match(ra, rb, mode, tolerance=tolerance) {
          hit_b = true
          break
        }
      }
      if hit_b { all_missing = false; break }
      fi = fi + 1
    }
    if all_missing { out.push(ra) }
  }
  VcfFile::{
    meta_headers: first.meta_headers,
    fixed_header: first.fixed_header,
    records: out,
  }
}

// ---------------------------------------------------------------------------
// 辅助：meta header 去重合并
// ---------------------------------------------------------------------------

fn merge_meta_headers(
  a : Array[VcfMetaHeader],
  b : Array[VcfMetaHeader],
) -> Array[VcfMetaHeader] {
  let out : Array[VcfMetaHeader] = a.copy()
  for mh in b {
    // 简单去重：相同 key+value 跳过
    let mut dup = false
    for existing in out {
      if existing.key == mh.key && existing.value == mh.value {
        dup = true
        break
      }
    }
    if !dup { out.push(mh) }
  }
  out
}

// ---------------------------------------------------------------------------
// 辅助：用 Map 加速的二元交集实现（输入较大时推荐使用）
// ---------------------------------------------------------------------------

///|
/// 二元交集的 Map 加速版：先对 `b` 按 pos 建索引，再遍历 `a` 的记录
/// 做精确查找。在 Exact / ExactWithRef 模式下效果显著。
pub fn vcf_isec_intersection_fast(
  a : VcfFile,
  b : VcfFile,
  mode~ : IsecMatchMode = Exact,
) -> VcfFile {
  let idx_b = build_pos_index(b)
  let out : Array[VcfRecord] = Array(capacity=0)
  for ra in a.records {
    let key = ra.chrom + "\t" + ra.pos.to_string()
    match idx_b.get(key) {
      Some(b_indices) => {
        let mut matched = false
        for bi in b_indices {
          let rb = b.records[bi]
          match mode {
            Exact => { matched = true; break }
            ExactWithRef => {
              if ra.ref_seq == rb.ref_seq { matched = true; break }
            }
            Nearby => { matched = true; break }
          }
        }
        if matched { out.push(ra) }
      }
      None => ()
    }
  }
  VcfFile::{
    meta_headers: a.meta_headers,
    fixed_header: a.fixed_header,
    records: out,
  }
}

// ---------------------------------------------------------------------------
// 统计辅助
// ---------------------------------------------------------------------------

///|
/// 返回 `(only_in_a, only_in_b, shared)` 三元组计数。
pub fn vcf_isec_summary(
  a : VcfFile,
  b : VcfFile,
  mode~ : IsecMatchMode = Exact,
  tolerance~ : Int32 = 0,
) -> (Int32, Int32, Int32) {
  let only_a = vcf_isec_unique_first(a, b, mode=mode, tolerance=tolerance).records.length()
  let shared = a.records.length() - only_a
  let only_b = vcf_isec_unique_second(a, b, mode=mode, tolerance=tolerance).records.length()
  (only_a, only_b, shared)
}

///|
/// Jaccard index（两 VCF 的相似度：|A∩B| / |A∪B|）。范围 0.0~1.0。
pub fn vcf_isec_jaccard(
  a : VcfFile,
  b : VcfFile,
  mode~ : IsecMatchMode = Exact,
) -> Double {
  let inter = vcf_isec_intersection_fast(a, b, mode=mode).records.length().to_double()
  let union_a = a.records.length().to_double()
  let union_b = b.records.length().to_double()
  let union = union_a + union_b - inter
  if union == 0.0 { return 0.0 }
  inter / union
}
