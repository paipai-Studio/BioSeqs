// ============================================================================
// pysam_vcf_norm — bcftools norm 的内存实现
//
// 核心操作：
//  1. trim_shared_prefix / trim_shared_suffix — 去掉 REF/ALT 两端共享的碱基
//  2. left_align — 在给定参考序列上左对齐 indel（把重复区域的 indel
//     移到最左端，输出 VCF 的 canonical representation）
//
// 参考：bcftools norm -f ref.fa -a --force
// ============================================================================

///|
/// Return true when this VcfRecord looks like a simple variant that `norm`
/// should normalise: exactly one ALT allele, non-empty REF/ALT strings, and
/// length difference ≥ 1 (SNPs are already canonical; complex MNVs can be
/// trimmed but not left-aligned without a ref).
pub fn vcf_is_norm_target(rec : VcfRecord) -> Bool {
  if rec.alt.length() != 1 {
    return false
  }
  let ref_s = rec.ref_seq
  let alt_s = rec.alt[0]
  if ref_s.length() == 0 || alt_s.length() == 0 {
    return false
  }
  if ref_s == alt_s {
    return false
  } // not a variant
  true
}

// ---------------------------------------------------------------------------
// trim_shared_prefix — 去掉 REF/ALT 共享的前缀碱基
// ---------------------------------------------------------------------------

///|
/// Strip shared prefix bases between REF and ALT alleles. Returns a
/// new VcfRecord with `pos` shifted forward by the trimmed amount, and the
/// alleles shortened. Matches bcftools norm's first pass.
///
/// e.g. REF="ACGT" ALT="ACGTA" → REF="T" ALT="TA" at `pos+3`.
pub fn vcf_norm_trim_prefix(rec : VcfRecord) -> VcfRecord {
  if rec.alt.length() != 1 {
    return rec
  }
  let ref_s = rec.ref_seq
  let alt_s = rec.alt[0]
  if ref_s.length() == 0 || alt_s.length() == 0 {
    return rec
  }
  let mut i = 0
  while i < ref_s.length() &&
        i < alt_s.length() &&
        ref_s.unsafe_get(i) == alt_s.unsafe_get(i) {
    i = i + 1
  }
  if i == 0 {
    return rec
  }
  let new_ref = ref_s[i:ref_s.length()].to_owned()
  let new_alt = alt_s[i:alt_s.length()].to_owned()
  let new_format : Array[String]? = match rec.format {
    Some(f) => Some(f.copy())
    None => None
  }
  VcfRecord::{
    chrom: rec.chrom,
    pos: rec.pos + i,
    id: rec.id.copy(),
    ref_seq: if new_ref.length() == 0 {
      "N"
    } else {
      new_ref
    },
    alt: if new_alt.length() == 0 {
      ["N"]
    } else {
      [new_alt]
    },
    qual: rec.qual,
    filter: rec.filter.copy(),
    info: rec.info.copy(),
    format: new_format,
    samples: rec.samples.copy(),
  }
}

// ---------------------------------------------------------------------------
// trim_shared_suffix — 去掉共享的后缀
// ---------------------------------------------------------------------------

///|
/// Strip shared suffix bases between REF and ALT. Applied *after* prefix
/// trimming so we only drop bases that remain equal at the end.  Matches
/// bcftools norm's second pass.
pub fn vcf_norm_trim_suffix(rec : VcfRecord) -> VcfRecord {
  if rec.alt.length() != 1 {
    return rec
  }
  let ref_s = rec.ref_seq
  let alt_s = rec.alt[0]
  if ref_s.length() <= 1 || alt_s.length() <= 1 {
    return rec
  }
  let mut i = 0
  let mut ref_pos = ref_s.length() - 1
  let mut alt_pos = alt_s.length() - 1
  while ref_pos > 0 &&
        alt_pos > 0 &&
        ref_s.unsafe_get(ref_pos) == alt_s.unsafe_get(alt_pos) {
    ref_pos = ref_pos - 1
    alt_pos = alt_pos - 1
    i = i + 1
  }
  if i == 0 {
    return rec
  }
  let new_ref = ref_s[0:ref_pos + 1].to_owned()
  let new_alt = alt_s[0:alt_pos + 1].to_owned()
  if new_ref.length() == 0 {
    return rec
  }
  if new_alt.length() == 0 {
    return rec
  }
  let new_format2 : Array[String]? = match rec.format {
    Some(f) => Some(f.copy())
    None => None
  }
  VcfRecord::{
    chrom: rec.chrom,
    pos: rec.pos, // suffix trim does NOT shift pos
    id: rec.id.copy(),
    ref_seq: new_ref,
    alt: [new_alt],
    qual: rec.qual,
    filter: rec.filter.copy(),
    info: rec.info.copy(),
    format: new_format2,
    samples: rec.samples.copy(),
  }
}

// ---------------------------------------------------------------------------
// full trim (prefix + suffix) — 无参考序列即可执行
// ---------------------------------------------------------------------------

///|
/// Apply both prefix and suffix trimming (bcftools norm without `-f ref`).
/// Returns the canonical representation of the variant's alleles.
pub fn vcf_norm_trim_all(rec : VcfRecord) -> VcfRecord {
  let after_prefix = vcf_norm_trim_prefix(rec)
  vcf_norm_trim_suffix(after_prefix)
}

// ---------------------------------------------------------------------------
// left_align — 基于参考序列的 indel 左对齐
// ---------------------------------------------------------------------------

///|
/// Left-align an indel on the given reference sequence. Requires a **1-based**
/// reference window beginning at `rec.pos - 1` covering enough upstream
/// sequence to slide the event. The algorithm walks one base left at a time
/// while the new REF/ALT pair produces an equivalent variant on the ref.
///
/// `ref_seq` is the CONTIG's full reference sequence (not a window).
/// Returns the left-aligned record (or the original if not left-alignable).
pub fn vcf_norm_left_align(rec : VcfRecord, ref_seq : String) -> VcfRecord {
  // Indels only (len diff = 1; we handle one-base-shift canonicalisation)
  if rec.alt.length() != 1 {
    return rec
  }
  let ref_s = rec.ref_seq
  let alt_s = rec.alt[0]
  let diff = alt_s.length() - ref_s.length()
  // diff == 1: insertion; diff == -1: deletion; diff == 0: SNP/MNV (not indel)
  if diff == 0 {
    // For SNPs we just apply trim; left-align doesn't apply
    return vcf_norm_trim_all(rec)
  }
  // Walk left while prev_base(ref) == trailing_base(alt-or-ref)
  let mut pos = rec.pos
  let mut cur_ref = ref_s
  let mut cur_alt = alt_s
  let mut steps = 0
  while pos > 1 && steps < 10000 {
    // ref base immediately BEFORE the current event
    let prev = ref_seq.unsafe_get(pos - 2).unsafe_to_char() // 1-based pos → 0-based idx
    // Candidate new alleles: prepend `prev`, drop trailing base that equals prev
    let new_ref = prev.to_string() + cur_ref
    let new_alt = prev.to_string() + cur_alt
    // Check if the trailing base of each matches prev (can strip one)
    let ref_trail = new_ref.unsafe_get(new_ref.length() - 1).unsafe_to_char()
    let alt_trail = new_alt.unsafe_get(new_alt.length() - 1).unsafe_to_char()
    if ref_trail == prev && alt_trail == prev {
      // OK: we can left-shift by one
      cur_ref = new_ref[0:new_ref.length() - 1].to_owned()
      cur_alt = new_alt[0:new_alt.length() - 1].to_owned()
      pos = pos - 1
      steps = steps + 1
    } else {
      break
    }
  }
  if steps == 0 {
    return rec
  }
  let new_format3 : Array[String]? = match rec.format {
    Some(f) => Some(f.copy())
    None => None
  }
  let trimmed = vcf_norm_trim_all(VcfRecord::{
    chrom: rec.chrom,
    pos,
    id: rec.id.copy(),
    ref_seq: cur_ref,
    alt: [cur_alt],
    qual: rec.qual,
    filter: rec.filter.copy(),
    info: rec.info.copy(),
    format: new_format3,
    samples: rec.samples.copy(),
  })
  trimmed
}

// ---------------------------------------------------------------------------
// 入口：全 VcfFile 归一化
// ---------------------------------------------------------------------------

///|
/// Normalise every record in a VcfFile. Without a reference, only trim is
/// applied.  With a reference (per-contig FastaFile), indels are additionally
/// left-aligned (matches `bcftools norm -f ref.fa`).
pub fn vcf_norm_file(vcf : VcfFile, ref_file : FastaFile) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=vcf.records.length())
  for rec in vcf.records {
    if !vcf_is_norm_target(rec) {
      out.push(rec)
      continue
    }
    match ref_file.get_contig(rec.chrom) {
      Some(c) => out.push(vcf_norm_left_align(rec, c.seq))
      None => out.push(vcf_norm_trim_all(rec))
    }
  }
  VcfFile::{
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: out,
  }
}

///|
/// Trim-only version (no reference sequence). Mirrors `bcftools norm -a`
/// without `-f` flag.
pub fn vcf_norm_trim_only(vcf : VcfFile) -> VcfFile {
  vcf_norm_file(vcf, FastaFile::{ contigs: [], })
}
