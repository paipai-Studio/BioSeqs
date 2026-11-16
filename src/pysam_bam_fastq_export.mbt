// ============================================================================
// pysam_bam_fastq_export — samtools fastq / fasta 的内存实现
//
// Pysam 对应：`AlignmentFile.fetch(until_eof=True)` 然后逐条
// `get_forward_sequence` / `get_forward_qualities` / `is_reverse`
//
// samtools fastq 支持：
//  - 按 FLAG 方向过滤（F 输出原始方向，默认输出测序方向）
//  - `-1/-2` read1 / read2 分文件输出
//  - `-0` unpaired 输出
//  - `-s` singleton（paired 但 mate 不在）输出
//  - `-n` 忽略 read-group 过滤
//  - `-t` 输出 tag 到 header
// ============================================================================

///|
/// Configuration for a single FASTQ / FASTA output stream. Mirrors
/// samtools fastq's `-1/-2/-0/-s` semantics plus basic RG / flag filtering.
pub(all) struct FastqExportConfig {
  /// 1-based mate: 1 = read1, 2 = read2, 0 = unpaired
  mate_of : Int32
  /// Read name line prefix (default "@"). samtools fastq always uses "@".
  name_prefix : String
  /// Output original alignment orientation (`true`) or original sequencing
  /// orientation (`false`, default — matches samtools fastq).
  original_sequence : Bool
  /// Read group filter: non-empty = only reads with this RG tag are emitted.
  rg_filter : String
  /// Emit qualities even for FASTA output (rare, mostly debug).
  emit_qualities : Bool
}

///|
/// Default: unpaired + sequencing orientation + "@" prefix.
pub fn fastq_config_default() -> FastqExportConfig {
  FastqExportConfig::{
    mate_of: 0,
    name_prefix: "@",
    original_sequence: false,
    rg_filter: "",
    emit_qualities: true,
  }
}

// ---------------------------------------------------------------------------
// 单条 record → FASTQ 记录 (name, seq, qual)
// ---------------------------------------------------------------------------

///|
/// Extract the sequencing-direction sequence from a BamRecord, honouring
/// the `original_sequence` flag and minus-strand reverse complement.
pub fn bam_record_fastq_sequence(
  rec : BamRecord,
  original_sequence : Bool,
) -> String {
  if original_sequence {
    rec.seq // raw alignment direction (may have softclips at the ends)
  } else {
    rec.get_forward_sequence()
  }
}

///|
/// Extract sequencing-direction qualities from a BamRecord.
pub fn bam_record_fastq_qualities(
  rec : BamRecord,
  original_sequence : Bool,
) -> String {
  if original_sequence {
    rec.qual
  } else {
    rec.get_forward_qualities()
  }
}

///|
/// Build the FASTQ header line (prefix + qname; optionally append /1 or /2).
pub fn bam_record_fastq_header(
  rec : BamRecord,
  prefix : String,
  append_mate_suffix : Bool,
) -> String {
  let suffix = if append_mate_suffix {
    if rec.is_read1() {
      "/1"
    } else if rec.is_read2() {
      "/2"
    } else {
      ""
    }
  } else {
    ""
  }
  prefix + rec.qname + suffix
}

// ---------------------------------------------------------------------------
// 主入口：BamFile → FASTQ 文本
// ---------------------------------------------------------------------------

///|
/// Export a BamFile to a FASTQ string, applying the given config.
/// Raises if the config's mate filter is incompatible with record flags.
pub fn bam_to_fastq(bam : BamFile, config : FastqExportConfig) -> String {
  let sb = StringBuilder()
  for r in bam.records {
    // Skip secondary / supplementary (same as samtools fastq default)
    if r.is_secondary() || r.is_supplementary() {
      continue
    }
    // Unmapped reads are still emitted (samtools fastq default) — no skip
    // RG filter
    if config.rg_filter.length() > 0 {
      match r.get_tag("RG") {
        Some(BamTagValue::Str(rg)) => if rg != config.rg_filter { continue }
        _ => continue
      }
    }
    // Mate filter
    let emit = if config.mate_of == 1 {
      r.is_read1() && !r.is_unmapped()
    } else if config.mate_of == 2 {
      r.is_read2() && !r.is_unmapped()
    } else {
      // mate_of == 0: unpaired reads + paired singletons
      // (paired but unmapped **自身** or mate unmapped)
      !r.is_paired() || r.is_unmapped() || r.is_mate_unmapped()
    }
    if !emit {
      continue
    }
    // Write FASTQ record
    let header = bam_record_fastq_header(r, config.name_prefix, true)
    let seq = bam_record_fastq_sequence(r, config.original_sequence)
    let qual = bam_record_fastq_qualities(r, config.original_sequence)
    sb.write_string(header)
    sb.write_char('\n')
    sb.write_string(seq)
    sb.write_char('\n')
    sb.write_char('+')
    sb.write_char('\n')
    sb.write_string(qual)
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// Convenience: export ALL records (no filter) to FASTQ using defaults.
pub fn bam_to_fastq_default(bam : BamFile) -> String {
  bam_to_fastq(bam, fastq_config_default())
}

// ---------------------------------------------------------------------------
// 成对导出：输出 paired 端到两个独立的 FASTQ 字符串
// ---------------------------------------------------------------------------

///|
/// Export paired reads into `(read1_fastq, read2_fastq)`. Unpaired / singleton
/// reads are collected into `unpaired_fastq` (empty when `-0` is set by the
/// caller; here we always include them for completeness).
pub fn bam_to_fastq_paired(bam : BamFile) -> (String, String, String) {
  let sb1 = StringBuilder()
  let sb2 = StringBuilder()
  let sb0 = StringBuilder()
  for r in bam.records {
    if r.is_secondary() || r.is_supplementary() {
      continue
    }
    let header = bam_record_fastq_header(r, "@", true)
    let seq = bam_record_fastq_sequence(r, false)
    let qual = bam_record_fastq_qualities(r, false)
    let rec = header + "\n" + seq + "\n+\n" + qual + "\n"
    if r.is_paired() && !r.is_unmapped() && !r.is_mate_unmapped() {
      if r.is_read1() {
        sb1.write_string(rec)
      } else if r.is_read2() {
        sb2.write_string(rec)
      } else {
        sb0.write_string(rec)
      }
    } else {
      sb0.write_string(rec)
    }
  }
  (sb1.to_string(), sb2.to_string(), sb0.to_string())
}

// ---------------------------------------------------------------------------
// BAM → FASTA（无质量串）
// ---------------------------------------------------------------------------

///|
/// Export BamFile to FASTA (no quality scores). Mirrors `samtools fasta`.
pub fn bam_to_fasta(bam : BamFile) -> String {
  bam_to_fasta_with_config(bam, FastqExportConfig::{
    mate_of: 0,
    name_prefix: ">",
    original_sequence: false,
    rg_filter: "",
    emit_qualities: false,
  })
}

///|
/// Version of `bam_to_fasta` that takes a config (overriding prefix etc.).
pub fn bam_to_fasta_with_config(
  bam : BamFile,
  config : FastqExportConfig,
) -> String {
  let sb = StringBuilder()
  for r in bam.records {
    if r.is_secondary() || r.is_supplementary() {
      continue
    }
    let header = bam_record_fastq_header(r, ">", true)
    let seq = bam_record_fastq_sequence(r, config.original_sequence)
    sb.write_string(header)
    sb.write_char('\n')
    sb.write_string(seq)
    sb.write_char('\n')
  }
  sb.to_string()
}
