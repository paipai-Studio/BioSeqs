///|
/// A single genomic interval (0-based, half-open [start, end]).
/// Ported from htslib synced_bcf_reader.c region1_t.
pub struct Region1 {
  start : Int
  end : Int
} derive(Eq)
pub extend Region1 with Eq::{equal, not_equal}

///|
/// Regions for a single chromosome/sequence.
/// Ported from htslib synced_bcf_reader.c region_t.
pub struct SeqRegions {
  seq_name : String
  regs : Array[Region1]
  creg : Int // current active region index, -1 = none
}

///|
/// A collection of genomic regions, parsed from strings like
/// "chr1:1000-2000,chr2:5000-6000".
/// Ported from htslib synced_bcf_reader.c bcf_sr_regions_t (in-memory part).
pub struct BcfSrRegions {
  seqs : Array[SeqRegions]
  iseq : Int // current sequence index
  start : Int // current region start (0-based)
  end : Int // current region end (0-based)
  prev_seq : Int
  prev_start : Int
  prev_end : Int
}

///|
/// Maximum coordinate value (analogous to MAX_CSI_COOR in htslib).
let max_csi_coor : Int = 1073741823 // (1 << 30) - 1

///|
/// Find or create a sequence index in the regions collection.
fn find_or_add_seq(
  seqs : Array[SeqRegions],
  name : String
) -> Int {
  for i = 0; i < seqs.length(); i = i + 1 {
    if seqs[i].seq_name == name {
      return i
    }
  }
  -1
}

///|
/// Add a region (chr, start, end) to the collection.
/// Coordinates are 1-based on input, stored as 0-based internally.
/// Ported from htslib _regions_add().
pub fn BcfSrRegions::add_region(
  self : BcfSrRegions,
  chr : String,
  start : Int,
  end : Int
) -> BcfSrRegions {
  let mut s = start
  let mut e = end
  if s == -1 && e == -1 {
    s = 0
    e = max_csi_coor
  } else {
    s = s - 1
    e = e - 1
  }
  let new_seqs = []
  let mut found_idx = -1
  for i = 0; i < self.seqs.length(); i = i + 1 {
    new_seqs.push(self.seqs[i])
    if self.seqs[i].seq_name == chr {
      found_idx = i
    }
  }
  if found_idx < 0 {
    found_idx = new_seqs.length()
    new_seqs.push(SeqRegions::{ seq_name: chr, regs: [], creg: -1 })
  }
  let seq = new_seqs[found_idx]
  let new_regs = []
  for i = 0; i < seq.regs.length(); i = i + 1 {
    new_regs.push(seq.regs[i])
  }
  new_regs.push(Region1::{ start: s, end: e })
  new_seqs[found_idx] = SeqRegions::{
    seq_name: seq.seq_name,
    regs: new_regs,
    creg: seq.creg,
  }
  BcfSrRegions::{
    seqs: new_seqs,
    iseq: self.iseq,
    start: self.start,
    end: self.end,
    prev_seq: self.prev_seq,
    prev_start: self.prev_start,
    prev_end: self.prev_end,
  }
}

///|
/// Parse a region string like "chr1:1000-2000,chr2,chr3:500".
/// Ported from htslib _regions_init_string().
pub fn bcf_sr_regions_init_string(regions_str : String) -> BcfSrRegions {
  let mut result = BcfSrRegions::{
    seqs: [],
    iseq: 0,
    start: -1,
    end: -1,
    prev_seq: -1,
    prev_start: -1,
    prev_end: -1,
  }
  let n = regions_str.length()
  let mut pos = 0
  while pos < n {
    // Skip leading comma
    if pos < n && regions_str[pos] == ',' {
      pos = pos + 1
    }
    if pos >= n {
      break
    }
    // Parse chromosome name (handle {weird-name} syntax)
    let mut chr = ""
    if regions_str[pos] == '{' {
      let mut end_brace = pos + 1
      while end_brace < n && regions_str[end_brace] != '}' {
        end_brace = end_brace + 1
      }
      if end_brace < n {
        chr = regions_str[pos + 1:end_brace].to_owned()
        pos = end_brace + 1
      } else {
        break
      }
    } else {
      let mut end_chr = pos
      while end_chr < n &&
        regions_str[end_chr] != ',' &&
        regions_str[end_chr] != ':' {
        end_chr = end_chr + 1
      }
      chr = regions_str[pos:end_chr].to_owned()
      pos = end_chr
    }
    if chr.length() == 0 {
      if pos < n && regions_str[pos] == ',' {
        pos = pos + 1
      }
      continue
    }
    // Parse position range
    if pos < n && regions_str[pos] == ':' {
      pos = pos + 1
      // Parse "from"
      let from_start = pos
      while pos < n &&
        regions_str[pos] != '-' &&
        regions_str[pos] != ',' {
        pos = pos + 1
      }
      let from = parse_int(regions_str[from_start:pos].to_owned())
      if pos < n && regions_str[pos] == '-' {
        pos = pos + 1
        // Parse "to"
        let to_start = pos
        while pos < n && regions_str[pos] != ',' {
          pos = pos + 1
        }
        if pos == to_start {
          result = result.add_region(chr, from, max_csi_coor)
        } else {
          result = result.add_region(chr, from, parse_int(regions_str[to_start:pos].to_owned()))
        }
      } else {
        // Single position: from == to
        result = result.add_region(chr, from, from)
      }
    } else {
      // Whole chromosome
      result = result.add_region(chr, -1, -1)
    }
  }
  result.sort_and_merge()
}

///|
/// Sort regions by start position and merge overlapping intervals.
/// Ported from htslib _regions_sort_and_merge() + regions_merge().
pub fn BcfSrRegions::sort_and_merge(self : BcfSrRegions) -> BcfSrRegions {
  let new_seqs = []
  for i = 0; i < self.seqs.length(); i = i + 1 {
    let seq = self.seqs[i]
    let regs = []
    for j = 0; j < seq.regs.length(); j = j + 1 {
      regs.push(seq.regs[j])
    }
    // Sort by (start, end)
    regs.sort_by(fn(a, b) {
      if a.start != b.start {
        a.start - b.start
      } else {
        a.end - b.end
      }
    })
    // Merge overlapping intervals
    let merged = []
    let mut i2 = 0
    while i2 < regs.length() {
      let mut current = regs[i2]
      let mut j = i2 + 1
      while j < regs.length() && current.end >= regs[j].start {
        if current.end < regs[j].end {
          current = Region1::{ start: current.start, end: regs[j].end }
        }
        j = j + 1
      }
      merged.push(current)
      i2 = j
    }
    new_seqs.push(SeqRegions::{
      seq_name: seq.seq_name,
      regs: merged,
      creg: seq.creg,
    })
  }
  BcfSrRegions::{
    seqs: new_seqs,
    iseq: self.iseq,
    start: self.start,
    end: self.end,
    prev_seq: self.prev_seq,
    prev_start: self.prev_start,
    prev_end: self.prev_end,
  }
}

///|
/// Check if a position/interval overlaps any region.
/// Returns 0 if overlap, -1 if no overlap, -2 if no more regions on this seq.
/// Ported from htslib _bcf_sr_regions_overlap().
pub fn BcfSrRegions::overlap(
  self : BcfSrRegions,
  seq : String,
  start : Int,
  end : Int
) -> Int {
  // Find sequence index
  let mut iseq = -1
  for i = 0; i < self.seqs.length(); i = i + 1 {
    if self.seqs[i].seq_name == seq {
      iseq = i
      break
    }
  }
  if iseq < 0 {
    return -1
  }
  // Check each region in this sequence
  let regs = self.seqs[iseq].regs
  for i = 0; i < regs.length(); i = i + 1 {
    // Region overlaps [start, end] if regs[i].start <= end && regs[i].end >= start
    if regs[i].start <= end && regs[i].end >= start {
      return 0
    }
    // Regions are sorted, if past the query, no overlap
    if regs[i].start > end {
      return -1
    }
  }
  -1
}

///|
/// Get the number of sequences.
pub fn BcfSrRegions::nseqs(self : BcfSrRegions) -> Int {
  self.seqs.length()
}

///|
/// Get sequence name by index.
pub fn BcfSrRegions::seq_name(self : BcfSrRegions, idx : Int) -> String? {
  if idx < 0 || idx >= self.seqs.length() {
    return None
  }
  Some(self.seqs[idx].seq_name)
}

///|
/// Get the number of regions for a sequence.
pub fn BcfSrRegions::nregs_for_seq(self : BcfSrRegions, seq : String) -> Int {
  for i = 0; i < self.seqs.length(); i = i + 1 {
    if self.seqs[i].seq_name == seq {
      return self.seqs[i].regs.length()
    }
  }
  0
}

///|
/// Get all regions for a sequence as (start, end) pairs (0-based).
pub fn BcfSrRegions::get_regs(
  self : BcfSrRegions,
  seq : String
) -> Array[(Int, Int)] {
  for i = 0; i < self.seqs.length(); i = i + 1 {
    if self.seqs[i].seq_name == seq {
      let result = []
      for j = 0; j < self.seqs[i].regs.length(); j = j + 1 {
        result.push((self.seqs[i].regs[j].start, self.seqs[i].regs[j].end))
      }
      return result
    }
  }
  []
}

// ==========================================================================
// Variant type classification and pair/collapse logic
// Ported from htslib bcf_sr_sort.c
// ==========================================================================

///|
/// Variant type flags (ported from bcf_sr_sort.c SR_* defines).
pub let sr_ref : Int = 1
pub let sr_snp : Int = 2
pub let sr_indel : Int = 4
pub let sr_other : Int = 8

///|
/// Pair logic flags (ported from htslib BCF_SR_PAIR_* defines).
pub let bcf_sr_pair_snps : Int = 1
pub let bcf_sr_pair_indels : Int = 2
pub let bcf_sr_pair_any : Int = 4
pub let bcf_sr_pair_some : Int = 8
pub let bcf_sr_pair_snp_ref : Int = 16
pub let bcf_sr_pair_indel_ref : Int = 32
pub let bcf_sr_pair_exact : Int = 64
pub let bcf_sr_pair_id : Int = 128
pub let bcf_sr_pair_both : Int = 3 // SNPS | INDELS
pub let bcf_sr_pair_both_ref : Int = 51 // SNPS|INDELS|SNP_REF|INDEL_REF

///|
/// Collapse flags (legacy, mapped to pair logic).
pub let collapse_none : Int = 0
pub let collapse_snps : Int = 1
pub let collapse_indels : Int = 2
pub let collapse_any : Int = 4
pub let collapse_some : Int = 8
pub let collapse_both : Int = 3

///|
/// Classify an allele as SNP, indel, ref, or other.
/// Returns one of sr_ref, sr_snp, sr_indel, sr_other.
/// Ported from bcf_sr_sort.c variant type classification.
pub fn classify_allele(ref_allele : String, alt_allele : String) -> Int {
  if alt_allele == "." || alt_allele == "" {
    return sr_ref
  }
  let ref_len = ref_allele.length()
  let alt_len = alt_allele.length()
  if ref_len == 1 && alt_len == 1 {
    return sr_snp
  }
  if ref_len == alt_len {
    // MNP, treat as SNP-like
    return sr_snp
  }
  // One is longer than the other → indel
  if ref_len == 0 || alt_len == 0 {
    return sr_indel
  }
  // Check if it's a pure insertion/deletion
  let shorter = if ref_len < alt_len { ref_allele } else { alt_allele }
  let longer = if ref_len < alt_len { alt_allele } else { ref_allele }
  // If shorter is a prefix of longer, it's a clean indel
  let mut is_prefix = true
  for i = 0; i < shorter.length(); i = i + 1 {
    if shorter[i] != longer[i] {
      is_prefix = false
      break
    }
  }
  if is_prefix {
    return sr_indel
  }
  sr_other
}

///|
/// Classify a variant given REF and all ALTs.
/// Returns the combined type (bitwise OR of all allele types).
pub fn classify_variant(ref_allele : String, alts : Array[String]) -> Int {
  let mut vtype = 0
  for i = 0; i < alts.length(); i = i + 1 {
    vtype = vtype | classify_allele(ref_allele, alts[i])
  }
  if vtype == 0 {
    vtype = sr_ref
  }
  vtype
}

///|
/// Check if two variant types are compatible given pair logic.
/// Ported from bcf_sr_sort.c bcf_sr_init_scores() scoring.
pub fn pair_compatible(type_a : Int, type_b : Int, pair_logic : Int) -> Bool {
  let mut pair = pair_logic
  if (pair & bcf_sr_pair_any) != 0 {
    pair = pair | bcf_sr_pair_snps | bcf_sr_pair_indels | bcf_sr_pair_snp_ref | bcf_sr_pair_indel_ref
  }
  if (pair & bcf_sr_pair_exact) != 0 {
    // Exact match required — only same type
    return type_a == type_b
  }
  // SNP-SNP
  if (type_a & sr_snp) != 0 && (type_b & sr_snp) != 0 && (pair & bcf_sr_pair_snps) != 0 {
    return true
  }
  // INDEL-INDEL
  if (type_a & sr_indel) != 0 && (type_b & sr_indel) != 0 && (pair & bcf_sr_pair_indels) != 0 {
    return true
  }
  // SNP-REF
  if (type_a & sr_snp) != 0 && (type_b & sr_ref) != 0 && (pair & bcf_sr_pair_snp_ref) != 0 {
    return true
  }
  if (type_a & sr_ref) != 0 && (type_b & sr_snp) != 0 && (pair & bcf_sr_pair_snp_ref) != 0 {
    return true
  }
  // INDEL-REF
  if (type_a & sr_indel) != 0 && (type_b & sr_ref) != 0 && (pair & bcf_sr_pair_indel_ref) != 0 {
    return true
  }
  if (type_a & sr_ref) != 0 && (type_b & sr_indel) != 0 && (pair & bcf_sr_pair_indel_ref) != 0 {
    return true
  }
  // ANY allows any combination
  if (pair_logic & bcf_sr_pair_any) != 0 {
    return true
  }
  false
}

// ==========================================================================
// Synced reader for in-memory VCF records
// Ported from htslib synced_bcf_reader.c next_line() algorithm
// ==========================================================================

///|
/// A minimal VCF variant record for synced reading.
pub struct SyncVariant {
  chrom : String
  pos : Int // 0-based position
  ref_allele : String
  alts : Array[String]
  id : String
} derive(Eq)
pub extend SyncVariant with Eq::{equal, not_equal}

///|
/// A single reader's buffer of variants.
pub struct SyncReader {
  name : String
  buffer : Array[SyncVariant]
  pos : Int // current position in buffer
}

///|
/// Synced reader state.
/// Ported from htslib bcf_srs_t (in-memory version).
pub struct BcfSrs {
  readers : Array[SyncReader]
  has_line : Array[Bool]
  pair_logic : Int
  errnum : Int
}

///|
/// Error codes (ported from bcf_sr_error enum).
pub let bcf_sr_ok : Int = 0
pub let bcf_sr_open_failed : Int = 1
pub let not_bgzf : Int = 2
pub let idx_load_failed : Int = 3
pub let file_type_error : Int = 4
pub let api_usage_error : Int = 5
pub let header_error : Int = 6
pub let vcf_parse_error : Int = 7
pub let bcf_read_error : Int = 8
pub let noidx_error : Int = 9

///|
/// Initialize a synced reader.
/// Ported from htslib bcf_sr_init().
pub fn bcf_sr_init() -> BcfSrs {
  BcfSrs::{
    readers: [],
    has_line: [],
    pair_logic: bcf_sr_pair_both_ref,
    errnum: bcf_sr_ok,
  }
}

///|
/// Set pair logic (collapse mode).
/// Ported from htslib bcf_sr_set_opt(BCF_SR_PAIR_LOGIC).
pub fn BcfSrs::set_pair_logic(self : BcfSrs, logic : Int) -> BcfSrs {
  BcfSrs::{ readers: self.readers, has_line: self.has_line, pair_logic: logic, errnum: self.errnum }
}

///|
/// Add a reader with a name and a list of variants.
/// Ported from htslib bcf_sr_add_reader() (in-memory version).
pub fn BcfSrs::add_reader(
  self : BcfSrs,
  name : String,
  variants : Array[SyncVariant]
) -> BcfSrs {
  let new_readers = []
  for i = 0; i < self.readers.length(); i = i + 1 {
    new_readers.push(self.readers[i])
  }
  new_readers.push(SyncReader::{ name, buffer: variants, pos: 0 })
  let new_has_line = []
  for i = 0; i < self.has_line.length(); i = i + 1 {
    new_has_line.push(self.has_line[i])
  }
  new_has_line.push(false)
  BcfSrs::{
    readers: new_readers,
    has_line: new_has_line,
    pair_logic: self.pair_logic,
    errnum: self.errnum,
  }
}

///|
/// Number of readers.
pub fn BcfSrs::nreaders(self : BcfSrs) -> Int {
  self.readers.length()
}

///|
/// Get the current line for a given reader, or None.
/// Ported from htslib bcf_sr_get_line().
pub fn BcfSrs::get_line(self : BcfSrs, idx : Int) -> SyncVariant? {
  if idx < 0 || idx >= self.has_line.length() {
    return None
  }
  if !self.has_line[idx] {
    return None
  }
  let reader = self.readers[idx]
  if reader.pos <= 0 || reader.pos > reader.buffer.length() {
    return None
  }
  Some(reader.buffer[reader.pos - 1])
}

///|
/// Check if a reader has a line at the current position.
/// Ported from htslib bcf_sr_has_line macro.
pub fn BcfSrs::has_line_at(self : BcfSrs, idx : Int) -> Bool {
  if idx < 0 || idx >= self.has_line.length() {
    return false
  }
  self.has_line[idx]
}

///|
/// Compare two variants for compatibility.
/// Returns true if they can be paired according to pair_logic.
fn variants_compatible(
  a : SyncVariant,
  b : SyncVariant,
  pair_logic : Int
) -> Bool {
  // Must be at same position
  if a.pos != b.pos {
    return false
  }
  // Must be on same chromosome
  if a.chrom != b.chrom {
    return false
  }
  if (pair_logic & bcf_sr_pair_exact) != 0 {
    // Exact: same REF and same ALT set
    if a.ref_allele != b.ref_allele {
      return false
    }
    if a.alts.length() != b.alts.length() {
      return false
    }
    for i = 0; i < a.alts.length(); i = i + 1 {
      if a.alts[i] != b.alts[i] {
        return false
      }
    }
    return true
  }
  if (pair_logic & bcf_sr_pair_id) != 0 {
    if a.id != b.id {
      return false
    }
  }
  let type_a = classify_variant(a.ref_allele, a.alts)
  let type_b = classify_variant(b.ref_allele, b.alts)
  pair_compatible(type_a, type_b, pair_logic)
}

///|
/// Advance to the next synchronized line.
/// Returns the number of readers with a line at this position, or 0 if done.
/// Ported from htslib bcf_sr_next_line() / next_line().
pub fn BcfSrs::next_line(self : BcfSrs) -> (BcfSrs, Int) {
  let n = self.readers.length()
  if n == 0 {
    return (self, 0)
  }
  // Find minimum (chrom, pos) across all readers
  let mut min_chrom = ""
  let mut min_pos = 2147483647 // INT_MAX
  let mut found = false
  for i = 0; i < n; i = i + 1 {
    let reader = self.readers[i]
    if reader.pos >= reader.buffer.length() {
      continue
    }
    let v = reader.buffer[reader.pos]
    if !found ||
      v.chrom < min_chrom ||
      (v.chrom == min_chrom && v.pos < min_pos) {
      min_chrom = v.chrom
      min_pos = v.pos
      found = true
    }
  }
  if !found {
    // All readers exhausted
    let final_has_line = []
    for i = 0; i < n; i = i + 1 {
      final_has_line.push(false)
    }
    return (
      BcfSrs::{
        readers: self.readers,
        has_line: final_has_line,
        pair_logic: self.pair_logic,
        errnum: self.errnum,
      },
      0,
    )
  }
  // Collect all variants at min (chrom, pos)
  // Group compatible variants together
  let has_line = []
  let new_readers = []
  for i = 0; i < n; i = i + 1 {
    has_line.push(false)
    new_readers.push(self.readers[i])
  }
  // First pass: mark readers that have a variant at (min_chrom, min_pos)
  let mut count = 0
  let mut ref_variant : SyncVariant? = None
  for i = 0; i < n; i = i + 1 {
    let reader = new_readers[i]
    if reader.pos < reader.buffer.length() {
      let v = reader.buffer[reader.pos]
      if v.chrom == min_chrom && v.pos == min_pos {
        match ref_variant {
          Some(rv) =>
            if variants_compatible(rv, v, self.pair_logic) {
              has_line[i] = true
              count = count + 1
            }
          None => {
            has_line[i] = true
            count = count + 1
            ref_variant = Some(v)
          }
        }
      }
    }
  }
  // Advance readers that have a line
  for i = 0; i < n; i = i + 1 {
    if has_line[i] {
      let reader = new_readers[i]
      new_readers[i] = SyncReader::{
        name: reader.name,
        buffer: reader.buffer,
        pos: reader.pos + 1,
      }
    }
  }
  // Also advance readers that are at a position before min (shouldn't happen
  // if sorted, but handle it)
  for i = 0; i < n; i = i + 1 {
    if !has_line[i] {
      let reader = new_readers[i]
      if reader.pos < reader.buffer.length() {
        let v = reader.buffer[reader.pos]
        if v.chrom == min_chrom && v.pos < min_pos {
          new_readers[i] = SyncReader::{
            name: reader.name,
            buffer: reader.buffer,
            pos: reader.pos + 1,
          }
        }
      }
    }
  }
  if count == 0 {
    // No compatible variants found, but there are variants at this position
    // from different readers. Advance all readers at this position and retry.
    for i = 0; i < n; i = i + 1 {
      let reader = new_readers[i]
      if reader.pos < reader.buffer.length() {
        let v = reader.buffer[reader.pos]
        if v.chrom == min_chrom && v.pos == min_pos {
          new_readers[i] = SyncReader::{
            name: reader.name,
            buffer: reader.buffer,
            pos: reader.pos + 1,
          }
        }
      }
    }
    // Recurse with updated state
    let updated = BcfSrs::{
      readers: new_readers,
      has_line: has_line,
      pair_logic: self.pair_logic,
      errnum: self.errnum,
    }
    return updated.next_line()
  }
  (
    BcfSrs::{
      readers: new_readers,
      has_line: has_line,
      pair_logic: self.pair_logic,
      errnum: self.errnum,
    },
    count,
  )
}

///|
/// Constructor for SyncVariant (for use from test package).
pub fn SyncVariant::make(
  chrom : String,
  pos : Int,
  ref_allele : String,
  alts : Array[String],
  id : String
) -> SyncVariant {
  SyncVariant::{ chrom, pos, ref_allele, alts, id }
}
