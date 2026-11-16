///|
/// Tabix (.tbi) index parsing and region querying - the random-access index
/// for coordinate-sorted, bgzip-compressed tab-delimited files (VCF, GFF,
/// BED, ...), as produced by `tabix`.
///
/// Port of the Rsamtools / htslib tabix reader. Tabix uses the same binning
/// scheme as the BAM index (BAI): a hierarchy of fixed-size bins plus a
/// linear index every 16384 bp. The per-reference binning structures are
/// therefore shared with `src/bam_index.mbt`.

///|
pub suberror TabixError {
  TabixError(String)
}

///|
fn[T] tbx_fail(msg : String) -> T raise TabixError {
  raise TabixError(msg)
}

///|
/// A parsed Tabix index.
///
/// `format` is 0 (generic), 1 (SAM) or 2 (VCF); `col_seq` / `col_beg` /
/// `col_end` are the 1-based file columns for sequence, start and end
/// (`col_end` is 0 when it defaults to `col_beg`); `meta` is the comment
/// character; `skip` the number of header lines; `names` the reference
/// sequence names. The per-reference binning is identical to BAI.
pub struct TabixIndex {
  n_ref : Int32
  format : Int32
  col_seq : Int32
  col_beg : Int32
  col_end : Int32
  meta : Int32
  skip : Int32
  names : Array[String]
  references : Array[BaiReference]
}

// ---------------------------------------------------------------------------
// Little-endian binary readers (mirror src/bam_index.mbt, raising TabixError)
// ---------------------------------------------------------------------------

///|
fn tbx_byte(data : String, pos : Int32) -> Int64 {
  (data.unsafe_get(pos).to_int() & 0xFF).to_int64()
}

///|
fn tbx_check_range(
  data : String,
  pos : Int32,
  len : Int32,
  what : String,
) -> Unit raise TabixError {
  if pos < 0 || pos + len > data.length() {
    tbx_fail(
      "Cannot read \{what} at offset \{pos}: premature end of data (length \{data.length()})",
    )
  }
}

///|
fn tbx_read_u32(data : String, pos : Int32) -> Int32 raise TabixError {
  tbx_check_range(data, pos, 4, "32-bit integer")
  (data.unsafe_get(pos).to_int() & 0xFF) |
  ((data.unsafe_get(pos + 1).to_int() & 0xFF) << 8) |
  ((data.unsafe_get(pos + 2).to_int() & 0xFF) << 16) |
  ((data.unsafe_get(pos + 3).to_int() & 0xFF) << 24)
}

///|
fn tbx_read_u64(data : String, pos : Int32) -> Int64 raise TabixError {
  tbx_check_range(data, pos, 8, "64-bit integer")
  tbx_byte(data, pos) |
  (tbx_byte(data, pos + 1) << 8) |
  (tbx_byte(data, pos + 2) << 16) |
  (tbx_byte(data, pos + 3) << 24) |
  (tbx_byte(data, pos + 4) << 32) |
  (tbx_byte(data, pos + 5) << 40) |
  (tbx_byte(data, pos + 6) << 48) |
  (tbx_byte(data, pos + 7) << 56)
}

///|
fn tbx_check_magic(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  data.unsafe_get(0).to_int() == 0x54 &&
  data.unsafe_get(1).to_int() == 0x42 &&
  data.unsafe_get(2).to_int() == 0x49 &&
  data.unsafe_get(3).to_int() == 0x01
}

///|
/// Split the NUL-terminated sequence names (total `l_nm` bytes).
fn tbx_parse_names(data : String, pos : Int32, l_nm : Int32) -> Array[String] {
  let names : Array[String] = []
  let mut cur = StringBuilder()
  for i = 0; i < l_nm; i = i + 1 {
    let byte = data.unsafe_get(pos + i).to_int() & 0xFF
    if byte == 0 {
      names.push(cur.to_string())
      cur = StringBuilder()
    } else {
      cur.write_char(data.unsafe_get(pos + i).unsafe_to_char())
    }
  }
  names
}

// ---------------------------------------------------------------------------
// Parsing and querying
// ---------------------------------------------------------------------------

///|
/// Parse a Tabix (.tbi) index from its raw bytes.
pub fn parse_tabix(data : String) -> TabixIndex raise TabixError {
  if !tbx_check_magic(data) {
    tbx_fail("Not a Tabix index: bad magic bytes")
  }
  let n_ref = tbx_read_u32(data, 4)
  let format = tbx_read_u32(data, 8)
  let col_seq = tbx_read_u32(data, 12)
  let col_beg = tbx_read_u32(data, 16)
  let col_end = tbx_read_u32(data, 20)
  let meta = tbx_read_u32(data, 24)
  let skip = tbx_read_u32(data, 28)
  let l_nm = tbx_read_u32(data, 32)
  let names = tbx_parse_names(data, 36, l_nm)
  let references : Array[BaiReference] = []
  let mut pos = 36 + l_nm
  for r = 0; r < n_ref; r = r + 1 {
    let n_bin = tbx_read_u32(data, pos)
    pos = pos + 4
    let bins : Map[Int32, Array[BaiChunk]] = Map([], capacity=n_bin)
    for b = 0; b < n_bin; b = b + 1 {
      let bin = tbx_read_u32(data, pos)
      pos = pos + 4
      let n_chunk = tbx_read_u32(data, pos)
      pos = pos + 4
      let chunks : Array[BaiChunk] = []
      for c = 0; c < n_chunk; c = c + 1 {
        let beg = tbx_read_u64(data, pos)
        let end = tbx_read_u64(data, pos + 8)
        pos = pos + 16
        chunks.push(BaiChunk::{ chunk_beg: beg, chunk_end: end, })
      }
      bins.set(bin, chunks)
    }
    let n_intv = tbx_read_u32(data, pos)
    pos = pos + 4
    let linear : Array[Int64] = []
    for iv = 0; iv < n_intv; iv = iv + 1 {
      linear.push(tbx_read_u64(data, pos))
      pos = pos + 8
    }
    references.push(BaiReference::{ bins, linear, })
  }
  TabixIndex::{
    n_ref,
    format,
    col_seq,
    col_beg,
    col_end,
    meta,
    skip,
    names,
    references,
  }
}

///|
/// Find the chunks of a Tabix index that overlap a half-open `[beg, end)`
/// interval on reference `ref_id`.
///
/// Returns the linear-index minimum offset together with every overlapping
/// chunk (BGZF virtual offsets), ready for `BgzfReader::seek`.
pub fn tabix_query(
  index : TabixIndex,
  ref_id : Int32,
  beg : Int32,
  end : Int32,
) -> BaiQueryResult raise TabixError {
  if ref_id < 0 || ref_id >= index.references.length() {
    tbx_fail("Reference id \{ref_id} is out of range (n_ref = \{index.n_ref})")
  }
  if beg < 0 || end < beg {
    tbx_fail("Invalid query interval [\{beg}, \{end})")
  }
  let reference = index.references[ref_id]
  let min_off = if reference.linear.length() > 0 {
    let iv = beg >> 14
    if iv < reference.linear.length() {
      reference.linear[iv]
    } else {
      reference.linear[reference.linear.length() - 1]
    }
  } else {
    0L
  }
  let chunks : Array[BaiChunk] = []
  let bins = bai_reg2bins(beg, end)
  for bin in bins {
    match reference.bins.get(bin) {
      Some(cs) =>
        for c in cs {
          chunks.push(c)
        }
      None => ()
    }
  }
  BaiQueryResult::{ min_off, chunks, }
}

// ============================================================================
// TabixFile — Pysam `pysam.TabixFile` style encapsulation
// ============================================================================

///|
/// In-memory TabixFile: holds a TabixIndex plus an array of parsed text lines
/// (already decompressed). The index maps region queries to chunks; the lines
/// are then scanned for overlap. Without a real BGZF reader this gives us
/// all the Tabix binning semantics for unit tests and small datasets.
pub struct TabixFile {
  index : TabixIndex
  /// Fully decompressed text lines (one entry per line, no newline).
  lines : Array[String]
  /// 1-based file columns for sequence, start, end (from index metadata).
  col_seq : Int32
  col_beg : Int32
  col_end : Int32
  /// Comment character (0 = none). Lines starting with this are skipped.
  meta : Char
  /// Number of header lines to skip before data (may be 0).
  skip : Int32
}

///|
/// Build a TabixFile from a TabixIndex plus decompressed text content. The
/// content is split on newlines and comment / header lines are recorded but
/// kept in the array (filtering is applied on fetch).
pub fn TabixFile::from_text(index : TabixIndex, content : String) -> TabixFile {
  let raw = content.split("\n").to_array()
  let lines : Array[String] = raw.map(fn(l) { l.to_owned() })
  TabixFile::{
    index,
    lines,
    col_seq: index.col_seq,
    col_beg: index.col_beg,
    col_end: index.col_end,
    meta: (index.meta & 0xFF).to_uint16().unsafe_to_char(),
    skip: index.skip,
  }
}

///|
/// Look up the contig id by name (matches Pysam `TabixFile.contigs`).
pub fn TabixFile::get_contig_id(self : TabixFile, name : String) -> Int32 {
  let mut i = 0
  while i < self.index.names.length() {
    if self.index.names[i] == name {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Contig names in declaration order. Matches Pysam `TabixFile.contigs`.
pub fn TabixFile::contigs(self : TabixFile) -> Array[String] {
  self.index.names.copy()
}

///|
/// Fetch all lines overlapping `[beg, end)` (0-based half-open) on
/// `contig_name`. Applies the Tabix bin overlap filter plus the per-line
/// start/end columns. Returns empty array on missing contig / no overlap.
/// Matches Pysam `TabixFile.fetch(contig, beg, end)`.
pub fn TabixFile::fetch(
  self : TabixFile,
  contig_name : String,
  beg : Int32,
  end : Int32,
) -> Array[String] raise TabixError {
  let ref_id = self.get_contig_id(contig_name)
  if ref_id < 0 {
    return []
  }
  // Run Tabix query to get candidate chunks (binning filter)
  let result = tabix_query(self.index, ref_id, beg, end)
  // For our in-memory model, chunks are synthetic so we just scan lines
  // within the header-skip zone and apply the per-line column filter.
  let out : Array[String] = Array(capacity=0)
  let mut line_idx = 0
  while line_idx < self.lines.length() {
    let line = self.lines[line_idx]
    // Skip header lines
    if line_idx < self.skip {
      line_idx = line_idx + 1
      continue
    }
    // Skip comment lines
    if line.length() > 0 && line.unsafe_get(0).unsafe_to_char() == self.meta {
      line_idx = line_idx + 1
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < self.col_beg {
      line_idx = line_idx + 1
      continue
    }
    // col_seq/col_beg/col_end are 1-based per Tabix spec
    let seq_name = fields[self.col_seq - 1].to_owned()
    if seq_name != contig_name {
      line_idx = line_idx + 1
      continue
    }
    let s = parse_int(fields[self.col_beg - 1].to_owned())
    let e = if self.col_end > 0 && self.col_end <= fields.length() {
      parse_int(fields[self.col_end - 1].to_owned())
    } else {
      s
    }
    if e > beg && s < end {
      out.push(line)
    }
    line_idx = line_idx + 1
  }
  out
}

///|
/// Fetch by samtools-style region string (1-based half-open). Matches Pysam
/// `TabixFile.fetch("chr1:100-200")`.
pub fn TabixFile::fetch_region(
  self : TabixFile,
  region : String,
) -> Array[String] raise TabixError {
  let (contig, start0, stop0) = parse_region(region)
  let zero : Int32 = 0
  let huge : Int32 = 1073741824
  let start = if start0 is Some(_) { start0.unwrap() } else { zero }
  let stop = if stop0 is Some(_) { stop0.unwrap() } else { huge }
  self.fetch(contig, start, stop)
}

// ============================================================================
// TBI serialization (write)
// ============================================================================

///|
/// Serialize a TabixIndex to its on-disk .tbi byte format (magic `TBI\x01` +
/// header + NUL-terminated names + per-reference binning + linear indices).
pub fn serialize_tabix(index : TabixIndex) -> String {
  let buf = StringBuilder()
  // Magic: 'T'=84, 'B'=66, 'I'=73, 1
  let t = 84
  let b = 66
  let i_ = 73
  let one = 1
  buf.write_char((t & 0xFF).to_uint16().unsafe_to_char())
  buf.write_char((b & 0xFF).to_uint16().unsafe_to_char())
  buf.write_char((i_ & 0xFF).to_uint16().unsafe_to_char())
  buf.write_char((one & 0xFF).to_uint16().unsafe_to_char())
  // Header u32 fields
  bai_write_u32(buf, index.n_ref)
  bai_write_u32(buf, index.format)
  bai_write_u32(buf, index.col_seq)
  bai_write_u32(buf, index.col_beg)
  bai_write_u32(buf, index.col_end)
  bai_write_u32(buf, index.meta)
  bai_write_u32(buf, index.skip)
  // NUL-terminated names blob
  let names_blob = StringBuilder()
  for n in index.names {
    names_blob.write_string(n)
    let nul = 0
    names_blob.write_char((nul & 0xFF).to_uint16().unsafe_to_char())
  }
  let names_str = names_blob.to_string()
  bai_write_u32(buf, names_str.length())
  buf.write_string(names_str)
  // Per-reference binning (identical layout to BAI)
  for ref_ in index.references {
    // Sorted bin keys
    let bin_keys : Array[Int32] = []
    for (k, _) in ref_.bins {
      bin_keys.push(k)
    }
    let n = bin_keys.length()
    let mut i = 0
    while i < n {
      let mut min_idx = i
      let mut j = i + 1
      while j < n {
        if bin_keys[j] < bin_keys[min_idx] {
          min_idx = j
        }
        j = j + 1
      }
      let tmp = bin_keys[i]
      bin_keys[i] = bin_keys[min_idx]
      bin_keys[min_idx] = tmp
      i = i + 1
    }
    bai_write_u32(buf, bin_keys.length())
    for bk in bin_keys {
      bai_write_u32(buf, bk)
      let chunks = ref_.bins[bk]
      bai_write_u32(buf, chunks.length())
      for c in chunks {
        bai_write_u64(buf, c.chunk_beg)
        bai_write_u64(buf, c.chunk_end)
      }
    }
    bai_write_u32(buf, ref_.linear.length())
    for off in ref_.linear {
      bai_write_u64(buf, off)
    }
  }
  buf.to_string()
}

// ============================================================================
// Memory-building a TabixIndex for a tab-delimited genomic file
// (samtools tabix -s -b -e equivalent, in-memory version)
// ============================================================================

///|
/// Build an in-memory TabixIndex from sorted tab-delimited lines. Mirrors
/// `samtools tabix -s col_seq -b col_beg -e col_end`. `col_seq/col_beg/col_end`
/// are **1-based**. The file must be sorted by (contig, start) for correct
/// binning.
pub fn build_tabix_index(
  lines : Array[String],
  col_seq : Int32,
  col_beg : Int32,
  col_end~ : Int32 = 0,
  meta~ : Char = '#',
  skip~ : Int32 = 0,
  format~ : Int32 = 0,
) -> TabixIndex {
  // Collect contig names in first-seen order
  let names : Array[String] = []
  let name_to_id : Map[String, Int32] = Map([], capacity=0)
  // Per-ref accumulators
  let ref_bins_acc : Array[Map[Int32, Array[BaiChunk]]] = []
  let ref_linear_last : Array[Int64] = []
  let ref_max_iv : Array[Int32] = []
  let mut next_voff : Int64 = 65536L
  let mut header_skip = skip
  let sentinel_v : Int64 = 9223372036854775807L
  for line in lines {
    // Skip headers / comments
    if header_skip > 0 ||
      (line.length() > 0 && line.unsafe_get(0).unsafe_to_char() == meta) {
      header_skip = header_skip - 1
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < col_beg {
      continue
    }
    let seq_name = fields[col_seq - 1].to_owned()
    let s = parse_int(fields[col_beg - 1].to_owned())
    let e = if col_end > 0 && col_end <= fields.length() {
      parse_int(fields[col_end - 1].to_owned())
    } else {
      s
    }
    // Resolve ref_id
    let ref_id = match name_to_id.get(seq_name) {
      Some(id) => id
      None => {
        let new_id = names.length()
        names.push(seq_name)
        name_to_id.set(seq_name, new_id)
        ref_bins_acc.push(Map([], capacity=0))
        ref_linear_last.push(sentinel_v)
        ref_max_iv.push(0)
        new_id
      }
    }
    let voff = next_voff
    let voff_end = voff + 256L
    next_voff = voff_end
    let chunk = BaiChunk::{ chunk_beg: voff, chunk_end: voff_end, }
    // Inline reg2bins here to avoid closure issue
    if s < e {
      let bins : Array[Int32] = []
      let e0 = e - 1
      bins.push(0)
      let mut k = 1 + (s >> 26)
      while k <= 1 + (e0 >> 26) {
        bins.push(k)
        k = k + 1
      }
      k = 9 + (s >> 23)
      while k <= 9 + (e0 >> 23) {
        bins.push(k)
        k = k + 1
      }
      k = 73 + (s >> 20)
      while k <= 73 + (e0 >> 20) {
        bins.push(k)
        k = k + 1
      }
      k = 585 + (s >> 17)
      while k <= 585 + (e0 >> 17) {
        bins.push(k)
        k = k + 1
      }
      k = 4681 + (s >> 14)
      while k <= 4681 + (e0 >> 14) {
        bins.push(k)
        k = k + 1
      }
      for bin in bins {
        let cur = ref_bins_acc[ref_id].get(bin)
        match cur {
          Some(arr) => {
            let copy = arr.copy()
            copy.push(chunk)
            ref_bins_acc[ref_id].set(bin, copy)
          }
          None => ref_bins_acc[ref_id].set(bin, [chunk])
        }
      }
    }
    let iv = s >> 14
    if iv > ref_max_iv[ref_id] {
      ref_max_iv[ref_id] = iv
    }
    if voff < ref_linear_last[ref_id] {
      ref_linear_last[ref_id] = voff
    }
  }
  // Assemble BaiReference per ref
  let references : Array[BaiReference] = Array(capacity=names.length())
  let mut i = 0
  while i < names.length() {
    let n_ivs = ref_max_iv[i] + 1
    let linear : Array[Int64] = []
    let mut last_good : Int64 = 0L
    let mut k = 0
    while k < n_ivs {
      if ref_linear_last[i] != sentinel_v {
        last_good = ref_linear_last[i]
      }
      linear.push(last_good)
      k = k + 1
    }
    references.push(BaiReference::{ bins: ref_bins_acc[i], linear, })
    i = i + 1
  }
  TabixIndex::{
    n_ref: names.length(),
    format,
    col_seq,
    col_beg,
    col_end,
    meta: meta.to_int(),
    skip,
    names,
    references,
  }
}
