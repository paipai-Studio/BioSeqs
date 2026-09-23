///|
/// SeqIO binary format unification adapters.
///
/// The unified `seqio_parse` entry point carries binary payloads as base64
/// text (spec 4.3.1): every `parse_*_records` function in this module first
/// decodes the base64 content with strict validation and then delegates to
/// the existing standalone decoder for its format, normalizing the
/// decoder-specific types into `SeqRecord`. Decoding logic is never
/// duplicated and the legacy lenient `base64_decode` (used by BAM) is left
/// untouched.
///
/// Formats (design 1.2.3, binary rows):
///   - 2bit:       magic-validated, then twobit_from_bytes +
///                 twobit_to_seq_records
///   - sff:        sff_parse (None -> "Malformed SFF data"); one record per
///                 read with per-base qualities in
///                 letter_annotations["phred_quality"]
///   - sff-trim:   same decoder, sequence/qualities restricted to the read
///                 clip intervals
///   - snapgene:   snapgene_from_bytes + snapgene_to_seqrecord
///   - gck:        "GCK_FILE" magic-validated, then gck_from_bytes +
///                 gck_to_seqrecord
///   - xdna:       xdna_from_bytes, one record per packed record
///   - nib:        raw packed payload (4 bases per byte); the raw format
///                 carries no identifier, so the synthetic id "nib" is used
///
/// Empty content yields an empty array, not an error (spec 5.1.3.4). The
/// `bytes_to_base64` encoder emits unpadded output, so the strict decoder
/// accepts both padded and unpadded input while still rejecting invalid
/// characters and misplaced padding.

///|
/// Map a base64 character code to its 6-bit value (-1 when invalid).
fn b64_char_value(code : Int) -> Int {
  if code >= 65 && code <= 90 {
    // 'A'..'Z'
    code - 65
  } else if code >= 97 && code <= 122 {
    // 'a'..'z'
    code - 97 + 26
  } else if code >= 48 && code <= 57 {
    // '0'..'9'
    code - 48 + 52
  } else if code == 43 {
    // '+'
    62
  } else if code == 47 {
    // '/'
    63
  } else {
    -1
  }
}

///|
/// Strict base64 core: validates every significant character and the
/// padding layout, raising a plain SeqIOError message that the wrappers
/// extend with the format name. ASCII whitespace is skipped between
/// groups.
fn b64_decode_impl(content : String) -> Array[Int] raise SeqIOError {
  let sig : Array[Int] = Array(capacity=0)
  let n = content.length()
  for i in 0..<n {
    let c = content.unsafe_get(i).to_int()
    if c == 32 || c == 9 || c == 10 || c == 13 {
      continue
    }
    sig.push(c)
  }
  let m = sig.length()
  // Padding: at most two trailing '=' characters.
  let mut pad = 0
  while pad < 2 && m - 1 - pad >= 0 && sig[m - 1 - pad] == 61 {
    pad = pad + 1
  }
  // '=' must not appear anywhere before the trailing run.
  for i in 0..<(m - pad) {
    if sig[i] == 61 {
      raise SeqIOError("invalid base64 padding")
    }
  }
  let alpha = m - pad
  // A single character can never encode a byte group.
  if alpha % 4 == 1 {
    raise SeqIOError("invalid base64 content")
  }
  let out : Array[Int] = Array(capacity=0)
  let mut i = 0
  while i < alpha {
    let chars = if i + 4 <= alpha { 4 } else { alpha - i }
    let vals : Array[Int] = Array::make(4, 0)
    for k in 0..<chars {
      let v = b64_char_value(sig[i + k])
      if v < 0 {
        raise SeqIOError("invalid base64 content")
      }
      vals[k] = v
    }
    out.push((vals[0] << 2) | (vals[1] >> 4))
    if chars >= 3 {
      out.push(((vals[1] & 0x0F) << 4) | (vals[2] >> 2))
    }
    if chars == 4 {
      out.push(((vals[2] & 0x03) << 6) | vals[3])
    }
    i = i + chars
  }
  out
}

///|
/// Decode a strictly validated base64 string into raw bytes. Invalid
/// characters or misplaced padding raise SeqIOError instead of silently
/// producing wrong bytes (spec 4.3.1; the legacy `base64_decode` used by
/// BAM is intentionally untouched).
pub fn base64_decode_strict(content : String) -> Array[Int] raise SeqIOError {
  b64_decode_impl(content)
}

///|
/// Strict base64 decode with the format name embedded in error messages.
fn unify_b64_decode(
  content : String,
  format : String,
) -> Array[Int] raise SeqIOError {
  match (try? b64_decode_impl(content)) {
    Ok(bytes) => bytes
    Err(SeqIOError(msg)) => raise SeqIOError(msg + " for " + format)
  }
}

///|
/// Parse base64-carried UCSC 2bit content into SeqRecords. The magic
/// number is validated before decoding (spec 5.1.3.3).
pub fn parse_twobit_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let bytes = unify_b64_decode(content, "2bit")
  if bytes.length() < 4 || twobit_read_u32_le(bytes, 0) != twobit_magic {
    raise SeqIOError("2bit: invalid magic number (not a 2bit file)")
  }
  twobit_to_seq_records(twobit_from_bytes(bytes))
}

///|
/// Shared SFF normalization: one SeqRecord per read with per-base Phred
/// qualities in letter_annotations["phred_quality"]. With `trim` the
/// sequence and qualities are restricted to the read clip intervals.
fn unify_sff_records(
  content : String,
  trim : Bool,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let bytes = unify_b64_decode(content, if trim { "sff-trim" } else { "sff" })
  let file = match sff_parse(bytes) {
    Some(f) => f
    None => raise SeqIOError("Malformed SFF data")
  }
  let records : Array[SeqRecord] = Array(capacity=0)
  for read in file.sff_reads() {
    let record = if trim {
      SeqRecord::new(
        Seq::new(read.sff_clipped_bases()),
        id=read.sff_name(),
        name=read.sff_name(),
      )
    } else {
      SeqRecord::new(
        Seq::new(read.sff_bases()),
        id=read.sff_name(),
        name=read.sff_name(),
      )
    }
    record.letter_annotations["phred_quality"] = if trim {
      read.sff_clipped_qualities()
    } else {
      read.sff_qualities()
    }
    records.push(record)
  }
  records
}

///|
/// Parse base64-carried SFF content into SeqRecords (full sequences).
pub fn parse_sff_records(content : String) -> Array[SeqRecord] raise SeqIOError {
  unify_sff_records(content, false)
}

///|
/// Parse base64-carried SFF content into SeqRecords, restricting each
/// read's sequence and qualities to its clip intervals.
pub fn parse_sff_trim_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  unify_sff_records(content, true)
}

///|
/// Parse base64-carried SnapGene content into a single SeqRecord (the
/// format stores one sequence with its feature table).
pub fn parse_snapgene_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let bytes = unify_b64_decode(content, "snapgene")
  [snapgene_to_seqrecord(snapgene_from_bytes(bytes))]
}

///|
/// Parse base64-carried GCK content into a single SeqRecord. The
/// "GCK_FILE" magic is validated before decoding (spec 5.1.3.3).
pub fn parse_gck_records(content : String) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let bytes = unify_b64_decode(content, "gck")
  if bytes.length() < 14 || gck_ascii_to_string(bytes, 0, 8) != "GCK_FILE" {
    raise SeqIOError("gck: invalid magic number (not a GCK file)")
  }
  [gck_to_seqrecord(gck_from_bytes(bytes))]
}

///|
/// Parse base64-carried SnapGene-style xDNA content into SeqRecords, one
/// per packed record with the sequence type mirrored into annotations.
pub fn parse_xdna_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let bytes = unify_b64_decode(content, "xdna")
  let file = xdna_from_bytes(bytes)
  let records : Array[SeqRecord] = Array(capacity=0)
  for rec in file.records() {
    let record = SeqRecord::new(
      Seq::new(rec.sequence()),
      id=rec.name(),
      name=rec.name(),
    )
    unify_set_annotation(
      record.annotations,
      "sequence_type",
      rec.seq_type().description(),
    )
    unify_set_annotation(record.annotations, "comment", rec.annotations())
    records.push(record)
  }
  records
}

///|
/// Parse base64-carried nib content into a single SeqRecord. The payload
/// is the raw packed nib data (4 bases per byte); raw nib files carry no
/// identifier, so the synthetic id "nib" is assigned. Sequences whose
/// length is not a multiple of 4 lose padding semantics on round-trip
/// (the last byte's unused low bits decode as extra bases), matching the
/// plain nib packed layout.
pub fn parse_nib_records(content : String) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let bytes = unify_b64_decode(content, "nib")
  let hex = StringBuilder::new()
  for b in bytes {
    hex.write_char(nib_hex_digit((b >> 4) & 0xF))
    hex.write_char(nib_hex_digit(b & 0xF))
  }
  let nib = nib_from_hex("nib", hex.to_string(), 4 * bytes.length())
  let sequence = nib_decode(nib.packed_data(), nib.length())
  [SeqRecord::new(Seq::new(sequence), id=nib.seq_id(), name=nib.seq_id())]
}
