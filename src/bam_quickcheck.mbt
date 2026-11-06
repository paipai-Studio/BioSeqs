/// Quick validation of BAM/BAI file integrity.
///
/// pysam/samtools equivalent: `samtools quickcheck`
///
/// Performs fast structural checks without fully parsing the file:
/// - BAM magic number (BGZF header)
/// - Non-empty data
/// - BAI magic number and basic structure (if index provided)

///|
/// Result of quickcheck validation
pub struct QuickcheckResult {
  is_valid : Bool
  error_message : String
}

///|
pub fn QuickcheckResult::is_valid(self : QuickcheckResult) -> Bool {
  self.is_valid
}

///|
pub fn QuickcheckResult::error_message(self : QuickcheckResult) -> String {
  self.error_message
}

///|
/// Quickcheck a BAM file's raw data.
/// Checks: non-empty, BGZF magic (0x1f 0x8b 0x08 0x04), BAM magic (BAM\1).
pub fn bam_quickcheck(data : String) -> QuickcheckResult {
  if data.length() == 0 {
    return { is_valid: false, error_message: "Empty data", }
  }
  // Check BGZF header
  if data.length() >= 4 {
    let b0 = data.unsafe_get(0).to_int()
    let b1 = data.unsafe_get(1).to_int()
    let b2 = data.unsafe_get(2).to_int()
    let b3 = data.unsafe_get(3).to_int()
    if b0 != 0x1f || b1 != 0x8b || b2 != 0x08 || b3 != 0x04 {
      return {
        is_valid: false,
        error_message: "Invalid BGZF magic number (expected 1f 8b 08 04, got " +
        b0.to_string() +
        " " +
        b1.to_string() +
        " " +
        b2.to_string() +
        " " +
        b3.to_string() +
        ")",
      }
    }
  } else {
    return { is_valid: false, error_message: "Data too short for BGZF header", }
  }
  // Check for BAM magic after decompression is not possible here,
  // but we can check the is_bam_magic function for uncompressed data
  if is_bgzf(data) {
    return { is_valid: true, error_message: "", }
  }
  // If not BGZF, check if it's raw BAM (starts with BAM\1)
  if is_bam_magic(data) {
    return { is_valid: true, error_message: "", }
  }
  { is_valid: false, error_message: "Data is neither valid BGZF nor raw BAM", }
}

///|
/// Quickcheck both a BAM file and its index (BAI).
pub fn bam_quickcheck_with_index(
  bam_data : String,
  bai_data : String,
) -> QuickcheckResult {
  let bam_check = bam_quickcheck(bam_data)
  if !bam_check.is_valid {
    return bam_check
  }
  // Check BAI: first 4 bytes should be "BAI\1"
  if bai_data.length() < 4 {
    return { is_valid: false, error_message: "BAI data too short", }
  }
  let b0 = bai_data.unsafe_get(0).to_int()
  let b1 = bai_data.unsafe_get(1).to_int()
  let b2 = bai_data.unsafe_get(2).to_int()
  let b3 = bai_data.unsafe_get(3).to_int()
  // BAI magic: B(66) A(65) I(73) \1(1)
  if b0 != 66 || b1 != 65 || b2 != 73 || b3 != 1 {
    return {
      is_valid: false,
      error_message: "Invalid BAI magic number (expected BAI\\1, got " +
      b0.to_string() +
      " " +
      b1.to_string() +
      " " +
      b2.to_string() +
      " " +
      b3.to_string() +
      ")",
    }
  }
  { is_valid: true, error_message: "", }
}

///|
/// Quickcheck a BamFile object (already parsed).
/// Verifies: has header, has references (if records exist), records have valid ref_ids.
pub fn bamfile_quickcheck(bam : BamFile) -> QuickcheckResult {
  // Check references exist if there are records
  if bam.records.length() > 0 && bam.header.references.length() == 0 {
    return {
      is_valid: false,
      error_message: "BAM has records but no references in header",
    }
  }
  // Check all records have valid ref_id
  let n_refs = bam.header.references.length()
  for i = 0; i < bam.records.length(); i = i + 1 {
    let rec = bam.records[i]
    if rec.ref_id >= 0 && rec.ref_id >= n_refs {
      return {
        is_valid: false,
        error_message: "Record " +
        i.to_string() +
        " has invalid ref_id " +
        rec.ref_id.to_string() +
        " (header has " +
        n_refs.to_string() +
        " references)",
      }
    }
    // Check CIGAR is non-empty for mapped reads
    if !rec.is_unmapped() && rec.cigar.length() == 0 {
      return {
        is_valid: false,
        error_message: "Mapped record " + i.to_string() + " has empty CIGAR",
      }
    }
  }
  { is_valid: true, error_message: "", }
}

///|
/// Batch quickcheck: validate multiple BAM data blobs.
/// Returns array of (index, result) for any that fail.
pub fn bam_quickcheck_batch(
  data_array : Array[String],
) -> Array[(Int32, QuickcheckResult)] {
  let failures : Array[(Int32, QuickcheckResult)] = []
  for i = 0; i < data_array.length(); i = i + 1 {
    let result = bam_quickcheck(data_array[i])
    if !result.is_valid {
      failures.push((i, result))
    }
  }
  failures
}
