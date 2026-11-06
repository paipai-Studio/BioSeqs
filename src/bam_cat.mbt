/// Concatenate BAM files (like `samtools cat`).
///
/// All input BAM files must have compatible headers (same references).
/// The output uses the first file's header and merges all records.

///|
/// Concatenate multiple BamFiles into one.
/// Uses the first file's header; all records are appended in order.
pub fn bam_cat(files : Array[BamFile]) -> BamFile {
  if files.length() == 0 {
    return { header: { text: "", references: [], }, records: [], }
  }
  let header = files[0].header
  let all_records : Array[BamRecord] = Array(capacity=0)
  for f in files {
    for r in f.records {
      all_records.push(r)
    }
  }
  { header, records: all_records, }
}

///|
/// Concatenate with header validation.
/// Returns an error if reference names/lengths don't match.
pub fn bam_cat_safe(files : Array[BamFile]) -> BamFile raise BamError {
  if files.length() == 0 {
    return { header: { text: "", references: [], }, records: [], }
  }
  let base_header = files[0].header
  let base_names = files[0].get_reference_names()
  for fi = 1; fi < files.length(); fi = fi + 1 {
    let other_names = files[fi].get_reference_names()
    if other_names.length() != base_names.length() {
      raise BamError("Reference count mismatch in file")
    }
    for j = 0; j < base_names.length(); j = j + 1 {
      if other_names[j] != base_names[j] {
        raise BamError("Reference name mismatch")
      }
      match files[0].get_reference_length(base_names[j]) {
        Some(l1) =>
          match files[fi].get_reference_length(base_names[j]) {
            Some(l2) =>
              if l1 != l2 {
                raise BamError("Reference length mismatch")
              }
            None => raise BamError("Missing reference in file")
          }
        None => ()
      }
    }
  }
  let all_records : Array[BamRecord] = Array(capacity=0)
  for f in files {
    for r in f.records {
      all_records.push(r)
    }
  }
  { header: base_header, records: all_records, }
}
