///|
/// SeqIO index extension (spec 5.4.3): the seqxml record scanner used by
/// the unified `seqio_index` entry point.
///
/// A seqxml record starts at the line carrying the `<entry` start tag and
/// ends at the line carrying the matching `</entry>` terminator. The key is
/// the entry's `id` attribute. A start tag without a terminator is a broken
/// record boundary and raises an error that pinpoints the offset (spec
/// 5.4.3.3).

///|
/// Scan seqxml content into index entries.
///
/// Records are delimited line-wise exactly like the other scanners: the raw
/// slice covers the `<entry ...>` start line through the `</entry>` line
/// (inclusive of its trailing newline when one follows). Duplicate
/// identifiers are resolved by the caller's position map (last one wins),
/// matching the fasta/genbank index behaviour.
fn seq_index_scan_seqxml(
  data : String,
  lines : Array[IndexedLine],
) -> Array[SeqIndexEntry] raise SeqIOError {
  let entries : Array[SeqIndexEntry] = Array(capacity=0)
  let n = lines.length()
  let mut i = 0
  while i < n {
    let line = lines[i]
    let trimmed = line.text.trim().to_owned()
    if trimmed.has_prefix("<entry") {
      let key = match seqxml_parse_entry_id(line.text) {
        Some(id) => id
        None =>
          raise SeqIOError(
            "seqxml record without an id attribute at offset " +
            line.start.to_string(),
          )
      }
      // Self-closing entry: the whole record lives on the start line.
      if trimmed.has_suffix("/>") {
        let end = if i + 1 < n { lines[i + 1].start } else { data.length() }
        entries.push(SeqIndexEntry::new(key, line.start, end - line.start))
        i = i + 1
        continue
      }
      // Find the `</entry>` terminator line.
      let mut j = i + 1
      let mut close_line = -1
      while j < n {
        let candidate = lines[j].text.trim().to_owned()
        if candidate.has_prefix("</entry>") {
          close_line = j
          break
        }
        // A nested or sibling start tag before the terminator means the
        // boundary is broken; report the offset of the unterminated record.
        if candidate.has_prefix("<entry") {
          break
        }
        j = j + 1
      }
      if close_line < 0 {
        raise SeqIOError(
          "Malformed seqxml record: unterminated <entry> (missing </entry>) at offset " +
          line.start.to_string(),
        )
      }
      let end = if close_line + 1 < n {
        lines[close_line + 1].start
      } else {
        data.length()
      }
      entries.push(SeqIndexEntry::new(key, line.start, end - line.start))
      i = close_line + 1
    } else {
      i = i + 1
    }
  }
  entries
}
