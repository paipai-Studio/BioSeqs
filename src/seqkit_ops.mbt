///|
/// SeqKit record-level operations: `head`/`range` (subset), `grep` (filter),
/// `sort` (by length/name), `rename` (renumber), `duplicate`/`rmdup`
/// (repeat/deduplicate), and `common` (shared IDs). All operate on parsed
/// `SeqRecord`s.

///|
/// Full FASTA header of a record (`id` plus `description` when present).
fn seqkit_full_name(rec : SeqRecord) -> String {
  if rec.description.length() == 0 {
    rec.id
  } else {
    rec.id + " " + rec.description
  }
}

///|
/// Naive substring test (returns true when `sub` is empty).
fn str_contains(s : String, sub : String) -> Bool {
  let m = sub.length()
  if m == 0 {
    return true
  }
  let n = s.length()
  let mut i = 0
  while i + m <= n {
    let mut j = 0
    while j < m && s.unsafe_get(i + j) == sub.unsafe_get(j) {
      j = j + 1
    }
    if j == m {
      return true
    }
    i = i + 1
  }
  false
}

///|
/// First `n` records (`seqkit head -n n`).
pub fn seqkit_head(records : Array[SeqRecord], n : Int) -> Array[SeqRecord] {
  let out : Array[SeqRecord] = Array::new()
  let limit = if n < records.length() { n } else { records.length() }
  let mut i = 0
  while i < limit {
    out.push(records[i])
    i = i + 1
  }
  out
}

///|
/// Records in the 1-based inclusive range `[start, end]` (`seqkit range`).
pub fn seqkit_range(
  records : Array[SeqRecord],
  start : Int,
  end : Int,
) -> Array[SeqRecord] {
  let out : Array[SeqRecord] = Array::new()
  let mut i = start - 1
  while i < end && i < records.length() {
    if i >= 0 {
      out.push(records[i])
    }
    i = i + 1
  }
  out
}

///|
/// Keep records whose full header contains `pattern`; `invert` keeps the
/// non-matching records instead (`seqkit grep`, substring match).
pub fn seqkit_grep(
  records : Array[SeqRecord],
  pattern : String,
  invert : Bool,
) -> Array[SeqRecord] {
  let out : Array[SeqRecord] = Array::new()
  for rec in records {
    let matched = str_contains(seqkit_full_name(rec), pattern)
    if matched != invert {
      out.push(rec)
    }
  }
  out
}

///|
/// Sort records by sequence length (`seqkit sort -l`; `descending=true` puts
/// the longest first, matching SeqKit's default).
pub fn seqkit_sort_length(
  records : Array[SeqRecord],
  descending : Bool,
) -> Array[SeqRecord] {
  let out = records.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 {
      let need_swap = if descending {
        out[j].seq.length() < key.seq.length()
      } else {
        out[j].seq.length() > key.seq.length()
      }
      if need_swap {
        out[j + 1] = out[j]
        j = j - 1
      } else {
        break
      }
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Sort records by `id` ascending (`seqkit sort -n`).
pub fn seqkit_sort_name(records : Array[SeqRecord]) -> Array[SeqRecord] {
  let out = records.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j].id > key.id {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Renumber records (`seqkit rename` default): IDs become `1`, `2`, … and the
/// description is cleared.
pub fn seqkit_renumber(records : Array[SeqRecord]) -> Array[SeqRecord] {
  let out : Array[SeqRecord] = Array::new()
  let mut i = 0
  while i < records.length() {
    let rec = records[i]
    let new_id = (i + 1).to_string()
    out.push(SeqRecord::{
      seq: rec.seq,
      id: new_id,
      name: new_id,
      description: "",
      letter_annotations: rec.letter_annotations,
      annotations: rec.annotations,
      features: rec.features,
    })
    i = i + 1
  }
  out
}

///|
/// Write each record `n` times (`seqkit duplicate -n n`).
pub fn seqkit_duplicate(records : Array[SeqRecord], n : Int) -> Array[SeqRecord] {
  let out : Array[SeqRecord] = Array::new()
  for rec in records {
    for _k = 0; _k < n; _k = _k + 1 {
      out.push(rec)
    }
  }
  out
}

///|
/// Remove records with duplicate sequence content, keeping the first occurrence
/// (`seqkit rmdup -s`).
pub fn seqkit_rmdup(records : Array[SeqRecord]) -> Array[SeqRecord] {
  let seen : Map[String, Bool] = Map([], capacity=0)
  let out : Array[SeqRecord] = Array::new()
  for rec in records {
    let s = rec.seq.to_string()
    match seen.get(s) {
      Some(_) => ()
      None => {
        seen.set(s, true)
        out.push(rec)
      }
    }
  }
  out
}

///|
/// Records of `a` whose full header also appears in `b` (`seqkit common`), in
/// `a`'s order.
pub fn seqkit_common(a : Array[SeqRecord], b : Array[SeqRecord]) -> Array[SeqRecord] {
  let b_names : Map[String, Bool] = Map([], capacity=0)
  for rec in b {
    b_names.set(seqkit_full_name(rec), true)
  }
  let out : Array[SeqRecord] = Array::new()
  for rec in a {
    match b_names.get(seqkit_full_name(rec)) {
      Some(_) => out.push(rec)
      None => ()
    }
  }
  out
}///|
/// Sliding-window subsequences (`seqkit sliding -W window -s step`,
/// non-greedy): windows `[i, i+window)` for i = 0, step, 2·step, … while they
/// fit; a partial trailing window is dropped.
pub fn seqkit_sliding(seq : String, window : Int, step : Int) -> Array[String] {
  let out : Array[String] = Array::new()
  let l = seq.length()
  let mut i = 0
  while i + window <= l {
    out.push(seq[i:i + window].to_owned())
    i = i + step
  }
  out
}

///|
/// Translate a DNA/RNA sequence to protein (`seqkit translate -f frame`),
/// forward frames 1/2/3 (frame 2/3 skip the leading 1/2 bases). Uses the
/// standard codon table via `Seq::translate`.
pub fn seqkit_translate(seq : String, frame : Int) -> String raise SeqError {
  let s = if frame <= 1 {
    seq
  } else {
    seq[frame - 1:seq.length()].to_owned()
  }
  Seq::new(s).translate().to_string()
}

///|
/// Point mutation: set the 1-based `pos`-th base to `base` (`seqkit mutate -p`).
pub fn seqkit_mutate_point(seq : String, pos : Int, base : String) -> String {
  let idx = pos - 1
  if idx < 0 || idx >= seq.length() {
    return seq
  }
  seq[0:idx].to_owned() + base + seq[idx + 1:seq.length()].to_owned()
}

///|
/// Deletion: remove the 1-based inclusive range `[start, end]` (`seqkit mutate
/// -d`).
pub fn seqkit_mutate_delete(seq : String, start : Int, end : Int) -> String {
  let s = start - 1
  if s < 0 || end > seq.length() || s >= end {
    return seq
  }
  seq[0:s].to_owned() + seq[end:seq.length()].to_owned()
}

///|
/// Insertion: insert `ins` before the 0-based index `pos` (`seqkit mutate -i`;
/// pos 0 = beginning, pos = length = end).
pub fn seqkit_mutate_insert(seq : String, pos : Int, ins : String) -> String {
  if pos < 0 || pos > seq.length() {
    return seq
  }
  seq[0:pos].to_owned() + ins + seq[pos:seq.length()].to_owned()
}

///|
fn seqkit_lcg(seed : Int) -> Int {
  (seed * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
/// Deterministically shuffle records (`seqkit shuffle -s seed`, LCG-based
/// Fisher–Yates).
pub fn seqkit_shuffle(records : Array[SeqRecord], seed : Int) -> Array[SeqRecord] {
  let out = records.copy()
  let mut s = seed
  let mut i = out.length() - 1
  while i > 0 {
    s = seqkit_lcg(s)
    let j = s % (i + 1)
    let tmp = out[i]
    out[i] = out[j]
    out[j] = tmp
    i = i - 1
  }
  out
}

///|
/// Randomly sample up to `n` records (`seqkit sample -n n`).
pub fn seqkit_sample(records : Array[SeqRecord], n : Int, seed : Int) -> Array[SeqRecord] {
  seqkit_head(seqkit_shuffle(records, seed), n)
}