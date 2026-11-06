/// Calculate MD and NM tags from alignment and reference sequence.
///
/// pysam/samtools equivalent: `samtools calmd` / `samtools fillmd`
///
/// The MD tag encodes the mismatching positions: a run of matches as a number,
/// followed by a mismatching base from the reference. The NM tag counts the
/// total number of mismatches (substitutions + insertions + deletions).

///|
/// Result of MD/NM tag calculation
pub struct MdNmResult {
  md_tag : String
  nm_count : Int32
}

///|
pub fn MdNmResult::md_tag(self : MdNmResult) -> String {
  self.md_tag
}

///|
pub fn MdNmResult::nm_count(self : MdNmResult) -> Int32 {
  self.nm_count
}

///|
/// Calculate MD and NM tags for a single BAM record given a reference sequence.
///
/// ref_seq: the reference sequence (ACGT string) covering the alignment region
/// ref_start: 0-based start position of ref_seq in the reference coordinate
///
/// Returns (MD_tag_string, NM_count).
pub fn bam_calmd(
  rec : BamRecord,
  ref_seq : String,
  ref_start : Int32,
) -> MdNmResult {
  let buf = StringBuilder()
  let mut nm : Int32 = 0
  let mut match_run : Int32 = 0
  let mut qpos : Int32 = 0
  let mut rpos = rec.pos - ref_start
  for el in rec.cigar {
    let oplen = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff =>
        for k = 0; k < oplen; k = k + 1 {
          if rpos >= 0 && rpos < ref_seq.length() && qpos < rec.seq.length() {
            let query_base = rec.seq.unsafe_get(qpos).to_int()
            let ref_base = ref_seq.unsafe_get(rpos).to_int()
            if el.op == BamCigarOp::Equal || query_base == ref_base {
              match_run = match_run + 1
            } else {
              if match_run > 0 {
                buf.write_string(match_run.to_string())
                match_run = 0
              }
              buf.write_char(Int32::unsafe_to_char(ref_base))
              nm = nm + 1
            }
          }
          qpos = qpos + 1
          rpos = rpos + 1
        }
      BamCigarOp::Insert | BamCigarOp::SoftClip => {
        if el.op == BamCigarOp::Insert {
          if match_run > 0 {
            buf.write_string(match_run.to_string())
            match_run = 0
          }
          buf.write_string("^")
          for k = 0; k < oplen; k = k + 1 {
            if qpos + k < rec.seq.length() {
              buf.write_char(
                Int32::unsafe_to_char(rec.seq.unsafe_get(qpos + k).to_int()),
              )
            }
          }
          nm = nm + oplen
        }
        qpos = qpos + oplen
      }
      BamCigarOp::Delete => {
        if match_run > 0 {
          buf.write_string(match_run.to_string())
          match_run = 0
        }
        buf.write_string("^")
        for k = 0; k < oplen; k = k + 1 {
          if rpos + k >= 0 && rpos + k < ref_seq.length() {
            buf.write_char(
              Int32::unsafe_to_char(ref_seq.unsafe_get(rpos + k).to_int()),
            )
          }
        }
        nm = nm + oplen
        rpos = rpos + oplen
      }
      BamCigarOp::Skip => rpos = rpos + oplen
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
    }
  }
  if match_run > 0 {
    buf.write_string(match_run.to_string())
  }
  { md_tag: buf.to_string(), nm_count: nm, }
}

///|
/// Add or update MD and NM tags on a BAM record, returning a new record.
pub fn bam_fillmd(
  rec : BamRecord,
  ref_seq : String,
  ref_start : Int32,
) -> BamRecord {
  let result = bam_calmd(rec, ref_seq, ref_start)
  let new_tags : Array[BamTag] = Array(capacity=0)
  let mut found_md = false
  let mut found_nm = false
  for t in rec.tags {
    if t.name == "MD" {
      new_tags.push({
        name: "MD",
        type_code: 'Z',
        value: BamTagValue::Str(result.md_tag),
      })
      found_md = true
    } else if t.name == "NM" {
      new_tags.push({
        name: "NM",
        type_code: 'i',
        value: BamTagValue::Int(result.nm_count),
      })
      found_nm = true
    } else {
      new_tags.push(t)
    }
  }
  if !found_md {
    new_tags.push({
      name: "MD",
      type_code: 'Z',
      value: BamTagValue::Str(result.md_tag),
    })
  }
  if !found_nm {
    new_tags.push({
      name: "NM",
      type_code: 'i',
      value: BamTagValue::Int(result.nm_count),
    })
  }
  { ..rec, tags: new_tags, }
}
