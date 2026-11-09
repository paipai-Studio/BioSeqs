/// pysam_bam_filter_expr.mbt — Filter reads by expression
/// Inspired by pysam/samtools view -e expression filtering

///|
/// Filter condition.
pub enum PysamFilterCond {
  MapqAtLeast(Int32) // mapq >= n
  MapqAtMost(Int32) // mapq <= n
  FlagHas(Int32) // flag & n != 0
  FlagNot(Int32) // flag & n == 0
  RnameIs(String) // rname == s
  RnameNot(String) // rname != s
  PosIn(Int32, Int32) // start <= pos < end
  SeqLenAtLeast(Int32) // seq.length >= n
  ContainsTag(String) // has tag
  NotContainsTag(String) // !has tag
  ProperPair // properly paired
  Primary // not secondary/supplementary
  Mapped // mapped
  Unmapped // unmapped
}

///|
/// Check if a SAM record matches a filter condition.
pub fn pysam_filter_match(
  fields : Array[StringView],
  cond : PysamFilterCond,
) -> Bool {
  if fields.length() < 11 {
    return false
  }
  let flag = parse_int(fields[1].to_owned())
  let rname = fields[2].to_owned()
  let pos = parse_int(fields[3].to_owned()) - 1
  let mapq = parse_int(fields[4].to_owned())
  let seq = fields[9].to_owned()
  match cond {
    PysamFilterCond::MapqAtLeast(n) => mapq >= n
    PysamFilterCond::MapqAtMost(n) => mapq <= n
    PysamFilterCond::FlagHas(n) => (flag & n) != 0
    PysamFilterCond::FlagNot(n) => (flag & n) == 0
    PysamFilterCond::RnameIs(s) => rname == s
    PysamFilterCond::RnameNot(s) => rname != s
    PysamFilterCond::PosIn(s, e) => pos >= s && pos < e
    PysamFilterCond::SeqLenAtLeast(n) => seq.length() >= n
    PysamFilterCond::ContainsTag(t) => {
      let mut found = false
      for i = 11; i < fields.length(); i = i + 1 {
        let tag = fields[i].to_owned()
        if tag.length() >= 2 && tag[0:2].to_owned() == t {
          found = true
          break
        }
      }
      found
    }
    PysamFilterCond::NotContainsTag(t) => {
      let mut found = false
      for i = 11; i < fields.length(); i = i + 1 {
        let tag = fields[i].to_owned()
        if tag.length() >= 2 && tag[0:2].to_owned() == t {
          found = true
          break
        }
      }
      !found
    }
    PysamFilterCond::ProperPair => (flag & 2) != 0
    PysamFilterCond::Primary => (flag & 256) == 0 && (flag & 2048) == 0
    PysamFilterCond::Mapped => (flag & 4) == 0
    PysamFilterCond::Unmapped => (flag & 4) != 0
  }
}

///|
/// Filter SAM by multiple conditions (AND logic).
pub fn pysam_filter_sam(
  sam_text : String,
  conds : Array[PysamFilterCond],
) -> String {
  let lines = sam_text.split("\n").to_array()
  let sb = StringBuilder()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      let _ = sb.write(line + "\n")
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let mut all_match = true
    for c = 0; c < conds.length(); c = c + 1 {
      if !pysam_filter_match(fields, conds[c]) {
        all_match = false
        break
      }
    }
    if all_match {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

///|
/// Filter SAM by any condition (OR logic).
pub fn pysam_filter_sam_any(
  sam_text : String,
  conds : Array[PysamFilterCond],
) -> String {
  let lines = sam_text.split("\n").to_array()
  let sb = StringBuilder()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      let _ = sb.write(line + "\n")
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let mut any_match = false
    for c = 0; c < conds.length(); c = c + 1 {
      if pysam_filter_match(fields, conds[c]) {
        any_match = true
        break
      }
    }
    if any_match {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

///|
/// Count reads matching conditions.
pub fn pysam_count_matching(
  sam_text : String,
  conds : Array[PysamFilterCond],
) -> Int32 {
  let filtered = pysam_filter_sam(sam_text, conds)
  let lines = filtered.split("\n").to_array()
  let mut count = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() > 0 && line[0].to_int() != 64 {
      count = count + 1
    }
  }
  count
}

///|
/// Suppress unused constructor warnings.
pub fn pysam_bam_unused_suppress_2() -> Unit {
  ignore(PysamFilterCond::MapqAtLeast(0))
  ignore(PysamFilterCond::MapqAtMost(0))
  ignore(PysamFilterCond::FlagHas(0))
  ignore(PysamFilterCond::FlagNot(0))
  ignore(PysamFilterCond::RnameIs(""))
  ignore(PysamFilterCond::RnameNot(""))
  ignore(PysamFilterCond::PosIn(0, 0))
  ignore(PysamFilterCond::SeqLenAtLeast(0))
  ignore(PysamFilterCond::ContainsTag(""))
  ignore(PysamFilterCond::NotContainsTag(""))
  ignore(PysamFilterCond::ProperPair)
  ignore(PysamFilterCond::Primary)
  ignore(PysamFilterCond::Mapped)
  ignore(PysamFilterCond::Unmapped)
}
