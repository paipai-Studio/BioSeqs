/// Add read group (RG) tags to BAM records.
///
/// pysam/samtools equivalent: `samtools addrg` / `samtools view -r`

///|
/// Add a read group tag to all records in a BamFile.
/// If overwrite is false, records that already have an RG tag are left unchanged.
/// If overwrite is true, existing RG tags are replaced.
pub fn bam_addrg(
  bam : BamFile,
  read_group : String,
  overwrite~ : Bool = false,
) -> BamFile {
  let new_records : Array[BamRecord] = Array(capacity=0)
  for rec in bam.records {
    let mut has_rg = false
    let new_tags : Array[BamTag] = Array(capacity=0)
    for t in rec.tags {
      if t.name == "RG" {
        has_rg = true
        if overwrite {
          new_tags.push({
            name: "RG",
            type_code: 'Z',
            value: BamTagValue::Str(read_group),
          })
        } else {
          new_tags.push(t)
        }
      } else {
        new_tags.push(t)
      }
    }
    if !has_rg {
      new_tags.push({
        name: "RG",
        type_code: 'Z',
        value: BamTagValue::Str(read_group),
      })
    }
    new_records.push({ ..rec, tags: new_tags, })
  }
  { header: bam.header, records: new_records, }
}

///|
/// Add read group to records and update the header with an @RG line.
pub fn bam_addrg_with_header(
  bam : BamFile,
  read_group : String,
  sample : String,
  overwrite~ : Bool = false,
) -> BamFile {
  let updated = bam_addrg(bam, read_group, overwrite~)
  let rg_line = "@RG\tID:" + read_group + "\tSM:" + sample
  let mut new_text = updated.header.text
  if !new_text.contains(rg_line) {
    if new_text.length() > 0 && !new_text.has_suffix("\n") {
      new_text = new_text + "\n"
    }
    new_text = new_text + rg_line + "\n"
  }
  {
    header: { text: new_text, references: updated.header.references, },
    records: updated.records,
  }
}

///|
/// Add different read groups to records based on a mapping of qname prefix -> RG.
pub fn bam_addrg_by_prefix(
  bam : BamFile,
  prefix_map : Map[String, String],
  default_rg~ : String = "",
) -> BamFile {
  let new_records : Array[BamRecord] = Array(capacity=0)
  for rec in bam.records {
    let mut rg = default_rg
    let keys = prefix_map.keys().to_array()
    for k = 0; k < keys.length(); k = k + 1 {
      let prefix = keys[k]
      if rec.qname.length() >= prefix.length() &&
        rec.qname[0:prefix.length()].to_string() == prefix {
        rg = prefix_map[prefix]
        break
      }
    }
    let new_tags : Array[BamTag] = Array(capacity=0)
    let mut found = false
    for t in rec.tags {
      if t.name == "RG" {
        found = true
        if rg != "" {
          new_tags.push({
            name: "RG",
            type_code: 'Z',
            value: BamTagValue::Str(rg),
          })
        } else {
          new_tags.push(t)
        }
      } else {
        new_tags.push(t)
      }
    }
    if !found && rg != "" {
      new_tags.push({ name: "RG", type_code: 'Z', value: BamTagValue::Str(rg), })
    }
    new_records.push({ ..rec, tags: new_tags, })
  }
  { header: bam.header, records: new_records, }
}
