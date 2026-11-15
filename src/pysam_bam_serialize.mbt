// ============================================================================
// pysam_bam_serialize — BamFile → SAM 文本 序列化（支持 round-trip）
//
// 核心函数：bam_to_sam_text(bam) → String
// 反向函数：sam_to_bam（已存在于 bam.mbt）
// ============================================================================

fn tag_value_to_string(v : BamTagValue) -> String {
  match v {
    BamTagValue::Char(c) => c.to_string()
    BamTagValue::Int(i) => i.to_string()
    BamTagValue::Float(d) => d.to_string()
    BamTagValue::Str(s) => s
    BamTagValue::IntArray(arr) => {
      let mut out = ""
      let mut i = 0
      while i < arr.length() {
        if i > 0 { out = out + "," }
        out = out + arr[i].to_string()
        i = i + 1
      }
      out
    }
    BamTagValue::FloatArray(arr) => {
      let mut out = ""
      let mut i = 0
      while i < arr.length() {
        if i > 0 { out = out + "," }
        out = out + arr[i].to_string()
        i = i + 1
      }
      out
    }
  }
}

fn refs_to_string(refs : FixedArray[BamReference]) -> String {
  let mut out = ""
  let mut i = 0
  while i < refs.length() {
    out = out + "@SQ\tSN:" + refs[i].name + "\tLN:" + refs[i].length.to_string() + "\n"
    i = i + 1
  }
  out
}

fn flags_to_sam_string(header_text : String, refs : FixedArray[BamReference]) -> String {
  let mut out = ""
  // 如果 header_text 已有内容（包含 @HD, @RG, @PG 等），直接输出
  if header_text.length() > 0 {
    out = header_text
    // 确保最后有换行
    if out.length() > 0 {
      let last : Char = out[out.length() - 1].unsafe_to_char()
      if last != '\n' { out = out + "\n" }
    }
  } else {
    // 空 header：用默认 @HD
    out = "@HD\tVN:1.6\tSO:unsorted\n"
  }
  // 总是输出 refs_to_string（保证 @SQ 完整，不依赖 header.text）
  if refs.length() > 0 {
    out = out + refs_to_string(refs)
  }
  out
}

///|
/// 序列化一条 BamRecord 为一行 SAM。
/// 列顺序：QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
pub fn bam_record_to_sam_line(rec : BamRecord, refs : FixedArray[BamReference]) -> String {
  // QNAME
  let qname = rec.qname
  // FLAG
  let flag_str = rec.flag.to_string()
  // RNAME / ref_id → 名字
  let rname = if rec.ref_id >= 0 && rec.ref_id < refs.length() {
    refs[rec.ref_id].name
  } else if rec.ref_id == -1 {
    "*"
  } else {
    "*"
  }
  // POS: BAM pos 是 0-based，SAM 输出 1-based
  let pos_str = (rec.pos + 1).to_string()
  // MAPQ
  let mapq_str = rec.mapq.to_string()
  // CIGAR
  let cigar_str = if rec.cigar.length() == 0 {
    "*"
  } else {
    (fn() {
        let out : Array[BamCigarElement] = Array(capacity=rec.cigar.length())
        let mut i = 0
        while i < rec.cigar.length() { out.push(rec.cigar[i]); i = i + 1 }
        cigar_to_string(out)
      })()
  }
  // RNEXT (next ref). SAM "=" means mate on same ref; bam_to_sam uses
  // same-ref-id detection as canonical signal (more robust than -2).
  let rnext = if rec.next_ref_id == -1 {
    "*"
  } else if rec.next_ref_id == rec.ref_id {
    "="
  } else if rec.next_ref_id >= 0 && rec.next_ref_id < refs.length() {
    refs[rec.next_ref_id].name
  } else {
    "*"
  }
  // PNEXT: BAM 可能是 0-based；SAM 要求 1-based；若为 0 输出 0
  let pnext = if rec.next_pos <= 0 {
    "0"
  } else {
    (rec.next_pos + 1).to_string()
  }
  // TLEN
  let tlen_str = rec.tlen.to_string()
  // SEQ
  let seq_str = if rec.seq.length() > 0 { rec.seq } else { "*" }
  // QUAL
  let qual_str = if rec.qual.length() > 0 { rec.qual } else { "*" }

  // 组装前 11 列
  let mut out = qname + "\t" +
                flag_str + "\t" +
                rname + "\t" +
                pos_str + "\t" +
                mapq_str + "\t" +
                cigar_str + "\t" +
                rnext + "\t" +
                pnext + "\t" +
                tlen_str + "\t" +
                seq_str + "\t" +
                qual_str

  // TAGS: 每行 tab 追加；格式为 TAG:TYPE:VALUE
  let mut ti = 0
  while ti < rec.tags.length() {
    let tag = rec.tags[ti]
    let val_str = tag_value_to_string(tag.value)
    let type_char = tag.type_code.to_string()
    out = out + "\t" + tag.name + ":" + type_char + ":" + val_str
    ti = ti + 1
  }

  out
}

///|
/// 把 BamFile 完整序列化为 SAM 文本。
pub fn bam_to_sam_text(bam : BamFile) -> String {
  let mut out = flags_to_sam_string(bam.header.text, bam.header.references)
  let mut ri = 0
  while ri < bam.records.length() {
    out = out + bam_record_to_sam_line(bam.records[ri], bam.header.references) + "\n"
    ri = ri + 1
  }
  out
}

///|
/// Convenience alias.
pub fn write_sam(bam : BamFile) -> String {
  bam_to_sam_text(bam)
}
