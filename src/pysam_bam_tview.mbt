///|
/// pysam_bam_tview — `samtools tview` 的内存版实现
///
/// 给定一个 BamFile + 可选 contig:position，返回浏览器风格的
/// read 对齐视图文本：参考序列行 + 每一行 read（带 CIGAR 展示）。
/// 这是 BioSeqs / Pysam 用户快速预览 region 内 read 布局的轻量工具。

///|
/// Render a tview window for `bam` between `[beg, end)` on `contig`.
///
/// - `bam`         : the BAM file (already parsed)
/// - `contig`      : target contig name; if "" all reads are shown
/// - `beg`, `end`  : half-open window on the contig. Auto-widened to 40+
/// - `ref_seq`     : optional reference sequence covering the window
/// - `max_reads`   : maximum reads to display (FIFO from window start)
pub fn bam_tview(
  bam : BamFile,
  contig : String,
  beg : Int32,
  end : Int32,
  ref_seq : String?,
  max_reads : Int32,
) -> String {
  let w_beg : Int32 = if end - beg < 40 { beg } else { beg }
  let w_end : Int32 = if end - beg < 40 { beg + 40 } else { end }
  let w_len = w_end - w_beg

  let sb = StringBuilder()

  // Ruler line: every 10th base, using decimal of the position
  let mut ruler_pos : Int32 = w_beg
  let mut first = true
  sb.write_string("    ")
  while ruler_pos < w_end {
    let mod100 = ruler_pos % 100
    let mod10 = ruler_pos % 10
    // Print tens digit of mod 100, or space for non-boundaries
    let digit : Char = if mod10 == 0 {
      let tens = (mod100 / 10).to_uint16().unsafe_to_char()
      tens
    } else {
      ' '
    }
    sb.write_char(digit)
    if !first && mod10 == 0 {
      sb.write_char('\n')
      sb.write_string("    ")
    }
    first = false
    ruler_pos = ruler_pos + 1
  }
  sb.write_char('\n')

  // REF row
  match ref_seq {
    Some(rs) => {
      sb.write_string("REF ")
      let mut i : Int32 = 0
      while i < w_len && i < rs.length() {
        sb.write_char(rs.unsafe_get(i).unsafe_to_char())
        i = i + 1
      }
      sb.write_char('\n')
    }
    None => ()
  }

  // Collect reads overlapping the window
  let reads : Array[BamRecord] = []
  let mut shown : Int32 = 0
  for rec in bam.records {
    if contig != "" {
      let rn = rec.get_reference_name(bam.header)
      if rn != contig { continue }
    }
    if rec.ref_id < 0 { continue }
    let rec_end = rec.pos + rec.infer_read_length()
    if rec.pos >= w_end || rec_end <= w_beg { continue }
    if max_reads > 0 && shown >= max_reads { break }
    reads.push(rec)
    shown = shown + 1
  }

  if reads.length() == 0 {
    sb.write_string("(no reads overlap this window)\n")
    return sb.to_string()
  }

  for rec in reads {
    render_tview_read(sb, rec, bam.header, w_beg, w_end)
    sb.write_char('\n')
  }

  sb.to_string()
}

///|
/// Render a single read as one tview row.
fn render_tview_read(
  sb : StringBuilder,
  rec : BamRecord,
  header : BamHeader,
  w_beg : Int32,
  w_end : Int32,
) -> Unit {
  // qname — pad/truncate to 10
  let qname : String = if rec.qname.length() > 10 {
    rec.qname[0:10].to_owned()
  } else {
    rec.qname
  }
  sb.write_string(qname)
  // pad
  let mut pad : Int32 = 10 - qname.length()
  while pad > 0 { sb.write_char(' '); pad = pad - 1 }

  let op_start : Int32 = match rec.query_alignment_start() {
    Some(s) => s
    None => 0
  }
  let op_end : Int32 = match rec.query_alignment_end() {
    Some(e) => e
    None => op_start
  }
  let read_seq : String = rec.query_alignment_sequence()

  // Skip CIGAR walking to bring ref_pos into [w_beg, w_end)
  let mut ref_pos : Int32 = rec.pos
  let mut q_pos : Int32 = 0
  while ref_pos < w_beg {
    let (sd, qd) = cigar_step(rec.cigar, q_pos)
    if sd == 0 && qd == 0 { break }
    ref_pos = ref_pos + sd
    q_pos = q_pos + qd
  }

  let mut col : Int32 = 0
  while ref_pos < w_end && col < (w_end - w_beg) {
    let (adv_ref, adv_q) = cigar_step(rec.cigar, q_pos)
    let ch : Char = if adv_ref > 0 && adv_q > 0 {
      let q_in_aligned = q_pos >= op_start && q_pos < op_end
      if q_in_aligned && (q_pos - op_start) < read_seq.length() {
        read_seq.unsafe_get(q_pos - op_start).unsafe_to_char()
      } else {
        '.'
      }
    } else if adv_ref > 0 {
      if cigar_is_skip(rec.cigar, q_pos) { 'N' } else { '-' }
    } else if adv_q > 0 {
      'I'
    } else {
      '.'
    }
    sb.write_char(ch)
    ref_pos = ref_pos + adv_ref
    q_pos = q_pos + adv_q
    col = col + 1
  }
}

///|
/// Advance through the CIGAR one element at a time, returning
/// (ref_delta, query_delta) at q_pos.
fn cigar_step(cigar : FixedArray[BamCigarElement], q_pos : Int32) -> (Int32, Int32) {
  let mut remaining : Int32 = q_pos
  for elem in cigar {
    let qlen = match elem.op {
      BamCigarOp::Match => elem.length
      BamCigarOp::Insert => elem.length
      BamCigarOp::SoftClip => elem.length
      BamCigarOp::Equal => elem.length
      BamCigarOp::Diff => elem.length
      _ => 0
    }
    let rlen = match elem.op {
      BamCigarOp::Match => elem.length
      BamCigarOp::Delete => elem.length
      BamCigarOp::Skip => elem.length
      BamCigarOp::Equal => elem.length
      BamCigarOp::Diff => elem.length
      _ => 0
    }
    if remaining < qlen {
      // We are inside this op. Return proportional advance.
      // For simplicity return either full step or zero.
      let left = if remaining < 0 { 0 } else { remaining }
      if qlen == 0 { return (rlen, 0) }
      return (
        ((left.to_double() / qlen.to_double()) * rlen.to_double()).to_int(),
        left,
      )
    }
    remaining = remaining - qlen
  }
  (0, 0)
}

///|
/// Is the operation covering q_pos a skip (N)?
fn cigar_is_skip(cigar : FixedArray[BamCigarElement], q_pos : Int32) -> Bool {
  let mut remaining : Int32 = q_pos
  for elem in cigar {
    let qlen = match elem.op {
      BamCigarOp::Match => elem.length
      BamCigarOp::Insert => elem.length
      BamCigarOp::SoftClip => elem.length
      BamCigarOp::Equal => elem.length
      BamCigarOp::Diff => elem.length
      _ => 0
    }
    if remaining < qlen {
      return elem.op == BamCigarOp::Skip
    }
    remaining = remaining - qlen
  }
  false
}
