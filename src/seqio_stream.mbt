///|
/// 流式(记录级惰性)序列解析.
///
/// 现有 parse_fasta/parse_fastq 会先切出全文行数组并构建完整
/// Array[SeqRecord], 大文件(如千万 reads 的 FASTQ)内存开销倍增.
/// 本模块以游标在原文上按偏移推进, 逐条构造记录, 每条仅复制自身内容,
/// 内存峰值 O(单条记录) 而非 O(全文件):
/// - `fasta_next_record` / `fastq_next_record`: 迭代式取下一条
/// - `fasta_each_record` / `fastq_each_record`: 回调式遍历
/// - `fasta_chunks` / `fastq_chunks`: 分批处理(批次间释放记录引用)
/// 解析语义与 parse_fasta/parse_fastq 对齐(含 CRLF 容错与相同报错).
pub struct SeqioCursor {
  mut pos : Int32
  mut line_no : Int32
}

///|
/// 新建游标(位于文本起始, 行号从 1 计).
pub fn SeqioCursor::new() -> SeqioCursor {
  { pos: 0, line_no: 0, }
}

///|
/// 读取下一行(不含换行符, 去行尾 \r), 到达文本末尾返回 None.
/// 返回的 view 指向原文, 不发生复制; 游标随之推进.
fn seqio_stream_next_line(
  content : String,
  cursor : SeqioCursor,
) -> StringView? {
  let n = content.length()
  if cursor.pos >= n {
    return None
  }
  let start = cursor.pos
  let rest = content[start:]
  let nl = rest.find("\n")
  let (line, next) = match nl {
    Some(i) => (rest[0:i], start + i + 1)
    None => (rest, n)
  }
  cursor.pos = next
  cursor.line_no = cursor.line_no + 1
  let ln = line.length()
  if ln > 0 && line.unsafe_get(ln - 1) == '\r' {
    Some(line[0:ln - 1])
  } else {
    Some(line)
  }
}

///|
/// 去除序列行中的空格与制表符(与 parse_fasta 一致).
fn seqio_stream_clean_seq_line(line : StringView) -> String {
  let sb = StringBuilder()
  for c in line {
    if c != ' ' && c != '\t' {
      sb.write_char(c)
    }
  }
  sb.to_string()
}

///|
/// 流式取下一条 FASTA 记录; 无更多记录返回 None.
/// 游标停在下一记录标题行行首, 可继续调用.
pub fn fasta_next_record(
  content : String,
  cursor : SeqioCursor,
) -> SeqRecord? raise SeqIOError {
  let mut title : String? = None
  let mut sb = StringBuilder()
  let mut finished = false
  while !finished {
    let line_start = cursor.pos
    let line_no = cursor.line_no
    match seqio_stream_next_line(content, cursor) {
      None => finished = true
      Some(line) => {
        if line.length() == 0 {
          continue
        }
        if line.unsafe_get(0) == '>' {
          if title is Some(_) {
            // 回退: 该标题行属于下一条记录
            cursor.pos = line_start
            cursor.line_no = line_no
            finished = true
          } else {
            title = Some(
              if line.length() > 1 {
                line[1:].to_owned()
              } else {
                ""
              },
            )
            sb = StringBuilder()
          }
        } else {
          sb.write_string(seqio_stream_clean_seq_line(line))
        }
      }
    }
  }
  match title {
    None => None
    Some(t) => Some(fasta_title_to_record(t, sb.to_string()))
  }
}

///|
/// 流式取下一条 FASTQ 记录(4 行一组); 无更多记录返回 None.
/// 报错携带行号; 序列与质量长度不一致时抛 SeqIOError.
pub fn fastq_next_record(
  content : String,
  cursor : SeqioCursor,
) -> SeqRecord? raise SeqIOError {
  // 跳过空行, 定位记录头
  let mut header : StringView? = None
  let mut exhausted = false
  while !exhausted {
    match seqio_stream_next_line(content, cursor) {
      None => exhausted = true
      Some(line) =>
        if line.length() > 0 {
          header = Some(line)
          exhausted = true
        }
    }
  }
  let header = match header {
    None => return None
    Some(h) => h
  }
  if header.unsafe_get(0) != '@' {
    raise SeqIOError(
      "Expected '@' at start of FASTQ header at line " +
      cursor.line_no.to_string() +
      ", got: " +
      header.to_owned(),
    )
  }
  // 序列行
  let seq_line = match seqio_stream_next_line(content, cursor) {
    None =>
      raise SeqIOError(
        "Truncated FASTQ record at line " + cursor.line_no.to_string(),
      )
    Some(l) => l
  }
  // 分隔行
  let plus_line = match seqio_stream_next_line(content, cursor) {
    None =>
      raise SeqIOError(
        "Truncated FASTQ record at line " + cursor.line_no.to_string(),
      )
    Some(l) => l
  }
  if plus_line.length() == 0 || plus_line.unsafe_get(0) != '+' {
    raise SeqIOError(
      "Expected '+' at start of FASTQ separator at line " +
      cursor.line_no.to_string() +
      ", got: " +
      plus_line.to_owned(),
    )
  }
  // 质量行
  let qual_line = match seqio_stream_next_line(content, cursor) {
    None =>
      raise SeqIOError(
        "Truncated FASTQ record at line " + cursor.line_no.to_string(),
      )
    Some(l) => l
  }
  if seq_line.length() != qual_line.length() {
    raise SeqIOError(
      "Sequence length " +
      seq_line.length().to_string() +
      " and quality length " +
      qual_line.length().to_string() +
      " don't match in FASTQ record at line " +
      cursor.line_no.to_string(),
    )
  }
  let title = if header.length() > 1 { header[1:].to_owned() } else { "" }
  let (id, description) = match title.find(" ") {
    Some(pos) => (title[0:pos].to_owned(), title)
    None => (title, title)
  }
  let quals : Array[Int32] = Array(capacity=qual_line.length())
  for c in qual_line {
    quals.push(c.to_int() - 33)
  }
  let letter_ann : Map[String, Array[Int32]] = Map([], capacity=1)
  letter_ann["phred_quality"] = quals
  Some(SeqRecord::{
    seq: Seq::new(seq_line.to_owned()),
    id,
    name: id,
    description,
    letter_annotations: letter_ann,
    annotations: Map([], capacity=0),
    features: [],
  })
}

///|
/// 回调式遍历全部 FASTA 记录(每条用毕即可被回收).
pub fn fasta_each_record(
  content : String,
  callback : (SeqRecord) -> Unit,
) -> Unit raise SeqIOError {
  let cursor = SeqioCursor::new()
  let mut finished = false
  while !finished {
    match fasta_next_record(content, cursor) {
      Some(record) => callback(record)
      None => finished = true
    }
  }
}

///|
/// 回调式遍历全部 FASTQ 记录.
pub fn fastq_each_record(
  content : String,
  callback : (SeqRecord) -> Unit,
) -> Unit raise SeqIOError {
  let cursor = SeqioCursor::new()
  let mut finished = false
  while !finished {
    match fastq_next_record(content, cursor) {
      Some(record) => callback(record)
      None => finished = true
    }
  }
}

///|
/// 分批遍历 FASTA 记录, 每批至多 chunk 条(批次间释放记录引用,
/// 内存峰值 O(chunk 条记录)). chunk <= 0 时按 1 条处理.
pub fn fasta_chunks(
  content : String,
  chunk : Int32,
  callback : (Array[SeqRecord]) -> Unit,
) -> Unit raise SeqIOError {
  let size = if chunk < 1 { 1 } else { chunk }
  let cursor = SeqioCursor::new()
  let batch : Array[SeqRecord] = Array(capacity=size)
  let mut finished = false
  while !finished {
    match fasta_next_record(content, cursor) {
      Some(record) => {
        batch.push(record)
        if batch.length() >= size {
          callback(batch.copy())
          batch.clear()
        }
      }
      None => {
        if batch.length() > 0 {
          callback(batch.copy())
          batch.clear()
        }
        finished = true
      }
    }
  }
}

///|
/// 分批遍历 FASTQ 记录, 语义同 fasta_chunks.
pub fn fastq_chunks(
  content : String,
  chunk : Int32,
  callback : (Array[SeqRecord]) -> Unit,
) -> Unit raise SeqIOError {
  let size = if chunk < 1 { 1 } else { chunk }
  let cursor = SeqioCursor::new()
  let batch : Array[SeqRecord] = Array(capacity=size)
  let mut finished = false
  while !finished {
    match fastq_next_record(content, cursor) {
      Some(record) => {
        batch.push(record)
        if batch.length() >= size {
          callback(batch.copy())
          batch.clear()
        }
      }
      None => {
        if batch.length() > 0 {
          callback(batch.copy())
          batch.clear()
        }
        finished = true
      }
    }
  }
}
