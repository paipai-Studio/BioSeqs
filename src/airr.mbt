///|
/// airr.mbt — AIRR 标准 Rearrangement TSV 读写.
///
/// 对齐 AIRR Community 的 MiAIRR / Rearrangement schema 核心列子集:
/// sequence_id / locus / v_call / d_call / j_call / cdr3 / cdr3_aa /
/// productive / v_identity / j_identity / junction_length.
/// 与 Change-O / Immcantation 的 Change-O `airr` 工具语义兼容,
/// 作为免疫组库(vdj.mbt)分析结果的交换格式. 数据不出本机.
pub suberror AirrError {
  AirrError(String)
}

///|
/// AIRR Rearrangement 记录(核心列子集).
pub(all) struct AirrRecord {
  sequence_id : String
  locus : String
  v_call : String
  d_call : String
  j_call : String
  cdr3 : String
  cdr3_aa : String
  productive : Bool
  v_identity : Double
  j_identity : Double
  junction_length : Int32
}

///|
/// AIRR TSV 的标准列序.
pub fn airr_columns() -> Array[String] {
  [
    "sequence_id", "locus", "v_call", "d_call", "j_call", "cdr3", "cdr3_aa", "productive",
    "v_identity", "j_identity", "junction_length",
  ]
}

///|
/// 由 VdjChain 转换为 AIRR 记录.
pub fn AirrRecord::from_vdj_chain(
  chain : VdjChain,
  locus~ : String = "",
) -> AirrRecord {
  {
    sequence_id: chain.sequence_id,
    locus,
    v_call: chain.v_call,
    d_call: chain.d_call,
    j_call: chain.j_call,
    cdr3: chain.cdr3_nt,
    cdr3_aa: chain.cdr3_aa,
    productive: chain.productive,
    v_identity: chain.v_identity,
    j_identity: chain.j_identity,
    junction_length: chain.cdr3_nt.length(),
  }
}

///|
/// 批量转换.
pub fn airr_from_vdj_chains(
  chains : Array[VdjChain],
  locus~ : String = "",
) -> Array[AirrRecord] {
  let out : Array[AirrRecord] = Array(capacity=0)
  for ch in chains {
    out.push(AirrRecord::from_vdj_chain(ch, locus~))
  }
  out
}

///|
fn airr_fmt_bool(b : Bool) -> String {
  if b {
    "T"
  } else {
    "F"
  }
}

///|
fn airr_fmt_double(x : Double) -> String {
  let s = x.to_string()
  if s.contains(".") ||
    s.contains("e") ||
    s.contains("E") ||
    s.contains("n") ||
    s.contains("a") {
    s
  } else {
    s + ".0"
  }
}

///|
fn airr_escape_field(s : String) -> String {
  if s.contains("\t") || s.contains("\n") || s.contains("\"") {
    "\"" + s.replace(old="\"", new="\"\"") + "\""
  } else {
    s
  }
}

///|
/// 写出 AIRR TSV 文本(含表头行).
pub fn airr_write(records : Array[AirrRecord]) -> String {
  let cols = airr_columns()
  let sb = StringBuilder()
  sb.write_string(cols.join("\t"))
  sb.write_char('\n')
  for r in records {
    let fields : Array[String] = [
      airr_escape_field(r.sequence_id),
      airr_escape_field(r.locus),
      airr_escape_field(r.v_call),
      airr_escape_field(r.d_call),
      airr_escape_field(r.j_call),
      airr_escape_field(r.cdr3),
      airr_escape_field(r.cdr3_aa),
      airr_fmt_bool(r.productive),
      airr_fmt_double(r.v_identity),
      airr_fmt_double(r.j_identity),
      r.junction_length.to_string(),
    ]
    sb.write_string(fields.join("\t"))
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
fn airr_unescape_field(s : String) -> String {
  if s.length() >= 2 && s.has_prefix("\"") && s.has_suffix("\"") {
    s[1:s.length() - 1].to_owned().replace(old="\"\"", new="\"")
  } else {
    s
  }
}

///|
fn airr_parse_bool(s : String) -> Bool {
  let t = trim(s)
  if t == "T" || t == "TRUE" || t == "true" || t == "1" {
    true
  } else {
    false
  }
}

///|
fn airr_parse_double(s : String) -> Double {
  parse_double(trim(s)).unwrap_or(0.0)
}

///|
fn airr_parse_int(s : String) -> Int32 {
  parse_int(trim(s))
}

///|
/// 解析 AIRR TSV 文本. 表头列序任意; 缺失列按空值容错.
pub fn airr_parse(text : String) -> Array[AirrRecord] raise AirrError {
  let content : String = text.replace(old="\r", new="")
  let raw_lines = content.split("\n").to_array()
  let lines : Array[String] = Array(capacity=0)
  for l in raw_lines {
    let t = trim(l.to_owned())
    if t.length() > 0 {
      lines.push(t)
    }
  }
  if lines.length() == 0 {
    return []
  }
  let header = airr_split_tsv(lines[0])
  let col_index : Map[String, Int32] = Map([], capacity=0)
  for i = 0; i < header.length(); i = i + 1 {
    col_index.set(trim(header[i]), i)
  }
  let out : Array[AirrRecord] = Array(capacity=0)
  for li = 1; li < lines.length(); li = li + 1 {
    let fields = airr_split_tsv(lines[li])
    if fields.length() == 0 ||
      (fields.length() == 1 && trim(fields[0]).length() == 0) {
      continue
    }
    let get = fn(col : String) -> String {
      match col_index.get(col) {
        Some(i) =>
          if i < fields.length() {
            airr_unescape_field(fields[i])
          } else {
            ""
          }
        None => ""
      }
    }
    let seq_id = get("sequence_id")
    if seq_id.length() == 0 {
      raise AirrError("row \{li}: missing sequence_id")
    }
    let productive_str = get("productive")
    let vi_str = get("v_identity")
    let ji_str = get("j_identity")
    let jl_str = get("junction_length")
    out.push({
      sequence_id: seq_id,
      locus: get("locus"),
      v_call: get("v_call"),
      d_call: get("d_call"),
      j_call: get("j_call"),
      cdr3: get("cdr3"),
      cdr3_aa: get("cdr3_aa"),
      productive: if productive_str.length() > 0 {
        airr_parse_bool(productive_str)
      } else {
        false
      },
      v_identity: if vi_str.length() > 0 {
        airr_parse_double(vi_str)
      } else {
        0.0
      },
      j_identity: if ji_str.length() > 0 {
        airr_parse_double(ji_str)
      } else {
        0.0
      },
      junction_length: if jl_str.length() > 0 {
        airr_parse_int(jl_str)
      } else {
        get("cdr3").length()
      },
    })
  }
  out
}

///|
/// 按制表符切分一行(带引号字段感知).
fn airr_split_tsv(line : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let n = line.length()
  let mut i = 0
  let mut sb = StringBuilder()
  let mut in_quotes = false
  while i < n {
    let c = line.unsafe_get(i).unsafe_to_char()
    if in_quotes {
      if c == '"' {
        if i + 1 < n && line.unsafe_get(i + 1).unsafe_to_char() == '"' {
          sb.write_char('"')
          i = i + 2
          continue
        }
        in_quotes = false
        i = i + 1
        continue
      }
      sb.write_char(c)
      i = i + 1
    } else if c == '"' && sb.to_string().length() == 0 {
      in_quotes = true
      i = i + 1
    } else if c == '\t' {
      out.push(sb.to_string())
      sb = StringBuilder()
      i = i + 1
    } else {
      sb.write_char(c)
      i = i + 1
    }
  }
  out.push(sb.to_string())
  out
}
