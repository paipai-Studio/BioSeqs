// ============================================================================
// pysam_vcf_query — bcftools query 的内存实现
// ============================================================================

// ---------------------------------------------------------------------------
// 模板 token
// ---------------------------------------------------------------------------

///|
pub enum QueryTemplateToken {
  Lit(String)
  Field(String, Int32?)
  InfoKey(String, Int32?)
  FormatKey(String, String?)
}

///|
fn parse_int_opt(s : String) -> Int32? {
  // MoonBit 顶层 parse_int 直接可用
  if s.length() == 0 {
    return None
  }
  Some(parse_int(s))
}

///|
fn char_index_of(arr : Array[Char], ch : Char) -> Int32 {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == ch {
      return i
    }
    i = i + 1
  }
  -1
}

///|
fn str_index_of(s : String, c : Char) -> Int32 {
  let mut i = 0
  while i < s.length() {
    let ch : Char = s[i].unsafe_to_char()
    if ch == c {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Compile a query template string into an array of tokens.
fn join_chars(arr : Array[Char]) -> String {
  let mut out = ""
  let mut i = 0
  while i < arr.length() {
    out = out + arr[i].to_string()
    i = i + 1
  }
  out
}

///|
pub fn vcf_query_compile(template : String) -> Array[QueryTemplateToken] {
  let tokens : Array[QueryTemplateToken] = Array(capacity=0)
  let buf : Array[Char] = Array(capacity=0)
  let n = template.length()
  let mut i : Int32 = 0

  let push_lit = fn() {
    let s = join_chars(buf)
    if s.length() > 0 {
      tokens.push(QueryTemplateToken::Lit(s))
      buf.clear()
    }
  }

  while i < n {
    let c : Char = template[i].unsafe_to_char()
    if c == '%' {
      push_lit()
      i = i + 1
      // 收集标识符
      let ident_chars : Array[Char] = Array(capacity=0)
      while i < n {
        let ch : Char = template[i].unsafe_to_char()
        let is_ok = (ch >= 'a' && ch <= 'z') ||
          (ch >= 'A' && ch <= 'Z') ||
          (ch >= '0' && ch <= '9') ||
          ch == '/' ||
          ch == '{' ||
          ch == '}'
        if is_ok {
          ident_chars.push(ch)
          i = i + 1
        } else {
          break
        }
      }
      let ident : String = join_chars(ident_chars)

      // 处理 {i} 下标 — 手动扫描
      let mut idx : Int32? = None
      let brace_open = char_index_of(ident_chars, '{')
      let brace_close = char_index_of(ident_chars, '}')
      let clean_chars : Array[Char] = Array(capacity=ident_chars.length())
      if brace_open >= 0 && brace_close > brace_open {
        // 提取 {..} 中间的数字
        let mut num_str : String = ""
        let mut k = brace_open + 1
        while k < brace_close {
          num_str = num_str + ident_chars[k].to_string()
          k = k + 1
        }
        idx = parse_int_opt(num_str)
        // { 之前的部分
        let mut k2 = 0
        while k2 < brace_open {
          clean_chars.push(ident_chars[k2])
          k2 = k2 + 1
        }
      } else {
        // 没有下标，整个 ident 干净
        let mut k = 0
        while k < ident_chars.length() {
          clean_chars.push(ident_chars[k])
          k = k + 1
        }
      }
      let clean_ident : String = join_chars(clean_chars)

      // 分类：找第一个 /
      let slash_idx = str_index_of(clean_ident, '/')
      if slash_idx >= 0 {
        let prefix = clean_ident.substring(start=0, end=slash_idx)
        let rest = clean_ident.substring(start=slash_idx + 1)
        if prefix == "INFO" {
          tokens.push(QueryTemplateToken::InfoKey(rest, idx))
        } else if prefix == "FORMAT" {
          let slash2 = str_index_of(rest, '/')
          if slash2 >= 0 {
            let key = rest.substring(start=0, end=slash2)
            let sample = rest.substring(start=slash2 + 1)
            tokens.push(QueryTemplateToken::FormatKey(key, Some(sample)))
          } else {
            tokens.push(QueryTemplateToken::FormatKey(rest, None))
          }
        } else {
          tokens.push(QueryTemplateToken::Lit("%" + clean_ident))
        }
      } else {
        tokens.push(QueryTemplateToken::Field(clean_ident, idx))
      }
    } else if c == '\\' {
      i = i + 1
      if i < n {
        let next : Char = template[i].unsafe_to_char()
        if next == 'n' {
          buf.push('\n')
        } else if next == 't' {
          buf.push('\t')
        } else if next == '%' {
          buf.push('%')
        } else {
          buf.push(next)
        }
        i = i + 1
      }
    } else {
      buf.push(c)
      i = i + 1
    }
  }
  push_lit()
  tokens
}

// ---------------------------------------------------------------------------
// 渲染辅助
// ---------------------------------------------------------------------------

///|
fn alt_join(alt : Array[String]) -> String {
  let mut out = ""
  let mut i = 0
  while i < alt.length() {
    if i > 0 {
      out = out + ","
    }
    out = out + alt[i]
    i = i + 1
  }
  out
}

///|
fn qual_string(q : Double?) -> String {
  match q {
    Some(v) => v.to_string()
    None => "."
  }
}

///|
fn filter_string(f : Array[String]) -> String {
  if f.length() == 0 {
    return "PASS"
  }
  let mut out = ""
  let mut i = 0
  while i < f.length() {
    if i > 0 {
      out = out + ";"
    }
    out = out + f[i]
    i = i + 1
  }
  out
}

///|
fn id_string(id : Array[String]) -> String {
  if id.length() == 0 {
    return "."
  }
  let mut out = ""
  let mut i = 0
  while i < id.length() {
    if i > 0 {
      out = out + ";"
    }
    out = out + id[i]
    i = i + 1
  }
  out
}

///|
fn sample_names(vcf : VcfFile) -> Array[String] {
  match vcf.fixed_header {
    Some(fh) => {
      let mut idx : Int32 = -1
      let mut i : Int32 = 0
      while i < fh.columns.length() {
        if fh.columns[i] == "FORMAT" {
          idx = i
          break
        }
        i = i + 1
      }
      if idx < 0 {
        return []
      }
      let out : Array[String] = Array(capacity=0)
      let mut j = idx + 1
      while j < fh.columns.length() {
        out.push(fh.columns[j])
        j = j + 1
      }
      out
    }
    None => []
  }
}

///|
fn split_to_strings(s : String, sep : String) -> Array[String] {
  let raw = s.split(sep).to_array()
  let out : Array[String] = Array(capacity=raw.length())
  let mut i = 0
  while i < raw.length() {
    out.push(raw[i].to_owned())
    i = i + 1
  }
  out
}

///|
fn has_prefix(s : String, p : String) -> Bool {
  if p.length() > s.length() {
    return false
  }
  let mut i = 0
  while i < p.length() {
    if s[i] != p[i] {
      return false
    }
    i = i + 1
  }
  true
}

// ---------------------------------------------------------------------------
// 单 record 渲染
// ---------------------------------------------------------------------------

///|
fn vcf_query_render_row(
  vcf : VcfFile,
  rec : VcfRecord,
  tokens : Array[QueryTemplateToken],
) -> String {
  let samples = sample_names(vcf)
  let mut out = ""
  let mut ti = 0
  while ti < tokens.length() {
    let tok = tokens[ti]
    match tok {
      QueryTemplateToken::Lit(s) => out = out + s
      QueryTemplateToken::Field(name, idx_opt) =>
        match name {
          "CHROM" => out = out + rec.chrom
          "POS" => out = out + rec.pos.to_string()
          "REF" => out = out + rec.ref_seq
          "ALT" =>
            match idx_opt {
              Some(i) if i >= 1 && i <= rec.alt.length() =>
                out = out + rec.alt[i - 1]
              _ => out = out + alt_join(rec.alt)
            }
          "QUAL" => out = out + qual_string(rec.qual)
          "FILTER" => out = out + filter_string(rec.filter)
          "ID" => out = out + id_string(rec.id)
          "END" => out = out + rec.pos_end().to_string()
          "LEN" => out = out + rec.ref_seq.length().to_string()
          "NALT" => out = out + rec.alt.length().to_string()
          "SAMPLE" => {
            let mut s = ""
            let mut i = 0
            while i < samples.length() {
              if i > 0 {
                s = s + ","
              }
              s = s + samples[i]
              i = i + 1
            }
            out = out + s
          }
          _ => out = out + "."
        }
      QueryTemplateToken::InfoKey(key, idx_opt) =>
        match rec.info_get(key) {
          Some(val) =>
            match idx_opt {
              Some(i) => {
                let parts = split_to_strings(val, ",")
                let ii = i - 1
                if ii >= 0 && ii < parts.length() {
                  out = out + parts[ii]
                } else {
                  out = out + "."
                }
              }
              None => out = out + val
            }
          None => out = out + "."
        }
      QueryTemplateToken::FormatKey(key, sample_opt) =>
        match rec.format {
          Some(fmt) => {
            let mut key_idx : Int32 = -1
            let mut i : Int32 = 0
            while i < fmt.length() {
              if fmt[i] == key {
                key_idx = i
                break
              }
              i = i + 1
            }
            if key_idx >= 0 {
              match sample_opt {
                Some(sample) =>
                  match rec.samples.get(sample) {
                    Some(arr) =>
                      if key_idx < arr.length() {
                        out = out + arr[key_idx]
                      } else {
                        out = out + "."
                      }
                    None => out = out + "."
                  }
                None => {
                  let mut joined = ""
                  let mut j = 0
                  while j < samples.length() {
                    if j > 0 {
                      joined = joined + ","
                    }
                    match rec.samples.get(samples[j]) {
                      Some(arr) =>
                        if key_idx < arr.length() {
                          joined = joined + arr[key_idx]
                        } else {
                          joined = joined + "."
                        }
                      None => joined = joined + "."
                    }
                    j = j + 1
                  }
                  out = out + joined
                }
              }
            } else {
              out = out + "."
            }
          }
          None => out = out + "."
        }
    }
    ti = ti + 1
  }
  out
}

// ---------------------------------------------------------------------------
// 顶层 API
// ---------------------------------------------------------------------------

///|
pub fn vcf_query(vcf : VcfFile, template : String) -> String {
  let tokens = vcf_query_compile(template)
  let mut out = ""
  let mut ri = 0
  while ri < vcf.records.length() {
    let row = vcf_query_render_row(vcf, vcf.records[ri], tokens)
    out = out + row
    ri = ri + 1
  }
  out
}

///|
pub fn vcf_query_with_tokens(
  vcf : VcfFile,
  tokens : Array[QueryTemplateToken],
) -> String {
  let mut out = ""
  let mut ri = 0
  while ri < vcf.records.length() {
    let row = vcf_query_render_row(vcf, vcf.records[ri], tokens)
    out = out + row
    ri = ri + 1
  }
  out
}

///|
pub fn vcf_query_fields(
  vcf : VcfFile,
  field_names : Array[String],
) -> Array[Array[String]] {
  let samples = sample_names(vcf)
  let rows : Array[Array[String]] = Array(capacity=0)
  let mut ri = 0
  while ri < vcf.records.length() {
    let rec = vcf.records[ri]
    let row : Array[String] = Array(capacity=field_names.length())
    let mut fi = 0
    while fi < field_names.length() {
      let name = field_names[fi]
      let val : String = match name {
        "CHROM" => rec.chrom
        "POS" => rec.pos.to_string()
        "REF" => rec.ref_seq
        "ALT" => alt_join(rec.alt)
        "QUAL" => qual_string(rec.qual)
        "FILTER" => filter_string(rec.filter)
        "ID" => id_string(rec.id)
        "END" => rec.pos_end().to_string()
        "LEN" => rec.ref_seq.length().to_string()
        "NALT" => rec.alt.length().to_string()
        _ => {
          let info_prefix = "INFO/"
          let fmt_prefix = "FORMAT/"
          if has_prefix(name, info_prefix) {
            let key = name.substring(start=info_prefix.length())
            match rec.info_get(key) {
              Some(v) => v
              None => "."
            }
          } else if has_prefix(name, fmt_prefix) {
            let key = name.substring(start=fmt_prefix.length())
            match rec.format {
              Some(fmt) => {
                let mut key_idx : Int32 = -1
                let mut i : Int32 = 0
                while i < fmt.length() {
                  if fmt[i] == key {
                    key_idx = i
                    break
                  }
                  i = i + 1
                }
                if key_idx >= 0 {
                  let mut joined = ""
                  let mut j = 0
                  while j < samples.length() {
                    if j > 0 {
                      joined = joined + ","
                    }
                    match rec.samples.get(samples[j]) {
                      Some(arr) =>
                        if key_idx < arr.length() {
                          joined = joined + arr[key_idx]
                        } else {
                          joined = joined + "."
                        }
                      None => joined = joined + "."
                    }
                    j = j + 1
                  }
                  joined
                } else {
                  "."
                }
              }
              None => "."
            }
          } else {
            "."
          }
        }
      }
      row.push(val)
      fi = fi + 1
    }
    rows.push(row)
    ri = ri + 1
  }
  rows
}

///|
pub fn vcf_query_with_header(
  vcf : VcfFile,
  header_fields : Array[String],
) -> String {
  let rows = vcf_query_fields(vcf, header_fields)
  let mut out = ""
  let mut fi = 0
  while fi < header_fields.length() {
    if fi > 0 {
      out = out + "\t"
    }
    out = out + header_fields[fi]
    fi = fi + 1
  }
  out = out + "\n"
  let mut ri = 0
  while ri < rows.length() {
    let row = rows[ri]
    let mut ci = 0
    while ci < row.length() {
      if ci > 0 {
        out = out + "\t"
      }
      out = out + row[ci]
      ci = ci + 1
    }
    out = out + "\n"
    ri = ri + 1
  }
  out
}
