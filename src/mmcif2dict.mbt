///|
/// MMCIF2Dict - turn an mmCIF file into a flat dictionary.
///
/// Port of Biopython's Bio.PDB.MMCIF2Dict: a raw tokenizer (single/double
/// quoted strings, `;`-delimited multi-line blocks, `#` comments, `loop_`
/// expansion) that produces a flat map of token -> values. Every value is a
/// list of strings; the leading `data_` block marker maps to a one-element
/// list (Biopython stores a bare string there).

///|
pub suberror Mmcif2DictError {
  Mmcif2DictError(String)
}

///|
fn[T] m2d_fail(msg : String) -> T raise Mmcif2DictError {
  raise Mmcif2DictError(msg)
}

///|
fn m2d_is_ws(cp : Int) -> Bool {
  cp == 32 || cp == 9
}

///|
fn m2d_is_quote(cp : Int) -> Bool {
  cp == 39 || cp == 34
}

///|
fn m2d_strip(s : String) -> String {
  let mut start = 0
  while start < s.length() && m2d_is_ws(s.unsafe_get(start).to_int()) {
    start = start + 1
  }
  let mut end = s.length()
  while end > start && m2d_is_ws(s.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  s[start:end].to_owned()
}

///|
fn m2d_rstrip(s : String) -> String {
  let mut end = s.length()
  while end > 0 && m2d_is_ws(s.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  s[0:end].to_owned()
}

///|
fn m2d_join(arr : Array[String], sep : String) -> String {
  let buf = StringBuilder()
  for i = 0; i < arr.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(sep)
    }
    buf.write_string(arr[i])
  }
  buf.to_string()
}

///|
/// Split a document into lines (drop a final trailing newline, strip \r).
fn m2d_lines(content : String) -> Array[String] {
  let raw = split_by_char(content, '\n'.to_int())
  let lines : Array[String] = []
  let n = raw.length()
  for i = 0; i < n; i = i + 1 {
    let mut line = raw[i]
    if line.length() > 0 && line.unsafe_get(line.length() - 1).to_int() == 13 {
      line = line[0:line.length() - 1].to_owned()
    }
    if i == n - 1 && line.length() == 0 {
      // drop the trailing empty piece produced by a final newline
    } else {
      lines.push(line)
    }
  }
  lines
}

///|
/// Split one line into whitespace-delimited tokens, honouring single/double
/// quoted strings and mid-line `#` comments (Port of _splitline).
fn m2d_split_line(line : String) -> Array[String] raise Mmcif2DictError {
  let tokens : Array[String] = []
  let mut in_token = false
  let mut quote_open = -1
  let mut start_i = 0
  let n = line.length()
  let mut i = 0
  while i < n {
    let cp = line.unsafe_get(i).to_int()
    if m2d_is_ws(cp) {
      if in_token && quote_open < 0 {
        in_token = false
        tokens.push(line[start_i:i].to_owned())
      }
    } else if m2d_is_quote(cp) {
      if quote_open < 0 && !in_token {
        quote_open = cp
        in_token = true
        start_i = i + 1
      } else if cp == quote_open &&
        (i + 1 == n || m2d_is_ws(line.unsafe_get(i + 1).to_int())) {
        quote_open = -1
        in_token = false
        tokens.push(line[start_i:i].to_owned())
      }
    } else if cp == 35 && !in_token {
      // '#' comment: skip the rest of the line
      return tokens
    } else if !in_token {
      in_token = true
      start_i = i
    } else {
      // inside a token, keep accumulating
    }
    i = i + 1
  }
  if in_token {
    tokens.push(line[start_i:n].to_owned())
  }
  if quote_open >= 0 {
    m2d_fail("Line ended with quote open: \{line}")
  }
  tokens
}

///|
/// Tokenize the whole document (Port of _tokenize): skip `#` comment lines
/// and `;`-multi-line blocks (emitting each block as a single token).
fn m2d_tokenize(content : String) -> Array[String] raise Mmcif2DictError {
  let lines = m2d_lines(content)
  let tokens : Array[String] = []
  let mut empty = true
  let mut i = 0
  while i < lines.length() {
    let mut line = lines[i]
    i = i + 1
    empty = false
    if starts_with(line, 0, "#") {
      continue
    }
    if starts_with(line, 0, ";") {
      let buffer : Array[String] = []
      buffer.push(m2d_rstrip(line[1:line.length()].to_owned()))
      let mut closed = false
      while i < lines.length() {
        let l2 = lines[i]
        i = i + 1
        let l2s = m2d_rstrip(l2)
        if starts_with(l2s, 0, ";") {
          tokens.push(m2d_join(buffer, "\n"))
          line = l2s[1:l2s.length()].to_owned()
          if line.length() > 0 && !m2d_is_ws(line.unsafe_get(0).to_int()) {
            m2d_fail("Missing whitespace")
          }
          closed = true
          break
        }
        buffer.push(l2s)
      }
      if !closed {
        m2d_fail("Missing closing semicolon")
      }
    }
    let splitted = m2d_split_line(m2d_strip(line))
    for t in splitted {
      tokens.push(t)
    }
  }
  if empty {
    m2d_fail("Empty file.")
  }
  tokens
}

///|
/// Parse an mmCIF document into a flat `{token: [values]}` dictionary.
///
/// The `data_` block marker maps to a one-element list holding the block
/// name; `loop_` columns and single key-value pairs map to their value list.
pub fn mmcif2dict_parse(
  content : String,
) -> Map[String, Array[String]] raise Mmcif2DictError {
  let tokens = m2d_tokenize(content)
  let result : Map[String, Array[String]] = Map([], capacity=tokens.length() + 1)
  if tokens.length() == 0 {
    return result
  }
  let first = tokens[0]
  if !starts_with(first, 0, "data_") {
    m2d_fail("The input mmCIF file must begin with a 'data_' directive.")
  }
  result.set("data_", [first[5:first.length()].to_owned()])

  let mut loop_flag = false
  let mut key : String? = None
  let mut keys : Array[String] = []
  let mut col = 0
  let mut n_cols = 0
  let mut idx = 1
  while idx < tokens.length() {
    let token = tokens[idx]
    if to_lowercase(token) == "loop_" {
      loop_flag = true
      keys = []
      col = 0
      n_cols = 0
      idx = idx + 1
      continue
    } else if loop_flag {
      if starts_with(token, 0, "_") && (n_cols == 0 || col % n_cols == 0) {
        if col > 0 {
          loop_flag = false
          // fall through to the key-value handling below
        } else {
          result.set(token, [])
          keys.push(token)
          n_cols = n_cols + 1
          idx = idx + 1
          continue
        }
      } else {
        let k = keys[col % n_cols]
        let arr : Array[String] = match result.get(k) {
          Some(a) => a
          None => []
        }
        arr.push(token)
        result.set(k, arr)
        col = col + 1
        idx = idx + 1
        continue
      }
    }
    match key {
      None => key = Some(token)
      Some(k) => {
        result.set(k, [token])
        key = None
      }
    }
    idx = idx + 1
  }
  result
}