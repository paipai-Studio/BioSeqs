///|
/// A small JSON parser/encoder. BIOM 1.0 (and several other modules) read JSON,
/// and this provides a shared `JsonValue` tree instead of per-module ad-hoc
/// string extraction.

///|
pub enum JsonValue {
  JsonNull
  JsonBool(Bool)
  JsonNum(Double)
  JsonStr(String)
  JsonArr(Array[JsonValue])
  JsonObj(Map[String, JsonValue])
} derive(Show)

///|
/// Skip JSON whitespace (space, tab, LF, CR).
fn json_ws(s : String, pos : Int) -> Int {
  let mut p = pos
  while p < s.length() {
    let c = s.unsafe_get(p).to_int()
    if c == 32 || c == 9 || c == 10 || c == 13 {
      p = p + 1
    } else {
      break
    }
  }
  p
}

///|
fn json_parse_literal(
  s : String,
  pos : Int,
  lit : String,
  value : JsonValue,
) -> (JsonValue?, Int) {
  if pos + lit.length() > s.length() {
    return (None, pos)
  }
  if s[pos:pos + lit.length()] != lit {
    return (None, pos)
  }
  (Some(value), pos + lit.length())
}

///|
fn json_hex(s : String, pos : Int) -> Int? {
  if pos + 4 > s.length() {
    return None
  }
  let mut v = 0
  for i = 0; i < 4; i = i + 1 {
    let c = s.unsafe_get(pos + i).to_int()
    let d = if c >= 48 && c <= 57 {
      c - 48
    } else if c >= 97 && c <= 102 {
      c - 87
    } else if c >= 65 && c <= 70 {
      c - 55
    } else {
      return None
    }
    v = (v << 4) | d
  }
  Some(v)
}

///|
fn json_parse_string(s : String, pos : Int) -> (JsonValue?, Int) {
  // pos 指向开头的 `"`, 前进到内容
  let mut p = pos + 1
  let sb = StringBuilder()
  while p < s.length() {
    let c = s.unsafe_get(p).to_int()
    if c == 34 {
      return (Some(JsonStr(sb.to_string())), p + 1)
    }
    if c == 92 {
      p = p + 1
      if p >= s.length() {
        return (None, p)
      }
      let e = s.unsafe_get(p).to_int()
      if e == 34 {
        sb.write_char('"')
      } else if e == 92 {
        sb.write_char('\\')
      } else if e == 47 {
        sb.write_char('/')
      } else if e == 98 {
        sb.write_char('\b')
      } else if e == 102 {
        sb.write_char('\f')
      } else if e == 110 {
        sb.write_char('\n')
      } else if e == 114 {
        sb.write_char('\r')
      } else if e == 116 {
        sb.write_char('\t')
      } else if e == 117 {
        let h = json_hex(s, p + 1)
        match h {
          None => return (None, p)
          Some(code) => {
            sb.write_char(code.to_uint16().unsafe_to_char())
            p = p + 4
          }
        }
      } else {
        return (None, p)
      }
    } else {
      sb.write_char(c.to_uint16().unsafe_to_char())
    }
    p = p + 1
  }
  (None, p)
}

///|
fn json_parse_number(s : String, pos : Int) -> (JsonValue?, Int) {
  let mut p = pos
  let start = p
  if p < s.length() && s.unsafe_get(p).to_int() == 45 {
    p = p + 1
  }
  while p < s.length() {
    let c = s.unsafe_get(p).to_int()
    if c >= 48 && c <= 57 {
      p = p + 1
    } else {
      break
    }
  }
  if p < s.length() && s.unsafe_get(p).to_int() == 46 {
    p = p + 1
    while p < s.length() {
      let c = s.unsafe_get(p).to_int()
      if c >= 48 && c <= 57 {
        p = p + 1
      } else {
        break
      }
    }
  }
  if p < s.length() {
    let c = s.unsafe_get(p).to_int()
    if c == 101 || c == 69 {
      p = p + 1
      if p < s.length() && (s.unsafe_get(p).to_int() == 43 || s.unsafe_get(p).to_int() == 45) {
        p = p + 1
      }
      while p < s.length() {
        let d = s.unsafe_get(p).to_int()
        if d >= 48 && d <= 57 {
          p = p + 1
        } else {
          break
        }
      }
    }
  }
  let text = s[start:p].to_string()
  let v = parse_double(text)
  match v {
    None => (None, pos)
    Some(x) => (Some(JsonNum(x)), p)
  }
}

///|
fn json_parse_object(s : String, pos : Int) -> (JsonValue?, Int) {
  let mut p = json_ws(s, pos + 1)
  let m : Map[String, JsonValue] = Map([], capacity=8)
  if p < s.length() && s.unsafe_get(p).to_int() == 125 {
    return (Some(JsonObj(m)), p + 1)
  }
  while true {
    p = json_ws(s, p)
    if p >= s.length() || s.unsafe_get(p).to_int() != 34 {
      return (None, p)
    }
    let (k, p1) = json_parse_string(s, p)
    let key = match k {
      Some(JsonStr(x)) => x
      _ => return (None, p)
    }
    let p2 = json_ws(s, p1)
    if p2 >= s.length() || s.unsafe_get(p2).to_int() != 58 {
      return (None, p2)
    }
    let (v, p3) = json_parse_value(s, json_ws(s, p2 + 1))
    let val = match v {
      Some(x) => x
      None => return (None, p3)
    }
    m.set(key, val)
    let p4 = json_ws(s, p3)
    if p4 >= s.length() {
      return (None, p4)
    }
    let c = s.unsafe_get(p4).to_int()
    if c == 44 {
      p = p4 + 1
    } else if c == 125 {
      return (Some(JsonObj(m)), p4 + 1)
    } else {
      return (None, p4)
    }
  }
  (None, pos)
}

///|
fn json_parse_array(s : String, pos : Int) -> (JsonValue?, Int) {
  let mut p = json_ws(s, pos + 1)
  let arr : Array[JsonValue] = []
  if p < s.length() && s.unsafe_get(p).to_int() == 93 {
    return (Some(JsonArr(arr)), p + 1)
  }
  while true {
    let (v, p1) = json_parse_value(s, json_ws(s, p))
    match v {
      None => return (None, p1)
      Some(x) => arr.push(x)
    }
    let p2 = json_ws(s, p1)
    if p2 >= s.length() {
      return (None, p2)
    }
    let c = s.unsafe_get(p2).to_int()
    if c == 44 {
      p = p2 + 1
    } else if c == 93 {
      return (Some(JsonArr(arr)), p2 + 1)
    } else {
      return (None, p2)
    }
  }
  (None, pos)
}

///|
fn json_parse_value(s : String, pos : Int) -> (JsonValue?, Int) {
  let p = json_ws(s, pos)
  if p >= s.length() {
    return (None, p)
  }
  let c = s.unsafe_get(p).to_int()
  if c == 123 {
    json_parse_object(s, p)
  } else if c == 91 {
    json_parse_array(s, p)
  } else if c == 34 {
    json_parse_string(s, p)
  } else if c == 116 {
    json_parse_literal(s, p, "true", JsonBool(true))
  } else if c == 102 {
    json_parse_literal(s, p, "false", JsonBool(false))
  } else if c == 110 {
    json_parse_literal(s, p, "null", JsonNull)
  } else {
    json_parse_number(s, p)
  }
}

///|
/// Parse a complete JSON document into a value tree, or `None` on error.
pub fn parse_json(text : String) -> JsonValue? {
  let (v, p) = json_parse_value(text, 0)
  match v {
    None => None
    Some(x) => {
      let tail = json_ws(text, p)
      if tail == text.length() {
        Some(x)
      } else {
        None
      }
    }
  }
}

///|
/// Serialize a `JsonValue` back to compact JSON.
pub fn json_to_string(value : JsonValue) -> String {
  match value {
    JsonNull => "null"
    JsonBool(b) => if b { "true" } else { "false" }
    JsonNum(d) => {
      if d == d.to_int().to_double() && d.abs() < 9.007199254740992e15 {
        d.to_int().to_string()
      } else {
        d.to_string()
      }
    }
    JsonStr(s) => "\"" + json_escape_str(s) + "\""
    JsonArr(a) => {
      let sb = StringBuilder()
      sb.write_char('[')
      for i = 0; i < a.length(); i = i + 1 {
        if i > 0 {
          sb.write_char(',')
        }
        sb.write_string(json_to_string(a[i]))
      }
      sb.write_char(']')
      sb.to_string()
    }
    JsonObj(m) => {
      let sb = StringBuilder()
      sb.write_char('{')
      let mut first = true
      for kv in m {
        if !(first) {
          sb.write_char(',')
        }
        first = false
        sb.write_string("\"" + json_escape_str(kv.0) + "\":")
        sb.write_string(json_to_string(kv.1))
      }
      sb.write_char('}')
      sb.to_string()
    }
  }
}

///|
fn json_escape_str(s : String) -> String {
  let sb = StringBuilder()
  for i = 0; i < s.length(); i = i + 1 {
    let c = s.unsafe_get(i).to_int()
    if c == 34 {
      sb.write_string("\\\"")
    } else if c == 92 {
      sb.write_string("\\\\")
    } else if c == 10 {
      sb.write_string("\\n")
    } else if c == 13 {
      sb.write_string("\\r")
    } else if c == 9 {
      sb.write_string("\\t")
    } else {
      sb.write_char(c.to_uint16().unsafe_to_char())
    }
  }
  sb.to_string()
}