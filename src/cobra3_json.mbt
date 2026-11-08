///|
/// cobra3_json.mbt — JSON 解析/序列化与模型 JSON/YAML/TSV 交换.
///
/// YAML 采用 JSON 超集形式（flow style），保证 roundtrip 可解析；
/// TSV 为两段式：代谢物表 + 反应表（计量用 id:coeff 串）。

///|
pub(all) enum Cobra3Json {
  JNull
  JBool(Bool)
  JNum(Double)
  JStr(String)
  JArr(Array[Cobra3Json])
  JObj(Map[String, Cobra3Json])
}

///|
fn cobra3_json_skip(s : String, pos : Int32) -> Int32 {
  let mut i = pos
  while i < s.length() {
    let c = s[i].to_int()
    if c == 32 || c == 9 || c == 10 || c == 13 {
      i = i + 1
    } else {
      break
    }
  }
  i
}

///|
fn cobra3_json_hex4(s : String, pos : Int32) -> Int32? {
  if pos + 4 > s.length() {
    return None
  }
  let mut value = 0
  for k in 0..<4 {
    let c = s[pos + k].to_int()
    let digit = if c >= 48 && c <= 57 {
      c - 48
    } else if c >= 65 && c <= 70 {
      c - 55
    } else if c >= 97 && c <= 102 {
      c - 87
    } else {
      return None
    }
    value = value * 16 + digit
  }
  Some(value)
}

///|
fn cobra3_json_parse_string(s : String, pos : Int32) -> (String?, Int32) {
  if pos >= s.length() || s[pos] != '"' {
    return (None, -1)
  }
  let sb = StringBuilder()
  let mut i = pos + 1
  while i < s.length() {
    let ch = s[i]
    if ch == '"' {
      return (Some(sb.to_string()), i + 1)
    }
    if ch == '\\' {
      if i + 1 >= s.length() {
        return (None, -1)
      }
      let esc = s[i + 1]
      if esc == '"' {
        let _ = sb.write_char('"')
      } else if esc == '\\' {
        let _ = sb.write_char('\\')
      } else if esc == '/' {
        let _ = sb.write_char('/')
      } else if esc == 'n' {
        let _ = sb.write_char('\n')
      } else if esc == 't' {
        let _ = sb.write_char('\t')
      } else if esc == 'r' {
        let _ = sb.write_char('\r')
      } else if esc == 'b' {
        let _ = sb.write_char((8).to_char().unwrap_or(' '))
      } else if esc == 'f' {
        let _ = sb.write_char((12).to_char().unwrap_or(' '))
      } else if esc == 'u' {
        match cobra3_json_hex4(s, i + 2) {
          None => return (None, -1)
          Some(hi) => {
            let mut code = hi
            let mut advance = 4
            // UTF-16 代理对：高代理 + 低代理组合成一个码点
            if hi >= 55296 &&
              hi <= 56319 &&
              i + 12 <= s.length() &&
              s[i + 6] == '\\' &&
              s[i + 7] == 'u' {
              match cobra3_json_hex4(s, i + 8) {
                Some(lo) =>
                  if lo >= 56320 && lo <= 57343 {
                    code = 65536 + ((hi - 55296) << 10) + (lo - 56320)
                    advance = 10
                  }
                None => ()
              }
            }
            match code.to_char() {
              Some(ch) => {
                let _ = sb.write_char(ch)
              }
              None => {
                let _ = sb.write_char((0xfffd).to_char().unwrap_or('?'))
              }
            }
            i = i + advance
          }
        }
      } else {
        let _ = sb.write(s[i + 1:i + 2].to_owned())
      }
      i = i + 2
      continue
    }
    let _ = sb.write(s[i:i + 1].to_owned())
    i = i + 1
  }
  (None, -1)
}

///|
fn cobra3_json_parse_number(s : String, pos : Int32) -> (Cobra3Json?, Int32) {
  let mut i = pos
  while i < s.length() {
    let c = s[i].to_int()
    if (c >= 48 && c <= 57) ||
      c == 45 ||
      c == 43 ||
      c == 46 ||
      c == 101 ||
      c == 69 {
      i = i + 1
    } else {
      break
    }
  }
  if i == pos {
    return (None, -1)
  }
  match parse_double(s[pos:i].to_owned()) {
    Some(v) => (Some(JNum(v)), i)
    None => (None, -1)
  }
}

///|
fn cobra3_json_parse_value(s : String, pos : Int32) -> (Cobra3Json?, Int32) {
  let i = cobra3_json_skip(s, pos)
  if i >= s.length() {
    return (None, -1)
  }
  let ch = s[i]
  if ch == '"' {
    let (str_opt, next) = cobra3_json_parse_string(s, i)
    match str_opt {
      Some(v) => (Some(JStr(v)), next)
      None => (None, -1)
    }
  } else if ch == '{' {
    let obj : Map[String, Cobra3Json] = Map([], capacity=8)
    let mut p = cobra3_json_skip(s, i + 1)
    if p < s.length() && s[p] == '}' {
      (Some(JObj(obj)), p + 1)
    } else {
      let mut ok = true
      while p < s.length() {
        let (key, next) = cobra3_json_parse_string(s, p)
        match key {
          None => {
            ok = false
            break
          }
          Some(k) => {
            p = cobra3_json_skip(s, next)
            if p >= s.length() || s[p] != ':' {
              ok = false
              break
            }
            let (val, next2) = cobra3_json_parse_value(s, p + 1)
            match val {
              None => {
                ok = false
                break
              }
              Some(v) => {
                obj[k] = v
                p = cobra3_json_skip(s, next2)
                if p < s.length() && s[p] == ',' {
                  p = cobra3_json_skip(s, p + 1)
                } else if p < s.length() && s[p] == '}' {
                  return (Some(JObj(obj)), p + 1)
                } else {
                  ok = false
                  break
                }
              }
            }
          }
        }
      }
      if ok {
        (None, -1)
      } else {
        (None, -1)
      }
    }
  } else if ch == '[' {
    let arr : Array[Cobra3Json] = Array(capacity=8)
    let mut p = cobra3_json_skip(s, i + 1)
    if p < s.length() && s[p] == ']' {
      return (Some(JArr(arr)), p + 1)
    }
    while p < s.length() {
      let (val, next) = cobra3_json_parse_value(s, p)
      match val {
        None => return (None, -1)
        Some(v) => {
          arr.push(v)
          p = cobra3_json_skip(s, next)
          if p < s.length() && s[p] == ',' {
            p = cobra3_json_skip(s, p + 1)
          } else if p < s.length() && s[p] == ']' {
            return (Some(JArr(arr)), p + 1)
          } else {
            return (None, -1)
          }
        }
      }
    }
    (None, -1)
  } else if ch == 't' && i + 4 <= s.length() && s[i:i + 4].to_owned() == "true" {
    (Some(JBool(true)), i + 4)
  } else if ch == 'f' && i + 5 <= s.length() && s[i:i + 5].to_owned() == "false" {
    (Some(JBool(false)), i + 5)
  } else if ch == 'n' && i + 4 <= s.length() && s[i:i + 4].to_owned() == "null" {
    (Some(JNull), i + 4)
  } else {
    cobra3_json_parse_number(s, i)
  }
}

///|
pub fn cobra3_json_parse(text : String) -> Cobra3Json? {
  let (value, next) = cobra3_json_parse_value(text, 0)
  match value {
    None => None
    Some(v) =>
      if next >= 0 && cobra3_json_skip(text, next) == text.length() {
        Some(v)
      } else {
        None
      }
  }
}

///|
fn cobra3_json_get(obj : Map[String, Cobra3Json], key : String) -> Cobra3Json? {
  obj.get(key)
}

///|
fn cobra3_json_string(j : Cobra3Json) -> String? {
  match j {
    JStr(s) => Some(s)
    _ => None
  }
}

///|
fn cobra3_json_number(j : Cobra3Json) -> Double? {
  match j {
    JNum(v) => Some(v)
    _ => None
  }
}

///|
fn cobra3_json_object(j : Cobra3Json) -> Map[String, Cobra3Json]? {
  match j {
    JObj(v) => Some(v)
    _ => None
  }
}

///|
fn cobra3_json_escape(s : String) -> String {
  let sb = StringBuilder()
  let hex = "0123456789abcdef"
  for i in 0..<s.length() {
    let code = s[i].to_int()
    if code == 34 {
      let _ = sb.write("\\\"")
    } else if code == 92 {
      let _ = sb.write("\\\\")
    } else if code == 10 {
      let _ = sb.write("\\n")
    } else if code == 9 {
      let _ = sb.write("\\t")
    } else if code == 13 {
      let _ = sb.write("\\r")
    } else if code < 32 {
      let _ = sb.write(
        "\\u00" +
        hex[code / 16:code / 16 + 1].to_owned() +
        hex[code % 16:code % 16 + 1].to_owned(),
      )
    } else {
      let _ = sb.write(s[i:i + 1].to_owned())
    }
  }
  sb.to_string()
}

///|
fn cobra3_json_str(s : String) -> String {
  "\"" + cobra3_json_escape(s) + "\""
}

///|
fn cobra3_annotation_json(annotation : Map[String, String]) -> String {
  let sb = StringBuilder()
  let _ = sb.write("{")
  let keys : Array[String] = Array(capacity=annotation.length())
  for k in annotation.keys() {
    keys.push(k)
  }
  cobra3_sort_strings(keys)
  for i in 0..<keys.length() {
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write(
      cobra3_json_str(keys[i]) + ":" + cobra3_json_str(annotation[keys[i]]),
    )
  }
  let _ = sb.write("}")
  sb.to_string()
}

///|
pub fn cobra3_model_to_json(model : Cobra3Model) -> String {
  let sb = StringBuilder()
  let _ = sb.write("{")
  let _ = sb.write("\"id\":" + cobra3_json_str(model.id) + ",")
  let _ = sb.write("\"name\":" + cobra3_json_str(model.name) + ",")
  let _ = sb.write("\"compartments\":{")
  let comp_keys : Array[String] = Array(capacity=model.compartments.length())
  for k in model.compartments.keys() {
    comp_keys.push(k)
  }
  cobra3_sort_strings(comp_keys)
  for i in 0..<comp_keys.length() {
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write(
      cobra3_json_str(comp_keys[i]) +
      ":" +
      cobra3_json_str(model.compartments[comp_keys[i]]),
    )
  }
  let _ = sb.write("},")
  let _ = sb.write("\"metabolites\":[")
  for i in 0..<model.metabolites.length() {
    let m = model.metabolites[i]
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write("{\"id\":" + cobra3_json_str(m.id))
    let _ = sb.write(",\"name\":" + cobra3_json_str(m.name))
    let _ = sb.write(",\"formula\":" + cobra3_json_str(m.formula))
    let _ = sb.write(",\"charge\":" + m.charge.to_string())
    let _ = sb.write(",\"compartment\":" + cobra3_json_str(m.compartment))
    let _ = sb.write(
      ",\"boundary\":" + (if m.boundary { "true" } else { "false" }),
    )
    let _ = sb.write(
      ",\"annotation\":" + cobra3_annotation_json(m.annotation) + "}",
    )
  }
  let _ = sb.write("],")
  let _ = sb.write("\"genes\":[")
  for i in 0..<model.genes.length() {
    let g = model.genes[i]
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write("{\"id\":" + cobra3_json_str(g.id))
    let _ = sb.write(",\"name\":" + cobra3_json_str(g.name))
    let _ = sb.write(
      ",\"functional\":" + (if g.functional { "true" } else { "false" }),
    )
    let _ = sb.write(
      ",\"annotation\":" + cobra3_annotation_json(g.annotation) + "}",
    )
  }
  let _ = sb.write("],")
  let _ = sb.write("\"reactions\":[")
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write("{\"id\":" + cobra3_json_str(r.id))
    let _ = sb.write(",\"name\":" + cobra3_json_str(r.name))
    let _ = sb.write(
      ",\"lower_bound\":" +
      r.lower_bound.to_string() +
      ",\"upper_bound\":" +
      r.upper_bound.to_string() +
      ",\"objective\":" +
      r.objective.to_string(),
    )
    let _ = sb.write(",\"subsystem\":" + cobra3_json_str(r.subsystem))
    let gpr = match r.gpr {
      Some(g) => cobra3_gpr_to_string(g)
      None => ""
    }
    let _ = sb.write(",\"gpr\":" + cobra3_json_str(gpr))
    let _ = sb.write(",\"metabolites\":{")
    let parts = cobra3_reaction_participants(r)
    for k in 0..<parts.length() {
      if k > 0 {
        let _ = sb.write(",")
      }
      let _ = sb.write(
        cobra3_json_str(parts[k].0) + ":" + parts[k].1.to_string(),
      )
    }
    let _ = sb.write("}")
    let _ = sb.write(
      ",\"annotation\":" + cobra3_annotation_json(r.annotation) + "}",
    )
  }
  let _ = sb.write("]}")
  sb.to_string()
}

///|
fn cobra3_json_to_annotation(j : Cobra3Json?) -> Map[String, String] {
  let out : Map[String, String] = Map([], capacity=2)
  match j {
    Some(JObj(obj)) =>
      for k in obj.keys() {
        match cobra3_json_string(obj[k]) {
          Some(v) => out[k] = v
          None => ()
        }
      }
    _ => ()
  }
  out
}

///|
pub fn cobra3_model_from_json(text : String) -> Cobra3Model? {
  let root = match cobra3_json_parse(text) {
    Some(v) => v
    None => return None
  }
  let obj = match cobra3_json_object(root) {
    Some(v) => v
    None => return None
  }
  let id = match cobra3_json_get(obj, "id") {
    Some(JStr(s)) => s
    _ => "model"
  }
  let name = match cobra3_json_get(obj, "name") {
    Some(JStr(s)) => s
    _ => id
  }
  let mut model = cobra3_model_new(id)
  model = {
    id: model.id,
    name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions: model.reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  // compartments
  match cobra3_json_get(obj, "compartments") {
    Some(JObj(cm)) =>
      for k in cm.keys() {
        match cobra3_json_string(cm[k]) {
          Some(v) => model.compartments[k] = v
          None => ()
        }
      }
    _ => ()
  }
  // metabolites
  match cobra3_json_get(obj, "metabolites") {
    Some(JArr(arr)) =>
      for item in arr {
        match cobra3_json_object(item) {
          None =>
            match item {
              JStr(s) =>
                if s.length() > 0 {
                  model = cobra3_model_add_metabolite(
                    model,
                    cobra3_metabolite_new(s),
                  )
                }
              _ => ()
            }
          Some(mo) => {
            let mid = match cobra3_json_get(mo, "id") {
              Some(JStr(s)) => s
              _ => ""
            }
            if mid.length() == 0 {
              continue
            }
            let mname = match cobra3_json_get(mo, "name") {
              Some(JStr(s)) => s
              _ => mid
            }
            let formula = match cobra3_json_get(mo, "formula") {
              Some(JStr(s)) => s
              _ => ""
            }
            let charge = match cobra3_json_get(mo, "charge") {
              Some(JNum(v)) => v
              _ => 0.0
            }
            let compartment = match cobra3_json_get(mo, "compartment") {
              Some(JStr(s)) => s
              _ => "c"
            }
            let boundary = match cobra3_json_get(mo, "boundary") {
              Some(JBool(v)) => v
              _ => false
            }
            let met : Cobra3Metabolite = {
              id: mid,
              name: mname,
              formula,
              charge,
              compartment,
              boundary,
              annotation: cobra3_json_to_annotation(
                cobra3_json_get(mo, "annotation"),
              ),
            }
            model = cobra3_model_add_metabolite(model, met)
          }
        }
      }
    _ => ()
  }
  // genes
  match cobra3_json_get(obj, "genes") {
    Some(JArr(arr)) =>
      for item in arr {
        match cobra3_json_object(item) {
          None => ()
          Some(go) => {
            let gid = match cobra3_json_get(go, "id") {
              Some(JStr(s)) => s
              _ => ""
            }
            if gid.length() == 0 {
              continue
            }
            let gname = match cobra3_json_get(go, "name") {
              Some(JStr(s)) => s
              _ => gid
            }
            let functional = match cobra3_json_get(go, "functional") {
              Some(JBool(v)) => v
              _ => true
            }
            let gene : Cobra3Gene = {
              id: gid,
              name: gname,
              annotation: cobra3_json_to_annotation(
                cobra3_json_get(go, "annotation"),
              ),
              functional,
            }
            model = cobra3_model_add_gene(model, gene)
          }
        }
      }
    _ => ()
  }
  // reactions
  match cobra3_json_get(obj, "reactions") {
    Some(JArr(arr)) =>
      for item in arr {
        match cobra3_json_object(item) {
          None => ()
          Some(ro) => {
            let rid = match cobra3_json_get(ro, "id") {
              Some(JStr(s)) => s
              _ => ""
            }
            if rid.length() == 0 {
              continue
            }
            let rname = match cobra3_json_get(ro, "name") {
              Some(JStr(s)) => s
              _ => rid
            }
            let l_bound = match cobra3_json_get(ro, "lower_bound") {
              Some(JNum(v)) => v
              _ =>
                match cobra3_json_get(ro, "lb") {
                  Some(JNum(v)) => v
                  _ => 0.0
                }
            }
            let u_bound = match cobra3_json_get(ro, "upper_bound") {
              Some(JNum(v)) => v
              _ =>
                match cobra3_json_get(ro, "ub") {
                  Some(JNum(v)) => v
                  _ => 1000.0
                }
            }
            let objective = match cobra3_json_get(ro, "objective") {
              Some(JNum(v)) => v
              _ =>
                match cobra3_json_get(ro, "obj") {
                  Some(JNum(v)) => v
                  _ => 0.0
                }
            }
            let subsystem = match cobra3_json_get(ro, "subsystem") {
              Some(JStr(s)) => s
              _ => ""
            }
            let gpr_text = match cobra3_json_get(ro, "gpr") {
              Some(JStr(s)) => s
              _ => ""
            }
            let mets : Map[String, Double] = Map([], capacity=4)
            match cobra3_json_get(ro, "metabolites") {
              Some(JObj(mo)) =>
                for k in mo.keys() {
                  match cobra3_json_number(mo[k]) {
                    Some(v) => mets[k] = v
                    None => ()
                  }
                }
              _ => ()
            }
            let gpr_from_genes : Cobra3Gpr? = match
              cobra3_json_get(ro, "genes") {
              Some(JArr(gs)) => {
                let nodes : Array[Cobra3Gpr] = Array(capacity=gs.length())
                for g in gs {
                  match g {
                    JStr(s) => nodes.push(Cobra3Gpr::Gene(s))
                    _ => ()
                  }
                }
                if nodes.length() == 0 {
                  None
                } else {
                  Some(Cobra3Gpr::Or(nodes))
                }
              }
              _ => None
            }
            let reaction : Cobra3Reaction = {
              id: rid,
              name: rname,
              metabolites: mets,
              lower_bound: l_bound,
              upper_bound: u_bound,
              objective,
              subsystem,
              gpr: if gpr_text.length() == 0 {
                gpr_from_genes
              } else {
                cobra3_gpr_parse(gpr_text)
              },
              annotation: cobra3_json_to_annotation(
                cobra3_json_get(ro, "annotation"),
              ),
            }
            model = cobra3_model_add_reaction(model, reaction)
          }
        }
      }
    _ => ()
  }
  Some(model)
}

///|
/// YAML 采用 JSON 超集输出，读取复用 JSON 解析器.
pub fn cobra3_model_to_yaml(model : Cobra3Model) -> String {
  cobra3_model_to_json(model)
}

///|
pub fn cobra3_model_from_yaml(text : String) -> Cobra3Model? {
  cobra3_model_from_json(text)
}

///|
fn cobra3_tsv_escape(s : String) -> String {
  let sb = StringBuilder()
  for i in 0..<s.length() {
    let code = s[i].to_int()
    if code == 9 || code == 10 || code == 13 {
      let _ = sb.write_char(' ')
    } else {
      let _ = sb.write(s[i:i + 1].to_owned())
    }
  }
  sb.to_string()
}

///|
pub fn cobra3_model_to_tsv(model : Cobra3Model) -> String {
  let sb = StringBuilder()
  let _ = sb.write(
    "#metabolites\tid\tname\tformula\tcharge\tcompartment\tboundary\n",
  )
  for m in model.metabolites {
    let _ = sb.write(
      "M\t" +
      cobra3_tsv_escape(m.id) +
      "\t" +
      cobra3_tsv_escape(m.name) +
      "\t" +
      cobra3_tsv_escape(m.formula) +
      "\t" +
      m.charge.to_string() +
      "\t" +
      cobra3_tsv_escape(m.compartment) +
      "\t" +
      (if m.boundary { "1" } else { "0" }) +
      "\n",
    )
  }
  let _ = sb.write(
    "#reactions\tid\tname\tlower\tupper\tobjective\tgpr\tmetabolites\n",
  )
  for r in model.reactions {
    let parts = cobra3_reaction_participants(r)
    let sb_stoich = StringBuilder()
    for i in 0..<parts.length() {
      if i > 0 {
        let _ = sb_stoich.write(";")
      }
      let _ = sb_stoich.write(parts[i].0 + ":" + parts[i].1.to_string())
    }
    let gpr = match r.gpr {
      Some(g) => cobra3_gpr_to_string(g)
      None => ""
    }
    let _ = sb.write(
      "R\t" +
      cobra3_tsv_escape(r.id) +
      "\t" +
      cobra3_tsv_escape(r.name) +
      "\t" +
      r.lower_bound.to_string() +
      "\t" +
      r.upper_bound.to_string() +
      "\t" +
      r.objective.to_string() +
      "\t" +
      cobra3_tsv_escape(gpr) +
      "\t" +
      sb_stoich.to_string() +
      "\n",
    )
  }
  sb.to_string()
}

///|
fn cobra3_split_lines(text : String) -> Array[String] {
  let out : Array[String] = Array(capacity=16)
  let sb = StringBuilder()
  for i in 0..<text.length() {
    let c = text[i]
    if c == '\n' {
      out.push(sb.to_string())
      let _ = sb.reset()
    } else if c != '\r' {
      let _ = sb.write(text[i:i + 1].to_owned())
    }
  }
  if sb.to_string().length() > 0 {
    out.push(sb.to_string())
  }
  out
}

///|
fn cobra3_split_fields(line : String) -> Array[String] {
  let out : Array[String] = Array(capacity=8)
  let sb = StringBuilder()
  for i in 0..<line.length() {
    let c = line[i]
    if c == '\t' {
      out.push(sb.to_string())
      let _ = sb.reset()
    } else {
      let _ = sb.write(line[i:i + 1].to_owned())
    }
  }
  out.push(sb.to_string())
  out
}

///|
pub fn cobra3_model_from_tsv(text : String) -> Cobra3Model? {
  let lines = cobra3_split_lines(text)
  let mut model = cobra3_model_new("model")
  for line in lines {
    if line.length() == 0 || line.has_prefix("#") {
      continue
    }
    let f = cobra3_split_fields(line)
    if f.length() < 2 {
      continue
    }
    let kind = f[0]
    if kind == "M" && f.length() >= 7 {
      let met : Cobra3Metabolite = {
        id: f[1],
        name: f[2],
        formula: f[3],
        charge: parse_double(f[4]).unwrap_or(0.0),
        compartment: f[5],
        boundary: f[6] == "1",
        annotation: Map([], capacity=1),
      }
      model = cobra3_model_add_metabolite(model, met)
    } else if kind == "R" && f.length() >= 8 {
      let mets : Map[String, Double] = Map([], capacity=4)
      let pairs = f[7].split(";")
      for p in pairs {
        let s = p.to_owned()
        let idx = s.find(":")
        match idx {
          Some(k) => {
            let mid = s[0:k].to_owned()
            let coeff = parse_double(s[k + 1:].to_owned()).unwrap_or(0.0)
            if mid.length() > 0 {
              mets[mid] = coeff
            }
          }
          None => ()
        }
      }
      let reaction : Cobra3Reaction = {
        id: f[1],
        name: f[2],
        metabolites: mets,
        lower_bound: parse_double(f[3]).unwrap_or(0.0),
        upper_bound: parse_double(f[4]).unwrap_or(1000.0),
        objective: parse_double(f[5]).unwrap_or(0.0),
        subsystem: "",
        gpr: if f[6].length() == 0 {
          None
        } else {
          cobra3_gpr_parse(f[6])
        },
        annotation: Map([], capacity=1),
      }
      model = cobra3_model_add_reaction(model, reaction)
    }
  }
  Some(model)
}
