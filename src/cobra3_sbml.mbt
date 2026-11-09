///|
/// cobra3_sbml.mbt — SBML Level 3 Version 2 + FBC v2 读写.
///
/// 写入包含 compartment/species/反应计量/FBC fluxBounds/FBC geneProduct/GPR/objective；
/// 读取兼容本模块输出与常见 FBC 风格（bounds 引用 parameter 或直接数值）.
fn cobra3_xml_escape(s : String) -> String {
  let sb = StringBuilder()
  for i in 0..<s.length() {
    let code = s[i].to_int()
    if code == 38 {
      let _ = sb.write("&amp;")
    } else if code == 60 {
      let _ = sb.write("&lt;")
    } else if code == 62 {
      let _ = sb.write("&gt;")
    } else if code == 34 {
      let _ = sb.write("&quot;")
    } else if code == 39 {
      let _ = sb.write("&apos;")
    } else {
      let _ = sb.write(s[i:i + 1].to_owned())
    }
  }
  sb.to_string()
}

///|
fn cobra3_xml_unescape(s : String) -> String {
  let mut out = s
  out = out.replace(old="&lt;", new="<")
  out = out.replace(old="&gt;", new=">")
  out = out.replace(old="&quot;", new="\"")
  out = out.replace(old="&apos;", new="'")
  out = out.replace(old="&amp;", new="&")
  out
}

///|
fn cobra3_xml_safe_id(id : String) -> String {
  let sb = StringBuilder()
  for i in 0..<id.length() {
    let code = id[i].to_int()
    if (code >= 97 && code <= 122) ||
      (code >= 65 && code <= 90) ||
      (code >= 48 && code <= 57) ||
      code == 95 ||
      code == 45 ||
      code == 46 {
      let _ = sb.write(id[i:i + 1].to_owned())
    } else {
      let _ = sb.write("_")
    }
  }
  let s = sb.to_string()
  if s.length() == 0 {
    "_"
  } else {
    s
  }
}

///|
fn cobra3_xml_attr(tag : String, name : String) -> String? {
  let needle = " " + name + "=\""
  match tag.find(needle) {
    None => None
    Some(idx) => {
      let start = idx + needle.length()
      match tag[start:].find("\"") {
        None => None
        Some(end) =>
          Some(cobra3_xml_unescape(tag[start:start + end].to_owned()))
      }
    }
  }
}

///|
fn cobra3_xml_gpr(gpr : Cobra3Gpr) -> String {
  match gpr {
    Gene(id) =>
      "<fbc:geneProductRef fbc:geneProduct=\"" + cobra3_xml_safe_id(id) + "\"/>"
    And(children) => {
      let sb = StringBuilder()
      let _ = sb.write("<fbc:and>")
      for c in children {
        let _ = sb.write(cobra3_xml_gpr(c))
      }
      let _ = sb.write("</fbc:and>")
      sb.to_string()
    }
    Or(children) => {
      let sb = StringBuilder()
      let _ = sb.write("<fbc:or>")
      for c in children {
        let _ = sb.write(cobra3_xml_gpr(c))
      }
      let _ = sb.write("</fbc:or>")
      sb.to_string()
    }
  }
}

///|
pub fn cobra3_model_to_sbml(model : Cobra3Model) -> String {
  let sb = StringBuilder()
  let _ = sb.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
  let _ = sb.write(
    "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version2/core\" " +
    "xmlns:fbc=\"http://www.sbml.org/sbml/level3/version1/fbc/version2\" " +
    "level=\"3\" version=\"2\" fbc:required=\"false\">\n",
  )
  let _ = sb.write(
    "  <model id=\"" +
    cobra3_xml_safe_id(model.id) +
    "\" name=\"" +
    cobra3_xml_escape(model.name) +
    "\">\n",
  )
  let _ = sb.write("    <listOfCompartments>\n")
  let comp_keys : Array[String] = Array(capacity=model.compartments.length())
  for k in model.compartments.keys() {
    comp_keys.push(k)
  }
  cobra3_sort_strings(comp_keys)
  for k in comp_keys {
    let _ = sb.write(
      "      <compartment id=\"" +
      cobra3_xml_safe_id(k) +
      "\" name=\"" +
      cobra3_xml_escape(model.compartments[k]) +
      "\" constant=\"true\"/>\n",
    )
  }
  let _ = sb.write("    </listOfCompartments>\n")
  let _ = sb.write("    <listOfSpecies>\n")
  for m in model.metabolites {
    let _ = sb.write(
      "      <species id=\"" +
      cobra3_xml_safe_id(m.id) +
      "\" name=\"" +
      cobra3_xml_escape(m.name) +
      "\" compartment=\"" +
      cobra3_xml_safe_id(m.compartment) +
      "\" boundaryCondition=\"" +
      (if m.boundary { "true" } else { "false" }) +
      "\" constant=\"false\"",
    )
    if m.formula.length() > 0 {
      let _ = sb.write(
        " fbc:chemicalFormula=\"" + cobra3_xml_escape(m.formula) + "\"",
      )
    }
    let _ = sb.write(" fbc:charge=\"" + m.charge.to_string() + "\"/>\n")
  }
  let _ = sb.write("    </listOfSpecies>\n")
  // gene products
  let genes : Array[String] = Array(capacity=0)
  for g in model.genes {
    genes.push(g.id)
  }
  for r in model.reactions {
    match r.gpr {
      Some(gpr) =>
        for g in cobra3_gpr_genes(gpr) {
          if !genes.contains(g) {
            genes.push(g)
          }
        }
      None => ()
    }
  }
  cobra3_sort_strings(genes)
  let _ = sb.write("    <fbc:listOfGeneProducts>\n")
  for g in genes {
    let _ = sb.write(
      "      <fbc:geneProduct fbc:id=\"" +
      cobra3_xml_safe_id(g) +
      "\" fbc:label=\"" +
      cobra3_xml_escape(g) +
      "\"/>\n",
    )
  }
  let _ = sb.write("    </fbc:listOfGeneProducts>\n")
  // parameters for bounds
  let _ = sb.write("    <listOfParameters>\n")
  for r in model.reactions {
    let safe = cobra3_xml_safe_id(r.id)
    let _ = sb.write(
      "      <parameter id=\"" +
      safe +
      "_lower\" value=\"" +
      r.lower_bound.to_string() +
      "\" constant=\"true\"/>\n",
    )
    let _ = sb.write(
      "      <parameter id=\"" +
      safe +
      "_upper\" value=\"" +
      r.upper_bound.to_string() +
      "\" constant=\"true\"/>\n",
    )
  }
  let _ = sb.write("    </listOfParameters>\n")
  let _ = sb.write("    <listOfReactions>\n")
  for r in model.reactions {
    let safe = cobra3_xml_safe_id(r.id)
    let reversible = r.lower_bound < 0.0
    let _ = sb.write(
      "      <reaction id=\"" +
      safe +
      "\" name=\"" +
      cobra3_xml_escape(r.name) +
      "\" reversible=\"" +
      (if reversible { "true" } else { "false" }) +
      "\" fbc:lowerFluxBound=\"" +
      safe +
      "_lower\" fbc:upperFluxBound=\"" +
      safe +
      "_upper\">\n",
    )
    let parts = cobra3_reaction_participants(r)
    let _ = sb.write("        <listOfReactants>\n")
    for p in parts {
      if p.1 < 0.0 {
        let _ = sb.write(
          "          <speciesReference species=\"" +
          cobra3_xml_safe_id(p.0) +
          "\" stoichiometry=\"" +
          (-p.1).to_string() +
          "\"/>\n",
        )
      }
    }
    let _ = sb.write("        </listOfReactants>\n")
    let _ = sb.write("        <listOfProducts>\n")
    for p in parts {
      if p.1 > 0.0 {
        let _ = sb.write(
          "          <speciesReference species=\"" +
          cobra3_xml_safe_id(p.0) +
          "\" stoichiometry=\"" +
          p.1.to_string() +
          "\"/>\n",
        )
      }
    }
    let _ = sb.write("        </listOfProducts>\n")
    match r.gpr {
      Some(gpr) => {
        let _ = sb.write("        <fbc:geneProductAssociation>\n          ")
        let _ = sb.write(cobra3_xml_gpr(gpr))
        let _ = sb.write("\n        </fbc:geneProductAssociation>\n")
      }
      None => ()
    }
    let _ = sb.write("      </reaction>\n")
  }
  let _ = sb.write("    </listOfReactions>\n")
  // objective
  let objectives : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    if r.objective.abs() > 1.0e-12 {
      objectives.push(r)
    }
  }
  if objectives.length() > 0 {
    let _ = sb.write("    <fbc:listOfObjectives>\n")
    let _ = sb.write(
      "      <fbc:objective fbc:id=\"obj1\" fbc:type=\"maximize\">\n",
    )
    let _ = sb.write("        <fbc:listOfFluxObjectives>\n")
    for r in objectives {
      let _ = sb.write(
        "          <fbc:fluxObjective fbc:reaction=\"" +
        cobra3_xml_safe_id(r.id) +
        "\" fbc:coefficient=\"" +
        r.objective.to_string() +
        "\"/>\n",
      )
    }
    let _ = sb.write("        </fbc:listOfFluxObjectives>\n")
    let _ = sb.write("      </fbc:objective>\n")
    let _ = sb.write("    </fbc:listOfObjectives>\n")
  }
  let _ = sb.write("  </model>\n")
  let _ = sb.write("</sbml>\n")
  sb.to_string()
}

///|
fn cobra3_collect_blocks(
  text : String,
  open_pattern : String,
  close_pattern : String,
) -> Array[String] {
  let out : Array[String] = Array(capacity=8)
  let mut pos = 0
  while pos < text.length() {
    match text[pos:].find(open_pattern) {
      None => break
      Some(i) => {
        let start = pos + i
        match text[start:].find(close_pattern) {
          None => break
          Some(j) => {
            out.push(text[start:start + j + close_pattern.length()].to_owned())
            pos = start + j + close_pattern.length()
          }
        }
      }
    }
  }
  out
}

///|
fn cobra3_gpr_from_xml(xml : String) -> Cobra3Gpr? {
  let stack : Array[Array[Cobra3Gpr]] = Array(capacity=8)
  let kinds : Array[Int32] = Array(capacity=8)
  stack.push(Array(capacity=4))
  kinds.push(0)
  let mut i = 0
  while i < xml.length() {
    let rest = xml[i:]
    if rest.has_prefix("<fbc:and>") {
      stack.push(Array(capacity=2))
      kinds.push(1)
      i = i + 9
      continue
    }
    if rest.has_prefix("<fbc:or>") {
      stack.push(Array(capacity=2))
      kinds.push(2)
      i = i + 8
      continue
    }
    if rest.has_prefix("</fbc:and>") || rest.has_prefix("</fbc:or>") {
      let is_and = rest.has_prefix("</fbc:and>")
      i = i + (if is_and { 10 } else { 9 })
      if stack.length() > 1 {
        let children = stack[stack.length() - 1]
        let _ = stack.pop()
        let kind = kinds[kinds.length() - 1]
        let _ = kinds.pop()
        let node = if kind == 1 {
          Cobra3Gpr::And(children)
        } else {
          Cobra3Gpr::Or(children)
        }
        stack[stack.length() - 1].push(node)
      }
      continue
    }
    if rest.has_prefix("<fbc:geneProductRef") {
      match rest.find("/>") {
        None => break
        Some(end) => {
          let tag = xml[i:i + end + 2].to_owned()
          match cobra3_xml_attr(tag, "fbc:geneProduct") {
            Some(g) => stack[stack.length() - 1].push(Cobra3Gpr::Gene(g))
            None => ()
          }
          i = i + end + 2
          continue
        }
      }
    }
    if rest.has_prefix("<fbc:geneProduct ") {
      match rest.find("/>") {
        None => break
        Some(end) => {
          let tag = xml[i:i + end + 2].to_owned()
          match cobra3_xml_attr(tag, "fbc:id") {
            Some(g) => stack[stack.length() - 1].push(Cobra3Gpr::Gene(g))
            None => ()
          }
          i = i + end + 2
          continue
        }
      }
    }
    i = i + 1
  }
  let root = stack[0]
  if root.length() == 0 {
    None
  } else if root.length() == 1 {
    Some(root[0])
  } else {
    Some(Cobra3Gpr::Or(root))
  }
}

///|
pub fn cobra3_read_sbml(text : String) -> Cobra3Model? {
  let mut model = cobra3_model_new("model")
  match cobra3_xml_attr(text, "model id") {
    Some(_) => ()
    None => ()
  }
  // model id/name：从 <model ...> 标签取
  match text.find("<model ") {
    Some(mi) =>
      match text[mi:].find(">") {
        Some(me) => {
          let tag = text[mi:mi + me].to_owned()
          let id = match cobra3_xml_attr(tag, "id") {
            Some(v) => v
            None => "model"
          }
          let name = match cobra3_xml_attr(tag, "name") {
            Some(v) => v
            None => id
          }
          model = {
            id,
            name,
            metabolites: model.metabolites,
            genes: model.genes,
            reactions: model.reactions,
            compartments: model.compartments,
            annotation: model.annotation,
            undo_stack: model.undo_stack,
          }
        }
        None => ()
      }
    None => ()
  }
  // compartments
  for block in cobra3_collect_blocks(text, "<compartment ", "/>") {
    match cobra3_xml_attr(block, "id") {
      Some(id) => {
        let name = match cobra3_xml_attr(block, "name") {
          Some(v) => v
          None => id
        }
        model.compartments[id] = name
      }
      None => ()
    }
  }
  // parameters
  let params : Map[String, Double] = Map([], capacity=16)
  for block in cobra3_collect_blocks(text, "<parameter ", "/>") {
    match cobra3_xml_attr(block, "id") {
      Some(id) =>
        match cobra3_xml_attr(block, "value") {
          Some(v) =>
            match parse_double(v) {
              Some(d) => params[id] = d
              None => ()
            }
          None => ()
        }
      None => ()
    }
  }
  // species
  for block in cobra3_collect_blocks(text, "<species ", "/>") {
    match cobra3_xml_attr(block, "id") {
      Some(id) => {
        let name = match cobra3_xml_attr(block, "name") {
          Some(v) => v
          None => id
        }
        let compartment = match cobra3_xml_attr(block, "compartment") {
          Some(v) => v
          None => "c"
        }
        let formula = match cobra3_xml_attr(block, "fbc:chemicalFormula") {
          Some(v) => v
          None => ""
        }
        let charge = match cobra3_xml_attr(block, "fbc:charge") {
          Some(v) => parse_double(v).unwrap_or(0.0)
          None => 0.0
        }
        let boundary = match cobra3_xml_attr(block, "boundaryCondition") {
          Some(v) => v == "true"
          None => false
        }
        let met : Cobra3Metabolite = {
          id,
          name,
          formula,
          charge,
          compartment,
          boundary,
          annotation: Map([], capacity=1),
        }
        model = cobra3_model_add_metabolite(model, met)
      }
      None => ()
    }
  }
  // objective 系数（model 级 listOfFluxObjectives）
  let objective_map : Map[String, Double] = Map([], capacity=4)
  for fo in cobra3_collect_blocks(text, "<fbc:fluxObjective ", "/>") {
    match cobra3_xml_attr(fo, "fbc:reaction") {
      Some(rid) => {
        let coeff = match cobra3_xml_attr(fo, "fbc:coefficient") {
          Some(v) => parse_double(v).unwrap_or(0.0)
          None => 0.0
        }
        objective_map[rid] = coeff
      }
      None => ()
    }
  }
  // reactions
  for block in cobra3_collect_blocks(text, "<reaction ", "</reaction>") {
    match cobra3_xml_attr(block, "id") {
      Some(id) => {
        let name = match cobra3_xml_attr(block, "name") {
          Some(v) => v
          None => id
        }
        let lb = match cobra3_xml_attr(block, "fbc:lowerFluxBound") {
          Some(v) =>
            match params.get(v) {
              Some(d) => d
              None => parse_double(v).unwrap_or(-1000.0)
            }
          None => -1000.0
        }
        let ub = match cobra3_xml_attr(block, "fbc:upperFluxBound") {
          Some(v) =>
            match params.get(v) {
              Some(d) => d
              None => parse_double(v).unwrap_or(1000.0)
            }
          None => 1000.0
        }
        // objective coefficient（来自 model 级 map）
        let objective = match objective_map.get(id) {
          Some(v) => v
          None => 0.0
        }
        let mets : Map[String, Double] = Map([], capacity=4)
        let react_start = match block.find("<listOfReactants>") {
          Some(x) => x
          None => 0
        }
        let react_end = match block.find("</listOfReactants>") {
          Some(x) => x
          None => react_start
        }
        let react_section = block[react_start:react_end].to_owned()
        for
          sr in cobra3_collect_blocks(react_section, "<speciesReference ", "/>") {
          match cobra3_xml_attr(sr, "species") {
            Some(sid) => {
              let coeff = match cobra3_xml_attr(sr, "stoichiometry") {
                Some(v) => parse_double(v).unwrap_or(1.0)
                None => 1.0
              }
              mets[sid] = -coeff
            }
            None => ()
          }
        }
        let prod_start = match block.find("<listOfProducts>") {
          Some(x) => x
          None => 0
        }
        let prod_end = match block.find("</listOfProducts>") {
          Some(x) => x
          None => prod_start
        }
        let prod_section = block[prod_start:prod_end].to_owned()
        for
          sr in cobra3_collect_blocks(prod_section, "<speciesReference ", "/>") {
          match cobra3_xml_attr(sr, "species") {
            Some(sid) => {
              let coeff = match cobra3_xml_attr(sr, "stoichiometry") {
                Some(v) => parse_double(v).unwrap_or(1.0)
                None => 1.0
              }
              let prev = match mets.get(sid) {
                Some(v) => v
                None => 0.0
              }
              mets[sid] = prev + coeff
            }
            None => ()
          }
        }
        // GPR
        let gpr = match block.find("<fbc:geneProductAssociation>") {
          Some(gs) =>
            match block.find("</fbc:geneProductAssociation>") {
              Some(ge) => {
                let section = block[gs:ge].to_owned()
                cobra3_gpr_from_xml(section)
              }
              None => None
            }
          None => None
        }
        let reaction : Cobra3Reaction = {
          id,
          name,
          metabolites: mets,
          lower_bound: lb,
          upper_bound: ub,
          objective,
          subsystem: "",
          gpr,
          annotation: Map([], capacity=1),
        }
        model = cobra3_model_add_reaction(model, reaction)
      }
      None => ()
    }
  }
  // gene products（补充未被反应引用的基因）
  for block in cobra3_collect_blocks(text, "<fbc:geneProduct ", "/>") {
    match cobra3_xml_attr(block, "fbc:id") {
      Some(id) => {
        let name = match cobra3_xml_attr(block, "fbc:label") {
          Some(v) => v
          None => id
        }
        if cobra3_model_gene_index(model, id) < 0 {
          let gene : Cobra3Gene = {
            id,
            name,
            annotation: Map([], capacity=1),
            functional: true,
          }
          model = cobra3_model_add_gene(model, gene)
        }
      }
      None => ()
    }
  }
  Some(model)
}
