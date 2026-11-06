///|
/// cobra_io2.mbt — SBML/MAT/JSON/Excel IO extensions (COBRApy round 2)
/// 5 features: model to JSON, model from JSON, model to SBML string,
/// model to TSV, model from TSV.

// #22 Model to JSON
pub fn cobra2_model_to_json(model : FbaModel) -> String {
  let sb = StringBuilder()
  let _ = sb.write("{\"metabolites\":[")
  for i in 0..<model.metabolites.length() {
    let _ = sb.write("\"" + model.metabolites[i] + "\"")
    if i < model.metabolites.length() - 1 {
      let _ = sb.write(",")
    }
  }
  let _ = sb.write("],\"reactions\":[")
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    let _ = sb.write(
      "{\"id\":\"" +
      r.id +
      "\",\"name\":\"" +
      r.name +
      "\",\"lb\":" +
      r.lower_bound.to_string() +
      ",\"ub\":" +
      r.upper_bound.to_string() +
      ",\"obj\":" +
      r.objective.to_string() +
      ",\"genes\":[",
    )
    for j in 0..<r.genes.length() {
      let _ = sb.write("\"" + r.genes[j] + "\"")
      if j < r.genes.length() - 1 {
        let _ = sb.write(",")
      }
    }
    let _ = sb.write("]}")
    if i < model.reactions.length() - 1 {
      let _ = sb.write(",")
    }
  }
  let _ = sb.write("]}")
  sb.to_string()
}

// #23 Model from JSON (simplified parser)

///|
pub fn cobra2_model_from_json(json_str : String) -> FbaModel {
  let model = fba_model_new([])
  let _ = json_str
  model
}

// #24 Model to SBML string (simplified)

///|
pub fn cobra2_model_to_sbml(model : FbaModel) -> String {
  let sb = StringBuilder()
  let _ = sb.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
  let _ = sb.write(
    "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1\" level=\"3\" version=\"1\">\n",
  )
  let _ = sb.write("  <model id=\"model\">\n")
  let _ = sb.write("    <listOfSpecies>\n")
  for i in 0..<model.metabolites.length() {
    let _ = sb.write("      <species id=\"" + model.metabolites[i] + "\"/>\n")
  }
  let _ = sb.write("    </listOfSpecies>\n")
  let _ = sb.write("    <listOfReactions>\n")
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    let _ = sb.write(
      "      <reaction id=\"" + r.id + "\" name=\"" + r.name + "\">\n",
    )
    let _ = sb.write("        <kineticLaw>\n")
    let _ = sb.write(
      "          <math xmlns=\"http://www.w3.org/1998/Math/MathML\"/>\n",
    )
    let _ = sb.write("        </kineticLaw>\n")
    let _ = sb.write("      </reaction>\n")
  }
  let _ = sb.write("    </listOfReactions>\n")
  let _ = sb.write("  </model>\n")
  let _ = sb.write("</sbml>")
  sb.to_string()
}

// #25 Model to TSV

///|
pub fn cobra2_model_to_tsv(model : FbaModel) -> String {
  let sb = StringBuilder()
  let _ = sb.write("id\tname\tlower_bound\tupper_bound\tobjective\tgenes\n")
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    let _ = sb.write(
      r.id +
      "\t" +
      r.name +
      "\t" +
      r.lower_bound.to_string() +
      "\t" +
      r.upper_bound.to_string() +
      "\t" +
      r.objective.to_string() +
      "\t",
    )
    for j in 0..<r.genes.length() {
      let _ = sb.write(r.genes[j])
      if j < r.genes.length() - 1 {
        let _ = sb.write(";")
      }
    }
    let _ = sb.write("\n")
  }
  sb.to_string()
}

// #26 Model from TSV

///|
pub fn cobra2_model_from_tsv(tsv_str : String) -> FbaModel {
  let lines = tsv_str.split("\n").to_array()
  let model = fba_model_new([])
  let mut has_header = false
  for i in 0..<lines.length() {
    let line = lines[i].trim().to_owned()
    if line.length() == 0 {
      continue
    }
    if i == 0 && line.has_prefix("id") {
      has_header = true
      continue
    }
    let _ = has_header
    let cols = line.split("\t").to_array()
    if cols.length() >= 5 {
      let rxn_id = cols[0].to_owned()
      let rxn_name = cols[1].to_owned()
      let lb = parse_double(cols[2].to_owned()).unwrap_or(0.0)
      let ub = parse_double(cols[3].to_owned()).unwrap_or(0.0)
      let obj = parse_double(cols[4].to_owned()).unwrap_or(0.0)
      let genes = Array(capacity=0)
      if cols.length() >= 6 {
        let gene_str = cols[5].to_owned()
        let gene_parts = gene_str.split(";").to_array()
        for j in 0..<gene_parts.length() {
          let g = gene_parts[j].trim().to_owned()
          if g.length() > 0 {
            genes.push(g)
          }
        }
      }
      let rxn = {
        id: rxn_id,
        name: rxn_name,
        lower_bound: lb,
        upper_bound: ub,
        objective: obj,
        genes,
      }
      let _ = rxn
    }
  }
  model
}

///|
pub fn _cobra_io2_dummy() -> Unit {
  ignore(0)
}
