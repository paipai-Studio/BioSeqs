/// cobra_sbml_io.mbt — SBML format read/write (COBRApy io.sbml)
/// Minimal SBML Level 3 Version 1 parser/writer for metabolic models.

///|
/// Write FBA model to SBML L3V1 XML string.
pub fn cobra_write_sbml(model : FbaModel) -> String {
  let mut sb = ""
  // XML header
  sb = sb + "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  sb = sb + "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\" "
  sb = sb + "level=\"3\" version=\"1\">\n"
  // Model
  sb = sb + "  <model id=\"FbaModel\">\n"
  // Compartments (single default compartment)
  sb = sb + "    <listOfCompartments>\n"
  sb = sb +
    "      <compartment id=\"c\" name=\"cytosol\" size=\"1\" constant=\"true\"/>\n"
  sb = sb + "    </listOfCompartments>\n"
  // Species (metabolites)
  sb = sb + "    <listOfSpecies>\n"
  for met in model.metabolites {
    sb = sb + "      <species id=\"" + met + "\" compartment=\"c\" "
    sb = sb + "boundaryCondition=\"false\" constant=\"false\"/>\n"
  }
  sb = sb + "    </listOfSpecies>\n"
  // Reactions
  sb = sb + "    <listOfReactions>\n"
  for j = 0; j < model.reactions.length(); j = j + 1 {
    let r = model.reactions[j]
    let reversible = r.lower_bound < 0.0
    sb = sb + "      <reaction id=\"" + r.id + "\" name=\"" + r.name + "\" "
    sb = sb +
      "reversible=\"" +
      (if reversible { "true" } else { "false" }) +
      "\">\n"
    // Reactants and products
    sb = sb + "        <listOfReactants>\n"
    for i = 0; i < model.metabolites.length(); i = i + 1 {
      let s = model.stoichiometry[i][j]
      if s < -1.0e-9 {
        sb = sb +
          "          <speciesReference species=\"" +
          model.metabolites[i] +
          "\" "
        sb = sb + "stoichiometry=\"" + (-s).to_string() + "\"/>\n"
      }
    }
    sb = sb + "        </listOfReactants>\n"
    sb = sb + "        <listOfProducts>\n"
    for i = 0; i < model.metabolites.length(); i = i + 1 {
      let s = model.stoichiometry[i][j]
      if s > 1.0e-9 {
        sb = sb +
          "          <speciesReference species=\"" +
          model.metabolites[i] +
          "\" "
        sb = sb + "stoichiometry=\"" + s.to_string() + "\"/>\n"
      }
    }
    sb = sb + "        </listOfProducts>\n"
    // Bounds via fbc extension (simplified as notes)
    sb = sb + "        <notes>\n"
    sb = sb +
      "          <p>lower_bound: " +
      r.lower_bound.to_string() +
      "</p>\n"
    sb = sb +
      "          <p>upper_bound: " +
      r.upper_bound.to_string() +
      "</p>\n"
    sb = sb + "          <p>objective: " + r.objective.to_string() + "</p>\n"
    sb = sb + "        </notes>\n"
    sb = sb + "      </reaction>\n"
  }
  sb = sb + "    </listOfReactions>\n"
  sb = sb + "  </model>\n"
  sb = sb + "</sbml>\n"
  sb
}

///|
/// Parse SBML XML string to FBA model；优先使用 FBC v2 解析，失败时回退旧版简化解析.
pub fn cobra_read_sbml(text : String) -> FbaModel? {
  match cobra3_read_sbml(text) {
    Some(m) => Some(cobra3_to_fba_model(m))
    None => cobra_read_sbml_legacy(text)
  }
}

///|
/// 旧版简化 SBML 解析（仅作为 FBC 解析失败时的回退）.
fn cobra_read_sbml_legacy(text : String) -> FbaModel? {
  // Extract species
  let metabolites : Array[String] = Array(capacity=0)
  let species_pattern = "species id=\""
  let mut pos = 0
  while pos < text.length() {
    let idx = text[pos:].find(species_pattern)
    match idx {
      None => break
      Some(i) => {
        let start = pos + i + species_pattern.length()
        let end = text[start:].find("\"")
        match end {
          None => break
          Some(e) => {
            metabolites.push(text[start:start + e].to_owned())
            pos = start + e + 1
          }
        }
      }
    }
  }
  // Extract reactions
  let reactions : Array[FbaReaction] = Array(capacity=0)
  let rxn_pattern = "reaction id=\""
  pos = 0
  while pos < text.length() {
    let idx = text[pos:].find(rxn_pattern)
    match idx {
      None => break
      Some(i) => {
        let start = pos + i + rxn_pattern.length()
        let end = text[start:].find("\"")
        match end {
          None => break
          Some(e) => {
            let rxn_id = text[start:start + e].to_owned()
            // Find the end of this reaction element
            let rxn_end = text[start:].find("</reaction>")
            let rxn_section = match rxn_end {
              None => text[start:].to_owned()
              Some(re) => text[start:start + re].to_owned()
            }
            // Parse bounds from notes
            let lb = cobra_extract_note(rxn_section, "lower_bound: ")
            let ub = cobra_extract_note(rxn_section, "upper_bound: ")
            let obj = cobra_extract_note(rxn_section, "objective: ")
            reactions.push({
              id: rxn_id,
              name: rxn_id,
              lower_bound: lb,
              upper_bound: ub,
              objective: obj,
              genes: [],
            })
            pos = start + e + 1
          }
        }
      }
    }
  }
  // Build stoichiometry from speciesReference entries
  let m = metabolites.length()
  let n_r = reactions.length()
  let stoich : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    stoich.push(Array::make(n_r, 0.0))
  }
  // Parse reactants and products for each reaction
  for j = 0; j < n_r; j = j + 1 {
    let rxn_id = reactions[j].id
    let rxn_start_pattern = "reaction id=\"" + rxn_id + "\""
    let ridx = text.find(rxn_start_pattern)
    match ridx {
      None => continue
      Some(rs) => {
        let rend = text[rs:].find("</reaction>")
        let section = match rend {
          None => text[rs:].to_owned()
          Some(re) => text[rs:rs + re].to_owned()
        }
        // Parse reactants (negative stoichiometry)
        let reactant_section = cobra_extract_section(section, "listOfReactants")
        for i = 0; i < m; i = i + 1 {
          let met = metabolites[i]
          let ref_pattern = "species=\"" + met + "\" stoichiometry=\""
          let ref_idx = reactant_section.find(ref_pattern)
          match ref_idx {
            None => ()
            Some(ri) => {
              let num_start = ri + ref_pattern.length()
              let num_end = reactant_section[num_start:].find("\"")
              match num_end {
                None => ()
                Some(ne) => {
                  let val = parse_double(
                    reactant_section[num_start:num_start + ne].to_owned(),
                  ).unwrap_or(1.0)
                  stoich[i][j] = stoich[i][j] - val
                }
              }
            }
          }
        }
        // Parse products (positive stoichiometry)
        let product_section = cobra_extract_section(section, "listOfProducts")
        for i = 0; i < m; i = i + 1 {
          let met = metabolites[i]
          let ref_pattern = "species=\"" + met + "\" stoichiometry=\""
          let ref_idx = product_section.find(ref_pattern)
          match ref_idx {
            None => ()
            Some(pi) => {
              let num_start = pi + ref_pattern.length()
              let num_end = product_section[num_start:].find("\"")
              match num_end {
                None => ()
                Some(ne) => {
                  let val = parse_double(
                    product_section[num_start:num_start + ne].to_owned(),
                  ).unwrap_or(1.0)
                  stoich[i][j] = stoich[i][j] + val
                }
              }
            }
          }
        }
      }
    }
  }
  if metabolites.length() == 0 && reactions.length() == 0 {
    return None
  }
  Some({ metabolites, reactions, stoichiometry: stoich, })
}

///|
/// Extract a numeric value from a notes section.
fn cobra_extract_note(text : String, key : String) -> Double {
  let idx = text.find(key)
  match idx {
    None => 0.0
    Some(i) => {
      let start = i + key.length()
      let end = text[start:].find("<")
      match end {
        None => parse_double(text[start:].to_owned()).unwrap_or(0.0)
        Some(e) => parse_double(text[start:start + e].to_owned()).unwrap_or(0.0)
      }
    }
  }
}

///|
/// Extract content between <listOfX>...</listOfX>
fn cobra_extract_section(text : String, tag : String) -> String {
  let open_tag = "<" + tag + ">"
  let close_tag = "</" + tag + ">"
  let start = text.find(open_tag)
  match start {
    None => ""
    Some(s) => {
      let content_start = s + open_tag.length()
      let end = text[content_start:].find(close_tag)
      match end {
        None => text[content_start:].to_owned()
        Some(e) => text[content_start:content_start + e].to_owned()
      }
    }
  }
}
