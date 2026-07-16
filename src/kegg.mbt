///|
/// KEGG database parsing, inspired by Biopython Bio.KEGG.

///|
pub struct KEGGGene {
  entry : String
  name : String
  definition : String
  orthology : String
  organism : String
  pathways : Array[String]
  position : String
  dblinks : Array[(String, String)]
}

///|
pub fn KEGGGene::new(
  entry : String,
  name : String,
  definition : String,
  orthology : String,
  organism : String,
  pathways : Array[String],
  position : String,
  dblinks : Array[(String, String)],
) -> KEGGGene {
  KEGGGene::{
    entry,
    name,
    definition,
    orthology,
    organism,
    pathways,
    position,
    dblinks,
  }
}

///|
pub struct KEGGPathway {
  entry : String
  name : String
  definition : String
  class : String
  genes : Array[(String, String)]
  compounds : Array[String]
  enzymes : Array[String]
  organisms : Array[String]
}

///|
pub fn KEGGPathway::new(
  entry : String,
  name : String,
  definition : String,
  class : String,
  genes : Array[(String, String)],
  compounds : Array[String],
  enzymes : Array[String],
  organisms : Array[String],
) -> KEGGPathway {
  KEGGPathway::{
    entry,
    name,
    definition,
    class,
    genes,
    compounds,
    enzymes,
    organisms,
  }
}

///|
pub struct KEGGCompound {
  entry : String
  name : String
  formula : String
  mass : Double
  pathways : Array[String]
  enzymes : Array[String]
  reactions : Array[String]
}

///|
pub fn KEGGCompound::new(
  entry : String,
  name : String,
  formula : String,
  mass : Double,
  pathways : Array[String],
  enzymes : Array[String],
  reactions : Array[String],
) -> KEGGCompound {
  KEGGCompound::{ entry, name, formula, mass, pathways, enzymes, reactions }
}

///|
pub struct KEGGEnzyme {
  entry : String
  name : String
  class : String
  genes : Array[(String, String)]
  pathways : Array[String]
  reactions : Array[String]
  cofactors : Array[String]
}

///|
pub fn KEGGEnzyme::new(
  entry : String,
  name : String,
  class : String,
  genes : Array[(String, String)],
  pathways : Array[String],
  reactions : Array[String],
  cofactors : Array[String],
) -> KEGGEnzyme {
  KEGGEnzyme::{ entry, name, class, genes, pathways, reactions, cofactors }
}

///|
pub fn parse_kegg_gene(text : String) -> KEGGGene? {
  let lines = split_by_char(text, '\n'.to_int())

  let mut entry = ""
  let mut name = ""
  let mut definition = ""
  let mut orthology = ""
  let mut organism = ""
  let pathways : Array[String] = Array::new()
  let mut position = ""
  let dblinks : Array[(String, String)] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_string()

    if starts_with(line, 0, "ENTRY") {
      entry = line[6:].trim().to_owned()
    } else if starts_with(line, 0, "NAME") {
      name = line[5:].trim().to_owned()
    } else if starts_with(line, 0, "DEFINITION") {
      definition = line[11:].trim().to_owned()
    } else if starts_with(line, 0, "ORTHOLOGY") {
      orthology = line[10:].trim().to_owned()
    } else if starts_with(line, 0, "ORGANISM") {
      organism = line[9:].trim().to_owned()
    } else if starts_with(line, 0, "PATHWAY") {
      let p = line[8:].trim().to_owned()
      pathways.push(p)
    } else if starts_with(line, 0, "POSITION") {
      position = line[9:].trim().to_owned()
    } else if starts_with(line, 0, "DBLINKS") {
      let parts = split_by_char(line[8:].trim().to_owned(), ':'.to_int())
      if parts.length() >= 2 {
        dblinks.push((parts[0], parts[1]))
      }
    }

    i = i + 1
  }

  if entry == "" {
    None
  } else {
    Some(
      KEGGGene::new(
        entry, name, definition, orthology, organism, pathways, position, dblinks,
      ),
    )
  }
}

///|
pub fn parse_kegg_pathway(text : String) -> KEGGPathway? {
  let lines = split_by_char(text, '\n'.to_int())

  let mut entry = ""
  let mut name = ""
  let mut definition = ""
  let mut class = ""
  let genes : Array[(String, String)] = Array::new()
  let compounds : Array[String] = Array::new()
  let enzymes : Array[String] = Array::new()
  let organisms : Array[String] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_string()

    if starts_with(line, 0, "ENTRY") {
      entry = line[6:].trim().to_owned()
    } else if starts_with(line, 0, "NAME") {
      name = line[5:].trim().to_owned()
    } else if starts_with(line, 0, "DEFINITION") {
      definition = line[11:].trim().to_owned()
    } else if starts_with(line, 0, "CLASS") {
      class = line[6:].trim().to_owned()
    } else if starts_with(line, 0, "GENE") {
      let parts = split_by_char(line[5:].trim().to_owned(), ' '.to_int())
      let mut j = 0
      while j < parts.length() {
        if parts[j] != "" {
          let mut k = j + 1
          while k < parts.length() && parts[k] == "" {
            k = k + 1
          }
          if k < parts.length() {
            genes.push((parts[j], parts[k]))
          }
          break
        }
        j = j + 1
      }
    } else if starts_with(line, 0, "COMPOUND") {
      let c = line[9:].trim().to_owned()
      compounds.push(c)
    } else if starts_with(line, 0, "ENZYME") {
      let e = line[7:].trim().to_owned()
      enzymes.push(e)
    } else if starts_with(line, 0, "ORGANISM") {
      let o = line[9:].trim().to_owned()
      organisms.push(o)
    }

    i = i + 1
  }

  if entry == "" {
    None
  } else {
    Some(
      KEGGPathway::new(
        entry, name, definition, class, genes, compounds, enzymes, organisms,
      ),
    )
  }
}

///|
pub fn parse_kegg_compound(text : String) -> KEGGCompound? {
  let lines = split_by_char(text, '\n'.to_int())

  let mut entry = ""
  let mut name = ""
  let mut formula = ""
  let mut mass = 0.0
  let pathways : Array[String] = Array::new()
  let enzymes : Array[String] = Array::new()
  let reactions : Array[String] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_string()

    if starts_with(line, 0, "ENTRY") {
      entry = line[6:].trim().to_owned()
    } else if starts_with(line, 0, "NAME") {
      name = line[5:].trim().to_owned()
    } else if starts_with(line, 0, "FORMULA") {
      formula = line[8:].trim().to_owned()
    } else if starts_with(line, 0, "MASS") {
      let m_str = line[5:].trim().to_owned()
      mass = parse_double(m_str).unwrap_or(0.0)
    } else if starts_with(line, 0, "PATHWAY") {
      let p = line[8:].trim().to_owned()
      pathways.push(p)
    } else if starts_with(line, 0, "ENZYME") {
      let e = line[7:].trim().to_owned()
      enzymes.push(e)
    } else if starts_with(line, 0, "REACTION") {
      let r = line[9:].trim().to_owned()
      reactions.push(r)
    }

    i = i + 1
  }

  if entry == "" {
    None
  } else {
    Some(
      KEGGCompound::new(
        entry, name, formula, mass, pathways, enzymes, reactions,
      ),
    )
  }
}

///|
pub fn parse_kegg_enzyme(text : String) -> KEGGEnzyme? {
  let lines = split_by_char(text, '\n'.to_int())

  let mut entry = ""
  let mut name = ""
  let mut class = ""
  let genes : Array[(String, String)] = Array::new()
  let pathways : Array[String] = Array::new()
  let reactions : Array[String] = Array::new()
  let cofactors : Array[String] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_string()

    if starts_with(line, 0, "ENTRY") {
      entry = line[6:].trim().to_owned()
    } else if starts_with(line, 0, "NAME") {
      name = line[5:].trim().to_owned()
    } else if starts_with(line, 0, "CLASS") {
      class = line[6:].trim().to_owned()
    } else if starts_with(line, 0, "GENE") {
      let parts = split_by_char(line[5:].trim().to_owned(), ' '.to_int())
      let mut j = 0
      while j < parts.length() {
        if parts[j] != "" {
          let mut k = j + 1
          while k < parts.length() && parts[k] == "" {
            k = k + 1
          }
          if k < parts.length() {
            genes.push((parts[j], parts[k]))
          }
          break
        }
        j = j + 1
      }
    } else if starts_with(line, 0, "PATHWAY") {
      let p = line[8:].trim().to_owned()
      pathways.push(p)
    } else if starts_with(line, 0, "REACTION") {
      let r = line[9:].trim().to_owned()
      reactions.push(r)
    } else if starts_with(line, 0, "COFACTOR") {
      let c = line[9:].trim().to_owned()
      cofactors.push(c)
    }

    i = i + 1
  }

  if entry == "" {
    None
  } else {
    Some(
      KEGGEnzyme::new(entry, name, class, genes, pathways, reactions, cofactors),
    )
  }
}

///|
pub fn kegg_pathway_gene_ids(pathway : KEGGPathway) -> Array[String] {
  let ids : Array[String] = Array::new()
  let mut i = 0
  while i < pathway.genes.length() {
    ids.push(pathway.genes[i].0)
    i = i + 1
  }
  ids
}

///|
pub fn kegg_gene_pathway_count(gene : KEGGGene) -> Int {
  gene.pathways.length()
}

///|
pub fn kegg_compound_in_pathway(
  compound : KEGGCompound,
  pathway_id : String,
) -> Bool {
  let len = compound.pathways.length()
  let id_len = pathway_id.length()
  let mut i = 0
  while i < len {
    let path = compound.pathways[i]
    let path_len = path.length()
    if path_len >= id_len {
      let sub = path[0:id_len]
      if sub.to_string() == pathway_id {
        return true
      }
    }
    i = i + 1
  }
  false
}
