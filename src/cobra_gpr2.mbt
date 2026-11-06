///|
/// cobra_gpr2.mbt — Advanced GPR parsing and evaluation (COBRApy round 2)
/// 5 features: GPR parse tree, evaluate with gene states, GPR to string,
/// extract genes, simplify GPR.

// #16 GPR parse tree
pub enum Cobra2GprNode {
  Gene(String)
  And(Array[Cobra2GprNode])
  Or(Array[Cobra2GprNode])
}

// #17 Evaluate GPR with gene activity states

///|
pub fn cobra2_eval_gpr(
  node : Cobra2GprNode,
  active_genes : Array[String],
) -> Bool {
  match node {
    Cobra2GprNode::Gene(g) => {
      let mut found = false
      for i in 0..<active_genes.length() {
        if active_genes[i] == g {
          found = true
          break
        }
      }
      found
    }
    Cobra2GprNode::And(children) => {
      let mut all_true = true
      for i in 0..<children.length() {
        if not(cobra2_eval_gpr(children[i], active_genes)) {
          all_true = false
          break
        }
      }
      all_true
    }
    Cobra2GprNode::Or(children) => {
      let mut any_true = false
      for i in 0..<children.length() {
        if cobra2_eval_gpr(children[i], active_genes) {
          any_true = true
          break
        }
      }
      any_true
    }
  }
}

// #18 GPR to string

///|
pub fn cobra2_gpr_to_string(node : Cobra2GprNode) -> String {
  match node {
    Cobra2GprNode::Gene(g) => g
    Cobra2GprNode::And(children) => {
      let parts = Array(capacity=children.length())
      for i in 0..<children.length() {
        parts.push(cobra2_gpr_to_string(children[i]))
      }
      let sb = StringBuilder()
      let _ = sb.write("(")
      for i in 0..<parts.length() {
        let _ = sb.write(parts[i])
        if i < parts.length() - 1 {
          let _ = sb.write(" and ")
        }
      }
      let _ = sb.write(")")
      sb.to_string()
    }
    Cobra2GprNode::Or(children) => {
      let parts = Array(capacity=children.length())
      for i in 0..<children.length() {
        parts.push(cobra2_gpr_to_string(children[i]))
      }
      let sb = StringBuilder()
      let _ = sb.write("(")
      for i in 0..<parts.length() {
        let _ = sb.write(parts[i])
        if i < parts.length() - 1 {
          let _ = sb.write(" or ")
        }
      }
      let _ = sb.write(")")
      sb.to_string()
    }
  }
}

// #19 Extract all genes from GPR tree

///|
pub fn cobra2_extract_genes(node : Cobra2GprNode) -> Array[String] {
  let genes = Array(capacity=0)
  cobra2_extract_genes_helper(node, genes)
  genes
}

///|
fn cobra2_extract_genes_helper(
  node : Cobra2GprNode,
  genes : Array[String],
) -> Unit {
  match node {
    Cobra2GprNode::Gene(g) => {
      let mut found = false
      for i in 0..<genes.length() {
        if genes[i] == g {
          found = true
          break
        }
      }
      if not(found) {
        genes.push(g)
      }
    }
    Cobra2GprNode::And(children) =>
      for i in 0..<children.length() {
        cobra2_extract_genes_helper(children[i], genes)
      }
    Cobra2GprNode::Or(children) =>
      for i in 0..<children.length() {
        cobra2_extract_genes_helper(children[i], genes)
      }
  }
}

// #20 Simplify GPR — remove duplicate genes, flatten nested same-type nodes

///|
pub fn cobra2_simplify_gpr(node : Cobra2GprNode) -> Cobra2GprNode {
  match node {
    Gene(_) => node
    Cobra2GprNode::And(children) => {
      let simplified = Array(capacity=0)
      let seen : Map[String, Bool] = Map([])
      for i in 0..<children.length() {
        let child = cobra2_simplify_gpr(children[i])
        let key = cobra2_gpr_to_string(child)
        if not(seen.contains(key)) {
          seen[key] = true
          simplified.push(child)
        }
      }
      if simplified.length() == 1 {
        simplified[0]
      } else {
        Cobra2GprNode::And(simplified)
      }
    }
    Cobra2GprNode::Or(children) => {
      let simplified = Array(capacity=0)
      let seen : Map[String, Bool] = Map([])
      for i in 0..<children.length() {
        let child = cobra2_simplify_gpr(children[i])
        let key = cobra2_gpr_to_string(child)
        if not(seen.contains(key)) {
          seen[key] = true
          simplified.push(child)
        }
      }
      if simplified.length() == 1 {
        simplified[0]
      } else {
        Cobra2GprNode::Or(simplified)
      }
    }
  }
}

// #21 Parse simple GPR string (e.g. "geneA and geneB or geneC")

///|
pub fn cobra2_parse_gpr_simple(gpr_str : String) -> Cobra2GprNode {
  let trimmed = gpr_str.trim().to_owned()
  if trimmed.length() == 0 {
    return Cobra2GprNode::Gene("")
  }
  let or_parts = trimmed.split(" or ").to_array()
  if or_parts.length() > 1 {
    let children = Array(capacity=or_parts.length())
    for i in 0..<or_parts.length() {
      children.push(cobra2_parse_and_part(or_parts[i].to_owned()))
    }
    return Cobra2GprNode::Or(children)
  }
  cobra2_parse_and_part(trimmed)
}

///|
fn cobra2_parse_and_part(s : String) -> Cobra2GprNode {
  let and_parts = s.split(" and ").to_array()
  if and_parts.length() > 1 {
    let children = Array(capacity=and_parts.length())
    for i in 0..<and_parts.length() {
      let g = and_parts[i].trim().to_owned()
      children.push(Cobra2GprNode::Gene(g))
    }
    return Cobra2GprNode::And(children)
  }
  Gene(s.trim().to_owned())
}

///|
pub fn _cobra_gpr2_dummy() -> Unit {
  let _ = Cobra2GprNode::Gene("")
  let _ = Cobra2GprNode::And([])
  let _ = Cobra2GprNode::Or([])
}
