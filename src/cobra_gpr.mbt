/// cobra_gpr.mbt — Gene-Protein-Reaction (GPR) boolean expression parsing
/// COBRApy core.GPR: parse, evaluate, and manipulate GPR rules.
/// GPR rules use boolean AND/OR: e.g. "gene1 and gene2 or gene3"

///|
/// GPR AST node types.
pub enum CobraGprNode {
  CobraGprGene(String)
  CobraGprAnd(Array[CobraGprNode])
  CobraGprOr(Array[CobraGprNode])
}

///|
/// GPR rule: parsed boolean expression.
pub struct CobraGpr {
  root : CobraGprNode?
  raw : String
}

///|
/// Parse a GPR string into an AST.
/// Supports: "gene1 and gene2", "gene1 or gene2", nested with parentheses.
pub fn cobra_gpr_parse(rule : String) -> CobraGpr {
  let trimmed = rule.trim().to_owned()
  if trimmed.length() == 0 || trimmed == "None" || trimmed == "" {
    return { root: None, raw: rule, }
  }
  let tokens = cobra_gpr_tokenize(trimmed)
  let (node, _) = cobra_gpr_parse_or(tokens, 0)
  { root: node, raw: rule, }
}

///|
/// Tokenize GPR string into tokens.
fn cobra_gpr_tokenize(s : String) -> Array[String] {
  let tokens : Array[String] = Array(capacity=0)
  let mut current = ""
  let mut i = 0
  let chars = s
  while i < chars.length() {
    let ch = chars[i]
    if ch == '(' || ch == ')' {
      if current.length() > 0 {
        tokens.push(current)
        current = ""
      }
      tokens.push(ch.to_string())
    } else if ch == ' ' || ch == '\t' {
      if current.length() > 0 {
        tokens.push(current)
        current = ""
      }
    } else {
      current = current + ch.to_string()
    }
    i = i + 1
  }
  if current.length() > 0 {
    tokens.push(current)
  }
  tokens
}

///|
/// Parse OR-level expression: term (or term)*
fn cobra_gpr_parse_or(
  tokens : Array[String],
  pos : Int32,
) -> (CobraGprNode?, Int32) {
  let (left, next) = cobra_gpr_parse_and(tokens, pos)
  if left is None {
    return (None, next)
  }
  let terms : Array[CobraGprNode] = Array(capacity=0)
  match left {
    Some(n) => terms.push(n)
    None => ()
  }
  let mut p = next
  while p < tokens.length() && tokens[p] == "or" {
    p = p + 1
    let (term, after) = cobra_gpr_parse_and(tokens, p)
    match term {
      Some(n) => terms.push(n)
      None => ()
    }
    p = after
  }
  if terms.length() == 1 {
    (Some(terms[0]), p)
  } else if terms.length() > 1 {
    (Some(CobraGprOr(terms)), p)
  } else {
    (None, p)
  }
}

///|
/// Parse AND-level expression: factor (and factor)*
fn cobra_gpr_parse_and(
  tokens : Array[String],
  pos : Int32,
) -> (CobraGprNode?, Int32) {
  let (left, next) = cobra_gpr_parse_factor(tokens, pos)
  if left is None {
    return (None, next)
  }
  let factors : Array[CobraGprNode] = Array(capacity=0)
  match left {
    Some(n) => factors.push(n)
    None => ()
  }
  let mut p = next
  while p < tokens.length() && tokens[p] == "and" {
    p = p + 1
    let (factor, after) = cobra_gpr_parse_factor(tokens, p)
    match factor {
      Some(n) => factors.push(n)
      None => ()
    }
    p = after
  }
  if factors.length() == 1 {
    (Some(factors[0]), p)
  } else if factors.length() > 1 {
    (Some(CobraGprAnd(factors)), p)
  } else {
    (None, p)
  }
}

///|
/// Parse factor: gene | (expr)
fn cobra_gpr_parse_factor(
  tokens : Array[String],
  pos : Int32,
) -> (CobraGprNode?, Int32) {
  if pos >= tokens.length() {
    return (None, pos)
  }
  let tok = tokens[pos]
  if tok == "(" {
    let (node, next) = cobra_gpr_parse_or(tokens, pos + 1)
    // Skip closing paren
    let after = if next < tokens.length() && tokens[next] == ")" {
      next + 1
    } else {
      next
    }
    return (node, after)
  }
  if tok == ")" {
    return (None, pos)
  }
  if tok == "and" || tok == "or" {
    return (None, pos)
  }
  (Some(CobraGprGene(tok)), pos + 1)
}

///|
/// Evaluate GPR with a set of knocked-out genes.
/// Returns true if the reaction is functional (all required genes active).
pub fn cobra_gpr_eval(gpr : CobraGpr, knockouts : Array[String]) -> Bool {
  match gpr.root {
    None => true // No GPR = always functional
    Some(node) => cobra_gpr_eval_node(node, knockouts)
  }
}

///|
fn cobra_gpr_eval_node(node : CobraGprNode, knockouts : Array[String]) -> Bool {
  match node {
    CobraGprGene(g) => !knockouts.contains(g)
    CobraGprAnd(children) => {
      let mut result = true
      for child in children {
        if !cobra_gpr_eval_node(child, knockouts) {
          result = false
          break
        }
      }
      result
    }
    CobraGprOr(children) => {
      let mut result = false
      for child in children {
        if cobra_gpr_eval_node(child, knockouts) {
          result = true
          break
        }
      }
      result
    }
  }
}

///|
/// Get all gene IDs mentioned in a GPR rule.
pub fn cobra_gpr_genes(gpr : CobraGpr) -> Array[String] {
  let genes : Array[String] = Array(capacity=0)
  match gpr.root {
    None => ()
    Some(node) => cobra_gpr_collect_genes(node, genes)
  }
  genes
}

///|
fn cobra_gpr_collect_genes(node : CobraGprNode, genes : Array[String]) -> Unit {
  match node {
    CobraGprGene(g) => if !genes.contains(g) { genes.push(g) }
    CobraGprAnd(children) =>
      for c in children {
        cobra_gpr_collect_genes(c, genes)
      }
    CobraGprOr(children) =>
      for c in children {
        cobra_gpr_collect_genes(c, genes)
      }
  }
}

///|
/// Convert GPR AST back to string.
pub fn cobra_gpr_to_string(gpr : CobraGpr) -> String {
  match gpr.root {
    None => ""
    Some(node) => cobra_gpr_node_to_string(node)
  }
}

///|
fn cobra_gpr_node_to_string(node : CobraGprNode) -> String {
  match node {
    CobraGprGene(g) => g
    CobraGprAnd(children) => {
      let parts : Array[String] = Array(capacity=0)
      for c in children {
        parts.push(cobra_gpr_node_to_string(c))
      }
      parts.join(" and ")
    }
    CobraGprOr(children) => {
      let parts : Array[String] = Array(capacity=0)
      for c in children {
        parts.push("(" + cobra_gpr_node_to_string(c) + ")")
      }
      parts.join(" or ")
    }
  }
}
