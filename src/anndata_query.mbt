/// anndata_query.mbt — pandas-style query expressions over obs/var.
///
/// Mirrors `adata.obs.query("n_genes > 500 & batch == 'b1'")` from the
/// anndata/pandas ecosystem: a small recursive-descent parser evaluates
/// comparison predicates combined with `&` (and), `|` (or) and `!` (not)
/// over annotation columns, returning matching row/column indices.

///|
/// Query comparison operators.
/// Tokens produced by the query lexer.
pub(all) enum QTok {
  QIdent(String) // column name
  QStr(String) // 'value' or "value"
  QNum(Double) // numeric literal
  QOp(String) // == != >= <= > < & | !
  QLP // (
  QRP // )
  QEnd
}

///|
/// Comparison operand: string literal or number.
pub(all) enum QVal {
  VStr(String)
  VNum(Double)
}

///|
/// Parsed query predicate tree.
pub(all) enum QPred {
  QCmp(String, String, QVal) // column, op, value
  QAnd(QPred, QPred)
  QOr(QPred, QPred)
  QNot(QPred)
}

///|
/// Tokenize a query expression. Returns `None` on stray characters.
pub fn anndata_query_tokenize(expr : String) -> Array[QTok]? {
  let toks : Array[QTok] = Array(capacity=0)
  let chars = expr.to_array()
  let n = chars.length()
  let mut i = 0
  while i < n {
    let c = chars[i]
    let code = c.to_int()
    if code == 32 || code == 9 {
      i = i + 1
    } else if code == 40 {
      // (
      toks.push(QLP)
      i = i + 1
    } else if code == 41 {
      // )
      toks.push(QRP)
      i = i + 1
    } else if code == 39 || code == 34 {
      // ' or "
      let quote = code
      let mut j = i + 1
      let sb = StringBuilder()
      while j < n && chars[j].to_int() != quote {
        let _ = sb.write_char(chars[j])
        j = j + 1
      }
      if j >= n {
        return None
      }
      toks.push(QStr(sb.to_string()))
      i = j + 1
    } else if (code >= 48 && code <= 57) || code == 45 || code == 46 {
      // digit, '-' or '.': number literal
      let sb = StringBuilder()
      let mut j = i
      while j < n {
        let cj = chars[j].to_int()
        if (cj >= 48 && cj <= 57) ||
          cj == 46 ||
          cj == 45 ||
          cj == 101 ||
          cj == 69 {
          // digits, '.', '-', 'e', 'E'
          let _ = sb.write_char(chars[j])
          j = j + 1
        } else {
          break
        }
      }
      match parse_double(sb.to_string()) {
        Some(d) => {
          toks.push(QNum(d))
          i = j
        }
        None => return None
      }
    } else if (code >= 65 && code <= 90) ||
      (code >= 97 && code <= 122) ||
      code == 95 {
      // letter or '_': identifier
      let sb = StringBuilder()
      let mut j = i
      while j < n {
        let cj = chars[j].to_int()
        if (cj >= 65 && cj <= 90) ||
          (cj >= 97 && cj <= 122) ||
          (cj >= 48 && cj <= 57) ||
          cj == 95 {
          let _ = sb.write_char(chars[j])
          j = j + 1
        } else {
          break
        }
      }
      toks.push(QIdent(sb.to_string()))
      i = j
    } else if code == 38 || code == 124 || code == 33 {
      // & | !
      toks.push(QOp(c.to_string()))
      i = i + 1
    } else if code == 61 || code == 60 || code == 62 {
      // = < : maybe two-char op
      if i + 1 < n {
        let c2 = chars[i + 1].to_int()
        if (code == 61 && c2 == 61) ||
          (code == 60 && c2 == 61) ||
          (code == 62 && c2 == 61) {
          toks.push(QOp(c.to_string() + chars[i + 1].to_string()))
          i = i + 2
        } else {
          toks.push(QOp(c.to_string()))
          i = i + 1
        }
      } else {
        toks.push(QOp(c.to_string()))
        i = i + 1
      }
    } else {
      return None
    }
  }
  toks.push(QEnd)
  Some(toks)
}

///|
/// Recursive-descent parser state.
priv struct QParser {
  toks : Array[QTok]
  mut pos : Int32
}

///|
fn QParser::peek(self : QParser) -> QTok {
  self.toks[self.pos]
}

///|
fn QParser::next_tok(self : QParser) -> QTok {
  let t = self.toks[self.pos]
  if self.pos < self.toks.length() - 1 {
    self.pos = self.pos + 1
  }
  t
}

///|
/// or_expr := and_expr ('|' and_expr)*
fn QParser::parse_or(self : QParser) -> QPred? {
  let left = match self.parse_and() {
    Some(p) => p
    None => return None
  }
  let mut acc = left
  while true {
    match self.peek() {
      QOp(op) =>
        if op == "|" {
          ignore(self.next_tok())
          match self.parse_and() {
            Some(right) => acc = QOr(acc, right)
            None => return None
          }
        } else {
          break
        }
      _ => break
    }
  }
  Some(acc)
}

///|
/// and_expr := unary ('&' unary)*
fn QParser::parse_and(self : QParser) -> QPred? {
  let mut acc = match self.parse_unary() {
    Some(p) => p
    None => return None
  }
  while true {
    match self.peek() {
      QOp(op) =>
        if op == "&" {
          ignore(self.next_tok())
          match self.parse_unary() {
            Some(right) => acc = QAnd(acc, right)
            None => return None
          }
        } else {
          break
        }
      _ => break
    }
  }
  Some(acc)
}

///|
/// unary := ('!' unary) | atom
fn QParser::parse_unary(self : QParser) -> QPred? {
  match self.peek() {
    QOp(op) =>
      if op == "!" {
        ignore(self.next_tok())
        match self.parse_unary() {
          Some(inner) => Some(QNot(inner))
          None => None
        }
      } else {
        None
      }
    _ => self.parse_atom()
  }
}

///|
/// atom := '(' or_expr ')' | comparison
fn QParser::parse_atom(self : QParser) -> QPred? {
  match self.peek() {
    QLP => {
      ignore(self.next_tok())
      match self.parse_or() {
        Some(inner) =>
          match self.peek() {
            QRP => {
              ignore(self.next_tok())
              Some(inner)
            }
            _ => None
          }
        None => None
      }
    }
    QIdent(col) => {
      ignore(self.next_tok())
      match self.peek() {
        QOp(op) => {
          ignore(self.next_tok())
          let value : QVal = match self.peek() {
            QStr(s) => {
              ignore(self.next_tok())
              VStr(s)
            }
            QNum(d) => {
              ignore(self.next_tok())
              VNum(d)
            }
            _ => return None
          }
          Some(QCmp(col, op, value))
        }
        _ => None
      }
    }
    _ => None
  }
}

///|
/// Parse a query expression into a predicate tree.
pub fn anndata_query_parse(expr : String) -> QPred? {
  match anndata_query_tokenize(expr) {
    Some(toks) => {
      let p = { toks, pos: 0, }
      match p.parse_or() {
        Some(pred) => if p.peek() is QEnd { Some(pred) } else { None }
        None => None
      }
    }
    None => None
  }
}

///|
/// Fetch an obs value as string ("" when the column is missing).
fn query_obs_value(adata : AnnData, col : String, i : Int32) -> String {
  match adata.obs.get(col) {
    Some(arr) => if i < arr.length() { arr[i] } else { "" }
    None => ""
  }
}

///|
/// Fetch a var value as string.
fn query_var_value(adata : AnnData, col : String, j : Int32) -> String {
  match adata.variables.get(col) {
    Some(arr) => if j < arr.length() { arr[j] } else { "" }
    None => ""
  }
}

///|
/// Compare a string cell against the predicate operand.
fn query_cmp_str(cell : String, op : String, value : QVal) -> Bool {
  match value {
    VStr(s) =>
      if op == "==" {
        cell == s
      } else if op == "!=" {
        cell != s
      } else if op == ">" {
        cell > s
      } else if op == ">=" {
        cell >= s
      } else if op == "<" {
        cell < s
      } else if op == "<=" {
        cell <= s
      } else {
        false
      }
    VNum(d) =>
      // numeric comparison on a string cell: parse or fail
      match parse_double(cell) {
        Some(c) => query_cmp_num(c, op, d)
        None => false
      }
  }
}

///|
/// Compare a numeric cell against a numeric operand.
fn query_cmp_num(cell : Double, op : String, value : Double) -> Bool {
  if op == "==" {
    cell == value
  } else if op == "!=" {
    cell != value
  } else if op == ">" {
    cell > value
  } else if op == ">=" {
    cell >= value
  } else if op == "<" {
    cell < value
  } else if op == "<=" {
    cell <= value
  } else {
    false
  }
}

///|
/// Evaluate the predicate against obs row `i`.
pub fn anndata_query_eval_obs(adata : AnnData, pred : QPred, i : Int32) -> Bool {
  match pred {
    QCmp(col, op, value) =>
      query_cmp_str(query_obs_value(adata, col, i), op, value)
    QAnd(a, b) =>
      anndata_query_eval_obs(adata, a, i) && anndata_query_eval_obs(adata, b, i)
    QOr(a, b) =>
      anndata_query_eval_obs(adata, a, i) || anndata_query_eval_obs(adata, b, i)
    QNot(a) => not(anndata_query_eval_obs(adata, a, i))
  }
}

///|
/// Evaluate the predicate against var column `j`.
pub fn anndata_query_eval_var(adata : AnnData, pred : QPred, j : Int32) -> Bool {
  match pred {
    QCmp(col, op, value) =>
      query_cmp_str(query_var_value(adata, col, j), op, value)
    QAnd(a, b) =>
      anndata_query_eval_var(adata, a, j) && anndata_query_eval_var(adata, b, j)
    QOr(a, b) =>
      anndata_query_eval_var(adata, a, j) || anndata_query_eval_var(adata, b, j)
    QNot(a) => not(anndata_query_eval_var(adata, a, j))
  }
}

///|
/// Query observations: indices of rows matching the expression
/// (`adata.obs.query(expr)` analogue). Returns `[]` on a parse error.
pub fn anndata_obs_query(adata : AnnData, expr : String) -> Array[Int32] {
  match anndata_query_parse(expr) {
    Some(pred) => {
      let out : Array[Int32] = Array(capacity=0)
      for i = 0; i < adata.n_obs(); i = i + 1 {
        if anndata_query_eval_obs(adata, pred, i) {
          out.push(i)
        }
      }
      out
    }
    None => []
  }
}

///|
/// Query variables: indices of columns matching the expression.
pub fn anndata_var_query(adata : AnnData, expr : String) -> Array[Int32] {
  match anndata_query_parse(expr) {
    Some(pred) => {
      let out : Array[Int32] = Array(capacity=0)
      for j = 0; j < adata.n_vars(); j = j + 1 {
        if anndata_query_eval_var(adata, pred, j) {
          out.push(j)
        }
      }
      out
    }
    None => []
  }
}

///|
/// Query + subset in one step: returns the AnnData restricted to the
/// observations matching the expression (an empty AnnData on parse error).
pub fn anndata_query_subset(adata : AnnData, expr : String) -> AnnData {
  adata.subset(anndata_obs_query(adata, expr), range(adata.n_vars()))
}

///|
/// Query + subset along the var axis.
pub fn anndata_query_subset_var(adata : AnnData, expr : String) -> AnnData {
  adata.subset(range(adata.n_obs()), anndata_var_query(adata, expr))
}
