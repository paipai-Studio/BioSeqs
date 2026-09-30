///|
/// htslib expression filter (hts_expr.c port).
///
/// Recursive descent parser/evaluator for filter expressions like
/// "QUAL>30 && DP>10", "MAPQ>=20", "FLAG&4==0", etc.
///
/// Supported syntax (by precedence, low to high):
///   && ||          boolean and/or
///   == != =~ !~    equality / regex match
///   > >= < <=      comparison
///   |              bitwise or
///   ^              bitwise xor
///   &              bitwise and
///   + -            addition / subtraction
///   * / %          multiplication / division / modulo
///   + - ! ~        unary plus / minus / not / bitwise not
///   ( )            grouping
///   func()         length, min, max, avg, sqrt, pow, exp, log, exists, default
///   "..."          string literal
///   number         integer or float literal
///   identifier     variable (looked up via symbol table)

///|
/// Expression value: can be numeric, string, or undefined.
pub struct ExprVal {
  is_str : Bool
  is_true : Int
  d : Double
  s : String
}

///|
/// Create a numeric value.
pub fn ExprVal::num(d : Double) -> ExprVal {
  ExprVal::{ is_str: false, is_true: if d != 0.0 { 1 } else { 0 }, d, s: "" }
}

///|
/// Create a string value.
pub fn ExprVal::str(s : String) -> ExprVal {
  ExprVal::{ is_str: true, is_true: if s.length() > 0 { 1 } else { 0 }, d: 0.0, s }
}

///|
/// Create an undefined value.
pub fn ExprVal::undef() -> ExprVal {
  ExprVal::{ is_str: false, is_true: -1, d: 0.0, s: "" }
}

///|
/// Check if value exists (is defined).
fn ExprVal::exists(self : ExprVal) -> Bool {
  self.is_true >= 0
}

///|
/// Check if value exists and is truthy.
fn ExprVal::exists_t(self : ExprVal) -> Bool {
  self.is_true > 0 || (self.is_true == 0 && self.is_str && self.s.length() > 0) || (self.is_true == 0 && not(self.is_str) && self.d != 0.0)
}

///|
/// Expression filter state.
pub struct HtsFilter {
  expr : String
}

///|
/// Create a filter from an expression string.
pub fn hts_filter_init(expr : String) -> HtsFilter {
  HtsFilter::{ expr }
}

///|
/// Parser state for recursive descent parsing.
struct Parser {
  chars : Array[Char]
  mut pos : Int
  symbols : Map[String, ExprVal]
  mut depth : Int
}

///|
/// Skip whitespace.
fn Parser::skip_ws(self : Parser) -> Unit {
  while self.pos < self.chars.length() {
    let c = self.chars[self.pos]
    if c == ' ' || c == '\t' {
      self.pos = self.pos + 1
    } else {
      break
    }
  }
}

///|
/// Peek current character (or '\uFFFF' if end).
fn Parser::peek(self : Parser) -> Char {
  if self.pos < self.chars.length() {
    self.chars[self.pos]
  } else {
    '\uFFFF'
  }
}

///|
/// Peek next character (or '\uFFFF' if end).
fn Parser::peek2(self : Parser) -> Char {
  if self.pos + 1 < self.chars.length() {
    self.chars[self.pos + 1]
  } else {
    '\uFFFF'
  }
}

///|
/// Parse a number (integer or float).
fn Parser::parse_number(self : Parser) -> ExprVal? {
  let start = self.pos
  let n = self.chars.length()
  // Optional sign handled by unary_expr
  // Digits
  while self.pos < n && self.chars[self.pos] >= '0' && self.chars[self.pos] <= '9' {
    self.pos = self.pos + 1
  }
  // Fractional part
  if self.pos < n && self.chars[self.pos] == '.' {
    self.pos = self.pos + 1
    while self.pos < n && self.chars[self.pos] >= '0' && self.chars[self.pos] <= '9' {
      self.pos = self.pos + 1
    }
  }
  // Exponent
  if self.pos < n && (self.chars[self.pos] == 'e' || self.chars[self.pos] == 'E') {
    self.pos = self.pos + 1
    if self.pos < n && (self.chars[self.pos] == '+' || self.chars[self.pos] == '-') {
      self.pos = self.pos + 1
    }
    while self.pos < n && self.chars[self.pos] >= '0' && self.chars[self.pos] <= '9' {
      self.pos = self.pos + 1
    }
  }
  if self.pos == start {
    return None
  }
  let sb = StringBuilder::new()
  for i = start; i < self.pos; i = i + 1 {
    sb.write_char(self.chars[i])
  }
  let num_str = sb.to_string()
  match parse_double(num_str) {
    Some(v) => Some(ExprVal::num(v))
    None => None
  }
}

///|
/// Parse a string literal (double-quoted).
fn Parser::parse_string(self : Parser) -> ExprVal? {
  if self.peek() != '"' {
    return None
  }
  self.pos = self.pos + 1
  let result = StringBuilder::new()
  while self.pos < self.chars.length() && self.chars[self.pos] != '"' {
    if self.chars[self.pos] == '\\' && self.pos + 1 < self.chars.length() {
      self.pos = self.pos + 1
      match self.chars[self.pos] {
        '"' => result.write_char('"')
        '\\' => result.write_char('\\')
        't' => result.write_char('\t')
        'n' => result.write_char('\n')
        'r' => result.write_char('\r')
        c => {
          result.write_char('\\')
          result.write_char(c)
        }
      }
    } else {
      result.write_char(self.chars[self.pos])
    }
    self.pos = self.pos + 1
  }
  if self.peek() != '"' {
    return None
  }
  self.pos = self.pos + 1
  Some(ExprVal::str(result.to_string()))
}

///|
/// Parse an identifier (variable name or function name).
fn Parser::parse_identifier(self : Parser) -> String {
  let start = self.pos
  let n = self.chars.length()
  while self.pos < n {
    let c = self.chars[self.pos]
    if (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' {
      self.pos = self.pos + 1
    } else {
      break
    }
  }
  {
    let sb2 = StringBuilder::new()
    for i = start; i < self.pos; i = i + 1 {
      sb2.write_char(self.chars[i])
    }
    sb2.to_string()
  }
}

///|
/// Look up a variable from the symbol table.
fn Parser::lookup(self : Parser, name : String) -> ExprVal? {
  match self.symbols.get(name) {
    Some(v) => Some(v)
    None => None
  }
}

///|
/// Parse a function call: func(expr) or func(expr, expr).
fn Parser::parse_func(self : Parser, name : String) -> ExprVal? {
  // Already consumed name and '('
  let arg1 = self.parse_expression()
  if arg1.is_empty() {
    return None
  }
  let mut result = arg1[0]
  // Check for comma (for two-arg functions like pow, default)
  let mut arg2 : ExprVal? = None
  self.skip_ws()
  if self.peek() == ',' {
    self.pos = self.pos + 1
    self.skip_ws()
    let a2 = self.parse_expression()
    if a2.is_empty() {
      return None
    }
    arg2 = Some(a2[0])
  }
  self.skip_ws()
  if self.peek() != ')' {
    return None
  }
  self.pos = self.pos + 1
  // Apply function
  match name {
    "length" => {
      if not(result.is_str) {
        return None
      }
      result = ExprVal::num(result.s.length().to_double())
    }
    "min" => {
      if not(result.is_str) {
        return None
      }
      let mut v = 999999
      for i = 0; i < result.s.length(); i = i + 1 {
        let c = result.s.unsafe_get(i).to_int()
        if c < v {
          v = c
        }
      }
      result = if v == 999999 { ExprVal::undef() } else { ExprVal::num(v.to_double()) }
    }
    "max" => {
      if not(result.is_str) {
        return None
      }
      let mut v = -1
      for i = 0; i < result.s.length(); i = i + 1 {
        let c = result.s.unsafe_get(i).to_int()
        if c > v {
          v = c
        }
      }
      result = if v == -1 { ExprVal::undef() } else { ExprVal::num(v.to_double()) }
    }
    "avg" => {
      if not(result.is_str) {
        return None
      }
      let mut sum = 0.0
      for i = 0; i < result.s.length(); i = i + 1 {
        sum = sum + result.s.unsafe_get(i).to_int().to_double()
      }
      let avg = if result.s.length() > 0 { sum / result.s.length().to_double() } else { 0.0 }
      result = ExprVal::num(avg)
    }
    "sqrt" => {
      if result.is_str {
        return None
      }
      result = ExprVal::num(result.d.sqrt())
    }
    "exp" => {
      if result.is_str {
        return None
      }
      result = ExprVal::num(@math.exp(result.d))
    }
    "log" => {
      if result.is_str {
        return None
      }
      result = ExprVal::num(@math.ln(result.d))
    }
    "pow" => {
      match arg2 {
        Some(v2) => {
          if result.is_str || v2.is_str {
            return None
          }
          result = ExprVal::num(@math.pow(result.d, v2.d))
        }
        None => return None
      }
    }
    "exists" => {
      let exists = result.exists()
      result = ExprVal::{ is_str: false, is_true: if exists { 1 } else { 0 }, d: if exists { 1.0 } else { 0.0 }, s: "" }
    }
    "default" => {
      match arg2 {
        Some(v2) => {
          if not(result.exists()) {
            result = v2
          }
        }
        None => return None
      }
    }
    _ => return None
  }
  Some(result)
}

///|
/// Parse simple expression: number, string, func(), (expr), variable.
fn Parser::parse_simple(self : Parser) -> ExprVal? {
  self.skip_ws()
  let c = self.peek()
  // Parenthesized expression
  if c == '(' {
    self.pos = self.pos + 1
    let result = self.parse_expression()
    if result.is_empty() {
      return None
    }
    self.skip_ws()
    if self.peek() != ')' {
      return None
    }
    self.pos = self.pos + 1
    return Some(result[0])
  }
  // String literal
  if c == '"' {
    return self.parse_string()
  }
  // Number
  if (c >= '0' && c <= '9') || c == '.' {
    return self.parse_number()
  }
  // Identifier (variable or function)
  if (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_' {
    let name = self.parse_identifier()
    self.skip_ws()
    if self.peek() == '(' {
      self.pos = self.pos + 1
      return self.parse_func(name)
    }
    // Variable lookup
    match self.lookup(name) {
      Some(v) => return Some(v)
      None => return Some(ExprVal::undef())
    }
  }
  None
}

///|
/// Parse unary expression: +x, -x, !x, ~x, or simple_expr.
fn Parser::parse_unary(self : Parser) -> ExprVal? {
  self.skip_ws()
  let c = self.peek()
  if c == '+' || c == '-' {
    self.pos = self.pos + 1
    match self.parse_simple() {
      Some(v0) => {
        let mut v = v0
        if not(v.exists()) {
          return Some(ExprVal::undef())
        }
        if v.is_str {
          return None
        }
        if c == '-' {
          v = ExprVal::num(-v.d)
        }
        Some(v)
      }
      None => None
    }
  } else if c == '!' {
    self.pos = self.pos + 1
    match self.parse_unary() {
      Some(v) => {
        if v.is_true > 0 {
          Some(ExprVal::num(0.0))
        } else if not(v.exists()) {
          Some(ExprVal::num(1.0))
        } else if v.is_str {
          Some(ExprVal::num(if v.s.length() == 0 { 1.0 } else { 0.0 }))
        } else {
          Some(ExprVal::num(if v.d == 0.0 { 1.0 } else { 0.0 }))
        }
      }
      None => None
    }
  } else if c == '~' {
    self.pos = self.pos + 1
    match self.parse_unary() {
      Some(v) => {
        if not(v.exists()) || v.is_str {
          return Some(ExprVal::undef())
        }
        Some(ExprVal::num((v.d.to_int().lnot()).to_double()))
      }
      None => None
    }
  } else {
    self.parse_simple()
  }
}

///|
/// Parse multiplication expression: unary (* | / | % unary)*
fn Parser::parse_mul(self : Parser) -> ExprVal? {
  match self.parse_unary() {
    Some(result0) => {
      let mut result = result0
      while true {
        self.skip_ws()
        let c = self.peek()
        if c == '*' || c == '/' || c == '%' {
          self.pos = self.pos + 1
          match self.parse_unary() {
            Some(v) => {
              if not(v.exists()) || not(result.exists()) || v.is_str || result.is_str {
                return Some(ExprVal::undef())
              }
              match c {
                '*' => result = ExprVal::num(result.d * v.d)
                '/' => result = ExprVal::num(result.d / v.d)
                '%' => {
                  if v.d == 0.0 {
                    return Some(ExprVal::undef())
                  }
                  result = ExprVal::num((result.d.to_int() % v.d.to_int()).to_double())
                }
                _ => break
              }
            }
            None => return None
          }
        } else {
          break
        }
      }
      Some(result)
    }
    None => None
  }
}

///|
/// Parse addition expression: mul (+ | - mul)*
fn Parser::parse_add(self : Parser) -> ExprVal? {
  match self.parse_mul() {
    Some(result0) => {
      let mut result = result0
      while true {
        self.skip_ws()
        let c = self.peek()
        if c == '+' || c == '-' {
          self.pos = self.pos + 1
          match self.parse_mul() {
            Some(v) => {
              if not(v.exists()) || not(result.exists()) || v.is_str || result.is_str {
                return Some(ExprVal::undef())
              }
              match c {
                '+' => result = ExprVal::num(result.d + v.d)
                '-' => result = ExprVal::num(result.d - v.d)
                _ => break
              }
            }
            None => return None
          }
        } else {
          break
        }
      }
      Some(result)
    }
    None => None
  }
}

///|
/// Parse bitwise and: add (& add)*  (but not &&)
fn Parser::parse_bitand(self : Parser) -> ExprVal? {
  match self.parse_add() {
    Some(result0) => {
      let mut result = result0
      while true {
        self.skip_ws()
        if self.peek() == '&' && self.peek2() != '&' {
          self.pos = self.pos + 1
          match self.parse_add() {
            Some(v) => {
              if not(v.exists()) || not(result.exists()) || v.is_str || result.is_str {
                return Some(ExprVal::undef())
              }
              result = ExprVal::num((result.d.to_int() & v.d.to_int()).to_double())
            }
            None => return None
          }
        } else {
          break
        }
      }
      Some(result)
    }
    None => None
  }
}

///|
/// Parse bitwise xor: bitand (^ bitand)*
fn Parser::parse_bitxor(self : Parser) -> ExprVal? {
  match self.parse_bitand() {
    Some(result0) => {
      let mut result = result0
      while true {
        self.skip_ws()
        if self.peek() == '^' {
          self.pos = self.pos + 1
          match self.parse_bitand() {
            Some(v) => {
              if not(v.exists()) || not(result.exists()) || v.is_str || result.is_str {
                return Some(ExprVal::undef())
              }
              result = ExprVal::num((result.d.to_int() ^ v.d.to_int()).to_double())
            }
            None => return None
          }
        } else {
          break
        }
      }
      Some(result)
    }
    None => None
  }
}

///|
/// Parse bitwise or: bitxor (| bitxor)*  (but not ||)
fn Parser::parse_bitor(self : Parser) -> ExprVal? {
  match self.parse_bitxor() {
    Some(result0) => {
      let mut result = result0
      while true {
        self.skip_ws()
        if self.peek() == '|' && self.peek2() != '|' {
          self.pos = self.pos + 1
          match self.parse_bitxor() {
            Some(v) => {
              if not(v.exists()) || not(result.exists()) || v.is_str || result.is_str {
                return Some(ExprVal::undef())
              }
              result = ExprVal::num((result.d.to_int() | v.d.to_int()).to_double())
            }
            None => return None
          }
        } else {
          break
        }
      }
      Some(result)
    }
    None => None
  }
}

///|
/// Parse comparison: bitor (< | <= | > | >= bitor)?
fn Parser::parse_cmp(self : Parser) -> ExprVal? {
  match self.parse_bitor() {
    Some(result) => {
      self.skip_ws()
      let c = self.peek()
      let c2 = self.peek2()
      // >=
      if c == '>' && c2 == '=' {
        self.pos = self.pos + 2
        match self.parse_bitor() {
          Some(v) => {
            if not(v.exists()) || not(result.exists()) {
              return Some(ExprVal::undef())
            }
            if result.is_str && v.is_str {
              let cmp = compare(result.s, v.s)
              return Some(ExprVal::num(if cmp >= 0 { 1.0 } else { 0.0 }))
            }
            if result.is_str || v.is_str {
              return Some(ExprVal::num(0.0))
            }
            Some(ExprVal::num(if result.d >= v.d { 1.0 } else { 0.0 }))
          }
          None => None
        }
      } else if c == '>' {
        self.pos = self.pos + 1
        match self.parse_bitor() {
          Some(v) => {
            if not(v.exists()) || not(result.exists()) {
              return Some(ExprVal::undef())
            }
            if result.is_str && v.is_str {
              let cmp = compare(result.s, v.s)
              return Some(ExprVal::num(if cmp > 0 { 1.0 } else { 0.0 }))
            }
            if result.is_str || v.is_str {
              return Some(ExprVal::num(0.0))
            }
            Some(ExprVal::num(if result.d > v.d { 1.0 } else { 0.0 }))
          }
          None => None
        }
      } else if c == '<' && c2 == '=' {
        self.pos = self.pos + 2
        match self.parse_bitor() {
          Some(v) => {
            if not(v.exists()) || not(result.exists()) {
              return Some(ExprVal::undef())
            }
            if result.is_str && v.is_str {
              let cmp = compare(result.s, v.s)
              return Some(ExprVal::num(if cmp <= 0 { 1.0 } else { 0.0 }))
            }
            if result.is_str || v.is_str {
              return Some(ExprVal::num(0.0))
            }
            Some(ExprVal::num(if result.d <= v.d { 1.0 } else { 0.0 }))
          }
          None => None
        }
      } else if c == '<' {
        self.pos = self.pos + 1
        match self.parse_bitor() {
          Some(v) => {
            if not(v.exists()) || not(result.exists()) {
              return Some(ExprVal::undef())
            }
            if result.is_str && v.is_str {
              let cmp = compare(result.s, v.s)
              return Some(ExprVal::num(if cmp < 0 { 1.0 } else { 0.0 }))
            }
            if result.is_str || v.is_str {
              return Some(ExprVal::num(0.0))
            }
            Some(ExprVal::num(if result.d < v.d { 1.0 } else { 0.0 }))
          }
          None => None
        }
      } else {
        Some(result)
      }
    }
    None => None
  }
}

///|
/// Parse equality: cmp (== | != cmp)?
fn Parser::parse_eq(self : Parser) -> ExprVal? {
  match self.parse_cmp() {
    Some(result) => {
      self.skip_ws()
      let c = self.peek()
      let c2 = self.peek2()
      if c == '=' && c2 == '=' {
        self.pos = self.pos + 2
        match self.parse_cmp() {
          Some(v) => {
            if not(v.exists()) || not(result.exists()) {
              return Some(ExprVal::undef())
            }
            if result.is_str && v.is_str {
              return Some(ExprVal::num(if result.s == v.s { 1.0 } else { 0.0 }))
            }
            if result.is_str || v.is_str {
              return Some(ExprVal::num(0.0))
            }
            Some(ExprVal::num(if result.d == v.d { 1.0 } else { 0.0 }))
          }
          None => None
        }
      } else if c == '!' && c2 == '=' {
        self.pos = self.pos + 2
        match self.parse_cmp() {
          Some(v) => {
            if not(v.exists()) || not(result.exists()) {
              return Some(ExprVal::undef())
            }
            if result.is_str && v.is_str {
              return Some(ExprVal::num(if result.s != v.s { 1.0 } else { 0.0 }))
            }
            if result.is_str || v.is_str {
              return Some(ExprVal::num(1.0))
            }
            Some(ExprVal::num(if result.d != v.d { 1.0 } else { 0.0 }))
          }
          None => None
        }
      } else {
        Some(result)
      }
    }
    None => None
  }
}

///|
/// Parse boolean and/or: eq (&& | || eq)*
fn Parser::parse_and(self : Parser) -> ExprVal? {
  match self.parse_eq() {
    Some(result0) => {
      let mut result = result0
      while true {
        self.skip_ws()
        let c = self.peek()
        let c2 = self.peek2()
        if c == '&' && c2 == '&' {
          self.pos = self.pos + 2
          match self.parse_eq() {
            Some(v) => {
              let lhs = result.exists_t()
              let rhs = v.exists_t()
              result = ExprVal::num(if lhs && rhs { 1.0 } else { 0.0 })
            }
            None => return None
          }
        } else if c == '|' && c2 == '|' {
          self.pos = self.pos + 2
          match self.parse_eq() {
            Some(v) => {
              let lhs = result.exists_t()
              let rhs = v.exists_t()
              result = ExprVal::num(if lhs || rhs { 1.0 } else { 0.0 })
            }
            None => return None
          }
        } else {
          break
        }
      }
      Some(result)
    }
    None => None
  }
}

///|
/// Parse full expression (top level).
fn Parser::parse_expression(self : Parser) -> Array[ExprVal] {
  self.depth = self.depth + 1
  if self.depth > 100 {
    return []
  }
  let result = match self.parse_and() {
    Some(v) => [v]
    None => []
  }
  self.depth = self.depth - 1
  result
}

///|
/// Evaluate a filter expression against a symbol table.
///
/// Returns the result as an ExprVal. The expression is parsed and
/// evaluated in one pass. Variables are looked up from the symbols map.
///
/// Example:
/// ```
/// let filter = hts_filter_init("QUAL > 30 && DP > 10")
/// let symbols = {"QUAL": ExprVal::num(40.0), "DP": ExprVal::num(15.0)}
/// let result = hts_filter_eval(filter, symbols)
/// // result.d == 1.0 (true)
/// ```
pub fn hts_filter_eval(
  filter : HtsFilter,
  symbols : Map[String, ExprVal],
) -> ExprVal {
  let chars = filter.expr.to_array()
  let parser = Parser::{ chars, pos: 0, symbols, depth: 0 }
  let result = parser.parse_expression()
  if result.is_empty() {
    return ExprVal::undef()
  }
  let v = result[0]
  // Final evaluation: strings are true if non-null, numbers are true if nonzero
  if v.is_str {
    let is_true = if v.is_true > 0 { 1 } else { 1 }  // strings are true
    ExprVal::{ is_str: true, is_true, d: is_true.to_double(), s: v.s }
  } else if v.exists() {
    let is_true = if v.d != 0.0 { 1 } else { 0 }
    ExprVal::{ is_str: false, is_true, d: v.d, s: "" }
  } else {
    v
  }
}

///|
/// Convenience: evaluate a filter expression and return true if it passes.
pub fn hts_filter_test(
  filter : HtsFilter,
  symbols : Map[String, ExprVal],
) -> Bool {
  let result = hts_filter_eval(filter, symbols)
  result.exists_t()
}
