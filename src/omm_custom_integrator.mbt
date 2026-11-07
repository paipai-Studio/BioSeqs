/// Custom integrator for user-defined integration schemes.
///
/// OpenMM equivalent: `CustomIntegrator`
///
/// Allows users to define integration steps using simple expressions,
/// enabling non-standard integration schemes.

///|
/// A single step in a custom integration scheme
pub struct OmmIntegratorStep {
  variable : String // target variable: "x", "v", or custom
  expression : String // expression to compute
}

///|
/// Custom integrator
pub struct OmmCustomIntegrator {
  dt : Double
  steps : Array[OmmIntegratorStep]
  global_variables : Array[(String, Double)]
}

///|
pub fn OmmCustomIntegrator::new(dt : Double) -> OmmCustomIntegrator {
  { dt, steps: [], global_variables: [], }
}

///|
pub fn OmmCustomIntegrator::add_step(
  self : OmmCustomIntegrator,
  variable~ : String,
  expression~ : String,
) -> Unit {
  self.steps.push({ variable, expression, })
}

///|
pub fn OmmCustomIntegrator::add_global_variable(
  self : OmmCustomIntegrator,
  name : String,
  value : Double,
) -> Unit {
  self.global_variables.push((name, value))
}

///|
/// Look up a variable name and return its value.
fn omm_lookup_var(
  name : String,
  x : Double,
  v : Double,
  f : Double,
  m : Double,
  dt : Double,
  globals : Array[(String, Double)],
) -> Double {
  if name == "x" {
    x
  } else if name == "v" {
    v
  } else if name == "f" {
    f
  } else if name == "m" {
    m
  } else if name == "dt" {
    dt
  } else {
    for g in globals {
      if g.0 == name {
        return g.1
      }
    }
    0.0
  }
}

///|
/// Tokenize expression into numbers, operators, and variable names
fn omm_tokenize(expr : String) -> Array[String] {
  let tokens : Array[String] = Array(capacity=0)
  let mut current = ""
  for c in expr.to_array() {
    if c == ' ' || c == '\t' {
      if current.length() > 0 {
        tokens.push(current)
        current = ""
      }
    } else if c == '+' ||
      c == '-' ||
      c == '*' ||
      c == '/' ||
      c == '(' ||
      c == ')' {
      if current.length() > 0 {
        tokens.push(current)
        current = ""
      }
      tokens.push(c.to_string())
    } else {
      current = current + c.to_string()
    }
  }
  if current.length() > 0 {
    tokens.push(current)
  }
  tokens
}

///|
/// Get operator precedence
fn omm_op_prec(op : String) -> Int32 {
  if op == "+" || op == "-" {
    1
  } else if op == "*" || op == "/" {
    2
  } else {
    0
  }
}

///|
/// Check if token is an operator
fn omm_is_operator(token : String) -> Bool {
  token == "+" || token == "-" || token == "*" || token == "/"
}

///|
/// Try to parse a string as a double; returns Option
fn omm_try_parse_double(s : String) -> Double? {
  let mut result = 0.0
  let mut sign = 1.0
  let mut i = 0
  let chars = s.to_array()
  let n = chars.length()
  if n == 0 {
    return None
  }
  if chars[0] == '-' {
    sign = -1.0
    i = 1
  } else if chars[0] == '+' {
    i = 1
  }
  let mut has_digit = false
  let mut has_dot = false
  while i < n {
    let c = chars[i]
    if c >= '0' && c <= '9' {
      result = result * 10.0 + (c.to_int() - 48).to_double()
      has_digit = true
    } else if c == '.' && !has_dot {
      has_dot = true
    } else {
      return None
    }
    i = i + 1
  }
  if !has_digit {
    return None
  }
  // Handle decimal part
  if has_dot {
    let mut int_part = 0.0
    let mut frac_part = 0.0
    let mut in_frac = false
    let mut frac_div = 1.0
    let mut j = 0
    if chars[0] == '-' || chars[0] == '+' {
      j = 1
    }
    while j < n {
      let c = chars[j]
      if c == '.' {
        in_frac = true
      } else if c >= '0' && c <= '9' {
        if in_frac {
          frac_div = frac_div * 10.0
          frac_part = frac_part + (c.to_int() - 48).to_double() / frac_div
        } else {
          int_part = int_part * 10.0 + (c.to_int() - 48).to_double()
        }
      }
      j = j + 1
    }
    result = int_part + frac_part
  }
  Some(sign * result)
}

///|
/// Evaluate tokenized expression using Shunting Yard + RPN eval
fn omm_eval_tokens(
  tokens : Array[String],
  x : Double,
  v : Double,
  f : Double,
  m : Double,
  dt : Double,
  globals : Array[(String, Double)],
) -> Double {
  // Convert to RPN using Shunting Yard
  let output : Array[String] = Array(capacity=0)
  let operators : Array[String] = Array(capacity=0)
  for token in tokens {
    if omm_is_operator(token) {
      while operators.length() > 0 &&
            operators[operators.length() - 1] != "(" &&
            omm_op_prec(operators[operators.length() - 1]) >= omm_op_prec(token) {
        output.push(operators.remove(operators.length() - 1))
      }
      operators.push(token)
    } else if token == "(" {
      operators.push(token)
    } else if token == ")" {
      while operators.length() > 0 && operators[operators.length() - 1] != "(" {
        output.push(operators.remove(operators.length() - 1))
      }
      if operators.length() > 0 {
        ignore(operators.remove(operators.length() - 1))
      }
    } else {
      output.push(token)
    }
  }
  while operators.length() > 0 {
    output.push(operators.remove(operators.length() - 1))
  }
  // Evaluate RPN
  let stack : Array[Double] = Array(capacity=0)
  for token in output {
    if omm_is_operator(token) {
      if stack.length() < 2 {
        return 0.0
      }
      let b = stack.remove(stack.length() - 1)
      let a = stack.remove(stack.length() - 1)
      if token == "+" {
        stack.push(a + b)
      } else if token == "-" {
        stack.push(a - b)
      } else if token == "*" {
        stack.push(a * b)
      } else if b != 0.0 {
        stack.push(a / b)
      } else {
        stack.push(0.0)
      }
    } else {
      let parsed = omm_try_parse_double(token)
      match parsed {
        Some(num) => stack.push(num)
        None => stack.push(omm_lookup_var(token, x, v, f, m, dt, globals))
      }
    }
  }
  if stack.length() > 0 {
    stack[0]
  } else {
    0.0
  }
}

///|
/// Evaluate expression for custom integrator.
fn omm_eval_expr(
  expr : String,
  x : Double,
  v : Double,
  f : Double,
  m : Double,
  dt : Double,
  globals : Array[(String, Double)],
) -> Double {
  let tokens = omm_tokenize(expr)
  omm_eval_tokens(tokens, x, v, f, m, dt, globals)
}

///|
/// Execute one custom integrator step for all particles.
pub fn OmmCustomIntegrator::step(
  self : OmmCustomIntegrator,
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  forces : Array[Array[Double]],
  masses : Array[Double],
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let n = positions.length()
  let pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  let vel : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    pos[i] = [positions[i][0], positions[i][1], positions[i][2]]
    vel[i] = [velocities[i][0], velocities[i][1], velocities[i][2]]
  }
  for s in self.steps {
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        let result = omm_eval_expr(
          s.expression,
          pos[i][d],
          vel[i][d],
          forces[i][d],
          masses[i],
          self.dt,
          self.global_variables,
        )
        if s.variable == "x" {
          pos[i][d] = result
        } else if s.variable == "v" {
          vel[i][d] = result
        }
      }
    }
  }
  (pos, vel)
}
