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
/// Supported math functions (OpenMM-style): unary sqrt/exp/log/sin/cos/
/// tan/abs/floor/ceil/square/step, binary min/max, ternary select.
fn omm_is_function(token : String) -> Bool {
  token == "sqrt" ||
  token == "exp" ||
  token == "log" ||
  token == "sin" ||
  token == "cos" ||
  token == "tan" ||
  token == "abs" ||
  token == "floor" ||
  token == "ceil" ||
  token == "square" ||
  token == "step" ||
  token == "min" ||
  token == "max" ||
  token == "select"
}

///|
/// Arity of a supported function.
fn omm_function_arity(token : String) -> Int32 {
  if token == "min" || token == "max" {
    2
  } else if token == "select" {
    3
  } else {
    1
  }
}

///|
/// Apply a function to popped stack values (in source order).
fn omm_apply_function(fn_name : String, args : Array[Double]) -> Double {
  if fn_name == "sqrt" {
    @math.pow(args[0], 0.5)
  } else if fn_name == "exp" {
    @math.exp(args[0])
  } else if fn_name == "log" {
    @math.ln(args[0])
  } else if fn_name == "sin" {
    @math.sin(args[0])
  } else if fn_name == "cos" {
    @math.cos(args[0])
  } else if fn_name == "tan" {
    @math.tan(args[0])
  } else if fn_name == "abs" {
    args[0].abs()
  } else if fn_name == "floor" {
    @math.floor(args[0])
  } else if fn_name == "ceil" {
    @math.ceil(args[0])
  } else if fn_name == "square" {
    args[0] * args[0]
  } else if fn_name == "step" {
    if args[0] >= 0.0 {
      1.0
    } else {
      0.0
    }
  } else if fn_name == "min" {
    if args[0] < args[1] {
      args[0]
    } else {
      args[1]
    }
  } else if fn_name == "max" {
    if args[0] > args[1] {
      args[0]
    } else {
      args[1]
    }
  } else if fn_name == "select" {
    // select(x, y, z) = x > 0 ? y : z
    if args[0] > 0.0 {
      args[1]
    } else {
      args[2]
    }
  } else {
    0.0
  }
}

///|
/// Tokenize expression into numbers, operators, function names, variables,
/// and parentheses/commas. Unary minus becomes the "neg" token.
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
      c == ')' ||
      c == ',' {
      if current.length() > 0 {
        tokens.push(current)
        current = ""
      }
      if c == '-' {
        // Unary minus: at the start, after an operator, "(", or ",".
        let unary = tokens.length() == 0 ||
          omm_is_operator(tokens[tokens.length() - 1]) ||
          tokens[tokens.length() - 1] == "(" ||
          tokens[tokens.length() - 1] == "," ||
          omm_is_function(tokens[tokens.length() - 1])
        tokens.push(if unary { "neg" } else { "-" })
      } else {
        tokens.push(c.to_string())
      }
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
    if token == "neg" {
      // Unary minus: right-binding, treated like a function call.
      operators.push(token)
    } else if omm_is_function(token) {
      // Function names wait for their argument list; highest precedence.
      operators.push(token)
    } else if omm_is_operator(token) {
      while operators.length() > 0 &&
            operators[operators.length() - 1] != "(" &&
            !omm_is_function(operators[operators.length() - 1]) &&
            operators[operators.length() - 1] != "neg" &&
            omm_op_prec(operators[operators.length() - 1]) >= omm_op_prec(token) {
        output.push(operators.remove(operators.length() - 1))
      }
      operators.push(token)
    } else if token == "(" {
      operators.push(token)
    } else if token == "," {
      // Argument separator: flush operators back to the opening "(".
      while operators.length() > 0 && operators[operators.length() - 1] != "(" {
        output.push(operators.remove(operators.length() - 1))
      }
    } else if token == ")" {
      while operators.length() > 0 && operators[operators.length() - 1] != "(" {
        output.push(operators.remove(operators.length() - 1))
      }
      if operators.length() > 0 {
        ignore(operators.remove(operators.length() - 1))
      }
      // The function itself goes to the output after its arguments.
      if operators.length() > 0 &&
        omm_is_function(operators[operators.length() - 1]) {
        output.push(operators.remove(operators.length() - 1))
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
    if token == "neg" {
      if stack.length() < 1 {
        return 0.0
      }
      let a = stack.remove(stack.length() - 1)
      stack.push(-a)
    } else if omm_is_operator(token) {
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
    } else if omm_is_function(token) {
      let arity = omm_function_arity(token)
      if stack.length() < arity {
        return 0.0
      }
      let args : Array[Double] = Array(capacity=0)
      for _k = 0; _k < arity; _k = _k + 1 {
        args.insert(0, stack.remove(stack.length() - 1))
      }
      stack.push(omm_apply_function(token, args))
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
