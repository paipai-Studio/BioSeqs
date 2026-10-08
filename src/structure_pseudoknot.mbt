///|
/// Pseudoknot order detection for RNA secondary structures.
/// Ported from biotite.structure.pseudoknot.
///
/// A pseudoknot occurs when base pairs cross: (i,j) and (k,l) with
/// i < k < j < l. The pseudoknot order is the minimum number of
/// nested sets of non-crossing base pairs needed to cover all pairs.

///|
/// Compute the pseudoknot order of a dot-bracket structure.
/// Returns the order (0 = no pseudoknots, 1 = simple PK, etc.).
///
/// Algorithm: iterative peeling — repeatedly extract the largest
/// non-crossing subset (order 0), remove those pairs, repeat.
/// The number of iterations is the pseudoknot order.
pub fn pseudoknot_order(dot_bracket : String) -> Int {
  let pairs = parse_dot_bracket_pairs(dot_bracket)
  if pairs.length() == 0 {
    return 0
  }
  let mut remaining = pairs
  let mut order = 0
  while remaining.length() > 0 {
    // Extract non-crossing subset (greedy: keep pairs that don't cross any other kept pair)
    let non_crossing = extract_non_crossing(remaining)
    if non_crossing.length() == 0 {
      // Shouldn't happen, but safety
      return order + remaining.length()
    }
    // Remove non_crossing from remaining
    let nc_set : Map[(Int, Int), Bool] = Map([], capacity=non_crossing.length())
    for p in non_crossing {
      nc_set[p] = true
    }
    let new_remaining : Array[(Int, Int)] = Array(capacity=0)
    for p in remaining {
      match nc_set.get(p) {
        Some(_) => () // skip (it's in non_crossing)
        None => new_remaining.push(p)
      }
    }
    remaining = new_remaining
    order = order + 1
  }
  order - 1
}

///|
/// Parse dot-bracket into list of base pairs (i, j) with i < j.
/// Supports single-level pseudoknots with [] and {}.
fn parse_dot_bracket_pairs(db : String) -> Array[(Int, Int)] {
  let pairs : Array[(Int, Int)] = Array(capacity=0)
  let stack_paren : Array[Int] = Array(capacity=0)
  let stack_bracket : Array[Int] = Array(capacity=0)
  let stack_brace : Array[Int] = Array(capacity=0)
  for i = 0; i < db.length(); i = i + 1 {
    let c = db.unsafe_get(i).unsafe_to_char()
    match c {
      '(' => stack_paren.push(i)
      ')' => {
        match stack_paren.pop() {
          Some(j) => pairs.push((j, i))
          None => ()
        }
      }
      '[' => stack_bracket.push(i)
      ']' => {
        match stack_bracket.pop() {
          Some(j) => pairs.push((j, i))
          None => ()
        }
      }
      '{' => stack_brace.push(i)
      '}' => {
        match stack_brace.pop() {
          Some(j) => pairs.push((j, i))
          None => ()
        }
      }
      _ => ()
    }
  }
  // Sort by first index
  for i = 1; i < pairs.length(); i = i + 1 {
    let key = pairs[i]
    let mut j = i - 1
    while j >= 0 && pairs[j].0 > key.0 {
      pairs[j + 1] = pairs[j]
      j = j - 1
    }
    pairs[j + 1] = key
  }
  pairs
}

///|
/// Extract the largest non-crossing subset of base pairs.
/// Two pairs (i,j) and (k,l) cross if i < k < j < l or k < i < l < j.
/// Uses a greedy algorithm: sort by j, then keep pairs that don't cross any kept pair.
fn extract_non_crossing(pairs : Array[(Int, Int)]) -> Array[(Int, Int)] {
  if pairs.length() == 0 {
    return Array(capacity=0)
  }
  // Sort by second element (j) ascending
  let sorted : Array[(Int, Int)] = Array(capacity=pairs.length())
  for p in pairs {
    sorted.push(p)
  }
  for i = 1; i < sorted.length(); i = i + 1 {
    let key = sorted[i]
    let mut j = i - 1
    while j >= 0 && sorted[j].1 > key.1 {
      sorted[j + 1] = sorted[j]
      j = j - 1
    }
    sorted[j + 1] = key
  }
  // Greedy: keep pairs that don't cross any already-kept pair
  let result : Array[(Int, Int)] = Array(capacity=0)
  for p in sorted {
    let mut crosses = false
    for q in result {
      // Check if p and q cross
      let (pi, pj) = p
      let (qi, qj) = q
      if (pi < qi && qi < pj && pj < qj) || (qi < pi && pi < qj && qj < pj) {
        crosses = true
      }
    }
    if !crosses {
      result.push(p)
    }
  }
  result
}

///|
/// Annotate each base pair with its pseudoknot order level.
/// Returns an array of order values, one per base pair (in input order).
pub fn pseudoknot_order_levels(dot_bracket : String) -> Array[Int] {
  let pairs = parse_dot_bracket_pairs(dot_bracket)
  let n = pairs.length()
  if n == 0 {
    return Array(capacity=0)
  }
  let levels : Array[Int] = Array::make(n, 0)
  let remaining : Array[Int] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    remaining.push(i)
  }
  let mut current_level = 0
  while remaining.length() > 0 {
    // Build remaining pairs
    let rem_pairs : Array[(Int, Int)] = Array(capacity=remaining.length())
    for idx in remaining {
      rem_pairs.push(pairs[idx])
    }
    let non_crossing = extract_non_crossing(rem_pairs)
    // Mark non_crossing pairs with current level
    let nc_set : Map[(Int, Int), Bool] = Map([], capacity=non_crossing.length())
    for p in non_crossing {
      nc_set[p] = true
    }
    let new_remaining : Array[Int] = Array(capacity=0)
    for idx in remaining {
      match nc_set.get(pairs[idx]) {
        Some(_) => levels[idx] = current_level
        None => new_remaining.push(idx)
      }
    }
    // Update remaining
    remaining.clear()
    for idx in new_remaining {
      remaining.push(idx)
    }
    current_level = current_level + 1
  }
  levels
}
