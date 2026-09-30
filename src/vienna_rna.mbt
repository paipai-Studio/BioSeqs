// ============================================================================
// ViennaRNA: partition function, base pair probability, MEA, and utilities
// ============================================================================

///|
/// Check if two RNA bases can form a canonical pair (AU, GC, GU wobble).
fn vienna_can_pair(b1 : Char, b2 : Char) -> Bool {
  let is_au = (b1 == 'A' || b1 == 'a') && (b2 == 'U' || b2 == 'u')
  let is_ua = (b2 == 'A' || b2 == 'a') && (b1 == 'U' || b1 == 'u')
  let is_gc = (b1 == 'G' || b1 == 'g') && (b2 == 'C' || b2 == 'c')
  let is_cg = (b2 == 'G' || b2 == 'g') && (b1 == 'C' || b1 == 'c')
  let is_gu = (b1 == 'G' || b1 == 'g') && (b2 == 'U' || b2 == 'u')
  let is_ug = (b2 == 'G' || b2 == 'g') && (b1 == 'U' || b1 == 'u')
  is_au || is_ua || is_gc || is_cg || is_gu || is_ug
}

///|
/// Simple base pair energy: AU=-1, GC=-2, GU=-0.5 (kcal/mol).
fn vienna_bp_energy(b1 : Char, b2 : Char) -> Double {
  let is_gc = ((b1 == 'G' || b1 == 'g') && (b2 == 'C' || b2 == 'c')) ||
    ((b2 == 'G' || b2 == 'g') && (b1 == 'C' || b1 == 'c'))
  let is_au = ((b1 == 'A' || b1 == 'a') && (b2 == 'U' || b2 == 'u')) ||
    ((b2 == 'A' || b2 == 'a') && (b1 == 'U' || b1 == 'u'))
  if is_gc { -2.0 } else if is_au { -1.0 } else { -0.5 }
}

///|
/// Compute the partition function for RNA secondary structure.
/// Uses a simplified energy model with Boltzmann weights.
/// Returns the n×n partition function matrix Q where Q[i][j] is the
/// partition function for the subsequence i..j.
pub fn vienna_partition_function(
  sequence : String,
  temperature : Double,
  min_loop : Int,
) -> Array[Array[Double]] {
  let seq = sequence.to_array()
  let n = seq.length()
  let q : Array[Array[Double]] = Array::makei(n, fn(_) { Array::make(n, 1.0) })
  let kt = 0.00198717 * temperature
  let mut length = 0
  while length <= n - 1 {
    let mut i = 0
    while i + length < n {
      let j = i + length
      let mut total = 1.0
      if length > min_loop {
        let mut k = i
        while k < j {
          if vienna_can_pair(seq[k], seq[j]) {
            let energy = vienna_bp_energy(seq[k], seq[j])
            let boltzmann = @math.exp(-energy / kt)
            let left = if k > i { q[i][k - 1] } else { 1.0 }
            let inner = if k + 1 <= j - 1 { q[k + 1][j - 1] } else { 1.0 }
            total = total + left * inner * boltzmann
          }
          k = k + 1
        }
      }
      q[i][j] = total
      i = i + 1
    }
    length = length + 1
  }
  q
}

///|
/// Compute base pair probability matrix from the partition function.
/// Returns an n×n matrix where P[i][j] = probability that i and j are paired.
pub fn vienna_bp_probability(
  sequence : String,
  temperature : Double,
  min_loop : Int,
) -> Array[Array[Double]] {
  let seq = sequence.to_array()
  let n = seq.length()
  if n < 2 {
    return []
  }
  let q = vienna_partition_function(sequence, temperature, min_loop)
  let z = q[0][n - 1]
  if z <= 0.0 {
    return Array::makei(n, fn(_) { Array::make(n, 0.0) })
  }
  let kt = 0.00198717 * temperature
  let prob : Array[Array[Double]] = Array::makei(n, fn(_) { Array::make(n, 0.0) })
  let mut i = 0
  while i < n {
    let mut j = i + min_loop + 1
    while j < n {
      if vienna_can_pair(seq[i], seq[j]) {
        let energy = vienna_bp_energy(seq[i], seq[j])
        let boltzmann = @math.exp(-energy / kt)
        let inner = if i + 1 <= j - 1 { q[i + 1][j - 1] } else { 1.0 }
        let outside_l = if i > 0 { q[0][i - 1] } else { 1.0 }
        let outside_r = if j < n - 1 { q[j + 1][n - 1] } else { 1.0 }
        prob[i][j] = outside_l * boltzmann * inner * outside_r / z
        prob[j][i] = prob[i][j]
      }
      j = j + 1
    }
    i = i + 1
  }
  prob
}

///|
/// Maximum Expected Accuracy (MEA) structure from base pair probabilities.
/// `gamma` is the pairing weight (higher → more pairs).
/// Returns dot-bracket notation.
pub fn vienna_mea(
  bp_prob : Array[Array[Double]],
  gamma : Double,
) -> String {
  let n = bp_prob.length()
  if n == 0 {
    return ""
  }
  let m : Array[Array[Double]] = Array::makei(n, fn(_) { Array::make(n, 0.0) })
  let bt : Array[Array[Int]] = Array::makei(n, fn(_) { Array::make(n, -1) })
  let mut length = 0
  while length < n {
    let mut i = 0
    while i + length < n {
      let j = i + length
      let mut best = if i + 1 <= j { m[i + 1][j] } else { 0.0 }
      let mut best_k = -1
      let mut k = i + 1
      while k <= j {
        let pp = bp_prob[i][k] * gamma
        let inner = if i + 1 <= k - 1 { m[i + 1][k - 1] } else { 0.0 }
        let right = if k + 1 <= j { m[k + 1][j] } else { 0.0 }
        let val = pp + inner + right
        if val > best {
          best = val
          best_k = k
        }
        k = k + 1
      }
      m[i][j] = best
      bt[i][j] = best_k
      i = i + 1
    }
    length = length + 1
  }
  // Backtrack
  let pairs : Array[Int] = Array::make(n, -1)
  let stack : Array[(Int, Int)] = [(0, n - 1)]
  while stack.length() > 0 {
    let (i, j) = stack[stack.length() - 1]
    ignore(stack.remove(stack.length() - 1))
    if i > j {
      continue
    }
    let k = bt[i][j]
    if k == -1 {
      if i + 1 <= j {
        stack.push((i + 1, j))
      }
    } else {
      pairs[i] = k
      pairs[k] = i
      if i + 1 <= k - 1 {
        stack.push((i + 1, k - 1))
      }
      if k + 1 <= j {
        stack.push((k + 1, j))
      }
    }
  }
  // Convert to dot-bracket
  let result : Array[Char] = Array::make(n, '.')
  for i = 0; i < n; i = i + 1 {
    if pairs[i] > i {
      result[i] = '('
      result[pairs[i]] = ')'
    }
  }
  let buf : Array[String] = []
  for c in result {
    buf.push(c.to_string())
  }
  buf.join("")
}

///|
/// Generate mountain plot values from dot-bracket notation.
/// Returns an array of (position, height) pairs.
/// Height at position i = number of base pairs enclosing position i.
pub fn vienna_mountain_plot(
  dot_bracket : String,
) -> Array[(Int, Int)] {
  let s = dot_bracket.to_array()
  let n = s.length()
  let result : Array[(Int, Int)] = Array(capacity=n)
  let mut depth = 0
  for i = 0; i < n; i = i + 1 {
    if s[i] == '(' {
      depth = depth + 1
      result.push((i + 1, depth))
    } else if s[i] == ')' {
      result.push((i + 1, depth))
      depth = depth - 1
    } else {
      result.push((i + 1, depth))
    }
  }
  result
}

///|
/// Compute base pair distance between two secondary structures.
/// d = |S1 \ S2| + |S2 \ S1| = number of pairs in one but not the other.
pub fn vienna_bp_distance(
  s1 : String,
  s2 : String,
) -> Int {
  let a1 = s1.to_array()
  let a2 = s2.to_array()
  let n = a1.length()
  if a2.length() != n {
    return -1
  }
  let pairs1 : Map[(Int, Int), Bool] = Map([])
  let stack1 : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    if a1[i] == '(' {
      stack1.push(i)
    } else if a1[i] == ')' && stack1.length() > 0 {
      let j = stack1[stack1.length() - 1]
      ignore(stack1.remove(stack1.length() - 1))
      pairs1[(j, i)] = true
    }
  }
  let pairs2 : Map[(Int, Int), Bool] = Map([])
  let stack2 : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    if a2[i] == '(' {
      stack2.push(i)
    } else if a2[i] == ')' && stack2.length() > 0 {
      let j = stack2[stack2.length() - 1]
      ignore(stack2.remove(stack2.length() - 1))
      pairs2[(j, i)] = true
    }
  }
  let mut dist = 0
  for k in pairs1.keys() {
    if not(pairs2.contains(k)) {
      dist = dist + 1
    }
  }
  for k in pairs2.keys() {
    if not(pairs1.contains(k)) {
      dist = dist + 1
    }
  }
  dist
}

///|
/// Classify loops in a secondary structure into hairpin, interior, bulge,
/// multiloop, and external (exterior) loops.
/// Returns (hairpin, interior, bulge, multiloop, external) counts.
pub fn vienna_loop_types(
  dot_bracket : String,
) -> (Int, Int, Int, Int, Int) {
  let s = dot_bracket.to_array()
  let n = s.length()
  let pair : Array[Int] = Array::make(n, -1)
  let stack : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    if s[i] == '(' {
      stack.push(i)
    } else if s[i] == ')' && stack.length() > 0 {
      let j = stack[stack.length() - 1]
      ignore(stack.remove(stack.length() - 1))
      pair[i] = j
      pair[j] = i
    }
  }
  let counts : Array[Int] = [0, 0, 0, 0, 0] // hairpin, interior, bulge, multiloop, external
  let mut i = 0
  while i < n {
    if pair[i] > i {
      let j = pair[i]
      // Check if external
      let mut is_ext = true
      let mut m = 0
      while m < n {
        if pair[m] > m && m < i && pair[m] > j {
          is_ext = false
          break
        }
        m = m + 1
      }
      if is_ext {
        counts[4] = counts[4] + 1
      }
      vienna_loop_classify(pair, i, j, counts)
      i = j + 1
    } else {
      i = i + 1
    }
  }
  (counts[0], counts[1], counts[2], counts[3], counts[4])
}

///|
/// Recursively classify loop enclosed by pair (i, j).
fn vienna_loop_classify(
  pair : Array[Int],
  i : Int,
  j : Int,
  counts : Array[Int],
) -> Unit {
  let inner_pairs : Array[(Int, Int)] = []
  let mut k = i + 1
  while k < j {
    if pair[k] > k && pair[k] < j {
      inner_pairs.push((k, pair[k]))
      k = pair[k] + 1
    } else {
      k = k + 1
    }
  }
  if inner_pairs.length() == 0 {
    counts[0] = counts[0] + 1 // hairpin
  } else if inner_pairs.length() == 1 {
    let (ip, jp) = inner_pairs[0]
    let unpaired_left = ip - i - 1
    let unpaired_right = j - jp - 1
    if unpaired_left == 0 && unpaired_right == 0 {
      counts[1] = counts[1] + 1 // interior (stack)
    } else if unpaired_left == 0 || unpaired_right == 0 {
      counts[2] = counts[2] + 1 // bulge
    } else {
      counts[1] = counts[1] + 1 // interior
    }
  } else {
    counts[3] = counts[3] + 1 // multiloop
  }
  for (ip, jp) in inner_pairs {
    vienna_loop_classify(pair, ip, jp, counts)
  }
}

///|
/// Compute simplified free energy of a given structure.
/// AU=-1, GC=-2, GU=-0.5 per pair, +0.5 per unpaired base in hairpin loop.
pub fn vienna_structure_energy(
  sequence : String,
  dot_bracket : String,
) -> Double {
  let seq = sequence.to_array()
  let s = dot_bracket.to_array()
  let n = s.length()
  let pair : Array[Int] = Array::make(n, -1)
  let stack : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    if s[i] == '(' {
      stack.push(i)
    } else if s[i] == ')' && stack.length() > 0 {
      let j = stack[stack.length() - 1]
      ignore(stack.remove(stack.length() - 1))
      pair[i] = j
      pair[j] = i
    }
  }
  let mut energy = 0.0
  let mut i = 0
  while i < n {
    if pair[i] > i {
      let j = pair[i]
      // Base pair energy
      energy = energy + vienna_bp_energy(seq[i], seq[j])
      // Check if this is a hairpin (no inner pairs)
      let mut has_inner = false
      let mut unpaired_count = 0
      let mut k = i + 1
      while k < j {
        if pair[k] > k && pair[k] < j {
          has_inner = true
          k = pair[k] + 1
        } else {
          if pair[k] == -1 {
            unpaired_count = unpaired_count + 1
          }
          k = k + 1
        }
      }
      if not(has_inner) {
        // Hairpin loop penalty
        energy = energy + 0.5 * unpaired_count.to_double()
      }
      i = j + 1
    } else {
      i = i + 1
    }
  }
  energy
}
///|
/// Convert dot-bracket to tree representation (ViennaRNA style).
/// Each base pair becomes an internal node, unpaired bases become leaves.
/// Returns a nested string like "((.)(...))" -> "((U)(UUU))"
pub fn vienna_tree_representation(
  dot_bracket : String,
) -> String {
  let s = dot_bracket.to_array()
  vienna_tree_build(s, 0, s.length())
}

fn vienna_tree_build(
  s : Array[Char],
  i : Int,
  n : Int,
) -> String {
  if i >= n {
    return ""
  }
  if s[i] == '.' {
    return "U" + vienna_tree_build(s, i + 1, n)
  }
  if s[i] == '(' {
    let mut depth = 1
    let mut j = i + 1
    while j < n && depth > 0 {
      if s[j] == '(' {
        depth = depth + 1
      } else if s[j] == ')' {
        depth = depth - 1
      }
      if depth > 0 {
        j = j + 1
      }
    }
    let inner = vienna_tree_build(s, i + 1, j)
    let rest = vienna_tree_build(s, j + 1, n)
    return "(" + inner + ")" + rest
  }
  vienna_tree_build(s, i + 1, n)
}
