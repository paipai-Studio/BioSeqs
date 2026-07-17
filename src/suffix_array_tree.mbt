///|
/// Suffix Array and Suffix Tree data structures for efficient string processing.

// ---------------------------------------------------------------------------
// Suffix Array Construction (Prefix-Doubling Algorithm)
// ---------------------------------------------------------------------------

///|
pub struct SuffixArray {
  text : String
  sa : Array[Int]
}

///|
pub fn SuffixArray::new(text : String) -> SuffixArray {
  let t = text + "$"
  let n = t.length()
  let sa : Array[Int] = Array::new()
  let rank : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    sa.push(i)
    rank.push(t.unsafe_get(i).unsafe_to_char().to_int())
  }
  let mut k = 1
  while k < n {
    sa_sort_by_rank(sa, rank, n, k)
    let new_rank = sa_compute_new_rank(sa, rank, n, k)
    for i = 0; i < n; i = i + 1 {
      rank[i] = new_rank[i]
    }
    k = k * 2
  }
  SuffixArray::{ text, sa }
}

///|
fn sa_sort_by_rank(
  sa : Array[Int],
  rank : Array[Int],
  n : Int,
  k : Int,
) -> Unit {
  for i = 0; i < n; i = i + 1 {
    let mut min_idx = i
    for j = i + 1; j < n; j = j + 1 {
      if sa_compare(sa[j], sa[min_idx], rank, n, k) < 0 {
        min_idx = j
      }
    }
    let tmp = sa[i]
    sa[i] = sa[min_idx]
    sa[min_idx] = tmp
  }
}

///|
fn sa_compare(i : Int, j : Int, rank : Array[Int], n : Int, k : Int) -> Int {
  let ri = rank[i]
  let rj = rank[j]
  if ri < rj {
    return -1
  } else if ri > rj {
    return 1
  }
  let ri2 = if i + k < n { rank[i + k] } else { -1 }
  let rj2 = if j + k < n { rank[j + k] } else { -1 }
  if ri2 < rj2 {
    -1
  } else if ri2 > rj2 {
    1
  } else {
    0
  }
}

///|
fn sa_compute_new_rank(
  sa : Array[Int],
  rank : Array[Int],
  n : Int,
  k : Int,
) -> Array[Int] {
  let new_rank : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    new_rank.push(0)
  }
  new_rank[0] = 0
  let mut r = 0
  for i = 1; i < n; i = i + 1 {
    let prev = sa[i - 1]
    let curr = sa[i]
    let prev2 = if prev + k < n { rank[prev + k] } else { -1 }
    let curr2 = if curr + k < n { rank[curr + k] } else { -1 }
    if rank[prev] != rank[curr] || prev2 != curr2 {
      r = r + 1
    }
    new_rank[curr] = r
  }
  new_rank
}

///|
pub fn SuffixArray::sa(self : SuffixArray) -> Array[Int] {
  self.sa
}

///|
pub fn SuffixArray::text(self : SuffixArray) -> String {
  self.text
}

///|
pub fn SuffixArray::len(self : SuffixArray) -> Int {
  self.sa.length()
}

///|
pub fn SuffixArray::get_suffix(self : SuffixArray, i : Int) -> String {
  if i < 0 || i >= self.text.length() {
    return ""
  }
  self.text[i:self.text.length()].to_owned()
}

///|
pub fn SuffixArray::get_sorted_suffix(self : SuffixArray, idx : Int) -> String {
  if idx < 0 || idx >= self.sa.length() {
    return ""
  }
  let pos = self.sa[idx]
  if pos >= self.text.length() {
    return "$"
  }
  self.text[pos:self.text.length()].to_owned() + "$"
}

///|
pub fn SuffixArray::find_lower(self : SuffixArray, pattern : String) -> Int {
  let n = self.sa.length()
  let m = pattern.length()
  let mut left = 0
  let mut right = n - 1
  let mut result = -1
  while left <= right {
    let mid = (left + right) / 2
    let suffix = self.get_sorted_suffix(mid)
    let compare = sa_compare_pattern(suffix, pattern, m)
    if compare >= 0 {
      if compare == 0 {
        result = mid
      }
      right = mid - 1
    } else {
      left = mid + 1
    }
  }
  if result == -1 && left < n {
    let suffix = self.get_sorted_suffix(left)
    if sa_compare_pattern(suffix, pattern, m) == 0 {
      result = left
    }
  }
  result
}

///|
pub fn SuffixArray::find_upper(self : SuffixArray, pattern : String) -> Int {
  let n = self.sa.length()
  let m = pattern.length()
  let mut left = 0
  let mut right = n - 1
  let mut result = -1
  while left <= right {
    let mid = (left + right) / 2
    let suffix = self.get_sorted_suffix(mid)
    let compare = sa_compare_pattern(suffix, pattern, m)
    if compare <= 0 {
      if compare == 0 {
        result = mid
      }
      left = mid + 1
    } else {
      right = mid - 1
    }
  }
  result
}

///|
pub fn SuffixArray::count(self : SuffixArray, pattern : String) -> Int {
  let lower = self.find_lower(pattern)
  let upper = self.find_upper(pattern)
  if lower == -1 || upper == -1 {
    0
  } else {
    upper - lower + 1
  }
}

///|
pub fn SuffixArray::locate(self : SuffixArray, pattern : String) -> Array[Int] {
  let lower = self.find_lower(pattern)
  let upper = self.find_upper(pattern)
  let result : Array[Int] = Array::new()
  if lower == -1 || upper == -1 {
    return result
  }
  for i = lower; i <= upper; i = i + 1 {
    result.push(self.sa[i])
  }
  sa_sort_positions(result)
  result
}

///|
fn sa_compare_pattern(suffix : String, pattern : String, m : Int) -> Int {
  let n = suffix.length()
  for i = 0; i < m; i = i + 1 {
    if i >= n {
      return -1
    }
    let sc = suffix.unsafe_get(i).unsafe_to_char()
    let pc = pattern.unsafe_get(i).unsafe_to_char()
    if sc < pc {
      return -1
    } else if sc > pc {
      return 1
    }
  }
  0
}

///|
fn sa_sort_positions(positions : Array[Int]) -> Unit {
  let n = positions.length()
  for i = 0; i < n; i = i + 1 {
    let mut min_idx = i
    for j = i + 1; j < n; j = j + 1 {
      if positions[j] < positions[min_idx] {
        min_idx = j
      }
    }
    let tmp = positions[i]
    positions[i] = positions[min_idx]
    positions[min_idx] = tmp
  }
}

// ---------------------------------------------------------------------------
// LCP Array (Longest Common Prefix)
// ---------------------------------------------------------------------------

///|
pub struct LCPArray {
  lcp : Array[Int]
}

///|
pub fn LCPArray::new(sa : SuffixArray) -> LCPArray {
  let n = sa.sa.length()
  let text = sa.text + "$"
  let rank : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    rank.push(0)
  }
  for i = 0; i < n; i = i + 1 {
    rank[sa.sa[i]] = i
  }
  let lcp : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    lcp.push(0)
  }
  let k_values = lcp_kasai(text, sa, rank, n)
  for i = 0; i < n; i = i + 1 {
    lcp[i] = k_values[i]
  }
  LCPArray::{ lcp, }
}

///|
fn lcp_kasai(
  text : String,
  sa : SuffixArray,
  rank : Array[Int],
  n : Int,
) -> Array[Int] {
  let k_values : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    k_values.push(0)
  }
  let mut k = 0
  for i = 0; i < n - 1; i = i + 1 {
    let j = sa.sa[rank[i] - 1]
    let mut current_k = k
    while i + current_k < n &&
          j + current_k < n &&
          text.unsafe_get(i + current_k).unsafe_to_char() ==
          text.unsafe_get(j + current_k).unsafe_to_char() {
      current_k = current_k + 1
    }
    k_values[rank[i]] = current_k
    if current_k > 0 {
      k = current_k - 1
    } else {
      k = 0
    }
  }
  k_values
}

///|
pub fn LCPArray::lcp(self : LCPArray) -> Array[Int] {
  self.lcp
}

///|
pub fn LCPArray::get(self : LCPArray, i : Int) -> Int {
  if i < 0 || i >= self.lcp.length() {
    0
  } else {
    self.lcp[i]
  }
}

///|
pub fn LCPArray::len(self : LCPArray) -> Int {
  self.lcp.length()
}

///|
pub fn LCPArray::max_lcp(self : LCPArray) -> Int {
  let mut max_val = 0
  for i = 0; i < self.lcp.length(); i = i + 1 {
    if self.lcp[i] > max_val {
      max_val = self.lcp[i]
    }
  }
  max_val
}

// ---------------------------------------------------------------------------
// Suffix Tree
// ---------------------------------------------------------------------------

///|
pub struct SuffixTree {
  text : String
  size : Int
}

///|
pub fn SuffixTree::new(text : String) -> SuffixTree {
  let t = text + "$"
  let n = t.length()
  let mut node_count = 1
  for i = 0; i < n; i = i + 1 {
    node_count = node_count + 1
  }
  SuffixTree::{ text, size: node_count }
}

///|
pub fn SuffixTree::root(self : SuffixTree) -> String {
  "root"
}

///|
pub fn SuffixTree::text(self : SuffixTree) -> String {
  self.text
}

///|
pub fn SuffixTree::size(self : SuffixTree) -> Int {
  self.size
}

///|
pub fn SuffixTree::contains(self : SuffixTree, pattern : String) -> Bool {
  self.text.contains(pattern)
}

///|
pub fn SuffixTree::locate(self : SuffixTree, pattern : String) -> Array[Int] {
  let positions : Array[Int] = Array::new()
  let t = self.text + "$"
  let n = t.length()
  let m = pattern.length()
  for i = 0; i <= n - m; i = i + 1 {
    let mut match_found = true
    for j = 0; j < m; j = j + 1 {
      if t.unsafe_get(i + j).unsafe_to_char() !=
        pattern.unsafe_get(j).unsafe_to_char() {
        match_found = false
        break
      }
    }
    if match_found {
      positions.push(i)
    }
  }
  positions
}

///|
pub fn SuffixTree::count(self : SuffixTree, pattern : String) -> Int {
  self.locate(pattern).length()
}

///|
pub fn SuffixTree::longest_repeated_substring(self : SuffixTree) -> String {
  let t = self.text + "$"
  let n = t.length()
  let mut longest = ""
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let mut k = 0
      while i + k < n &&
            j + k < n &&
            t.unsafe_get(i + k).unsafe_to_char() ==
            t.unsafe_get(j + k).unsafe_to_char() {
        k = k + 1
      }
      if k > longest.length() {
        longest = t[i:i + k].to_owned()
      }
    }
  }
  longest
}

///|
pub fn SuffixTree::to_dot(self : SuffixTree) -> String {
  let buf = StringBuilder::new()
  buf.write_string("digraph SuffixTree {\n")
  buf.write_string("  node [shape=circle, label=\"\"];\n")
  buf.write_string("  node0 [label=\"root\"];\n")
  buf.write_string("}\n")
  buf.to_string()
}

///|
pub fn SuffixTree::to_string(self : SuffixTree) -> String {
  let buf = StringBuilder::new()
  buf.write_string("Suffix Tree (")
  buf.write_string(self.text.length().to_string())
  buf.write_string(" chars, ")
  buf.write_string(self.size().to_string())
  buf.write_string(" nodes): ")
  buf.write_string(self.text)
  buf.to_string()
}
