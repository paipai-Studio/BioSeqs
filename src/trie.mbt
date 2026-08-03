///|
/// Bio.trie - Trie (prefix tree) data structure for efficient prefix matching.
/// Inspired by Biopython's Bio.trie and Bio.triefind modules.
///
/// A trie stores key-value pairs where keys are strings, enabling:
///   - O(m) lookup (m = key length)
///   - Prefix-based search (find all keys with a given prefix)
///   - Longest-prefix matching (match beginning of a string)
///   - Substring search (find keys anywhere in a string)
///
/// Common bioinformatics use cases:
///   - Restriction enzyme site detection
///   - Primer matching
///   - k-mer lookup
///   - Motif scanning

///|
/// A node in the trie.
pub struct TrieNode {
  mut children : Map[String, TrieNode]
  mut value : TrieValue
}

///|
/// Value stored in a trie node (None = not a key endpoint, Some(v) = key with value v).
pub enum TrieValue {
  Empty
  Value(String)
}

///|
/// Helper constructors for cross-package access.
pub fn trie_value_empty() -> TrieValue {
  TrieValue::Empty
}

///|
pub fn trie_value_of(v : String) -> TrieValue {
  TrieValue::Value(v)
}

///|
/// The trie data structure.
pub struct Trie {
  root : TrieNode
  mut count : Int
}

///|
/// Create a new empty trie.
pub fn Trie::new() -> Trie {
  Trie::{
    root: TrieNode::{
      children: Map([], capacity=4),
      value: TrieValue::Empty,
    },
    count: 0,
  }
}

///|
/// Create a new TrieNode.
pub fn TrieNode::new() -> TrieNode {
  TrieNode::{
    children: Map([], capacity=4),
    value: TrieValue::Empty,
  }
}

///|
/// Insert a key-value pair into the trie.
/// If the key already exists, its value is updated.
pub fn trie_insert(trie : Trie, key : String, value : String) -> Unit {
  let node = trie.root
  let mut current = node
  let mut i = 0
  while i < key.length() {
    let ch = key[i:i + 1].to_owned()
    let child = current.children.get(ch)
    match child {
      Some(c) => current = c
      None => {
        let new_node = TrieNode::new()
        current.children[ch] = new_node
        current = new_node
      }
    }
    i = i + 1
  }
  let was_empty = match current.value {
    TrieValue::Empty => true
    _ => false
  }
  current.value = TrieValue::Value(value)
  if was_empty {
    trie.count = trie.count + 1
  }
}

///|
/// Lookup the value for a key. Returns None if not found.
pub fn trie_lookup(trie : Trie, key : String) -> String? {
  let node = trie_find_node(trie.root, key)
  match node {
    Some(n) =>
      match n.value {
        TrieValue::Value(v) => Some(v)
        TrieValue::Empty => None
      }
    None => None
  }
}

///|
/// Check if a key exists in the trie.
pub fn trie_has_key(trie : Trie, key : String) -> Bool {
  let node = trie_find_node(trie.root, key)
  match node {
    Some(n) =>
      match n.value {
        TrieValue::Empty => false
        _ => true
      }
    None => false
  }
}

///|
/// Check if any key in the trie has the given prefix.
pub fn trie_has_prefix(trie : Trie, prefix : String) -> Bool {
  let node = trie_find_node(trie.root, prefix)
  node is Some(_)
}

///|
/// Get the number of keys in the trie.
pub fn trie_size(trie : Trie) -> Int {
  trie.count
}

///|
/// Find all keys with the given prefix.
/// Returns an array of (key, value) tuples.
pub fn trie_with_prefix(trie : Trie, prefix : String) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array::new()
  let start_node = trie_find_node(trie.root, prefix)
  match start_node {
    Some(n) => trie_collect_pairs(n, prefix, result)
    None => ()
  }
  result
}

///|
/// Get all keys in the trie.
pub fn trie_keys(trie : Trie) -> Array[String] {
  let result : Array[String] = Array::new()
  trie_collect_keys(trie.root, "", result)
  result
}

///|
/// Get all (key, value) pairs in the trie.
pub fn trie_items(trie : Trie) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array::new()
  trie_collect_pairs(trie.root, "", result)
  result
}

///|
/// Remove a key from the trie. Returns true if the key was found and removed.
pub fn trie_remove(trie : Trie, key : String) -> Bool {
  let node = trie_find_node(trie.root, key)
  match node {
    Some(n) =>
      match n.value {
        TrieValue::Empty => false
        TrieValue::Value(_) => {
          n.value = TrieValue::Empty
          trie.count = trie.count - 1
          true
        }
      }
    None => false
  }
}

///|
/// Clear all keys from the trie.
pub fn trie_clear(trie : Trie) -> Unit {
  trie.root.children = Map([], capacity=4)
  trie.root.value = TrieValue::Empty
  trie.count = 0
}

// ============================================================================
// Bio.triefind - Finding occurrences of trie keys in strings
// ============================================================================

///|
/// Find the longest key in the trie that matches the beginning of the string.
/// Returns the key and its value, or None if no match.
pub fn triefind_match(trie : Trie, text : String) -> (String, String)? {
  let mut current = trie.root
  let mut best_key : String? = None
  let mut best_val : String? = None
  let mut prefix = ""
  let mut i = 0
  while i < text.length() {
    let ch = text[i:i + 1].to_owned()
    let child = current.children.get(ch)
    match child {
      Some(c) => {
        current = c
        prefix = prefix + ch
        match c.value {
          TrieValue::Value(v) => {
            best_key = Some(prefix)
            best_val = Some(v)
          }
          TrieValue::Empty => ()
        }
      }
      None => break
    }
    i = i + 1
  }
  match (best_key, best_val) {
    (Some(k), Some(v)) => Some((k, v))
    _ => None
  }
}

///|
/// Find all keys in the trie that match the beginning of the string.
/// Returns a list of (key, value) tuples, ordered by increasing key length.
pub fn triefind_match_all(trie : Trie, text : String) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array::new()
  let mut current = trie.root
  let mut prefix = ""
  let mut i = 0
  while i < text.length() {
    let ch = text[i:i + 1].to_owned()
    let child = current.children.get(ch)
    match child {
      Some(c) => {
        current = c
        prefix = prefix + ch
        match c.value {
          TrieValue::Value(v) => result.push((prefix, v))
          TrieValue::Empty => ()
        }
      }
      None => break
    }
    i = i + 1
  }
  result
}

///|
/// Find all keys in the trie that match anywhere in the string.
/// Returns a list of (key, value, start, end) tuples.
pub fn triefind_find(trie : Trie, text : String) -> Array[(String, String, Int, Int)] {
  let result : Array[(String, String, Int, Int)] = Array::new()
  let mut start = 0
  while start < text.length() {
    let matches = triefind_match_all(trie, text[start:text.length()].to_owned())
    let mut j = 0
    while j < matches.length() {
      let pair = matches[j]
      let key = pair.0
      let val = pair.1
      let end = start + key.length()
      result.push((key, val, start, end))
      j = j + 1
    }
    start = start + 1
  }
  result
}

///|
/// Find all keys in the trie that match whole words in the string.
/// Word boundaries are defined as any punctuation, whitespace, or string start/end.
/// Returns a list of (key, value, start, end) tuples.
pub fn triefind_find_words(trie : Trie, text : String) -> Array[(String, String, Int, Int)] {
  let result : Array[(String, String, Int, Int)] = Array::new()
  let mut start = 0
  while start < text.length() {
    // Check if this position is at a word boundary
    if start == 0 || trie_is_boundary(text[start - 1:start].to_owned()) {
      let matches = triefind_match_all(trie, text[start:text.length()].to_owned())
      let mut j = 0
      while j < matches.length() {
        let pair = matches[j]
        let key = pair.0
        let val = pair.1
        let end = start + key.length()
        // Check that the character after the match is also a boundary
        if end >= text.length() || trie_is_boundary(text[end:end + 1].to_owned()) {
          result.push((key, val, start, end))
        }
        j = j + 1
      }
    }
    start = start + 1
  }
  result
}

// ============================================================================
// Internal helpers
// ============================================================================

///|
/// Find the node corresponding to a key (or prefix). Returns None if not found.
fn trie_find_node(root : TrieNode, key : String) -> TrieNode? {
  let mut current = root
  let mut i = 0
  while i < key.length() {
    let ch = key[i:i + 1].to_owned()
    let child = current.children.get(ch)
    match child {
      Some(c) => current = c
      None => return None
    }
    i = i + 1
  }
  Some(current)
}

///|
/// Recursively collect all keys from a node.
fn trie_collect_keys(node : TrieNode, prefix : String, result : Array[String]) -> Unit {
  match node.value {
    TrieValue::Value(_) => result.push(prefix)
    TrieValue::Empty => ()
  }
  let keys = node.children.keys()
  let key_arr : Array[String] = keys.collect()
  let mut i = 0
  while i < key_arr.length() {
    let ch = key_arr[i]
    let child = node.children[ch]
    trie_collect_keys(child, prefix + ch, result)
    i = i + 1
  }
}

///|
/// Recursively collect all (key, value) pairs from a node.
fn trie_collect_pairs(
  node : TrieNode,
  prefix : String,
  result : Array[(String, String)]
) -> Unit {
  match node.value {
    TrieValue::Value(v) => result.push((prefix, v))
    TrieValue::Empty => ()
  }
  let keys = node.children.keys()
  let key_arr : Array[String] = keys.collect()
  let mut i = 0
  while i < key_arr.length() {
    let ch = key_arr[i]
    let child = node.children[ch]
    trie_collect_pairs(child, prefix + ch, result)
    i = i + 1
  }
}

///|
/// Check if a character is a word boundary (whitespace or punctuation).
fn trie_is_boundary(ch : String) -> Bool {
  match ch {
    " " => true
    "\t" => true
    "\n" => true
    "\r" => true
    "." => true
    "," => true
    ";" => true
    ":" => true
    "!" => true
    "?" => true
    "(" => true
    ")" => true
    "[" => true
    "]" => true
    "{" => true
    "}" => true
    "\"" => true
    "'" => true
    "/" => true
    "\\" => true
    "|" => true
    "-" => true
    "_" => true
    _ => false
  }
}

///|
/// Create a sample trie with restriction enzyme recognition sites.
pub fn trie_sample_restriction_sites() -> Trie {
  let trie = Trie::new()
  trie_insert(trie, "GAATTC", "EcoRI")
  trie_insert(trie, "GGATCC", "BamHI")
  trie_insert(trie, "AAGCTT", "HindIII")
  trie_insert(trie, "CTGCAG", "PstI")
  trie_insert(trie, "GTCGAC", "SalI")
  trie_insert(trie, "CCCGGG", "SmaI")
  trie_insert(trie, "GCGGCCGC", "NotI")
  trie
}

///|
/// Create a sample trie with DNA primers.
pub fn trie_sample_primers() -> Trie {
  let trie = Trie::new()
  trie_insert(trie, "ATCGATCGAT", "Primer_Fwd_1")
  trie_insert(trie, "ATCGATCGAC", "Primer_Fwd_2")
  trie_insert(trie, "GCTAGCTAGC", "Primer_Rev_1")
  trie_insert(trie, "TTTTGGGGCC", "Primer_Rev_2")
  trie
}
