///|
/// Bio.Seq.Approximate module - Approximate string matching for sequences.
/// Inspired by Biopython's Bio.Seq.Approximate package.
/// Provides fuzzy search with substitutions, insertions, and deletions.

///|
/// Approximate match result.
pub struct ApproxMatch {
  pattern : String
  query : String
  start : Int
  end : Int
  mismatches : Int
  score : Int
  is_match : Bool
}

///|
/// Get the pattern string of the match.
pub fn ApproxMatch::get_pattern(self : ApproxMatch) -> String {
  self.pattern
}

///|
/// Get the query string of the match.
pub fn ApproxMatch::get_query(self : ApproxMatch) -> String {
  self.query
}

///|
/// Get the start position of the match (0-based).
pub fn ApproxMatch::get_start(self : ApproxMatch) -> Int {
  self.start
}

///|
/// Get the end position of the match (exclusive, 0-based).
pub fn ApproxMatch::get_end(self : ApproxMatch) -> Int {
  self.end
}

///|
/// Get the number of mismatches in the alignment.
pub fn ApproxMatch::get_mismatches(self : ApproxMatch) -> Int {
  self.mismatches
}

///|
/// Get the score of the match (positive = match, negative = mismatch).
pub fn ApproxMatch::get_score(self : ApproxMatch) -> Int {
  self.score
}

///|
/// Check if the match is valid.
pub fn ApproxMatch::is_match(self : ApproxMatch) -> Bool {
  self.is_match
}

///|
/// Get a string representation of the match.
pub fn ApproxMatch::to_string(self : ApproxMatch) -> String {
  let mut s = "ApproxMatch[".to_string()
  s = s + "pattern=\{self.pattern}, query=\{self.query}"
  s = s + ", start=\{self.start}, end=\{self.end}"
  s = s + ", mismatches=\{self.mismatches}"
  s = s + ", score=\{self.score}, valid=\{self.is_match}]"
  s
}

///|
/// Create a new ApproxMatch.
pub fn approx_match_new(
  pattern : String,
  query : String,
  start : Int,
  end : Int,
  mismatches : Int,
  score : Int,
  is_match : Bool,
) -> ApproxMatch {
  ApproxMatch::{ pattern, query, start, end, mismatches, score, is_match }
}

///|
/// Count mismatches between two equal-length strings.
pub fn count_mismatches(s1 : String, s2 : String) -> Int {
  let len = if s1.length() < s2.length() { s1.length() } else { s2.length() }
  let mut count = 0
  let mut i = 0
  while i < len {
    let c1 = s1.substring(start=i, end=i + 1).to_upper()
    let c2 = s2.substring(start=i, end=i + 1).to_upper()
    if c1 != c2 {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Approximate search with a maximum number of mismatches allowed.
/// Returns the first match found, or a non-matching ApproxMatch if not found.
pub fn approx_search(
  pattern : String,
  query : String,
  max_mismatches : Int,
) -> ApproxMatch {
  let pat_len = pattern.length()
  let qry_len = query.length()
  if pat_len == 0 || qry_len < pat_len {
    return ApproxMatch::{
      pattern,
      query,
      start: 0,
      end: 0,
      mismatches: pat_len,
      score: 0,
      is_match: false,
    }
  }

  let mut best : Option[ApproxMatch] = None
  let mut pos = 0
  let max_pos = qry_len - pat_len

  while pos <= max_pos {
    let substr = query.substring(start=pos, end=pos + pat_len)
    let mm = count_mismatches(pattern, substr)

    if mm <= max_mismatches {
      let score = pat_len - mm
      let m = ApproxMatch::{
        pattern,
        query,
        start: pos,
        end: pos + pat_len,
        mismatches: mm,
        score,
        is_match: true,
      }
      best = match best {
        Some(prev) => {
          if mm < prev.mismatches { Some(m) } else { Some(prev) }
        }
        None => Some(m)
      }
      if mm == 0 {
        return m
      }
    }
    pos = pos + 1
  }

  match best {
    Some(m) => m
    None => ApproxMatch::{
      pattern,
      query,
      start: 0,
      end: 0,
      mismatches: pat_len,
      score: 0,
      is_match: false,
    }
  }
}

///|
/// Approximate search allowing insertions and deletions (indels).
/// Uses a simple dynamic programming approach for edit distance.
pub fn approx_search_with_indels(
  pattern : String,
  query : String,
  max_errors : Int,
) -> ApproxMatch {
  let pat_len = pattern.length()
  let qry_len = query.length()

  if pat_len == 0 {
    return ApproxMatch::{
      pattern,
      query,
      start: 0,
      end: 0,
      mismatches: 0,
      score: 0,
      is_match: true,
    }
  }

  if qry_len == 0 {
    return ApproxMatch::{
      pattern,
      query,
      start: 0,
      end: 0,
      mismatches: pat_len,
      score: 0,
      is_match: false,
    }
  }

  let mut best_errors = max_errors + 1
  let mut best_start = 0
  let mut best_end = 0
  let mut found = false

  let mut start = 0
  while start < qry_len {
    // Allow segment to be shorter (deletions) or longer (insertions) than pattern
    let diff = if pat_len > max_errors { pat_len - max_errors } else { 0 }
    let min_end = start + diff
    let max_end_limit = start + pat_len + max_errors
    let max_end = if max_end_limit < qry_len { max_end_limit } else { qry_len }
    let init_end = if min_end > start { min_end } else { start + 1 }
    let mut end = init_end
    while end <= max_end {
      let seg = query.substring(start=start, end=end)
      let errors = levenshtein_distance(pattern, seg)
      if errors <= max_errors && errors < best_errors {
        best_errors = errors
        best_start = start
        best_end = end
        found = true
        if errors == 0 {
          return ApproxMatch::{
            pattern,
            query,
            start,
            end,
            mismatches: 0,
            score: pat_len,
            is_match: true,
          }
        }
      }
      end = end + 1
    }
    start = start + 1
  }

  if found {
    let seg = query.substring(start=best_start, end=best_end)
    let score = pat_len - best_errors
    ApproxMatch::{
      pattern,
      query,
      start: best_start,
      end: best_end,
      mismatches: best_errors,
      score,
      is_match: true,
    }
  } else {
    ApproxMatch::{
      pattern,
      query,
      start: 0,
      end: 0,
      mismatches: max_errors + 1,
      score: 0,
      is_match: false,
    }
  }
}

///|
/// Calculate Levenshtein edit distance between two strings.
/// Allowed operations: substitution, insertion, deletion (all cost 1).
pub fn levenshtein_distance(s1 : String, s2 : String) -> Int {
  let len1 = s1.length()
  let len2 = s2.length()

  if len1 == 0 { return len2 }
  if len2 == 0 { return len1 }

  let mut prev : Array[Int] = Array::new()
  let mut i = 0
  while i <= len2 {
    prev.push(i)
    i = i + 1
  }

  let mut row : Array[Int] = Array::new()
  i = 1
  while i <= len1 {
    row = Array::new()
    row.push(i)
    let c1 = s1.substring(start=i - 1, end=i).to_upper()

    let mut j = 1
    while j <= len2 {
      let c2 = s2.substring(start=j - 1, end=j).to_upper()
      let cost = if c1 == c2 { 0 } else { 1 }
      let deletion = prev[j] + 1
      let insertion = row[j - 1] + 1
      let substitution = prev[j - 1] + cost
      let min_val = if deletion < insertion { if deletion < substitution { deletion } else { substitution } } else { if insertion < substitution { insertion } else { substitution } }
      row.push(min_val)
      j = j + 1
    }

    prev = row
    i = i + 1
  }

  prev[len2]
}

///|
/// Find all approximate matches of a pattern in a query string.
pub fn approx_find_all(
  pattern : String,
  query : String,
  max_mismatches : Int,
) -> Array[ApproxMatch] {
  let pat_len = pattern.length()
  let qry_len = query.length()
  let results : Array[ApproxMatch] = Array::new()

  if pat_len == 0 || qry_len < pat_len {
    return results
  }

  let mut pos = 0
  let max_pos = qry_len - pat_len

  while pos <= max_pos {
    let substr = query.substring(start=pos, end=pos + pat_len)
    let mm = count_mismatches(pattern, substr)

    if mm <= max_mismatches {
      let score = pat_len - mm
      results.push(ApproxMatch::{
        pattern,
        query,
        start: pos,
        end: pos + pat_len,
        mismatches: mm,
        score,
        is_match: true,
      })
    }
    pos = pos + 1
  }

  results
}

///|
/// Approximate word search: split query into words and find approximate matches for each.
pub struct ApproxWordResult {
  word : String
  match_start : Int
  match_end : Int
  edit_distance : Int
  is_found : Bool
}

///|
/// Perform approximate word search for a single word in a text.
pub fn approx_word_search(
  word : String,
  text : String,
  max_errors : Int,
) -> ApproxWordResult {
  let len = word.length()
  let text_len = text.length()

  if len == 0 {
    return ApproxWordResult::{ word, match_start: 0, match_end: 0, edit_distance: 0, is_found: false }
  }

  let mut best_dist = max_errors + 1
  let mut best_start = 0
  let mut best_end = 0
  let mut found = false

  let mut start = 0
  while start < text_len {
    let max_end = if start + len + max_errors < text_len { start + len + max_errors } else { text_len }
    let mut end = start + len
    while end <= max_end {
      let seg = text.substring(start=start, end=end)
      let dist = levenshtein_distance(word, seg)
      if dist <= max_errors && dist < best_dist {
        best_dist = dist
        best_start = start
        best_end = end
        found = true
        if dist == 0 {
          return ApproxWordResult::{
            word,
            match_start: start,
            match_end: end,
            edit_distance: 0,
            is_found: true,
          }
        }
      }
      end = end + 1
    }
    start = start + 1
  }

  ApproxWordResult::{
    word,
    match_start: best_start,
    match_end: best_end,
    edit_distance: best_dist,
    is_found: found,
  }
}

///|
/// Find the best approximate match using a sliding window approach.
/// Optimized for longer sequences with a given error rate.
pub fn approx_best_match(
  pattern : String,
  query : String,
  max_error_rate : Double,
) -> ApproxMatch {
  let pat_len = pattern.length()
  let max_errors = (pat_len.to_double() * max_error_rate).to_int()

  approx_search_with_indels(pattern, query, max_errors)
}
