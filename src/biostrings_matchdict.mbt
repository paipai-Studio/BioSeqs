///|
/// PDict - Pattern dictionary for dictionary-based sequence matching.
///
/// Inspired by Bioconductor Biostrings::PDict.
/// Stores a collection of patterns that can be matched against a subject
/// sequence with support for mismatches and indels.
pub struct PDict {
  patterns : Array[String]
  max_mismatches : Int
  with_indels : Bool
}

///|
/// MatchPDictHit - A single hit from dictionary matching.
///
/// Records a successful alignment of a pattern from the PDict
/// against a region of the subject sequence.
pub struct MatchPDictHit {
  pattern_idx : Int
  pattern : String
  start : Int
  end : Int
  width : Int
  mismatches : Int
}

///|
/// MatchPDictResult - Result of dictionary matching against a subject.
///
/// Contains all hits found, the source PDict, and the subject length
/// for downstream analysis.
pub struct MatchPDictResult {
  hits : Array[MatchPDictHit]
  pdict : PDict
  subject_length : Int
}

///|
/// CountPatternResult - Count of pattern occurrences across subjects.
///
/// Stores per-subject counts and the total number of occurrences
/// found by vectorized pattern matching.
pub struct CountPatternResult {
  pattern : String
  counts : Array[Int]
  total : Int
}

///|
/// Create a pattern dictionary (PDict) from an array of patterns.
///
/// @param patterns Array of pattern strings to store in the dictionary.
/// @param max_mismatches Maximum number of mismatches allowed during matching.
///     Default is 0 (exact match).
/// @param with_indels Whether to allow insertions and deletions during matching.
///     Default is false.
/// @return A PDict structure ready for dictionary-based matching.
pub fn bmd_create_pdict(
  patterns~ : Array[String],
  max_mismatches? : Int = 0,
  with_indels? : Bool = false,
) -> PDict {
  PDict::{
    patterns,
    max_mismatches,
    with_indels,
  }
}

///|
/// Match all patterns in a PDict against a subject sequence.
///
/// Scans the subject sequence with every pattern in the dictionary,
/// returning all hits that satisfy the mismatch/indel constraints.
///
/// @param pdict The pattern dictionary containing patterns to search.
/// @param subject The subject sequence string to search against.
/// @return A MatchPDictResult containing all detected hits.
pub fn bmd_match_pdict(pdict~ : PDict, subject~ : String) -> MatchPDictResult {
  let hits : Array[MatchPDictHit] = Array::new()

  if subject.length() == 0 || pdict.patterns.length() == 0 {
    return MatchPDictResult::{
      hits,
      pdict,
      subject_length: subject.length(),
    }
  }

  let mut p_idx = 0
  while p_idx < pdict.patterns.length() {
    let pattern = pdict.patterns[p_idx]
    let plen = pattern.length()

    if plen == 0 || plen > subject.length() {
      p_idx = p_idx + 1
      continue
    }

    let mut i = 0
    while i <= subject.length() - plen {
      let mut mismatches = 0
      let mut j = 0
      while j < plen {
        let p_ch = pattern.unsafe_get(j)
        let s_ch = subject.unsafe_get(i + j)

        if p_ch != s_ch {
          mismatches = mismatches + 1
          if mismatches > pdict.max_mismatches {
            break
          }
        }
        j = j + 1
      }

      if j == plen && mismatches <= pdict.max_mismatches {
        hits.push(MatchPDictHit::{
          pattern_idx: p_idx,
          pattern,
          start: i + 1,
          end: i + plen,
          width: plen,
          mismatches,
        })
      }
      i = i + 1
    }
    p_idx = p_idx + 1
  }

  MatchPDictResult::{
    hits,
    pdict,
    subject_length: subject.length(),
  }
}

///|
/// Vectorized count of pattern occurrences across multiple subjects.
///
/// Counts how many times a pattern appears in each subject sequence,
/// allowing up to max_mismatches mismatches.
///
/// @param pattern The pattern string to search for.
/// @param subjects Array of subject sequences to search.
/// @param max_mismatches Maximum allowed mismatches (default 0).
/// @return A CountPatternResult with per-subject counts and total.
pub fn bmd_vcount_pattern(
  pattern~ : String,
  subjects~ : Array[String],
  max_mismatches? : Int = 0,
) -> CountPatternResult {
  let counts : Array[Int] = Array::new()
  let plen = pattern.length()
  let mut total = 0

  if plen == 0 {
    let mut i = 0
    while i < subjects.length() {
      counts.push(0)
      i = i + 1
    }
    return CountPatternResult::{
      pattern,
      counts,
      total: 0,
    }
  }

  let mut i = 0
  while i < subjects.length() {
    let subj = subjects[i]
    let mut count = 0

    if subj.length() >= plen {
      let mut j = 0
      while j <= subj.length() - plen {
        let mut mismatches = 0
        let mut k = 0
        while k < plen {
          let p_ch = pattern.unsafe_get(k)
          let s_ch = subj.unsafe_get(j + k)

          if p_ch != s_ch {
            mismatches = mismatches + 1
            if mismatches > max_mismatches {
              break
            }
          }
          k = k + 1
        }

        if k == plen && mismatches <= max_mismatches {
          count = count + 1
        }
        j = j + 1
      }
    }

    counts.push(count)
    total = total + count
    i = i + 1
  }

  CountPatternResult::{
    pattern,
    counts,
    total,
  }
}

///|
/// Vectorized matching of multiple patterns against multiple subjects.
///
/// Each pattern is counted across all subjects, returning an array
/// of CountPatternResult (one per pattern).
///
/// @param patterns Array of pattern strings to search for.
/// @param subjects Array of subject sequences to search.
/// @param max_mismatches Maximum allowed mismatches (default 0).
/// @return An array of CountPatternResult, one per input pattern.
pub fn bmd_vmatch_pattern(
  patterns~ : Array[String],
  subjects~ : Array[String],
  max_mismatches? : Int = 0,
) -> Array[CountPatternResult] {
  let results : Array[CountPatternResult] = Array::new()

  let mut i = 0
  while i < patterns.length() {
    let result = bmd_vcount_pattern(
      pattern=patterns[i],
      subjects=subjects,
      max_mismatches=max_mismatches,
    )
    results.push(result)
    i = i + 1
  }

  results
}

///|
/// Determine which subjects have at least one match to any pattern in the PDict.
///
/// @param pdict The pattern dictionary to search with.
/// @param subjects Array of subject sequences to test.
/// @return An array of Bool where true indicates the subject has at least one match.
pub fn bmd_which(pdict~ : PDict, subjects~ : Array[String]) -> Array[Bool] {
  let result : Array[Bool] = Array::new()

  let mut i = 0
  while i < subjects.length() {
    let subj = subjects[i]
    let mut found = false

    if subj.length() > 0 && pdict.patterns.length() > 0 {
      let mut p_idx = 0
      while p_idx < pdict.patterns.length() && !found {
        let pattern = pdict.patterns[p_idx]
        let plen = pattern.length()

        if plen > 0 && plen <= subj.length() {
          let mut j = 0
          while j <= subj.length() - plen && !found {
            let mut mismatches = 0
            let mut k = 0
            while k < plen {
              let p_ch = pattern.unsafe_get(k)
              let s_ch = subj.unsafe_get(j + k)

              if p_ch != s_ch {
                mismatches = mismatches + 1
                if mismatches > pdict.max_mismatches {
                  break
                }
              }
              k = k + 1
            }

            if k == plen && mismatches <= pdict.max_mismatches {
              found = true
            }
            j = j + 1
          }
        }
        p_idx = p_idx + 1
      }
    }

    result.push(found)
    i = i + 1
  }

  result
}

///|
/// Return the indices of subjects that have at least one match.
///
/// @param pdict The pattern dictionary to search with.
/// @param subjects Array of subject sequences to test.
/// @return An array of Int indices indicating matching subjects.
pub fn bmd_which_index(pdict~ : PDict, subjects~ : Array[String]) -> Array[Int] {
  let which_result = bmd_which(pdict=pdict, subjects=subjects)
  let indices : Array[Int] = Array::new()

  let mut i = 0
  while i < which_result.length() {
    if which_result[i] {
      indices.push(i)
    }
    i = i + 1
  }

  indices
}

///|
/// Count the total number of pattern occurrences across the subject.
///
/// Every occurrence of every pattern in the PDict is counted,
/// including overlapping matches.
///
/// @param pdict The pattern dictionary containing patterns to search.
/// @param subject The subject sequence string to search.
/// @return Total number of pattern occurrences found.
pub fn bmd_count_occurrences(pdict~ : PDict, subject~ : String) -> Int {
  let result = bmd_match_pdict(pdict=pdict, subject=subject)
  result.hits.length()
}

///|
/// Find the best (fewest mismatches) hit from dictionary matching.
///
/// Returns the single best MatchPDictHit across all patterns,
/// or None if no hits are found.
///
/// @param pdict The pattern dictionary to search with.
/// @param subject The subject sequence string to search.
/// @return The best MatchPDictHit, or None if no matches found.
pub fn bmd_find_best_match(pdict~ : PDict, subject~ : String) -> MatchPDictHit? {
  let result = bmd_match_pdict(pdict=pdict, subject=subject)

  if result.hits.length() == 0 {
    return None
  }

  let mut best = result.hits[0]
  let mut i = 1
  while i < result.hits.length() {
    if result.hits[i].mismatches < best.mismatches {
      best = result.hits[i]
    }
    i = i + 1
  }

  Some(best)
}

// ============================================================================
// Inline Tests
// ============================================================================

///|
test "bmd_create_pdict_basic" {
  let pdict = bmd_create_pdict(patterns=["ATG", "CCC", "GGG"])
  assert_eq(pdict.patterns.length(), 3)
  assert_eq(pdict.max_mismatches, 0)
  assert_eq(pdict.with_indels, false)
  assert_eq(pdict.patterns[0], "ATG")
  assert_eq(pdict.patterns[1], "CCC")
  assert_eq(pdict.patterns[2], "GGG")
}

///|
test "bmd_create_pdict_with_mismatches" {
  let pdict = bmd_create_pdict(patterns=["ATG", "TTT"], max_mismatches=2)
  assert_eq(pdict.max_mismatches, 2)
  assert_eq(pdict.patterns.length(), 2)
}

///|
test "bmd_create_pdict_empty" {
  let pdict = bmd_create_pdict(patterns=[])
  assert_eq(pdict.patterns.length(), 0)
}

///|
test "bmd_create_pdict_with_indels" {
  let pdict = bmd_create_pdict(patterns=["ATG"], max_mismatches=1, with_indels=true)
  assert_eq(pdict.with_indels, true)
  assert_eq(pdict.max_mismatches, 1)
}

///|
test "bmd_match_pdict_single_pattern" {
  let pdict = bmd_create_pdict(patterns=["ATG"])
  let result = bmd_match_pdict(pdict=pdict, subject="AATGCTAG")
  assert_eq(result.hits.length(), 1)
  assert_eq(result.hits[0].pattern, "ATG")
  assert_eq(result.hits[0].start, 2)
  assert_eq(result.hits[0].end, 4)
  assert_eq(result.subject_length, 8)
}

///|
test "bmd_match_pdict_multiple_patterns" {
  let pdict = bmd_create_pdict(patterns=["ATG", "GCT"])
  let result = bmd_match_pdict(pdict=pdict, subject="ATGGCT")
  assert_eq(result.hits.length(), 2)
  assert_eq(result.hits[0].pattern, "ATG")
  assert_eq(result.hits[1].pattern, "GCT")
}

///|
test "bmd_match_pdict_overlapping_hits" {
  let pdict = bmd_create_pdict(patterns=["AAA"])
  let result = bmd_match_pdict(pdict=pdict, subject="AAAA")
  assert_eq(result.hits.length(), 2)
}

///|
test "bmd_match_pdict_no_match" {
  let pdict = bmd_create_pdict(patterns=["XYZ"])
  let result = bmd_match_pdict(pdict=pdict, subject="ATGC")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_match_pdict_empty_subject" {
  let pdict = bmd_create_pdict(patterns=["ATG"])
  let result = bmd_match_pdict(pdict=pdict, subject="")
  assert_eq(result.hits.length(), 0)
  assert_eq(result.subject_length, 0)
}

///|
test "bmd_match_pdict_empty_patterns" {
  let pdict = bmd_create_pdict(patterns=[])
  let result = bmd_match_pdict(pdict=pdict, subject="ATGC")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_match_pdict_with_mismatches" {
  let pdict = bmd_create_pdict(patterns=["ATG"], max_mismatches=1)
  let result = bmd_match_pdict(pdict=pdict, subject="AXG")
  assert_eq(result.hits.length(), 1)
  assert_eq(result.hits[0].mismatches, 1)
}

///|
test "bmd_match_pdict_mismatch_too_many" {
  let pdict = bmd_create_pdict(patterns=["ATG"], max_mismatches=1)
  let result = bmd_match_pdict(pdict=pdict, subject="XYZ")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_match_pdict_pattern_longer_than_subject" {
  let pdict = bmd_create_pdict(patterns=["ATGCT"])
  let result = bmd_match_pdict(pdict=pdict, subject="AT")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_vcount_pattern_basic" {
  let subjects = ["ATGATG", "ATGCAT", "CCCCCC"]
  let result = bmd_vcount_pattern(pattern="ATG", subjects=subjects)
  assert_eq(result.pattern, "ATG")
  assert_eq(result.counts.length(), 3)
  assert_eq(result.counts[0], 2)
  assert_eq(result.counts[1], 1)
  assert_eq(result.counts[2], 0)
  assert_eq(result.total, 3)
}

///|
test "bmd_vcount_pattern_no_match" {
  let subjects = ["CCCC", "GGGG", "TTTT"]
  let result = bmd_vcount_pattern(pattern="AAAA", subjects=subjects)
  assert_eq(result.total, 0)
  assert_eq(result.counts[0], 0)
  assert_eq(result.counts[1], 0)
  assert_eq(result.counts[2], 0)
}

///|
test "bmd_vcount_pattern_with_mismatches" {
  let subjects = ["AXG", "AYG", "AZG", "ATG"]
  let result = bmd_vcount_pattern(pattern="ATG", subjects=subjects, max_mismatches=1)
  assert_eq(result.counts[0], 1)
  assert_eq(result.counts[1], 1)
  assert_eq(result.counts[2], 1)
  assert_eq(result.counts[3], 1)
  assert_eq(result.total, 4)
}

///|
test "bmd_vcount_pattern_empty_pattern" {
  let subjects = ["ATGC", "GCAT"]
  let result = bmd_vcount_pattern(pattern="", subjects=subjects)
  assert_eq(result.total, 0)
  assert_eq(result.counts.length(), 2)
}

///|
test "bmd_vcount_pattern_empty_subjects" {
  let subjects : Array[String] = []
  let result = bmd_vcount_pattern(pattern="ATG", subjects=subjects)
  assert_eq(result.total, 0)
  assert_eq(result.counts.length(), 0)
}

///|
test "bmd_vmatch_pattern_basic" {
  let patterns = ["ATG", "CCC"]
  let subjects = ["ATGATG", "CCCCCC"]
  let results = bmd_vmatch_pattern(patterns=patterns, subjects=subjects)
  assert_eq(results.length(), 2)
  assert_eq(results[0].pattern, "ATG")
  assert_eq(results[0].total, 2)
  assert_eq(results[1].pattern, "CCC")
  assert_eq(results[1].total, 4)
}

///|
test "bmd_vmatch_pattern_multiple_subjects" {
  let patterns = ["ATG"]
  let subjects = ["ATGATG", "ATXATG", "TTTTTT"]
  let results = bmd_vmatch_pattern(patterns=patterns, subjects=subjects)
  assert_eq(results.length(), 1)
  assert_eq(results[0].counts[0], 2)
  assert_eq(results[0].counts[2], 0)
}

///|
test "bmd_vmatch_pattern_empty" {
  let results = bmd_vmatch_pattern(patterns=[], subjects=["ATGC"])
  assert_eq(results.length(), 0)
}

///|
test "bmd_which_basic" {
  let pdict = bmd_create_pdict(patterns=["ATG"])
  let subjects = ["ATGCTA", "TTTTTT", "ATGCAT"]
  let result = bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result.length(), 3)
  assert_eq(result[0], true)
  assert_eq(result[1], false)
  assert_eq(result[2], true)
}

///|
test "bmd_which_no_match" {
  let pdict = bmd_create_pdict(patterns=["ZZZZZ"])
  let subjects = ["ATGC", "GCAT"]
  let result = bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result[0], false)
  assert_eq(result[1], false)
}

///|
test "bmd_which_empty_pdict" {
  let pdict = bmd_create_pdict(patterns=[])
  let subjects = ["ATGC", "GCAT"]
  let result = bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result[0], false)
  assert_eq(result[1], false)
}

///|
test "bmd_which_index_basic" {
  let pdict = bmd_create_pdict(patterns=["ATG"])
  let subjects = ["TTTT", "ATGC", "GGGG", "ATGG"]
  let indices = bmd_which_index(pdict=pdict, subjects=subjects)
  assert_eq(indices.length(), 2)
  assert_eq(indices[0], 1)
  assert_eq(indices[1], 3)
}

///|
test "bmd_which_index_no_match" {
  let pdict = bmd_create_pdict(patterns=["ZZZZ"])
  let subjects = ["ATGC", "GCAT"]
  let indices = bmd_which_index(pdict=pdict, subjects=subjects)
  assert_eq(indices.length(), 0)
}

///|
test "bmd_count_occurrences_basic" {
  let pdict = bmd_create_pdict(patterns=["ATG"])
  let count = bmd_count_occurrences(pdict=pdict, subject="ATGATGATG")
  assert_eq(count, 3)
}

///|
test "bmd_count_occurrences_multiple_patterns" {
  let pdict = bmd_create_pdict(patterns=["AT", "TG"])
  let count = bmd_count_occurrences(pdict=pdict, subject="ATGATG")
  assert_eq(count, 4)
}

///|
test "bmd_count_occurrences_no_match" {
  let pdict = bmd_create_pdict(patterns=["ZZZZ"])
  let count = bmd_count_occurrences(pdict=pdict, subject="ATGC")
  assert_eq(count, 0)
}

///|
test "bmd_count_occurrences_empty" {
  let pdict = bmd_create_pdict(patterns=[""])
  let count = bmd_count_occurrences(pdict=pdict, subject="ATGC")
  assert_eq(count, 0)
}

///|
test "bmd_find_best_match_exact" {
  let pdict = bmd_create_pdict(patterns=["ATG"])
  let best = bmd_find_best_match(pdict=pdict, subject="AATGCT")
  assert_true(best is Some(_))
  assert_eq(best.unwrap().pattern, "ATG")
  assert_eq(best.unwrap().start, 2)
}

///|
test "bmd_find_best_match_with_mismatches" {
  let pdict = bmd_create_pdict(patterns=["ATG", "AXG", "AYG"], max_mismatches=1)
  let best = bmd_find_best_match(pdict=pdict, subject="AAAG")
  assert_true(best is Some(_))
  assert_eq(best.unwrap().mismatches, 1)
  assert_eq(best.unwrap().pattern, "ATG")
}

///|
test "bmd_find_best_match_no_hit" {
  let pdict = bmd_create_pdict(patterns=["ZZZZZ"])
  let best = bmd_find_best_match(pdict=pdict, subject="ATGC")
  assert_true(best is None)
}

///|
test "bmd_find_best_match_multiple_hits" {
  let pdict = bmd_create_pdict(patterns=["ATG"], max_mismatches=2)
  let best = bmd_find_best_match(pdict=pdict, subject="AATGCTXG")
  assert_true(best is Some(_))
  assert_eq(best.unwrap().pattern, "ATG")
  assert_eq(best.unwrap().mismatches, 0)
}

///|
test "bmd_which_with_mismatches" {
  let pdict = bmd_create_pdict(patterns=["ATG"], max_mismatches=1)
  let subjects = ["AXGCTA", "TTTTTT", "AYGCAT"]
  let result = bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result[0], true)
  assert_eq(result[1], false)
  assert_eq(result[2], true)
}

///|
test "bmd_vcount_pattern_single_subject" {
  let subjects = ["AAAA"]
  let result = bmd_vcount_pattern(pattern="AA", subjects=subjects)
  assert_eq(result.counts.length(), 1)
  assert_eq(result.counts[0], 3)
  assert_eq(result.total, 3)
}

///|
test "bmd_match_pdict_hit_properties" {
  let pdict = bmd_create_pdict(patterns=["ATG", "TGC"])
  let result = bmd_match_pdict(pdict=pdict, subject="ATGC")
  assert_eq(result.hits.length(), 2)
  assert_eq(result.hits[0].pattern_idx, 0)
  assert_eq(result.hits[0].width, 3)
  assert_eq(result.hits[0].mismatches, 0)
  assert_eq(result.hits[1].pattern_idx, 1)
  assert_eq(result.hits[1].start, 2)
  assert_eq(result.hits[1].end, 4)
}

///|
test "bmd_count_occurrences_overlapping" {
  let pdict = bmd_create_pdict(patterns=["AA"])
  let count = bmd_count_occurrences(pdict=pdict, subject="AAAA")
  assert_eq(count, 3)
}

///|
test "bmd_vmatch_pattern_multiple_mismatches" {
  let patterns = ["ATG", "CCC"]
  let subjects = ["AXGATG", "CXCCCX"]
  let results = bmd_vmatch_pattern(patterns=patterns, subjects=subjects, max_mismatches=1)
  assert_eq(results.length(), 2)
  assert_eq(results[0].counts[0], 2)
  assert_eq(results[1].counts[1], 4)
}
