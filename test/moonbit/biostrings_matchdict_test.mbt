///|
test "bmd_create_pdict_basic" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG", "CCC", "GGG"])
  assert_eq(pdict.patterns.length(), 3)
  assert_eq(pdict.max_mismatches, 0)
  assert_eq(pdict.with_indels, false)
  assert_eq(pdict.patterns[0], "ATG")
  assert_eq(pdict.patterns[1], "CCC")
  assert_eq(pdict.patterns[2], "GGG")
}

///|
test "bmd_create_pdict_with_mismatches" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG", "TTT"], max_mismatches=2)
  assert_eq(pdict.max_mismatches, 2)
  assert_eq(pdict.patterns.length(), 2)
}

///|
test "bmd_create_pdict_empty" {
  let pdict = @src.bmd_create_pdict(patterns=[])
  assert_eq(pdict.patterns.length(), 0)
}

///|
test "bmd_create_pdict_with_indels" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"], max_mismatches=1, with_indels=true)
  assert_eq(pdict.with_indels, true)
  assert_eq(pdict.max_mismatches, 1)
}

///|
test "bmd_match_pdict_single_pattern" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="AATGCTAG")
  assert_eq(result.hits.length(), 1)
  assert_eq(result.hits[0].pattern, "ATG")
  assert_eq(result.hits[0].start, 2)
  assert_eq(result.hits[0].end, 4)
  assert_eq(result.subject_length, 8)
}

///|
test "bmd_match_pdict_multiple_patterns" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG", "GCT"])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="ATGGCT")
  assert_eq(result.hits.length(), 2)
  assert_eq(result.hits[0].pattern, "ATG")
  assert_eq(result.hits[1].pattern, "GCT")
}

///|
test "bmd_match_pdict_overlapping_hits" {
  let pdict = @src.bmd_create_pdict(patterns=["AAA"])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="AAAA")
  assert_eq(result.hits.length(), 2)
}

///|
test "bmd_match_pdict_no_match" {
  let pdict = @src.bmd_create_pdict(patterns=["XYZ"])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="ATGC")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_match_pdict_empty_subject" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="")
  assert_eq(result.hits.length(), 0)
  assert_eq(result.subject_length, 0)
}

///|
test "bmd_match_pdict_empty_patterns" {
  let pdict = @src.bmd_create_pdict(patterns=[])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="ATGC")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_match_pdict_with_mismatches" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"], max_mismatches=1)
  let result = @src.bmd_match_pdict(pdict=pdict, subject="AXG")
  assert_eq(result.hits.length(), 1)
  assert_eq(result.hits[0].mismatches, 1)
}

///|
test "bmd_match_pdict_mismatch_too_many" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"], max_mismatches=1)
  let result = @src.bmd_match_pdict(pdict=pdict, subject="XYZ")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_match_pdict_pattern_longer_than_subject" {
  let pdict = @src.bmd_create_pdict(patterns=["ATGCT"])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="AT")
  assert_eq(result.hits.length(), 0)
}

///|
test "bmd_vcount_pattern_basic" {
  let subjects = ["ATGATG", "ATGCAT", "CCCCCC"]
  let result = @src.bmd_vcount_pattern(pattern="ATG", subjects=subjects)
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
  let result = @src.bmd_vcount_pattern(pattern="AAAA", subjects=subjects)
  assert_eq(result.total, 0)
  assert_eq(result.counts[0], 0)
  assert_eq(result.counts[1], 0)
  assert_eq(result.counts[2], 0)
}

///|
test "bmd_vcount_pattern_with_mismatches" {
  let subjects = ["AXG", "AYG", "AZG", "ATG"]
  let result = @src.bmd_vcount_pattern(pattern="ATG", subjects=subjects, max_mismatches=1)
  assert_eq(result.counts[0], 1)
  assert_eq(result.counts[1], 1)
  assert_eq(result.counts[2], 1)
  assert_eq(result.counts[3], 1)
  assert_eq(result.total, 4)
}

///|
test "bmd_vcount_pattern_empty_pattern" {
  let subjects = ["ATGC", "GCAT"]
  let result = @src.bmd_vcount_pattern(pattern="", subjects=subjects)
  assert_eq(result.total, 0)
  assert_eq(result.counts.length(), 2)
}

///|
test "bmd_vcount_pattern_empty_subjects" {
  let subjects : Array[String] = []
  let result = @src.bmd_vcount_pattern(pattern="ATG", subjects=subjects)
  assert_eq(result.total, 0)
  assert_eq(result.counts.length(), 0)
}

///|
test "bmd_vmatch_pattern_basic" {
  let patterns = ["ATG", "CCC"]
  let subjects = ["ATGATG", "CCCCCC"]
  let results = @src.bmd_vmatch_pattern(patterns=patterns, subjects=subjects)
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
  let results = @src.bmd_vmatch_pattern(patterns=patterns, subjects=subjects)
  assert_eq(results.length(), 1)
  assert_eq(results[0].counts[0], 2)
  assert_eq(results[0].counts[2], 0)
}

///|
test "bmd_vmatch_pattern_empty" {
  let results = @src.bmd_vmatch_pattern(patterns=[], subjects=["ATGC"])
  assert_eq(results.length(), 0)
}

///|
test "bmd_which_basic" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"])
  let subjects = ["ATGCTA", "TTTTTT", "ATGCAT"]
  let result = @src.bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result.length(), 3)
  assert_eq(result[0], true)
  assert_eq(result[1], false)
  assert_eq(result[2], true)
}

///|
test "bmd_which_no_match" {
  let pdict = @src.bmd_create_pdict(patterns=["ZZZZZ"])
  let subjects = ["ATGC", "GCAT"]
  let result = @src.bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result[0], false)
  assert_eq(result[1], false)
}

///|
test "bmd_which_empty_pdict" {
  let pdict = @src.bmd_create_pdict(patterns=[])
  let subjects = ["ATGC", "GCAT"]
  let result = @src.bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result[0], false)
  assert_eq(result[1], false)
}

///|
test "bmd_which_index_basic" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"])
  let subjects = ["TTTT", "ATGC", "GGGG", "ATGG"]
  let indices = @src.bmd_which_index(pdict=pdict, subjects=subjects)
  assert_eq(indices.length(), 2)
  assert_eq(indices[0], 1)
  assert_eq(indices[1], 3)
}

///|
test "bmd_which_index_no_match" {
  let pdict = @src.bmd_create_pdict(patterns=["ZZZZ"])
  let subjects = ["ATGC", "GCAT"]
  let indices = @src.bmd_which_index(pdict=pdict, subjects=subjects)
  assert_eq(indices.length(), 0)
}

///|
test "bmd_count_occurrences_basic" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"])
  let count = @src.bmd_count_occurrences(pdict=pdict, subject="ATGATGATG")
  assert_eq(count, 3)
}

///|
test "bmd_count_occurrences_multiple_patterns" {
  let pdict = @src.bmd_create_pdict(patterns=["AT", "TG"])
  let count = @src.bmd_count_occurrences(pdict=pdict, subject="ATGATG")
  assert_eq(count, 4)
}

///|
test "bmd_count_occurrences_no_match" {
  let pdict = @src.bmd_create_pdict(patterns=["ZZZZ"])
  let count = @src.bmd_count_occurrences(pdict=pdict, subject="ATGC")
  assert_eq(count, 0)
}

///|
test "bmd_count_occurrences_empty" {
  let pdict = @src.bmd_create_pdict(patterns=[""])
  let count = @src.bmd_count_occurrences(pdict=pdict, subject="ATGC")
  assert_eq(count, 0)
}

///|
test "bmd_find_best_match_exact" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"])
  let best = @src.bmd_find_best_match(pdict=pdict, subject="AATGCT")
  assert_true(best is Some(_))
  assert_eq(best.unwrap().pattern, "ATG")
  assert_eq(best.unwrap().start, 2)
}

///|
test "bmd_find_best_match_with_mismatches" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG", "AXG", "AYG"], max_mismatches=1)
  let best = @src.bmd_find_best_match(pdict=pdict, subject="AAAG")
  assert_true(best is Some(_))
  assert_eq(best.unwrap().mismatches, 1)
  assert_eq(best.unwrap().pattern, "ATG")
}

///|
test "bmd_find_best_match_no_hit" {
  let pdict = @src.bmd_create_pdict(patterns=["ZZZZZ"])
  let best = @src.bmd_find_best_match(pdict=pdict, subject="ATGC")
  assert_true(best is None)
}

///|
test "bmd_find_best_match_multiple_hits" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"], max_mismatches=2)
  let best = @src.bmd_find_best_match(pdict=pdict, subject="AATGCTXG")
  assert_true(best is Some(_))
  assert_eq(best.unwrap().pattern, "ATG")
  assert_eq(best.unwrap().mismatches, 0)
}

///|
test "bmd_which_with_mismatches" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG"], max_mismatches=1)
  let subjects = ["AXGCTA", "TTTTTT", "AYGCAT"]
  let result = @src.bmd_which(pdict=pdict, subjects=subjects)
  assert_eq(result[0], true)
  assert_eq(result[1], false)
  assert_eq(result[2], true)
}

///|
test "bmd_vcount_pattern_single_subject" {
  let subjects = ["AAAA"]
  let result = @src.bmd_vcount_pattern(pattern="AA", subjects=subjects)
  assert_eq(result.counts.length(), 1)
  assert_eq(result.counts[0], 3)
  assert_eq(result.total, 3)
}

///|
test "bmd_match_pdict_hit_properties" {
  let pdict = @src.bmd_create_pdict(patterns=["ATG", "TGC"])
  let result = @src.bmd_match_pdict(pdict=pdict, subject="ATGC")
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
  let pdict = @src.bmd_create_pdict(patterns=["AA"])
  let count = @src.bmd_count_occurrences(pdict=pdict, subject="AAAA")
  assert_eq(count, 3)
}

///|
test "bmd_vmatch_pattern_multiple_mismatches" {
  let patterns = ["ATG", "CCC"]
  let subjects = ["AXGATG", "CXCCCX"]
  let results = @src.bmd_vmatch_pattern(patterns=patterns, subjects=subjects, max_mismatches=1)
  assert_eq(results.length(), 2)
  assert_eq(results[0].counts[0], 2)
  assert_eq(results[1].counts[1], 4)
}
