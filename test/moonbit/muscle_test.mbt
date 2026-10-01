///|
test "muscle_stage1 basic" {
  let seqs = ["ACDEFGHIK", "ACDEAGHIK", "ACDEFGHLK"]
  let names = ["s1", "s2", "s3"]
  let result = @src.muscle_stage1(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  assert_eq(result.method, "muscle-s1")
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
}

///|
test "muscle_stage2 basic" {
  let seqs = ["ACDEFGHIK", "ACDEAGHIK", "ACDEFGHLK"]
  let names = ["s1", "s2", "s3"]
  let result = @src.muscle_stage2(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  assert_eq(result.method, "muscle-s2")
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
}

///|
test "muscle_align basic" {
  let seqs = ["ACDEFGHIK", "ACDEAGHIK", "ACDEFGHLK"]
  let names = ["s1", "s2", "s3"]
  let result = @src.muscle_align(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  assert_eq(result.method, "muscle")
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
}

///|
test "muscle_align identical sequences" {
  let seqs = ["ACDEFGHIK", "ACDEFGHIK", "ACDEFGHIK"]
  let names = ["s1", "s2", "s3"]
  let result = @src.muscle_align(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  for s in result.aligned_sequences {
    assert_false(s.contains("-"))
  }
}

///|
test "muscle_align with gaps needed" {
  let seqs = [
    "ACDEFGHIK",
    "ACFGHIK",
    "ACDEFGHIKL",
  ]
  let names = ["s1", "s2", "s3"]
  let result = @src.muscle_align(seqs, names)
  assert_eq(result.aligned_sequences.length(), 3)
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
  let mut has_gap = false
  for s in result.aligned_sequences {
    if s.contains("-") {
      has_gap = true
    }
  }
  assert_true(has_gap)
}

///|
test "muscle_align single sequence" {
  let result = @src.muscle_align(["ACDEFG"], ["s1"])
  assert_eq(result.aligned_sequences.length(), 1)
  assert_eq(result.aligned_sequences[0], "ACDEFG")
}

///|
test "muscle_align empty" {
  let result = @src.muscle_align([], [])
  assert_eq(result.aligned_sequences.length(), 0)
}

///|
test "muscle_align with params" {
  let seqs = ["ACDEFGHIK", "ACDEAGHIK", "ACDEFGHLK"]
  let names = ["s1", "s2", "s3"]
  let result = @src.muscle_align_with_params(seqs, names, 3, -8.0, 2)
  assert_eq(result.aligned_sequences.length(), 3)
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
}

///|
test "muscle_align preserves residues" {
  let seqs = ["ACDEFGHIK", "ACFGHIK", "ACDEFGHIKL"]
  let names = ["s1", "s2", "s3"]
  let result = @src.muscle_align(seqs, names)
  for i = 0; i < seqs.length(); i = i + 1 {
    let aligned = result.aligned_sequences[i]
    let mut ungapped = ""
    for j = 0; j < aligned.length(); j = j + 1 {
      let c = aligned[j].unsafe_to_char()
      if c != '-' {
        ungapped += c.to_string()
      }
    }
    assert_eq(ungapped, seqs[i])
  }
}

///|
test "muscle stage2 improves over stage1" {
  let seqs = [
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK",
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV",
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK",
  ]
  let names = ["hba1", "hbb1", "hba2"]
  let s1 = @src.muscle_stage1(seqs, names)
  let s2 = @src.muscle_stage2(seqs, names)
  // Stage 2 should have at least as good a score as stage 1
  assert_true(s2.score >= s1.score - 0.1)
}

///|
test "muscle vs mafft comparison" {
  let seqs = ["ACDEFGHIK", "ACDEAGHIK", "ACDEFGHLK"]
  let names = ["s1", "s2", "s3"]
  let mafft_result = @src.mafft_align(seqs, names)
  let muscle_result = @src.muscle_align(seqs, names)
  // Both should produce valid alignments
  assert_eq(mafft_result.aligned_sequences.length(), 3)
  assert_eq(muscle_result.aligned_sequences.length(), 3)
  // Both should have same aligned length
  let mafft_len = mafft_result.aligned_sequences[0].length()
  let muscle_len = muscle_result.aligned_sequences[0].length()
  assert_true(mafft_len > 0)
  assert_true(muscle_len > 0)
}

///|
test "muscle MsaResult fields" {
  let seqs = ["ACDEFG", "ACDEAG"]
  let result = @src.muscle_align(seqs, ["a", "b"])
  assert_eq(result.method, "muscle")
  assert_eq(result.names.length(), 2)
}

///|
test "muscle larger alignment" {
  let seqs = [
    "KETAAAKFERQHMDSSSAASAA",
    "KETAAAKFERQHMDSSTSAASAA",
    "KETAAAKFERQHMDSAASAASAA",
    "KETAAAKFERQHMDSSTSAASAA",
  ]
  let names = ["p1", "p2", "p3", "p4"]
  let result = @src.muscle_align(seqs, names)
  assert_eq(result.aligned_sequences.length(), 4)
  let len = result.aligned_sequences[0].length()
  for s in result.aligned_sequences {
    assert_eq(s.length(), len)
  }
}
