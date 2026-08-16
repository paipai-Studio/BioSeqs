///|
fn shared_reference_test_sequence(
  name : String,
  aligned : String,
) -> @src.SharedReferenceSequence {
  @src.shared_reference_sequence(name, aligned) catch {
    _ => abort("valid shared-reference sequence should build")
  }
}

///|
fn shared_reference_test_input(
  aligned_reference : String,
  queries : Array[@src.SharedReferenceSequence],
) -> @src.SharedReferenceInput {
  @src.shared_reference_input(
    "ACGT",
    aligned_reference,
    queries,
    reference_name="reference",
  ) catch {
    _ => abort("valid shared-reference input should build")
  }
}

///|
fn shared_reference_test_sample() -> @src.SharedReferenceAlignment {
  @src.shared_reference_alignment_sample() catch {
    _ => abort("shared-reference sample should merge")
  }
}

///|
test "shared reference: query constructor derives ungapped sequence" {
  let query = shared_reference_test_sequence("read", "AC-GT")
  assert_eq(query.name, "read")
  assert_eq(query.sequence, "ACGT")
  assert_eq(query.aligned_sequence, "AC-GT")
  assert_eq(query.start, 0)
  assert_eq(query.end, 4)
}

///|
test "shared reference: query constructor preserves local coordinates" {
  let query = @src.shared_reference_sequence(
    "read",
    "AC-GT",
    sequence="NNACGTNN",
    start=2,
  ) catch {
    _ => abort("valid local query should build")
  }
  assert_eq(query.sequence, "NNACGTNN")
  assert_eq(query.start, 2)
  assert_eq(query.end, 6)
}

///|
test "shared reference: query constructor rejects negative start" {
  let failed = try {
    ignore(@src.shared_reference_sequence("read", "ACGT", start=-1))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: query constructor rejects sequence mismatch" {
  let failed = try {
    ignore(
      @src.shared_reference_sequence("read", "ACGT", sequence="AACCAA", start=1),
    )
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: query constructor rejects whitespace" {
  let failed = try {
    ignore(@src.shared_reference_sequence("read", "AC GT"))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: input validates and records reference span" {
  let query = shared_reference_test_sequence("read", "AC-GT")
  let input = shared_reference_test_input("AC-GT", [query])
  assert_eq(input.reference, "ACGT")
  assert_eq(input.reference_start, 0)
  assert_eq(input.reference_end, 4)
  assert_eq(input.queries.length(), 1)
}

///|
test "shared reference: input preserves reference metadata" {
  let query = shared_reference_test_sequence("read", "ACGT")
  let input = @src.shared_reference_input(
    "ACGT",
    "ACGT",
    [query],
    reference_name="plasmid",
    reference_description="circular construct",
  ) catch {
    _ => abort("valid input should build")
  }
  assert_eq(input.reference_name, "plasmid")
  assert_eq(input.reference_description, "circular construct")
}

///|
test "shared reference: input copies caller query array" {
  let queries = [shared_reference_test_sequence("read", "ACGT")]
  let input = shared_reference_test_input("ACGT", queries)
  queries.push(shared_reference_test_sequence("later", "ACGT"))
  assert_eq(input.queries.length(), 1)
}

///|
test "shared reference: input rejects empty reference" {
  let query = shared_reference_test_sequence("read", "")
  let failed = try {
    ignore(@src.shared_reference_input("", "", [query]))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: input rejects missing queries" {
  let failed = try {
    ignore(@src.shared_reference_input("ACGT", "ACGT", []))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: input rejects unequal row lengths" {
  let query = shared_reference_test_sequence("read", "ACGT")
  let failed = try {
    ignore(@src.shared_reference_input("ACGT", "AC-GT", [query]))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: input rejects aligned reference mismatch" {
  let query = shared_reference_test_sequence("read", "AGGT")
  let failed = try {
    ignore(@src.shared_reference_input("ACGT", "AGGT", [query]))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: official mixed PWA MSA example" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.num_sequences(), 4)
  assert_eq(alignment.num_queries(), 3)
  assert_eq(alignment.alignment_length(), 5)
  assert_eq(alignment.row(0), Some("ACG-T"))
  assert_eq(alignment.row(1), Some("AC--T"))
  assert_eq(alignment.row(2), Some("ACGGT"))
  assert_eq(alignment.row(3), Some("A---T"))
}

///|
test "shared reference: sample synchronizes insertion slots" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.insertion_widths, [0, 0, 0, 1, 0])
  assert_eq(alignment.aligned_reference, "ACG-T")
}

///|
test "shared reference: query ordering and metadata are stable" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.queries[0].name, "seq1")
  assert_eq(alignment.queries[1].name, "seq2")
  assert_eq(alignment.queries[2].name, "seq3")
  assert_eq(alignment.queries[1].description, "one insertion")
}

///|
test "shared reference: leading insertion is synchronized" {
  let first = shared_reference_test_input("-ACGT", [
    shared_reference_test_sequence("inserted", "TACGT"),
  ])
  let second = shared_reference_test_input("ACGT", [
    shared_reference_test_sequence("plain", "ACGT"),
  ])
  let alignment = @src.alignments_with_same_reference([first, second]) catch {
    _ => abort("leading insertion should merge")
  }
  assert_eq(alignment.aligned_reference, "-ACGT")
  assert_eq(alignment.queries[0].aligned_sequence, "TACGT")
  assert_eq(alignment.queries[1].aligned_sequence, "-ACGT")
}

///|
test "shared reference: trailing insertion is synchronized" {
  let first = shared_reference_test_input("ACGT--", [
    shared_reference_test_sequence("inserted", "ACGTAA"),
  ])
  let second = shared_reference_test_input("ACGT", [
    shared_reference_test_sequence("plain", "ACGT"),
  ])
  let alignment = @src.alignments_with_same_reference([first, second]) catch {
    _ => abort("trailing insertion should merge")
  }
  assert_eq(alignment.aligned_reference, "ACGT--")
  assert_eq(alignment.queries[1].aligned_sequence, "ACGT--")
}

///|
test "shared reference: widest insertion controls slot width" {
  let short = shared_reference_test_input("AC-GT", [
    shared_reference_test_sequence("short", "ACTGT"),
  ])
  let long = shared_reference_test_input("AC---GT", [
    shared_reference_test_sequence("long", "ACTTTGT"),
  ])
  let alignment = @src.alignments_with_same_reference([short, long]) catch {
    _ => abort("different insertion widths should merge")
  }
  assert_eq(alignment.aligned_reference, "AC---GT")
  assert_eq(alignment.queries[0].aligned_sequence, "ACT--GT")
  assert_eq(alignment.queries[1].aligned_sequence, "ACTTTGT")
}

///|
test "shared reference: insertions at distinct boundaries are retained" {
  let left = shared_reference_test_input("A-CGT", [
    shared_reference_test_sequence("left", "ATCGT"),
  ])
  let right = shared_reference_test_input("ACG-T", [
    shared_reference_test_sequence("right", "ACGGT"),
  ])
  let alignment = @src.alignments_with_same_reference([left, right]) catch {
    _ => abort("distinct insertion slots should merge")
  }
  assert_eq(alignment.aligned_reference, "A-CG-T")
  assert_eq(alignment.queries[0].aligned_sequence, "ATCG-T")
  assert_eq(alignment.queries[1].aligned_sequence, "A-CGGT")
}

///|
test "shared reference: multi-query input preserves internal gap columns" {
  let first = shared_reference_test_input("AC--GT", [
    shared_reference_test_sequence("one", "ACT-GT"),
    shared_reference_test_sequence("two", "AC-TGT"),
  ])
  let second = shared_reference_test_input("ACGT", [
    shared_reference_test_sequence("plain", "ACGT"),
  ])
  let alignment = @src.alignments_with_same_reference([first, second]) catch {
    _ => abort("multi-query input should merge")
  }
  assert_eq(alignment.queries[0].aligned_sequence, "ACT-GT")
  assert_eq(alignment.queries[1].aligned_sequence, "AC-TGT")
  assert_eq(alignment.queries[2].aligned_sequence, "AC--GT")
}

///|
test "shared reference: reference comparison is case insensitive" {
  let upper = shared_reference_test_input("ACGT", [
    shared_reference_test_sequence("upper", "ACGT"),
  ])
  let lower_query = shared_reference_test_sequence("lower", "acgt")
  let lower = @src.shared_reference_input("acgt", "acgt", [lower_query]) catch {
    _ => abort("lowercase input should build")
  }
  let alignment = @src.alignments_with_same_reference([upper, lower]) catch {
    _ => abort("case-insensitive references should merge")
  }
  assert_eq(alignment.num_queries(), 2)
}

///|
test "shared reference: merge rejects reference length mismatch" {
  let first = shared_reference_test_input("ACGT", [
    shared_reference_test_sequence("one", "ACGT"),
  ])
  let query = shared_reference_test_sequence("two", "ACGTA")
  let second = @src.shared_reference_input("ACGTA", "ACGTA", [query]) catch {
    _ => abort("second input should build")
  }
  let failed = try {
    ignore(@src.alignments_with_same_reference([first, second]))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: merge rejects reference content mismatch" {
  let first = shared_reference_test_input("ACGT", [
    shared_reference_test_sequence("one", "ACGT"),
  ])
  let query = shared_reference_test_sequence("two", "AGGT")
  let second = @src.shared_reference_input("AGGT", "AGGT", [query]) catch {
    _ => abort("second input should build")
  }
  let failed = try {
    ignore(@src.alignments_with_same_reference([first, second]))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: merge rejects inconsistent reference coordinates" {
  let query1 = shared_reference_test_sequence("one", "ACGT")
  let first = @src.shared_reference_input(
    "AACGTT",
    "ACGT",
    [query1],
    reference_start=1,
  ) catch {
    _ => abort("first local input should build")
  }
  let query2 = shared_reference_test_sequence("two", "AACG")
  let second = @src.shared_reference_input(
    "AACGTT",
    "AACG",
    [query2],
    reference_start=0,
  ) catch {
    _ => abort("second local input should build")
  }
  let failed = try {
    ignore(@src.alignments_with_same_reference([first, second]))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: merge rejects empty input list" {
  let failed = try {
    ignore(@src.alignments_with_same_reference([]))
    false
  } catch {
    SharedReferenceAlignmentError(_) => true
  }
  assert_true(failed)
}

///|
test "shared reference: PairwiseAlignment adapter preserves global rows" {
  let pairwise = @src.pairaligner_align("ACGT", "ACGGT")
  let input = @src.shared_reference_input_from_pairwise(
    pairwise,
    reference_name="ref",
    query_name="read",
  ) catch {
    _ => abort("global pairwise alignment should adapt")
  }
  assert_eq(input.reference, "ACGT")
  assert_eq(input.reference_name, "ref")
  assert_eq(input.queries[0].name, "read")
  assert_eq(input.reference_start, 0)
  assert_eq(input.reference_end, 4)
}

///|
test "shared reference: PairwiseAlignment adapter preserves local offsets" {
  let config = @src.PairwiseAlignerConfig::default_dna().set_mode(
    @src.pairaligner_local(),
  )
  let pairwise = @src.pairaligner_align("NNACGTNN", "TTACGTTT", config~)
  let input = @src.shared_reference_input_from_pairwise(
    pairwise,
    query_name="local",
  ) catch {
    _ => abort("local pairwise alignment should adapt")
  }
  assert_eq(input.reference_start, 2)
  assert_eq(input.reference_end, 6)
  assert_eq(input.queries[0].start, 2)
  assert_eq(input.queries[0].end, 6)
}

///|
test "shared reference: row access includes reference and bounds checks" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.row_name(0), Some("reference"))
  assert_eq(alignment.row_name(2), Some("seq2"))
  assert_eq(alignment.row(-1), None)
  assert_eq(alignment.row(4), None)
  assert_eq(alignment.row_name(4), None)
}

///|
test "shared reference: column access returns all rows" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.column(0), Some("AAAA"))
  assert_eq(alignment.column(2), Some("G-G-"))
  assert_eq(alignment.column(3), Some("--G-"))
  assert_eq(alignment.column(-1), None)
  assert_eq(alignment.column(5), None)
}

///|
test "shared reference: reference column mapping skips insertion columns" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.reference_to_column(0), Some(0))
  assert_eq(alignment.reference_to_column(2), Some(2))
  assert_eq(alignment.reference_to_column(3), Some(4))
  assert_eq(alignment.column_to_reference(3), None)
  assert_eq(alignment.column_to_reference(4), Some(3))
}

///|
test "shared reference: query column mapping tracks insertions" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.query_to_column(1, 2), Some(2))
  assert_eq(alignment.query_to_column(1, 3), Some(3))
  assert_eq(alignment.query_to_column(1, 4), Some(4))
  assert_eq(alignment.column_to_query(1, 3), Some(3))
  assert_eq(alignment.column_to_query(0, 3), None)
}

///|
test "shared reference: query to reference maps insertion to none" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.query_to_reference(1, 2), Some(2))
  assert_eq(alignment.query_to_reference(1, 3), None)
  assert_eq(alignment.query_to_reference(1, 4), Some(3))
}

///|
test "shared reference: reference to query maps deletion to none" {
  let alignment = shared_reference_test_sample()
  assert_eq(alignment.reference_to_query(0, 0), Some(0))
  assert_eq(alignment.reference_to_query(0, 2), None)
  assert_eq(alignment.reference_to_query(0, 3), Some(2))
}

///|
test "shared reference: counts classify deletion" {
  let alignment = shared_reference_test_sample()
  let counts = alignment.counts(0).unwrap()
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.insertions, 0)
  assert_eq(counts.deletions, 1)
  assert_eq(counts.identity_percent(), 100.0)
}

///|
test "shared reference: counts classify insertion" {
  let alignment = shared_reference_test_sample()
  let counts = alignment.counts(1).unwrap()
  assert_eq(counts.identities, 4)
  assert_eq(counts.insertions, 1)
  assert_eq(counts.deletions, 0)
  assert_eq(counts.aligned_pairs, 4)
  assert_eq(alignment.counts(-1), None)
}

///|
test "shared reference: find query returns metadata" {
  let alignment = shared_reference_test_sample()
  match alignment.find_query("seq2") {
    Some(query) => assert_eq(query.description, "one insertion")
    None => abort("seq2 should be present")
  }
  assert_eq(alignment.find_query("missing"), None)
}

///|
test "shared reference: converts to MultipleSeqAlignment" {
  let alignment = shared_reference_test_sample()
  let msa = alignment.to_msa() catch {
    _ => abort("valid shared-reference alignment should convert")
  }
  assert_eq(msa.num_records(), 4)
  assert_eq(msa.get_alignment_length(), 5)
  assert_eq(msa.records[0].id, "reference")
  assert_eq(msa.records[1].description, "one deletion")
  assert_eq(msa.records[2].seq.to_string(), "ACGGT")
}

///|
test "shared reference: aligned FASTA preserves metadata" {
  let alignment = shared_reference_test_sample()
  assert_eq(
    alignment.to_fasta(),
    ">reference shared reference\nACG-T\n" +
    ">seq1 one deletion\nAC--T\n" +
    ">seq2 one insertion\nACGGT\n" +
    ">seq3 two deletions\nA---T\n",
  )
}

///|
test "shared reference: summary reports dimensions and span" {
  let alignment = shared_reference_test_sample()
  assert_eq(
    alignment.summary(),
    "Shared-reference alignment: 4 sequences, 5 columns, reference 0:4",
  )
}

///|
test "shared reference: local alignment uses absolute coordinates" {
  let first_query = @src.shared_reference_sequence(
    "one",
    "AC-GT",
    sequence="NNACGTNN",
    start=2,
  ) catch {
    _ => abort("local query should build")
  }
  let first = @src.shared_reference_input(
    "NNAACGTNN",
    "AC-GT",
    [first_query],
    reference_start=3,
  ) catch {
    _ => abort("local input should build")
  }
  let second_query = @src.shared_reference_sequence(
    "two",
    "ACGGT",
    sequence="TTACGGTTT",
    start=2,
  ) catch {
    _ => abort("inserted local query should build")
  }
  let second = @src.shared_reference_input(
    "NNAACGTNN",
    "ACG-T",
    [second_query],
    reference_start=3,
  ) catch {
    _ => abort("second local input should build")
  }
  let alignment = @src.alignments_with_same_reference([first, second]) catch {
    _ => abort("local inputs should merge")
  }
  assert_eq(alignment.reference_start, 3)
  assert_eq(alignment.reference_end, 7)
  assert_eq(alignment.reference_to_column(6), Some(5))
  assert_eq(alignment.query_to_reference(1, 5), None)
  assert_eq(alignment.query_to_reference(1, 6), Some(6))
}

///|
test "shared reference: supports protein residues and mismatches" {
  let first_query = shared_reference_test_sequence("one", "MKTL")
  let first = @src.shared_reference_input("MKTL", "MKTL", [first_query]) catch {
    _ => abort("protein input should build")
  }
  let second_query = shared_reference_test_sequence("two", "MRATL")
  let second = @src.shared_reference_input("MKTL", "MK-TL", [second_query]) catch {
    _ => abort("inserted protein input should build")
  }
  let alignment = @src.alignments_with_same_reference([first, second]) catch {
    _ => abort("protein alignments should merge")
  }
  let counts = alignment.counts(1).unwrap()
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 1)
  assert_eq(counts.insertions, 1)
}
