///|
fn alignment_map_pair(
  target : String,
  query : String,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> @src.CoordinatePairwiseAlignment {
  @src.coordinate_pairwise_alignment(
    "target", target, "query", query, target_coordinates, query_coordinates,
  ) catch {
    _ => abort("valid coordinate alignment should build")
  }
}

///|
fn alignment_map_tutorial_outer() -> @src.CoordinatePairwiseAlignment {
  @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA",
    "transcript",
    "CCCCCCCGGGGGG",
    [8, 15, 26, 32],
    [0, 7, 7, 13],
  ) catch {
    _ => abort("valid tutorial outer alignment should build")
  }
}

///|
fn alignment_map_tutorial_read() -> @src.CoordinatePairwiseAlignment {
  @src.coordinate_pairwise_alignment(
    "transcript",
    "CCCCCCCGGGGGG",
    "read",
    "CCCCGGGG",
    [3, 11],
    [0, 8],
  ) catch {
    _ => abort("valid tutorial read alignment should build")
  }
}

///|
test "alignment map: official chromosome transcript read coordinates" {
  let mapped = alignment_map_tutorial_outer().map(alignment_map_tutorial_read()) catch {
    _ => abort("official map should succeed")
  }
  assert_eq(mapped.target_coordinates, [11, 15, 26, 30])
  assert_eq(mapped.query_coordinates, [0, 4, 4, 8])
  assert_eq(mapped.target_name, "chromosome")
  assert_eq(mapped.query_name, "read")
}

///|
test "alignment map: official example renders exon gap" {
  let mapped = @src.alignment_map_example() catch {
    _ => abort("example should map")
  }
  match mapped.aligned_rows() {
    Some((target, query)) => {
      assert_eq(target, "CCCCAAAAAAAAAAAGGGG")
      assert_eq(query, "CCCC-----------GGGG")
    }
    None => abort("concrete sequences should render")
  }
}

///|
test "alignment map: official example PSL" {
  let mapped = @src.alignment_map_example() catch {
    _ => abort("example should map")
  }
  let psl = mapped.to_psl() catch {
    _ => abort("non-empty map should serialize")
  }
  assert_eq(
    psl, "8\t0\t0\t0\t0\t0\t1\t11\t+\tread\t8\t0\t8\tchromosome\t40\t11\t30\t2\t4,4,\t0,4,\t11,26,\n",
  )
}

///|
test "alignment map: internal overlap" {
  let outer = @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA",
    "transcript",
    "GGGGGGGCCCCCGGGGGGA",
    [12, 31],
    [0, 19],
  ) catch {
    _ => abort("outer should build")
  }
  let inner = @src.coordinate_pairwise_alignment(
    "transcript",
    "GGGGGGGCCCCCGGGGGGA",
    "sequence",
    "GGCCCCCGGG",
    [5, 15],
    [0, 10],
  ) catch {
    _ => abort("inner should build")
  }
  let mapped = outer.map(inner) catch {
    _ => abort("internal map should succeed")
  }
  assert_eq(mapped.target_coordinates, [17, 27])
  assert_eq(mapped.query_coordinates, [0, 10])
}

///|
test "alignment map: clips left overhang" {
  let outer = @src.coordinate_pairwise_alignment(
    "chromosome",
    "GGGCCCCCGGGGGGAAAAAAAAAA",
    "transcript",
    "AGGGGGCCCCCGGGGGGA",
    [0, 15],
    [3, 18],
  ) catch {
    _ => abort("outer should build")
  }
  let inner = @src.coordinate_pairwise_alignment(
    "transcript",
    "AGGGGGCCCCCGGGGGGA",
    "sequence",
    "GGGGGCCCCCGGG",
    [1, 14],
    [0, 13],
  ) catch {
    _ => abort("inner should build")
  }
  let mapped = outer.map(inner) catch { _ => abort("left clipping should map") }
  assert_eq(mapped.target_coordinates, [0, 11])
  assert_eq(mapped.query_coordinates, [2, 13])
}

///|
test "alignment map: clips right overhang" {
  let outer = @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAAAAAAGGGGGGGCCCCCGGG",
    "transcript",
    "GGGGGGGCCCCCGGGGGGA",
    [12, 27],
    [0, 15],
  ) catch {
    _ => abort("outer should build")
  }
  let inner = @src.coordinate_pairwise_alignment(
    "transcript",
    "GGGGGGGCCCCCGGGGGGA",
    "sequence",
    "GGCCCCCGGGGG",
    [5, 17],
    [0, 12],
  ) catch {
    _ => abort("inner should build")
  }
  let mapped = outer.map(inner) catch {
    _ => abort("right clipping should map")
  }
  assert_eq(mapped.target_coordinates, [17, 27])
  assert_eq(mapped.query_coordinates, [0, 10])
}

///|
test "alignment map: reverse middle and reverse query yield forward result" {
  let outer = @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA",
    "transcript",
    "TCCCCCCGGGGGCCCCCCC",
    [12, 31],
    [19, 0],
  ) catch {
    _ => abort("reverse outer should build")
  }
  let inner = @src.coordinate_pairwise_alignment(
    "transcript",
    "TCCCCCCGGGGGCCCCCCC",
    "sequence",
    "GGCCCCCGGG",
    [4, 14],
    [10, 0],
  ) catch {
    _ => abort("reverse inner should build")
  }
  let mapped = outer.map(inner) catch {
    _ => abort("double reverse should map")
  }
  assert_eq(mapped.target_coordinates, [17, 27])
  assert_eq(mapped.query_coordinates, [0, 10])
  assert_false(mapped.is_reverse())
}

///|
test "alignment map: reverse final query is preserved" {
  let outer = @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA",
    "transcript",
    "GGGGGGGCCCCCGGGGGGA",
    [12, 31],
    [0, 19],
  ) catch {
    _ => abort("outer should build")
  }
  let inner = @src.coordinate_pairwise_alignment(
    "transcript",
    "GGGGGGGCCCCCGGGGGGA",
    "sequence",
    "CCCGGGGGCC",
    [5, 15],
    [10, 0],
  ) catch {
    _ => abort("reverse sequence should build")
  }
  let mapped = outer.map(inner) catch { _ => abort("reverse query should map") }
  assert_eq(mapped.target_coordinates, [17, 27])
  assert_eq(mapped.query_coordinates, [10, 0])
  assert_true(mapped.is_reverse())
}

///|
test "alignment map: reverse middle flips forward query" {
  let outer = @src.coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA",
    "transcript",
    "TCCCCCCGGGGGCCCCCCC",
    [12, 31],
    [19, 0],
  ) catch {
    _ => abort("reverse outer should build")
  }
  let inner = @src.coordinate_pairwise_alignment(
    "transcript",
    "TCCCCCCGGGGGCCCCCCC",
    "sequence",
    "CCCGGGGGCC",
    [4, 14],
    [0, 10],
  ) catch {
    _ => abort("inner should build")
  }
  let mapped = outer.map(inner) catch {
    _ => abort("reverse middle should map")
  }
  assert_eq(mapped.target_coordinates, [17, 27])
  assert_eq(mapped.query_coordinates, [10, 0])
}

///|
test "alignment map: coordinate blocks and counts" {
  let mapped = @src.alignment_map_example() catch {
    _ => abort("example should map")
  }
  let blocks = mapped.blocks()
  assert_eq(blocks.length(), 2)
  assert_eq(blocks[0].target_start, 11)
  assert_eq(blocks[0].query_start, 0)
  assert_eq(blocks[1].target_start, 26)
  let counts = mapped.counts()
  assert_eq(counts.aligned, 8)
  assert_eq(counts.query_gap_events, 1)
  assert_eq(counts.query_gap_bases, 11)
  assert_eq(counts.target_gap_events, 0)
}

///|
test "alignment map: target to query coordinate lookup" {
  let mapped = @src.alignment_map_example() catch {
    _ => abort("example should map")
  }
  assert_eq(mapped.target_to_query(11), Some(0))
  assert_eq(mapped.target_to_query(14), Some(3))
  assert_eq(mapped.target_to_query(20), None)
  assert_eq(mapped.target_to_query(29), Some(7))
}

///|
test "alignment map: query to target coordinate lookup" {
  let mapped = @src.alignment_map_example() catch {
    _ => abort("example should map")
  }
  assert_eq(mapped.query_to_target(0), Some(11))
  assert_eq(mapped.query_to_target(3), Some(14))
  assert_eq(mapped.query_to_target(4), Some(26))
  assert_eq(mapped.query_to_target(7), Some(29))
}

///|
test "alignment map: reverse coordinate lookup" {
  let alignment = alignment_map_pair("AACCGGTT", "AACCGGTT", [1, 5], [7, 3])
  assert_eq(alignment.target_to_query(1), Some(6))
  assert_eq(alignment.target_to_query(4), Some(3))
  assert_eq(alignment.query_to_target(6), Some(1))
  assert_eq(alignment.query_to_target(3), Some(4))
}

///|
test "alignment map: reverse rows use reverse complement" {
  let alignment = alignment_map_pair("AAAAAA", "AAGC", [1, 5], [4, 0])
  match alignment.aligned_rows() {
    Some((target, query)) => {
      assert_eq(target, "AAAA")
      assert_eq(query, "GCTT")
    }
    None => abort("rows should render")
  }
}

///|
test "alignment map: reverse PSL strand and query starts" {
  let alignment = alignment_map_pair("AAAAAA", "AAGC", [1, 5], [4, 0])
  let psl = alignment.to_psl() catch {
    _ => abort("reverse alignment should serialize")
  }
  assert_true(psl.contains("\t-\tquery\t4\t0\t4\ttarget\t6\t1\t5\t"))
  assert_true(psl.has_suffix("\t4,\t0,\t1,\n"))
}

///|
test "alignment map: summary reports direction" {
  let alignment = alignment_map_pair("AAAAAA", "AAGC", [1, 5], [4, 0])
  assert_eq(
    alignment.summary(),
    "CoordinatePairwiseAlignment(target <- query, blocks=1, aligned=4, strand=-)",
  )
}

///|
test "alignment map: undefined sequence contents still map" {
  let outer = @src.coordinate_pairwise_alignment_with_lengths(
    "genome",
    100,
    "transcript",
    20,
    [10, 20, 30, 40],
    [0, 10, 10, 20],
  ) catch {
    _ => abort("undefined outer should build")
  }
  let inner = @src.coordinate_pairwise_alignment_with_lengths(
    "transcript",
    20,
    "read",
    8,
    [6, 14],
    [0, 8],
  ) catch {
    _ => abort("undefined inner should build")
  }
  let mapped = outer.map(inner) catch {
    _ => abort("undefined alignments should map")
  }
  assert_eq(mapped.target_coordinates, [16, 20, 30, 34])
  assert_eq(mapped.aligned_rows(), None)
}

///|
test "alignment map: map many preserves order" {
  let outer = alignment_map_tutorial_outer()
  let first = alignment_map_tutorial_read()
  let second = @src.coordinate_pairwise_alignment(
    "transcript",
    "CCCCCCCGGGGGG",
    "other",
    "CCGG",
    [5, 9],
    [0, 4],
  ) catch {
    _ => abort("second read should build")
  }
  let mapped = outer.map_many([first, second]) catch {
    _ => abort("map many should succeed")
  }
  assert_eq(mapped.length(), 2)
  assert_eq(mapped[0].query_name, "read")
  assert_eq(mapped[1].query_name, "other")
}

///|
test "alignment map: no coordinate overlap returns empty" {
  let outer = alignment_map_pair("AAAAAAAAAA", "CCCCCCCCCC", [0, 4], [0, 4])
  let inner = alignment_map_pair("CCCCCCCCCC", "GG", [7, 9], [0, 2])
  let mapped = outer.map(inner) catch {
    _ => abort("disjoint paths should return an empty map")
  }
  assert_true(mapped.is_empty())
  assert_eq(mapped.blocks().length(), 0)
}

///|
test "alignment map: empty input remains empty" {
  let outer = @src.coordinate_pairwise_alignment_with_lengths(
    "outer",
    10,
    "middle",
    10,
    [],
    [],
  ) catch {
    _ => abort("empty coordinate map should build")
  }
  let inner = alignment_map_pair("AAAAAAAAAA", "AA", [0, 2], [0, 2])
  let mapped = outer.map(inner) catch { _ => abort("empty map should compose") }
  assert_true(mapped.is_empty())
}

///|
test "alignment map: constructor copies coordinate arrays" {
  let target_coordinates = [0, 4]
  let query_coordinates = [0, 4]
  let alignment = alignment_map_pair(
    "AAAA", "AAAA", target_coordinates, query_coordinates,
  )
  target_coordinates[1] = 2
  query_coordinates[1] = 3
  assert_eq(alignment.target_coordinates, [0, 4])
  assert_eq(alignment.query_coordinates, [0, 4])
}

///|
test "alignment map: adapts existing PairwiseAlignment" {
  let pairwise = @src.pairaligner_align("AACCGG", "ACG")
  let alignment = @src.coordinate_alignment_from_pairwise(
    pairwise,
    target_name="reference",
    query_name="read",
  ) catch {
    _ => abort("pairwise adapter should succeed")
  }
  assert_eq(
    alignment.target_coordinates.length(),
    pairwise.aligned_target.length() + 1,
  )
  assert_eq(alignment.target_coordinates[0], pairwise.target_start)
  assert_eq(
    alignment.target_coordinates[alignment.target_coordinates.length() - 1],
    pairwise.target_end,
  )
  assert_eq(alignment.target_name, "reference")
}

///|
test "alignment map: pairwise adapter preserves rendered gaps" {
  let pairwise = @src.pairaligner_align("AACCGG", "ACG")
  let alignment = @src.coordinate_alignment_from_pairwise(pairwise) catch {
    _ => abort("pairwise adapter should succeed")
  }
  match alignment.aligned_rows() {
    Some((target, query)) => {
      assert_eq(target, pairwise.aligned_target)
      assert_eq(query, pairwise.aligned_query)
    }
    None => abort("adapted concrete rows should render")
  }
}

///|
test "alignment map: rejects empty name" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment("", "AAAA", "query", "AAAA", [0, 4], [
        0, 4,
      ]),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects gapped raw sequence" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment(
        "target",
        "AA-A",
        "query",
        "AAAA",
        [0, 4],
        [0, 4],
      ),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects unequal coordinate arrays" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment(
        "target",
        "AAAA",
        "query",
        "AAAA",
        [0, 4],
        [0],
      ),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects one-point path" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment(
        "target",
        "AAAA",
        "query",
        "AAAA",
        [0],
        [0],
      ),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects coordinate out of bounds" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment(
        "target",
        "AAAA",
        "query",
        "AAAA",
        [0, 5],
        [0, 4],
      ),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects decreasing target coordinates" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment(
        "target",
        "AAAA",
        "query",
        "AAAA",
        [4, 0],
        [4, 0],
      ),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects query direction changes" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment(
        "target",
        "AAAA",
        "query",
        "AAAA",
        [0, 2, 4],
        [0, 2, 1],
      ),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects repeated coordinate point" {
  let failed = try {
    ignore(
      @src.coordinate_pairwise_alignment(
        "target",
        "AAAA",
        "query",
        "AAAA",
        [0, 0],
        [0, 0],
      ),
    )
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects incompatible middle lengths" {
  let outer = alignment_map_pair("AAAAAAAA", "CCCC", [0, 4], [0, 4])
  let inner = alignment_map_pair("CCCCCC", "GG", [0, 2], [0, 2])
  let failed = try {
    ignore(outer.map(inner))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: rejects unequal aligned steps" {
  let outer = alignment_map_pair("AAAAAAAA", "CCCC", [0, 4], [0, 2])
  let inner = alignment_map_pair("CCCC", "GGGG", [0, 4], [0, 4])
  let failed = try {
    ignore(outer.map(inner))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment map: empty alignment cannot serialize PSL" {
  let empty = @src.coordinate_pairwise_alignment_with_lengths(
    "target",
    4,
    "query",
    4,
    [],
    [],
  ) catch {
    _ => abort("empty alignment should build")
  }
  let failed = try {
    ignore(empty.to_psl())
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment mapall: builds coordinate MSA from gapped rows" {
  let alignment = @src.coordinate_multiple_alignment(["p1", "p2"], [
    "M-K", "MTK",
  ]) catch {
    _ => abort("valid MSA should build")
  }
  assert_eq(alignment.sequences, ["MK", "MTK"])
  assert_eq(alignment.coordinates[0], [0, 1, 1, 2])
  assert_eq(alignment.coordinates[1], [0, 1, 2, 3])
  assert_eq(alignment.alignment_length(), 3)
}

///|
test "alignment mapall: protein MSA projects to codon MSA" {
  let alignment = @src.coordinate_multiple_alignment(["p1", "p2"], [
    "M-K", "MTK",
  ]) catch {
    _ => abort("valid protein MSA should build")
  }
  let first = @src.coordinate_pairwise_alignment(
    "p1",
    "MK",
    "dna1",
    "ATGAAA",
    [0, 2],
    [0, 6],
  ) catch {
    _ => abort("valid codon mapping should build")
  }
  let second = @src.coordinate_pairwise_alignment(
    "p2",
    "MTK",
    "dna2",
    "ATGACCAAA",
    [0, 3],
    [0, 9],
  ) catch {
    _ => abort("valid codon mapping should build")
  }
  let mapped = alignment.mapall([first, second]) catch {
    _ => abort("codon mapall should succeed")
  }
  assert_eq(mapped.num_sequences(), 2)
  assert_eq(mapped.alignment_length(), 9)
  assert_eq(mapped.row(0), Some("ATG---AAA"))
  assert_eq(mapped.row(1), Some("ATGACCAAA"))
}

///|
test "alignment mapall: nucleotide one-to-one projection" {
  let alignment = @src.coordinate_multiple_alignment(["r1", "r2"], [
    "AC-G", "ACTG",
  ]) catch {
    _ => abort("valid nucleotide MSA should build")
  }
  let first = alignment_map_pair("ACG", "ACG", [0, 3], [0, 3])
  let second = alignment_map_pair("ACTG", "ACTG", [0, 4], [0, 4])
  let mapped = alignment.mapall([first, second]) catch {
    _ => abort("one-to-one mapall should succeed")
  }
  assert_eq(mapped.row(0), Some("AC-G"))
  assert_eq(mapped.row(1), Some("ACTG"))
}

///|
test "alignment mapall: FASTA serialization" {
  let alignment = @src.coordinate_multiple_alignment(["p1", "p2"], [
    "M-K", "MTK",
  ]) catch {
    _ => abort("valid MSA should build")
  }
  assert_eq(alignment.to_fasta(), ">p1\nM-K\n>p2\nMTK\n")
  assert_eq(
    alignment.summary(),
    "CoordinateMultipleAlignment(rows=2, columns=3)",
  )
}

///|
test "alignment mapall: rejects unequal MSA row width" {
  let failed = try {
    ignore(@src.coordinate_multiple_alignment(["p1", "p2"], ["MK", "M-K"]))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment mapall: rejects missing mapping" {
  let alignment = @src.coordinate_multiple_alignment(["p1", "p2"], [
    "M-K", "MTK",
  ]) catch {
    _ => abort("valid MSA should build")
  }
  let first = @src.coordinate_pairwise_alignment(
    "p1",
    "MK",
    "dna1",
    "ATGAAA",
    [0, 2],
    [0, 6],
  ) catch {
    _ => abort("valid mapping should build")
  }
  let failed = try {
    ignore(alignment.mapall([first]))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment mapall: rejects mapping target length mismatch" {
  let alignment = @src.coordinate_multiple_alignment(["p1"], ["MK"]) catch {
    _ => abort("valid MSA should build")
  }
  let mapping = @src.coordinate_pairwise_alignment(
    "p1",
    "MKK",
    "dna1",
    "ATGAAAAAA",
    [0, 3],
    [0, 9],
  ) catch {
    _ => abort("valid mapping should build")
  }
  let failed = try {
    ignore(alignment.mapall([mapping]))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment mapall: rejects inconsistent factors" {
  let alignment = @src.coordinate_multiple_alignment(["p1", "p2"], ["MK", "MK"]) catch {
    _ => abort("valid MSA should build")
  }
  let nucleotide = alignment_map_pair("MK", "MK", [0, 2], [0, 2])
  let codon = @src.coordinate_pairwise_alignment(
    "p2",
    "MK",
    "dna2",
    "ATGAAA",
    [0, 2],
    [0, 6],
  ) catch {
    _ => abort("valid codon mapping should build")
  }
  let failed = try {
    ignore(alignment.mapall([nucleotide, codon]))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment mapall: rejects unsupported two-to-one factor" {
  let alignment = @src.coordinate_multiple_alignment(["p1"], ["MK"]) catch {
    _ => abort("valid MSA should build")
  }
  let mapping = alignment_map_pair("MK", "ACGT", [0, 2], [0, 4])
  let failed = try {
    ignore(alignment.mapall([mapping]))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment mapall: rejects mapping without aligned block" {
  let alignment = @src.coordinate_multiple_alignment(["p1"], ["MK"]) catch {
    _ => abort("valid MSA should build")
  }
  let mapping = alignment_map_pair("MK", "AA", [0, 2, 2], [0, 0, 2])
  let failed = try {
    ignore(alignment.mapall([mapping]))
    false
  } catch {
    AlignmentMapError(_) => true
  }
  assert_true(failed)
}
