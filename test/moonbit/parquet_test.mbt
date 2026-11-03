///|
test "parquet construction and schema" {
  let t = @src.ParquetTable::new([
    ("gene", @src.PStr),
    ("fpkm", @src.PDouble),
    ("expressed", @src.PBool),
  ])
  assert_eq(t.num_row_groups(), 0)
  assert_eq(t.num_rows(), 0)
  assert_eq(t.pending_rows(), 0)
  let names = t.column_names()
  assert_eq(names, ["gene", "fpkm", "expressed"])
  assert_eq(t.schema().length(), 3)
}

///|
test "parquet duplicate column raises" {
  let msg = try {
    ignore(@src.ParquetTable::new([("a", @src.PInt32), ("a", @src.PDouble)]))
    "NO_ERROR"
  } catch {
    ParquetError(m) => m
  }
  assert_true(msg.contains("duplicate column"))
}

///|
test "parquet add row width mismatch raises" {
  let t = @src.ParquetTable::new([("a", @src.PInt32), ("b", @src.PDouble)])
  let msg = try {
    ignore(t.add_row([@src.PvInt(1)]))
    "NO_ERROR"
  } catch {
    ParquetError(m) => m
  }
  assert_true(msg.contains("row width"))
}

///|
test "parquet add row type mismatch raises" {
  let t = @src.ParquetTable::new([("a", @src.PInt32)])
  let msg = try {
    ignore(t.add_row([@src.PvDouble(1.0)]))
    "NO_ERROR"
  } catch {
    ParquetError(m) => m
  }
  assert_true(msg.contains("type mismatch"))
  ignore(t.add_row([@src.PvNull]))
  ignore(t.add_row([@src.PvInt(7)]))
  assert_eq(t.pending_rows(), 2)
}

///|
test "parquet flush and read column across groups" {
  let t = @src.ParquetTable::new([("a", @src.PInt32)])
  for i in 0..<3 {
    ignore(t.add_row([@src.PvInt(i)]))
  }
  assert_eq(t.flush_row_group(), 3)
  for i in 3..<5 {
    ignore(t.add_row([@src.PvInt(i)]))
  }
  assert_eq(t.flush_row_group(), 2)
  assert_eq(t.flush_row_group(), 0)
  assert_eq(t.num_row_groups(), 2)
  assert_eq(t.num_rows(), 5)
  let col = t.read_column("a")
  assert_eq(col.length(), 5)
  assert_eq(col[0], @src.PvInt(0))
  assert_eq(col[4], @src.PvInt(4))
  let msg = try {
    ignore(t.read_column("zzz"))
    "NO_ERROR"
  } catch {
    ParquetError(m) => m
  }
  assert_true(msg.contains("unknown column"))
}

///|
test "parquet to rows transpose" {
  let t = @src.ParquetTable::new([("g", @src.PStr), ("v", @src.PDouble)])
  ignore(t.add_row([@src.PvStr("x"), @src.PvDouble(1.5)]))
  ignore(t.add_row([@src.PvStr("y"), @src.PvDouble(2.5)]))
  ignore(t.flush_row_group())
  let rows = t.to_rows()
  assert_eq(rows.length(), 2)
  assert_eq(rows[0], [@src.PvStr("x"), @src.PvDouble(1.5)])
  assert_eq(rows[1], [@src.PvStr("y"), @src.PvDouble(2.5)])
}

///|
test "parquet chunk stats" {
  let t = @src.ParquetTable::new([("v", @src.PInt32)])
  ignore(t.add_row([@src.PvInt(5)]))
  ignore(t.add_row([@src.PvNull]))
  ignore(t.add_row([@src.PvInt(-2)]))
  ignore(t.add_row([@src.PvInt(9)]))
  ignore(t.flush_row_group())
  let st = t.chunk_stats(0, "v")
  assert_eq(st.min, Some(@src.PvInt(-2)))
  assert_eq(st.max, Some(@src.PvInt(9)))
  assert_eq(st.null_count, 1)
}

///|
test "parquet scan projection" {
  let t = @src.ParquetTable::new([
    ("g", @src.PStr),
    ("v", @src.PInt32),
    ("w", @src.PDouble),
  ])
  ignore(t.add_row([@src.PvStr("a"), @src.PvInt(1), @src.PvDouble(0.1)]))
  ignore(t.add_row([@src.PvStr("b"), @src.PvInt(2), @src.PvDouble(0.2)]))
  ignore(t.flush_row_group())
  let res = t.scan(columns=["g", "v"])
  assert_eq(res.groups_scanned, 1)
  assert_eq(res.groups_skipped, 0)
  assert_eq(res.rows.length(), 2)
  assert_eq(res.rows[0].length(), 2)
  assert_eq(res.rows[0][0], @src.PvStr("a"))
  assert_eq(res.rows[1][1], @src.PvInt(2))
}

///|
test "parquet scan numeric filter with pushdown skip" {
  let t = @src.ParquetTable::new([("v", @src.PInt32)])
  // 行组1: [1,3], 行组2: [10,30] — 谓词 v >= 10 应整组跳过行组1
  for i in 0..<3 {
    ignore(t.add_row([@src.PvInt(i + 1)]))
  }
  ignore(t.flush_row_group())
  for i in 0..<3 {
    ignore(t.add_row([@src.PvInt((i + 1) * 10)]))
  }
  ignore(t.flush_row_group())
  let res = t.scan(filters=[
    { column: "v", op: @src.OpGe, value: @src.PvInt(10), },
  ])
  assert_eq(res.groups_scanned, 1)
  assert_eq(res.groups_skipped, 1)
  assert_eq(res.rows.length(), 3)
  assert_eq(res.rows[0][0], @src.PvInt(10))
  assert_eq(res.rows[2][0], @src.PvInt(30))
}

///|
test "parquet scan string filter" {
  let t = @src.ParquetTable::new([("g", @src.PStr)])
  ignore(t.add_row([@src.PvStr("geneA")]))
  ignore(t.add_row([@src.PvStr("geneB")]))
  ignore(t.add_row([@src.PvStr("tfC")]))
  ignore(t.flush_row_group())
  let res = t.scan(filters=[
    { column: "g", op: @src.OpGe, value: @src.PvStr("gene"), },
  ])
  assert_eq(res.rows.length(), 2)
  assert_eq(res.rows[0][0], @src.PvStr("geneA"))
}

///|
test "parquet scan unknown column raises" {
  let t = @src.ParquetTable::new([("a", @src.PInt32)])
  ignore(t.add_row([@src.PvInt(1)]))
  ignore(t.flush_row_group())
  let msg = try {
    ignore(t.scan(columns=["nope"]))
    "NO_ERROR"
  } catch {
    ParquetError(m) => m
  }
  assert_true(msg.contains("unknown column"))
  let msg2 = try {
    ignore(
      t.scan(filters=[{ column: "nope", op: @src.OpEq, value: @src.PvInt(1), }]),
    )
    "NO_ERROR"
  } catch {
    ParquetError(m) => m
  }
  assert_true(msg2.contains("unknown filter column"))
}

///|
test "parquet scan null never matches" {
  let t = @src.ParquetTable::new([("v", @src.PInt32)])
  ignore(t.add_row([@src.PvNull]))
  ignore(t.add_row([@src.PvInt(5)]))
  ignore(t.flush_row_group())
  let res = t.scan(filters=[
    { column: "v", op: @src.OpGt, value: @src.PvInt(-1), },
  ])
  assert_eq(res.rows.length(), 1)
  assert_eq(res.rows[0][0], @src.PvInt(5))
}

///|
test "parquet scan empty table" {
  let t = @src.ParquetTable::new([("a", @src.PInt32)])
  let res = t.scan()
  assert_eq(res.rows.length(), 0)
  assert_eq(res.groups_scanned, 0)
  assert_eq(res.groups_skipped, 0)
}

///|
test "parquet pending rows not scanned" {
  let t = @src.ParquetTable::new([("a", @src.PInt32)])
  ignore(t.add_row([@src.PvInt(1)]))
  let res = t.scan()
  assert_eq(res.rows.length(), 0)
  assert_eq(t.num_rows(), 0)
  assert_eq(t.pending_rows(), 1)
}

///|
test "parquet dict encode" {
  let (dict, idx) = @src.parquet_dict_encode([
    "chr1", "chr2", "chr1", "chr1", "chr3", "chr2",
  ])
  assert_eq(dict, ["chr1", "chr2", "chr3"])
  assert_eq(idx, [0, 1, 0, 0, 2, 1])
}

///|
test "parquet json roundtrip" {
  let t = @src.ParquetTable::new([
    ("g", @src.PStr),
    ("v", @src.PInt32),
    ("l", @src.PInt64),
    ("d", @src.PDouble),
    ("b", @src.PBool),
  ])
  ignore(
    t.add_row([
      @src.PvStr("a"),
      @src.PvInt(1),
      @src.PvLong((42).to_int64()),
      @src.PvDouble(0.5),
      @src.PvBool(true),
    ]),
  )
  ignore(
    t.add_row([
      @src.PvNull,
      @src.PvInt(2),
      @src.PvLong((7).to_int64()),
      @src.PvDouble(1.5),
      @src.PvBool(false),
    ]),
  )
  ignore(t.flush_row_group())
  let back = @src.parquet_from_json(t.to_json()).unwrap()
  assert_eq(back.num_rows(), 2)
  assert_eq(back.num_row_groups(), 1)
  assert_eq(back.read_column("g")[0], @src.PvStr("a"))
  assert_true(back.read_column("g")[1] is @src.PvNull)
  assert_eq(back.read_column("v")[1], @src.PvInt(2))
  assert_eq(back.read_column("l")[0], @src.PvLong((42).to_int64()))
  assert_eq(back.read_column("d")[1], @src.PvDouble(1.5))
  assert_eq(back.read_column("b")[0], @src.PvBool(true))
  assert_eq(back.column_names(), ["g", "v", "l", "d", "b"])
}

///|
test "parquet json invalid returns none" {
  assert_true(@src.parquet_from_json("not json") is None)
  assert_true(@src.parquet_from_json("{\"row_groups\":[]}") is None)
  assert_true(
    @src.parquet_from_json("{\"schema\":[{\"name\":\"a\",\"type\":\"nope\"}]}")
    is None,
  )
}

///|
test "parquet filter combines pushdown and row filter" {
  // 行组1 v∈[1,2] padj∈[0.5,0.9]; 行组2 v∈[10,20] padj∈[0.01,0.04]
  // 谓词 v>=5 AND padj<0.05: 行组1 被 v 谓词整组跳过; 行组2 扫描后全过
  let t = @src.ParquetTable::new([("v", @src.PInt32), ("padj", @src.PDouble)])
  ignore(t.add_row([@src.PvInt(1), @src.PvDouble(0.5)]))
  ignore(t.add_row([@src.PvInt(2), @src.PvDouble(0.9)]))
  ignore(t.flush_row_group())
  ignore(t.add_row([@src.PvInt(10), @src.PvDouble(0.01)]))
  ignore(t.add_row([@src.PvInt(20), @src.PvDouble(0.04)]))
  ignore(t.flush_row_group())
  let res = t.scan(filters=[
    { column: "v", op: @src.OpGe, value: @src.PvInt(5), },
    { column: "padj", op: @src.OpLt, value: @src.PvDouble(0.05), },
  ])
  assert_eq(res.groups_skipped, 1)
  assert_eq(res.groups_scanned, 1)
  assert_eq(res.rows.length(), 2)
  assert_eq(res.rows[0][0], @src.PvInt(10))
}
