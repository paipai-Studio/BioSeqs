///|
/// Tests for BIOM 1.0 parse + write round-trip.
test "biom parse dense" {
  // 最小稠密 BIOM 1.0：2 observations × 3 samples
  let json =
    #|{
    #|  "id": "test_table",
    #|  "format": "Biological Observation Matrix 1.0.0",
    #|  "matrix_type": "dense",
    #|  "matrix_element_type": "float",
    #|  "shape": [2, 3],
    #|  "rows": [
    #|    {"id": "OTU1", "metadata": {"taxonomy": ["Bacteria", "Firmicutes"]}},
    #|    {"id": "OTU2", "metadata": {"taxonomy": ["Bacteria", "Bacteroidetes"]}}
    #|  ],
    #|  "columns": [
    #|    {"id": "S1", "metadata": {"body_site": "gut"}},
    #|    {"id": "S2", "metadata": {"body_site": "skin"}},
    #|    {"id": "S3", "metadata": {}}
    #|  ],
    #|  "data": [[5.0, 0.0, 3.0], [0.0, 7.0, 1.0]]
    #|}

  let t = @src.parse_biom(json)
  assert_true(t is Some(_))

  let table = t.unwrap()
  assert_eq(table.table_id, "test_table")
  assert_eq(table.matrix_type, "dense")
  assert_eq(table.observation_ids.length(), 2)
  assert_eq(table.sample_ids.length(), 3)
  assert_eq(table.observation_ids[0], "OTU1")
  assert_eq(table.observation_ids[1], "OTU2")
  assert_eq(table.sample_ids[0], "S1")
  assert_eq(table.sample_ids[1], "S2")
  assert_eq(table.sample_ids[2], "S3")
  // data: [[5,0,3],[0,7,1]]
  assert_eq(table.data[0][0], 5.0)
  assert_eq(table.data[0][1], 0.0)
  assert_eq(table.data[0][2], 3.0)
  assert_eq(table.data[1][0], 0.0)
  assert_eq(table.data[1][1], 7.0)
  assert_eq(table.data[1][2], 1.0)
}

test "biom parse sparse" {
  // 稀疏格式：同上数据但用 [row, col, value] 三元组
  let json =
    #|{
    #|  "id": "sparse_table",
    #|  "format": "Biological Observation Matrix 1.0.0",
    #|  "matrix_type": "sparse",
    #|  "matrix_element_type": "float",
    #|  "shape": [2, 3],
    #|  "rows": [
    #|    {"id": "OTU1", "metadata": null},
    #|    {"id": "OTU2", "metadata": null}
    #|  ],
    #|  "columns": [
    #|    {"id": "S1", "metadata": null},
    #|    {"id": "S2", "metadata": null},
    #|    {"id": "S3", "metadata": null}
    #|  ],
    #|  "data": [[0, 0, 5.0], [0, 2, 3.0], [1, 1, 7.0], [1, 2, 1.0]]
    #|}

  let t = @src.parse_biom(json)
  assert_true(t is Some(_))

  let table = t.unwrap()
  assert_eq(table.matrix_type, "sparse")
  assert_eq(table.observation_ids.length(), 2)
  assert_eq(table.sample_ids.length(), 3)
  // 稀疏应展开为同样的稠密矩阵
  assert_eq(table.data[0][0], 5.0)
  assert_eq(table.data[0][1], 0.0)
  assert_eq(table.data[0][2], 3.0)
  assert_eq(table.data[1][0], 0.0)
  assert_eq(table.data[1][1], 7.0)
  assert_eq(table.data[1][2], 1.0)
}

test "biom roundtrip dense" {
  let original =
    #|{
    #|  "id": "rt_dense",
    #|  "format": "Biological Observation Matrix 1.0.0",
    #|  "matrix_type": "dense",
    #|  "matrix_element_type": "float",
    #|  "shape": [2, 2],
    #|  "rows": [
    #|    {"id": "O1", "metadata": {"phylum": "Firmicutes"}},
    #|    {"id": "O2", "metadata": {"phylum": "Bacteroidetes"}}
    #|  ],
    #|  "columns": [
    #|    {"id": "C1", "metadata": {"site": "gut"}},
    #|    {"id": "C2", "metadata": {"site": "oral"}}
    #|  ],
    #|  "data": [[1.0, 2.0], [3.0, 4.0]]
    #|}

  let t1 = @src.parse_biom(original).unwrap()
  let json_out = t1.to_json()
  let t2 = @src.parse_biom(json_out).unwrap()

  assert_eq(t2.table_id, t1.table_id)
  assert_eq(t2.matrix_type, t1.matrix_type)
  assert_eq(t2.observation_ids, t1.observation_ids)
  assert_eq(t2.sample_ids, t1.sample_ids)
  assert_eq(t2.data.length(), t1.data.length())
  assert_eq(t2.data[0][0], 1.0)
  assert_eq(t2.data[0][1], 2.0)
  assert_eq(t2.data[1][0], 3.0)
  assert_eq(t2.data[1][1], 4.0)
}

test "biom roundtrip sparse" {
  // parse sparse → to_json (should emit sparse) → parse again
  let original =
    #|{
    #|  "id": "rt_sparse",
    #|  "format": "Biological Observation Matrix 1.0.0",
    #|  "matrix_type": "sparse",
    #|  "matrix_element_type": "float",
    #|  "shape": [3, 2],
    #|  "rows": [
    #|    {"id": "A", "metadata": null},
    #|    {"id": "B", "metadata": null},
    #|    {"id": "C", "metadata": null}
    #|  ],
    #|  "columns": [
    #|    {"id": "X", "metadata": null},
    #|    {"id": "Y", "metadata": null}
    #|  ],
    #|  "data": [[0, 1, 10.0], [2, 0, 20.0]]
    #|}

  let t1 = @src.parse_biom(original).unwrap()
  let json_out = t1.to_json()
  let t2 = @src.parse_biom(json_out).unwrap()

  assert_eq(t2.matrix_type, "sparse")
  assert_eq(t2.observation_ids.length(), 3)
  assert_eq(t2.sample_ids.length(), 2)
  assert_eq(t2.data[0][1], 10.0)
  assert_eq(t2.data[2][0], 20.0)
  assert_eq(t2.data[0][0], 0.0)
  assert_eq(t2.data[1][0], 0.0)
}

test "biom empty metadata" {
  // metadata 为 null 或缺失时应正常解析
  let json =
    #|{
    #|  "id": "no_meta",
    #|  "format": "Biological Observation Matrix 1.0.0",
    #|  "matrix_type": "dense",
    #|  "matrix_element_type": "float",
    #|  "shape": [1, 1],
    #|  "rows": [{"id": "O1", "metadata": null}],
    #|  "columns": [{"id": "C1", "metadata": null}],
    #|  "data": [[42.0]]
    #|}

  let t = @src.parse_biom(json).unwrap()
  assert_eq(t.observation_ids[0], "O1")
  assert_eq(t.sample_ids[0], "C1")
  assert_eq(t.data[0][0], 42.0)
}

test "json parser basic" {
  // 验证通用 JSON 解析器
  let v = @src.parse_json("{\"a\": 1, \"b\": [true, false, null], \"c\": \"hello\"}")
  assert_true(v is Some(_))
  let j = v.unwrap()
  assert_true(j is @src.JsonObj(_))

  let s = @src.json_to_string(j)
  // 重新解析应一致
  let v2 = @src.parse_json(s)
  assert_true(v2 is Some(_))
}

// ============================================================
// Tests for BiomTable operations (new + parse + TSV + merge)
// ============================================================

test "biom new and accessors" {
  let table = @src.BiomTable::new(
    ["OTU1", "OTU2"],
    ["S1", "S2", "S3"],
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    table_id="test",
  )
  assert_eq(table.table_id, "test")
  assert_eq(table.n_observations(), 2)
  assert_eq(table.n_samples(), 3)
  assert_eq(table.observation_ids[0], "OTU1")
  assert_eq(table.sample_ids[2], "S3")
}

test "biom get_value" {
  let table = @src.BiomTable::new(
    ["A", "B"],
    ["X", "Y"],
    [[10.0, 20.0], [30.0, 40.0]],
    table_id="t",
  )
  assert_eq(table.get_value(0, 0).unwrap(), 10.0)
  assert_eq(table.get_value(0, 1).unwrap(), 20.0)
  assert_eq(table.get_value(1, 1).unwrap(), 40.0)
  // Out of bounds
  match table.get_value(5, 5) {
    Some(_) => fail("should be None")
    None => ()
  }
}

test "biom observation_index" {
  let table = @src.BiomTable::new(
    ["OTU1", "OTU2", "OTU3"],
    ["S1"],
    [[1.0], [2.0], [3.0]],
    table_id="t",
  )
  assert_eq(table.observation_index("OTU2"), 1)
  assert_eq(table.observation_index("OTU3"), 2)
  assert_eq(table.observation_index("MISSING"), -1)
}

test "biom sample_index" {
  let table = @src.BiomTable::new(
    ["O1"],
    ["SA", "SB", "SC"],
    [[1.0, 2.0, 3.0]],
    table_id="t",
  )
  assert_eq(table.sample_index("SB"), 1)
  assert_eq(table.sample_index("SC"), 2)
  assert_eq(table.sample_index("MISSING"), -1)
}

test "biom sample_sums" {
  let table = @src.BiomTable::new(
    ["O1", "O2"],
    ["S1", "S2"],
    [[1.0, 3.0], [2.0, 4.0]],
    table_id="t",
  )
  // S1: 1+2=3, S2: 3+4=7
  let sums = table.sample_sums()
  assert_eq(sums.length(), 2)
  assert_eq(sums[0], 3.0)
  assert_eq(sums[1], 7.0)
}

test "biom observation_sums" {
  let table = @src.BiomTable::new(
    ["O1", "O2"],
    ["S1", "S2"],
    [[1.0, 3.0], [2.0, 4.0]],
    table_id="t",
  )
  // O1: 1+3=4, O2: 2+4=6
  let sums = table.observation_sums()
  assert_eq(sums.length(), 2)
  assert_eq(sums[0], 4.0)
  assert_eq(sums[1], 6.0)
}

test "biom transform" {
  let table = @src.BiomTable::new(
    ["O1"],
    ["S1", "S2"],
    [[2.0, 4.0]],
    table_id="t",
  )
  let doubled = table.transform(fn(v, _) { v * 2.0 })
  assert_eq(doubled.data[0][0], 4.0)
  assert_eq(doubled.data[0][1], 8.0)
  // Original unchanged
  assert_eq(table.data[0][0], 2.0)
}

test "biom normalize_samples" {
  let table = @src.BiomTable::new(
    ["O1", "O2"],
    ["S1", "S2"],
    [[1.0, 3.0], [3.0, 9.0]],
    table_id="t",
  )
  // S1 total: 1+3=4, S2 total: 3+9=12
  let norm = table.normalize_samples()
  // S1: 1/4=0.25, 3/4=0.75
  assert_eq(norm.data[0][0], 0.25)
  assert_eq(norm.data[1][0], 0.75)
  // S2: 3/12=0.25, 9/12=0.75
  assert_eq(norm.data[0][1], 0.25)
  assert_eq(norm.data[1][1], 0.75)
}

test "biom filter_observations" {
  let table = @src.BiomTable::new(
    ["A", "B", "C"],
    ["S1", "S2"],
    [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
    table_id="t",
  )
  let filtered = table.filter_observations(["A", "C"])
  assert_eq(filtered.n_observations(), 2)
  assert_eq(filtered.observation_ids[0], "A")
  assert_eq(filtered.observation_ids[1], "C")
  assert_eq(filtered.data[0][0], 1.0)
  assert_eq(filtered.data[1][0], 5.0)
}

test "biom filter_samples" {
  let table = @src.BiomTable::new(
    ["O1", "O2"],
    ["SA", "SB", "SC"],
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    table_id="t",
  )
  let filtered = table.filter_samples(["SB", "SC"])
  assert_eq(filtered.n_samples(), 2)
  assert_eq(filtered.sample_ids[0], "SB")
  assert_eq(filtered.sample_ids[1], "SC")
  assert_eq(filtered.data[0][0], 2.0)
  assert_eq(filtered.data[1][1], 6.0)
}

test "biom to_tsv" {
  let table = @src.BiomTable::new(
    ["OTU1", "OTU2"],
    ["S1", "S2"],
    [[5.0, 0.0], [3.0, 7.0]],
    table_id="test_tsv",
  )
  let tsv = table.to_tsv()
  // Should contain header
  assert_true(tsv.contains("#OTU ID"))
  assert_true(tsv.contains("S1"))
  assert_true(tsv.contains("S2"))
  assert_true(tsv.contains("OTU1"))
  assert_true(tsv.contains("OTU2"))
  assert_true(tsv.contains("Observation count: 2"))
  assert_true(tsv.contains("Sample count: 2"))
}

test "biom from_tsv" {
  let tsv =
    #|#OTU ID	S1	S2	S3
    #|OTU1	5	0	3
    #|OTU2	0	7	1
    #|# Observation count: 2; Sample count: 3
  match @src.parse_biom_tsv(tsv) {
    Some(table) => {
      assert_eq(table.n_observations(), 2)
      assert_eq(table.n_samples(), 3)
      assert_eq(table.observation_ids[0], "OTU1")
      assert_eq(table.sample_ids[0], "S1")
      assert_eq(table.data[0][0], 5.0)
      assert_eq(table.data[0][2], 3.0)
      assert_eq(table.data[1][1], 7.0)
    }
    None => fail("should parse TSV")
  }
}

test "biom tsv roundtrip" {
  let original = @src.BiomTable::new(
    ["O1", "O2", "O3"],
    ["A", "B"],
    [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
    table_id="roundtrip",
  )
  let tsv = original.to_tsv()
  match @src.parse_biom_tsv(tsv) {
    Some(parsed) => {
      assert_eq(parsed.n_observations(), 3)
      assert_eq(parsed.n_samples(), 2)
      assert_eq(parsed.observation_ids[0], "O1")
      assert_eq(parsed.sample_ids[1], "B")
      assert_eq(parsed.data[0][0], 1.0)
      assert_eq(parsed.data[2][1], 6.0)
    }
    None => fail("TSV roundtrip should work")
  }
}

test "biom merge" {
  let table_a = @src.BiomTable::new(
    ["O1", "O2"],
    ["S1", "S2"],
    [[1.0, 2.0], [3.0, 4.0]],
    table_id="a",
  )
  let table_b = @src.BiomTable::new(
    ["O2", "O3"],
    ["S2", "S3"],
    [[5.0, 6.0], [7.0, 8.0]],
    table_id="b",
  )
  let merged = table_a.merge(table_b)
  // Union: 3 observations (O1, O2, O3), 3 samples (S1, S2, S3)
  assert_eq(merged.n_observations(), 3)
  assert_eq(merged.n_samples(), 3)
  // O1 only in a: [1, 2, 0]
  assert_eq(merged.data[0][0], 1.0)
  assert_eq(merged.data[0][1], 2.0)
  assert_eq(merged.data[0][2], 0.0)
  // O2 in both: [3, 4+5=9, 0+6=6]
  assert_eq(merged.data[1][0], 3.0)
  assert_eq(merged.data[1][1], 9.0)
  assert_eq(merged.data[1][2], 6.0)
  // O3 only in b: [0, 0, 8]
  assert_eq(merged.data[2][0], 0.0)
  assert_eq(merged.data[2][2], 8.0)
}
