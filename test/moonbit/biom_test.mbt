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