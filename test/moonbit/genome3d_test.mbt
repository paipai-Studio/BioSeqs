// ===========================================================================
// genome3d_test.mbt — 3D 基因组建模测试(接触图 -> cMDS -> PDB)
// ===========================================================================

///|
/// 构造模拟接触图: 线性基因组邻近衰减 + 一个域内强接触.
fn genome3d_test_contacts(n : Int32) -> Array[Array[Double]] {
  let c : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    c.push(row)
  }
  // 接触强度随 |i-j| 衰减
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let sep = if j - i > 0 { j - i } else { i - j }
      let v = if sep <= 5 {
        10.0 / sep.to_double()
      } else if sep <= 15 {
        2.0 / sep.to_double()
      } else {
        0.5 / sep.to_double()
      }
      c[i][j] = v
      c[j][i] = v
    }
  }
  // 域内强接触(第二个域, 边界保护)
  let dom_end = if n / 2 + 5 < n { n / 2 + 5 } else { n }
  for i = n / 2; i < dom_end; i = i + 1 {
    for j = n / 2; j < dom_end; j = j + 1 {
      if i != j {
        c[i][j] = c[i][j] + 5.0
      }
    }
  }
  c
}

///|
test "genome3d contact to distance" {
  let contacts = genome3d_test_contacts(10)
  let d = @src.genome3d_contact_to_distance(contacts)
  // 对角 0
  assert_eq(d[0][0], 0.0)
  // 非零接触倒数: c[0][1]=10 -> d=0.1
  assert_true(d[0][1] > 0.099 && d[0][1] < 0.101)
  // 对称
  assert_eq(d[0][1], d[1][0])
  // 单调性: 更远接触 -> 更大距离
  assert_true(d[0][1] < d[0][9])
}

///|
test "genome3d contact to distance zero contact uses dmax" {
  let contacts : Array[Array[Double]] = [
    [0.0, 2.0, 0.0],
    [2.0, 0.0, 4.0],
    [0.0, 4.0, 0.0],
  ]
  let d = @src.genome3d_contact_to_distance(contacts, d_max=5.0)
  assert_eq(d[0][2], 5.0)
  assert_eq(d[0][1], 0.5)
  assert_eq(d[1][2], 0.25)
}

///|
test "genome3d contact to distance validation" {
  let mut raised = false
  try {
    let d = @src.genome3d_contact_to_distance([])
    ignore(d)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
  raised = false
  try {
    let d = @src.genome3d_contact_to_distance([[1.0, 2.0], [2.0]])
    ignore(d)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "genome3d cmds recovers linear structure" {
  // 完美线性距离(1D 基因组): d_ij = |i-j|
  let n = 12
  let d : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      let s = i - j
      row.push(if s >= 0 { s.to_double() } else { (-s).to_double() })
    }
    d.push(row)
  }
  let (coords, evals) = @src.genome3d_cmds(d)
  assert_eq(coords.length(), n)
  // 1D 结构: 第一维特征值占绝对主导
  assert_true(evals[0] > 0.0)
  if evals.length() > 1 {
    assert_true(evals[0] > evals[1] * 5.0)
  }
  // 第一维坐标应单调(线性排序)
  let mut monotonic = true
  for i = 1; i < n; i = i + 1 {
    if coords[i][0] < coords[i - 1][0] - 1.0e-9 {
      // 允许符号翻转: 单调递减也算
      monotonic = monotonic && false
      break
    }
  }
  // 单调递增或递减均可(特征向量符号任意)
  let increasing = {
    let mut ok = true
    for i = 1; i < n; i = i + 1 {
      if coords[i][0] < coords[i - 1][0] - 1.0e-9 {
        ok = false
        break
      }
    }
    ok
  }
  let decreasing = {
    let mut ok = true
    for i = 1; i < n; i = i + 1 {
      if coords[i][0] > coords[i - 1][0] + 1.0e-9 {
        ok = false
        break
      }
    }
    ok
  }
  assert_true(increasing || decreasing)
  ignore(monotonic)
}

///|
test "genome3d model full flow" {
  let contacts = genome3d_test_contacts(20)
  let model = @src.genome3d_model(contacts)
  assert_eq(model.n_beads, 20)
  // 邻近 bead 的模型距离应小于远端 bead
  let near = model.model_distance(5, 6)
  let far = model.model_distance(5, 18)
  assert_true(near < far)
  // 距离相关(cMDS 对大动态范围距离存在 3D 嵌入压缩, 阈值 0.5)
  assert_true(model.stress_correlation > 0.5)
}

///|
test "genome3d model summary" {
  let contacts = genome3d_test_contacts(10)
  let model = @src.genome3d_model(contacts)
  let s = model.summary()
  assert_true(s.contains("genome3d"))
  assert_true(s.contains("10 beads"))
  assert_true(s.contains("distance correlation"))
}

///|
test "genome3d model distance out of range raises" {
  let contacts = genome3d_test_contacts(5)
  let model = @src.genome3d_model(contacts)
  let mut raised = false
  try {
    let d = model.model_distance(0, 99)
    ignore(d)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "genome3d pdb output format" {
  let contacts = genome3d_test_contacts(8)
  let model = @src.genome3d_model(contacts)
  let pdb = model.to_pdb()
  let lines = pdb.split("\n").to_array()
  // ATOM 行数 = bead 数
  let mut n_atoms = 0
  for l in lines {
    if l.to_owned().has_prefix("ATOM") {
      n_atoms = n_atoms + 1
    }
  }
  assert_eq(n_atoms, 8)
  assert_true(pdb.contains("TER"))
  assert_true(pdb.contains("END"))
  assert_true(pdb.contains("CA"))
  // 固定列: ATOM 行长 80(含尾 C 元素)
  let atom_line = {
    let mut found = ""
    for l in lines {
      let s = l.to_owned()
      if s.has_prefix("ATOM") {
        found = s
        break
      }
    }
    found
  }
  assert_true(atom_line.length() >= 66)
  // 坐标列可解析
  assert_true(atom_line[30:38].to_owned().trim().length() > 0)
}

///|
test "genome3d pdb roundtrip" {
  let contacts = genome3d_test_contacts(8)
  let model = @src.genome3d_model(contacts)
  let pdb = model.to_pdb()
  // 解析回坐标
  let parsed = @src.genome3d_parse_pdb_ca(pdb)
  assert_eq(parsed.length(), 8)
  // 坐标往返一致(3 位小数精度内)
  for i = 0; i < 8; i = i + 1 {
    let dx = parsed[i][0] - model.coordinates[i][0]
    let dy = parsed[i][1] - model.coordinates[i][1]
    let dz = parsed[i][2] - model.coordinates[i][2]
    assert_true(dx > -0.002 && dx < 0.002)
    assert_true(dy > -0.002 && dy < 0.002)
    assert_true(dz > -0.002 && dz < 0.002)
  }
}

///|
test "genome3d parse pdb empty raises" {
  let mut raised = false
  try {
    let p = @src.genome3d_parse_pdb_ca("no atoms here\n")
    ignore(p)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "genome3d domain compaction" {
  // 带强域的接触图: 域内 bead 应比域外更近
  let contacts = genome3d_test_contacts(24)
  let model = @src.genome3d_model(contacts)
  // 域 [12,17) 内部距离 vs 与远处 bead 的距离
  let intra = model.model_distance(13, 15)
  let inter = model.model_distance(13, 22)
  assert_true(intra < inter)
}
