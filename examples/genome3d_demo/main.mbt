///|
/// 3D 基因组建模演示: 接触图 -> 距离转换 -> cMDS -> 3D 坐标 -> PDB.
fn main raise {
  println("=== 3D 基因组建模 (接触图 -> cMDS -> PDB, ShRec 风格) ===")

  // 模拟接触图: 线性邻近衰减 + 第二域内强接触
  let n = 20
  let contacts : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    contacts.push(Array::make(n, 0.0))
  }
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
      contacts[i][j] = v
      contacts[j][i] = v
    }
  }
  let dom_end = if n / 2 + 5 < n { n / 2 + 5 } else { n }
  for i = n / 2; i < dom_end; i = i + 1 {
    for j = n / 2; j < dom_end; j = j + 1 {
      if i != j {
        contacts[i][j] = contacts[i][j] + 5.0
      }
    }
  }
  println("接触图 \{n}x\{n}(含 [10,15) 强接触域)")

  // 1. 距离转换
  println("\n1. 距离转换(ShRec3D 倒数变换)")
  let d = @src.genome3d_contact_to_distance(contacts)
  println("  d[0][1]=\{d[0][1]} (接触 10 -> 距离 0.1)")
  println("  d[0][19]=\{d[0][19]} (弱接触 -> 大距离)")

  // 2. cMDS 建模
  println("\n2. 经典 MDS 建模")
  let model = @src.genome3d_model(contacts)
  println("  " + model.summary())

  // 3. 结构验证
  println("\n3. 结构验证")
  let near = model.model_distance(5, 6)
  let far = model.model_distance(5, 18)
  println("  邻近 bead(5,6) 距离: \{near}")
  println("  远端 bead(5,18) 距离: \{far}")
  println("  邻近 < 远端: \{near < far}")
  // 域内压缩
  let intra = model.model_distance(13, 15)
  let inter = model.model_distance(13, 19)
  println(
    "  域内 bead(13,15): \{intra} < 域外(13,19): \{inter}: \{intra < inter}",
  )

  // 4. 坐标输出
  println("\n4. 3D 坐标(前 5 个 bead)")
  for i = 0; i < 5; i = i + 1 {
    let p = model.coordinates[i]
    println("  bead\{i}: (\{p[0]}, \{p[1]}, \{p[2]})")
  }

  // 5. PDB 输出与解析闭环
  println("\n5. PDB 输出(与结构栈闭环)")
  let pdb = model.to_pdb()
  let lines = pdb.split("\n").to_array()
  println("  ATOM 行:")
  for li = 0; li < lines.length() && li < 5; li = li + 1 {
    let s = lines[li].to_owned()
    if s.has_prefix("ATOM") || s.has_prefix("REMARK") {
      println("  | \{s}")
    }
  }
  let parsed = @src.genome3d_parse_pdb_ca(pdb)
  println("  解析回 \{parsed.length()} 个 CA 坐标(roundtrip 闭环)")
  println("  bead0 往返: (\{parsed[0][0]}, \{parsed[0][1]}, \{parsed[0][2]})")

  println(
    "\n完成: 接触图 -> 距离 -> cMDS 3D 坐标 -> PDB 全链路本地完成.",
  )
}
