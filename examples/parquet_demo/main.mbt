///|
/// parquet 列存内存模型演示: 基因表达表 -> 行组统计 -> 谓词下推 ->
/// 投影下推 -> 字典编码 -> JSON 持久化.
fn main raise {
  // 1. 构建基因表达列存表(3 个行组, fpkm 量级按行组分离)
  let t = @src.ParquetTable::new([
    ("gene", @src.PStr),
    ("chrom", @src.PStr),
    ("fpkm", @src.PDouble),
    ("log2fc", @src.PDouble),
    ("padj", @src.PDouble),
    ("expressed", @src.PBool),
  ])
  let rng = @src.Rng::from_int(42)
  let bases = [0.0, 50.0, 200.0]
  let chroms = ["chr1", "chr2", "chr3"]
  for g in 0..<3 {
    for i in 0..<40 {
      let gi = g * 40 + i
      let fpkm = bases[g] + rng.next_int(40).to_double()
      let padj = (rng.next_int(100).to_double() + 1.0) / 400.0
      ignore(
        t.add_row([
          @src.PvStr("gene" + (gi + 1).to_string()),
          @src.PvStr(chroms[rng.next_int(3)]),
          @src.PvDouble(fpkm),
          @src.PvDouble((rng.next_int(200).to_double() - 100.0) / 100.0),
          @src.PvDouble(padj),
          @src.PvBool(fpkm > 1.0),
        ]),
      )
    }
    ignore(t.flush_row_group())
  }
  println("=== parquet 列存表 (120 行 x 6 列, 3 行组) ===")
  println("  列: " + t.column_names().join(", "))
  println(
    "  行组: " +
    t.num_row_groups().to_string() +
    ", 总行数: " +
    t.num_rows().to_string(),
  )
  // 2. 行组统计(min/max 元数据)
  println("=== 行组统计 (fpkm 列 min/max) ===")
  for g in 0..<3 {
    let st = t.chunk_stats(g, "fpkm")
    let lo = match st.min {
      Some(@src.PvDouble(v)) => v.to_string()
      _ => "?"
    }
    let hi = match st.max {
      Some(@src.PvDouble(v)) => v.to_string()
      _ => "?"
    }
    println("  行组 " + g.to_string() + ": [" + lo + ", " + hi + "]")
  }
  // 3. 谓词下推: fpkm >= 100 只应扫描行组 3
  let hot = t.scan(filters=[
    { column: "fpkm", op: @src.OpGe, value: @src.PvDouble(100.0), },
  ])
  println("=== 谓词下推 scan (fpkm >= 100) ===")
  println(
    "  命中 " +
    hot.rows.length().to_string() +
    " 行; 扫描行组 " +
    hot.groups_scanned.to_string() +
    " / 跳过 " +
    hot.groups_skipped.to_string() +
    " (统计区间裁剪, 未读行组数据)",
  )
  // 4. 投影下推 + 组合谓词
  let sig = t.scan(columns=["gene", "fpkm", "padj"], filters=[
    { column: "padj", op: @src.OpLt, value: @src.PvDouble(0.05), },
    { column: "expressed", op: @src.OpEq, value: @src.PvBool(true), },
  ])
  println("=== 投影 + 谓词 scan (padj < 0.05 AND expressed) ===")
  println(
    "  显著差异基因: " +
    sig.rows.length().to_string() +
    " 个 (只物化 3/6 列)",
  )
  if sig.rows.length() > 0 {
    let row = sig.rows[0]
    let gn = match row[0] {
      @src.PvStr(s) => s
      _ => "?"
    }
    let fv = match row[1] {
      @src.PvDouble(v) => v.to_string()
      _ => "?"
    }
    let pv = match row[2] {
      @src.PvDouble(v) => v.to_string()
      _ => "?"
    }
    println("  示例: " + gn + " fpkm=" + fv + " padj=" + pv)
  }
  // 5. 字典编码(chrom 列高度重复)
  let chrom_col : Array[@src.ParquetValue] = t.read_column("chrom")
  let chrom_strs : Array[String] = []
  for v in chrom_col {
    match v {
      @src.PvStr(s) => chrom_strs.push(s)
      _ => ()
    }
  }
  let (dict, idx) = @src.parquet_dict_encode(chrom_strs)
  println("=== 字典编码 (chrom 列, DICT 编码语义) ===")
  println(
    "  " +
    chrom_strs.length().to_string() +
    " 个值 → 字典 " +
    dict.length().to_string() +
    " 项 (压缩比 " +
    (chrom_strs.length().to_double() / dict.length().to_double()).to_string() +
    "x); 首行索引: " +
    idx[0].to_string(),
  )
  // 6. JSON 持久化往返
  let json = t.to_json()
  let back = @src.parquet_from_json(json).unwrap()
  println("=== JSON 持久化往返 ===")
  println(
    "  行数 " +
    back.num_rows().to_string() +
    " == " +
    t.num_rows().to_string() +
    ", 行组 " +
    back.num_row_groups().to_string() +
    " == " +
    t.num_row_groups().to_string() +
    ", JSON 字节数: " +
    json.length().to_string(),
  )
  println("")
  println(
    "注: 本模块为 parquet 列存内存模型(JSON 持久化), 非 .parquet 二进制",
  )
  println(
    "    容器; 与 .parquet 互操作建议经 JSON/CSV 导入导出衔接.",
  )
}
