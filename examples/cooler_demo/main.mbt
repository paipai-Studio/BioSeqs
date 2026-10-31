///|
/// cooler 接触矩阵演示: 构建矩阵 -> ICE 平衡 -> mcool 多分辨率 -> HiCDC 衔接.
fn main raise {
  // 1. 合成 Hi-C 矩阵: 衰减背景 + 一个局部环(强接触)
  let c = @src.Cooler::new(10000)
  ignore(c.add_chrom("chr1", 2000000))
  let rng = @src.Rng::from_int(42)
  let nbins = 200
  // 背景: 随距离幂律衰减的随机接触
  for i in 0..<nbins {
    for j in i..<nbins {
      let d = j - i
      let intensity = 20.0 / (d + 1).to_double().sqrt()
      let k = rng.next_int(3) +
        (if intensity > 1.0 { intensity.to_int() } else { 1 })
      if k > 0 {
        ignore(c.add_pixel(i, j, k))
      }
    }
  }
  // 局部环: bin 50-58 之间强接触
  for i in 50..<59 {
    for j in i..<59 {
      ignore(c.add_pixel(i, j, 30))
    }
  }
  println("=== cooler 构建 (bin_size=10kb, chr1: 2Mb) ===")
  println("  bins: " + c.bin_count().to_string())
  println("  非零像素: " + c.pixel_count().to_string())
  println("  总计数: " + c.total_counts().to_string())
  // 2. ICE 平衡
  let weights = c.balance()
  let mut w_lo = 1.0
  let mut w_hi = 0.0
  for w in weights {
    if w < w_lo {
      w_lo = w
    }
    if w > w_hi {
      w_hi = w
    }
  }
  println("=== ICE 平衡(10 轮) ===")
  println("  bin 权重范围: " + w_lo.to_string() + " ~ " + w_hi.to_string())
  // 3. mcool 多分辨率
  let mid = c.coarsen(2)
  let coarse = c.coarsen(4)
  println("=== mcool 多分辨率 ===")
  println(
    "  10kb -> 20kb -> 40kb, 总计数守恒: " +
    c.total_counts().to_string() +
    " == " +
    mid.total_counts().to_string() +
    " == " +
    coarse.total_counts().to_string(),
  )
  // 4. HiCDC 衔接: 差异分析输入
  let contacts = c.to_hicdc_contacts()
  println("=== HiCDC 衔接 ===")
  println("  顺式接触记录: " + contacts.length().to_string() + " 条")
  // 5. JSON 持久化往返
  let json = c.to_json()
  let back = @src.cooler_from_json(json).unwrap()
  println(
    "=== JSON 往返 === 总计数 " +
    back.total_counts().to_string() +
    " == " +
    c.total_counts().to_string() +
    ", 像素数 " +
    back.pixel_count().to_string() +
    " == " +
    c.pixel_count().to_string(),
  )
  // 6. BED 输出
  let bed = c.bins_as_bed("chr1", nonzero_only=true)
  let nlines = {
    let mut n = 0
    for ch in bed {
      if ch == '\n' {
        n = n + 1
      }
    }
    n
  }
  println(
    "=== bins -> BED(非零) === 非零 bin: " + nlines.to_string() + " 行",
  )
  println("")
  println(
    "注: 本模块为 cooler 内存数据模型(JSON 持久化), 非 HDF5 二进制容器;",
  )
  println("    与 .cool 文件互操作建议经 CSV/JSON 导入导出衔接.")
}
