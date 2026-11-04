///|
/// cooler 接触矩阵内存模型.
///
/// cooler 是 Hi-C 接触矩阵的社区标准(HDF5 容器: chroms/bins/pixels 三表).
/// 本模块实现其**内存数据模型与分析语义**(bin 网格/像素聚合/ICE 平衡/
/// 多分辨率粗化/HiCDC 衔接), 持久化采用 JSON 而非 HDF5 二进制
/// (MoonBit 生态无 HDF5 读写; rhdf5 为内存模拟). 与 cooler 文件的
/// 互操作建议经 CSV/JSON 导入导出衔接.
///
/// mcool = 同一数据的多个分辨率(bin_size 递增)集合, 见 cooler_coarsen.
pub suberror CoolerError {
  CoolerError(String)
}

///|
pub struct Cooler {
  bin_size : Int32
  chrom_names : Array[String]
  chrom_lengths : Map[String, Int32]
  pixels : Map[Int32, Map[Int32, Int32]]
}

///|
/// 一条像素(接触)记录.
pub struct CoolerPixel {
  bin1 : Int32
  bin2 : Int32
  count : Int32
}

///|
/// bin1 全局索引.
pub fn CoolerPixel::bin1(self : CoolerPixel) -> Int32 {
  self.bin1
}

///|
/// bin2 全局索引.
pub fn CoolerPixel::bin2(self : CoolerPixel) -> Int32 {
  self.bin2
}

///|
/// 接触计数.
pub fn CoolerPixel::count(self : CoolerPixel) -> Int32 {
  self.count
}

///|
/// 新建 cooler, bin_size 为分箱大小(bp).
pub fn Cooler::new(bin_size : Int32) -> Cooler {
  if bin_size <= 0 {
    abort("Cooler bin_size must be positive")
  }
  {
    bin_size,
    chrom_names: [],
    chrom_lengths: Map([], capacity=0),
    pixels: Map([], capacity=0),
  }
}

///|
/// bin 大小(bp).
pub fn Cooler::bin_size(self : Cooler) -> Int32 {
  self.bin_size
}

///|
/// 注册染色体(重复注册按首次生效).
pub fn Cooler::add_chrom(
  self : Cooler,
  name : String,
  length : Int32,
) -> Cooler {
  if !self.chrom_lengths.contains(name) && length > 0 {
    self.chrom_names.push(name)
    self.chrom_lengths[name] = length
  }
  self
}

///|
/// 染色体名列表(注册顺序).
pub fn Cooler::chrom_names(self : Cooler) -> Array[String] {
  self.chrom_names
}

///|
/// 染色体长度(bp).
pub fn Cooler::chrom_length(self : Cooler, name : String) -> Int32? {
  self.chrom_lengths.get(name)
}

///|
/// 全局 bin 总数(按各染色体长度向上取整求和).
pub fn Cooler::bin_count(self : Cooler) -> Int32 {
  let mut total = 0
  for name in self.chrom_names {
    let len = self.chrom_lengths[name]
    total = total + (len + self.bin_size - 1) / self.bin_size
  }
  total
}

///|
/// 染色体首个全局 bin 索引.
fn cooler_chrom_offset(c : Cooler, name : String) -> Int32? {
  let mut offset = 0
  for cn in c.chrom_names {
    if cn == name {
      return Some(offset)
    }
    let len = c.chrom_lengths[cn]
    offset = offset + (len + c.bin_size - 1) / c.bin_size
  }
  None
}

///|
/// 基因组坐标 → 全局 bin 索引; 未知染色体或越界返回 None.
pub fn Cooler::bin_index(self : Cooler, chrom : String, pos : Int32) -> Int32? {
  match self.chrom_lengths.get(chrom) {
    None => None
    Some(len) =>
      if pos < 0 || pos >= len {
        None
      } else {
        match cooler_chrom_offset(self, chrom) {
          None => None
          Some(off) => Some(off + pos / self.bin_size)
        }
      }
  }
}

///|
/// 全局 bin 索引 → (染色体, bin 起始 bp).
pub fn Cooler::bin_location(self : Cooler, bin : Int32) -> (String, Int32)? {
  let mut offset = 0
  for cn in self.chrom_names {
    let len = self.chrom_lengths[cn]
    let nbins = (len + self.bin_size - 1) / self.bin_size
    if bin < offset + nbins {
      return Some((cn, (bin - offset) * self.bin_size))
    }
    offset = offset + nbins
  }
  None
}

///|
/// 添加像素(全局 bin 索引, 半开上三角 bin1 <= bin2), 同格计数累加.
/// 索引倒置自动交换; 越界抛 CoolerError.
pub fn Cooler::add_pixel(
  self : Cooler,
  bin1 : Int32,
  bin2 : Int32,
  count : Int32,
) -> Cooler raise CoolerError {
  let nbins = self.bin_count()
  let (a, b) = if bin1 <= bin2 { (bin1, bin2) } else { (bin2, bin1) }
  if a < 0 || b >= nbins {
    raise CoolerError(
      "pixel bin index out of range: " +
      bin1.to_string() +
      "," +
      bin2.to_string() +
      " (total bins " +
      nbins.to_string() +
      ")",
    )
  }
  match self.pixels.get(a) {
    Some(row) =>
      match row.get(b) {
        Some(v) => row[b] = v + count
        None => row[b] = count
      }
    None => {
      let row : Map[Int32, Int32] = Map([], capacity=1)
      row[b] = count
      self.pixels[a] = row
    }
  }
  self
}

///|
/// 以基因组坐标添加接触(未知染色体/坐标越界抛 CoolerError).
pub fn Cooler::add_contact(
  self : Cooler,
  chrom : String,
  pos1 : Int32,
  pos2 : Int32,
  count : Int32,
) -> Cooler raise CoolerError {
  let b1 = match self.bin_index(chrom, pos1) {
    None =>
      raise CoolerError(
        "invalid coordinate " + chrom + ":" + pos1.to_string() + " for cooler",
      )
    Some(b) => b
  }
  let b2 = match self.bin_index(chrom, pos2) {
    None =>
      raise CoolerError(
        "invalid coordinate " + chrom + ":" + pos2.to_string() + " for cooler",
      )
    Some(b) => b
  }
  self.add_pixel(b1, b2, count)
}

///|
/// 全部像素(按 bin1 升序、bin2 升序).
pub fn Cooler::pixels(self : Cooler) -> Array[CoolerPixel] {
  let out : Array[CoolerPixel] = Array(capacity=0)
  let mut a = 0
  let n = self.bin_count()
  while a < n {
    match self.pixels.get(a) {
      Some(row) => {
        let keys : Array[Int32] = []
        for entry in row {
          keys.push(entry.0)
        }
        let mut i = 0
        while i < keys.length() {
          let mut j = i + 1
          while j < keys.length() {
            if keys[j] < keys[i] {
              let t = keys[i]
              keys[i] = keys[j]
              keys[j] = t
            }
            j = j + 1
          }
          i = i + 1
        }
        for k in keys {
          out.push({ bin1: a, bin2: k, count: row[k], })
        }
      }
      None => ()
    }
    a = a + 1
  }
  out
}

///|
/// 非零像素数.
pub fn Cooler::pixel_count(self : Cooler) -> Int32 {
  let mut total = 0
  for entry in self.pixels {
    total = total + entry.1.length()
  }
  total
}

///|
/// 总接触计数(全部像素 count 之和).
pub fn Cooler::total_counts(self : Cooler) -> Int32 {
  let mut total = 0
  for outer in self.pixels {
    for inner in outer.1 {
      total = total + inner.1
    }
  }
  total
}

///|
/// 读取一格像素计数(无记录返回 None).
pub fn Cooler::get(self : Cooler, bin1 : Int32, bin2 : Int32) -> Int32? {
  let (a, b) = if bin1 <= bin2 { (bin1, bin2) } else { (bin2, bin1) }
  match self.pixels.get(a) {
    Some(row) => row.get(b)
    None => None
  }
}

///|
/// ICE 风格平衡(迭代校正): 逐轮按 bin 行和调整权重, 返回每个 bin 的
/// 权重(行和均匀化后). 低覆盖 bin(行和 0)权重为 0.
/// iterations <= 0 时按 10 轮处理.
pub fn Cooler::balance(
  self : Cooler,
  iterations~ : Int32 = 10,
) -> Array[Double] {
  let n = self.bin_count()
  let weights = Array::make(n, 1.0)
  let rounds = if iterations <= 0 { 10 } else { iterations }
  let mut round = 0
  while round < rounds {
    // 行和(以当前权重加权)
    let row_sum = Array::make(n, 0.0)
    for outer in self.pixels {
      let i = outer.0
      if i < n {
        for inner in outer.1 {
          let j = inner.0
          if j < n {
            let v = inner.1.to_double() * weights[i] * weights[j]
            row_sum[i] = row_sum[i] + v
            row_sum[j] = row_sum[j] + v
          }
        }
      }
    }
    // 更新权重
    let mut norm = 0.0
    let mut cnt = 0
    for i in 0..<n {
      if row_sum[i] > 0.0 {
        weights[i] = weights[i] / row_sum[i].sqrt()
        norm = norm + weights[i]
        cnt = cnt + 1
      } else {
        weights[i] = 0.0
      }
    }
    // 归一化均值 1
    if cnt > 0 {
      let mean = norm / cnt.to_double()
      for i in 0..<n {
        weights[i] = weights[i] / mean
      }
    }
    round = round + 1
  }
  weights
}

///|
/// 粗化(factor 倍大 bin): 相邻 bin 合并计数, 染色体保留,
/// 生成 mcool 风格的低分辨率视图. factor <= 1 时原样复制浅层结构.
pub fn Cooler::coarsen(
  self : Cooler,
  factor : Int32,
) -> Cooler raise CoolerError {
  if factor <= 1 {
    let out = Cooler::new(self.bin_size)
    for name in self.chrom_names {
      ignore(out.add_chrom(name, self.chrom_lengths[name]))
    }
    for p in self.pixels() {
      ignore(out.add_pixel(p.bin1, p.bin2, p.count))
    }
    return out
  }
  let new_bin = self.bin_size * factor
  let out = Cooler::new(new_bin)
  for name in self.chrom_names {
    ignore(out.add_chrom(name, self.chrom_lengths[name]))
  }
  for p in self.pixels() {
    // 全局索引直接除以 factor 聚合(染色体 bin 数对齐按全局序保持)
    ignore(out.add_pixel(p.bin1 / factor, p.bin2 / factor, p.count))
  }
  out
}

///|
/// 转为 HiCDC 分析输入(HiCContact 数组): 跨染色体的像素跳过
/// (HiCDC 为顺式模型), gc/mappability 留默认由调用方补充.
pub fn Cooler::to_hicdc_contacts(self : Cooler) -> Array[HiCContact] {
  let out : Array[HiCContact] = Array(capacity=0)
  for p in self.pixels() {
    match (self.bin_location(p.bin1), self.bin_location(p.bin2)) {
      (Some((c1, _)), Some((c2, _))) =>
        if c1 == c2 {
          out.push(HiCContact::new(c1, p.bin1, p.bin2, p.count, self.bin_size))
        }
      _ => ()
    }
  }
  out
}

///|
/// bins 转 BED 文本(指定染色体, 全部 bin 或仅非零 bin).
pub fn Cooler::bins_as_bed(
  self : Cooler,
  chrom : String,
  nonzero_only~ : Bool = false,
) -> String {
  match self.chrom_lengths.get(chrom) {
    None => ""
    Some(len) => {
      let sb = StringBuilder()
      let nbins = (len + self.bin_size - 1) / self.bin_size
      let offset = match cooler_chrom_offset(self, chrom) {
        Some(o) => o
        None => return ""
      }
      for i in 0..<nbins {
        let start = i * self.bin_size
        let end = if start + self.bin_size < len {
          start + self.bin_size
        } else {
          len
        }
        if nonzero_only {
          let mut nonzero = false
          for j in 0..<nbins {
            if self.get(offset + i, offset + j) is Some(_) ||
              self.get(offset + j, offset + i) is Some(_) {
              nonzero = true
              break
            }
          }
          if !nonzero {
            continue
          }
        }
        sb.write_string(
          chrom + "\t" + start.to_string() + "\t" + end.to_string() + "\n",
        )
      }
      sb.to_string()
    }
  }
}

///|
/// 序列化为 JSON(持久化交换格式, 语义等价 cooler 的三表).
pub fn Cooler::to_json(self : Cooler) -> String {
  let chroms : Array[JsonValue] = Array(capacity=self.chrom_names.length())
  for name in self.chrom_names {
    let m : Map[String, JsonValue] = Map([], capacity=2)
    m["name"] = JsonStr(name)
    m["length"] = JsonNum(self.chrom_lengths[name].to_double())
    chroms.push(JsonObj(m))
  }
  let pxs : Array[JsonValue] = Array(capacity=0)
  for p in self.pixels() {
    let m : Map[String, JsonValue] = Map([], capacity=3)
    m["b1"] = JsonNum(p.bin1.to_double())
    m["b2"] = JsonNum(p.bin2.to_double())
    m["count"] = JsonNum(p.count.to_double())
    pxs.push(JsonObj(m))
  }
  let root : Map[String, JsonValue] = Map([], capacity=3)
  root["bin_size"] = JsonNum(self.bin_size.to_double())
  root["chroms"] = JsonArr(chroms)
  root["pixels"] = JsonArr(pxs)
  json_to_string(JsonObj(root))
}

///|
/// 从 JSON 反序列化(cooler_to_json 的逆).
pub fn cooler_from_json(text : String) -> Cooler? raise CoolerError {
  match parse_json(text) {
    None => None
    Some(JsonObj(root)) => {
      let bin_size = match root.get("bin_size") {
        Some(JsonNum(v)) => v.to_int()
        _ => 0
      }
      if bin_size <= 0 {
        return None
      }
      let cooler = Cooler::new(bin_size)
      match root.get("chroms") {
        Some(JsonArr(chroms)) =>
          for c in chroms {
            match c {
              JsonObj(m) =>
                match (m.get("name"), m.get("length")) {
                  (Some(JsonStr(n)), Some(JsonNum(l))) =>
                    ignore(cooler.add_chrom(n, l.to_int()))
                  _ => ()
                }
              _ => ()
            }
          }
        _ => ()
      }
      match root.get("pixels") {
        Some(JsonArr(pxs)) =>
          for p in pxs {
            match p {
              JsonObj(m) =>
                match (m.get("b1"), m.get("b2"), m.get("count")) {
                  (Some(JsonNum(a)), Some(JsonNum(b)), Some(JsonNum(c))) =>
                    ignore(cooler.add_pixel(a.to_int(), b.to_int(), c.to_int()))
                  _ => ()
                }
              _ => ()
            }
          }
        _ => ()
      }
      Some(cooler)
    }
    _ => None
  }
}
