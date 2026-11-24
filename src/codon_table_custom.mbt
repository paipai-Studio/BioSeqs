///|
/// codon_table_custom.mbt — 统一遗传密码表对象
/// 对应 Biotite `biotite.sequence.codon.CodonTable` / Biopython `Bio.Data.CodonTable`
///
/// 特性:
///   - 对象化封装 27 张 NCBI 标准表 (1-6, 9-16, 21-33)
///   - 密码子 → 氨基酸双向索引 (forward_table / backward_table)
///   - 起始密码子 (start_codons) 与终止密码子 (stop_codons)
///   - IUPAC 简并密码子解析 (展开一致才判定, 否则 'X')
///   - 统一翻译入口 translate_with_table / translate_with_table_id
///
/// 命名说明: 本模块是全库唯一的"任意遗传密码表翻译"权威入口,
/// 原先散落在 skbio3_*/skbio_* 模块中的翻译实现现已收敛到本模块。

///|
/// 遗传密码表相关错误, 对应 Biopython 的 TranslationError。
pub suberror CodonTableError {
  CodonTableError(String)
}

///|
/// 不可变的遗传密码表对象。
///
/// 字段:
///   - table_id: NCBI 密码子表编号 (自定义表为 0)
///   - table_name: 表的官方名称
///   - forward: 密码子 → 氨基酸符号 (含 "*" 表示终止)
///   - backward: 氨基酸符号 → 密码子列表 (不含终止密码子)
///   - starts: 起始密码子列表
///   - stops: 终止密码子列表
pub struct CodonTable {
  table_id : Int32
  table_name : String
  forward : Map[String, String]
  backward : Map[String, Array[String]]
  starts : Array[String]
  stops : Array[String]
}

// ============================================================
// 构造
// ============================================================

///|
/// 复制密码子映射。
///
/// codon_table 返回的是模块级共享 Map, 直接持有会让调用方有机会改到全局数据,
/// 因此所有 CodonTable 都持有一份独立副本。
fn _ct_copy_map(src : Map[String, String]) -> Map[String, String] {
  let dst : Map[String, String] = Map([], capacity=64)
  for codon in src.keys() {
    match src.get(codon) {
      Some(aa) => dst[codon] = aa
      None => ()
    }
  }
  dst
}

///|
/// 构造反向索引: 氨基酸符号 → 密码子列表 (终止密码子 "*" 不计入)。
fn _ct_build_backward(
  forward : Map[String, String],
) -> Map[String, Array[String]] {
  let backward : Map[String, Array[String]] = Map([], capacity=20)
  let keys = forward.keys().to_array()
  let mut i = 0
  while i < keys.length() {
    let codon = keys[i]
    match forward.get(codon) {
      Some(aa) =>
        if aa != "*" {
          match backward.get(aa) {
            Some(list) => list.push(codon)
            None => backward[aa] = [codon]
          }
        }
      None => ()
    }
    i = i + 1
  }
  backward
}

///|
/// 由 NCBI 表编号构造密码子表对象 (不校验编号有效性)。
fn _ct_from_id(id : Int32) -> CodonTable {
  let forward = _ct_copy_map(codon_table(id))
  CodonTable::{
    table_id: id,
    table_name: codon_table_name(id),
    forward,
    backward: _ct_build_backward(forward),
    starts: codon_table_start_codons(id),
    stops: codon_table_stop_codons(id),
  }
}

///|
/// 按 NCBI 表编号构造密码子表对象。
///
/// 参数:
///   - id: NCBI 密码子表编号 (1-6, 9-16, 21-33)
///
/// 异常: 编号不在有效集合内时抛出 CodonTableError。
pub fn codon_table_by_id(id : Int32) -> CodonTable raise CodonTableError {
  if !codon_table_is_valid(id) {
    raise CodonTableError("Invalid NCBI codon table id: " + id.to_string())
  }
  _ct_from_id(id)
}

///|
/// 构造标准遗传密码表 (NCBI Table 1) 对象。
pub fn codon_table_default() -> CodonTable {
  _ct_from_id(1)
}

///|
/// 按 NCBI 表编号构造密码子表对象, 编号无效时回退到标准表 (Table 1)。
///
/// 与 codon_table_by_id 的区别: 本函数不会抛出异常, 行为与
/// codon_table(id) 一致 (非法编号静默回退到表 1)。
pub fn codon_table_or_default(id : Int32) -> CodonTable {
  if codon_table_is_valid(id) {
    _ct_from_id(id)
  } else {
    _ct_from_id(1)
  }
}

///|
/// 构造细菌/古菌/植物质体密码子表 (NCBI Table 11) 对象。
pub fn codon_table_bacterial() -> CodonTable {
  _ct_from_id(11)
}

///|
/// 构造交替酵母核密码子表 (NCBI Table 12) 对象。
pub fn codon_table_alt_yeast() -> CodonTable {
  _ct_from_id(12)
}

///|
/// 构造纤毛虫核密码子表 (NCBI Table 6) 对象。
pub fn codon_table_ciliate() -> CodonTable {
  _ct_from_id(6)
}

///|
/// 构造脊椎动物线粒体密码子表 (NCBI Table 2) 对象。
pub fn codon_table_vertebrate_mito() -> CodonTable {
  _ct_from_id(2)
}

///|
/// 由自定义 codon → 氨基酸 映射构造密码子表对象。
///
/// 参数:
///   - forward: 密码子 → 氨基酸符号 的映射
///   - name: 表名称 (默认 "Custom")
///   - start_codons: 起始密码子列表 (默认空)
///   - stop_codons: 终止密码子列表 (默认空)
pub fn codon_table_from_map(
  forward : Map[String, String],
  name~ : String = "Custom",
  start_codons~ : Array[String] = [],
  stop_codons~ : Array[String] = [],
) -> CodonTable {
  let own = _ct_copy_map(forward)
  CodonTable::{
    table_id: 0,
    table_name: name,
    forward: own,
    backward: _ct_build_backward(own),
    starts: start_codons.copy(),
    stops: stop_codons.copy(),
  }
}

// ============================================================
// 密码子顺序
// ============================================================

///|
/// 64 个密码子的标准顺序 (TCAG 三重循环), 与 codon_fixed_table 的索引一致。
///
/// 用于需要按索引排列密码子 (如 RSCU / CAI 权重数组) 的场景。
pub fn canonical_codon_order() -> Array[String] {
  let bases = "TCAG"
  let order : Array[String] = Array(capacity=64)
  for i in 0..<4 {
    for j in 0..<4 {
      for k in 0..<4 {
        let sb = StringBuilder()
        sb.write_char(bases.unsafe_get(i).unsafe_to_char())
        sb.write_char(bases.unsafe_get(j).unsafe_to_char())
        sb.write_char(bases.unsafe_get(k).unsafe_to_char())
        order.push(sb.to_string())
      }
    }
  }
  order
}

// ============================================================
// 访问器
// ============================================================

///|
/// NCBI 表编号 (自定义表为 0)。
pub fn CodonTable::id(self : CodonTable) -> Int32 {
  self.table_id
}

///|
/// 表名称。
pub fn CodonTable::name(self : CodonTable) -> String {
  self.table_name
}

///|
/// 正向索引: 密码子 → 氨基酸符号 (返回副本, 修改不影响表对象)。
pub fn CodonTable::forward_table(self : CodonTable) -> Map[String, String] {
  _ct_copy_map(self.forward)
}

///|
/// 反向索引: 氨基酸符号 → 密码子列表 (不含终止密码子, 返回副本)。
pub fn CodonTable::backward_table(
  self : CodonTable,
) -> Map[String, Array[String]] {
  let out : Map[String, Array[String]] = Map([], capacity=20)
  for aa in self.backward.keys() {
    match self.backward.get(aa) {
      Some(list) => out[aa] = list.copy()
      None => ()
    }
  }
  out
}

///|
/// 起始密码子列表。
pub fn CodonTable::start_codons(self : CodonTable) -> Array[String] {
  self.starts.copy()
}

///|
/// 终止密码子列表。
pub fn CodonTable::stop_codons(self : CodonTable) -> Array[String] {
  self.stops.copy()
}

///|
/// 某氨基酸对应的密码子列表。
pub fn CodonTable::aa_to_codons(
  self : CodonTable,
  aa : String,
) -> Array[String] {
  match self.backward.get(aa) {
    Some(list) => list.copy()
    None => []
  }
}

// ============================================================
// 密码子解析
// ============================================================

///|
/// IUPAC 掩码展开为碱基字符串列表。
fn _ct_bases(mask : Int32) -> Array[String] {
  let bases : Array[String] = []
  if (mask & 1) != 0 {
    bases.push("A")
  }
  if (mask & 2) != 0 {
    bases.push("C")
  }
  if (mask & 4) != 0 {
    bases.push("G")
  }
  if (mask & 8) != 0 {
    bases.push("T")
  }
  bases
}

///|
/// 展开简并密码子为全部非简并密码子。非法字符或长度不为 3 时返回空数组。
fn _ct_expand_codon(codon : String) -> Array[String] {
  let up = codon.to_upper()
  if up.length() != 3 {
    return []
  }
  let m0 = iupac_base_mask[up.unsafe_get(0).to_int()]
  let m1 = iupac_base_mask[up.unsafe_get(1).to_int()]
  let m2 = iupac_base_mask[up.unsafe_get(2).to_int()]
  if m0 == 0 || m1 == 0 || m2 == 0 {
    return []
  }
  let b0 = _ct_bases(m0)
  let b1 = _ct_bases(m1)
  let b2 = _ct_bases(m2)
  let result : Array[String] = []
  let mut i = 0
  while i < b0.length() {
    let mut j = 0
    while j < b1.length() {
      let mut k = 0
      while k < b2.length() {
        result.push(b0[i] + b1[j] + b2[k])
        k = k + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// 密码子 → 氨基酸字符 (支持 IUPAC 简并)。
///
/// 规则:
///   - 完全确定的密码子直接查表
///   - 简并密码子展开后, 若所有展开结果一致则取该氨基酸
///   - 全部展开均为终止密码子时返回 '*'
///   - 无法判定 (展开结果不一致或含非法碱基) 返回 'X'
pub fn CodonTable::codon_to_aa(self : CodonTable, codon : String) -> Char {
  let up = codon.to_upper()
  match self.forward.get(up) {
    Some(aa) =>
      if aa.length() == 1 {
        aa.unsafe_get(0).unsafe_to_char()
      } else {
        'X'
      }
    None => {
      let expansions = _ct_expand_codon(up)
      if expansions.length() == 0 {
        return 'X'
      }
      let mut result = 'X'
      let mut first = true
      let mut i = 0
      while i < expansions.length() {
        match self.forward.get(expansions[i]) {
          Some(aa) => {
            let ch = if aa.length() == 1 {
              aa.unsafe_get(0).unsafe_to_char()
            } else {
              'X'
            }
            if first {
              result = ch
              first = false
            } else if ch != result {
              return 'X'
            }
          }
          None => return 'X'
        }
        i = i + 1
      }
      result
    }
  }
}

///|
/// 判断密码子是否为该表的起始密码子。
pub fn CodonTable::is_start_codon(self : CodonTable, codon : String) -> Bool {
  let up = codon.to_upper()
  let mut i = 0
  while i < self.starts.length() {
    if self.starts[i] == up {
      return true
    }
    i = i + 1
  }
  false
}

///|
/// 判断密码子是否为该表的终止密码子。
pub fn CodonTable::is_stop_codon(self : CodonTable, codon : String) -> Bool {
  let up = codon.to_upper()
  let mut i = 0
  while i < self.stops.length() {
    if self.stops[i] == up {
      return true
    }
    i = i + 1
  }
  false
}

// ============================================================
// 翻译
// ============================================================

///|
/// 按指定阅读框翻译 DNA 序列。
///
/// 参数:
///   - seq: DNA 序列
///   - frame: 阅读框偏移 (0/1/2), 序列开头被跳过的碱基数
///
/// 终止密码子翻译为 '*', 无法判定翻译为 'X'。
pub fn CodonTable::translate_frame(
  self : CodonTable,
  seq : String,
  frame : Int32,
) -> String {
  let up = seq.to_upper()
  let sb = StringBuilder()
  let mut i = frame
  while i + 3 <= up.length() {
    let codon = up[i:i + 3].to_owned()
    let ch = self.codon_to_aa(codon)
    sb.write_char(ch)
    i = i + 3
  }
  sb.to_string()
}

///|
/// 翻译 DNA 序列 (第 0 阅读框)。
pub fn CodonTable::translate(self : CodonTable, seq : String) -> String {
  self.translate_frame(seq, 0)
}

///|
/// 翻译 CDS 序列: 起始密码子统一翻译为 M, 并在首个终止密码子处结束。
///
/// 若序列首密码子不是该表的起始密码子, 则按普通翻译处理。
pub fn CodonTable::translate_cds(self : CodonTable, seq : String) -> String {
  let up = seq.to_upper()
  let sb = StringBuilder()
  let mut i = 0
  let mut is_first = true
  while i + 3 <= up.length() {
    let codon = up[i:i + 3].to_owned()
    let ch = if is_first && self.is_start_codon(codon) {
      'M'
    } else {
      self.codon_to_aa(codon)
    }
    sb.write_char(ch)
    is_first = false
    if ch == '*' {
      break
    }
    i = i + 3
  }
  sb.to_string()
}

///|
/// 统一翻译入口: 使用指定密码子表翻译 DNA 序列。
///
/// 参数:
///   - seq: DNA 序列
///   - table: 密码子表对象
pub fn translate_with_table(seq : String, table : CodonTable) -> String {
  table.translate(seq)
}

///|
/// 统一翻译入口: 按 NCBI 表编号翻译 DNA 序列。
///
/// 参数:
///   - seq: DNA 序列
///   - table_id: NCBI 密码子表编号
///
/// 异常: 编号无效时抛出 CodonTableError。
pub fn translate_with_table_id(
  seq : String,
  table_id : Int32,
) -> String raise CodonTableError {
  let table = codon_table_by_id(table_id)
  table.translate(seq)
}

///|
/// 统一翻译入口: 使用原始 codon → 氨基酸 映射翻译 DNA 序列。
///
/// 自动构造临时密码子表对象, 支持 IUPAC 简并密码子与大小写归一。
///
/// 参数:
///   - seq: DNA 序列
///   - codon_map: 密码子 → 氨基酸 映射
pub fn translate_with_codon_map(
  seq : String,
  codon_map : Map[String, String],
) -> String {
  codon_table_from_map(codon_map).translate(seq)
}

///|
/// 统一翻译入口: 按 NCBI 表编号翻译 CDS 序列。
///
/// 起始密码子统一翻译为 M, 遇到首个终止密码子即结束。
///
/// 参数:
///   - seq: CDS 序列
///   - table_id: NCBI 密码子表编号
///
/// 异常: 编号无效时抛出 CodonTableError。
pub fn translate_cds_with_table_id(
  seq : String,
  table_id : Int32,
) -> String raise CodonTableError {
  let table = codon_table_by_id(table_id)
  table.translate_cds(seq)
}
