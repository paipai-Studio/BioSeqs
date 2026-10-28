///|
/// parquet 列式存储内存模型与读取语义.
///
/// parquet 是分析型列存事实标准(行组/列块/统计元数据/谓词下推/投影).
/// 本模块实现其**内存数据模型与高效读取语义**(行组划分、列扫描、
/// min/max 统计驱动的谓词下推跳行组、投影下推、字典编码),
/// 持久化采用 JSON 而非 .parquet 二进制容器(MoonBit 生态无
/// thrift/snappy/zstd 解码栈). 与 .parquet 文件的互操作建议经
/// JSON/CSV 导入导出衔接.
///
/// 读取语义对齐 parquet: NULL 不满足任何比较谓词; 统计区间
/// [min, max] 与谓词不相交的行组整组跳过(谓词下推核心收益).
pub suberror ParquetError {
  ParquetError(String)
}

///|
/// 列物理类型(parquet 常用子集).
pub(all) enum ParquetType {
  PBool
  PInt32
  PInt64
  PDouble
  PStr
} derive(Eq)

///|
pub extend ParquetType with Eq::{equal, not_equal}

///|
/// 列值(NULL 语义对齐 parquet 可空列).
pub(all) enum ParquetValue {
  PvNull
  PvBool(Bool)
  PvInt(Int32)
  PvLong(Int64)
  PvDouble(Double)
  PvStr(String)
} derive(Eq, Debug)

///|
pub extend ParquetValue with Eq::{equal, not_equal}

///|
pub extend ParquetValue with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// 谓词算子.
pub(all) enum ParquetOp {
  OpEq
  OpNe
  OpLt
  OpLe
  OpGt
  OpGe
}

///|
/// 谓词(列名 + 算子 + 值); 类型须与列类型一致才可命中.
pub(all) struct ParquetFilter {
  column : String
  op : ParquetOp
  value : ParquetValue
}

///|
/// 列块(一个行组内某列的全部值, 含 NULL 槽位).
pub struct ColumnChunk {
  column : String
  ptype : ParquetType
  values : Array[ParquetValue]
}

///|
/// 行组(水平分片: 每列一个 chunk, 各 chunk 等长同序).
pub struct RowGroup {
  num_rows : Int32
  columns : Array[ColumnChunk]
}

///|
/// 列存表(schema + 行组集合; 未 flush 的行不可查询, 对齐
/// parquet 行组落盘后才可读的语义).
pub struct ParquetTable {
  schema : Array[(String, ParquetType)]
  row_groups : Array[RowGroup]
  pending : Array[Array[ParquetValue]]
}

///|
/// 列块统计(parquet 元数据语义: 谓词下推依据).
pub struct ParquetStats {
  min : ParquetValue?
  max : ParquetValue?
  null_count : Int32
}

///|
/// scan 结果(含谓词下推命中情况: 跳过的行组数可直接断言下推生效).
pub struct ParquetScanResult {
  rows : Array[Array[ParquetValue]]
  groups_scanned : Int32
  groups_skipped : Int32
}

///|
/// 新建列存表(schema 为 (列名, 类型) 数组; 列名重复抛 ParquetError).
pub fn ParquetTable::new(
  schema : Array[(String, ParquetType)],
) -> ParquetTable raise ParquetError {
  let seen : Map[String, Bool] = Map([], capacity=0)
  for pair in schema {
    if seen.contains(pair.0) {
      raise ParquetError("duplicate column in schema: " + pair.0)
    }
    seen[pair.0] = true
  }
  let s : Array[(String, ParquetType)] = Array(capacity=schema.length())
  for pair in schema {
    s.push(pair)
  }
  { schema: s, row_groups: [], pending: [], }
}

///|
/// schema(列名+类型, 构建顺序).
pub fn ParquetTable::schema(
  self : ParquetTable,
) -> Array[(String, ParquetType)] {
  self.schema
}

///|
/// 列名列表.
pub fn ParquetTable::column_names(self : ParquetTable) -> Array[String] {
  let out : Array[String] = Array(capacity=self.schema.length())
  for pair in self.schema {
    out.push(pair.0)
  }
  out
}

///|
/// 已 flush 行数(不含 pending).
pub fn ParquetTable::num_rows(self : ParquetTable) -> Int32 {
  let mut n = 0
  for rg in self.row_groups {
    n = n + rg.num_rows
  }
  n
}

///|
/// 行组数.
pub fn ParquetTable::num_row_groups(self : ParquetTable) -> Int32 {
  self.row_groups.length()
}

///|
/// 未 flush 的待切行数.
pub fn ParquetTable::pending_rows(self : ParquetTable) -> Int32 {
  self.pending.length()
}

///|
/// 追加一行(宽度与类型须匹配 schema; PvNull 任意列可空).
/// 行先进入 pending, 调用 flush_row_group 后才可被 scan 读取.
pub fn ParquetTable::add_row(
  self : ParquetTable,
  values : Array[ParquetValue],
) -> ParquetTable raise ParquetError {
  if values.length() != self.schema.length() {
    raise ParquetError(
      "row width " +
      values.length().to_string() +
      " != schema " +
      self.schema.length().to_string(),
    )
  }
  for i in 0..<values.length() {
    match parquet_type_of(values[i]) {
      None => ()
      Some(t) =>
        if t != self.schema[i].1 {
          raise ParquetError(
            "value type mismatch at column " +
            self.schema[i].0 +
            ": " +
            parquet_type_tag(t) +
            " != " +
            parquet_type_tag(self.schema[i].1),
          )
        }
    }
  }
  self.pending.push(values)
  self
}

///|
/// 将 pending 行切成一个行组(空 pending 不产生行组), 返回该组行数.
pub fn ParquetTable::flush_row_group(self : ParquetTable) -> Int32 {
  let n = self.pending.length()
  if n == 0 {
    return 0
  }
  let cols : Array[ColumnChunk] = Array(capacity=self.schema.length())
  for i in 0..<self.schema.length() {
    let pair = self.schema[i]
    let vals : Array[ParquetValue] = Array(capacity=n)
    for j in 0..<n {
      vals.push(self.pending[j][i])
    }
    cols.push({ column: pair.0, ptype: pair.1, values: vals, })
  }
  self.row_groups.push({ num_rows: n, columns: cols, })
  let mut k = n
  while k > 0 {
    ignore(self.pending.pop())
    k = k - 1
  }
  n
}

///|
/// 全表列扫描(跨行组按行序拼接; 未知列名抛 ParquetError).
pub fn ParquetTable::read_column(
  self : ParquetTable,
  name : String,
) -> Array[ParquetValue] raise ParquetError {
  let idx = match parquet_schema_index(self).get(name) {
    None => raise ParquetError("unknown column: " + name)
    Some(i) => i
  }
  let out : Array[ParquetValue] = []
  for rg in self.row_groups {
    for v in rg.columns[idx].values {
      out.push(v)
    }
  }
  out
}

///|
/// 指定行组指定列的统计(min/max 忽略 NULL; Bool 列无 min/max).
pub fn ParquetTable::chunk_stats(
  self : ParquetTable,
  rg_idx : Int32,
  name : String,
) -> ParquetStats raise ParquetError {
  if rg_idx < 0 || rg_idx >= self.row_groups.length() {
    raise ParquetError("row group index out of range: " + rg_idx.to_string())
  }
  let idx = match parquet_schema_index(self).get(name) {
    None => raise ParquetError("unknown column: " + name)
    Some(i) => i
  }
  parquet_chunk_stats(self.row_groups[rg_idx].columns[idx])
}

///|
/// 扫描: 投影下推(只读所需列) + 谓词下推(统计区间不相交的行组整组跳过).
/// columns 为空取全列; NULL 值不满足任何比较谓词(对齐 parquet 语义).
pub fn ParquetTable::scan(
  self : ParquetTable,
  columns~ : Array[String] = [],
  filters~ : Array[ParquetFilter] = [],
) -> ParquetScanResult raise ParquetError {
  let sidx = parquet_schema_index(self)
  let proj : Array[Int32] = []
  if columns.is_empty() {
    for i in 0..<self.schema.length() {
      proj.push(i)
    }
  } else {
    for cn in columns {
      match sidx.get(cn) {
        None => raise ParquetError("unknown column: " + cn)
        Some(i) => proj.push(i)
      }
    }
  }
  let fidx : Array[Int32] = Array(capacity=filters.length())
  for f in filters {
    match sidx.get(f.column) {
      None => raise ParquetError("unknown filter column: " + f.column)
      Some(i) => fidx.push(i)
    }
  }
  let rows : Array[Array[ParquetValue]] = []
  let mut scanned = 0
  let mut skipped = 0
  for rg in self.row_groups {
    // 谓词下推: 逐谓词用行组统计区间裁剪
    let mut may = true
    let mut fi = 0
    while fi < filters.length() {
      let st = parquet_chunk_stats(rg.columns[fidx[fi]])
      if !parquet_range_may_match(st.min, st.max, filters[fi]) {
        may = false
        break
      }
      fi = fi + 1
    }
    if may {
      scanned = scanned + 1
      for r in 0..<rg.num_rows {
        let mut ok = true
        let mut fj = 0
        while fj < filters.length() {
          if !parquet_value_matches(rg.columns[fidx[fj]].values[r], filters[fj]) {
            ok = false
            break
          }
          fj = fj + 1
        }
        if ok {
          let row : Array[ParquetValue] = Array(capacity=proj.length())
          for pi in proj {
            row.push(rg.columns[pi].values[r])
          }
          rows.push(row)
        }
      }
    } else {
      skipped = skipped + 1
    }
  }
  { rows, groups_scanned: scanned, groups_skipped: skipped, }
}

///|
/// 全表行式视图(列转行; 含全部已 flush 行组, 不含 pending).
pub fn ParquetTable::to_rows(self : ParquetTable) -> Array[Array[ParquetValue]] {
  let out : Array[Array[ParquetValue]] = []
  for rg in self.row_groups {
    for r in 0..<rg.num_rows {
      let row : Array[ParquetValue] = Array(capacity=rg.columns.length())
      for chunk in rg.columns {
        row.push(chunk.values[r])
      }
      out.push(row)
    }
  }
  out
}

///|
/// 字典编码(parquet DICT 编码语义: 重复字符串列 → 去重字典 + 索引).
/// 返回 (字典, 索引数组); 压缩比 = values 数 / 字典大小.
pub fn parquet_dict_encode(
  values : Array[String],
) -> (Array[String], Array[Int32]) {
  let dict : Array[String] = []
  let idx : Map[String, Int32] = Map([], capacity=0)
  let indices : Array[Int32] = Array(capacity=values.length())
  for v in values {
    match idx.get(v) {
      Some(i) => indices.push(i)
      None => {
        idx[v] = dict.length()
        dict.push(v)
        indices.push(dict.length() - 1)
      }
    }
  }
  (dict, indices)
}

///|
/// 序列化为 JSON(持久化交换格式: schema + 行组列块).
/// Int64 以浮点承载, |v| >= 2^53 有精度损失.
pub fn ParquetTable::to_json(self : ParquetTable) -> String {
  let sch : Array[JsonValue] = Array(capacity=self.schema.length())
  for pair in self.schema {
    let m : Map[String, JsonValue] = Map([], capacity=2)
    m["name"] = JsonStr(pair.0)
    m["type"] = JsonStr(parquet_type_tag(pair.1))
    sch.push(JsonObj(m))
  }
  let rgs : Array[JsonValue] = Array(capacity=self.row_groups.length())
  for rg in self.row_groups {
    let cols : Array[JsonValue] = Array(capacity=rg.columns.length())
    for chunk in rg.columns {
      let m : Map[String, JsonValue] = Map([], capacity=2)
      m["name"] = JsonStr(chunk.column)
      let vs : Array[JsonValue] = Array(capacity=chunk.values.length())
      for v in chunk.values {
        vs.push(parquet_value_to_json(v))
      }
      m["values"] = JsonArr(vs)
      cols.push(JsonObj(m))
    }
    let g : Map[String, JsonValue] = Map([], capacity=2)
    g["num_rows"] = JsonNum(rg.num_rows.to_double())
    g["columns"] = JsonArr(cols)
    rgs.push(JsonObj(g))
  }
  let root : Map[String, JsonValue] = Map([], capacity=2)
  root["schema"] = JsonArr(sch)
  root["row_groups"] = JsonArr(rgs)
  json_to_string(JsonObj(root))
}

///|
/// 从 JSON 反序列化(to_json 的逆; 解析失败/schema 为空/类型标签非法
/// 返回 None; 行组缺列或列长不齐则跳过该行组; JSON 值类型与 schema
/// 不符的位置置 NULL).
pub fn parquet_from_json(text : String) -> ParquetTable? raise ParquetError {
  match parse_json(text) {
    None => None
    Some(JsonObj(root)) => {
      let schema : Array[(String, ParquetType)] = []
      match root.get("schema") {
        Some(JsonArr(items)) =>
          for it in items {
            match it {
              JsonObj(m) =>
                match (m.get("name"), m.get("type")) {
                  (Some(JsonStr(n)), Some(JsonStr(t))) =>
                    match parquet_type_from_tag(t) {
                      None => return None
                      Some(pt) => schema.push((n, pt))
                    }
                  _ => ()
                }
              _ => ()
            }
          }
        _ => ()
      }
      if schema.is_empty() {
        return None
      }
      let table = ParquetTable::new(schema)
      match root.get("row_groups") {
        Some(JsonArr(groups)) =>
          for g in groups {
            match g {
              JsonObj(gm) => {
                let cols_json : Map[String, Array[JsonValue]] = Map(
                  [],
                  capacity=0,
                )
                match gm.get("columns") {
                  Some(JsonArr(cs)) =>
                    for c in cs {
                      match c {
                        JsonObj(cm) =>
                          match (cm.get("name"), cm.get("values")) {
                            (Some(JsonStr(n)), Some(JsonArr(vs))) =>
                              cols_json[n] = vs
                            _ => ()
                          }
                        _ => ()
                      }
                    }
                  _ => ()
                }
                let mut nrow = 0
                match gm.get("num_rows") {
                  Some(JsonNum(v)) => nrow = v.to_int()
                  _ => ()
                }
                let chunks : Array[ColumnChunk] = Array(
                  capacity=table.schema.length(),
                )
                let mut ok = true
                for pair in table.schema {
                  match cols_json.get(pair.0) {
                    None => {
                      ok = false
                      break
                    }
                    Some(vs) => {
                      if nrow == 0 {
                        nrow = vs.length()
                      }
                      if vs.length() != nrow {
                        ok = false
                        break
                      }
                      let vals : Array[ParquetValue] = Array(
                        capacity=vs.length(),
                      )
                      for v in vs {
                        vals.push(parquet_value_from_json(v, pair.1))
                      }
                      chunks.push({
                        column: pair.0,
                        ptype: pair.1,
                        values: vals,
                      })
                    }
                  }
                }
                if ok && nrow > 0 {
                  table.row_groups.push({ num_rows: nrow, columns: chunks, })
                }
              }
              _ => ()
            }
          }
        _ => ()
      }
      Some(table)
    }
    _ => None
  }
}

///|
/// 值 → 类型标签(仅非 NULL 值).
fn parquet_type_of(v : ParquetValue) -> ParquetType? {
  match v {
    PvNull => None
    PvBool(_) => Some(PBool)
    PvInt(_) => Some(PInt32)
    PvLong(_) => Some(PInt64)
    PvDouble(_) => Some(PDouble)
    PvStr(_) => Some(PStr)
  }
}

///|
/// 类型 → 标签.
fn parquet_type_tag(t : ParquetType) -> String {
  match t {
    PBool => "bool"
    PInt32 => "int32"
    PInt64 => "int64"
    PDouble => "double"
    PStr => "string"
  }
}

///|
/// 标签 → 类型.
fn parquet_type_from_tag(s : String) -> ParquetType? {
  if s == "bool" {
    Some(PBool)
  } else if s == "int32" {
    Some(PInt32)
  } else if s == "int64" {
    Some(PInt64)
  } else if s == "double" {
    Some(PDouble)
  } else if s == "string" {
    Some(PStr)
  } else {
    None
  }
}

///|
/// schema 列名 → 下标.
fn parquet_schema_index(t : ParquetTable) -> Map[String, Int32] {
  let idx : Map[String, Int32] = Map([], capacity=t.schema.length())
  for i in 0..<t.schema.length() {
    idx[t.schema[i].0] = i
  }
  idx
}

///|
/// 同类型全序比较(Bool 按 false < true, 字符串按字典序);
/// NULL 或类型不同返回 None.
fn parquet_value_cmp(a : ParquetValue, b : ParquetValue) -> Int32? {
  match (a, b) {
    (PvBool(x), PvBool(y)) => Some(if x == y { 0 } else if x { 1 } else { -1 })
    (PvInt(x), PvInt(y)) =>
      Some(if x == y { 0 } else if x < y { -1 } else { 1 })
    (PvLong(x), PvLong(y)) =>
      Some(if x == y { 0 } else if x < y { -1 } else { 1 })
    (PvDouble(x), PvDouble(y)) =>
      Some(if x == y { 0 } else if x < y { -1 } else { 1 })
    (PvStr(x), PvStr(y)) =>
      Some(if x == y { 0 } else if x < y { -1 } else { 1 })
    _ => None
  }
}

///|
/// NULL 不满足任何比较谓词; 类型不匹配同样不命中.
fn parquet_value_matches(v : ParquetValue, f : ParquetFilter) -> Bool {
  match parquet_value_cmp(v, f.value) {
    None => false
    Some(c) =>
      match f.op {
        OpEq => c == 0
        OpNe => c != 0
        OpLt => c < 0
        OpLe => c <= 0
        OpGt => c > 0
        OpGe => c >= 0
      }
  }
}

///|
/// 统计区间 [min, max] 与谓词是否可能相交(不可能 → 整组跳过).
/// 无统计(全 NULL/Bool 列)时保守返回 true.
fn parquet_range_may_match(
  min_v : ParquetValue?,
  max_v : ParquetValue?,
  f : ParquetFilter,
) -> Bool {
  match (min_v, max_v) {
    (Some(lo), Some(hi)) => {
      let lo_c = match parquet_value_cmp(lo, f.value) {
        None => return true
        Some(c) => c
      }
      let hi_c = match parquet_value_cmp(hi, f.value) {
        None => return true
        Some(c) => c
      }
      match f.op {
        OpEq => lo_c <= 0 && hi_c >= 0
        OpNe => !(lo_c == 0 && hi_c == 0)
        OpLt => lo_c < 0
        OpLe => lo_c <= 0
        OpGt => hi_c > 0
        OpGe => hi_c >= 0
      }
    }
    _ => true
  }
}

///|
/// 列块统计(遍历值; min/max 忽略 NULL).
fn parquet_chunk_stats(chunk : ColumnChunk) -> ParquetStats {
  let mut min_v : ParquetValue? = None
  let mut max_v : ParquetValue? = None
  let mut nulls = 0
  for v in chunk.values {
    match v {
      PvNull => nulls = nulls + 1
      _ => {
        min_v = parquet_min_update(min_v, v)
        max_v = parquet_max_update(max_v, v)
      }
    }
  }
  { min: min_v, max: max_v, null_count: nulls, }
}

///|
fn parquet_min_update(cur : ParquetValue?, v : ParquetValue) -> ParquetValue? {
  match cur {
    None => Some(v)
    Some(c) =>
      match parquet_value_cmp(v, c) {
        Some(k) => if k < 0 { Some(v) } else { Some(c) }
        None => Some(c)
      }
  }
}

///|
fn parquet_max_update(cur : ParquetValue?, v : ParquetValue) -> ParquetValue? {
  match cur {
    None => Some(v)
    Some(c) =>
      match parquet_value_cmp(v, c) {
        Some(k) => if k > 0 { Some(v) } else { Some(c) }
        None => Some(c)
      }
  }
}

///|
fn parquet_value_to_json(v : ParquetValue) -> JsonValue {
  match v {
    PvNull => JsonNull
    PvBool(b) => JsonBool(b)
    PvInt(i) => JsonNum(i.to_double())
    PvLong(l) => JsonNum(l.to_double())
    PvDouble(d) => JsonNum(d)
    PvStr(s) => JsonStr(s)
  }
}

///|
/// JSON 值按 schema 类型转换(不符置 NULL).
fn parquet_value_from_json(v : JsonValue, t : ParquetType) -> ParquetValue {
  match v {
    JsonNull => PvNull
    JsonBool(b) => if t is PBool { PvBool(b) } else { PvNull }
    JsonNum(d) =>
      match t {
        PInt32 => PvInt(d.to_int())
        PInt64 => PvLong(d.to_int64())
        PDouble => PvDouble(d)
        PBool => PvNull
        PStr => PvNull
      }
    JsonStr(s) => if t is PStr { PvStr(s) } else { PvNull }
    _ => PvNull
  }
}
