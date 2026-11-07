/// anndata_typed.mbt — typed / nullable obs & var columns.
///
/// Python anndata stores obs/var as pandas DataFrames with real dtypes:
/// string, integer, double, boolean and categorical columns, including
/// pandas' nullable-integer / nullable-boolean / masked-categorical columns
/// (missing values via `pd.NA`). The MoonBit AnnData stores obs/var as
/// string maps; this module adds a typed column layer on top: columns are
/// written back to the string maps (values stringified) while the dtype,
/// the NA mask and category labels are kept in `uns` metadata so the typed
/// view survives round-trips.

///|
/// Column dtype kinds (mirrors pandas dtypes used by anndata).
pub(all) enum AnnDataColKind {
  KStr // object / string
  KInt // nullable-integer
  KDouble // float64
  KBool // nullable-boolean
  KCat // categorical
}

///|
/// A typed annotation column with optional missing-value mask.
pub struct AnnDataTypedColumn {
  name : String
  kind : AnnDataColKind
  strs : Array[String]
  ints : Array[Int32]
  doubles : Array[Double]
  bools : Array[Bool]
  categories : Array[String] // for KCat
  codes : Array[Int32] // for KCat
  na_mask : Array[Bool] // true = missing (pd.NA)
  nullable : Bool
}

///|
/// String column (non-nullable).
pub fn anndata_col_str(
  name : String,
  values : Array[String],
) -> AnnDataTypedColumn {
  let mask : Array[Bool] = Array(capacity=values.length())
  for _v in values {
    mask.push(false)
  }
  {
    name,
    kind: KStr,
    strs: values.copy(),
    ints: [],
    doubles: [],
    bools: [],
    categories: [],
    codes: [],
    na_mask: mask,
    nullable: false,
  }
}

///|
/// Nullable string column; `na[i] = true` marks a missing entry.
pub fn anndata_col_str_nullable(
  name : String,
  values : Array[String],
  na : Array[Bool],
) -> AnnDataTypedColumn {
  {
    name,
    kind: KStr,
    strs: values.copy(),
    ints: [],
    doubles: [],
    bools: [],
    categories: [],
    codes: [],
    na_mask: na.copy(),
    nullable: true,
  }
}

///|
/// Integer column (non-nullable).
pub fn anndata_col_int(
  name : String,
  values : Array[Int32],
) -> AnnDataTypedColumn {
  let mask : Array[Bool] = Array(capacity=values.length())
  for _v in values {
    mask.push(false)
  }
  {
    name,
    kind: KInt,
    strs: [],
    ints: values.copy(),
    doubles: [],
    bools: [],
    categories: [],
    codes: [],
    na_mask: mask,
    nullable: false,
  }
}

///|
/// Nullable integer column (pandas nullable-integer analogue).
pub fn anndata_col_int_nullable(
  name : String,
  values : Array[Int32],
  na : Array[Bool],
) -> AnnDataTypedColumn {
  {
    name,
    kind: KInt,
    strs: [],
    ints: values.copy(),
    doubles: [],
    bools: [],
    categories: [],
    codes: [],
    na_mask: na.copy(),
    nullable: true,
  }
}

///|
/// Double column (non-nullable).
pub fn anndata_col_double(
  name : String,
  values : Array[Double],
) -> AnnDataTypedColumn {
  let mask : Array[Bool] = Array(capacity=values.length())
  for _v in values {
    mask.push(false)
  }
  {
    name,
    kind: KDouble,
    strs: [],
    ints: [],
    doubles: values.copy(),
    bools: [],
    categories: [],
    codes: [],
    na_mask: mask,
    nullable: false,
  }
}

///|
/// Nullable double column.
pub fn anndata_col_double_nullable(
  name : String,
  values : Array[Double],
  na : Array[Bool],
) -> AnnDataTypedColumn {
  {
    name,
    kind: KDouble,
    strs: [],
    ints: [],
    doubles: values.copy(),
    bools: [],
    categories: [],
    codes: [],
    na_mask: na.copy(),
    nullable: true,
  }
}

///|
/// Boolean column (non-nullable).
pub fn anndata_col_bool(
  name : String,
  values : Array[Bool],
) -> AnnDataTypedColumn {
  let mask : Array[Bool] = Array(capacity=values.length())
  for _v in values {
    mask.push(false)
  }
  {
    name,
    kind: KBool,
    strs: [],
    ints: [],
    doubles: [],
    bools: values.copy(),
    categories: [],
    codes: [],
    na_mask: mask,
    nullable: false,
  }
}

///|
/// Nullable boolean column (pandas nullable-boolean analogue).
pub fn anndata_col_bool_nullable(
  name : String,
  values : Array[Bool],
  na : Array[Bool],
) -> AnnDataTypedColumn {
  {
    name,
    kind: KBool,
    strs: [],
    ints: [],
    doubles: [],
    bools: values.copy(),
    categories: [],
    codes: [],
    na_mask: na.copy(),
    nullable: true,
  }
}

///|
/// Categorical column from explicit categories + codes.
pub fn anndata_col_cat(
  name : String,
  categories : Array[String],
  codes : Array[Int32],
) -> AnnDataTypedColumn {
  let mask : Array[Bool] = Array(capacity=codes.length())
  for c in codes {
    mask.push(c < 0)
  }
  {
    name,
    kind: KCat,
    strs: [],
    ints: [],
    doubles: [],
    bools: [],
    categories: categories.copy(),
    codes: codes.copy(),
    na_mask: mask,
    nullable: false,
  }
}

///|
/// Column length.
pub fn AnnDataTypedColumn::len(self : AnnDataTypedColumn) -> Int32 {
  match self.kind {
    KStr => self.strs.length()
    KInt => self.ints.length()
    KDouble => self.doubles.length()
    KBool => self.bools.length()
    KCat => self.codes.length()
  }
}

///|
/// Dtype name (matches the h5ad encoding-type family).
pub fn AnnDataTypedColumn::kind_name(self : AnnDataTypedColumn) -> String {
  match self.kind {
    KStr => "string"
    KInt => "integer"
    KDouble => "double"
    KBool => "boolean"
    KCat => "categorical"
  }
}

///|
/// True when entry `i` is missing (pd.NA analogue).
pub fn AnnDataTypedColumn::is_na(self : AnnDataTypedColumn, i : Int32) -> Bool {
  if i >= 0 && i < self.na_mask.length() {
    self.na_mask[i]
  } else {
    false
  }
}

///|
/// Number of missing entries.
pub fn AnnDataTypedColumn::n_missing(self : AnnDataTypedColumn) -> Int32 {
  let mut n = 0
  for m in self.na_mask {
    if m {
      n = n + 1
    }
  }
  n
}

///|
/// Categorical value at `i` ("" when missing / out of range).
pub fn AnnDataTypedColumn::cat_value(
  self : AnnDataTypedColumn,
  i : Int32,
) -> String {
  if i >= 0 &&
    i < self.codes.length() &&
    self.codes[i] >= 0 &&
    self.codes[i] < self.categories.length() {
    self.categories[self.codes[i]]
  } else {
    ""
  }
}

///|
/// Stringify the column (missing entries become ""); this is what gets
/// stored in the AnnData string maps.
pub fn AnnDataTypedColumn::to_strings(
  self : AnnDataTypedColumn,
) -> Array[String] {
  let out : Array[String] = Array(capacity=self.len())
  for i = 0; i < self.len(); i = i + 1 {
    if self.is_na(i) {
      out.push("")
    } else {
      match self.kind {
        KStr => out.push(self.strs[i])
        KInt => out.push(self.ints[i].to_string())
        KDouble => out.push(self.doubles[i].to_string())
        KBool => out.push(if self.bools[i] { "true" } else { "false" })
        KCat => out.push(self.cat_value(i))
      }
    }
  }
  out
}

///|
/// Subset the column by indices (keeps dtype, NA mask and categories).
pub fn AnnDataTypedColumn::take(
  self : AnnDataTypedColumn,
  idx : Array[Int32],
) -> AnnDataTypedColumn {
  let strs : Array[String] = Array(capacity=idx.length())
  let ints : Array[Int32] = Array(capacity=idx.length())
  let doubles : Array[Double] = Array(capacity=idx.length())
  let bools : Array[Bool] = Array(capacity=idx.length())
  let codes : Array[Int32] = Array(capacity=idx.length())
  let mask : Array[Bool] = Array(capacity=idx.length())
  for i in idx {
    let valid = i >= 0 && i < self.len()
    mask.push(valid && self.is_na(i))
    match self.kind {
      KStr => strs.push(if valid { self.strs[i] } else { "" })
      KInt => ints.push(if valid { self.ints[i] } else { 0 })
      KDouble => doubles.push(if valid { self.doubles[i] } else { 0.0 })
      KBool => bools.push(valid && self.bools[i])
      KCat => codes.push(if valid { self.codes[i] } else { -1 })
    }
  }
  {
    name: self.name,
    kind: self.kind,
    strs,
    ints,
    doubles,
    bools,
    categories: self.categories.copy(),
    codes,
    na_mask: mask,
    nullable: self.nullable,
  }
}

///|
/// Serialize the column metadata (dtype / NA mask / categories) to uns.
fn typed_col_meta_json(col : AnnDataTypedColumn) -> JsonValue {
  let m : Map[String, JsonValue] = Map([], capacity=4)
  m["kind"] = JsonStr(col.kind_name())
  m["nullable"] = JsonBool(col.nullable)
  let na : Array[JsonValue] = Array(capacity=col.na_mask.length())
  for b in col.na_mask {
    na.push(JsonBool(b))
  }
  m["na"] = JsonArr(na)
  if col.kind is KCat {
    let cats : Array[JsonValue] = Array(capacity=col.categories.length())
    for s in col.categories {
      cats.push(JsonStr(s))
    }
    m["categories"] = JsonArr(cats)
  }
  JsonObj(m)
}

///|
/// Metadata key under which a typed obs column is recorded in uns.
fn typed_obs_meta_key(name : String) -> String {
  "__typed__/obs/" + name
}

///|
/// Metadata key under which a typed var column is recorded in uns.
fn typed_var_meta_key(name : String) -> String {
  "__typed__/var/" + name
}

///|
/// Write a typed column into `obs` (values stringified; dtype + NA mask
/// recorded in uns). Returns `self` for chaining.
pub fn AnnData::set_obs_typed(
  self : AnnData,
  col : AnnDataTypedColumn,
) -> AnnData {
  self.obs[col.name] = col.to_strings()
  self.uns[typed_obs_meta_key(col.name)] = typed_col_meta_json(col)
  self
}

///|
/// Write a typed column into `var`.
pub fn AnnData::set_var_typed(
  self : AnnData,
  col : AnnDataTypedColumn,
) -> AnnData {
  self.variables[col.name] = col.to_strings()
  self.uns[typed_var_meta_key(col.name)] = typed_col_meta_json(col)
  self
}

///|
/// Rebuild a typed column from its stringified obs storage + uns metadata.
/// Without metadata the column is recovered as plain strings (KStr).
pub fn AnnData::get_obs_typed(
  self : AnnData,
  name : String,
) -> AnnDataTypedColumn? {
  match self.obs.get(name) {
    Some(strs) =>
      Some(
        typed_from_strings(name, strs, self.uns.get(typed_obs_meta_key(name))),
      )
    None => None
  }
}

///|
/// Rebuild a typed var column.
pub fn AnnData::get_var_typed(
  self : AnnData,
  name : String,
) -> AnnDataTypedColumn? {
  match self.variables.get(name) {
    Some(strs) =>
      Some(
        typed_from_strings(name, strs, self.uns.get(typed_var_meta_key(name))),
      )
    None => None
  }
}

///|
/// Reconstruct a typed column from stringified values and optional
/// uns metadata.
fn typed_from_strings(
  name : String,
  strs : Array[String],
  meta : JsonValue?,
) -> AnnDataTypedColumn {
  let mut kind = KStr
  let mut nullable = false
  let meta_na : Array[Bool] = Array(capacity=0)
  let categories : Array[String] = Array(capacity=0)
  match meta {
    Some(JsonObj(m)) => {
      match m.get("kind") {
        Some(JsonStr(s)) =>
          if s == "integer" {
            kind = KInt
          } else if s == "double" {
            kind = KDouble
          } else if s == "boolean" {
            kind = KBool
          } else if s == "categorical" {
            kind = KCat
          }
        _ => ()
      }
      match m.get("nullable") {
        Some(JsonBool(b)) => nullable = b
        _ => ()
      }
      match m.get("na") {
        Some(JsonArr(arr)) =>
          for v in arr {
            match v {
              JsonBool(b) => meta_na.push(b)
              _ => meta_na.push(false)
            }
          }
        _ => ()
      }
      match m.get("categories") {
        Some(JsonArr(arr)) =>
          for v in arr {
            match v {
              JsonStr(s) => categories.push(s)
              _ => categories.push("")
            }
          }
        _ => ()
      }
    }
    _ => ()
  }
  let n = strs.length()
  let mask : Array[Bool] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    mask.push(nullable && i < meta_na.length() && meta_na[i])
  }
  let base = {
    name,
    kind,
    strs: strs.copy(),
    ints: Array(capacity=0),
    doubles: Array(capacity=0),
    bools: Array(capacity=0),
    categories,
    codes: Array(capacity=0),
    na_mask: mask,
    nullable,
  }
  if kind is KInt {
    let vals : Array[Int32] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      if base.is_na(i) {
        vals.push(0)
      } else {
        vals.push(parse_int(strs[i]))
      }
    }
    {
      name: base.name,
      kind: base.kind,
      strs: [],
      ints: vals,
      doubles: [],
      bools: [],
      categories: base.categories,
      codes: [],
      na_mask: base.na_mask,
      nullable: base.nullable,
    }
  } else if kind is KDouble {
    let vals : Array[Double] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      if base.is_na(i) {
        vals.push(0.0)
      } else {
        vals.push(parse_double(strs[i]).unwrap_or(0.0))
      }
    }
    {
      name: base.name,
      kind: base.kind,
      strs: [],
      ints: [],
      doubles: vals,
      bools: [],
      categories: base.categories,
      codes: [],
      na_mask: base.na_mask,
      nullable: base.nullable,
    }
  } else if kind is KBool {
    let vals : Array[Bool] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      vals.push(not(base.is_na(i)) && strs[i] == "true")
    }
    {
      name: base.name,
      kind: base.kind,
      strs: [],
      ints: [],
      doubles: [],
      bools: vals,
      categories: base.categories,
      codes: [],
      na_mask: base.na_mask,
      nullable: base.nullable,
    }
  } else if kind is KCat {
    let codes : Array[Int32] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      if base.is_na(i) {
        codes.push(-1)
      } else {
        let mut code = -1
        for k = 0; k < base.categories.length(); k = k + 1 {
          if base.categories[k] == strs[i] {
            code = k
            break
          }
        }
        codes.push(code)
      }
    }
    {
      name: base.name,
      kind: base.kind,
      strs: [],
      ints: [],
      doubles: [],
      bools: [],
      categories: base.categories,
      codes,
      na_mask: base.na_mask,
      nullable: base.nullable,
    }
  } else {
    base
  }
}
