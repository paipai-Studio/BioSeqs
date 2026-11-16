/// AnnData categorical encoding for obs/var columns.
///
/// Encodes string columns as categorical (categories + codes),
/// saving memory and enabling groupby operations.

///|
/// Categorical encoding result
pub struct AnnDataCategorical {
  categories : Array[String] // unique category labels
  codes : Array[Int32] // integer codes (index into categories)
}

///|
/// Encode a string array as categorical.
/// Returns categories (unique values in order of first appearance)
/// and codes (integer index for each element).
pub fn anndata_encode_categorical(values : Array[String]) -> AnnDataCategorical {
  let categories : Array[String] = Array(capacity=0)
  let cat_map : Map[String, Int32] = Map([], capacity=values.length())
  let codes : Array[Int32] = Array(capacity=values.length())
  for s in values {
    match cat_map.get(s) {
      Some(idx) => codes.push(idx)
      None => {
        cat_map[s] = categories.length()
        categories.push(s)
        codes.push(categories.length() - 1)
      }
    }
  }
  { categories, codes, }
}

///|
/// Decode categorical codes back to strings.
pub fn anndata_decode_categorical(cat : AnnDataCategorical) -> Array[String] {
  let out : Array[String] = Array(capacity=cat.codes.length())
  for c in cat.codes {
    if c >= 0 && c < cat.categories.length() {
      out.push(cat.categories[c])
    } else {
      out.push("")
    }
  }
  out
}

///|
/// Number of categories.
pub fn AnnDataCategorical::n_categories(self : AnnDataCategorical) -> Int32 {
  self.categories.length()
}

///|
/// Get category name by code.
pub fn AnnDataCategorical::category(
  self : AnnDataCategorical,
  code : Int32,
) -> String {
  if code >= 0 && code < self.categories.length() {
    self.categories[code]
  } else {
    ""
  }
}

///|
/// Get code for a category name.
pub fn AnnDataCategorical::code(
  self : AnnDataCategorical,
  name : String,
) -> Int32 {
  for i = 0; i < self.categories.length(); i = i + 1 {
    if self.categories[i] == name {
      return i
    }
  }
  -1
}

///|
/// Count elements in each category.
pub fn AnnDataCategorical::counts(self : AnnDataCategorical) -> Array[Int32] {
  let out : Array[Int32] = Array::make(self.categories.length(), 0)
  for c in self.codes {
    if c >= 0 && c < out.length() {
      out[c] = out[c] + 1
    }
  }
  out
}

///|
/// Encode an obs column as categorical and store in uns.
pub fn anndata_obs_to_categorical(
  adata : AnnData,
  key : String,
) -> AnnDataCategorical? {
  match adata.obs_column(key) {
    Some(col) => Some(anndata_encode_categorical(col))
    None => None
  }
}

///|
/// Encode a var column as categorical.
pub fn anndata_var_to_categorical(
  adata : AnnData,
  key : String,
) -> AnnDataCategorical? {
  match adata.var_column(key) {
    Some(col) => Some(anndata_encode_categorical(col))
    None => None
  }
}

///|
/// Replace an obs column with its categorical codes (as strings).
/// Stores the category mapping in uns for later decoding.
pub fn anndata_obs_factorize(adata : AnnData, key : String) -> Unit {
  match adata.obs_column(key) {
    Some(col) => {
      let cat = anndata_encode_categorical(col)
      let code_strs : Array[String] = Array(capacity=cat.codes.length())
      for c in cat.codes {
        code_strs.push(c.to_string())
      }
      ignore(adata.set_obs(key, code_strs))
      // Store categories in uns
      let cat_arr : Array[JsonValue] = Array(capacity=cat.categories.length())
      for s in cat.categories {
        cat_arr.push(JsonStr(s))
      }
      ignore(adata.set_uns(key + "_categories", JsonArr(cat_arr)))
    }
    None => ()
  }
}

///|
/// Replace a var column with its categorical codes.
pub fn anndata_var_factorize(adata : AnnData, key : String) -> Unit {
  match adata.var_column(key) {
    Some(col) => {
      let cat = anndata_encode_categorical(col)
      let code_strs : Array[String] = Array(capacity=cat.codes.length())
      for c in cat.codes {
        code_strs.push(c.to_string())
      }
      ignore(adata.set_var(key, code_strs))
      let cat_arr : Array[JsonValue] = Array(capacity=cat.categories.length())
      for s in cat.categories {
        cat_arr.push(JsonStr(s))
      }
      ignore(adata.set_uns(key + "_categories", JsonArr(cat_arr)))
    }
    None => ()
  }
}

///|
/// Convert all string obs/var columns with repeated values into
/// categorical (code) storage, mirroring `AnnData.strings_to_categoricals`.
/// Each converted column is replaced by its integer codes (as strings) and
/// the category labels are stored in `uns` under `<key>_categories`.
/// Returns the number of columns converted.
pub fn anndata_strings_to_categoricals(adata : AnnData) -> Int32 {
  let mut converted = 0
  // obs columns (snapshot keys: factorize mutates the maps)
  let obs_keys = adata.obs.keys().collect()
  for k in obs_keys {
    match adata.obs.get(k) {
      Some(col) =>
        if col.length() > 0 {
          let cat = anndata_encode_categorical(col)
          // Only convert when there are actual repeats (true categorical)
          if cat.categories.length() < col.length() {
            anndata_obs_factorize(adata, k)
            converted = converted + 1
          }
        }
      None => ()
    }
  }
  // var columns
  let var_keys = adata.variables.keys().collect()
  for k in var_keys {
    match adata.variables.get(k) {
      Some(col) =>
        if col.length() > 0 {
          let cat = anndata_encode_categorical(col)
          if cat.categories.length() < col.length() {
            anndata_var_factorize(adata, k)
            converted = converted + 1
          }
        }
      None => ()
    }
  }
  converted
}

// ============================================================================
// AnnData.rename_categories — 官方 API
// ============================================================================

///|
/// Rename categories of an obs or var column.
///
/// `mapping` 是 Array[(old_name, new_name)] 映射。
///
/// 如果该列之前被 factorize 过（categories 存在 uns[key + "_categories"]），
/// 则重映射 obs 列里的 codes 并更新 uns；否则重建列值。
///
/// 返回 true 表示成功。
pub fn AnnData::rename_categories(
  self : AnnData,
  key : String,
  mapping : Array[(String, String)],
  axis~ : String = "obs",
) -> Bool {
  let uns_cat_key = key + "_categories"

  // 1. 查原 obs/var 列
  let codes_str : Array[String] = if axis == "var" {
    match self.var_column(key) {
      Some(c) => c
      None => return false
    }
  } else {
    match self.obs_column(key) {
      Some(c) => c
      None => return false
    }
  }

  // 2. 查 uns 里存的 categories（如果之前 factorize 过）
  let old_cats : Array[String] = match
    anndata_uns_get_str_array(self, uns_cat_key) {
    Some(a) => a
    None => {
      // 没有预存 categories → 从 mapping 的 old 部分构造
      let out : Array[String] = Array(capacity=mapping.length())
      for pair in mapping {
        let (old_name, _new_name) = pair
        out.push(old_name)
      }
      out
    }
  }

  // 3. 建立 new_cats：对每个 old_cat，找 mapping 里的 new（没映射就保持原名）
  let new_cats : Array[String] = Array(capacity=old_cats.length())
  for c in old_cats {
    let mut mapped = c
    for pair in mapping {
      let (old_n, new_n) = pair
      if old_n == c {
        mapped = new_n
        break
      }
    }
    new_cats.push(mapped)
  }

  // 4. 重建列值：codes 字符串 → 按 new_cats 解码 → category 名
  let new_col : Array[String] = Array(capacity=codes_str.length())
  for cs in codes_str {
    let name = match parse_int_string_local(cs) {
      Some(code_idx) =>
        if code_idx >= 0 && code_idx < new_cats.length() {
          new_cats[code_idx]
        } else {
          cs
        }
      None => cs
    }
    new_col.push(name)
  }

  // 5. 写回 obs/var
  if axis == "var" {
    ignore(self.set_var(key, new_col))
  } else {
    ignore(self.set_obs(key, new_col))
  }

  // 6. 更新 uns 里的 categories
  let cat_arr : Array[JsonValue] = Array(capacity=new_cats.length())
  for s in new_cats {
    cat_arr.push(JsonStr(s))
  }
  ignore(self.set_uns(uns_cat_key, JsonArr(cat_arr)))

  true
}

///|
/// 直接用新 categories 列表替换整个 categories。
/// `new_categories` 长度必须匹配原 categories 数。
pub fn AnnData::rename_categories_list(
  self : AnnData,
  key : String,
  new_categories : Array[String],
  axis~ : String = "obs",
) -> Bool {
  let uns_cat_key = key + "_categories"
  match anndata_uns_get_str_array(self, uns_cat_key) {
    None => false
    Some(old) =>
      if old.length() != new_categories.length() {
        false
      } else {
        let mapping : Array[(String, String)] = Array(capacity=old.length())
        for i = 0; i < old.length(); i = i + 1 {
          mapping.push((old[i], new_categories[i]))
        }
        self.rename_categories(key, mapping, axis~)
      }
  }
}

///|
/// AnnDataCategorical 重命名方法
pub fn AnnDataCategorical::rename(
  self : AnnDataCategorical,
  mapping : Array[(String, String)],
) -> AnnDataCategorical {
  let new_cats : Array[String] = Array(capacity=self.categories.length())
  for c in self.categories {
    let mut mapped = c
    for pair in mapping {
      let (old_name, new_name) = pair
      if old_name == c {
        mapped = new_name
        break
      }
    }
    new_cats.push(mapped)
  }
  anndata_categorical_make(new_cats, self.codes)
}

///|
/// Helper: 从两个 Array 构造 AnnDataCategorical（文件内所以能 access literal）
/// MoonBit struct literal 只支持 shorthand: { field1, field2, }
fn anndata_categorical_make(
  categories : Array[String],
  codes : Array[Int32],
) -> AnnDataCategorical {
  { categories, codes, }
}

///|
/// 本地 parse_int_string（因为 anndata.mbt 里的是私有）
fn parse_int_string_local(s : String) -> Int32? {
  if s.length() == 0 {
    return None
  }
  let first_i : Int32 = s.unsafe_get(0).unsafe_to_char().to_int()
  let mut is_num = (first_i >= 48 && first_i <= 57) || first_i == 45
  let mut i : Int32 = 1
  while i < s.length() && is_num {
    let ci : Int32 = s.unsafe_get(i).unsafe_to_char().to_int()
    if !(ci >= 48 && ci <= 57) {
      is_num = false
    }
    i = i + 1
  }
  if !is_num {
    return None
  }
  let neg = first_i == 45
  let start = if neg { 1 } else { 0 }
  let mut val : Int32 = 0
  let mut j = start
  while j < s.length() {
    let ci : Int32 = s.unsafe_get(j).unsafe_to_char().to_int()
    val = val * 10 + (ci - 48)
    j = j + 1
  }
  if neg {
    val = -val
  }
  Some(val)
}
