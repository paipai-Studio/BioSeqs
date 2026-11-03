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
