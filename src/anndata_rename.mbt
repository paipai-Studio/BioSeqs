/// AnnData rename columns and make unique names.
///
/// `adata.obs.rename(columns=...)` and
/// `adata.obs_names_make_unique()` equivalents.

///|
/// Rename an obs column key.
pub fn anndata_rename_obs_key(
  adata : AnnData,
  old_key : String,
  new_key : String,
) -> Unit {
  match adata.obs_column(old_key) {
    Some(col) => {
      ignore(adata.set_obs(new_key, col))
      adata.obs.remove(old_key)
    }
    None => ()
  }
}

///|
/// Rename a var column key.
pub fn anndata_rename_var_key(
  adata : AnnData,
  old_key : String,
  new_key : String,
) -> Unit {
  match adata.var_column(old_key) {
    Some(col) => {
      ignore(adata.set_var(new_key, col))
      adata.variables.remove(old_key)
    }
    None => ()
  }
}

///|
/// Rename multiple obs columns at once.
/// pairs: array of (old_key, new_key).
pub fn anndata_rename_obs_keys(
  adata : AnnData,
  pairs : Array[(String, String)],
) -> Unit {
  for p in pairs {
    anndata_rename_obs_key(adata, p.0, p.1)
  }
}

///|
/// Rename multiple var columns at once.
pub fn anndata_rename_var_keys(
  adata : AnnData,
  pairs : Array[(String, String)],
) -> Unit {
  for p in pairs {
    anndata_rename_var_key(adata, p.0, p.1)
  }
}

///|
/// Make obs_names unique by appending suffix to duplicates.
/// e.g. ["cell", "cell"] -> ["cell", "cell-1"]
pub fn anndata_obs_names_make_unique(adata : AnnData) -> Unit {
  let names = adata.get_obs_names()
  let unique = anndata_make_unique_names(names)
  // Update obs_names by setting a new array
  // Since obs_names is a field, we need to clear and repopulate
  for i = 0; i < names.length(); i = i + 1 {
    if i < unique.length() {
      names[i] = unique[i]
    }
  }
}

///|
/// Make var_names unique by appending suffix to duplicates.
pub fn anndata_var_names_make_unique(adata : AnnData) -> Unit {
  let names = adata.get_var_names()
  let unique = anndata_make_unique_names(names)
  for i = 0; i < names.length(); i = i + 1 {
    if i < unique.length() {
      names[i] = unique[i]
    }
  }
}

///|
/// Make a string array unique by appending "-N" to duplicates.
pub fn anndata_make_unique_names(names : Array[String]) -> Array[String] {
  let out : Array[String] = Array(capacity=names.length())
  let counts : Map[String, Int32] = Map([], capacity=names.length())
  for s in names {
    match counts.get(s) {
      None => {
        out.push(s)
        counts[s] = 1
      }
      Some(c) => {
        let new_name = s + "-" + c.to_string()
        out.push(new_name)
        counts[s] = c + 1
      }
    }
  }
  out
}

///|
/// Check if obs_names are unique.
pub fn anndata_obs_names_is_unique(adata : AnnData) -> Bool {
  let names = adata.get_obs_names()
  let seen : Map[String, Bool] = Map([], capacity=names.length())
  for s in names {
    match seen.get(s) {
      Some(_) => return false
      None => seen[s] = true
    }
  }
  true
}

///|
/// Check if var_names are unique.
pub fn anndata_var_names_is_unique(adata : AnnData) -> Bool {
  let names = adata.get_var_names()
  let seen : Map[String, Bool] = Map([], capacity=names.length())
  for s in names {
    match seen.get(s) {
      Some(_) => return false
      None => seen[s] = true
    }
  }
  true
}

///|
/// Set obs_names to a new array.
pub fn anndata_set_obs_names(adata : AnnData, names : Array[String]) -> Unit {
  let current = adata.get_obs_names()
  for i = 0; i < current.length() && i < names.length(); i = i + 1 {
    current[i] = names[i]
  }
}

///|
/// Set var_names to a new array.
pub fn anndata_set_var_names(adata : AnnData, names : Array[String]) -> Unit {
  let current = adata.get_var_names()
  for i = 0; i < current.length() && i < names.length(); i = i + 1 {
    current[i] = names[i]
  }
}
