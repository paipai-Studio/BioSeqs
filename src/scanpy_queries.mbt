/// Scanpy-inspired: queries.mitochondrial_genes, queries.ribo_genes
/// Gene query utilities for identifying mitochondrial, ribosomal, and other gene categories.

///|
/// Identify mitochondrial genes by name pattern.
///
/// Parameters:
/// - `var_names`: gene names
/// - `species`: "human" or "mouse" (default "human")
///
/// Returns boolean mask (true = mitochondrial gene)
pub fn scanpy_queries_mito_genes(
  var_names : Array[String],
  species~ : String = "human",
) -> Array[Bool] {
  let prefixes = if species == "mouse" {
    ["mt-", "mt-", "MT-"]
  } else {
    ["MT-", "mt-"]
  }
  let result = Array(capacity=var_names.length())
  for name in var_names {
    let mut is_mito = false
    for p in prefixes {
      if name.has_prefix(p) {
        is_mito = true
      }
    }
    result.push(is_mito)
  }
  result
}

///|
/// Identify ribosomal genes by name pattern.
/// Ribosomal protein genes start with RPS or RPL (human) or Rps or Rpl (mouse).
///
/// Parameters:
/// - `var_names`: gene names
/// - `species`: "human" or "mouse" (default "human")
pub fn scanpy_queries_ribo_genes(
  var_names : Array[String],
  species~ : String = "human",
) -> Array[Bool] {
  let prefixes = if species == "mouse" {
    ["Rps", "Rpl", "rps", "rpl"]
  } else {
    ["RPS", "RPL", "rps", "rpl"]
  }
  let result = Array(capacity=var_names.length())
  for name in var_names {
    let mut is_ribo = false
    for p in prefixes {
      if name.has_prefix(p) {
        is_ribo = true
      }
    }
    result.push(is_ribo)
  }
  result
}

///|
/// Identify hemoglobin genes.
pub fn scanpy_queries_hemo_genes(
  var_names : Array[String],
  species~ : String = "human",
) -> Array[Bool] {
  let prefixes = if species == "mouse" {
    ["Hba", "Hbb", "Hba-", "Hbb-"]
  } else {
    [
      "HBA", "HBB", "HBA1", "HBA2", "HBB1", "HBD", "HBG1", "HBG2", "HBE1", "HBM",
      "HBQ1", "HBZ",
    ]
  }
  let result = Array(capacity=var_names.length())
  for name in var_names {
    let mut is_hemo = false
    for p in prefixes {
      if name.has_prefix(p) {
        is_hemo = true
      }
    }
    result.push(is_hemo)
  }
  result
}

///|
/// Query gene names by pattern (simple wildcard matching).
/// Supports prefix matching with '*' suffix.
///
/// Parameters:
/// - `var_names`: gene names
/// - `pattern`: pattern (e.g., "MT-*" for all MT- genes)
pub fn scanpy_queries_gene_pattern(
  var_names : Array[String],
  pattern : String,
) -> Array[Bool] {
  let result = Array(capacity=var_names.length())
  // Handle wildcard pattern
  if pattern.has_suffix("*") {
    let prefix = pattern[0:pattern.length() - 1].to_owned()
    for name in var_names {
      result.push(name.has_prefix(prefix))
    }
  } else {
    for name in var_names {
      result.push(name == pattern)
    }
  }
  result
}

///|
/// Get indices of genes matching a boolean mask.
pub fn scanpy_queries_mask_to_indices(mask : Array[Bool]) -> Array[Int32] {
  let result = Array(capacity=0)
  for i in 0..<mask.length() {
    if mask[i] {
      result.push(i)
    }
  }
  result
}
