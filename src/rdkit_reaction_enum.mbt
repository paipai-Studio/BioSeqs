/// rdkit_reaction_enum.mbt — reaction SMARTS parsing + application, combinatorial library enumeration
/// Items 13-14

///|
/// A reaction definition: reactant templates -> product templates
pub struct RdkitReaction {
  reactants : Array[String] // SMARTS patterns for reactants
  products : Array[String] // SMARTS patterns for products
  name : String
}

///|
/// Parse a reaction SMARTS string (reactants>>products format).
pub fn rdkit_parse_reaction(rxn_smarts : String) -> RdkitReaction {
  // Split on ">>"
  let parts = rdkit_smarts_split(rxn_smarts, ">>")
  let reactant_str = if parts.length() > 0 { parts[0] } else { "" }
  let product_str = if parts.length() > 1 { parts[1] } else { "" }
  // Split reactants on "."
  let reactants = rdkit_smarts_split(reactant_str, ".")
  // Split products on "."
  let products = rdkit_smarts_split(product_str, ".")
  RdkitReaction::{ reactants, products, name: "", }
}

///|
fn rdkit_smarts_split(s : String, sep : String) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let sep_len = sep.length()
  let mut start = 0
  let s_len = s.length()
  let mut i = 0
  while i <= s_len - sep_len {
    let mut is_match = true
    for j = 0; j < sep_len; j = j + 1 {
      if s[i + j] != sep[j] {
        is_match = false
      }
    }
    if is_match {
      result.push(s[start:i].to_owned())
      start = i + sep_len
      i = i + sep_len
    } else {
      i = i + 1
    }
  }
  if start <= s_len {
    result.push(s[start:s_len].to_owned())
  }
  // Filter empty
  let filtered : Array[String] = Array(capacity=0)
  for k = 0; k < result.length(); k = k + 1 {
    if result[k].length() > 0 {
      filtered.push(result[k])
    }
  }
  filtered
}

///|
/// Apply a reaction to a set of reactant molecules.
/// Returns the product SMILES (simplified: combines matched atoms).
pub fn rdkit_apply_reaction(
  rxn : RdkitReaction,
  reactant_smiles : Array[String],
) -> Array[String] {
  let products : Array[String] = Array(capacity=0)
  // Check that we have enough reactants
  if reactant_smiles.length() < rxn.reactants.length() {
    return products
  }
  // For each reactant, check if it matches the template
  let mut all_is_match = true
  for i = 0; i < rxn.reactants.length(); i = i + 1 {
    if i < reactant_smiles.length() {
      // Simplified: check substructure match
      if !rdkit_rxn_template_match(reactant_smiles[i], rxn.reactants[i]) {
        all_is_match = false
      }
    }
  }
  if !all_is_match {
    return products
  }
  // Generate product by combining reactants
  // Simplified: concatenate non-overlapping parts
  let combined = rdkit_rxn_combine(reactant_smiles, rxn)
  products.push(combined)
  products
}

///|
fn rdkit_rxn_template_match(smiles : String, template : String) -> Bool {
  // Simplified: if template is generic (contains *), accept
  if template.contains("*") {
    return true
  }
  // Try substructure match
  rdkit_has_substruct(smiles, template)
}

///|
fn rdkit_rxn_combine(reactants : Array[String], rxn : RdkitReaction) -> String {
  // Simplified: join with "." if multiple products
  if rxn.products.length() == 0 {
    // Default: join reactants
    let sb = StringBuilder()
    for i = 0; i < reactants.length(); i = i + 1 {
      if i > 0 {
        let _ = sb.write(".")
      }
      let _ = sb.write(reactants[i])
    }
    return sb.to_string()
  }
  // Use product template as scaffold
  let sb = StringBuilder()
  for i = 0; i < rxn.products.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write(".")
    }
    let _ = sb.write(rxn.products[i])
  }
  sb.to_string()
}

///|
/// Reaction application result with atom mapping
pub struct RdkitRxnResult {
  products : Array[String]
  reactant_indices : Array[Int32] // which reactants were used
  success : Bool
}

///|
/// Apply reaction with detailed result.
pub fn rdkit_apply_reaction_detailed(
  rxn : RdkitReaction,
  reactant_smiles : Array[String],
) -> RdkitRxnResult {
  let products = rdkit_apply_reaction(rxn, reactant_smiles)
  let indices : Array[Int32] = Array(capacity=0)
  for i = 0; i < rxn.reactants.length(); i = i + 1 {
    indices.push(i)
  }
  RdkitRxnResult::{
    products,
    reactant_indices: indices,
    success: products.length() > 0,
  }
}

///|
/// Combinatorial library enumeration: apply a reaction to all combinations of reactants.
pub fn rdkit_enumerate_library(
  rxn : RdkitReaction,
  reactant_pools : Array[Array[String]], // one pool per reactant slot
  max_products~ : Int32 = 1000,
) -> Array[String] {
  let products : Array[String] = Array(capacity=0)
  if reactant_pools.length() == 0 {
    return products
  }
  // Generate all combinations
  let indices : Array[Int32] = Array(capacity=reactant_pools.length())
  for i = 0; i < reactant_pools.length(); i = i + 1 {
    indices.push(0)
  }
  rdkit_enum_backtrack(0, indices, reactant_pools, rxn, products, max_products)
  products
}

///|
fn rdkit_enum_backtrack(
  depth : Int32,
  indices : Array[Int32],
  pools : Array[Array[String]],
  rxn : RdkitReaction,
  products : Array[String],
  max_products : Int32,
) -> Unit {
  if products.length() >= max_products {
    return
  }
  if depth == pools.length() {
    // Build reactant list from current indices
    let reactants : Array[String] = Array(capacity=depth)
    for i = 0; i < depth; i = i + 1 {
      reactants.push(pools[i][indices[i]])
    }
    let result = rdkit_apply_reaction(rxn, reactants)
    for i = 0; i < result.length(); i = i + 1 {
      if products.length() < max_products {
        products.push(result[i])
      }
    }
    return
  }
  for i = 0; i < pools[depth].length(); i = i + 1 {
    indices[depth] = i
    rdkit_enum_backtrack(depth + 1, indices, pools, rxn, products, max_products)
    if products.length() >= max_products {
      return
    }
  }
}

///|
/// Count unique products from enumeration (deduplicate).
pub fn rdkit_count_unique_products(products : Array[String]) -> Int32 {
  let seen : Set[String] = Set([])
  for i = 0; i < products.length(); i = i + 1 {
    seen.add(products[i])
  }
  seen.length()
}

///|
/// Enumerate library with deduplication.
pub fn rdkit_enumerate_library_unique(
  rxn : RdkitReaction,
  reactant_pools : Array[Array[String]],
  max_products~ : Int32 = 1000,
) -> Array[String] {
  let raw = rdkit_enumerate_library(rxn, reactant_pools, max_products~)
  let seen : Set[String] = Set([])
  let unique : Array[String] = Array(capacity=0)
  for i = 0; i < raw.length(); i = i + 1 {
    if !seen.contains(raw[i]) {
      seen.add(raw[i])
      unique.push(raw[i])
    }
  }
  unique
}
