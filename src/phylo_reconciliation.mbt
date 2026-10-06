// =========================================================================
// Phylogenetic Reconciliation: LCA Mapping (ete3)
//   Given a gene tree and a species tree, map each gene tree node to the
//   species tree node that is the most recent common ancestor (LCA) of all
//   species in the gene node's subtree.
//   This is the foundation for duplication/loss inference.
// =========================================================================

///|
/// Result of LCA mapping for one gene tree node.
pub(all) struct LcaMapResult {
  gene_leaves : String       // sorted gene leaf names in this subtree
  species_leaves : String    // sorted species leaf names under the LCA node
  is_duplication : Bool      // true if LCA map equals a child's LCA map
  species_names : Array[String] // species names found in this subtree
}

///|
/// Extract species name from gene name using a separator.
/// E.g., "Human_001" with sep "_" → "Human".
/// If sep is empty, the gene name itself is the species name.
pub fn phylo_extract_species(gene_name : String, sep : String) -> String {
  if sep == "" {
    return gene_name
  }
  // Find the first occurrence of sep
  let mut idx = -1
  for i in 0..<(gene_name.length() - sep.length() + 1) {
    let mut matched = true
    for j in 0..<sep.length() {
      if gene_name[i + j] != sep[j] {
        matched = false
        break
      }
    }
    if matched {
      idx = i
      break
    }
  }
  if idx > 0 {
    gene_name.substring(start=0, end=idx)
  } else {
    gene_name
  }
}

///|
/// Collect sorted leaf names from a Clade.
fn recon_collect_leaf_names(clade : Clade) -> Array[String] {
  let names : Array[String] = Array(capacity=0)
  recon_collect_leaf_names_rec(clade, names)
  names
}

///|
fn recon_collect_leaf_names_rec(
  clade : Clade,
  out : Array[String]
) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => out.push(n)
      None => out.push("?")
    }
    return
  }
  for ch in clade.clades {
    recon_collect_leaf_names_rec(ch, out)
  }
}

///|
/// Sort an array of strings in place.
fn recon_str_lt(a : String, b : String) -> Bool {
  let la = a.length()
  let lb = b.length()
  let min_len = if la < lb { la } else { lb }
  for k in 0..<min_len {
    if a[k] < b[k] {
      return true
    }
    if a[k] > b[k] {
      return false
    }
  }
  la < lb
}

///|
fn recon_sort_strings(arr : Array[String]) -> Unit {
  for i in 0..<arr.length() {
    let mut j = i + 1
    while j < arr.length() {
      if recon_str_lt(arr[j], arr[i]) {
        let tmp = arr[i]
        arr[i] = arr[j]
        arr[j] = tmp
      }
      j = j + 1
    }
  }
}

///|
/// Join sorted leaf names into a canonical key.
fn recon_leaf_key(names : Array[String]) -> String {
  let sorted = Array(capacity=names.length())
  for n in names {
    sorted.push(n)
  }
  recon_sort_strings(sorted)
  let sb = StringBuilder()
  for n in sorted {
    sb.write_string(n)
    sb.write_char(',')
  }
  sb.to_string()
}

///|
/// Check if a clade's subtree contains all given species names.
fn recon_contains_all(clade : Clade, names : Array[String]) -> Bool {
  let leaf_set : Map[String, Bool] = Map([], capacity=0)
  let leaves = recon_collect_leaf_names(clade)
  for l in leaves {
    leaf_set.set(l, true)
  }
  for n in names {
    if not(leaf_set.contains(n)) {
      return false
    }
  }
  true
}

///|
/// Find the LCA of a set of species names in the species tree.
/// Uses subtree membership: the deepest node whose subtree contains all names.
/// Returns the sorted leaf names of the LCA's subtree, or "" if not found.
fn recon_find_lca_key(
  species_tree : Clade,
  species_names : Array[String]
) -> String {
  if species_names.length() == 0 {
    return ""
  }
  if species_names.length() == 1 {
    return species_names[0] + ","
  }
  // Find the deepest node whose subtree contains all species names
  let lca = recon_find_lca_node(species_tree, species_names)
  match lca {
    Some(node) => {
      let lca_leaves = recon_collect_leaf_names(node)
      recon_leaf_key(lca_leaves)
    }
    None => ""
  }
}

///|
/// Find the deepest node whose subtree contains all given names.
fn recon_find_lca_node(
  clade : Clade,
  names : Array[String]
) -> Clade? {
  if not(recon_contains_all(clade, names)) {
    return None
  }
  // Check if any child contains all names (deeper LCA)
  for ch in clade.clades {
    if recon_contains_all(ch, names) {
      return recon_find_lca_node(ch, names)
    }
  }
  // This is the deepest node containing all names
  Some(clade)
}

///|
/// LCA mapping from gene tree to species tree.
/// Species names are extracted from gene leaf names using `sep`
/// (e.g., "Human_001" with sep "_" → species "Human").
/// If sep is empty, gene leaf names are used directly as species names.
/// Returns one `LcaMapResult` per gene tree node (internal + leaf).
pub fn phylo_lca_mapping(
  gene_tree : Clade,
  species_tree : Clade,
  sep : String
) -> Array[LcaMapResult] {
  let gene_to_species : Map[String, String] = Map([], capacity=0)
  phylo_lca_mapping_explicit(gene_tree, species_tree, gene_to_species, sep)
}

///|
/// LCA mapping with explicit gene-to-species map.
/// If a gene name is in `gene_to_species`, use that mapping.
/// Otherwise, extract using `sep`.
pub fn phylo_lca_mapping_explicit(
  gene_tree : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> Array[LcaMapResult] {
  let results : Array[LcaMapResult] = Array(capacity=0)
  ignore(phylo_lca_mapping_rec(
    gene_tree,
    species_tree,
    gene_to_species,
    sep,
    results
  ))
  results
}

///|
/// Recursive helper for LCA mapping.
/// Returns (gene_leaf_key, species_leaf_key, species_names) for the current node.
fn phylo_lca_mapping_rec(
  node : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String,
  out : Array[LcaMapResult]
) -> (String, String, Array[String]) {
  if node.is_terminal() {
    let gene_name = match node.name {
      Some(n) => n
      None => "?"
    }
    let species_name = match gene_to_species.get(gene_name) {
      Some(s) => s
      None => phylo_extract_species(gene_name, sep)
    }
    let gene_key = gene_name + ","
    let species_key = species_name + ","
    let species_names = [species_name]
    out.push(LcaMapResult::{
      gene_leaves: gene_key,
      species_leaves: species_key,
      is_duplication: false,
      species_names,
    })
    return (gene_key, species_key, species_names)
  }
  // Internal node: recurse into children
  let child_results : Array[(String, String, Array[String])] = Array(
    capacity=0
  )
  for ch in node.clades {
    child_results.push(
      phylo_lca_mapping_rec(ch, species_tree, gene_to_species, sep, out)
    )
  }
  // Collect all gene leaf names and species names
  let all_gene_names = recon_collect_leaf_names(node)
  let gene_key = recon_leaf_key(all_gene_names)
  // Deduplicate species names
  let unique_species : Map[String, Bool] = Map([], capacity=0)
  for cr in child_results {
    for sn in cr.2 {
      unique_species.set(sn, true)
    }
  }
  let dedup_species : Array[String] = Array(capacity=0)
  for (sn, _) in unique_species {
    dedup_species.push(sn)
  }
  // Find LCA in species tree
  let species_key = recon_find_lca_key(species_tree, dedup_species)
  // Check for duplication: LCA map equals any child's LCA map
  let mut is_dup = false
  for cr in child_results {
    if cr.1 == species_key && child_results.length() > 1 {
      is_dup = true
    }
  }
  out.push(LcaMapResult::{
    gene_leaves: gene_key,
    species_leaves: species_key,
    is_duplication: is_dup,
    species_names: dedup_species,
  })
  (gene_key, species_key, dedup_species)
}

///|
/// Count the number of duplications in an LCA mapping result.
pub fn phylo_count_duplications(lca_maps : Array[LcaMapResult]) -> Int {
  let mut count = 0
  for m in lca_maps {
    if m.is_duplication {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of losses implied by an LCA mapping.
/// For each gene tree edge (parent → child):
///   - Compute depth difference between parent's LCA and child's LCA in species tree
///   - Speciation node: losses = diff - 1 (if diff > 0)
///   - Duplication node: losses = diff
pub fn phylo_count_losses(
  gene_tree : Clade,
  species_tree : Clade,
  sep : String
) -> Int {
  let gene_to_species : Map[String, String] = Map([], capacity=0)
  phylo_count_losses_explicit(gene_tree, species_tree, gene_to_species, sep)
}

///|
/// Count losses with explicit gene-to-species map.
pub fn phylo_count_losses_explicit(
  gene_tree : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> Int {
  let (_, _, losses) = phylo_count_losses_rec(
    gene_tree,
    species_tree,
    gene_to_species,
    sep
  )
  losses
}

///|
/// Recursive helper for counting losses.
/// Returns (species_key, species_names, loss_count).
fn phylo_count_losses_rec(
  node : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> (String, Array[String], Int) {
  if node.is_terminal() {
    let gene_name = match node.name {
      Some(n) => n
      None => "?"
    }
    let species_name = match gene_to_species.get(gene_name) {
      Some(s) => s
      None => phylo_extract_species(gene_name, sep)
    }
    return (species_name + ",", [species_name], 0)
  }
  // Internal node: recurse into children
  let child_results : Array[(String, Array[String], Int)] = Array(
    capacity=0
  )
  for ch in node.clades {
    child_results.push(
      phylo_count_losses_rec(ch, species_tree, gene_to_species, sep)
    )
  }
  // Collect all species names
  let unique_species : Map[String, Bool] = Map([], capacity=0)
  for cr in child_results {
    for sn in cr.1 {
      unique_species.set(sn, true)
    }
  }
  let dedup : Array[String] = Array(capacity=0)
  for (sn, _) in unique_species {
    dedup.push(sn)
  }
  let parent_key = recon_find_lca_key(species_tree, dedup)
  // Check if this node is a duplication
  let mut is_dup = false
  for cr in child_results {
    if cr.0 == parent_key && child_results.length() > 1 {
      is_dup = true
    }
  }
  // Count losses on each child edge
  let mut total_losses = 0
  for cr in child_results {
    let child_key = cr.0
    if child_key != parent_key && child_key != "" && parent_key != "" {
      // Count edges between parent LCA and child LCA in species tree
      let p_depth = recon_depth(species_tree, parent_key)
      let c_depth = recon_depth(species_tree, child_key)
      let diff = c_depth - p_depth
      if diff > 0 {
        if is_dup {
          total_losses = total_losses + diff
        } else {
          total_losses = total_losses + diff - 1
        }
      }
    }
    // Add child's accumulated losses
    total_losses = total_losses + cr.2
  }
  (parent_key, dedup, total_losses)
}

///|
/// Parse a leaf key back to names.
fn recon_key_to_names(key : String) -> Array[String] {
  let names : Array[String] = Array(capacity=0)
  let mut start = 0
  for i in 0..<key.length() {
    if key[i] == 44 { // ','
      if i > start {
        names.push(key.substring(start=start, end=i))
      }
      start = i + 1
    }
  }
  names
}

///|
/// Compute the depth of a node (identified by leaf key) in the species tree.
fn recon_depth(species_tree : Clade, leaf_key : String) -> Int {
  recon_depth_rec(species_tree, leaf_key, 0)
}

///|
fn recon_depth_rec(
  clade : Clade,
  leaf_key : String,
  depth : Int
) -> Int {
  let my_leaves = recon_collect_leaf_names(clade)
  let my_key = recon_leaf_key(my_leaves)
  if my_key == leaf_key {
    return depth
  }
  for ch in clade.clades {
    let result = recon_depth_rec(ch, leaf_key, depth + 1)
    if result >= 0 {
      return result
    }
  }
  -1
}
