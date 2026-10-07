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
  recon_sort_strings(dedup_species)
  // Find LCA in species tree
  let species_key = recon_find_lca_key(species_tree, dedup_species)
  // Check for duplication: LCA map equals any child's LCA map
  let mut is_dup = false
  for cr in child_results {
    if cr.1 == species_key && child_results.length() > 1 && species_key != "" {
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
  recon_sort_strings(dedup)
  let parent_key = recon_find_lca_key(species_tree, dedup)
  // Check if this node is a duplication
  let mut is_dup = false
  for cr in child_results {
    if cr.0 == parent_key && child_results.length() > 1 && parent_key != "" {
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

// =========================================================================
// Dup/Loss Inference: comprehensive reconciliation result (ete3)
//   Combines LCA mapping, duplication detection, and loss counting into
//   a structured result with per-node information and summary output.
// =========================================================================

///|
/// Per-node reconciliation information.
pub(all) struct NodeReconciliation {
  gene_leaves : String         // sorted gene leaf names in this subtree
  species_leaves : String      // sorted species leaf names under LCA
  is_duplication : Bool        // true if this node is a duplication
  is_speciation : Bool         // true if this node is a speciation
  edge_losses : Int            // losses on edges from this node to its children
  species_names : Array[String]
}

///|
/// Comprehensive reconciliation result.
pub(all) struct ReconciliationResult {
  duplication_count : Int
  loss_count : Int
  speciation_count : Int
  total_cost : Int             // duplication_count + loss_count
  node_count : Int             // total gene tree nodes
  leaf_count : Int             // gene tree leaves
  nodes : Array[NodeReconciliation]
}

///|
/// Reconcile a gene tree with a species tree.
/// Returns a comprehensive result with per-node duplication/loss information.
/// Species names are extracted from gene leaf names using `sep`.
pub fn phylo_reconcile(
  gene_tree : Clade,
  species_tree : Clade,
  sep : String
) -> ReconciliationResult {
  let gene_to_species : Map[String, String] = Map([], capacity=0)
  phylo_reconcile_explicit(gene_tree, species_tree, gene_to_species, sep)
}

///|
/// Reconcile with explicit gene-to-species map.
pub fn phylo_reconcile_explicit(
  gene_tree : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> ReconciliationResult {
  // Build per-node reconciliation with edge losses
  let nodes : Array[NodeReconciliation] = Array(capacity=0)
  let (_, _, edge_loss_map) = phylo_reconcile_rec(
    gene_tree,
    species_tree,
    gene_to_species,
    sep,
    nodes
  )
  // Count summary statistics
  let mut dup_count = 0
  let mut spec_count = 0
  let mut total_losses = 0
  for node in nodes {
    if node.is_duplication {
      dup_count = dup_count + 1
    }
    if node.is_speciation {
      spec_count = spec_count + 1
    }
    total_losses = total_losses + node.edge_losses
  }
  let leaf_count = gene_tree.count_terminals()
  ReconciliationResult::{
    duplication_count: dup_count,
    loss_count: total_losses,
    speciation_count: spec_count,
    total_cost: dup_count + total_losses,
    node_count: nodes.length(),
    leaf_count,
    nodes,
  }
}

///|
/// Recursive helper for reconciliation.
/// Returns (species_key, species_names, edge_loss_map).
/// Populates `nodes` array with per-node information.
fn phylo_reconcile_rec(
  node : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String,
  out : Array[NodeReconciliation]
) -> (String, Array[String], Map[String, Int]) {
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
    out.push(NodeReconciliation::{
      gene_leaves: gene_key,
      species_leaves: species_key,
      is_duplication: false,
      is_speciation: false,
      edge_losses: 0,
      species_names: [species_name],
    })
    let m : Map[String, Int] = Map([], capacity=0)
    return (species_key, [species_name], m)
  }
  // Internal node: recurse into children
  let child_results : Array[(String, Array[String], Map[String, Int])] = Array(
    capacity=0
  )
  for ch in node.clades {
    child_results.push(
      phylo_reconcile_rec(ch, species_tree, gene_to_species, sep, out)
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
  recon_sort_strings(dedup)
  let parent_key = recon_find_lca_key(species_tree, dedup)
  // Check duplication
  let mut is_dup = false
  for cr in child_results {
    if cr.0 == parent_key && child_results.length() > 1 && parent_key != "" {
      is_dup = true
    }
  }
  let is_spec = not(is_dup) && child_results.length() > 1
  // Count edge losses
  let mut edge_losses = 0
  for cr in child_results {
    let child_key = cr.0
    if child_key != parent_key && child_key != "" && parent_key != "" {
      let p_depth = recon_depth(species_tree, parent_key)
      let c_depth = recon_depth(species_tree, child_key)
      let diff = c_depth - p_depth
      if diff > 0 {
        if is_dup {
          edge_losses = edge_losses + diff
        } else {
          edge_losses = edge_losses + diff - 1
        }
      }
    }
  }
  // Build gene leaf key
  let all_gene_names = recon_collect_leaf_names(node)
  let gene_key = recon_leaf_key(all_gene_names)
  out.push(NodeReconciliation::{
    gene_leaves: gene_key,
    species_leaves: parent_key,
    is_duplication: is_dup,
    is_speciation: is_spec,
    edge_losses,
    species_names: dedup,
  })
  // Build edge loss map (gene_key → edge_losses for this node)
  let loss_map : Map[String, Int] = Map([], capacity=0)
  loss_map.set(gene_key, edge_losses)
  (parent_key, dedup, loss_map)
}

///|
/// Format a reconciliation result as a human-readable summary.
pub fn phylo_reconciliation_summary(result : ReconciliationResult) -> String {
  let sb = StringBuilder()
  sb.write_string("Reconciliation Summary:\n")
  sb.write_string("  Gene tree nodes: ")
  sb.write_string(result.node_count.to_string())
  sb.write_string(" (")
  sb.write_string(result.leaf_count.to_string())
  sb.write_string(" leaves)\n")
  sb.write_string("  Duplications: ")
  sb.write_string(result.duplication_count.to_string())
  sb.write_string("\n")
  sb.write_string("  Speciations:  ")
  sb.write_string(result.speciation_count.to_string())
  sb.write_string("\n")
  sb.write_string("  Losses:        ")
  sb.write_string(result.loss_count.to_string())
  sb.write_string("\n")
  sb.write_string("  Total cost:    ")
  sb.write_string(result.total_cost.to_string())
  sb.write_string(" (D + L)\n")
  // Per-node details
  sb.write_string("  Per-node details:\n")
  for i in 0..<result.nodes.length() {
    let n = result.nodes[i]
    sb.write_string("    [")
    sb.write_string(i.to_string())
    sb.write_string("] ")
    if n.is_duplication {
      sb.write_string("DUP   ")
    } else if n.is_speciation {
      sb.write_string("SPEC  ")
    } else {
      sb.write_string("LEAF  ")
    }
    sb.write_string("losses=")
    sb.write_string(n.edge_losses.to_string())
    sb.write_string("  species={")
    for j in 0..<n.species_names.length() {
      if j > 0 {
        sb.write_char(',')
      }
      sb.write_string(n.species_names[j])
    }
    sb.write_string("}\n")
  }
  sb.to_string()
}

// =========================================================================
// Ortholog/Paralog classification (ete3)
//   For each pair of genes, find their LCA in the gene tree.
//   If the LCA is a speciation → orthologs (diverged via speciation).
//   If the LCA is a duplication → paralogs (diverged via duplication).
//   Super-paralog: duplication within the same species (in-paralog).
//   Ultra-paralog: super-paralog with no intervening speciation.
// =========================================================================

///|
/// Gene relationship type.
pub(all) enum GeneRelationship {
  Ortholog       // diverged via speciation
  Paralog        // diverged via duplication (out-paralog)
  SuperParalog   // in-paralog: duplication within same species
  UltraParalog   // super-paralog with no intervening speciation
}

///|
/// Classification result for a pair of genes.
pub(all) struct PairClassification {
  gene1 : String
  gene2 : String
  species1 : String
  species2 : String
  relationship : GeneRelationship
  lca_species : String   // species LCA of the gene-tree LCA node
}

///|
/// Classify a pair of genes as ortholog or paralog.
/// Returns None if either gene is not found.
pub fn phylo_classify_pair(
  gene_tree : Clade,
  species_tree : Clade,
  gene1 : String,
  gene2 : String,
  sep : String
) -> PairClassification? {
  let gene_to_species : Map[String, String] = Map([], capacity=0)
  phylo_classify_pair_explicit(
    gene_tree, species_tree, gene_to_species, gene1, gene2, sep
  )
}

///|
/// Classify a pair with explicit gene-to-species map.
pub fn phylo_classify_pair_explicit(
  gene_tree : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  gene1 : String,
  gene2 : String,
  sep : String
) -> PairClassification? {
  // Find LCA of gene1 and gene2 in the gene tree
  let lca = recon_find_lca_node(gene_tree, [gene1, gene2])
  match lca {
    None => None
    Some(lca_node) => {
      // Get species names for both genes
      let s1 = match gene_to_species.get(gene1) {
        Some(s) => s
        None => phylo_extract_species(gene1, sep)
      }
      let s2 = match gene_to_species.get(gene2) {
        Some(s) => s
        None => phylo_extract_species(gene2, sep)
      }
      // Determine if LCA is a duplication
      let is_dup = phylo_is_duplication_node(
        lca_node, species_tree, gene_to_species, sep
      )
      // Get LCA's species map
      let lca_gene_leaves = recon_collect_leaf_names(lca_node)
      let lca_species_names : Array[String] = Array(capacity=0)
      for g in lca_gene_leaves {
        let sn = match gene_to_species.get(g) {
          Some(s) => s
          None => phylo_extract_species(g, sep)
        }
        lca_species_names.push(sn)
      }
      let lca_species_key = recon_find_lca_key(species_tree, lca_species_names)
      // Classify
      let rel = if not(is_dup) {
        Ortholog
      } else if s1 == s2 {
        // Duplication within same species → super-paralog
        // Check for ultra-paralog: no speciation between dup and LCA of species
        if phylo_is_ultra_paralog(
          gene_tree, lca_node, species_tree, gene_to_species, sep
        ) {
          UltraParalog
        } else {
          SuperParalog
        }
      } else {
        Paralog
      }
      Some(PairClassification::{
        gene1,
        gene2,
        species1: s1,
        species2: s2,
        relationship: rel,
        lca_species: lca_species_key,
      })
    }
  }
}

///|
/// Check if an internal node in the gene tree is a duplication.
/// A node is a duplication if its LCA map equals any child's LCA map.
fn phylo_is_duplication_node(
  node : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> Bool {
  if node.is_terminal() || node.clades.length() < 2 {
    return false
  }
  // Get species names for this node
  let all_leaves = recon_collect_leaf_names(node)
  let species_names : Array[String] = Array(capacity=0)
  for g in all_leaves {
    let sn = match gene_to_species.get(g) {
      Some(s) => s
      None => phylo_extract_species(g, sep)
    }
    species_names.push(sn)
  }
  let parent_key = recon_find_lca_key(species_tree, species_names)
  // Check each child
  for ch in node.clades {
    let ch_leaves = recon_collect_leaf_names(ch)
    let ch_species : Array[String] = Array(capacity=0)
    for g in ch_leaves {
      let sn = match gene_to_species.get(g) {
        Some(s) => s
        None => phylo_extract_species(g, sep)
      }
      ch_species.push(sn)
    }
    let ch_key = recon_find_lca_key(species_tree, ch_species)
    if ch_key == parent_key {
      return true
    }
  }
  false
}

///|
/// Check if a super-paralog is an ultra-paralog.
/// Ultra-paralog: all ancestors up to the root are also duplications
/// (no intervening speciation).
fn phylo_is_ultra_paralog(
  gene_tree : Clade,
  dup_node : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> Bool {
  // Find path from root to dup_node, check all ancestors are duplications
  let dup_key = recon_leaf_key(recon_collect_leaf_names(dup_node))
  phylo_check_ancestors(
    gene_tree, dup_key, species_tree, gene_to_species, sep, true
  )
}

///|
/// Recursively check if all ancestors of the node with `target_key` are
/// duplications. The first call starts with `is_ancestor=true` meaning
/// we haven't found the target yet.
fn phylo_check_ancestors(
  node : Clade,
  target_key : String,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String,
  found_target : Bool
) -> Bool {
  let my_key = recon_leaf_key(recon_collect_leaf_names(node))
  if my_key == target_key {
    // Found the target node — all ancestors were duplications
    return true
  }
  if found_target {
    // We're above the target — check if this node is a duplication
    if node.clades.length() >= 2 {
      if not(phylo_is_duplication_node(
        node, species_tree, gene_to_species, sep
      )) {
        // Found a speciation ancestor → not ultra-paralog
        return false
      }
    }
  }
  // Recurse into children
  for ch in node.clades {
    let ch_key = recon_leaf_key(recon_collect_leaf_names(ch))
    // Check if target is in this subtree
    if phylo_key_contains(ch_key, target_key) {
      return phylo_check_ancestors(
        ch, target_key, species_tree, gene_to_species, sep, true
      )
    }
  }
  false
}

///|
/// Check if `outer` key contains `inner` key (all inner species in outer).
fn phylo_key_contains(outer : String, inner : String) -> Bool {
  let outer_names = recon_key_to_names(outer)
  let inner_names = recon_key_to_names(inner)
  let outer_set : Map[String, Bool] = Map([], capacity=0)
  for n in outer_names {
    outer_set.set(n, true)
  }
  for n in inner_names {
    if not(outer_set.contains(n)) {
      return false
    }
  }
  true
}

///|
/// Classify all pairs of gene leaves in the gene tree.
pub fn phylo_classify_all_pairs(
  gene_tree : Clade,
  species_tree : Clade,
  sep : String
) -> Array[PairClassification] {
  let gene_to_species : Map[String, String] = Map([], capacity=0)
  phylo_classify_all_pairs_explicit(gene_tree, species_tree, gene_to_species, sep)
}

///|
/// Classify all pairs with explicit gene-to-species map.
pub fn phylo_classify_all_pairs_explicit(
  gene_tree : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> Array[PairClassification] {
  let leaves = recon_collect_leaf_names(gene_tree)
  let results : Array[PairClassification] = Array(capacity=0)
  for i in 0..<leaves.length() {
    for j in (i + 1)..<leaves.length() {
      match phylo_classify_pair_explicit(
        gene_tree, species_tree, gene_to_species, leaves[i], leaves[j], sep
      ) {
        Some(cls) => results.push(cls)
        None => ()
      }
    }
  }
  results
}

///|
/// Extract only ortholog pairs from classifications.
pub fn phylo_orthologs(
  classifications : Array[PairClassification]
) -> Array[PairClassification] {
  let result : Array[PairClassification] = Array(capacity=0)
  for c in classifications {
    match c.relationship {
      Ortholog => result.push(c)
      _ => ()
    }
  }
  result
}

///|
/// Extract only paralog pairs (including super/ultra) from classifications.
pub fn phylo_paralogs(
  classifications : Array[PairClassification]
) -> Array[PairClassification] {
  let result : Array[PairClassification] = Array(capacity=0)
  for c in classifications {
    match c.relationship {
      Paralog => result.push(c)
      SuperParalog => result.push(c)
      UltraParalog => result.push(c)
      _ => ()
    }
  }
  result
}

///|
/// Format a pair classification as a string.
pub fn phylo_pair_to_string(cls : PairClassification) -> String {
  let rel_str = match cls.relationship {
    Ortholog => "Ortholog"
    Paralog => "Paralog"
    SuperParalog => "SuperParalog"
    UltraParalog => "UltraParalog"
  }
  let sb = StringBuilder()
  sb.write_string(cls.gene1)
  sb.write_string(" -- ")
  sb.write_string(cls.gene2)
  sb.write_string(": ")
  sb.write_string(rel_str)
  sb.write_string(" (")
  sb.write_string(cls.species1)
  sb.write_string(", ")
  sb.write_string(cls.species2)
  sb.write_string(")")
  sb.to_string()
}

///|
/// Format all pair classifications as a summary.
pub fn phylo_pairs_summary(
  classifications : Array[PairClassification]
) -> String {
  let mut ortho_count = 0
  let mut para_count = 0
  let mut super_count = 0
  let mut ultra_count = 0
  for c in classifications {
    match c.relationship {
      Ortholog => ortho_count = ortho_count + 1
      Paralog => para_count = para_count + 1
      SuperParalog => super_count = super_count + 1
      UltraParalog => ultra_count = ultra_count + 1
    }
  }
  let sb = StringBuilder()
  sb.write_string("Ortholog/Paralog Classification:\n")
  sb.write_string("  Total pairs: ")
  sb.write_string(classifications.length().to_string())
  sb.write_string("\n")
  sb.write_string("  Orthologs:     ")
  sb.write_string(ortho_count.to_string())
  sb.write_string("\n")
  sb.write_string("  Paralogs:      ")
  sb.write_string(para_count.to_string())
  sb.write_string("\n")
  sb.write_string("  SuperParalogs: ")
  sb.write_string(super_count.to_string())
  sb.write_string("\n")
  sb.write_string("  UltraParalogs: ")
  sb.write_string(ultra_count.to_string())
  sb.write_string("\n")
  for c in classifications {
    sb.write_string("  ")
    sb.write_string(phylo_pair_to_string(c))
    sb.write_string("\n")
  }
  sb.to_string()
}

// =========================================================================
// Reconciliation Context Object (P0 optimization)
//   Inspired by scvi-tools model pattern: compute once, query many times.
//   Pre-computes:
//     1. LCA mapping for all gene tree nodes (one traversal)
//     2. Species tree depth table (Map[leaf_key → depth], one traversal)
//   Then exposes count_duplications / count_losses / reconcile /
//   classify_all_pairs / summary as methods that reuse cached data.
//   This eliminates the repeated LCA recomputation that occurred when
//   calling phylo_lca_mapping + phylo_count_losses + phylo_reconcile
//   separately (each redid the full mapping from scratch).
// =========================================================================

///|
/// Reconciliation context: pre-computes LCA mapping and species depth
/// table once, then supports multiple queries without recomputation.
///
/// Usage:
///   let recon = Reconciliation::new(gene_tree, species_tree, "_")
///   let dups = recon.count_duplications()
///   let losses = recon.count_losses()
///   let pairs = recon.classify_all_pairs()
///   let summary = recon.summary()
pub struct Reconciliation {
  gene_tree : Clade
  species_tree : Clade
  sep : String
  lca_maps : Array[LcaMapResult]           // pre-computed LCA mapping
  depth_table : Map[String, Int]           // species leaf_key → depth
  gene_to_species : Map[String, String]
}

///|
/// Construct a Reconciliation context by pre-computing the LCA mapping
/// and species depth table in a single pass each.
pub fn Reconciliation::new(
  gene_tree : Clade,
  species_tree : Clade,
  sep : String
) -> Reconciliation {
  let gene_to_species : Map[String, String] = Map([], capacity=0)
  Reconciliation::new_explicit(gene_tree, species_tree, gene_to_species, sep)
}

///|
/// Construct with an explicit gene-to-species map.
pub fn Reconciliation::new_explicit(
  gene_tree : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String
) -> Reconciliation {
  // Pre-compute LCA mapping (one traversal)
  let lca_maps = phylo_lca_mapping_explicit(
    gene_tree, species_tree, gene_to_species, sep
  )
  // Pre-compute species depth table (one traversal)
  let depth_table : Map[String, Int] = Map([], capacity=0)
  recon_build_depth_table(species_tree, 0, depth_table)
  Reconciliation::{
    gene_tree,
    species_tree,
    sep,
    lca_maps,
    depth_table,
    gene_to_species,
  }
}

///|
/// Build a depth table: maps each node's leaf_key to its depth.
/// One traversal of the species tree, replaces O(n) per-depth-query
/// lookups with O(1) map lookups.
fn recon_build_depth_table(
  clade : Clade,
  depth : Int,
  table : Map[String, Int]
) -> Unit {
  let leaves = recon_collect_leaf_names(clade)
  let key = recon_leaf_key(leaves)
  table.set(key, depth)
  for ch in clade.clades {
    recon_build_depth_table(ch, depth + 1, table)
  }
}

///|
/// Return the pre-computed LCA mapping.
pub fn Reconciliation::lca_maps(self : Reconciliation) -> Array[LcaMapResult] {
  self.lca_maps
}

///|
/// Count duplications from the pre-computed LCA mapping.
pub fn Reconciliation::count_duplications(self : Reconciliation) -> Int {
  phylo_count_duplications(self.lca_maps)
}

///|
/// Count losses using the pre-computed depth table (O(1) per depth query
/// instead of O(n) tree traversal).
pub fn Reconciliation::count_losses(self : Reconciliation) -> Int {
  let (_, _, losses) = recon_count_losses_cached(
    self.gene_tree,
    self.species_tree,
    self.gene_to_species,
    self.sep,
    self.depth_table
  )
  losses
}

///|
/// Recursive helper for loss counting with cached depth lookups.
fn recon_count_losses_cached(
  node : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String,
  depth_table : Map[String, Int]
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
  let child_results : Array[(String, Array[String], Int)] = Array(capacity=0)
  for ch in node.clades {
    child_results.push(
      recon_count_losses_cached(
        ch, species_tree, gene_to_species, sep, depth_table
      )
    )
  }
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
  recon_sort_strings(dedup)
  let parent_key = recon_find_lca_key(species_tree, dedup)
  let mut is_dup = false
  for cr in child_results {
    if cr.0 == parent_key && child_results.length() > 1 && parent_key != "" {
      is_dup = true
    }
  }
  let mut total_losses = 0
  for cr in child_results {
    let child_key = cr.0
    if child_key != parent_key && child_key != "" && parent_key != "" {
      // O(1) depth lookup from pre-computed table
      let p_depth = match depth_table.get(parent_key) {
        Some(d) => d
        None => -1
      }
      let c_depth = match depth_table.get(child_key) {
        Some(d) => d
        None => -1
      }
      if p_depth >= 0 && c_depth >= 0 {
        let diff = c_depth - p_depth
        if diff > 0 {
          if is_dup {
            total_losses = total_losses + diff
          } else {
            total_losses = total_losses + diff - 1
          }
        }
      }
    }
    total_losses = total_losses + cr.2
  }
  (parent_key, dedup, total_losses)
}

///|
/// Full reconciliation with per-node detail, using cached depth table.
pub fn Reconciliation::reconcile(self : Reconciliation) -> ReconciliationResult {
  let nodes : Array[NodeReconciliation] = Array(capacity=0)
  ignore(recon_reconcile_cached(
    self.gene_tree,
    self.species_tree,
    self.gene_to_species,
    self.sep,
    self.depth_table,
    nodes
  ))
  let mut dup_count = 0
  let mut spec_count = 0
  let mut total_losses = 0
  for node in nodes {
    if node.is_duplication {
      dup_count = dup_count + 1
    }
    if node.is_speciation {
      spec_count = spec_count + 1
    }
    total_losses = total_losses + node.edge_losses
  }
  let leaf_count = self.gene_tree.count_terminals()
  ReconciliationResult::{
    duplication_count: dup_count,
    loss_count: total_losses,
    speciation_count: spec_count,
    total_cost: dup_count + total_losses,
    node_count: nodes.length(),
    leaf_count,
    nodes,
  }
}

///|
/// Recursive helper for reconciliation with cached depth lookups.
fn recon_reconcile_cached(
  node : Clade,
  species_tree : Clade,
  gene_to_species : Map[String, String],
  sep : String,
  depth_table : Map[String, Int],
  out : Array[NodeReconciliation]
) -> (String, Array[String], Map[String, Int]) {
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
    out.push(NodeReconciliation::{
      gene_leaves: gene_key,
      species_leaves: species_key,
      is_duplication: false,
      is_speciation: false,
      edge_losses: 0,
      species_names: [species_name],
    })
    let m : Map[String, Int] = Map([], capacity=0)
    return (species_key, [species_name], m)
  }
  let child_results : Array[(String, Array[String], Map[String, Int])] = Array(
    capacity=0
  )
  for ch in node.clades {
    child_results.push(
      recon_reconcile_cached(
        ch, species_tree, gene_to_species, sep, depth_table, out
      )
    )
  }
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
  recon_sort_strings(dedup)
  let parent_key = recon_find_lca_key(species_tree, dedup)
  let mut is_dup = false
  for cr in child_results {
    if cr.0 == parent_key && child_results.length() > 1 && parent_key != "" {
      is_dup = true
    }
  }
  let is_spec = not(is_dup) && child_results.length() > 1
  let mut edge_loss = 0
  for cr in child_results {
    let child_key = cr.0
    if child_key != parent_key && child_key != "" && parent_key != "" {
      let p_depth = match depth_table.get(parent_key) {
        Some(d) => d
        None => -1
      }
      let c_depth = match depth_table.get(child_key) {
        Some(d) => d
        None => -1
      }
      if p_depth >= 0 && c_depth >= 0 {
        let diff = c_depth - p_depth
        if diff > 0 {
          if is_dup {
            edge_loss = edge_loss + diff
          } else {
            edge_loss = edge_loss + diff - 1
          }
        }
      }
    }
  }
  let all_gene_names = recon_collect_leaf_names(node)
  let gene_key = recon_leaf_key(all_gene_names)
  out.push(NodeReconciliation::{
    gene_leaves: gene_key,
    species_leaves: parent_key,
    is_duplication: is_dup,
    is_speciation: is_spec,
    edge_losses: edge_loss,
    species_names: dedup,
  })
  let loss_map : Map[String, Int] = Map([], capacity=0)
  loss_map.set(gene_key, edge_loss)
  (parent_key, dedup, loss_map)
}

///|
/// Classify all gene pairs as ortholog/paralog.
/// Reuses the pre-computed LCA mapping for duplication detection.
pub fn Reconciliation::classify_all_pairs(
  self : Reconciliation
) -> Array[PairClassification] {
  phylo_classify_all_pairs_explicit(
    self.gene_tree,
    self.species_tree,
    self.gene_to_species,
    self.sep
  )
}

///|
/// Classify a single gene pair.
pub fn Reconciliation::classify_pair(
  self : Reconciliation,
  gene1 : String,
  gene2 : String
) -> PairClassification? {
  phylo_classify_pair_explicit(
    self.gene_tree,
    self.species_tree,
    self.gene_to_species,
    gene1,
    gene2,
    self.sep
  )
}

///|
/// Formatted summary: duplication/loss/speciation counts + total cost.
pub fn Reconciliation::summary(self : Reconciliation) -> String {
  let result = self.reconcile()
  phylo_reconciliation_summary(result)
}
