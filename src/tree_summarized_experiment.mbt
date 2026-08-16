///|
/// TreeSummarizedExperiment - rectangular assays linked to row and column trees.
/// Inspired by the Bioconductor TreeSummarizedExperiment package.

///|
/// Link between one assay dimension entry and a node in a tree.
pub struct TseLink {
  node_label : String
  node_alias : String
  node_number : Int
  is_leaf : Bool
  tree_name : String
}

///|
/// Construct a tree link.
pub fn TseLink::new(
  node_label~ : String,
  node_alias~ : String,
  node_number~ : Int,
  is_leaf~ : Bool,
  tree_name? : String = "phylo",
) -> TseLink {
  TseLink::{ node_label, node_alias, node_number, is_leaf, tree_name }
}

///|
/// Return the node label.
pub fn TseLink::node_label(self : TseLink) -> String {
  self.node_label
}

///|
/// Return the stable node alias.
pub fn TseLink::node_alias(self : TseLink) -> String {
  self.node_alias
}

///|
/// Return the one-based node number.
pub fn TseLink::node_number(self : TseLink) -> Int {
  self.node_number
}

///|
/// Return whether the linked node is a leaf.
pub fn TseLink::is_leaf(self : TseLink) -> Bool {
  self.is_leaf
}

///|
/// Return the tree name used by the link.
pub fn TseLink::tree_name(self : TseLink) -> String {
  self.tree_name
}

///|
/// Supported aggregation functions for tree-level summaries.
pub enum TseAggregation {
  Sum
  Mean
  Min
  Max
} derive(Eq)

///|
/// Create a sum aggregation strategy.
pub fn tse_aggregation_sum() -> TseAggregation {
  TseAggregation::Sum
}

///|
/// Create a mean aggregation strategy.
pub fn tse_aggregation_mean() -> TseAggregation {
  TseAggregation::Mean
}

///|
/// Create a minimum aggregation strategy.
pub fn tse_aggregation_min() -> TseAggregation {
  TseAggregation::Min
}

///|
/// Create a maximum aggregation strategy.
pub fn tse_aggregation_max() -> TseAggregation {
  TseAggregation::Max
}

///|
/// A SummarizedExperiment with optional hierarchical structures on both axes.
pub struct TreeSummarizedExperiment {
  experiment : SummarizedExperiment
  row_tree : Tree?
  col_tree : Tree?
  row_links : Array[TseLink]
  col_links : Array[TseLink]
  reference_sequences : Array[String]
}

///|
/// Construct a TreeSummarizedExperiment.
///
/// `row_node_labels` and `col_node_labels` map assay entries to named tree
/// nodes. Repeated labels are allowed, matching the Bioconductor link model.
pub fn TreeSummarizedExperiment::new(
  experiment~ : SummarizedExperiment,
  row_tree? : Tree? = None,
  row_node_labels? : Array[String] = [],
  col_tree? : Tree? = None,
  col_node_labels? : Array[String] = [],
  reference_sequences? : Array[String] = [],
) -> TreeSummarizedExperiment {
  TreeSummarizedExperiment::{
    experiment,
    row_tree,
    col_tree,
    row_links: tse_make_links(row_tree, row_node_labels),
    col_links: tse_make_links(col_tree, col_node_labels),
    reference_sequences,
  }
}

///|
/// Return the underlying SummarizedExperiment.
pub fn TreeSummarizedExperiment::experiment(
  self : TreeSummarizedExperiment,
) -> SummarizedExperiment {
  self.experiment
}

///|
/// Return an assay by name.
pub fn TreeSummarizedExperiment::assay(
  self : TreeSummarizedExperiment,
  name : String,
) -> Array[Array[Double]]? {
  se_assay(self.experiment, name)
}

///|
/// Return all assay names.
pub fn TreeSummarizedExperiment::assay_names(
  self : TreeSummarizedExperiment,
) -> Array[String] {
  se_assay_names(self.experiment)
}

///|
/// Return the number of assay rows.
pub fn TreeSummarizedExperiment::nrow(self : TreeSummarizedExperiment) -> Int {
  se_nrow(self.experiment)
}

///|
/// Return the number of assay columns.
pub fn TreeSummarizedExperiment::ncol(self : TreeSummarizedExperiment) -> Int {
  se_ncol(self.experiment)
}

///|
/// Return the row tree.
pub fn TreeSummarizedExperiment::row_tree(
  self : TreeSummarizedExperiment,
) -> Tree? {
  self.row_tree
}

///|
/// Return the column tree.
pub fn TreeSummarizedExperiment::col_tree(
  self : TreeSummarizedExperiment,
) -> Tree? {
  self.col_tree
}

///|
/// Return links between assay rows and row-tree nodes.
pub fn TreeSummarizedExperiment::row_links(
  self : TreeSummarizedExperiment,
) -> Array[TseLink] {
  self.row_links
}

///|
/// Return links between assay columns and column-tree nodes.
pub fn TreeSummarizedExperiment::col_links(
  self : TreeSummarizedExperiment,
) -> Array[TseLink] {
  self.col_links
}

///|
/// Return optional per-row reference sequences.
pub fn TreeSummarizedExperiment::reference_sequences(
  self : TreeSummarizedExperiment,
) -> Array[String] {
  self.reference_sequences
}

///|
/// Check dimensions, tree links, reference sequences, and assay shapes.
pub fn TreeSummarizedExperiment::is_valid(
  self : TreeSummarizedExperiment,
) -> Bool {
  let nrow = self.nrow()
  let ncol = self.ncol()
  match self.row_tree {
    Some(_) => if self.row_links.length() != nrow { return false }
    None => if self.row_links.length() != 0 { return false }
  }
  match self.col_tree {
    Some(_) => if self.col_links.length() != ncol { return false }
    None => if self.col_links.length() != 0 { return false }
  }
  if self.reference_sequences.length() != 0 &&
    self.reference_sequences.length() != nrow {
    return false
  }
  for assay_name in self.assay_names() {
    match self.assay(assay_name) {
      Some(assay) => {
        if assay.length() != nrow {
          return false
        }
        for row in assay {
          if row.length() != ncol {
            return false
          }
        }
      }
      None => return false
    }
  }
  true
}

///|
/// Subset assay rows while preserving the full tree and updating row links.
pub fn TreeSummarizedExperiment::subset_rows(
  self : TreeSummarizedExperiment,
  indices : Array[Int],
) -> TreeSummarizedExperiment {
  let links : Array[TseLink] = []
  let sequences : Array[String] = []
  for index in indices {
    if index >= 0 && index < self.row_links.length() {
      links.push(self.row_links[index])
    }
    if index >= 0 && index < self.reference_sequences.length() {
      sequences.push(self.reference_sequences[index])
    }
  }
  TreeSummarizedExperiment::{
    experiment: se_subset_rows(self.experiment, indices),
    row_tree: self.row_tree,
    col_tree: self.col_tree,
    row_links: links,
    col_links: self.col_links,
    reference_sequences: sequences,
  }
}

///|
/// Subset assay columns while preserving the full tree and updating col links.
pub fn TreeSummarizedExperiment::subset_cols(
  self : TreeSummarizedExperiment,
  indices : Array[Int],
) -> TreeSummarizedExperiment {
  let links : Array[TseLink] = []
  for index in indices {
    if index >= 0 && index < self.col_links.length() {
      links.push(self.col_links[index])
    }
  }
  TreeSummarizedExperiment::{
    experiment: se_subset_cols(self.experiment, indices),
    row_tree: self.row_tree,
    col_tree: self.col_tree,
    row_links: self.row_links,
    col_links: links,
    reference_sequences: self.reference_sequences,
  }
}

///|
/// Subset rows linked to the selected nodes or any of their descendants.
pub fn TreeSummarizedExperiment::subset_by_row_nodes(
  self : TreeSummarizedExperiment,
  node_labels : Array[String],
) -> TreeSummarizedExperiment? {
  match self.row_tree {
    None => None
    Some(tree) => {
      let indices = tse_link_indices_under(tree, self.row_links, node_labels)
      Some(self.subset_rows(indices))
    }
  }
}

///|
/// Subset columns linked to the selected nodes or any of their descendants.
pub fn TreeSummarizedExperiment::subset_by_col_nodes(
  self : TreeSummarizedExperiment,
  node_labels : Array[String],
) -> TreeSummarizedExperiment? {
  match self.col_tree {
    None => None
    Some(tree) => {
      let indices = tse_link_indices_under(tree, self.col_links, node_labels)
      Some(self.subset_cols(indices))
    }
  }
}

///|
/// Aggregate assay rows to arbitrary nodes in the row tree.
pub fn TreeSummarizedExperiment::aggregate_rows(
  self : TreeSummarizedExperiment,
  node_labels : Array[String],
  aggregation? : TseAggregation = TseAggregation::Sum,
) -> TreeSummarizedExperiment? {
  match self.row_tree {
    None => None
    Some(tree) => {
      let targets = tse_existing_nodes(tree, node_labels)
      if targets.length() == 0 {
        return None
      }
      let groups = tse_link_groups_under(tree, self.row_links, targets)
      let assays : Map[String, Array[Array[Double]]] = Map([])
      for assay_name in self.assay_names() {
        match self.assay(assay_name) {
          None => ()
          Some(assay) => {
            let result : Array[Array[Double]] = []
            for group in groups {
              let row : Array[Double] = []
              for col = 0; col < self.ncol(); col = col + 1 {
                let values : Array[Double] = []
                for source_row in group {
                  values.push(assay[source_row][col])
                }
                row.push(tse_aggregate_values(values, aggregation))
              }
              result.push(row)
            }
            assays[assay_name] = result
          }
        }
      }
      let row_ranges = targets.map(fn(label) { (label, 0, 0) })
      let experiment = summarized_experiment(
        assays,
        row_ranges,
        se_col_data(self.experiment),
        se_metadata(self.experiment),
      )
      Some(TreeSummarizedExperiment::{
        experiment,
        row_tree: self.row_tree,
        col_tree: self.col_tree,
        row_links: tse_make_links(self.row_tree, targets),
        col_links: self.col_links,
        reference_sequences: [],
      })
    }
  }
}

///|
/// Aggregate assay columns to arbitrary nodes in the column tree.
pub fn TreeSummarizedExperiment::aggregate_cols(
  self : TreeSummarizedExperiment,
  node_labels : Array[String],
  aggregation? : TseAggregation = TseAggregation::Sum,
) -> TreeSummarizedExperiment? {
  match self.col_tree {
    None => None
    Some(tree) => {
      let targets = tse_existing_nodes(tree, node_labels)
      if targets.length() == 0 {
        return None
      }
      let groups = tse_link_groups_under(tree, self.col_links, targets)
      let assays : Map[String, Array[Array[Double]]] = Map([])
      for assay_name in self.assay_names() {
        match self.assay(assay_name) {
          None => ()
          Some(assay) => {
            let result : Array[Array[Double]] = []
            for source_row = 0
                source_row < self.nrow()
                source_row = source_row + 1 {
              let row : Array[Double] = []
              for group in groups {
                let values : Array[Double] = []
                for source_col in group {
                  values.push(assay[source_row][source_col])
                }
                row.push(tse_aggregate_values(values, aggregation))
              }
              result.push(row)
            }
            assays[assay_name] = result
          }
        }
      }
      let col_data = targets.map(fn(label) { Map([("tree_node", label)]) })
      let experiment = summarized_experiment(
        assays,
        se_row_ranges(self.experiment),
        col_data,
        se_metadata(self.experiment),
      )
      Some(TreeSummarizedExperiment::{
        experiment,
        row_tree: self.row_tree,
        col_tree: self.col_tree,
        row_links: self.row_links,
        col_links: tse_make_links(self.col_tree, targets),
        reference_sequences: self.reference_sequences,
      })
    }
  }
}

///|
/// Return descendant labels of a node.
pub fn tse_find_descendants(
  tree : Tree,
  node_label : String,
  leaves_only? : Bool = false,
) -> Array[String] {
  match tree.find_any(name=Some(node_label)) {
    None => []
    Some(node) => {
      let result : Array[String] = []
      tse_collect_descendant_labels(node, leaves_only, result)
      result
    }
  }
}

///|
/// Return ancestor labels from root to the parent of `node_label`.
pub fn tse_find_ancestors(tree : Tree, node_label : String) -> Array[String] {
  match tree.get_path(node_label) {
    None => []
    Some(path) => {
      let result : Array[String] = []
      for index = 0; index + 1 < path.length(); index = index + 1 {
        match path[index].name {
          Some(name) => result.push(name)
          None => ()
        }
      }
      result
    }
  }
}

///|
/// Return whether a named tree node is a leaf.
pub fn tse_is_leaf(tree : Tree, node_label : String) -> Bool {
  match tree.find_any(name=Some(node_label)) {
    Some(node) => node.is_terminal()
    None => false
  }
}

///|
fn tse_tree_name(tree : Tree) -> String {
  tree.name.unwrap_or("phylo")
}

///|
fn tse_ordered_nodes(tree : Tree) -> Array[Clade] {
  let nodes = tree.get_terminals()
  for node in tree.get_nonterminals() {
    nodes.push(node)
  }
  nodes
}

///|
fn tse_make_links(tree : Tree?, labels : Array[String]) -> Array[TseLink] {
  match tree {
    None => []
    Some(value) => {
      let nodes = tse_ordered_nodes(value)
      let result : Array[TseLink] = []
      for label in labels {
        let mut found = false
        for index = 0; index < nodes.length(); index = index + 1 {
          if !found && nodes[index].name == Some(label) {
            result.push(
              TseLink::new(
                node_label=label,
                node_alias="alias_" + (index + 1).to_string(),
                node_number=index + 1,
                is_leaf=nodes[index].is_terminal(),
                tree_name=tse_tree_name(value),
              ),
            )
            found = true
          }
        }
      }
      result
    }
  }
}

///|
fn tse_node_is_under(
  tree : Tree,
  node_label : String,
  ancestor : String,
) -> Bool {
  if node_label == ancestor {
    return true
  }
  match tree.get_path(node_label) {
    None => false
    Some(path) => {
      for node in path {
        if node.name == Some(ancestor) {
          return true
        }
      }
      false
    }
  }
}

///|
fn tse_link_indices_under(
  tree : Tree,
  links : Array[TseLink],
  node_labels : Array[String],
) -> Array[Int] {
  let result : Array[Int] = []
  for index = 0; index < links.length(); index = index + 1 {
    let mut selected = false
    for node_label in node_labels {
      if tse_node_is_under(tree, links[index].node_label, node_label) {
        selected = true
      }
    }
    if selected {
      result.push(index)
    }
  }
  result
}

///|
fn tse_existing_nodes(tree : Tree, labels : Array[String]) -> Array[String] {
  let result : Array[String] = []
  for label in labels {
    if tree.find_any(name=Some(label)) is Some(_) {
      let mut duplicate = false
      for existing in result {
        if existing == label {
          duplicate = true
        }
      }
      if !duplicate {
        result.push(label)
      }
    }
  }
  result
}

///|
fn tse_link_groups_under(
  tree : Tree,
  links : Array[TseLink],
  node_labels : Array[String],
) -> Array[Array[Int]] {
  node_labels.map(fn(label) { tse_link_indices_under(tree, links, [label]) })
}

///|
fn tse_aggregate_values(
  values : Array[Double],
  aggregation : TseAggregation,
) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  match aggregation {
    TseAggregation::Sum => {
      let mut total = 0.0
      for value in values {
        total = total + value
      }
      total
    }
    TseAggregation::Mean => {
      let mut total = 0.0
      for value in values {
        total = total + value
      }
      total / values.length().to_double()
    }
    TseAggregation::Min => {
      let mut result = values[0]
      for value in values {
        if value < result {
          result = value
        }
      }
      result
    }
    TseAggregation::Max => {
      let mut result = values[0]
      for value in values {
        if value > result {
          result = value
        }
      }
      result
    }
  }
}

///|
fn tse_collect_descendant_labels(
  node : Clade,
  leaves_only : Bool,
  result : Array[String],
) -> Unit {
  for child in node.clades {
    if !leaves_only || child.is_terminal() {
      match child.name {
        Some(name) => result.push(name)
        None => ()
      }
    }
    tse_collect_descendant_labels(child, leaves_only, result)
  }
}
