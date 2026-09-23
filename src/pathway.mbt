///|
/// Bio.Pathway - Biochemical pathway analysis.
/// Ported from Biopython 1.88 Bio.Pathway module.
///
/// Core objects:
/// - **Species**  — a biochemical species (substrate, product, enzyme, …).
/// - **Reaction** — a (potentially reversible) biochemical transformation.
///   Species with negative stoichiometric coefficient are consumed
///   (substrates); species with positive coefficient are produced.
/// - **System**   — a collection of reactions; provides the
///   ``stochiometry`` matrix (species × reactions).
/// - **Network**  — a directed multigraph whose nodes are Species and
///   whose edges are Interaction references.
///
/// Graph algorithms (``Graph`` / ``MultiGraph``, BFS / DFS) live in
/// ``pathway_graph.mbt`` alongside this module.

// ===========================================================================
// Species
// ===========================================================================

///|
/// A biochemical species.  The same species participates in many
/// reactions as substrate, product, or catalyst.
pub struct Species {
  id : String
  name : String
}

///|
/// Create a Species with the given id and name (falls back to id when
/// no name is supplied).
pub fn Species::new(id : String, name? : String) -> Species {
  let nm = match name {
    Some(n) => n
    None => id
  }
  Species::{ id, name: nm }
}

///|
/// Species equality by id (Biopython convention).
pub fn Species::eq(self : Species, other : Species) -> Bool {
  self.id == other.id
}

// ===========================================================================
// Reaction
// ===========================================================================

///|
/// A biochemical reaction.
///
/// The ``reactants`` map uses signed stoichiometric coefficients so that
/// the same map represents both substrates and products:
///
/// ```
/// // 2 glucose + 2 ATP -> 2 G6P + 2 ADP
/// Reaction::new({"glucose": -2.0, "atp": -2.0, "g6p": 2.0, "adp": 2.0}, catalysts, reversible)
/// ```
///
/// - **negative** coefficient ⇒ substrate / consumed
/// - **positive** coefficient ⇒ product / produced
pub struct Reaction {
  id : String
  reactants : Map[String, Double]
  catalysts : Array[String]
  reversible : Bool
}

///|
/// Create a new Reaction.
///
/// - ``reactants``: signed stoichiometric map (neg=substrate, pos=product).
///   Entries with coefficient == 0 are silently dropped.
/// - ``catalysts`` : enzyme names (empty by default).
/// - ``reversible``: whether the reaction can run in both directions.
pub fn Reaction::new(
  id : String,
  reactants : Map[String, Double],
  catalysts? : Array[String],
  reversible? : Bool,
) -> Reaction {
  let clean : Map[String, Double] = Map([])
  let keys = reactants.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let k = keys[i]
    let v = reactants[k]
    if v != 0.0 {
      clean[k] = v
    }
    i = i + 1
  }
  let cat = match catalysts {
    Some(c) => c
    None => []
  }
  let rev = match reversible {
    Some(r) => r
    None => false
  }
  Reaction::{ id, reactants: clean, catalysts: cat, reversible: rev }
}

///|
/// Return a reversed version of this reaction (all stoichiometric
/// signs flipped; catalysts preserved).
pub fn Reaction::reverse(self : Reaction) -> Reaction {
  let flipped : Map[String, Double] = Map([])
  let keys = self.reactants.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let k = keys[i]
    flipped[k] = -self.reactants[k]
    i = i + 1
  }
  Reaction::{
    id: self.id + "::rev",
    reactants: flipped,
    catalysts: self.catalysts,
    reversible: self.reversible,
  }
}

///|
/// Species involved in this reaction (substrates + products).
pub fn Reaction::species(self : Reaction) -> Array[String] {
  self.reactants.keys().collect()
}

///|
/// Human-readable string: ``S1 + S2 --> P1 + P2`` or with ``<->`` for
/// reversible reactions.
pub fn Reaction::to_string(self : Reaction) -> String {
  let substrates : Array[String] = Array(capacity=0)
  let products : Array[String] = Array(capacity=0)
  let keys = self.reactants.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let k = keys[i]
    let c = self.reactants[k]
    if c < 0.0 {
      substrates.push((-c).to_string() + " " + k)
    } else if c > 0.0 {
      products.push(c.to_string() + " " + k)
    }
    i = i + 1
  }
  let arrow = if self.reversible { " <-> " } else { " --> " }
  let lhs = if substrates.length() == 0 { "..." } else { substrates.join(" + ") }
  let rhs = if products.length() == 0 { "..." } else { products.join(" + ") }
  lhs + arrow + rhs
}

// ===========================================================================
// System  (a collection of reactions + stoichiometry matrix)
// ===========================================================================

///|
/// A collection of reactions without explicit species/edge connectivity.
pub struct System {
  id : String
  reactions : Array[Reaction]
}

///|
/// Create a System from a slice of reactions.
pub fn System::new(id : String, reactions : Array[Reaction]) -> System {
  System::{ id, reactions }
}

///|
/// All unique species referenced by any reaction in this system.
pub fn System::species(self : System) -> Array[String] {
  let seen : Map[String, Bool] = Map([])
  let mut i = 0
  while i < self.reactions.length() {
    let sp = self.reactions[i].species()
    let mut j = 0
    while j < sp.length() {
      seen[sp[j]] = true
      j = j + 1
    }
    i = i + 1
  }
  seen.keys().collect()
}

///|
/// The reactions registered in this system.
pub fn System::reactions(self : System) -> Array[Reaction] {
  self.reactions
}

///|
/// Compute the stoichiometry matrix for this system.
///
/// Returns a triple ``(species, reactions, stoch)`` where
/// ``stoch[row_idx][col_idx]`` is the signed coefficient of
/// species ``species[row_idx]`` in reaction ``reactions[col_idx]``.
/// Negative ⇒ substrate, positive ⇒ product.
pub fn System::stochiometry(self : System) -> (Array[String], Array[Reaction], Array[Array[Double]]) {
  let species = self.species()
  let reactions = self.reactions
  let n_sp = species.length()
  let n_rx = reactions.length()

  // Initialize matrix with zeros.
  let stoch : Array[Array[Double]] = Array(capacity=0)
  let mut i = 0
  while i < n_sp {
    let row : Array[Double] = Array::make(n_rx, 0.0)
    stoch.push(row)
    i = i + 1
  }

  // Fill in coefficients.
  let mut j = 0
  while j < n_rx {
    let rxn = reactions[j]
    let keys = rxn.reactants.keys().collect()
    let mut k = 0
    while k < keys.length() {
      let sp = keys[k]
      let coef = rxn.reactants[sp]
      // Find row index for this species.
      let mut ri = 0
      while ri < n_sp {
        if species[ri] == sp {
          stoch[ri][j] = coef
          break
        }
        ri = ri + 1
      }
      k = k + 1
    }
    j = j + 1
  }

  (species, reactions, stoch)
}

// ===========================================================================
// Interaction (biochemical edge abstraction)
// ===========================================================================

///|
/// An arbitrary interaction between species.  In Bio.Pathway this is
/// just a wrapper; here we store ``source`` and ``sink`` species and
/// an optional label (typically a reference to the Reaction that
/// accounts for the edge).
pub struct Interaction {
  source : String
  sink : String
  label : String
  data : String?
}

///|
/// Build an Interaction.
pub fn Interaction::new(
  source : String,
  sink : String,
  label? : String,
  data? : String,
) -> Interaction {
  let lbl = match label {
    Some(l) => l
    None => ""
  }
  Interaction::{ source, sink, label: lbl, data }
}

// ===========================================================================
// Network  (directed multigraph with labelled edges)
// ===========================================================================

///|
/// A directed multigraph whose nodes are species and whose labelled
/// edges are Interactions.
///
/// This is the Bio.Pathway.Network class, re-implemented here so that
/// it does not depend on Bio.Pathway.Rep.{Graph,MultiGraph}.
pub struct Network {
  nodes : Map[String, Bool]
  /// Outgoing adjacency: node -> list of (neighbor, label)
  out_adj : Map[String, Array[(String, String)]]
  /// Incoming adjacency: node -> list of (neighbor, label)
  in_adj : Map[String, Array[(String, String)]]
}

///|
/// Create an empty Network.
pub fn Network::new() -> Network {
  Network::{ nodes: Map([]), out_adj: Map([]), in_adj: Map([]) }
}

///|
/// Register a species in this network (no-op if already present).
pub fn Network::add_species(self : Network, species : String) -> Unit {
  self.nodes[species] = true
  if !self.out_adj.contains(species) {
    self.out_adj[species] = []
  }
  if !self.in_adj.contains(species) {
    self.in_adj[species] = []
  }
}

fn _pathway_get_or_empty(m : Map[String, Array[(String, String)]], key : String) -> Array[(String, String)] {
  match m.get(key) {
    Some(arr) => arr
    None => []
  }
}

///|
/// Add a directed edge ``source -> sink`` labelled by ``interaction_id``.
///
/// Both endpoints are auto-registered as nodes.
pub fn Network::add_interaction(
  self : Network,
  source : String,
  sink : String,
  interaction_id : String,
) -> Unit {
  self.add_species(source)
  self.add_species(sink)

  // Append to outgoing.
  let outs = _pathway_get_or_empty(self.out_adj, source)
  outs.push((sink, interaction_id))
  self.out_adj[source] = outs

  // Append to incoming.
  let ins = _pathway_get_or_empty(self.in_adj, sink)
  ins.push((source, interaction_id))
  self.in_adj[sink] = ins
}

///|
/// The set of species in this network.
pub fn Network::species(self : Network) -> Array[String] {
  self.nodes.keys().collect()
}

///|
/// Species that have a directed edge towards the given node
/// (i.e. who *produce* it).
pub fn Network::source(self : Network, species : String) -> Array[String] {
  let ins = _pathway_get_or_empty(self.in_adj, species)
  let result : Array[String] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([])
  let mut i = 0
  while i < ins.length() {
    let s = ins[i].0
    if !seen.contains(s) {
      seen[s] = true
      result.push(s)
    }
    i = i + 1
  }
  result
}

///|
/// Species that this node has a directed edge towards
/// (i.e. that it *produces*).
pub fn Network::sink(self : Network, species : String) -> Array[String] {
  let outs = _pathway_get_or_empty(self.out_adj, species)
  let result : Array[String] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([])
  let mut i = 0
  while i < outs.length() {
    let s = outs[i].0
    if !seen.contains(s) {
      seen[s] = true
      result.push(s)
    }
    i = i + 1
  }
  result
}

///|
/// ``(source, label)`` pairs whose edge produces the given species.
pub fn Network::source_interactions(
  self : Network,
  species : String,
) -> Array[(String, String)] {
  _pathway_get_or_empty(self.in_adj, species)
}

///|
/// ``(sink, label)`` pairs whose edge is produced by the given species.
pub fn Network::sink_interactions(
  self : Network,
  species : String,
) -> Array[(String, String)] {
  _pathway_get_or_empty(self.out_adj, species)
}

// ===========================================================================
// Example glycolysis system (for demos / tests)
// ===========================================================================

///|
/// Build a tiny 3-reaction system reminiscent of glycolysis.
/// Reactions (all irreversible, stoichiometry 1:1):
///
/// ```
///   hexokinase   : glucose + ATP -> G6P + ADP
///   gpi          : G6P           -> F6P
///   pfk          : F6P      + ATP -> F16BP + ADP
/// ```
pub fn create_example_pathway() -> System {
  let r_hk = Reaction::new(
    "hexokinase",
    Map([("glucose", -1.0), ("atp", -1.0), ("g6p", 1.0), ("adp", 1.0)]),
    catalysts=["Hexokinase"],
    reversible=false,
  )
  let r_gpi = Reaction::new(
    "gpi",
    Map([("g6p", -1.0), ("f6p", 1.0)]),
    catalysts=["G6P isomerase"],
    reversible=true,
  )
  let r_pfk = Reaction::new(
    "pfk",
    Map([("f6p", -1.0), ("atp", -1.0), ("f16bp", 1.0), ("adp", 1.0)]),
    catalysts=["Phosphofructokinase"],
    reversible=false,
  )
  System::new("glycolysis_core", [r_hk, r_gpi, r_pfk])
}

///|
/// Build the corresponding Network (S -> A -> B -> P linear chain).
pub fn create_example_network() -> Network {
  let net = Network::new()
  net.add_species("S")
  net.add_species("A")
  net.add_species("B")
  net.add_species("P")
  net.add_interaction("S", "A", "r1")
  net.add_interaction("A", "B", "r2")
  net.add_interaction("B", "P", "r3")
  net
}
