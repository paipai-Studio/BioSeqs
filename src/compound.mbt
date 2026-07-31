///|
/// Bio.Compound module - Chemical compound representation.
/// Inspired by Biopython's Bio.KEGG.Compound package.
/// Provides data structures for metabolites and reactions in biological pathways.

///|
/// Represents a chemical compound/metabolite in a biological pathway.
pub struct Compound {
  id : String
  name : String
  formula : String
  charge : Int
  smiles : String
  inchikey : String
  inchi : String
  aliases : Array[String]
  pathways : Array[String]
  reactions : Array[String]
  comment : String
}

///|
/// Create a new compound with minimal information.
pub fn Compound::new(id : String, name : String) -> Compound {
  Compound::{
    id,
    name,
    formula: "",
    charge: 0,
    smiles: "",
    inchikey: "",
    inchi: "",
    aliases: [],
    pathways: [],
    reactions: [],
    comment: "",
  }
}

///|
/// Create a compound with full chemical annotations.
pub fn Compound::with_chemical(
  id : String,
  name : String,
  formula : String,
  charge : Int,
  smiles : String,
) -> Compound {
  Compound::{
    id,
    name,
    formula,
    charge,
    smiles,
    inchikey: "",
    inchi: "",
    aliases: [],
    pathways: [],
    reactions: [],
    comment: "",
  }
}

///|
/// Get compound ID.
pub fn Compound::get_id(self : Compound) -> String {
  self.id
}

///|
/// Get compound name.
pub fn Compound::get_name(self : Compound) -> String {
  self.name
}

///|
/// Get molecular formula.
pub fn Compound::get_formula(self : Compound) -> String {
  self.formula
}

///|
/// Set molecular formula.
pub fn Compound::set_formula(self : Compound, formula : String) -> Compound {
  let c = self
  Compound::{ id: c.id, name: c.name, formula, charge: c.charge, smiles: c.smiles, inchikey: c.inchikey, inchi: c.inchi, aliases: c.aliases, pathways: c.pathways, reactions: c.reactions, comment: c.comment }
}

///|
/// Get charge.
pub fn Compound::get_charge(self : Compound) -> Int {
  self.charge
}

///|
/// Set charge.
pub fn Compound::set_charge(self : Compound, charge : Int) -> Compound {
  let c = self
  Compound::{ id: c.id, name: c.name, formula: c.formula, charge, smiles: c.smiles, inchikey: c.inchikey, inchi: c.inchi, aliases: c.aliases, pathways: c.pathways, reactions: c.reactions, comment: c.comment }
}

///|
/// Get SMILES string.
pub fn Compound::get_smiles(self : Compound) -> String {
  self.smiles
}

///|
/// Set SMILES string.
pub fn Compound::set_smiles(self : Compound, smiles : String) -> Compound {
  let c = self
  Compound::{ id: c.id, name: c.name, formula: c.formula, charge: c.charge, smiles, inchikey: c.inchikey, inchi: c.inchi, aliases: c.aliases, pathways: c.pathways, reactions: c.reactions, comment: c.comment }
}

///|
/// Get aliases (synonyms).
pub fn Compound::get_aliases(self : Compound) -> Array[String] {
  self.aliases
}

///|
/// Add an alias.
pub fn Compound::add_alias(self : Compound, alias : String) -> Compound {
  let c = self
  let new_aliases = c.aliases.copy()
  new_aliases.push(alias)
  Compound::{ id: c.id, name: c.name, formula: c.formula, charge: c.charge, smiles: c.smiles, inchikey: c.inchikey, inchi: c.inchi, aliases: new_aliases, pathways: c.pathways, reactions: c.reactions, comment: c.comment }
}

///|
/// Check if compound has a specific alias.
pub fn Compound::has_alias(self : Compound, name : String) -> Bool {
  let mut i = 0
  let target = name.to_upper()
  while i < self.aliases.length() {
    if self.aliases[i].to_upper() == target {
      return true
    }
    i = i + 1
  }
  self.name.to_upper() == target
}

///|
/// Get pathways this compound participates in.
pub fn Compound::get_pathways(self : Compound) -> Array[String] {
  self.pathways
}

///|
/// Add a pathway.
pub fn Compound::add_pathway(self : Compound, pathway : String) -> Compound {
  let c = self
  let new_pathways = c.pathways.copy()
  new_pathways.push(pathway)
  Compound::{ id: c.id, name: c.name, formula: c.formula, charge: c.charge, smiles: c.smiles, inchikey: c.inchikey, inchi: c.inchi, aliases: c.aliases, pathways: new_pathways, reactions: c.reactions, comment: c.comment }
}

///|
/// Get reactions involving this compound.
pub fn Compound::get_reactions(self : Compound) -> Array[String] {
  self.reactions
}

///|
/// Add a reaction.
pub fn Compound::add_reaction(self : Compound, reaction : String) -> Compound {
  let c = self
  let new_reactions = c.reactions.copy()
  new_reactions.push(reaction)
  Compound::{ id: c.id, name: c.name, formula: c.formula, charge: c.charge, smiles: c.smiles, inchikey: c.inchikey, inchi: c.inchi, aliases: c.aliases, pathways: c.pathways, reactions: new_reactions, comment: c.comment }
}

///|
/// Get a string summary of the compound.
pub fn Compound::to_string(self : Compound) -> String {
  let mut result = "Compound{".to_string()
  result = result + "id: '\{self.id}'"
  result = result + ", name: '\{self.name}'"
  if self.formula.length() > 0 {
    result = result + ", formula: '\{self.formula}'"
  }
  if self.charge != 0 {
    result = result + ", charge: \{self.charge.to_string()}"
  }
  if self.smiles.length() > 0 {
    result = result + ", smiles: '\{self.smiles}'"
  }
  if self.aliases.length() > 0 {
    result = result + ", aliases: ["
    let mut i = 0
    while i < self.aliases.length() {
      if i > 0 { result = result + ", " }
      result = result + "'\{self.aliases[i]}'"
      i = i + 1
    }
    result = result + "]"
  }
  result = result + "}"
  result
}

///|
/// Calculate molecular formula from a string representation.
/// Supports formats like "C6H12O6" (glucose).
pub fn parse_formula(formula : String) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=20)
  let mut i = 0
  let len = formula.length()

  while i < len {
    let c = formula.substring(start=i, end=i + 1)

    if c >= "A" && c <= "Z" {
      let mut atom = c
      i = i + 1
      if i < len {
        let next_c = formula.substring(start=i, end=i + 1)
        if next_c >= "a" && next_c <= "z" {
          atom = atom + next_c
          i = i + 1
        }
      }

      let mut count = ""
      while i < len {
        let digit = formula.substring(start=i, end=i + 1)
        if digit >= "0" && digit <= "9" {
          count = count + digit
          i = i + 1
        } else {
          break
        }
      }

      let n = if count.length() > 0 { parse_int(count) } else { 1 }
      let prev = match counts.get(atom) {
        Some(v) => v
        None => 0
      }
      counts[atom] = prev + n
    } else {
      i = i + 1
    }
  }

  counts
}

///|
/// Get molecular formula as a string from parsed counts.
pub fn formula_to_string(counts : Map[String, Int]) -> String {
  let mut result = ""
  let elements = ["C", "H", "N", "O", "S", "P", "Fe", "Mg", "Ca", "K", "Na", "Cl", "Br", "I"]

  let mut i = 0
  while i < elements.length() {
    let el = elements[i]
    match counts.get(el) {
      Some(n) => {
        if n > 0 {
          result = result + el
          if n > 1 {
            result = result + n.to_string()
          }
        }
      }
      None => ()
    }
    i = i + 1
  }

  // Add any remaining elements not in the predefined list
  let counted = elements.copy()
  for el in counts.keys() {
    let mut found = false
    let mut j = 0
    while j < counted.length() {
      if counted[j] == el {
        found = true
        break
      }
      j = j + 1
    }
    if !found {
      match counts.get(el) {
        Some(n) => {
          if n > 0 {
            result = result + el
            if n > 1 {
              result = result + n.to_string()
            }
          }
        }
        None => ()
      }
    }
  }

  result
}

///|
/// Calculate approximate molecular weight from a formula string.
pub fn compound_molecular_weight(formula : String) -> Double {
  let counts = parse_formula(formula)
  let weights : Map[String, Double] = Map([], capacity=15)
  weights["H"] = 1.008
  weights["C"] = 12.011
  weights["N"] = 14.007
  weights["O"] = 15.999
  weights["F"] = 18.998
  weights["Na"] = 22.990
  weights["P"] = 30.974
  weights["S"] = 32.065
  weights["Cl"] = 35.453
  weights["K"] = 39.098
  weights["Ca"] = 40.078
  weights["Fe"] = 55.845
  weights["Br"] = 79.904
  weights["I"] = 126.904
  weights["Mg"] = 24.305

  let mut total = 0.0
  for el in counts.keys() {
    match counts.get(el) {
      Some(count) => {
        match weights.get(el) {
          Some(w) => total = total + w * count.to_double()
          None => ()
        }
      }
      None => ()
    }
  }

  total
}

///|
/// Represents a biochemical reaction (to avoid conflict with pathway.Reaction).
pub struct ChemicalReaction {
  id : String
  name : String
  substrates : Array[String]
  products : Array[String]
  enzymes : Array[String]
  direction : String
  pathway : String
  comment : String
}

///|
/// Create a new chemical reaction.
pub fn ChemicalReaction::new(id : String, name : String) -> ChemicalReaction {
  ChemicalReaction::{
    id,
    name,
    substrates: [],
    products: [],
    enzymes: [],
    direction: "forward",
    pathway: "",
    comment: "",
  }
}

///|
/// Get reaction ID.
pub fn ChemicalReaction::get_id(self : ChemicalReaction) -> String {
  self.id
}

///|
/// Get reaction name.
pub fn ChemicalReaction::get_name(self : ChemicalReaction) -> String {
  self.name
}

///|
/// Get substrates (reactants).
pub fn ChemicalReaction::get_substrates(self : ChemicalReaction) -> Array[String] {
  self.substrates
}

///|
/// Add a substrate.
pub fn ChemicalReaction::add_substrate(self : ChemicalReaction, substrate : String) -> ChemicalReaction {
  let r = self
  let new_subs = r.substrates.copy()
  new_subs.push(substrate)
  ChemicalReaction::{ id: r.id, name: r.name, substrates: new_subs, products: r.products, enzymes: r.enzymes, direction: r.direction, pathway: r.pathway, comment: r.comment }
}

///|
/// Get products.
pub fn ChemicalReaction::get_products(self : ChemicalReaction) -> Array[String] {
  self.products
}

///|
/// Add a product.
pub fn ChemicalReaction::add_product(self : ChemicalReaction, product : String) -> ChemicalReaction {
  let r = self
  let new_prods = r.products.copy()
  new_prods.push(product)
  ChemicalReaction::{ id: r.id, name: r.name, substrates: r.substrates, products: new_prods, enzymes: r.enzymes, direction: r.direction, pathway: r.pathway, comment: r.comment }
}

///|
/// Get enzymes catalyzing this reaction.
pub fn ChemicalReaction::get_enzymes(self : ChemicalReaction) -> Array[String] {
  self.enzymes
}

///|
/// Add an enzyme.
pub fn ChemicalReaction::add_enzyme(self : ChemicalReaction, enzyme : String) -> ChemicalReaction {
  let r = self
  let new_enzymes = r.enzymes.copy()
  new_enzymes.push(enzyme)
  ChemicalReaction::{ id: r.id, name: r.name, substrates: r.substrates, products: r.products, enzymes: new_enzymes, direction: r.direction, pathway: r.pathway, comment: r.comment }
}

///|
/// Get reaction direction.
pub fn ChemicalReaction::get_direction(self : ChemicalReaction) -> String {
  self.direction
}

///|
/// Set reaction direction ("forward", "reverse", "reversible").
pub fn ChemicalReaction::set_direction(self : ChemicalReaction, direction : String) -> ChemicalReaction {
  let r = self
  ChemicalReaction::{ id: r.id, name: r.name, substrates: r.substrates, products: r.products, enzymes: r.enzymes, direction, pathway: r.pathway, comment: r.comment }
}

///|
/// Get reaction equation as a string.
pub fn ChemicalReaction::get_equation(self : ChemicalReaction) -> String {
  let mut lhs = ""
  let mut i = 0
  while i < self.substrates.length() {
    if i > 0 { lhs = lhs + " + " }
    lhs = lhs + self.substrates[i]
    i = i + 1
  }

  let mut rhs = ""
  i = 0
  while i < self.products.length() {
    if i > 0 { rhs = rhs + " + " }
    rhs = rhs + self.products[i]
    i = i + 1
  }

  let arrow = match self.direction {
    "forward" => " -> "
    "reverse" => " <- "
    _ => " <=> "
  }

  lhs + arrow + rhs
}

///|
/// Get a string summary of the reaction.
pub fn ChemicalReaction::to_string(self : ChemicalReaction) -> String {
  let mut result = "ChemicalReaction{\{self.id}: \{self.name}\n".to_string()
  result = result + "  Equation: \{self.get_equation()}\n"
  if self.enzymes.length() > 0 {
    result = result + "  Enzymes: "
    let mut i = 0
    while i < self.enzymes.length() {
      if i > 0 { result = result + ", " }
      result = result + self.enzymes[i]
      i = i + 1
    }
    result = result + "\n"
  }
  if self.pathway.length() > 0 {
    result = result + "  Pathway: \{self.pathway}\n"
  }
  result = result + "}"
  result
}

///|
/// A collection of compounds and reactions representing a biochemical pathway.
pub struct CompoundPathwayMap {
  name : String
  compounds : Array[Compound]
  reactions : Array[ChemicalReaction]
  description : String
}

///|
/// Create a new pathway map.
pub fn CompoundPathwayMap::new(name : String) -> CompoundPathwayMap {
  CompoundPathwayMap::{ name, compounds: [], reactions: [], description: "" }
}

///|
/// Get pathway name.
pub fn CompoundPathwayMap::get_name(self : CompoundPathwayMap) -> String {
  self.name
}

///|
/// Add a compound to the pathway.
pub fn CompoundPathwayMap::add_compound(self : CompoundPathwayMap, compound : Compound) -> CompoundPathwayMap {
  let p = self
  let new_compounds = p.compounds.copy()
  new_compounds.push(compound)
  CompoundPathwayMap::{ name: p.name, compounds: new_compounds, reactions: p.reactions, description: p.description }
}

///|
/// Add a reaction to the pathway.
pub fn CompoundPathwayMap::add_reaction(self : CompoundPathwayMap, reaction : ChemicalReaction) -> CompoundPathwayMap {
  let p = self
  let new_reactions = p.reactions.copy()
  new_reactions.push(reaction)
  CompoundPathwayMap::{ name: p.name, compounds: p.compounds, reactions: new_reactions, description: p.description }
}

///|
/// Get all compounds.
pub fn CompoundPathwayMap::get_compounds(self : CompoundPathwayMap) -> Array[Compound] {
  self.compounds
}

///|
/// Get all reactions.
pub fn CompoundPathwayMap::get_reactions(self : CompoundPathwayMap) -> Array[ChemicalReaction] {
  self.reactions
}

///|
/// Find a compound by ID.
pub fn CompoundPathwayMap::find_compound(self : CompoundPathwayMap, id : String) -> Option[Compound] {
  let mut i = 0
  while i < self.compounds.length() {
    if self.compounds[i].id == id {
      return Some(self.compounds[i])
    }
    i = i + 1
  }
  None
}

///|
/// Find a reaction by ID.
pub fn CompoundPathwayMap::find_reaction(self : CompoundPathwayMap, id : String) -> Option[ChemicalReaction] {
  let mut i = 0
  while i < self.reactions.length() {
    if self.reactions[i].id == id {
      return Some(self.reactions[i])
    }
    i = i + 1
  }
  None
}

///|
/// Generate a pathway summary.
pub fn CompoundPathwayMap::summary(self : CompoundPathwayMap) -> String {
  let mut result = "Pathway: \{self.name}\n".to_string()
  result = result + "Compounds: \{self.compounds.length()}\n"
  result = result + "Reactions: \{self.reactions.length()}\n"

  if self.compounds.length() > 0 {
    result = result + "\nCompounds:\n"
    let mut i = 0
    while i < self.compounds.length() {
      result = result + "  - \{self.compounds[i].id}: \{self.compounds[i].name}"
      if self.compounds[i].formula.length() > 0 {
        result = result + " (\{self.compounds[i].formula})"
      }
      result = result + "\n"
      i = i + 1
    }
  }

  if self.reactions.length() > 0 {
    result = result + "\nReactions:\n"
    let mut i = 0
    while i < self.reactions.length() {
      result = result + "  - \{self.reactions[i].id}: \{self.reactions[i].get_equation()}\n"
      i = i + 1
    }
  }

  result
}

///|
/// Create common known compound templates.
pub fn compound_glucose() -> Compound {
  Compound::with_chemical("C0003", "D-Glucose", "C6H12O6", 0, "C(C1C(C(C(C(O1)O)O)O)O)O")
}

///|
/// Create fructose compound.
pub fn compound_fructose() -> Compound {
  Compound::with_chemical("C0005", "D-Fructose", "C6H12O6", 0, "C(=O)C1(C(C(C(O1)O)O)O)O")
}

///|
/// Create pyruvate compound.
pub fn compound_pyruvate() -> Compound {
  Compound::with_chemical("C0024", "Pyruvate", "C3H3O3", -1, "CC(=O)C([O-])=O")
}

///|
/// Create acetate compound.
pub fn compound_acetate() -> Compound {
  Compound::with_chemical("C0009", "Acetate", "C2H3O2", -1, "CC(=O)[O-]")
}

///|
/// Create citrate compound.
pub fn compound_citrate() -> Compound {
  Compound::with_chemical("C0015", "Citrate", "C6H5O7", -3, "OC(=O)C(CC(=O)[O-])(O)CC(=O)[O-]")
}
