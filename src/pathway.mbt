///|
/// Bio.Pathway - Biochemical pathway analysis.

///|
pub struct Species {
  id : String
  name : String
}

///|
pub fn Species::new(id : String, name : String) -> Species {
  Species::{ id, name }
}

///|
pub struct Reaction {
  id : String
  name : String
  reactants : Array[String]
  products : Array[String]
  reversible : Bool
}

///|
pub fn Reaction::new(id : String, name : String, reactants : Array[String], products : Array[String], reversible : Bool) -> Reaction {
  Reaction::{ id, name, reactants, products, reversible }
}

///|
pub fn Reaction::to_string(self : Reaction) -> String {
  let arrow = if self.reversible { "<->" } else { "->" }
  let reactants_str = self.reactants.join(", ")
  let products_str = self.products.join(", ")
  reactants_str + arrow + products_str
}

///|
pub struct Pathway {
  id : String
  name : String
  species : Array[Species]
  reactions : Array[Reaction]
}

///|
pub fn Pathway::new(id : String, name : String, species : Array[Species], reactions : Array[Reaction]) -> Pathway {
  Pathway::{ id, name, species, reactions }
}

///|
pub fn pathway_find_species(pathway : Pathway, id : String) -> Species? {
  let mut i = 0
  while i < pathway.species.length() {
    if pathway.species[i].id == id {
      return Some(pathway.species[i])
    }
    i = i + 1
  }
  None
}

///|
pub fn pathway_find_reactions_with_species(pathway : Pathway, species_id : String) -> Array[Reaction] {
  let result : Array[Reaction] = Array::new()
  
  let mut i = 0
  while i < pathway.reactions.length() {
    let reaction = pathway.reactions[i]
    
    let mut j = 0
    while j < reaction.reactants.length() {
      if reaction.reactants[j] == species_id {
        result.push(reaction)
        break
      }
      j = j + 1
    }
    
    if j == reaction.reactants.length() {
      j = 0
      while j < reaction.products.length() {
        if reaction.products[j] == species_id {
          result.push(reaction)
          break
        }
        j = j + 1
      }
    }
    
    i = i + 1
  }
  
  result
}

///|
pub fn species_degree(pathway : Pathway, species_id : String) -> Int {
  let reactions = pathway_find_reactions_with_species(pathway, species_id)
  reactions.length()
}

///|
pub fn pathway_count_species(pathway : Pathway) -> Int {
  pathway.species.length()
}

///|
pub fn pathway_count_reactions(pathway : Pathway) -> Int {
  pathway.reactions.length()
}

///|
pub fn create_example_pathway() -> Pathway {
  let glucose = Species::new("glucose", "Glucose")
  let g6p = Species::new("g6p", "Glucose-6-phosphate")
  let f6p = Species::new("f6p", "Fructose-6-phosphate")
  let f16bp = Species::new("f16bp", "Fructose-1,6-bisphosphate")
  let dhap = Species::new("dhap", "Dihydroxyacetone phosphate")
  let g3p = Species::new("g3p", "Glyceraldehyde-3-phosphate")
  
  let species : Array[Species] = Array::new()
  species.push(glucose)
  species.push(g6p)
  species.push(f6p)
  species.push(f16bp)
  species.push(dhap)
  species.push(g3p)
  
  let reactants1 : Array[String] = Array::new()
  reactants1.push("glucose")
  let products1 : Array[String] = Array::new()
  products1.push("g6p")
  let r1 = Reaction::new("r1", "Hexokinase", reactants1, products1, false)
  
  let reactants2 : Array[String] = Array::new()
  reactants2.push("g6p")
  let products2 : Array[String] = Array::new()
  products2.push("f6p")
  let r2 = Reaction::new("r2", "Isomerase", reactants2, products2, true)
  
  let reactants3 : Array[String] = Array::new()
  reactants3.push("f6p")
  let products3 : Array[String] = Array::new()
  products3.push("f16bp")
  let r3 = Reaction::new("r3", "Phosphofructokinase", reactants3, products3, false)
  
  let reactants4 : Array[String] = Array::new()
  reactants4.push("f16bp")
  let products4 : Array[String] = Array::new()
  products4.push("dhap")
  products4.push("g3p")
  let r4 = Reaction::new("r4", "Aldolase", reactants4, products4, true)
  
  let reactants5 : Array[String] = Array::new()
  reactants5.push("dhap")
  let products5 : Array[String] = Array::new()
  products5.push("g3p")
  let r5 = Reaction::new("r5", "Triose-phosphate isomerase", reactants5, products5, true)
  
  let reactions : Array[Reaction] = Array::new()
  reactions.push(r1)
  reactions.push(r2)
  reactions.push(r3)
  reactions.push(r4)
  reactions.push(r5)
  
  Pathway::new("glycolysis", "Glycolysis", species, reactions)
}