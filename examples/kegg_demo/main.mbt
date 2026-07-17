///|
fn main {
  println("=== KEGG Database Demo ===")

  let gene_text = "ENTRY       hsa:12345\nNAME        GeneA\nDEFINITION  Description of gene A\nORTHOLOGY   K00001\nORGANISM    Homo sapiens\nPATHWAY     path:hsa00010 Glycolysis / Gluconeogenesis\nPATHWAY     path:hsa00020 Citrate cycle\nPOSITION    1p36.11\nDBLINKS     NCBI-GeneID: 12345\nDBLINKS     UniProt: P12345"

  println("\n1. Parse KEGG Gene Record:")
  println("==========================")
  let gene = @bio.parse_kegg_gene(gene_text)
  if gene.is_some() {
    let g = gene.unwrap()
    println("Entry: " + g.entry)
    println("Name: " + g.name)
    println("Definition: " + g.definition)
    println("Orthology: " + g.orthology)
    println("Organism: " + g.organism)
    println("Position: " + g.position)
    println("Pathways: " + g.pathways.length().to_string())

    let mut i = 0
    while i < g.pathways.length() {
      println("  - " + g.pathways[i])
      i = i + 1
    }

    println("DB Links: " + g.dblinks.length().to_string())
    i = 0
    while i < g.dblinks.length() {
      println("  - " + g.dblinks[i].0 + ": " + g.dblinks[i].1)
      i = i + 1
    }

    let pathway_count = @bio.kegg_gene_pathway_count(g)
    println("\nPathway count: " + pathway_count.to_string())
  }

  let pathway_text = "ENTRY       path:hsa00010\nNAME        Glycolysis / Gluconeogenesis\nDEFINITION  Glycolysis is the metabolic pathway that converts glucose...\nCLASS       Carbohydrate metabolism\nGENE        hsa:12345  GeneA\nGENE        hsa:67890  GeneB\nGENE        hsa:54321  GeneC\nCOMPOUND    C00033 Glucose\nCOMPOUND    C00084 Pyruvate\nENZYME      EC:2.7.1.1\nENZYME      EC:1.2.1.12"

  println("\n2. Parse KEGG Pathway Record:")
  println("==============================")
  let pathway = @bio.parse_kegg_pathway(pathway_text)
  if pathway.is_some() {
    let p = pathway.unwrap()
    println("Entry: " + p.entry)
    println("Name: " + p.name)
    println("Definition: " + p.definition)
    println("Class: " + p.class)
    println("Genes: " + p.genes.length().to_string())

    let mut i = 0
    while i < p.genes.length() {
      println("  - " + p.genes[i].0 + " (" + p.genes[i].1 + ")")
      i = i + 1
    }

    println("Compounds: " + p.compounds.length().to_string())
    i = 0
    while i < p.compounds.length() {
      println("  - " + p.compounds[i])
      i = i + 1
    }

    println("Enzymes: " + p.enzymes.length().to_string())
    i = 0
    while i < p.enzymes.length() {
      println("  - " + p.enzymes[i])
      i = i + 1
    }

    let gene_ids = @bio.kegg_pathway_gene_ids(p)
    println("\nGene IDs extracted:")
    i = 0
    while i < gene_ids.length() {
      println("  " + gene_ids[i])
      i = i + 1
    }
  }

  let compound_text = "ENTRY       C00033\nNAME        Glucose\nFORMULA     C6H12O6\nMASS        180.156\nPATHWAY     path:hsa00010 Glycolysis / Gluconeogenesis\nPATHWAY     path:hsa00052 Galactose metabolism\nENZYME      EC:2.7.1.1\nENZYME      EC:3.2.1.20\nREACTION    R01786\nREACTION    R02102"

  println("\n3. Parse KEGG Compound Record:")
  println("===============================")
  let compound = @bio.parse_kegg_compound(compound_text)
  if compound.is_some() {
    let c = compound.unwrap()
    println("Entry: " + c.entry)
    println("Name: " + c.name)
    println("Formula: " + c.formula)
    println("Mass: " + c.mass.to_string())
    println("Pathways: " + c.pathways.length().to_string())

    let mut i = 0
    while i < c.pathways.length() {
      println("  - " + c.pathways[i])
      i = i + 1
    }

    println("Enzymes: " + c.enzymes.length().to_string())
    i = 0
    while i < c.enzymes.length() {
      println("  - " + c.enzymes[i])
      i = i + 1
    }

    println("Reactions: " + c.reactions.length().to_string())
    i = 0
    while i < c.reactions.length() {
      println("  - " + c.reactions[i])
      i = i + 1
    }

    let mut found = false
    i = 0
    while i < c.pathways.length() {
      if c.pathways[i].has_prefix("path:hsa00010") {
        found = true
      }
      i = i + 1
    }
    println("\nGlucose is in pathway hsa00010: " + found.to_string())
  }

  let enzyme_text = "ENTRY       EC:2.7.1.1\nNAME        Hexokinase\nCLASS       Transferases; Transferring phosphorus-containing groups\nGENE        hsa:1234  HK1\nGENE        hsa:1235  HK2\nGENE        hsa:1236  HK3\nPATHWAY     path:hsa00010\nPATHWAY     path:hsa00051\nREACTION    R01786\nCOFACTOR    ATP\nCOFACTOR    Mg2+"

  println("\n4. Parse KEGG Enzyme Record:")
  println("=============================")
  let enzyme = @bio.parse_kegg_enzyme(enzyme_text)
  if enzyme.is_some() {
    let e = enzyme.unwrap()
    println("Entry: " + e.entry)
    println("Name: " + e.name)
    println("Class: " + e.class)
    println("Genes: " + e.genes.length().to_string())

    let mut i = 0
    while i < e.genes.length() {
      println("  - " + e.genes[i].0 + " (" + e.genes[i].1 + ")")
      i = i + 1
    }

    println("Pathways: " + e.pathways.length().to_string())
    i = 0
    while i < e.pathways.length() {
      println("  - " + e.pathways[i])
      i = i + 1
    }

    println("Reactions: " + e.reactions.length().to_string())
    i = 0
    while i < e.reactions.length() {
      println("  - " + e.reactions[i])
      i = i + 1
    }

    println("Cofactors: " + e.cofactors.length().to_string())
    i = 0
    while i < e.cofactors.length() {
      println("  - " + e.cofactors[i])
      i = i + 1
    }
  }

  println("\n=== Demo Complete ===")
}
