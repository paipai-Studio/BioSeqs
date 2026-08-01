///|
fn main {
  println("=== KGML (KEGG Markup Language) Demo ===")
  println("Parsing a KEGG pathway map in KGML XML format.")
  println("")

  // 1. Parse sample KGML pathway
  println("1. Parse KGML Pathway:")
  println("======================")
  let xml = @bio.kgml_sample_pathway()
  let pw_opt = @bio.parse_kgml(xml)
  if pw_opt.is_some() {
    let pw = pw_opt.unwrap()
    println("Name: " + pw.name())
    println("Organism: " + pw.org())
    println("Number: " + pw.number())
    println("Title: " + pw.title())
    println("Image: " + pw.image())
    println("Link: " + pw.link())
    println("Entries: " + pw.entries().length().to_string())
    println("Relations: " + pw.relations().length().to_string())
    println("Reactions: " + pw.reactions().length().to_string())
  }

  // 2. Browse entries
  println("\n2. Browse Pathway Entries:")
  println("==========================")
  if pw_opt.is_some() {
    let pw = pw_opt.unwrap()
    for e in pw.entries() {
      println("Entry " + e.id().to_string() + ": " + e.name() + " (type=" + e.etype() + ")")
      if e.reaction().length() > 0 {
        println("  Reaction: " + e.reaction())
      }
      let g = e.graphics()
      if g.is_some() {
        let gfx = g.unwrap()
        println("  Graphics: type=" + gfx.gtype() + " x=" + gfx.x().to_string() + " y=" + gfx.y().to_string())
        println("  Colors: fg=" + gfx.fgcolor() + " bg=" + gfx.bgcolor())
      }
    }
  }

  // 3. Browse relations
  println("\n3. Browse Relations:")
  println("====================")
  if pw_opt.is_some() {
    let pw = pw_opt.unwrap()
    for r in pw.relations() {
      println("Relation: " + r.entry1().to_string() + " -> " + r.entry2().to_string() + " (type=" + r.rtype() + ")")
      for st in r.subtypes() {
        println("  Subtype: " + st.name() + " = " + st.value())
      }
    }
  }

  // 4. Browse reactions
  println("\n4. Browse Reactions:")
  println("====================")
  if pw_opt.is_some() {
    let pw = pw_opt.unwrap()
    for rxn in pw.reactions() {
      println("Reaction: " + rxn.name() + " (type=" + rxn.rtype() + ")")
      for s in rxn.substrates() {
        println("  Substrate: " + s)
      }
      for p in rxn.products() {
        println("  Product: " + p)
      }
    }
  }

  // 5. Query entries by type
  println("\n5. Query Entries by Type:")
  println("=========================")
  if pw_opt.is_some() {
    let pw = pw_opt.unwrap()
    let genes = pw.get_entries_by_type("gene")
    println("Gene entries: " + genes.length().to_string())
    let compounds = pw.get_entries_by_type("compound")
    println("Compound entries: " + compounds.length().to_string())
    let maps = pw.get_entries_by_type("map")
    println("Map entries: " + maps.length().to_string())
  }

  // 6. Serialize back to KGML
  println("\n6. Serialize Pathway to KGML:")
  println("=============================")
  if pw_opt.is_some() {
    let pw = pw_opt.unwrap()
    let output = @bio.kgml_to_string(pw)
    println("Serialized length: " + output.length().to_string() + " chars")
    // Print first 200 chars
    let preview_len = if output.length() > 200 { 200 } else { output.length() }
    println("Preview: " + output[0:preview_len].to_owned() + "...")
  }

  println("\n=== Demo Complete ===")
}
