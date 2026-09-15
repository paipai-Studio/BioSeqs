///|
/// abstract_property_map demo - Unified residue property access.
/// Inspired by Biopython Bio.PDB.AbstractPropertyMap.

fn main {
  println("=== abstract_property_map Demo ===")
  println("")

  // Step 1: Build property map
  println("--- Step 1: Build Property Map ---")
  let pm = @src.abstract_property_map_example()
  println("Residues: " + pm.n_residues.to_string())
  println("Properties: " + pm.property_names.to_string())
  println("")

  // Step 2: Query properties
  println("--- Step 2: Query Properties ---")
  let ss = pm.get_property("A", 1, "SS")
  let asa = pm.get_property("A", 1, "ASA")
  println("A:1 SS = " + ss.unwrap_or("NA"))
  println("A:1 ASA = " + asa.unwrap_or("NA"))
  println("Has A:1 = " + pm.has_residue("A", 1).to_string())
  println("Has B:1 = " + pm.has_residue("B", 1).to_string())
  println("")

  // Step 3: Filter
  println("--- Step 3: Filter by Property ---")
  let helix = pm.filter_by_property("SS", fn(s) { s == "H" })
  println("Helix residues: " + helix.n_residues.to_string())
  let sheet = pm.filter_by_property("SS", fn(s) { s == "E" })
  println("Sheet residues: " + sheet.n_residues.to_string())
  println("")

  // Step 4: Statistics
  println("--- Step 4: Statistics ---")
  let (count, min_v, max_v, mean, sum) = pm.stats_property("ASA")
  println("ASA: count=" + count.to_string() + " min=" + min_v.to_string() + " max=" + max_v.to_string())
  println("     mean=" + mean.to_string() + " sum=" + sum.to_string())
  println("")

  // Step 5: Export TSV
  println("--- Step 5: TSV Export ---")
  println(pm.to_tsv())
  println("")

  // Step 6: Merge
  println("--- Step 6: Merge Maps ---")
  let pm2 = @src.AbstractPropertyMap::new(["Phi", "Psi"])
  let props : Map[String, String] = {
    let m : Map[String, String] = Map([], capacity=2)
    m["Phi"] = "-60.0"
    m["Psi"] = "-45.0"
    m
  }
  pm2.add("A", 1, props)
  let merged = @src.AbstractPropertyMap::merge(pm, pm2)
  let phi = merged.get_property("A", 1, "Phi")
  println("Merged A:1 Phi = " + phi.unwrap_or("NA"))
  println("Merged residues: " + merged.n_residues.to_string())
  println("")

  println("=== abstract_property_map Demo Complete ===")
}