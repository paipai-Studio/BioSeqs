/// Demo: Bio.Pathway — a tiny glycolysis-like system.
///
/// Run with: moon run examples/pathway_demo

fn main {
  // ---- System with 3 reactions -------------------------------------
  let sys = @src.create_example_pathway()
  println("=== System: " + sys.id + " ===")
  let mut i = 0
  while i < sys.reactions.length() {
    let rxn = sys.reactions[i]
    let rev = if rxn.reversible { " (reversible)" } else { "" }
    println("  [" + rxn.id + "] " + rxn.to_string() + rev)
    i = i + 1
  }

  // ---- Stoichiometry matrix ---------------------------------------
  let (species, _rxns, stoch) = sys.stochiometry()
  println("")
  println("=== Stoichiometry matrix ===")
  let mut hdr = "  species"
  let mut j = 0
  while j < sys.reactions.length() {
    hdr = hdr + " | " + sys.reactions[j].id
    j = j + 1
  }
  println(hdr)

  let mut k = 0
  while k < species.length() {
    let mut row = "  " + species[k]
    let mut j2 = 0
    while j2 < stoch[k].length() {
      let v = stoch[k][j2]
      let cell = v.to_string()
      row = row + " | " + cell
      j2 = j2 + 1
    }
    println(row)
    k = k + 1
  }

  // ---- Network: source / sink -------------------------------------
  println("")
  println("=== Network (linear S -> A -> B -> P chain) ===")
  let net = @src.create_example_network()
  let nodes = net.species()
  nodes.sort_by(fn(a : String, b : String) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  let mut m = 0
  while m < nodes.length() {
    let sp = nodes[m]
    let src = net.source(sp)
    let snk = net.sink(sp)
    let src_s = if src.length() == 0 { "none" } else { src.join(",") }
    let snk_s = if snk.length() == 0 { "none" } else { snk.join(",") }
    println("  source(" + sp + ") = [" + src_s + "]   sink(" + sp + ") = [" + snk_s + "]")
    m = m + 1
  }

  // ---- BFS / DFS on PathwayMultiGraph -----------------------------
  println("")
  println("=== Graph traversals (MultiGraph A->B->C->D->E + shortcut A->C) ===")
  let mg = @src.PathwayMultiGraph::new()
  mg.add_node("A"); mg.add_node("B"); mg.add_node("C"); mg.add_node("D"); mg.add_node("E")
  mg.add_edge("A", "B", "e1")
  mg.add_edge("B", "C", "e2")
  mg.add_edge("C", "D", "e3")
  mg.add_edge("D", "E", "e4")
  mg.add_edge("A", "C", "e5")

  let bfs_order = @src.bf_search(mg, root="A")
  let dfs_order = @src.df_search(mg, root="A")
  println("  BFS from A: " + bfs_order.join(" -> "))
  println("  DFS from A: " + dfs_order.join(" -> "))

  println("")
  println("=== Done ===")
}
