///|
/// Phylogenetic tree example: parse Newick, traverse, and analyze.
///
/// Run: moon run examples/phylo_demo/main.mbt
fn main {
  // === Parse Newick tree ===
  let newick = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
  let tree = @bio.parse_newick(newick)

  println("=== Newick Tree ===")
  println("  Input: \{newick}")
  println("  Rooted: \{tree.is_rooted()}")

  // === Tree statistics ===
  let terminals = tree.count_terminals()
  println("\n=== Tree Statistics ===")
  println("  Terminal nodes (leaves): \{terminals}")

  let nonterminals = tree.get_nonterminals()
  println("  Non-terminal nodes: \{nonterminals.length()}")

  let is_bifurcating = tree.is_bifurcating()
  println("  Is bifurcating: \{is_bifurcating}")

  // === Get terminal clades ===
  let leaf_clades = tree.get_terminals()
  println("\n=== Terminal Clades ===")
  for i = 0; i < leaf_clades.length(); i = i + 1 {
    let name = leaf_clades[i].name
    match name {
      Some(n) => println("  Leaf \{i}: \{n}")
      None => println("  Leaf \{i}: (unnamed)")
    }
  }

  // === Find clade by name ===
  let found = tree.find_any(name=Some("A"))
  match found {
    Some(c) => {
      let name = c.name
      match name {
        Some(n) => println("\n=== Find Clade 'A' ===\n  Found: \{n}")
        None => println("\n=== Find Clade 'A' ===\n  Found (unnamed)")
      }
    }
    None => println("\n=== Find Clade 'A' ===\n  Not found")
  }

  // === Get path between two taxa ===
  let path = tree.get_path("C")
  match path {
    Some(p) => {
      println("\n=== Path to C ===")
      for i = 0; i < p.length(); i = i + 1 {
        let name = p[i].name
        match name {
          Some(n) => println("  [\{i}] \{n}")
          None => println("  [\{i}] (internal)")
        }
      }
    }
    None => println("\n=== Path to C ===\n  No path found")
  }

  // === Calculate tree depths ===
  let depths = tree.depths()
  println("\n=== Clade Depths ===")
  for i = 0; i < depths.length(); i = i + 1 {
    let (name, depth) = depths[i]
    println("  \{name}: \{depth}")
  }

  // === Distance between two taxa ===
  let dist = tree.distance("A", "B")
  match dist {
    Some(d) => println("\n=== Distance A <-> B ===\n  \{d}")
    None => println("\n=== Distance A <-> B ===\n  N/A")
  }

  // === Common ancestor ===
  let ancestor = tree.common_ancestor(["A", "B"])
  match ancestor {
    Some(c) => {
      let name = c.name
      match name {
        Some(n) => println("\n=== Common Ancestor of A, B ===\n  \{n}")
        None => println("\n=== Common Ancestor of A, B ===\n  (internal node)")
      }
    }
    None => println("\n=== Common Ancestor of A, B ===\n  Not found")
  }

  // === Write Newick ===
  let newick_out = @bio.write_newick(tree)
  println("\n=== Serialized Newick ===")
  println("  \{newick_out}")

  // === Draw ASCII tree ===
  let ascii = tree.draw_ascii()
  println("\n=== ASCII Tree ===")
  println(ascii)
}
