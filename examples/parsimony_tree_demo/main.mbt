///|
/// Maximum parsimony tree construction demo (Biopython
/// ``Bio.Phylo.TreeConstruction`` equivalents): ParsimonyScorer (Fitch and
/// Sankoff), NNITreeSearcher hill climbing, ParsimonyTreeConstructor with
/// the default identity-UPGMA start, and Tree.root_at_midpoint.
///
/// Run: moon run examples/parsimony_tree_demo
fn main raise {
  println("=== Bio.Phylo.TreeConstruction Parsimony Demo ===")
  println("")

  // -------------------------------------------------------------------------
  // 1. Five-taxon nucleotide alignment
  // -------------------------------------------------------------------------
  let pairs = [
    ("Alpha", "AACGTGGCCACAT"),
    ("Beta", "AAGGTCGCCACAC"),
    ("Gamma", "CAGTTCGCCACAA"),
    ("Delta", "GAGATTTCCGCCT"),
    ("Epsilon", "GAGATCTCCGCCC"),
  ]
  let alignment = make_map(pairs)
  println("--- 1. Input alignment (5 taxa, 13 columns) ---")
  for p in pairs {
    println("  \{p.0}: \{p.1}")
  }
  println("")

  // Neighbor-joining start tree (unrooted), as built by Biopython.
  let calculator = @bio.DistanceCalculator::new(model="identity")
  let dm = calculator.get_distance(pairs)
  let names = pairs.map(fn(kv) { kv.0 })
  let nj = @bio.Tree::new(
    @bio.TreeConstructor::new(tree_method="nj").build_tree(dm, names),
    rooted=false,
  )
  println("--- 2. Unrooted NJ start tree ---")
  println("  " + nj.root.to_newick())
  println("")

  // Fitch parsimony: the scorer roots unrooted trees at their midpoint.
  let fitch = @bio.ParsimonyScorer::new()
  let midpoint = nj.root_at_midpoint()
  println("--- 3. Automatic midpoint rooting before scoring ---")
  println("  " + midpoint.root.to_newick())
  println("  Fitch parsimony score: \{fitch.get_score(nj, alignment)}")
  println("")

  // NNI neighbors and the hill-climbing search.
  println("--- 4. NNI hill climbing (Fitch) ---")
  let neighbors = @bio.nni_get_neighbors(midpoint)
  let neighbor_scores : Array[String] = []
  for n in neighbors {
    neighbor_scores.push(fitch.get_score(n, alignment).to_string())
  }
  println("  neighbor scores: [" + neighbor_scores.join(",") + "]")
  let result = @bio.NNITreeSearcher::new(fitch).search(nj, alignment)
  println("  final tree: " + result.root.to_newick())
  println("  final score: \{fitch.get_score(result, alignment)}")
  println("")

  // Default constructor: identity-UPGMA start, then NNI search.
  println("--- 5. ParsimonyTreeConstructor (default UPGMA start) ---")
  let ctor = @bio.ParsimonyTreeConstructor::new(
    @bio.NNITreeSearcher::new(@bio.ParsimonyScorer::new()),
  )
  let built = ctor.build_tree(alignment)
  println("  " + built.root.to_newick())
  println("  Fitch score: \{fitch.get_score(built, alignment)}")
  println("")

  // -------------------------------------------------------------------------
  // 2. Six-taxon Sankoff parsimony (transitions cost 1, transversions 5)
  // -------------------------------------------------------------------------
  println("--- 6. Sankoff NNI search with a custom cost matrix ---")
  let sankoff_alignment = make_map([
    ("S1", "AAAAA"),
    ("S2", "AAGAA"),
    ("S3", "GTATT"),
    ("S4", "GTAAT"),
    ("S5", "TTTTT"),
    ("S6", "TTTGT"),
  ])
  let start = @bio.Tree::new(
    node([
      node([
        node([
          node([leaf("S1", 0.11), leaf("S3", 0.13)], 0.31),
          leaf("S4", 0.14),
        ], 0.32),
        leaf("S5", 0.15),
      ], 0.33),
      node([leaf("S6", 0.16), leaf("S2", 0.12)], 0.34),
    ], 0.02),
    rooted=true,
  )
  let matrix = @bio.ParsimonyMatrix::new(
    ["A", "C", "G", "T"],
    [
      [0.0, 5.0, 1.0, 5.0],
      [5.0, 0.0, 5.0, 1.0],
      [1.0, 5.0, 0.0, 5.0],
      [5.0, 1.0, 5.0, 0.0],
    ],
  )
  let sankoff = @bio.ParsimonyScorer::new(matrix=Some(matrix))
  println("  start: " + start.root.to_newick())
  println("  Fitch score:   \{fitch.get_score(start, sankoff_alignment)}")
  println("  Sankoff score: \{sankoff.get_score(start, sankoff_alignment)}")

  let mut work = start
  let mut round = 1
  let mut done = false
  while !done {
    let current = sankoff.get_score(work, sankoff_alignment)
    let round_neighbors = @bio.nni_get_neighbors(work)
    let scores : Array[String] = []
    let mut best_score = current
    let mut accepted = -1
    for i in 0..<round_neighbors.length() {
      let s = sankoff.get_score(round_neighbors[i], sankoff_alignment)
      scores.push(s.to_string())
      if s < best_score {
        best_score = s
        accepted = i
        work = round_neighbors[i]
      }
    }
    println("  round \{round}: scores=[" + scores.join(",") + "]")
    if accepted < 0 {
      println("  no strictly better neighbor: search stops")
      done = true
    } else {
      println("    accept neighbor \{accepted}: " + work.root.to_newick())
    }
    round = round + 1
  }
  println("  final Sankoff score: \{sankoff.get_score(work, sankoff_alignment)}")
  println("  final Fitch score:   \{fitch.get_score(work, sankoff_alignment)}")
  println("")

  // -------------------------------------------------------------------------
  // 3. Midpoint rooting of a simple unrooted star
  // -------------------------------------------------------------------------
  println("--- 7. Midpoint rooting of an unrooted 3-taxon tree ---")
  let star = @bio.Tree::new(
    @bio.Clade::new(clades=[
      leaf("A", 1.0),
      leaf("B", 2.0),
      leaf("C", 4.0),
    ]),
    rooted=false,
  )
  println("  unrooted: " + star.root.to_newick())
  println("  midpoint: " + star.root_at_midpoint().root.to_newick())
}

fn make_map(pairs : Array[(String, String)]) -> Map[String, String] {
  let m : Map[String, String] = Map([])
  for p in pairs {
    m[p.0] = p.1
  }
  m
}

fn leaf(name : String, bl : Double) -> @bio.Clade {
  @bio.Clade::new(name=Some(name), branch_length=Some(bl))
}

fn node(children : Array[@bio.Clade], bl : Double) -> @bio.Clade {
  @bio.Clade::new(branch_length=Some(bl), clades=children)
}
