///|
/// DECIPHER demo - Biological sequence analysis toolkit.

fn main {
  println("=== DECIPHER Demo: Biological Sequence Analysis ===")
  println("")

  // Step 1: Sequence Alignment
  println("--- Step 1: Pairwise Alignment ---")
  let seq1 = "ACGTACGTACGTACGTACGT"
  let seq2 = "ACGTACGTACGTACGTACGA"
  let aln = @src.align_seqs(seq1, seq2)
  println("Sequence 1: " + seq1)
  println("Sequence 2: " + seq2)
  println("Aligned 1: " + aln.aligned_a)
  println("Aligned 2: " + aln.aligned_b)
  println("Score: " + aln.score.to_string())
  println("Identity: " + (aln.identity * 100.0).to_string() + "%")
  println("Mismatches: " + aln.mismatches.to_string())
  println("Gaps: " + aln.gaps.to_string())
  println("")

  // Step 2: Distance Matrix
  println("--- Step 2: Distance Matrix ---")
  let seqs = @src.create_example_sequences()
  let dist = @src.decipher_distance_matrix(seqs)
  println("Sequences: " + seqs.length().to_string())
  println("Distance matrix (" + dist.length().to_string() + "x" + dist.length().to_string() + "):")
  let mut i = 0
  while i < dist.length() {
    let sb = StringBuilder()
    sb.write_string("  [")
    let mut j = 0
    while j < dist[i].length() {
      if j > 0 {
        sb.write_string(", ")
      }
      sb.write_string(dist[i][j].to_string())
      j = j + 1
    }
    sb.write_string("]")
    println(sb.to_string())
    i = i + 1
  }
  println("")

  // Step 3: Sequence Clustering
  println("--- Step 3: Identity Clustering (threshold=0.9) ---")
  let clusters = @src.id_clusters(seqs, threshold=0.9)
  println("Number of clusters: " + clusters.n_clusters.to_string())
  let mut ci = 0
  while ci < clusters.n_clusters {
    let members = clusters.clusters[ci]
    let sb = StringBuilder()
    sb.write_string("  Cluster " + (ci + 1).to_string() + " (size=" + members.length().to_string() + "): ")
    let mut mi = 0
    while mi < members.length() {
      if mi > 0 {
        sb.write_string(", ")
      }
      sb.write_string("Seq" + members[mi].to_string())
      mi = mi + 1
    }
    println(sb.to_string())
    ci = ci + 1
  }
  println("")

  // Step 4: Phylogenetic Tree
  println("--- Step 4: Neighbor-Joining Tree ---")
  let names = ["Seq0", "Seq1", "Seq2", "Seq3", "Seq4", "Seq5"]
  let tree = @src.neighbor_joining(dist, names)
  let newick = @src.tree_to_newick(tree)
  println("Newick format:")
  println("  " + newick)
  println("")

  // Step 5: Sequence Database
  println("--- Step 5: Sequence Database ---")
  let db = @src.create_example_database()
  println("Database size: " + @src.db_size(db).to_string())
  let mut di = 0
  while di < @src.db_size(db) {
    let entry = @src.db_get(db, di)
    match entry {
      Some(e) => println("  [" + e.id.to_string() + "] " + e.description + ": " + e.sequence)
      None => println("  [empty]")
    }
    di = di + 1
  }
  println("")

  // Step 6: Database Search
  println("--- Step 6: Database Search ---")
  let query = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  let (idx, search_aln) = @src.db_search(db, query)
  println("Query: " + query)
  println("Best match index: " + idx.to_string())
  println("Identity: " + (search_aln.identity * 100.0).to_string() + "%")
  println("")

  // Step 7: Consensus Sequence
  println("--- Step 7: Consensus Sequence ---")
  let aligned_seqs = [
    "ACGTACGTACGT",
    "ACGTACGTACGT",
    "ACGTACGTACGA",
    "ACGTACGTACGT",
  ]
  let cons = @src.decipher_consensus(aligned_seqs)
  println("Aligned sequences:")
  let mut si = 0
  while si < aligned_seqs.length() {
    println("  " + aligned_seqs[si])
    si = si + 1
  }
  println("Consensus: " + cons)
  println("")

  // Step 8: Chimera Detection
  println("--- Step 8: Chimera Detection ---")
  let parent_a = "AAAAAAATTTTTTT"
  let parent_b = "GGGGGGGCCCCCCC"
  let chimera = "AAAAAAACCCCCCC"
  let all_seqs = [parent_a, parent_b, chimera]
  let abundances = [100, 100, 10]
  let chimera_flags = @src.find_chimeras(all_seqs, abundances, threshold=0.8)
  println("Sequences:")
  let mut chi = 0
  while chi < all_seqs.length() {
    let status = if chimera_flags[chi] { "CHIMERA" } else { "clean" }
    println("  " + all_seqs[chi] + " (abund=" + abundances[chi].to_string() + ") -> " + status)
    chi = chi + 1
  }
  println("")

  // Step 9: GC Content and Complexity
  println("--- Step 9: Sequence Statistics ---")
  let test_seqs = ["ACGTACGT", "GGCCGGCC", "AAAATTTT", "ATATATAT"]
  let mut ti = 0
  while ti < test_seqs.length() {
    let s = test_seqs[ti]
    let gc = @src.decipher_gc_content(s)
    let complexity = @src.decipher_complexity(s)
    println("  " + s + ": GC=" + (gc * 100.0).to_string() + "%, complexity=" + complexity.to_string())
    ti = ti + 1
  }
  println("")

  // Step 10: Browse Sequences
  println("--- Step 10: Browse Sequences ---")
  let display = @src.browse_seqs(
    ["ACGTACGTACGT", "TTGCGCGCATAT", "GCGCATATATAT"],
    names=["E.coli", "L.plantarum", "B.longum"],
  )
  println(display)

  println("=== DECIPHER Demo Complete ===")
}