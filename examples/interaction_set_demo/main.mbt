///|
/// InteractionSet Demo - Chromatin interaction data analysis.
fn main {
  println("=== InteractionSet Demo ===\n")

  println("1. Creating anchors...")
  let a1 = @src.Anchor::new(chrom="chr1", start=10000, end=20000)
  let a2 = @src.Anchor::new(chrom="chr1", start=50000, end=60000)
  let a3 = @src.Anchor::new(chrom="chr2", start=10000, end=20000)
  println(
    "   Anchor1: \{a1.chrom}:\{a1.start}-\{a1.end} (width=\{@src.Anchor::width(a1)})",
  )
  println("   Anchor2: \{a2.chrom}:\{a2.start}-\{a2.end}")

  println("\n2. Checking anchor overlaps...")
  let ov1 = @src.anchors_overlap(a1, a2)
  let ov2 = @src.anchors_overlap(a1, a3)
  println("   Anchor1 vs Anchor2 overlap: \{ov1}")
  println("   Anchor1 vs Anchor3 overlap: \{ov2}")

  println("\n3. Creating genomic interactions...")
  let i1 = @src.GInteraction::new(anchor1=a1, anchor2=a2, count=25)
  let i2 = @src.GInteraction::new(anchor1=a1, anchor2=a3, count=10)
  let i3 = @src.GInteraction::new(anchor1=a2, anchor2=a3, count=5)
  let interactions = [i1, i2, i3]
  let iset = @src.InteractionSet::new(interactions~)
  println("   \{@src.InteractionSet::summary(iset)}")

  println("\n4. Computing interaction distances...")
  match @src.GInteraction::distance(i1) {
    Some(d) => println("   Interaction1 distance: \{d} bp")
    None => println("   No distance (inter-chromosomal)")
  }
  match @src.GInteraction::distance(i2) {
    Some(d) => println("   Interaction2 distance: \{d} bp")
    None => println("   Interaction2: No distance (inter-chromosomal)")
  }

  println("\n5. Filtering by count (min_count=8)...")
  let filtered = @src.InteractionSet::filter_by_count(iset, 8)
  println(
    "   \{@src.InteractionSet::num_interactions(filtered)} interactions remaining",
  )

  println("\n6. Intra vs inter-chromosomal interactions...")
  let intra = @src.InteractionSet::intra_chromosomal(iset)
  let inter = @src.InteractionSet::inter_chromosomal(iset)
  println(
    "   Intra-chromosomal: \{@src.InteractionSet::num_interactions(intra)}",
  )
  println(
    "   Inter-chromosomal: \{@src.InteractionSet::num_interactions(inter)}",
  )
  println("   Inter fraction: \{@src.InteractionSet::inter_fraction(iset)}")

  println("\n7. Interactions on specific chromosome...")
  let on_chr1 = @src.InteractionSet::interactions_on_chrom(iset, "chr1")
  let on_chr2 = @src.InteractionSet::interactions_on_chrom(iset, "chr2")
  println("   On chr1: \{@src.InteractionSet::num_interactions(on_chr1)}")
  println("   On chr2: \{@src.InteractionSet::num_interactions(on_chr2)}")

  println("\n8. Building interaction matrix...")
  let matrix = @src.InteractionSet::interaction_matrix(iset, ["chr1", "chr2"])
  println("   chr1 x chr1: \{matrix[0][0]}")
  println("   chr1 x chr2: \{matrix[0][1]}")
  println("   chr2 x chr1: \{matrix[1][0]}")
  println("   chr2 x chr2: \{matrix[1][1]}")

  println("\n9. Top interactions...")
  let top = @src.InteractionSet::top_interactions(iset, 2)
  let mut i = 0
  while i < top.length() {
    let t = top[i]
    println(
      "   #\{i + 1}: count=\{t.count}, \{t.anchor1.chrom}:\{t.anchor1.start} <-> \{t.anchor2.chrom}:\{t.anchor2.start}",
    )
    i = i + 1
  }

  println("\n10. Creating binned interactions from Hi-C matrix...")
  let bins = [(0, 1000), (1000, 2000), (2000, 3000), (3000, 4000)]
  let counts = [
    [50, 30, 10, 5],
    [30, 80, 25, 8],
    [10, 25, 60, 15],
    [5, 8, 15, 40],
  ]
  let binned = @src.create_binned_interactions("chr1", bins, counts)
  println("   \{@src.InteractionSet::num_interactions(binned)} interactions")
  println("   Total count: \{@src.InteractionSet::total_count(binned)}")

  println("\n11. Distance distribution...")
  let dist = @src.InteractionSet::distance_distribution(binned, 500)
  let mut j = 0
  while j < dist.length() {
    println("   \{dist[j].0}-\{dist[j].0 + 500} bp: \{dist[j].1} interactions")
    j = j + 1
  }

  println("\n=== InteractionSet Demo Complete ===")
}
