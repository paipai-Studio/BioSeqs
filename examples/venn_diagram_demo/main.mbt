///|
/// VennDiagram module example - Venn diagram set operations and region calculations.
///
/// Run: moon run examples/venn_diagram_demo
fn main {
  println("=== VennDiagram Module Demo ===")
  println("Venn Diagram Set Operations and Region Calculations")
  println("")

  // Load sample data
  println("=== Sample Data ===")
  let sets = @bio.venn_sample_gene_sets()
  let names = @bio.venn_sample_gene_names()
  let n_sets = sets.length()
  println("Number of sets: \{n_sets}")
  let mut i = 0
  while i < n_sets {
    println("  \{names[i]}: \{sets[i].length()} elements")
    i = i + 1
  }
  println("")

  // Create Venn diagram
  println("=== Venn Diagram ===")
  let result = @bio.venn_diagram(sets, names=names)
  println("Total elements in universe: \{result.total_elements}")
  println("Number of regions: \{result.regions.length()}")
  println("")

  // Show only-in-set elements
  println("=== Elements Only in Each Set ===")
  i = 0
  while i < n_sets {
    let only = @bio.venn_only_in_set(result, i)
    println("  Only in \{names[i]}: \{only.length()} elements")
    if only.length() > 0 {
      println("    First 5: \{preview(only, 5)}")
    }
    i = i + 1
  }
  println("")

  // Show intersection of all sets
  println("=== Intersection of All Sets ===")
  let all_intersect = @bio.venn_all_intersect(result)
  println("Elements in all \{n_sets} sets: \{all_intersect.length()}")
  if all_intersect.length() > 0 {
    println("  Elements: \{preview(all_intersect, 10)}")
  }
  println("")

  // Show unique to each set
  println("=== Unique to Each Set ===")
  let unique = @bio.venn_unique_to_each(result)
  i = 0
  while i < unique.length() {
    println("  Unique to \{names[i]}: \{unique[i].length()} elements")
    i = i + 1
  }
  println("")

  // Pairwise overlap
  println("=== Pairwise Overlap ===")
  i = 0
  while i < n_sets {
    let mut j = i + 1
    while j < n_sets {
      if j < n_sets {
        let overlap = @bio.venn_pairwise_overlap(result, i, j)
        println("  \{names[i]} vs \{names[j]}:")
        println("    Intersection: \{overlap["intersection_count"]}")
        println("    Union: \{overlap["union_count"]}")
        println("    Jaccard: \{format_stat(overlap["jaccard_index"])}")
        println("    Overlap coefficient: \{format_stat(overlap["overlap_coefficient"])}")
        println("    Dice coefficient: \{format_stat(overlap["dice_coefficient"])}")
      }
      j = j + 1
    }
    i = i + 1
  }
  println("")

  // Set operations
  println("=== Set Operations ===")
  let a = sets[0]
  let b = sets[1]
  println("Set A (\{names[0]}): \{a.length()} elements")
  println("Set B (\{names[1]}): \{b.length()} elements")

  let inter = @bio.venn_intersection(a, b)
  println("Intersection (A ∩ B): \{inter.length()} elements")
  if inter.length() > 0 {
    println("  Elements: \{preview(inter, 5)}")
  }

  let uni = @bio.venn_union_two(a, b)
  println("Union (A ∪ B): \{uni.length()} elements")

  let diff = @bio.venn_difference(a, b)
  println("Difference (A - B): \{diff.length()} elements")
  if diff.length() > 0 {
    println("  Elements: \{preview(diff, 5)}")
  }
  println("")

  // Similarity metrics
  println("=== Similarity Metrics ===")
  let jacc = @bio.venn_jaccard(a, b)
  println("Jaccard index (A vs B): \{format_stat(jacc)}")
  let dice = @bio.venn_dice(a, b)
  println("Dice coefficient (A vs B): \{format_stat(dice)}")
  println("")

  // Summary
  println("=== Summary ===")
  let summary = @bio.venn_summary(result)
  println(summary)
}

fn preview(arr : Array[String], n : Int) -> String {
  let mut s = "["
  let limit = if arr.length() < n { arr.length() } else { n }
  let mut i = 0
  while i < limit {
    if i > 0 {
      s = s + ", "
    }
    s = s + arr[i]
    i = i + 1
  }
  if arr.length() > n {
    s = s + ", ..."
  }
  s = s + "]"
  s
}

fn format_stat(value : Double) -> String {
  let rounded = (value * 1000.0).round() / 1000.0
  "\{rounded}"
}
