///|
fn main {
  println("=== motifStack Demo ===")

  let motifs = @src.motifstack_sample_motifs()
  println("\n1. Motif collection")
  for m in motifs {
    println("   " + m.name + " (width=" + m.len().to_string() + ")")
  }

  println("\n2. Pairwise distances (1 - Pearson)")
  let dm = @src.motifstack_distance_matrix(motifs)
  for i in 0..<dm.length() {
    for j in (i + 1)..<dm.length() {
      println(
        "   " + motifs[i].name + " vs " + motifs[j].name +
        ": " + dm[i][j].to_string(),
      )
    }
  }

  println("\n3. Hierarchical clustering")
  match @src.motifstack_hclust(dm) {
    Some(root) => {
      let leaves = root.leaves()
      println("   leaf order: " + leaves.to_string())
      println("   merge height: " + root.height().to_string())
    }
    None => println("   no clustering")
  }

  println("\n4. Stacking order (by clustering)")
  let order = @src.motifstack_stack_order(motifs)
  println("   order: " + order.to_string())

  println("\n5. Per-column information content")
  for m in motifs {
    let ics = @src.motifstack_column_ic(m)
    let mut line = m.name + ": "
    for ic in ics {
      line = line + ic.to_string() + " "
    }
    println("   " + line)
  }

  println("\n6. ASCII logo for GATA1")
  println(@src.motifstack_ascii_logo(motifs[0], rows=5))

  println("\n7. Stacked logos (clustering order)")
  println(@src.motifstack_stack_ascii(motifs, rows=4))

  println("\n8. Flattened PWM vector")
  let flat = @src.motifstack_flatten(motifs[0])
  println("   " + motifs[0].name + ": " + flat.length().to_string() + " values")

  println("\n9. Node traversal")
  let leaf = @src.MotifStackNode::Leaf(0)
  println("   leaf height: " + leaf.height().to_string())
  println("   leaf leaves: " + leaf.leaves().to_string())
  let merge = @src.MotifStackNode::Merge(
    @src.MotifStackNode::Leaf(0),
    @src.MotifStackNode::Leaf(1),
    0.5,
  )
  println("   merge height: " + merge.height().to_string())
  println("   merge leaves: " + merge.leaves().to_string())
}