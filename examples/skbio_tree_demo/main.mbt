///|
/// scikit-bio tree + UniFrac demo: parse a rooted Newick tree, list its tips,
/// and compute unweighted/weighted UniFrac between two samples.

///|
fn main {
  let tree = @src.skbio_tree_read("(((a:1,b:2):0.5,(c:1,d:1):0.5):1.0)root;")

  println("=== tree tips ===")
  println(tree.tips())

  println("\n=== UniFrac between two samples (scikit-bio doctest tree) ===")
  let t = @src.skbio_tree_read(
    "(((((U1:0.5,U2:0.5):0.5,U3:1.0):1.0):0.0," +
    "(U4:0.75,(U5:0.5,((U6:0.33,U7:0.62):0.5,U8:0.5):0.5):0.5):1.25):0.0)root;",
  )
  let taxa = ["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8"]
  let u = [1.0, 0.0, 0.0, 4.0, 1.0, 2.0, 3.0, 0.0]
  let v = [0.0, 1.0, 1.0, 6.0, 0.0, 1.0, 0.0, 0.0]
  println("unweighted UniFrac = \{@src.skbio_unweighted_unifrac(u, v, taxa, t)}")
  println("weighted   UniFrac = \{@src.skbio_weighted_unifrac(u, v, taxa, t)}")
}