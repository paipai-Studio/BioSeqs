///|
/// Align Applications Demo - Demonstrates multiple sequence alignment tool wrappers.

fn main {
  println("=== Align Applications Demo ===")
  
  let clustalw = @bio.create_example_clustalw()
  println("ClustalW command:")
  println(clustalw.build())
  println("")
  
  let clustalo = @bio.create_example_clustalomega()
  println("Clustal Omega command:")
  println(clustalo.build())
  println("")
  
  let muscle = @bio.create_example_muscle()
  println("Muscle command:")
  println(muscle.build())
  println("")
  
  let mafft = @bio.create_example_mafft()
  println("MAFFT command:")
  println(mafft.build())
}