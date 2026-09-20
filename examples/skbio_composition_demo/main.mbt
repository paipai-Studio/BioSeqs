///|
/// scikit-bio composition (CoDA) demo: closure, clr/ilr/alr transforms and the
/// Aitchison arithmetic operations.

///|
fn main {
  let x = [0.1, 0.3, 0.4, 0.2]

  println("=== closure([2, 2, 6]) ===")
  println(@src.skbio_closure([2.0, 2.0, 6.0]))

  println("\n=== clr([0.1, 0.3, 0.4, 0.2]) ===")
  println(@src.skbio_clr(x))

  println("\n=== ilr (default Gram-Schmidt basis) ===")
  let basis = @src.skbio_gram_schmidt_basis(4)
  let ilr = @src.skbio_ilr(x, basis)
  println(ilr)
  println("ilr_inv(ilr(x)) round-trips to closure(x):")
  println(@src.skbio_ilr_inv(ilr, basis))

  println("\n=== alr (reference = component 0) ===")
  println(@src.skbio_alr([0.1, 0.3, 0.6], 0))

  println("\n=== Aitchison arithmetic ===")
  println("perturb([1/3,1/3,1/3], [1/2,1/2,1]) = \{@src.skbio_perturb([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], [0.5, 0.5, 1.0])}")
  println("inner([.1,.3,.4,.2], [.2,.4,.2,.2])  = \{@src.skbio_inner(x, [0.2, 0.4, 0.2, 0.2])}")

  println("\n=== multi_replace (zero imputation) ===")
  println(@src.skbio_multi_replace([0.2, 0.4, 0.4, 0.0]))
}