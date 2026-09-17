///|
/// scikit-bio alpha diversity demo.
///
/// Computes a spread of scikit-bio alpha diversity metrics on a small taxon
/// count vector.

///|
fn main {
  let c = [1.0, 1.0, 2.0, 2.0, 5.0, 8.0]
  println("counts          : \{c}")
  println("observed species: \{@src.skbio_alpha_sobs(c)}")
  println("singletons      : \{@src.skbio_alpha_singles(c)}")
  println("doubletons      : \{@src.skbio_alpha_doubles(c)}")
  println("dominance (D)   : \{@src.skbio_alpha_dominance(c)}")
  println("inv_simpson     : \{@src.skbio_alpha_inv_simpson(c)}")
  println("berger-parker d : \{@src.skbio_alpha_berger_parker_d(c)}")
  println("margalef        : \{@src.skbio_alpha_margalef(c)}")
  println("mcintosh D      : \{@src.skbio_alpha_mcintosh_d(c)}")
  println("brillouin d     : \{@src.skbio_alpha_brillouin_d(c)}")
  println("strong          : \{@src.skbio_alpha_strong(c)}")
  println("renyi(2)        : \{@src.skbio_alpha_renyi(c, 2.0)}")
  println("tsallis(2)      : \{@src.skbio_alpha_tsallis(c, 2.0)}")
  println("hill(2)         : \{@src.skbio_alpha_hill(c, 2.0)}")
}