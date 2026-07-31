///|
/// Example: reduced_demo
/// Demonstrates the Bio.Align.Reduced module (reduced amino acid alphabets).

fn main {
  println("=== Reduced Amino Acid Alphabets Demo ===")
  println("")
  
  // Example 1: RAD alphabet (6 groups)
  println("=== RAD Alphabet (6 groups) ===")
  let rad = @bio.rad_alphabet()
  let rad_name = rad.get_name()
  let rad_ng = rad.get_n_groups()
  println("Name: \{rad_name}")
  println("Groups: \{rad_ng}")
  println("Groups:")
  for g in rad.get_groups() {
    println("  \{g}")
  }
  println("")
  
  // Reduce a sequence
  let seq1 = "MVLSPADKTNVKAAWGCK"
  let reduced1 = rad.reduce_sequence(seq1)
  println("Original: \{seq1}")
  println("Reduced:  \{reduced1}")
  println("")
  
  // Example 2: Dayhoff alphabet
  println("=== Dayhoff Alphabet (6 groups) ===")
  let dh = @bio.dayhoff_alphabet()
  let dh_name = dh.get_name()
  println("Name: \{dh_name}")
  println("Groups:")
  for g in dh.get_groups() {
    println("  \{g}")
  }
  println("")
  
  let reduced2 = dh.reduce_sequence(seq1)
  println("Original: \{seq1}")
  println("Reduced:  \{reduced2}")
  println("")
  
  // Example 3: CHARM alphabet
  println("=== CHARM Alphabet (4 groups) ===")
  let charm = @bio.charm_alphabet()
  let charm_name = charm.get_name()
  println("Name: \{charm_name}")
  println("Groups:")
  for g in charm.get_groups() {
    println("  \{g}")
  }
  println("")
  
  let reduced3 = charm.reduce_sequence(seq1)
  println("Original: \{seq1}")
  println("Reduced:  \{reduced3}")
  println("")
  
  // Example 4: SDM12 alphabet
  println("=== SDM12 Alphabet (12 groups) ===")
  let sdm = @bio.sdm12_alphabet()
  let sdm_name = sdm.get_name()
  let sdm_ng = sdm.get_n_groups()
  println("Name: \{sdm_name}")
  println("Groups: \{sdm_ng}")
  println("")
  
  let reduced4 = sdm.reduce_sequence(seq1)
  println("Original: \{seq1}")
  println("Reduced:  \{reduced4}")
  println("")
  
  // Example 5: Compare sequences using reduced alphabet
  println("=== Sequence Comparison ===")
  let seq_a = "MVLSPADKTNVKAAWGCK"
  let seq_b = "MVLSPADKTNVKAAWGCK"
  let seq_c = "AAAAAAAAAAAAAAAARR"
  
  println("Comparing sequences with RAD:")
  let id_ab = @bio.reduced_identity(seq_a, seq_b, rad)
  let id_ac = @bio.reduced_identity(seq_a, seq_c, rad)
  println("  A vs B: \{id_ab}")
  println("  A vs C: \{id_ac}")
  
  println("")
  println("Comparing sequences with Dayhoff:")
  let id_ab_dh = @bio.reduced_identity(seq_a, seq_b, dh)
  let id_ac_dh = @bio.reduced_identity(seq_a, seq_c, dh)
  println("  A vs B: \{id_ab_dh}")
  println("  A vs C: \{id_ac_dh}")
  println("")
  
  // Example 6: Get alphabet by name
  println("=== Get Alphabet by Name ===")
  match @bio.get_reduced_alphabet("rad") {
    Some(a) => {
      let a_name = a.get_name()
      let a_ng = a.get_n_groups()
      println("Found: \{a_name} (\{a_ng} groups)")
    }
    None => println("Not found")
  }
  match @bio.get_reduced_alphabet("dayhoff") {
    Some(a) => {
      let a_name = a.get_name()
      let a_ng = a.get_n_groups()
      println("Found: \{a_name} (\{a_ng} groups)")
    }
    None => println("Not found")
  }
  match @bio.get_reduced_alphabet("unknown") {
    Some(_) => println("Found: unexpected")
    None => println("Not found: expected")
  }
  println("")
  
  // Example 7: Summary
  println("=== RAD Summary ===")
  println(rad.summary())
}
