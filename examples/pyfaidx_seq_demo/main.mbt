///|
/// pyfaidx sequence-utility demo.
///
/// Demonstrates IUPAC complement / reverse / reverse-complement, GC content,
/// FASTA line wrapping, and BED / UCSC region-string parsing.

///|
fn main raise {
  let seq = "ATCGTA"
  println("complement          : \{@src.pyfaidx_complement(seq)}")
  println("reverse             : \{@src.pyfaidx_reverse(seq)}")
  println("reverse complement  : \{@src.pyfaidx_reverse_complement(seq)}")

  let amb = "NMRATCGTA"
  println("gc        \{amb} = \{@src.pyfaidx_gc(amb)}")
  println("gc_strict \{amb} = \{@src.pyfaidx_gc_strict(amb)}")
  println("gc_iupac  \{amb} = \{@src.pyfaidx_gc_iupac(amb)}")

  println("\nwrapped (n=3):")
  println(@src.pyfaidx_wrap_sequence(3, "ACGTACGTAC"))

  let (r1, s1, e1) = @src.pyfaidx_bed_split("chr1\t10\t20")
  println("bed_split  : \{r1} \{s1} \{e1}")

  let (r2, s2, e2) = @src.pyfaidx_ucsc_split("chr1:100-200")
  println("ucsc_split : \{r2} \{s2} \{e2}")
}