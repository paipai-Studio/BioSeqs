///|
/// Fasta.get_spliced_seq demo: splice multiple 1-based inclusive intervals.

///|
fn main raise {
  let fasta = ">chr1\nACGTACGTACGTACGT\n>chr2\nAAACCCGGGTTT\n"
  let fa = @src.Fasta::from_content(fasta)

  let fwd = fa.get_spliced_seq("chr1", [(1, 5), (9, 13)], false)
  let rev = fa.get_spliced_seq("chr1", [(1, 5), (9, 13)], true)
  println("chr1 [1,5]+[9,13]  fwd = \{fwd.unwrap().to_string()}")
  println("chr1 [1,5]+[9,13]  rc  = \{rev.unwrap().to_string()}")

  let fwd2 = fa.get_spliced_seq("chr2", [(1, 4), (5, 7), (10, 13)], false)
  println("chr2 [1,4]+[5,7]+[10,13] fwd = \{fwd2.unwrap().to_string()}")
}