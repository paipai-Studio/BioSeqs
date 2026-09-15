fn main {
  // No Kozak features
  let seq = "TTTTTTATGT"
  let score = @src.orff_kozak_score(seq, 6)
  println("Kozak score for TTTTTTATGT at pos 6: " + score.to_string())
  // pos -3 = T (not A or G) -> 0
  // pos +4 = T (not G) -> 0
  // pos -2 = T (not C or T) -> 0.15
  // pos -1 = T (not C or T) -> 0.15
  // Expected: 0.3
}
