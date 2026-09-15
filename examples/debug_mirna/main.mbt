fn main {
  // Simple test with DNA
  let mirna = "ATGCATGC"
  let seed = @src.mir_extract_seed(mirna)
  println("mirna: " + mirna + " seed: " + seed)
  
  // seed = "TGCATGC", rc = "GCATGCA"
  // UTR with "GCATGCA" preceded by A -> 8mer
  let utr = "AAGCATGCA"
  let sites = @src.mir_find_seed_matches(mirna, utr, min_seed_type=2)
  println("utr: " + utr + " sites: " + sites.length().to_string())
  for s in sites {
    println("  type=" + @src.mir_seed_type_label(s.seed_type))
  }
}
