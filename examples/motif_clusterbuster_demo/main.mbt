///|
/// motif_clusterbuster demo - ClusterBuster clustered motif format.
/// Inspired by Biopython Bio.motifs (clusterbuster support).

fn main {
  println("=== motif_clusterbuster Demo: Clustered Motifs ===")
  println("")

  // Step 1: Parse ClusterBuster format
  println("--- Step 1: Parse ClusterBuster Format ---")
  let content = ">cluster1 2.5\nmotifA\n[0.7 0.1 0.1 0.1]\n[0.1 0.1 0.7 0.1]\n>cluster2 1.8\nmotifB\n[0.3 0.2 0.2 0.3]\n[0.2 0.3 0.1 0.4]\n"
  let clusters = @src.clusterbuster_parse(content)
  println("Clusters: " + clusters.length().to_string())
  let mut c = 0
  while c < clusters.length() {
    let cl = clusters[c]
    println("  " + cl.name + " (score=" + cl.score.to_string() + ") has " + cl.motifs.length().to_string() + " motifs")
    c = c + 1
  }
  println("")

  // Step 2: Inspect motifs
  println("--- Step 2: Motif Details ---")
  let mut c2 = 0
  while c2 < clusters.length() {
    let cl = clusters[c2]
    let mut m = 0
    while m < cl.motifs.length() {
      let motif = cl.motifs[m]
      println("  " + motif.name + " (" + motif.pfm.length().to_string() + " positions)")
      let mut p = 0
      while p < motif.pfm.length() {
        println("    pos " + p.to_string() + ": [" + motif.pfm[p][0].to_string() + " " + motif.pfm[p][1].to_string() + " " + motif.pfm[p][2].to_string() + " " + motif.pfm[p][3].to_string() + "]")
        p = p + 1
      }
      m = m + 1
    }
    c2 = c2 + 1
  }
  println("")

  // Step 3: Cluster score
  println("--- Step 3: Cluster Score ---")
  let mut c3 = 0
  while c3 < clusters.length() {
    let cl = clusters[c3]
    let score = @src.clusterbuster_cluster_score(cl)
    println("  " + cl.name + " cluster score: " + score.to_string())
    c3 = c3 + 1
  }
  println("")

  // Step 4: Round-trip
  println("--- Step 4: Round-trip Serialization ---")
  let output = @src.clusterbuster_to_string(clusters)
  println("Serialized (first 100 chars):")
  println("  " + output[0:100].to_string() + "...")
  let reparsed = @src.clusterbuster_parse(output)
  println("Re-parsed clusters: " + reparsed.length().to_string())
  println("")

  // Step 5: Convert to MotifPFM
  println("--- Step 5: Convert to MotifPFM ---")
  let pfm = @src.clusterbuster_to_motif_pfm(clusters[0].motifs[0])
  println(@src.motif_pfm_to_string(pfm))
  println("")

  println("=== motif_clusterbuster Demo Complete ===")
}