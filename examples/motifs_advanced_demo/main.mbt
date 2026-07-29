///|
/// Advanced Motifs Demo - JASPAR, TRANSFAC, alignment, and clustering.
///
/// Run: moon run examples/motifs_advanced_demo/main.mbt
fn main {
  println("=== Advanced Motifs Analysis Demo ===\n")

  // 1. JASPAR PFM parsing
  println("1. JASPAR PFM parsing:")
  let jaspar_pfm =
    ">MA0001.1 E2F1\n" +
    "A  [ 45   0 100   0   0   0   0   0  10   0  10  35  90 ]\n" +
    "C  [  1   0   0   0   0   0   0  90  80  10   0  10   0 ]\n" +
    "G  [  0  50   0   0 100 100  10   0  10  80  90  45  10 ]\n" +
    "T  [ 50   0   0 100   0   0  90   0   0  10   0  10   0 ]\n"
  let jaspar_motifs = @bio.parse_jaspar_pfm(jaspar_pfm)
  println("   Number of motifs: \{jaspar_motifs.length()}")
  if jaspar_motifs.length() > 0 {
    let jaspar_motif = jaspar_motifs[0]
    println("   Motif name: \{jaspar_motif.name}")
    println("   Motif width: \{jaspar_motif.pwm.length()}")
  }

  // 2. TRANSFAC parsing
  println("\n2. TRANSFAC format parsing:")
  let transfac =
    "ID  MA0002.1\n" +
    "NA  E2F1\n" +
    "PO    A   C   G   T\n" +
    "01  45   1   0  50\n" +
    "02   0   0  50   0\n" +
    "03 100   0   0   0\n" +
    "04   0   0   0 100\n" +
    "05   0   0 100   0\n" +
    "06   0   0 100   0\n" +
    "07   0   0  10  90\n" +
    "08   0  90   0   0\n" +
    "09  10  80  10   0\n" +
    "10   0  10  80  10\n" +
    "11  10   0  90   0\n" +
    "12  35  10  45  10\n" +
    "13  90   0  10   0\n" +
    "XX\n" +
    "//\n"
  let transfac_motifs = @bio.parse_transfac(transfac)
  println("   Number of motifs: \{transfac_motifs.length()}")
  if transfac_motifs.length() > 0 {
    let transfac_motif = transfac_motifs[0]
    println("   Motif name: \{transfac_motif.name}")
    println("   Motif width: \{transfac_motif.pwm.length()}")
  }

  // 3. Motif alignment
  println("\n3. Optimal motif alignment:")
  let pwm1 = [
    [0.9, 0.1, 0.0, 0.0],
    [0.1, 0.0, 0.0, 0.9],
    [0.0, 0.9, 0.1, 0.0],
    [0.8, 0.1, 0.1, 0.0],
  ]
  let pwm2 = [
    [0.85, 0.1, 0.05, 0.0],
    [0.1, 0.0, 0.0, 0.9],
    [0.0, 0.85, 0.15, 0.0],
  ]
  let alignment = @bio.optimal_motif_alignment(pwm1, pwm2)
  println("   Best score: \{alignment.score}")
  println("   Offset: \{alignment.offset}")
  println("   Strand: \{alignment.strand}")
  println("   Aligned length: \{alignment.aligned_length}")

  // 4. Motif comparison
  println("\n4. Motif comparison (KL divergence):")
  let kl_div = @bio.motif_kl_divergence(pwm1, pwm2)
  println("   KL divergence: \{kl_div}")

  println("\n5. Motif comparison (JS divergence):")
  let js_div = @bio.motif_js_divergence(pwm1, pwm2)
  println("   JS divergence: \{js_div}")

  // 6. Motif clustering
  println("\n6. Motif clustering:")
  let motif_names : Array[String] = ["motif1", "motif2"]
  let motif_pwms = [pwm1, pwm2]
  let clusters = @bio.cluster_motifs(motif_names, motif_pwms, 0.5)
  println("   Number of clusters: \{clusters.length()}")
  for i = 0; i < clusters.length(); i = i + 1 {
    let cluster = clusters[i]
    println("   Cluster \{cluster.cluster_id}: \{cluster.motifs.length()} motifs")
    println("     Representative: \{cluster.representative}")
    println("     Avg similarity: \{cluster.avg_similarity}")
  }

  println("\n=== Demo Complete ===")
}
