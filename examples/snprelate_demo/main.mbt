///|
/// Bioconductor SNPRelate Demo - SNP genetic relationship analysis.
///
/// Demonstrates the port of the Bioconductor `SNPRelate` package
/// (Zheng & Weir, 2016): per-SNP statistics, KING-robust kinship,
/// KING homo IBD (k0/k1), IBS distances, EIGMIX-style PCA and greedy
/// LD pruning with a sliding window. Rows are SNPs, columns are
/// samples; -1 marks a missing call.

///|
/// 8 SNPs x 6 samples drawn from two ancestry clusters (samples 0-2 vs
/// 3-5) with a missing call at SNP7/Sample4.
fn demo_dosages() -> Array[Array[Int]] {
  [
    [0, 1, 0, 2, 2, 2],
    [1, 1, 0, 2, 1, 2],
    [0, 0, 1, 1, 2, 2],
    [2, 2, 2, 0, 0, 1],
    [1, 2, 2, 0, 0, 0],
    [2, 1, 2, 0, 1, 0],
    [0, 1, 2, 2, 1, 1],
    [1, 1, 1, 0, -1, 2],
  ]
}

///|
fn main raise {
  let dosages = demo_dosages()

  println("=== snpgdsSnpStats ===")
  let stats = @src.snpgds_snp_stats(dosages)
  for k in 0..<stats.allele_frequency.length() {
    println(
      "SNP\{k}: alt freq=\{stats.allele_frequency[k]}, maf=\{stats.maf[k]}, missing=\{stats.missing_rate[k]}",
    )
  }

  println("\n=== snpgdsIBDKing (robust kinship) ===")
  let king = @src.snpgds_ibd_king(dosages)
  let pairs = [(0, 1), (0, 3), (1, 3), (2, 5)]
  for pair in pairs {
    let (i, j) = pair
    println(
      "pair (\{i},\{j}): kinship=\{king.kinship[i][j]}, IBS0=\{king.ibs0[i][j]}",
    )
  }

  println("\n=== snpgdsIBDKing homo (k0 / k1) ===")
  let homo = @src.snpgds_ibd_king_homo(dosages, stats.allele_frequency)
  for pair in pairs {
    let (i, j) = pair
    println(
      "pair (\{i},\{j}): k0=\{homo.k0[i][j]}, k1=\{homo.k1[i][j]}",
    )
  }

  println("\n=== snpgdsIBS / snpgdsDist ===")
  let ibs = @src.snpgds_ibs(dosages)
  for pair in pairs {
    let (i, j) = pair
    println(
      "pair (\{i},\{j}): IBS=\{ibs.ibs[i][j]}, distance=\{ibs.distance[i][j]}",
    )
  }

  println("\n=== snpgdsPCA (EIGMIX) ===")
  let pca = @src.snpgds_pca(dosages, 3)
  for k in 0..<pca.eigenvalues.length() {
    println(
      "PC\{k + 1}: eigenvalue=\{pca.eigenvalues[k]}, var explained=\{pca.var_explained[k]}",
    )
  }
  let pc1 = pca.eigenvectors[0]
  println("PC1 loadings: \{pc1[0]}, \{pc1[1]}, \{pc1[2]}, \{pc1[3]}, \{pc1[4]}, \{pc1[5]}")

  println("\n=== snpgdsLD (corr) ===")
  let ld_pairs = [(0, 1), (0, 3), (1, 3), (3, 4)]
  for pair in ld_pairs {
    let (i, j) = pair
    let r = @src.snpgds_ld(dosages, i, j)
    println("r(SNP\{i}, SNP\{j}) = \{r}")
  }

  println("\n=== snpgdsLDpruning (threshold=0.7, maf=0.05, slide=4) ===")
  let pruned = @src.snpgds_ld_pruning(dosages, 0.7, 0.05, 4)
  println("kept   : \{pruned.kept_indices}")
  println("removed: \{pruned.removed_indices}")

  println("\n=== snpgdsLDpruning (threshold=0.5, maf=0.05, slide=2) ===")
  let pruned2 = @src.snpgds_ld_pruning(dosages, 0.5, 0.05, 2)
  println("kept   : \{pruned2.kept_indices}")
  println("removed: \{pruned2.removed_indices}")
}
