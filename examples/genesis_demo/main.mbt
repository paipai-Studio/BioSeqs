///|
/// GENESIS demo.
///
/// Demonstrates population structure and kinship analysis.

pub fn main() {
  println("=== GENESIS Demo ===")
  println("")
  
  let genotypes = [
    [0.0, 1.0, 2.0, 1.0, 0.0],
    [0.0, 1.0, 2.0, 1.0, 0.0],
    [1.0, 2.0, 0.0, 0.0, 1.0],
    [1.0, 2.0, 0.0, 0.0, 1.0],
    [2.0, 0.0, 1.0, 2.0, 2.0],
    [2.0, 0.0, 1.0, 2.0, 2.0],
  ]
  
  println("Input genotypes: 6 samples, 5 SNPs")
  println("")
  
  let kinship = bio_genesis_estimate_kinship(genotypes)
  
  println("Kinship matrix:")
  let mut i = 0
  while i < kinship.kinship_matrix.length() {
    let mut row_str = ""
    let mut j = 0
    while j < kinship.kinship_matrix[i].length() {
      row_str = row_str.concat("{:.4f} ".format(kinship.kinship_matrix[i][j]))
      j = j + 1
    }
    println("  {}", row_str)
    i = i + 1
  }
  println("")
  
  let pca = bio_genesis_pca(genotypes, 3)
  
  println("PCA result:")
  println("  Eigenvalues: {:?}", pca.eigenvalues)
  println("  Variance explained: {:?}", pca.var_explained)
  println("")
  
  println("PCA coordinates:")
  let mut k = 0
  while k < pca.eigenvectors.length() {
    let mut coord_str = ""
    let mut l = 0
    while l < pca.eigenvectors[k].length() {
      coord_str = coord_str.concat("{:.4f} ".format(pca.eigenvectors[k][l]))
      l = l + 1
    }
    println("  {}: {}", pca.sample_ids[k], coord_str)
    k = k + 1
  }
  println("")
  
  let rel = bio_genesis_relatedness(genotypes[0], genotypes[1])
  println("Relatedness between sample 0 and 1: {:.4f}", rel)
  println("")
  
  let distance = bio_genesis_genetic_distance(genotypes)
  
  println("Genetic distance matrix:")
  let mut m = 0
  while m < distance.length() {
    let mut row_str = ""
    let mut n = 0
    while n < distance[m].length() {
      row_str = row_str.concat("{:.4f} ".format(distance[m][n]))
      n = n + 1
    }
    println("  {}", row_str)
    m = m + 1
  }
  
  println("")
  println("=== Demo Complete ===")
}