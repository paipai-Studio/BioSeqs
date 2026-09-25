///|
/// Bio.Cluster Demo.
/// Demonstrates the numerical clustering core ported from Biopython's
/// Bio.Cluster (the C Clustering Library): distance metrics, hierarchical
/// clustering, k-means/k-medoids, self-organizing maps, PCA, and the
/// Eisen Cluster/TreeView Record file format.

/// Format a Double to a String with a given number of decimal places.
fn fmt(value : Double, decimals : Int) -> String {
  let s = value.to_string()
  let parts = s.split(".").to_array()
  let int_part = if parts.length() > 0 { parts[0].to_owned() } else { "0" }
  let frac_part = if parts.length() > 1 { parts[1].to_owned() } else { "" }
  let mut padded_frac = frac_part
  if padded_frac.length() < decimals {
    padded_frac = padded_frac + "0".repeat(decimals - padded_frac.length())
  } else if padded_frac.length() > decimals {
    padded_frac = padded_frac[0:decimals].to_owned()
  }
  int_part + "." + padded_frac
}

fn print_vector(name : String, v : Array[Double]) -> Unit {
  let parts : Array[String] = Array(capacity=0)
  for x in v {
    parts.push(fmt(x, 4))
  }
  println("\{name}: [\{parts.join(", ")}]")
}

fn print_matrix(name : String, m : Array[Array[Double]]) -> Unit {
  println("\{name}:")
  for row in m {
    let parts : Array[String] = Array(capacity=0)
    for x in row {
      parts.push(fmt(x, 4))
    }
    println("  [\{parts.join(", ")}]")
  }
}

fn print_int_array(name : String, v : Array[Int]) -> Unit {
  let parts : Array[String] = Array(capacity=0)
  for x in v {
    parts.push("\{x}")
  }
  println("\{name}: [\{parts.join(", ")}]")
}

///|
fn main raise {
  // Sample gene expression data: 4 genes x 3 experiments.
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0],
    [2.0, 4.0, 6.0],
    [7.0, 3.0, 0.0],
    [4.0, 2.0, 6.0],
  ]

  println("=== 1. Distance metrics (clusterdistance) ===")
  // Note: like the C Clustering Library, the Euclidean distance is the
  // weighted mean of SQUARED differences (no square root).
  let d_euclid = @src.clusterdistance(data, index1=[0], index2=[1], dist="e")
  let d_city = @src.clusterdistance(data, index1=[0], index2=[1], dist="b")
  let d_corr = @src.clusterdistance(data, index1=[0], index2=[1], dist="c")
  let d_spear = @src.clusterdistance(data, index1=[0], index2=[1], dist="s")
  println("euclid(g0, g1)     = \{fmt(d_euclid, 4)}")
  println("cityblock(g0, g1)  = \{fmt(d_city, 4)}")
  println("correlation(g0,g1) = \{fmt(d_corr, 4)}")
  println("spearman(g0, g1)   = \{fmt(d_spear, 4)}")

  println("\n=== 2. Hierarchical clustering (treecluster) ===")
  let tree = @src.treecluster(data, method_="a", dist="e")
  println("number of nodes = \{tree.length()}")
  for i = 0; i < tree.length(); i = i + 1 {
    let n = tree.node(i)
    println(
      "node \{i}: left=\{n.left} right=\{n.right} distance=\{fmt(n.distance, 4)}",
    )
  }
  print_int_array("cut into 2 clusters", tree.cut(nclusters=2))

  println("\n=== 3. Distance matrix + hierarchical clustering from it ===")
  let dm = @src.distancematrix(data)
  println("lower-triangular distance matrix (\{dm.length()} rows):")
  for i = 0; i < dm.length(); i = i + 1 {
    let parts : Array[String] = Array(capacity=0)
    for j = 0; j < dm[i].length(); j = j + 1 {
      parts.push(fmt(dm[i][j], 4))
    }
    println("  [\{parts.join(", ")}]")
  }
  let tree2 = @src.treecluster_from_distance(dm, method_="s")
  print_int_array("single-linkage cut into 2 clusters", tree2.cut(nclusters=2))

  println("\n=== 4. k-means / k-medoids clustering ===")
  @src.cluster_seed(12345)
  let kc = @src.kcluster(data, nclusters=2, npass=10, dist="e")
  print_int_array("kcluster assignments", kc.clusterid())
  println("kcluster within-cluster sum of distances = \{fmt(kc.error(), 4)}")
  println("kcluster found \{kc.ifound} time(s)")
  let kmed = @src.kmedoids(dm, nclusters=2, npass=10)
  print_int_array("kmedoids assignments (medioids)", kmed.clusterid())
  println("kmedoids total distance = \{fmt(kmed.error(), 4)}")

  println("\n=== 5. Cluster centroids ===")
  let cent = @src.clustercentroids(data, clusterid=[0, 0, 1, 1], method_="a")
  print_matrix("cluster means", cent.cdata())

  println("\n=== 6. Self-organizing map (somcluster) ===")
  @src.cluster_seed(12345)
  let som = @src.somcluster(data, nxgrid=2, nygrid=1, inittau=0.02, niter=100)
  for i = 0; i < som.clusterid().length(); i = i + 1 {
    let (ix, iy) = som.clusterid()[i]
    println("gene \{i} -> cell (\{ix}, \{iy})")
  }

  println("\n=== 7. Principal component analysis (pca) ===")
  let pca = @src.cluster_pca(data)
  print_vector("column means", pca.columnmean())
  print_vector("eigenvalues", pca.eigenvalues())
  print_matrix("coordinates", pca.coordinates())
  print_matrix("principal components", pca.pc())

  println("\n=== 8. Record: Eisen Cluster/TreeView file format ===")
  let content = "UID\tNAME\texp1\texp2\texp3\n" +
    "g1\talpha\t1.0\t2.0\t3.0\n" +
    "g2\tbeta\t2.0\t4.0\t6.0\n" +
    "g3\tgamma\t7.0\t3.0\t0.0\n" +
    "g4\tdelta\t4.0\t2.0\t6.0\n" +
    "EWEIGHT\t\t1.0\t1.0\t1.0\n"
  let rec = @src.bio_cluster_read(content)
  println("uniqid = \{rec.uniqid()}")
  println("geneid = [\{rec.geneid().join(", ")}]")
  let names = match rec.genename() {
    Some(ns) => ns.join(", ")
    None => "(none)"
  }
  println("genename = [\{names}]")
  let rt = rec.treecluster(method_="s")
  print_int_array("Record.treecluster(method=\"s\") cut into 2",
    rt.cut(nclusters=2))
  let rk = rec.kcluster(nclusters=2, initialid=[0, 0, 1, 1])
  print_int_array("Record.kcluster(nclusters=2)", rk.clusterid())
  let rcd = rec.clusterdistance(index1=[0], index2=[3], dist="e")
  println("clusterdistance(g0, g3) = \{fmt(rcd, 4)}")
  let rdm = rec.distancematrix()
  println(
    "distancematrix row 3 = [\{fmt(rdm[3][0], 4)}, \{fmt(rdm[3][1], 4)}, \{fmt(rdm[3][2], 4)}]",
  )

  println("\nDone.")
}
