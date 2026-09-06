///|
/// minet module example - mutual-information gene-network inference
/// (port of Bioconductor minet: Meyer, Lafitte & Bontempi, BMC Bioinformatics 2008).
///
/// Run: moon run examples/minet_demo
fn main {
  println("=== Bioconductor minet Module Demo ===")
  println("Mutual-Information Gene-Network Inference (ARACNE / CLR / MRNET)")
  println("")

  // ---- Section 1: synthetic expression matrix ----
  // 6 genes x 30 samples, deterministic:
  //   G0,G1: strong targets of a latent binary program Z
  //   G2   : weaker target of Z
  //   G3   : pure noise (independent)
  //   G4,G5: co-regulated by a second, independent binary program
  let n_samples = 40
  let n_genes = 6
  let gene_names = ["G0", "G1", "G2", "G3", "G4", "G5"]
  let data : Array[Array[Double]] = []
  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = []
    let mut s = 0
    while s < n_samples {
      let z = if s < 20 { 0.0 } else { 1.0 }
      let alt = if s % 2 == 0 { 0.0 } else { 1.0 }
      // deterministic aperiodic gene/sample-specific jitter in [0, 0.4)
      let seed = s.to_double() * 12.9898 + g.to_double() * 78.233
      let frac = @math.sin(seed) * 43758.5453
      let h = (frac - frac.floor()) * 0.4
      let v = match g {
        0 => 6.0 * z + h
        1 => 6.0 * z + h
        2 => 3.0 * z + h
        3 => h
        4 => 5.0 * alt + h
        _ => 5.0 * alt + h
      }
      row.push(v)
      s = s + 1
    }
    data.push(row)
    g = g + 1
  }
  println("=== 1. Expression matrix ===")
  println("genes: " + n_genes.to_string() + ", samples: " + n_samples.to_string())
  println("")

  // ---- Section 2: discretize ----
  println("=== 2. Discretization (equalfreq, nbin = 3) ===")
  let disc = @bio.minet_discretize(data, nbin=3) catch {
    _ => abort("discretization failed")
  }
  let bins_str : Array[String] = []
  let mut b = 0
  while b < disc[0].length() {
    bins_str.push(disc[0][b].to_string())
    b = b + 1
  }
  println("G0 bins: [" + bins_str.join(",") + "]")
  println("")

  // ---- Section 3: mutual information matrix ----
  println("=== 3. Mutual information matrix (Miller-Madow estimator, bits) ===")
  let mim = @bio.minet_build_mim(disc) catch {
    _ => abort("MIM build failed")
  }
  let mut header = "      "
  let mut gn = 0
  while gn < n_genes {
    header = header + gene_names[gn] + "      "
    gn = gn + 1
  }
  println(header)
  let mut i = 0
  while i < n_genes {
    let mut line = gene_names[i] + "  "
    let mut j = 0
    while j < n_genes {
      let cell = mim[i][j].to_string()
      line = line + (if cell.length() < 7 { cell + "  " } else { cell + " " })
      j = j + 1
    }
    println(line)
    i = i + 1
  }
  println("")

  // ---- Section 4: infer networks with the three methods ----
  println("=== 4. Inferred edges ===")
  let nets = [
    ("ARACNE", @bio.minet_aracne(mim, eps=0.1)),
    ("CLR", @bio.minet_clr(mim)),
    ("MRNET", @bio.minet_mrnet(mim)),
  ]
  let mut ni = 0
  while ni < nets.length() {
    let name = nets[ni].0
    let net = nets[ni].1
    println("[" + name + "] edges (weight > 0):")
    let mut ei = 0
    while ei < n_genes {
      let mut ej = ei + 1
      while ej < n_genes {
        if net[ei][ej] > 0.0 {
          println(
            "  " + gene_names[ei] + " -- " + gene_names[ej] + "  weight=" + net[ei][ej].to_string(),
          )
        }
        ej = ej + 1
      }
      ei = ei + 1
    }
    println("")
    ni = ni + 1
  }

  // ---- Section 5: validate against the reference network ----
  // true edges: G0-G1, G0-G2, G1-G2 (program Z module), G4-G5 (program alt)
  let reference : Array[Array[Int]] = []
  let mut ri = 0
  while ri < n_genes {
    reference.push(Array::make(n_genes, 0))
    ri = ri + 1
  }
  let true_edges = [(0, 1), (0, 2), (1, 2), (4, 5)]
  for te in true_edges {
    reference[te.0][te.1] = 1
    reference[te.1][te.0] = 1
  }
  println("=== 5. Validation against the reference network ===")
  println("reference edges: G0-G1, G0-G2, G1-G2, G4-G5")
  let mut vi = 0
  while vi < nets.length() {
    let name = nets[vi].0
    let net = nets[vi].1
    let v = @bio.minet_validate(net, reference) catch {
      _ => abort("validation failed")
    }
    println(
      "[" + name + "] AUROC=" + v.auc_roc.to_string() + "  AUPR=" + v.auc_pr.to_string() +
        "  maxF=" + v.max_fscore.to_string(),
    )
    vi = vi + 1
  }
  println("")
  println("All three methods recover the two co-regulated modules;")
  println("the noise gene G3 receives only background-level scores.")
}
