///|
/// ropls demo: PCA, PLS(-DA), OPLS-DA, VIP, cross-validation statistics
/// and response permutation testing, following the Bioconductor ropls API.
///
/// Run with: moon run examples/ropls_demo

/// Round to 3 decimals for compact printed output.
fn r3(v : Double) -> String {
  ((v * 1000.0).round() / 1000.0).to_string()
}

fn fmt_vec(v : Array[Double]) -> String {
  let s = StringBuilder()
  s.write_char('[')
  for i in 0..<v.length() {
    if i > 0 {
      s.write_string(", ")
    }
    s.write_string(r3(v[i]))
  }
  s.write_char(']')
  s.to_string()
}

fn fmt_strs(v : Array[String]) -> String {
  let s = StringBuilder()
  s.write_char('[')
  for i in 0..<v.length() {
    if i > 0 {
      s.write_string(", ")
    }
    s.write_string(v[i])
  }
  s.write_char(']')
  s.to_string()
}

/// Deterministic two-group dataset (12 samples, 5 variables). Variables 0-2
/// separate the groups; variables 3-4 are pure noise.
fn make_x() -> Array[Array[Double]] {
  let x : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<12 {
    let shift = if i < 6 { -2.0 } else { 2.0 }
    let row = [
      shift + 0.35 * @math.sin((i + 1).to_double()),
      -shift + 0.35 * @math.cos((2 * i + 1).to_double()),
      shift + 0.20 * (i.to_double() - 5.5),
      0.30 * @math.sin((3 * i).to_double()),
      0.30 * @math.cos((2 * i + 3).to_double()),
    ]
    x.push(row)
  }
  x
}

fn make_y(x : Array[Array[Double]]) -> Array[Array[Double]] {
  // Continuous response mostly driven by variables 0 and 2.
  let y : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<x.length() {
    y.push([1.5 * x[i][0] - 0.8 * x[i][2] + 0.2 * @math.sin(i.to_double())])
  }
  y
}

fn main raise {
  let x = make_x()
  let labels : Array[String] = Array(capacity=0)
  for i in 0..<x.length() {
    labels.push(if i < 6 { "ctrl" } else { "treat" })
  }
  let y = make_y(x)

  println("== PCA (2 components, unit-variance scaling) ==")
  let pca = @src.ropls_pca(x, ncomp=2)
  println("R2X per component: " + fmt_vec(pca.r2x))
  println("R2X cumulative:    " + fmt_vec(pca.r2x_cum))
  println("NIPALS iterations: " + pca.iters[0].to_string() + ", " +
    pca.iters[1].to_string())

  println("")
  println("== PLS regression (2 predictive components) ==")
  let pls = @src.ropls_fit(x, y=@src.RoplsY::Continuous(y), ncomp=2)
  println(
    "R2X(cum)=" + fmt_vec(pls.r2x_cum) + " R2Y(cum)=" + fmt_vec(pls.r2y_cum) +
      " Q2(cum)=" + fmt_vec(pls.q2_cum),
  )
  println("RMSEE = " + r3(pls.rmsee))
  println("VIP   = " + fmt_vec(pls.vip))
  println("Significance: " + fmt_strs(pls.significance))

  println("")
  println("== PLS-DA (two classes, 2 predictive components) ==")
  let plsda = @src.ropls_fit(x, y=@src.RoplsY::Classes(labels), ncomp=2)
  println("classes: " + fmt_strs(plsda.classes))
  println(
    "R2Y(cum)=" + fmt_vec(plsda.r2y_cum) + " Q2(cum)=" + fmt_vec(plsda.q2_cum),
  )
  println("fitted classes: " + fmt_strs(plsda.fitted_classes))
  // Class probabilities / dummy responses for two representative rows.
  let probe = [x[0], x[11]]
  let pred_cls = @src.ropls_predict_classes(plsda, probe)
  println("predicted for samples 1 and 12: " + fmt_strs(pred_cls))

  println("")
  println("== OPLS-DA (1 predictive + 1 orthogonal component) ==")
  let oplsda =
    @src.ropls_fit(x, y=@src.RoplsY::Classes(labels), ncomp=1, ortho=1)
  println(
    "R2X(cum)=" + fmt_vec(oplsda.r2x_cum) + " R2Y(cum)=" +
      fmt_vec(oplsda.r2y_cum) + " Q2(cum)=" + fmt_vec(oplsda.q2_cum),
  )
  println("orthogonal VIP: " + fmt_vec(oplsda.ortho_vip))

  println("")
  println("== Response permutation test (10 deterministic permutations) ==")
  let n = x.length()
  let perms : Array[Array[Int]] = Array(capacity=0)
  for k in 0..<10 {
    let p : Array[Int] = Array(capacity=0)
    for i in 0..<n {
      p.push((i + 2 * k + 1) % n)
    }
    perms.push(p)
  }
  let perm =
    @src.ropls_permutation(x, y=@src.RoplsY::Continuous(y), perms=perms, ncomp=2)
  println("observed R2Y=" + r3(perm.obs_r2y) +
    " Q2=" + r3(perm.obs_q2))
  println("p_R2Y=" + r3(perm.p_r2y) +
    " p_Q2=" + r3(perm.p_q2))
}
