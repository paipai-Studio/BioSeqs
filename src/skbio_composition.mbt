///|
/// scikit-bio `skbio.stats.composition` port — compositional data analysis
/// (CoDA) in the Aitchison simplex. A composition is a vector whose components
/// sum to 1; the single-composition functions below take/return an
/// `Array[Double]` (one row), and `skbio_centralize` works on a matrix of rows.
///
/// Semantics cross-checked against scikit-bio 0.7.3 (`stats/composition/_base.py`).

///|
fn comp_ln(x : Double) -> Double {
  @math.ln(x)
}

///|
fn comp_exp(x : Double) -> Double {
  @math.exp(x)
}

///|
fn comp_sum(x : Array[Double]) -> Double {
  let mut s = 0.0
  for v in x {
    s = s + v
  }
  s
}

///|
fn comp_max(x : Array[Double]) -> Double {
  let mut m = x[0]
  for v in x {
    if v > m {
      m = v
    }
  }
  m
}

///|
/// Close a composition to sum to 1: x_i / Σx (`skbio.stats.composition.closure`).
pub fn skbio_closure(x : Array[Double]) -> Array[Double] {
  let s = comp_sum(x)
  let out : Array[Double] = Array(capacity=0)
  for v in x {
    out.push(v / s)
  }
  out
}

///|
/// Replace zeros with small positive values, keeping compositions summing to
/// 1 (`multi_replace`, formerly `multiplicative_replacement`). Zero entries
/// become `delta`, non-zero entries are scaled by `(1 - z·delta)` where `z` is
/// the number of zeros. `delta` defaults to `(1/D)²` (D = number of components).
fn multi_replace_impl(x : Array[Double], delta : Double) -> Array[Double] {
  let closed = skbio_closure(x)
  let mut z = 0
  for v in closed {
    if v == 0.0 {
      z = z + 1
    }
  }
  let zcnt = 1.0 - z.to_double() * delta
  let out : Array[Double] = Array(capacity=0)
  for v in closed {
    if v == 0.0 {
      out.push(delta)
    } else {
      out.push(zcnt * v)
    }
  }
  out
}

///|
/// Replace zeros using the default `delta = (1/D)²`.
pub fn skbio_multi_replace(x : Array[Double]) -> Array[Double] {
  let d = 1.0 / (x.length().to_double() * x.length().to_double())
  multi_replace_impl(x, d)
}

///|
/// Replace zeros using an explicit `delta`.
pub fn skbio_multi_replace_delta(x : Array[Double], delta : Double) -> Array[Double] {
  multi_replace_impl(x, delta)
}

///|
/// Aitchison perturbation: x ⊕ y = closure(x·y) (`perturb`).
pub fn skbio_perturb(x : Array[Double], y : Array[Double]) -> Array[Double] {
  let cx = skbio_closure(x)
  let cy = skbio_closure(y)
  let prod : Array[Double] = Array(capacity=0)
  for i = 0; i < x.length(); i = i + 1 {
    prod.push(cx[i] * cy[i])
  }
  skbio_closure(prod)
}

///|
/// Aitchison inverse perturbation: x ⊖ y = closure(x / y) (`perturb_inv`).
pub fn skbio_perturb_inv(x : Array[Double], y : Array[Double]) -> Array[Double] {
  let cx = skbio_closure(x)
  let cy = skbio_closure(y)
  let q : Array[Double] = Array(capacity=0)
  for i = 0; i < x.length(); i = i + 1 {
    q.push(cx[i] / cy[i])
  }
  skbio_closure(q)
}

///|
/// Aitchison power: x ⊙ a = closure(x^a) (`power`).
pub fn skbio_power(x : Array[Double], a : Double) -> Array[Double] {
  let cx = skbio_closure(x)
  let p : Array[Double] = Array(capacity=0)
  for v in cx {
    p.push(comp_exp(a * comp_ln(v)))
  }
  skbio_closure(p)
}

///|
/// Centered log-ratio transform: clr(x) = ln(x) - mean(ln(x)) (`clr`).
pub fn skbio_clr(x : Array[Double]) -> Array[Double] {
  let logs : Array[Double] = Array(capacity=0)
  let mut s = 0.0
  for v in x {
    let l = comp_ln(v)
    logs.push(l)
    s = s + l
  }
  let mean = s / x.length().to_double()
  let out : Array[Double] = Array(capacity=0)
  for l in logs {
    out.push(l - mean)
  }
  out
}

///|
/// Inverse clr: closure(exp(x - max(x))) (`clr_inv`, numerically stable).
pub fn skbio_clr_inv(x : Array[Double]) -> Array[Double] {
  let mx = comp_max(x)
  let e : Array[Double] = Array(capacity=0)
  for v in x {
    e.push(comp_exp(v - mx))
  }
  skbio_closure(e)
}

///|
/// Aitchison inner product: ⟨x, y⟩ = Σ_d clr(x)_d · clr(y)_d (`inner`).
pub fn skbio_inner(x : Array[Double], y : Array[Double]) -> Double {
  let cx = skbio_clr(x)
  let cy = skbio_clr(y)
  let mut s = 0.0
  for i = 0; i < cx.length(); i = i + 1 {
    s = s + cx[i] * cy[i]
  }
  s
}

///|
/// Additive log-ratio transform: alr(x) = ln(x_i / x_ref) for i ≠ ref
/// (`alr` with `ref_idx`, default 0 = first component is the denominator).
pub fn skbio_alr(x : Array[Double], ref_idx : Int) -> Array[Double] {
  let ref_l = comp_ln(x[ref_idx])
  let out : Array[Double] = Array(capacity=0)
  for i = 0; i < x.length(); i = i + 1 {
    if i != ref_idx {
      out.push(comp_ln(x[i]) - ref_l)
    }
  }
  out
}

///|
/// Inverse alr: insert 0 at `ref_idx`, then closure(exp(v - max(v))) (`alr_inv`).
pub fn skbio_alr_inv(x : Array[Double], ref_idx : Int) -> Array[Double] {
  let comp : Array[Double] = Array(capacity=0)
  for i = 0; i < ref_idx; i = i + 1 {
    comp.push(x[i])
  }
  comp.push(0.0)
  for i = ref_idx; i < x.length(); i = i + 1 {
    comp.push(x[i])
  }
  skbio_clr_inv(comp)
}

///|
/// Default orthonormal ILR basis (shape `(n-1) × n`), the J. J. Egozcue basis
/// built by Gram–Schmidt orthogonalization (scikit-bio `_gram_schmidt_basis`).
pub fn skbio_gram_schmidt_basis(n : Int) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < n - 1; j = j + 1 {
    let i = j + 1
    let id = i.to_double()
    let s = (id / (id + 1.0)).sqrt()
    let row : Array[Double] = Array(capacity=0)
    for d = 0; d < n; d = d + 1 {
      let v = if d <= j {
        1.0 / id * s
      } else if d == i {
        0.0 - s
      } else {
        0.0
      }
      row.push(v)
    }
    out.push(row)
  }
  out
}

///|
/// Isometric log-ratio transform: ilr(x) = clr(x) · basisᵀ. `basis` is
/// `(n-1) × n`; use `skbio_gram_schmidt_basis(n)` for the default basis.
pub fn skbio_ilr(x : Array[Double], basis : Array[Array[Double]]) -> Array[Double] {
  let c = skbio_clr(x)
  let out : Array[Double] = Array(capacity=0)
  for j = 0; j < basis.length(); j = j + 1 {
    let row = basis[j]
    let mut s = 0.0
    for d = 0; d < c.length(); d = d + 1 {
      s = s + c[d] * row[d]
    }
    out.push(s)
  }
  out
}

///|
/// Inverse isometric log-ratio: ilr_inv(x) = clr_inv(x · basis). `x` is
/// `(n-1)`-dimensional and `basis` is `(n-1) × n`.
pub fn skbio_ilr_inv(x : Array[Double], basis : Array[Array[Double]]) -> Array[Double] {
  let n = basis[0].length()
  let comp : Array[Double] = Array(capacity=0)
  for d = 0; d < n; d = d + 1 {
    let mut s = 0.0
    for j = 0; j < x.length(); j = j + 1 {
      s = s + x[j] * basis[j][d]
    }
    comp.push(s)
  }
  skbio_clr_inv(comp)
}

///|
fn comp_col_gmean(mat : Array[Array[Double]], col : Int) -> Double {
  let mut s = 0.0
  let mut n = 0.0
  for row in mat {
    s = s + comp_ln(row[col])
    n = n + 1.0
  }
  comp_exp(s / n)
}

///|
/// Center each composition around the per-component geometric mean across rows:
/// perturb_inv(closure(matᵢ), gmean(closure(mat), axis=0)) (`centralize`).
/// `mat` is a list of compositions (rows); returns the centered matrix.
pub fn skbio_centralize(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let closed : Array[Array[Double]] = Array(capacity=0)
  for row in mat {
    closed.push(skbio_closure(row))
  }
  let d = mat[0].length()
  let cen : Array[Double] = Array(capacity=0)
  for col = 0; col < d; col = col + 1 {
    cen.push(comp_col_gmean(closed, col))
  }
  let out : Array[Array[Double]] = Array(capacity=0)
  for row in closed {
    out.push(skbio_perturb_inv(row, cen))
  }
  out
}