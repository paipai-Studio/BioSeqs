///|
/// scikit-allel `allel.stats.sf` — the (joint) site frequency spectrum and its
/// folded / scaled variants.

///|
fn sfs_bincount(values : Array[Int], minlength : Int) -> Array[Int] {
  let out = Array::make(minlength, 0)
  for v in values {
    if v >= 0 && v < minlength {
      out[v] = out[v] + 1
    }
  }
  out
}

///|
fn sfs_max(values : Array[Int]) -> Int {
  let mut m = 0
  for v in values {
    if v > m {
      m = v
    }
  }
  m
}

///|
/// Site frequency spectrum from derived allele counts (`allel.sfs`).
pub fn skallel_sfs(dac : Array[Int], n : Int?) -> Array[Int] {
  let nn = match n {
    Some(x) => x
    None => sfs_max(dac)
  }
  sfs_bincount(dac, nn + 1)
}

///|
/// Folded site frequency spectrum from biallelic allele counts (`allel.sfs_folded`).
pub fn skallel_sfs_folded(ac : Array[Array[Int]], n : Int?) -> Array[Int] {
  let mac : Array[Int] = Array(capacity=0)
  let mut mx = 0
  for aci in ac {
    let m = if aci[0] < aci[1] { aci[0] } else { aci[1] }
    mac.push(m)
    if aci[0] + aci[1] > mx {
      mx = aci[0] + aci[1]
    }
  }
  let nn = match n {
    Some(x) => x
    None => mx
  }
  sfs_bincount(mac, nn / 2 + 1)
}

///|
/// Scale an unfolded site frequency spectrum by the derived-allele count
/// (`allel.scale_sfs`).
pub fn skallel_scale_sfs(s : Array[Int]) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for k = 0; k < s.length(); k = k + 1 {
    out.push(s[k].to_double() * k.to_double())
  }
  out
}

///|
/// Scaled site frequency spectrum (`allel.sfs_scaled`).
pub fn skallel_sfs_scaled(dac : Array[Int], n : Int?) -> Array[Double] {
  let s = skallel_sfs(dac, n)
  skallel_scale_sfs(s)
}

///|
/// Scale a folded site frequency spectrum (`allel.scale_sfs_folded`).
pub fn skallel_scale_sfs_folded(s : Array[Int], n : Int) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for k = 0; k < s.length(); k = k + 1 {
    out.push(
      s[k].to_double() * k.to_double() * (n - k).to_double() / n.to_double(),
    )
  }
  out
}

///|
/// Scaled folded site frequency spectrum (`allel.sfs_folded_scaled`).
pub fn skallel_sfs_folded_scaled(ac : Array[Array[Int]], n : Int?) -> Array[Double] {
  let s = skallel_sfs_folded(ac, n)
  let mx = match n {
    Some(x) => x
    None => {
      let mut m = 0
      for aci in ac {
        if aci[0] + aci[1] > m {
          m = aci[0] + aci[1]
        }
      }
      m
    }
  }
  skallel_scale_sfs_folded(s, mx)
}

///|
/// Joint site frequency spectrum between two populations (`allel.joint_sfs`).
pub fn skallel_joint_sfs(
  dac1 : Array[Int],
  dac2 : Array[Int],
  n1 : Int?,
  n2 : Int?,
) -> Array[Array[Int]] {
  let nn1 = match n1 {
    Some(x) => x
    None => sfs_max(dac1)
  }
  let nn2 = match n2 {
    Some(x) => x
    None => sfs_max(dac2)
  }
  let x = nn1 + 1
  let y = nn2 + 1
  let tmp : Array[Int] = Array(capacity=0)
  let nt = if dac1.length() < dac2.length() { dac1.length() } else { dac2.length() }
  for i = 0; i < nt; i = i + 1 {
    tmp.push(dac1[i] * y + dac2[i])
  }
  let flat = sfs_bincount(tmp, x * y)
  let out : Array[Array[Int]] = Array(capacity=0)
  for i = 0; i < x; i = i + 1 {
    let row : Array[Int] = Array(capacity=0)
    for j = 0; j < y; j = j + 1 {
      row.push(flat[i * y + j])
    }
    out.push(row)
  }
  out
}

///|
/// Folded joint site frequency spectrum (`allel.joint_sfs_folded`).
pub fn skallel_joint_sfs_folded(
  ac1 : Array[Array[Int]],
  ac2 : Array[Array[Int]],
  n1 : Int?,
  n2 : Int?,
) -> Array[Array[Int]] {
  let mac1 : Array[Int] = Array(capacity=0)
  let mac2 : Array[Int] = Array(capacity=0)
  let nt = if ac1.length() < ac2.length() { ac1.length() } else { ac2.length() }
  for i = 0; i < nt; i = i + 1 {
    mac1.push(if ac1[i][0] < ac1[i][1] { ac1[i][0] } else { ac1[i][1] })
    mac2.push(if ac2[i][0] < ac2[i][1] { ac2[i][0] } else { ac2[i][1] })
  }
  let sum1 : Array[Int] = Array(capacity=0)
  let sum2 : Array[Int] = Array(capacity=0)
  for i = 0; i < nt; i = i + 1 {
    sum1.push(ac1[i][0] + ac1[i][1])
    sum2.push(ac2[i][0] + ac2[i][1])
  }
  let mx1 = match n1 { Some(x) => x; None => sfs_max(sum1) }
  let mx2 = match n2 { Some(x) => x; None => sfs_max(sum2) }
  let x = mx1 / 2 + 1
  let y = mx2 / 2 + 1
  let tmp : Array[Int] = Array(capacity=0)
  for i = 0; i < nt; i = i + 1 {
    tmp.push(mac1[i] * y + mac2[i])
  }
  let flat = sfs_bincount(tmp, x * y)
  let out : Array[Array[Int]] = Array(capacity=0)
  for i = 0; i < x; i = i + 1 {
    let row : Array[Int] = Array(capacity=0)
    for j = 0; j < y; j = j + 1 {
      row.push(flat[i * y + j])
    }
    out.push(row)
  }
  out
}

///|
/// Fold a site frequency spectrum (`allel.fold_sfs`): `o[k] = s[k] + s[n−k]`.
pub fn skallel_fold_sfs(s : Array[Int], n : Int) -> Array[Int] {
  let padded : Array[Int] = Array(capacity=0)
  for k = 0; k < n + 1; k = k + 1 {
    padded.push(if k < s.length() { s[k] } else { 0 })
  }
  let nf = (n + 1) / 2
  let nn = nf * 2
  let out : Array[Int] = Array(capacity=0)
  for k = 0; k < nf; k = k + 1 {
    out.push(padded[k] + padded[nn - 1 - k])
  }
  out
}

///|
/// Fold a joint site frequency spectrum (`allel.fold_joint_sfs`).
pub fn skallel_fold_joint_sfs(
  s : Array[Array[Int]],
  n1 : Int,
  n2 : Int,
) -> Array[Array[Int]] {
  let mf = (n1 + 1) / 2
  let nf = (n2 + 1) / 2
  let r1 = mf * 2
  let r2 = nf * 2
  let out : Array[Array[Int]] = Array(capacity=0)
  for i = 0; i < mf; i = i + 1 {
    let row : Array[Int] = Array(capacity=0)
    for j = 0; j < nf; j = j + 1 {
      let get = fn(r : Int, c : Int) -> Int {
        if r < s.length() && c < s[r].length() {
          s[r][c]
        } else {
          0
        }
      }
      let v = get(i, j) + get(r1 - 1 - i, j) + get(i, r2 - 1 - j) + get(r1 - 1 - i, r2 - 1 - j)
      row.push(v)
    }
    out.push(row)
  }
  out
}