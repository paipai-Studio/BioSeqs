///|
/// genome3d.mbt — 3D 基因组建模(接触图 -> 经典 MDS -> 3D 坐标 -> PDB).
///
/// 对齐 ShRec3D / 3DMax 的核心建模语义(本地实现):
/// - 接触图: cooler/hicdc 产出的 bin x bin 接触矩阵(对称)
/// - 距离转换: d_ij = 1/contact_ij(ShRec3D 倒数变换, 零接触取 d_max)
/// - 经典 MDS: 距离平方双中心化 B = -0.5·J·D²·J, Jacobi 特征分解
///   取前 3 大特征值/向量得 3D 坐标 X = V·Λ^0.5
/// - 质量评估: 恢复欧氏距离与目标距离的 Pearson 相关(应力指标)
/// - 输出: PDB ATOM 记录(CA 链, 每 bin 一 bead), 与结构栈(tm_score/
///   cealign/dssp)闭环
///
/// 数据不出本机.
pub suberror Genome3dError {
  Genome3dError(String)
}

///|
/// 3D 基因组模型.
pub(all) struct Genome3dModel {
  n_beads : Int32
  /// 每 bead 的 3D 坐标 [x, y, z]
  coordinates : Array[Array[Double]]
  /// 前 3 特征值(MDS 各维方差)
  eigenvalues : Array[Double]
  stress_correlation : Double
  target_distances : Array[Array[Double]]
}

///|
/// 接触矩阵 -> 目标距离(ShRec3D 倒数变换).
/// 零接触(含对角线语义)取 d_max(默认 1.0); 对角距离 0.
pub fn genome3d_contact_to_distance(
  contacts : Array[Array[Double]],
  d_max~ : Double = 1.0,
) -> Array[Array[Double]] raise Genome3dError {
  let n = contacts.length()
  if n == 0 {
    raise Genome3dError("empty contact map")
  }
  for row in contacts {
    if row.length() != n {
      raise Genome3dError(
        "contact map must be square (got \{row.length()} != \{n})",
      )
    }
  }
  let d : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row[j] = 0.0
      } else if contacts[i][j] > 0.0 {
        row[j] = 1.0 / contacts[i][j]
      } else {
        row[j] = d_max
      }
    }
    d.push(row)
  }
  d
}

///|
/// 经典 MDS: 目标距离矩阵 -> 3D 坐标.
pub fn genome3d_cmds(
  distances : Array[Array[Double]],
  dims~ : Int32 = 3,
) -> (Array[Array[Double]], Array[Double]) raise Genome3dError {
  let n = distances.length()
  if n == 0 {
    raise Genome3dError("empty distance matrix")
  }
  for row in distances {
    if row.length() != n {
      raise Genome3dError("distance matrix must be square")
    }
  }
  // D² 双中心化: B = -0.5 J D² J, J = I - (1/n)11^T
  let d2 : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      row.push(distances[i][j] * distances[i][j])
    }
    d2.push(row)
  }
  // 行/列均值
  let row_mean : Array[Double] = Array::make(n, 0.0)
  let col_mean : Array[Double] = Array::make(n, 0.0)
  let mut grand = 0.0
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n; j = j + 1 {
      s = s + d2[i][j]
    }
    row_mean[i] = s / n.to_double()
    grand = grand + s
  }
  grand = grand / (n.to_double() * n.to_double())
  for j = 0; j < n; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < n; i = i + 1 {
      s = s + d2[i][j]
    }
    col_mean[j] = s / n.to_double()
  }
  let b : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      row.push(-0.5 * (d2[i][j] - row_mean[i] - col_mean[j] + grand))
    }
    b.push(row)
  }
  // Jacobi 特征分解
  let max_iter = 200
  let (evals_raw, evecs_raw) = linalg_jacobi(b, n, max_iter)
  let (evals, evecs) = linalg_sort_eigen(evals_raw, evecs_raw, n)
  // 取前 dims 维(特征值为正才有效)
  let k = if dims > 0 && dims < n { dims } else { 3 }
  let coords : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let p : Array[Double] = Array::make(k, 0.0)
    for c = 0; c < k; c = c + 1 {
      if evals[c] > 0.0 {
        // evecs[i][c] 为第 c 特征向量在 i 的分量
        p[c] = evals[c].sqrt() * evecs[i][c]
      }
    }
    coords.push(p)
  }
  let top : Array[Double] = Array(capacity=0)
  for c = 0; c < k; c = c + 1 {
    top.push(evals[c])
  }
  (coords, top)
}

///|
/// 恢复坐标的欧氏距离 vs 目标距离 Pearson 相关.
fn genome3d_eval_correlation(
  coordinates : Array[Array[Double]],
  target : Array[Array[Double]],
) -> Double {
  let n = coordinates.length()
  let xs : Array[Double] = Array(capacity=0)
  let ys : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let mut s = 0.0
      for c = 0; c < coordinates[i].length(); c = c + 1 {
        let d = coordinates[i][c] - coordinates[j][c]
        s = s + d * d
      }
      xs.push(s.sqrt())
      ys.push(target[i][j])
    }
  }
  pearson(xs, ys)
}

///|
/// Pearson 相关系数(复用库内实现语义, 内联以避免依赖).
fn pearson(xs : Array[Double], ys : Array[Double]) -> Double {
  let n = xs.length()
  if n == 0 {
    return 0.0
  }
  let mut sx = 0.0
  let mut sy = 0.0
  for i = 0; i < n; i = i + 1 {
    sx = sx + xs[i]
    sy = sy + ys[i]
  }
  let mx = sx / n.to_double()
  let my = sy / n.to_double()
  let mut sxy = 0.0
  let mut sxx = 0.0
  let mut syy = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = xs[i] - mx
    let dy = ys[i] - my
    sxy = sxy + dx * dy
    sxx = sxx + dx * dx
    syy = syy + dy * dy
  }
  let denom = (sxx * syy).sqrt()
  if denom <= 0.0 {
    return 0.0
  }
  sxy / denom
}

///|
/// 全流程建模: 接触矩阵 -> 3D 模型.
pub fn genome3d_model(
  contacts : Array[Array[Double]],
  d_max~ : Double = 1.0,
) -> Genome3dModel raise Genome3dError {
  let target = genome3d_contact_to_distance(contacts, d_max~)
  let (coords, evals) = genome3d_cmds(target)
  let corr = genome3d_eval_correlation(coords, target)
  {
    n_beads: coords.length(),
    coordinates: coords,
    eigenvalues: evals,
    stress_correlation: corr,
    target_distances: target,
  }
}

///|
/// 模型摘要.
pub fn Genome3dModel::summary(self : Genome3dModel) -> String {
  let ev : Array[String] = []
  for e in self.eigenvalues {
    ev.push(e.to_string())
  }
  "genome3d: \{self.n_beads} beads, top eigenvalues [" +
  ev.join(", ") +
  "], distance correlation = \{self.stress_correlation}"
}

///|
/// bead 间恢复距离(模型坐标).
pub fn Genome3dModel::model_distance(
  self : Genome3dModel,
  i : Int32,
  j : Int32,
) -> Double raise Genome3dError {
  if i < 0 || i >= self.n_beads || j < 0 || j >= self.n_beads {
    raise Genome3dError("bead index out of range")
  }
  let mut s = 0.0
  for c = 0; c < self.coordinates[i].length(); c = c + 1 {
    let d = self.coordinates[i][c] - self.coordinates[j][c]
    s = s + d * d
  }
  s.sqrt()
}

///|
/// PDB ATOM 记录输出(CA 链, 每 bin 一 bead, 残基序 = bin+1).
pub fn Genome3dModel::to_pdb(
  self : Genome3dModel,
  chain~ : Char = 'A',
) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "REMARK   1 BIOSEQS GENOME3D CONTACT MAP TO 3D VIA CLASSICAL MDS\n",
  )
  sb.write_string(
    "REMARK   2 DISTANCE CORRELATION \{self.stress_correlation}\n",
  )
  for i = 0; i < self.n_beads; i = i + 1 {
    let p = self.coordinates[i]
    let x = if p.length() > 0 { p[0] } else { 0.0 }
    let y = if p.length() > 1 { p[1] } else { 0.0 }
    let z = if p.length() > 2 { p[2] } else { 0.0 }
    // ATOM 行(列位对齐 PDB 固定格式)
    let serial = (i + 1) % 100000
    let resseq = (i + 1) % 10000
    sb.write_string(
      "ATOM  \{pdb_pad(serial.to_string(), 5)}  CA  BIN \{chain}\{pdb_pad(resseq.to_string(), 5)}   ",
    )
    sb.write_string(pdb_coord(x))
    sb.write_string(pdb_coord(y))
    sb.write_string(pdb_coord(z))
    sb.write_string("  1.00  0.00           C  \n")
  }
  sb.write_string("TER\nEND\n")
  sb.to_string()
}

///|
fn pdb_pad(s : String, width : Int32) -> String {
  let mut out = s
  while out.length() < width {
    out = " " + out
  }
  out
}

///|
fn pdb_coord(x : Double) -> String {
  // %8.3f 定点语义(避免科学计数法)
  let neg = x < 0.0
  let a = if neg { -x } else { x }
  // 千分位舍入
  let scaled = (a * 1000.0).round()
  let int_part = (scaled / 1000.0).trunc().to_int()
  let frac_val = (scaled - int_part.to_double() * 1000.0).to_int()
  let mut s = int_part.to_string() + "."
  let f = frac_val.to_string()
  let mut pad = f.length()
  while pad < 3 {
    s = s + "0"
    pad = pad + 1
  }
  s = s + f
  if neg {
    s = "-" + s
  }
  while s.length() < 8 {
    s = " " + s
  }
  if s.length() > 8 {
    s = s[s.length() - 8:].to_owned()
  }
  s
}

///|
/// 读取 PDB 模型的 CA 坐标(与结构栈 parse 闭环的轻量入口).
pub fn genome3d_parse_pdb_ca(
  content : String,
) -> Array[Array[Double]] raise Genome3dError {
  let out : Array[Array[Double]] = Array(capacity=0)
  let lines = content.split("\n").to_array()
  for l in lines {
    let line = l.to_owned()
    if line.has_prefix("ATOM") && line.length() >= 54 {
      let atom_name = line[12:16].to_owned().trim()
      if atom_name == "CA" {
        let xs = line[30:38].to_owned().trim()
        let ys = line[38:46].to_owned().trim()
        let zs = line[46:54].to_owned().trim()
        match
          (
            parse_double(xs.to_owned()),
            parse_double(ys.to_owned()),
            parse_double(zs.to_owned()),
          ) {
          (Some(x), Some(y), Some(z)) => out.push([x, y, z])
          _ => ()
        }
      }
    }
  }
  if out.length() == 0 {
    raise Genome3dError("no CA atoms found")
  }
  out
}
