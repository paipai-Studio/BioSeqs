///|
/// Bio.PDB.Polypose：多构象 PDB 结构（NMR 结构集合 / altLoc 交替构象）。
///
/// 提供构象级操作：残基与原子访问、Cα 坐标提取、Kabsch 叠合 pairwise
/// RMSD 矩阵、构象收敛统计、平均结构、altLoc 构象检测与 TSV 导出。
/// 叠合复用 quaternion_superimposer 的四元数刚体叠合器。

/// 3D 坐标点。
pub(all) struct PpyCoord {
  x : Double
  y : Double
  z : Double
}

/// 原子记录。
pub(all) struct PpyAtom {
  atom_name : String
  altloc : String
  element : String
  x : Double
  y : Double
  z : Double
}

/// 残基记录。
pub(all) struct PpyResidue {
  chain_id : String
  resname : String
  resseq : Int
  atoms : Array[PpyAtom]
}

/// 单个构象（一个 model / 一套坐标）。
pub(all) struct PpyConformer {
  conformer_id : String
  residues : Array[PpyResidue]
}

/// 多构象结构集合（NMR bundle）。
pub(all) struct PpyStructure {
  name : String
  conformers : Array[PpyConformer]
}

/// 构象收敛统计。
pub(all) struct PpyConvergence {
  n_conformers : Int
  n_calpha : Int
  mean_rmsd : Double
  min_rmsd : Double
  max_rmsd : Double
  median_rmsd : Double
  max_min_rmsd : Double
  mean_to_reference : Double
}

/// 构造一个原子。
pub fn PpyAtom::new(
  atom_name : String,
  x : Double,
  y : Double,
  z : Double,
  altloc? : String = "",
  element? : String = "C",
) -> PpyAtom {
  PpyAtom::{ atom_name: atom_name, altloc: altloc, element: element, x: x, y: y, z: z }
}

/// 构造一个残基。
pub fn PpyResidue::new(
  resname : String,
  resseq : Int,
  atoms : Array[PpyAtom],
  chain_id? : String = "A",
) -> PpyResidue {
  PpyResidue::{ chain_id: chain_id, resname: resname, resseq: resseq, atoms: atoms }
}

/// 构造一个构象。
pub fn PpyConformer::new(conformer_id : String, residues : Array[PpyResidue]) -> PpyConformer {
  PpyConformer::{ conformer_id: conformer_id, residues: residues }
}

/// 构造一个坐标点。
pub fn PpyCoord::new(x : Double, y : Double, z : Double) -> PpyCoord {
  PpyCoord::{ x: x, y: y, z: z }
}

/// 创建一个空的多构象结构集合。
pub fn ppy_new(name : String) -> PpyStructure {
  PpyStructure::{ name: name, conformers: [] }
}

/// 向集合添加一个构象（保持原有顺序）。
pub fn ppy_add_conformer(st : PpyStructure, conformer : PpyConformer) -> PpyStructure {
  let confs = st.conformers.copy()
  confs.push(conformer)
  PpyStructure::{ name: st.name, conformers: confs }
}

/// 构象数量。
pub fn ppy_n_conformers(st : PpyStructure) -> Int {
  st.conformers.length()
}

/// 按构象 ID 查找构象。
pub fn ppy_get_conformer(st : PpyStructure, conformer_id : String) -> PpyConformer? {
  for c in st.conformers {
    if c.conformer_id == conformer_id {
      return Some(c)
    }
  }
  None
}

/// 取出某一构象的全部残基。
pub fn ppy_conformer_residues(c : PpyConformer) -> Array[PpyResidue] {
  c.residues
}

/// 某一构象中指定残基号（及链）的残基。
pub fn ppy_find_residue(
  c : PpyConformer,
  resseq : Int,
  chain_id? : String = "A",
) -> PpyResidue? {
  for r in c.residues {
    if r.resseq == resseq && r.chain_id == chain_id {
      return Some(r)
    }
  }
  None
}

/// 残基内按原子名取原子（忽略 altLoc：优先取空 altLoc 或第一个匹配）。
pub fn ppy_atom(r : PpyResidue, atom_name : String) -> PpyAtom? {
  let mut fallback : PpyAtom? = None
  for a in r.atoms {
    if a.atom_name == atom_name {
      if a.altloc == "" {
        return Some(a)
      }
      if fallback.is_none() {
        fallback = Some(a)
      }
    }
  }
  fallback
}

/// 残基内按原子名取全部原子（含 altLoc 变体）。
pub fn ppy_atoms(r : PpyResidue, atom_name : String) -> Array[PpyAtom] {
  let out : Array[PpyAtom] = Array(capacity=0)
  for a in r.atoms {
    if a.atom_name == atom_name {
      out.push(a)
    }
  }
  out
}

/// 残基内 altLoc 标识列表（去重）。
pub fn ppy_altloc_ids(r : PpyResidue) -> Array[String] {
  let seen : Map[String, Int] = Map([], capacity=4)
  for a in r.atoms {
    if a.altloc != "" {
      seen.set(a.altloc, 1)
    }
  }
  let out : Array[String] = Array(capacity=0)
  for k in seen.keys() {
    out.push(k)
  }
  out
}

/// 构象中含 altLoc 交替构象的残基数。
pub fn ppy_altloc_residue_count(c : PpyConformer) -> Int {
  let mut n = 0
  for r in c.residues {
    let ids = ppy_altloc_ids(r)
    if ids.length() > 0 {
      n = n + 1
    }
  }
  n
}

/// 按原子名提取构象的坐标（跨残基，顺序为残基顺序）。
/// 跳过找不到该原子的残基。
pub fn ppy_coords(c : PpyConformer, atom_name : String) -> Array[PpyCoord] {
  let out : Array[PpyCoord] = Array(capacity=0)
  for r in c.residues {
    match ppy_atom(r, atom_name) {
      Some(a) => out.push(PpyCoord::{ x: a.x, y: a.y, z: a.z })
      None => continue
    }
  }
  out
}

/// 构象对应残基标签（chain:resseq），与 ppy_coords 顺序对齐。
pub fn ppy_residue_labels(c : PpyConformer, atom_name : String) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for r in c.residues {
    match ppy_atom(r, atom_name) {
      Some(_) => out.push(r.chain_id + ":" + r.resseq.to_string())
      None => continue
    }
  }
  out
}

/// 多构象集合中沿固定次序的坐标数组。
/// 返回失败条件：某构象的坐标数与其他构象不一致（返回 None）。
pub fn ppy_structure_coords(
  st : PpyStructure,
  atom_name : String,
) -> Array[Array[PpyCoord]]? {
  if st.conformers.length() == 0 {
    return Some([])
  }
  let first = ppy_coords(st.conformers[0], atom_name)
  let n = first.length()
  let out : Array[Array[PpyCoord]] = Array(capacity=0)
  out.push(first)
  for i in 1..<st.conformers.length() {
    let coords = ppy_coords(st.conformers[i], atom_name)
    if coords.length() != n {
      return None
    }
    out.push(coords)
  }
  Some(out)
}

/// 将坐标数组转换为叠合器可用的 QsPoint 数组。
fn ppy_to_qspoint(vals : Array[PpyCoord]) -> Array[QsPoint] {
  let out : Array[QsPoint] = Array(capacity=0)
  for p in vals {
    out.push(QsPoint::{ x: p.x, y: p.y, z: p.z })
  }
  out
}

/// 3x3 旋转矩阵作用于向量。
fn ppy_apply_rotation(
  rot : Array[Array[Double]],
  x : Double,
  y : Double,
  z : Double,
) -> (Double, Double, Double) {
  let rx = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z
  let ry = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z
  let rz = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z
  (rx, ry, rz)
}

/// Kabsch / 四元数刚体叠合的 RMSD（复用 quaternion_superimposer）。
/// 要求两坐标数组长度一致且残基一一对应。
pub fn ppy_kabsch_rmsd(a : Array[PpyCoord], b : Array[PpyCoord]) -> Double {
  if a.length() == 0 || a.length() != b.length() {
    return 0.0
  }
  let res = qs_superimpose(ppy_to_qspoint(a), ppy_to_qspoint(b))
  res.rmsd
}

/// 将 b 叠合到 a（返回叠合后的坐标）。
pub fn ppy_superpose(a : Array[PpyCoord], b : Array[PpyCoord]) -> Array[PpyCoord] {
  let out : Array[PpyCoord] = Array(capacity=0)
  if a.length() != b.length() {
    return out
  }
  if a.length() == 0 {
    return out
  }
  let res = qs_superimpose(ppy_to_qspoint(a), ppy_to_qspoint(b))
  let rot = res.rotation
  let t = res.translation
  for p in b {
    let (rx, ry, rz) = ppy_apply_rotation(rot, p.x, p.y, p.z)
    out.push(PpyCoord::{ x: rx + t[0], y: ry + t[1], z: rz + t[2] })
  }
  out
}

/// 构象两两 RMSD 矩阵（Cα），矩阵对称、对角为 0。
/// 若某构象坐标数与首个构象不一致，返回 None。
pub fn ppy_pairwise_rmsds(
  st : PpyStructure,
  atom_name? : String = "CA",
) -> Array[Array[Double]]? {
  let coords_opt = ppy_structure_coords(st, atom_name)
  match coords_opt {
    None => None
    Some(coords) => {
      let n = coords.length()
      let mat : Array[Array[Double]] = Array(capacity=0)
      for i in 0..<n {
        let row : Array[Double] = Array::make(n, 0.0)
        for j in 0..<n {
          if i != j {
            row[j] = ppy_kabsch_rmsd(coords[i], coords[j])
          }
        }
        mat.push(row)
      }
      Some(mat)
    }
  }
}

/// 前 k 个输入中最小的值。
fn ppy_qsort_asc(vals : Array[Double]) -> Array[Double] {
  let sorted = vals.copy()
  let n = sorted.length()
  for i in 0..<n {
    let mut min_idx = i
    let j_start = i + 1
    for j in j_start..<n {
      if sorted[j] < sorted[min_idx] {
        min_idx = j
      }
    }
    if min_idx != i {
      let tmp = sorted[i]
      sorted[i] = sorted[min_idx]
      sorted[min_idx] = tmp
    }
  }
  sorted
}

/// 数组的中位数。
fn ppy_median(vals : Array[Double]) -> Double {
  if vals.length() == 0 {
    return 0.0
  }
  let sorted = ppy_qsort_asc(vals)
  let n = sorted.length()
  if (n % 2) == 1 {
    sorted[(n - 1) / 2]
  } else {
    let a = sorted[n / 2 - 1]
    let b = sorted[n / 2]
    (a + b) / 2.0
  }
}

/// 构象收敛统计。
pub fn ppy_convergence(
  st : PpyStructure,
  atom_name? : String = "CA",
) -> PpyConvergence? {
  let coords_opt = ppy_structure_coords(st, atom_name)
  match coords_opt {
    None => None
    Some(coords) => {
      let n = coords.length()
      if n < 2 {
        return Some(PpyConvergence::{
          n_conformers: n, n_calpha: if n > 0 { coords[0].length() } else { 0 },
          mean_rmsd: 0.0, min_rmsd: 0.0, max_rmsd: 0.0, median_rmsd: 0.0,
          max_min_rmsd: 0.0, mean_to_reference: 0.0,
        })
      }
      let rmsds : Array[Double] = Array(capacity=0)
      let reference_rmsds : Array[Double] = Array(capacity=0)
      let mut max_min = 0.0
      for i in 0..<n {
        let mut min_for_i = 1.0e18
        for j in 0..<n {
          if i == j {
            continue
          }
          let d = ppy_kabsch_rmsd(coords[i], coords[j])
          if i < j {
            rmsds.push(d)
            if i == 0 {
              reference_rmsds.push(d)
            }
          }
          if d < min_for_i {
            min_for_i = d
          }
        }
        if min_for_i > max_min {
          max_min = min_for_i
        }
      }
      let mut total = 0.0
      for d in rmsds {
        total = total + d
      }
      let mut ref_total = 0.0
      for d in reference_rmsds {
        ref_total = ref_total + d
      }
      let n_r = rmsds.length()
      let n_ref = reference_rmsds.length()
      let mut mn = 1.0e18
      let mut mx = 0.0
      for d in rmsds {
        if d < mn {
          mn = d
        }
        if d > mx {
          mx = d
        }
      }
      Some(PpyConvergence::{
        n_conformers: n,
        n_calpha: coords[0].length(),
        mean_rmsd: if n_r > 0 { total / n_r.to_double() } else { 0.0 },
        min_rmsd: mn,
        max_rmsd: mx,
        median_rmsd: ppy_median(rmsds),
        max_min_rmsd: max_min,
        mean_to_reference: if n_ref > 0 { ref_total / n_ref.to_double() } else { 0.0 },
      })
    }
  }
}

/// 平均结构：将所有构象叠合到第 0 个构象后，逐残基平均坐标。
pub fn ppy_mean_structure(
  st : PpyStructure,
  atom_name? : String = "CA",
) -> Array[PpyCoord]? {
  let coords_opt = ppy_structure_coords(st, atom_name)
  match coords_opt {
    None => None
    Some(coords) => {
      let n_conf = coords.length()
      if n_conf == 0 {
        return Some([])
      }
      let n_atoms = coords[0].length()
      let ref_coords = coords[0]
      let sums : Array[(Double, Double, Double)] = Array::make(n_atoms, (0.0, 0.0, 0.0))
      for c in 0..<n_conf {
        let aligned = ppy_superpose(ref_coords, coords[c])
        for i in 0..<n_atoms {
          let (sx, sy2, sz) = sums[i]
          sums[i] = (sx + aligned[i].x, sy2 + aligned[i].y, sz + aligned[i].z)
        }
      }
      let out : Array[PpyCoord] = Array(capacity=0)
      let nf = n_conf.to_double()
      for i in 0..<n_atoms {
        let (sx, sy2, sz) = sums[i]
        out.push(PpyCoord::{ x: sx / nf, y: sy2 / nf, z: sz / nf })
      }
      Some(out)
    }
  }
}

/// RMSD 矩阵导出 TSV。
pub fn ppy_rmsd_tsv(
  matrix : Array[Array[Double]],
  ids : Array[String],
) -> String {
  let sb = StringBuilder()
  sb.write_string("conformer")
  for id in ids {
    sb.write_string("\t" + id)
  }
  sb.write_string("\n")
  for i in 0..<matrix.length() {
    let mut line = ids[i]
    for j in 0..<matrix[i].length() {
      line = line + "\t" + matrix[i][j].to_string()
    }
    sb.write_string(line + "\n")
  }
  sb.to_string()
}

/// 构象内部残基间同名单原子数（原子完整性检查用）。
pub fn ppy_n_atoms(r : PpyResidue) -> Int {
  r.atoms.length()
}

/// 便捷：从行列坐标数组构造结构的 Cα 坐标（测试用）。
pub fn ppy_coords_from_residues(
  residues : Array[PpyResidue],
  atom_name : String,
) -> Array[PpyCoord] {
  let out : Array[PpyCoord] = Array(capacity=0)
  for r in residues {
    match ppy_atom(r, atom_name) {
      Some(a) => out.push(PpyCoord::{ x: a.x, y: a.y, z: a.z })
      None => continue
    }
  }
  out
}