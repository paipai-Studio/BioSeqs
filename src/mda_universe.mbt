// ============================================================================
// MDAnalysis 风格 Universe / AtomGroup 核心类型
// ============================================================================

///|
/// 扩展选择语言用的原子信息结构体。
///
/// `mda_selection_ext.mbt` 里的选择表达式会引用这些字段。
pub struct UnivAtom {
  name : String
  resname : String
  chain : String
  segid : String
  element : String
  atom_type : String
  resid : Int32
  mass : Double
  charge : Double
  is_protein : Bool
  is_water : Bool
  is_backbone : Bool
}

///|
/// 构建 UnivAtom 默认值。
pub fn UnivAtom::new(
  name~ : String = "",
  resname~ : String = "UNK",
  chain~ : String = "",
  segid~ : String = "",
  element~ : String = "",
  atom_type~ : String = "",
  resid~ : Int32 = 1,
  mass~ : Double = 12.0,
  charge~ : Double = 0.0,
  is_protein~ : Bool = false,
  is_water~ : Bool = false,
  is_backbone~ : Bool = false,
) -> UnivAtom {
  UnivAtom::{
    name,
    resname,
    chain,
    segid,
    element,
    atom_type,
    resid,
    mass,
    charge,
    is_protein,
    is_water,
    is_backbone,
  }
}

// ============================================================================
// PBCBox
// ============================================================================

///|
/// 周期性边界条件（PBC）盒向量。
pub struct PBCBox {
  v1 : Vector3
  v2 : Vector3
  v3 : Vector3
  alpha : Double
  beta : Double
  gamma : Double
}

///|
/// 空 box（无 PBC）。
pub fn PBCBox::empty() -> PBCBox {
  PBCBox::{
    v1: Vector3::new(0.0, 0.0, 0.0),
    v2: Vector3::new(0.0, 0.0, 0.0),
    v3: Vector3::new(0.0, 0.0, 0.0),
    alpha: 90.0,
    beta: 90.0,
    gamma: 90.0,
  }
}

///|
/// 从正交盒长度构建。
pub fn PBCBox::from_lengths(a : Double, b : Double, c : Double) -> PBCBox {
  PBCBox::{
    v1: Vector3::new(a, 0.0, 0.0),
    v2: Vector3::new(0.0, b, 0.0),
    v3: Vector3::new(0.0, 0.0, c),
    alpha: 90.0,
    beta: 90.0,
    gamma: 90.0,
  }
}

// ============================================================================
// Universe
// ============================================================================

///|
/// MDAnalysis Universe — 轨迹 + 拓扑的统一容器。
pub struct Universe {
  atoms : Array[UnivAtom]
  trajectory : Array[Array[Vector3]]
  times : Array[Double]
  box : PBCBox
  mut current_frame : Int32
}

///|
/// 构建 Universe。
pub fn Universe::new(
  atoms : Array[UnivAtom],
  trajectory~ : Array[Array[Vector3]] = Array(capacity=0),
  times~ : Array[Double] = Array(capacity=0),
  box~ : PBCBox = PBCBox::empty(),
) -> Universe {
  Universe::{ atoms, trajectory, times, box, current_frame: 0, }
}

///|
/// 从 `Array[Array[Array[Double]]]`（DCD/Trajectory 风格 frames）构建。
///
/// 每个 frame 是 `Array[Array[Double]]`，即 `[atom_idx] -> [x, y, z]`。
pub fn Universe::from_dcd_frames(
  atoms : Array[UnivAtom],
  frames : Array[Array[Array[Double]]],
  box~ : PBCBox = PBCBox::empty(),
) -> Universe {
  let traj : Array[Array[Vector3]] = Array(capacity=frames.length())
  let mut f = 0
  while f < frames.length() {
    let frame = frames[f]
    let v3s : Array[Vector3] = Array(capacity=frame.length())
    let mut a = 0
    while a < frame.length() {
      let c = frame[a]
      if c.length() >= 3 {
        v3s.push(Vector3::new(c[0], c[1], c[2]))
      } else {
        v3s.push(Vector3::new(0.0, 0.0, 0.0))
      }
      a = a + 1
    }
    traj.push(v3s)
    f = f + 1
  }
  Universe::new(atoms, trajectory=traj, box~)
}

///|
/// 原子总数。
pub fn Universe::n_atoms(self : Universe) -> Int32 {
  self.atoms.length()
}

///|
/// 轨迹帧数。
pub fn Universe::n_frames(self : Universe) -> Int32 {
  self.trajectory.length()
}

///|
/// 某一帧的坐标。
pub fn Universe::frame(self : Universe, i : Int32) -> Array[Vector3] {
  if i < 0 || i >= self.trajectory.length() {
    Array(capacity=0)
  } else {
    self.trajectory[i]
  }
}

///|
/// 当前帧索引。
pub fn Universe::frame_index(self : Universe) -> Int32 {
  self.current_frame
}

///|
/// 切换当前帧（钳制到合法范围）。
pub fn Universe::set_frame(self : Universe, i : Int32) -> Unit {
  let n = self.trajectory.length()
  if n <= 0 {
    self.current_frame = 0
    return
  }
  let clamped = if i < 0 { 0 } else if i >= n { n - 1 } else { i }
  self.current_frame = clamped
}

///|
/// 当前帧坐标（返回 `Array[Array[Double]]` 兼容格式）。
pub fn Universe::positions(self : Universe) -> Array[Array[Double]] {
  self.positions_at(self.current_frame)
}

///|
/// 指定帧坐标（返回 `Array[Array[Double]]` 兼容格式）。
pub fn Universe::positions_at(
  self : Universe,
  i : Int32,
) -> Array[Array[Double]] {
  let f = self.frame(i)
  let out : Array[Array[Double]] = Array(capacity=f.length())
  let mut k = 0
  while k < f.length() {
    out.push([f[k].x, f[k].y, f[k].z])
    k = k + 1
  }
  out
}

///|
/// 基于选择表达式筛选原子（复用 mda_selection_ext 的选择器）。
pub fn Universe::select_atoms(self : Universe, expr : String) -> AtomGroup {
  let mask = mda_select_universe(
    expr,
    self.atoms,
    positions=Some(self.positions()),
  )
  AtomGroup::{ indices: mask, universe: self, }
}

///|
/// 拓扑元数据（供 omm_bridge 等外部调用）。
pub fn Universe::atom_info(self : Universe) -> Array[UnivAtom] {
  self.atoms
}

///|
/// 单帧 Universe。
pub fn Universe::from_frame(
  atoms : Array[UnivAtom],
  coords : Array[Vector3],
  box~ : PBCBox = PBCBox::empty(),
) -> Universe {
  let traj : Array[Array[Vector3]] = Array(capacity=1)
  traj.push(coords)
  Universe::new(atoms, trajectory=traj, box~)
}

// ============================================================================
// AtomGroup
// ============================================================================

///|
/// 原子组 — Universe 的子选择。
pub struct AtomGroup {
  indices : Array[Int32]
  universe : Universe
}

///|
/// 选所有原子。
pub fn Universe::all_atoms(self : Universe) -> AtomGroup {
  let indices : Array[Int32] = Array(capacity=self.atoms.length())
  let mut i = 0
  while i < self.atoms.length() {
    indices.push(i)
    i = i + 1
  }
  AtomGroup::{ indices, universe: self, }
}

///|
/// 按索引子集选原子。
pub fn Universe::select_atoms_by_indices(
  self : Universe,
  indices : Array[Int32],
) -> AtomGroup {
  AtomGroup::{ indices, universe: self, }
}

///|
/// AtomGroup 原子数。
pub fn AtomGroup::n_atoms(self : AtomGroup) -> Int32 {
  self.indices.length()
}

///|
/// 按名字过滤。
pub fn AtomGroup::select_by_name(self : AtomGroup, name : String) -> AtomGroup {
  let indices : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    if self.universe.atoms[idx].name == name {
      indices.push(idx)
    }
    i = i + 1
  }
  AtomGroup::{ indices, universe: self.universe, }
}

///|
/// 按残基名过滤。
pub fn AtomGroup::select_by_resname(
  self : AtomGroup,
  resname : String,
) -> AtomGroup {
  let indices : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    if self.universe.atoms[idx].resname == resname {
      indices.push(idx)
    }
    i = i + 1
  }
  AtomGroup::{ indices, universe: self.universe, }
}

///|
/// 按链 ID 过滤。
pub fn AtomGroup::select_by_chain(
  self : AtomGroup,
  chain : String,
) -> AtomGroup {
  let indices : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    if self.universe.atoms[idx].chain == chain {
      indices.push(idx)
    }
    i = i + 1
  }
  AtomGroup::{ indices, universe: self.universe, }
}

///|
/// 按残基序号范围过滤。
pub fn AtomGroup::select_by_resid_range(
  self : AtomGroup,
  start : Int32,
  end : Int32,
) -> AtomGroup {
  let indices : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    let r = self.universe.atoms[idx].resid
    if r >= start && r <= end {
      indices.push(idx)
    }
    i = i + 1
  }
  AtomGroup::{ indices, universe: self.universe, }
}

///|
/// 选择重原子（非 H）。
pub fn AtomGroup::heavy_atoms(self : AtomGroup) -> AtomGroup {
  let indices : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    let sym = self.universe.atoms[idx].element
    let is_h = if sym.length() > 0 {
      sym == "H"
    } else {
      let name = self.universe.atoms[idx].name
      name.length() > 0 && name.unsafe_get(0) == 72 // 'H' ASCII
    }
    if !is_h {
      indices.push(idx)
    }
    i = i + 1
  }
  AtomGroup::{ indices, universe: self.universe, }
}

///|
/// AtomGroup 的原子索引数组（0-based，指向 Universe.atoms）。
pub fn AtomGroup::indices(self : AtomGroup) -> Array[Int32] {
  self.indices
}

///|
/// 所有原子的 mass 之和。
pub fn AtomGroup::total_mass(self : AtomGroup) -> Double {
  let mut total = 0.0
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    if idx >= 0 && idx < self.universe.atoms.length() {
      let m = self.universe.atoms[idx].mass
      total = total + (if m > 0.0 { m } else { 1.0 })
    }
    i = i + 1
  }
  total
}

///|
/// 所有原子的 name 数组。
pub fn AtomGroup::names(self : AtomGroup) -> Array[String] {
  let out : Array[String] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    if idx >= 0 && idx < self.universe.atoms.length() {
      out.push(self.universe.atoms[idx].name)
    }
    i = i + 1
  }
  out
}

///|
/// 所有原子的 resid 数组。
pub fn AtomGroup::resids(self : AtomGroup) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    if idx >= 0 && idx < self.universe.atoms.length() {
      out.push(self.universe.atoms[idx].resid)
    }
    i = i + 1
  }
  out
}

///|
/// 当前帧下的坐标（返回 `Array[Array[Double]]` 兼容格式）。
pub fn AtomGroup::positions(self : AtomGroup) -> Array[Array[Double]] {
  let v3s = self.coords(self.universe.current_frame)
  let out : Array[Array[Double]] = Array(capacity=v3s.length())
  let mut i = 0
  while i < v3s.length() {
    out.push([v3s[i].x, v3s[i].y, v3s[i].z])
    i = i + 1
  }
  out
}

///|
/// 子组内按选择表达式再筛选。
pub fn AtomGroup::select(self : AtomGroup, expr : String) -> AtomGroup {
  let selected = mda_select_universe(
    expr,
    self.universe.atoms,
    positions=Some(self.universe.positions_at(self.universe.current_frame)),
  )
  // selected 是 Array[Int32]（选中的 indices 列表）
  let out : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    // 检查 idx 是否在 selected 里
    let mut found = false
    let mut j = 0
    while j < selected.length() {
      if selected[j] == idx {
        found = true
        j = selected.length() // break
      }
      j = j + 1
    }
    if found {
      out.push(idx)
    }
    i = i + 1
  }
  AtomGroup::{ indices: out, universe: self.universe, }
}

///|
/// 选择 backbone 原子（PDB 标准：N, CA, C, O）。
pub fn AtomGroup::backbone(self : AtomGroup) -> AtomGroup {
  let names = ["N", "CA", "C", "O"]
  let indices : Array[Int32] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    let n = self.universe.atoms[idx].name
    let mut j = 0
    while j < names.length() {
      if n == names[j] {
        indices.push(idx)
        j = names.length() // break
      }
      j = j + 1
    }
    i = i + 1
  }
  AtomGroup::{ indices, universe: self.universe, }
}

///|
/// AtomGroup 在某帧的坐标。
pub fn AtomGroup::coords(self : AtomGroup, frame : Int32) -> Array[Vector3] {
  let fr = self.universe.frame(frame)
  if fr.length() == 0 {
    return Array(capacity=0)
  }
  let result : Array[Vector3] = Array(capacity=self.indices.length())
  let mut i = 0
  while i < self.indices.length() {
    let idx = self.indices[i]
    if idx >= 0 && idx < fr.length() {
      result.push(fr[idx])
    }
    i = i + 1
  }
  result
}

///|
/// 几何中心（无权重，默认 frame 0）。
pub fn AtomGroup::center_of_geometry(
  self : AtomGroup,
  frame~ : Int32 = 0,
) -> Vector3 {
  let cs = self.coords(frame)
  let n = cs.length()
  if n == 0 {
    return Vector3::new(0.0, 0.0, 0.0)
  }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  let mut i = 0
  while i < n {
    sx = sx + cs[i].x
    sy = sy + cs[i].y
    sz = sz + cs[i].z
    i = i + 1
  }
  Vector3::new(sx / n.to_double(), sy / n.to_double(), sz / n.to_double())
}

///|
/// 质心（按 mass 权重，默认 frame 0）。
pub fn AtomGroup::center_of_mass(
  self : AtomGroup,
  frame~ : Int32 = 0,
) -> Vector3 {
  let cs = self.coords(frame)
  let n = cs.length()
  if n == 0 {
    return Vector3::new(0.0, 0.0, 0.0)
  }
  let mut total_mass = 0.0
  let mut wx = 0.0
  let mut wy = 0.0
  let mut wz = 0.0
  let mut i = 0
  while i < n {
    let idx = self.indices[i]
    let m = if idx >= 0 && idx < self.universe.atoms.length() {
      let am = self.universe.atoms[idx].mass
      if am > 0.0 {
        am
      } else {
        1.0
      }
    } else {
      1.0
    }
    total_mass = total_mass + m
    wx = wx + m * cs[i].x
    wy = wy + m * cs[i].y
    wz = wz + m * cs[i].z
    i = i + 1
  }
  if total_mass > 0.0 {
    Vector3::new(wx / total_mass, wy / total_mass, wz / total_mass)
  } else {
    Vector3::new(0.0, 0.0, 0.0)
  }
}

///|
/// 回转半径（按 mass 权重，默认 frame 0）。
pub fn AtomGroup::radius_of_gyration(
  self : AtomGroup,
  frame~ : Int32 = 0,
) -> Double {
  let cs = self.coords(frame)
  let n = cs.length()
  if n == 0 {
    return 0.0
  }
  let com = self.center_of_mass(frame~)
  let mut total_mass = 0.0
  let mut rg_sum = 0.0
  let mut i = 0
  while i < n {
    let idx = self.indices[i]
    let m = if idx >= 0 && idx < self.universe.atoms.length() {
      let am = self.universe.atoms[idx].mass
      if am > 0.0 {
        am
      } else {
        1.0
      }
    } else {
      1.0
    }
    total_mass = total_mass + m
    let dx = cs[i].x - com.x
    let dy = cs[i].y - com.y
    let dz = cs[i].z - com.z
    rg_sum = rg_sum + m * (dx * dx + dy * dy + dz * dz)
    i = i + 1
  }
  if total_mass > 0.0 {
    (rg_sum / total_mass).sqrt()
  } else {
    0.0
  }
}

// ============================================================================
// GRO → Universe 便捷构建
// ============================================================================

///|
/// 从 GRO 文本构建 Universe（单帧）。
///
/// 复用已有的 `gro_parse`（structure_gro.mbt）。
/// GroAtom 坐标单位是 nm，这里转成 Å。
pub fn Universe::from_gro(text : String) -> Universe {
  let gro = gro_parse(text)
  let natoms = gro.atoms.length()
  let atoms : Array[UnivAtom] = Array(capacity=natoms)
  let coords : Array[Vector3] = Array(capacity=natoms)
  let mut i = 0
  while i < natoms {
    let ga = gro.atoms[i]
    let name = ga.atom_name.trim().to_owned()
    let resname = ga.res_name.trim().to_owned()
    let element = if name.length() > 0 {
      let c = name.unsafe_get(0)
      if c >= 65 && c <= 90 {
        // uppercase letter -> first char as element (e.g. "CA" -> "C")
        let sb = StringBuilder()
        sb.write_char(c.unsafe_to_char())
        sb.to_string()
      } else if c >= 48 && c <= 57 {
        // digit prefix, take second char
        if name.length() > 1 {
          let c2 = name.unsafe_get(1)
          let sb = StringBuilder()
          sb.write_char(c2.unsafe_to_char())
          sb.to_string()
        } else {
          ""
        }
      } else {
        ""
      }
    } else {
      ""
    }
    // nm → Å
    let x = ga.x * 10.0
    let y = ga.y * 10.0
    let z = ga.z * 10.0
    atoms.push(
      UnivAtom::new(
        name~,
        resname~,
        resid=ga.res_id,
        chain="A",
        segid="SYS",
        element~,
        is_protein=mda_universe_is_protein(resname),
        is_water=resname == "WAT" || resname == "HOH" || resname == "SOL",
        is_backbone=name == "N" || name == "CA" || name == "C" || name == "O",
      ),
    )
    coords.push(Vector3::new(x, y, z))
    i = i + 1
  }
  let traj : Array[Array[Vector3]] = Array(capacity=1)
  traj.push(coords)
  let box = if gro.box.length() >= 3 && gro.box[0].length() >= 3 {
    PBCBox::from_lengths(
      gro.box[0][0] * 10.0,
      gro.box[1][1] * 10.0,
      gro.box[2][2] * 10.0,
    )
  } else {
    PBCBox::empty()
  }
  Universe::new(atoms, trajectory=traj, box~)
}

///|
/// 判断是否标准蛋白质残基名。
fn mda_universe_is_protein(resname : String) -> Bool {
  match resname {
    "ALA"
    | "ARG"
    | "ASN"
    | "ASP"
    | "CYS"
    | "GLN"
    | "GLU"
    | "GLY"
    | "HIS"
    | "ILE"
    | "LEU"
    | "LYS"
    | "MET"
    | "PHE"
    | "PRO"
    | "SER"
    | "THR"
    | "TRP"
    | "TYR"
    | "VAL" => true
    _ => false
  }
}

///|
/// 从 LAMMPS Dump text 构建 Universe（支持多帧）。
///
/// 直接调用 BioSeqs 已有的 `lammps_dump_read()`，
/// 把 LammpsFrame → UnivAtom + Vector3 trajectory。
pub fn Universe::from_lammps(text : String) -> Universe raise LammpsDumpError {
  let frames = lammps_dump_read(text)
  let nf = frames.length()
  if nf == 0 {
    return Universe::new(Array(capacity=0), trajectory=Array(capacity=0))
  }
  let first = frames[0]
  let natoms = first.natoms
  let atoms : Array[UnivAtom] = Array(capacity=natoms)
  let mut i = 0
  while i < natoms {
    let la = first.atoms[i]
    let tn = la.atom_type.to_string()
    atoms.push(
      UnivAtom::new(
        name="A" + tn,
        resname="LAM",
        resid=la.mol,
        chain="A",
        segid="SYS",
        element="X",
        mass=1.0,
        charge=0.0,
        is_protein=false,
        is_backbone=false,
      ),
    )
    i = i + 1
  }
  let traj : Array[Array[Vector3]] = Array(capacity=nf)
  let mut f = 0
  while f < nf {
    let fr = frames[f]
    let fv : Array[Vector3] = Array(capacity=fr.atoms.length())
    let mut k = 0
    while k < fr.atoms.length() {
      fv.push(Vector3::new(fr.atoms[k].x, fr.atoms[k].y, fr.atoms[k].z))
      k = k + 1
    }
    traj.push(fv)
    f = f + 1
  }
  let box = if first.box_lo.length() >= 3 && first.box_hi.length() >= 3 {
    let dx = first.box_hi[0] - first.box_lo[0]
    let dy = first.box_hi[1] - first.box_lo[1]
    let dz = first.box_hi[2] - first.box_lo[2]
    PBCBox::from_lengths(dx, dy, dz)
  } else {
    PBCBox::empty()
  }
  Universe::new(atoms, trajectory=traj, box~)
}
