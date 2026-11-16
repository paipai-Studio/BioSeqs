// ============================================================================
// MDAnalysis 辅助函数（test 文件引用）
// ============================================================================

///|
/// 从原子名推断元素符号。
///
/// 规则：前 1-2 个字母，数字前缀跳过。
/// 例子："CA" -> "C", "N" -> "N", "1HG2" -> "H", "CL1" -> "Cl", "FE" -> "Fe"
pub fn mda_infer_element(name : String) -> String {
  let n = name.length()
  if n == 0 { return "" }
  let mut start = 0
  // 跳过开头数字
  while start < n {
    let c = name.unsafe_get(start)
    if c >= 48 && c <= 57 { // '0'-'9'
      start = start + 1
    } else {
      break
    }
  }
  if start >= n { return "" }
  let first = name.unsafe_get(start)
  // 大写字母开头
  if first >= 65 && first <= 90 { // 'A'-'Z'
    if start + 1 < n {
      let second = name.unsafe_get(start + 1)
      let sec_is_alpha = (second >= 65 && second <= 90) || (second >= 97 && second <= 122)
      if sec_is_alpha {
        let sb = StringBuilder()
        sb.write_char(first.unsafe_to_char())
        // second 转小写（如果是大写 +32）
        let sec = if second >= 65 && second <= 90 { second + 32 } else { second }
        sb.write_char(sec.unsafe_to_char())
        let two = sb.to_string()
        // 第二位小写 → 标准两字母元素
        // 第二位大写 → 白名单
        if second >= 97 && second <= 122 {
          return two
        }
        if mda_infer_two_letter_uppercase(two) {
          return two
        }
      }
    }
    // 默认取一位
    let sb = StringBuilder()
    sb.write_char(first.unsafe_to_char())
    return sb.to_string()
  }
  return ""
}

///|
/// 当第二位是大写字母时，判断是否属于两字母元素白名单。
fn mda_infer_two_letter_uppercase(two : String) -> Bool {
  match two.to_lower() {
    // 只保留在 PDB 原子名中不太可能歧义的两字母元素
    "fe" | "zn" | "na" | "mg" | "al" | "si" | "cl" | "br" | "se" |
    "cu" | "ag" | "be" => true
    _ => false
  }
}

///|
/// 元素符号 → 标准原子量（amu）。
///
/// 未知元素回退 12.011（碳）。
pub fn mda_element_mass(element : String) -> Double {
  match element.to_upper() {
    "H"   => 1.008
    "HE"  => 4.0026
    "LI"  => 6.941
    "BE"  => 9.0122
    "B"   => 10.811
    "C"   => 12.011
    "N"   => 14.007
    "O"   => 15.999
    "F"   => 18.998
    "NE"  => 20.180
    "NA"  => 22.990
    "MG"  => 24.305
    "AL"  => 26.982
    "SI"  => 28.086
    "P"   => 30.974
    "S"   => 32.065
    "CL"  => 35.453
    "AR"  => 39.948
    "K"   => 39.098
    "CA"  => 40.078
    "FE"  => 55.845
    "CO"  => 58.933
    "NI"  => 58.693
    "CU"  => 63.546
    "ZN"  => 65.38
    "SE"  => 78.971
    "BR"  => 79.904
    "AG"  => 107.868
    "I"   => 126.904
    _     => 12.011
  }
}

///|
/// 判断水残基名（大小写不敏感）。
pub fn mda_is_water_resname(resname : String) -> Bool {
  match resname.to_upper() {
    "WAT" | "HOH" | "SOL" | "TIP3" | "TIP4" | "SPC" | "SPCE" | "WATER" | "H2O" => true
    _ => false
  }
}

// ============================================================================
// PDB → Universe
// ============================================================================

///|
/// 从 PDB 文本（支持多 MODEL）构建 Universe。
///
/// 每个 MODEL/ENDMDL 块是一帧；拓扑从第一帧提取。
pub fn universe_from_pdb_text(pdb : String) -> Universe {
  let structure = parse_pdb(pdb)
  match PDBTrajectory::from_structure(structure) {
    Some(traj) => {
      let natoms = traj.n_atoms
      let atoms : Array[UnivAtom] = Array(capacity=natoms)
      let mut i = 0
      while i < natoms {
        let name = traj.atom_names[i]
        let resname = traj.resnames[i]
        let element = if traj.elements[i].length() > 0 {
          traj.elements[i]
        } else {
          mda_infer_element(name)
        }
        let chain_sb = StringBuilder()
        let cid = traj.chainids[i]
        if cid != ' ' { chain_sb.write_char(cid) }
        atoms.push(UnivAtom::new(
          name=name,
          resname=resname,
          chain=chain_sb.to_string(),
          segid="SYS",
          element=element,
          atom_type="",
          resid=traj.resseqs[i],
          mass=mda_element_mass(element),
          is_protein=mda_universe_is_protein(resname),
          is_water=mda_is_water_resname(resname),
          is_backbone=(name == "N" || name == "CA" || name == "C" || name == "O"),
        ))
        i = i + 1
      }
      let out_frames : Array[Array[Vector3]] = Array(capacity=traj.n_frames)
      let mut f = 0
      while f < traj.n_frames {
        out_frames.push(traj.frames[f])
        f = f + 1
      }
      Universe::new(atoms, trajectory=out_frames)
    }
    None => {
      // 空 PDB
      Universe::new(Array(capacity=0))
    }
  }
}

// ============================================================================
// PQR → Universe
// ============================================================================

///|
/// 从 PQR 文本构建 Universe（单帧）。
///
/// PQR 格式类似 PDB ATOM 记录，但最后两个字段是 charge 和 radius。
/// 列宽：serial(7-11), name(13-16), resname(18-20), chain(22), resid(23-26),
///       x(31-38), y(39-46), z(47-54), charge(55-62), radius(63-70)
pub fn universe_from_pqr_text(pqr : String) -> Universe {
  let lines = pqr.split("\n").to_array()
  let atoms : Array[UnivAtom] = Array(capacity=lines.length())
  let coords : Array[Vector3] = Array(capacity=lines.length())
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].to_owned()
    // 只处理 ATOM / HETATM
    let rec = if line.length() >= 6 {
      line[:6].to_owned().trim().to_owned().to_upper()
    } else { "" }
    if rec != "ATOM" && rec != "HETATM" {
      i = i + 1
      continue
    }
    // PQR 固定列宽（PDB 标准），但用空格 split 更健壮
    let fields : Array[String] = Array(capacity=12)
    for w in line.split(" ") {
      let t = w.to_owned().trim().to_owned()
      if t.length() > 0 { fields.push(t) }
    }
    // fields: [record, serial, name, resname, chain, resid, x, y, z, charge, radius]
    // 但 chain 可能跟 resid 合并（"A1"），所以直接按名字/列取
    let name = if line.length() > 16 { line[12:16].to_owned().trim().to_owned() } else { "" }
    let resname = if line.length() > 20 { line[17:20].to_owned().trim().to_owned() } else { "" }
    let chain = if line.length() > 21 {
      let c = line.unsafe_get(21).to_int()
      if c == 32 { "" } else {
        let sb = StringBuilder()
        sb.write_char(c.unsafe_to_char())
        sb.to_string()
      }
    } else { "" }
    let resid = if line.length() > 26 {
      let v = parse_int(line[22:26].to_owned().trim().to_owned())
      if v > 0 { v } else { 1 }
    } else { 1 }
    let x = if line.length() > 38 {
      match parse_double(line[30:38].to_owned().trim().to_owned()) {
        Some(v) => v
        None => 0.0
      }
    } else { 0.0 }
    let y = if line.length() > 46 {
      match parse_double(line[38:46].to_owned().trim().to_owned()) {
        Some(v) => v
        None => 0.0
      }
    } else { 0.0 }
    let z = if line.length() > 54 {
      match parse_double(line[46:54].to_owned().trim().to_owned()) {
        Some(v) => v
        None => 0.0
      }
    } else { 0.0 }
    // PQR charge / radius: 最后两个空格分隔字段更可靠
    let charge = if fields.length() >= 10 {
      match parse_double(fields[9]) {
        Some(v) => v
        None => 0.0
      }
    } else { 0.0 }
    let element = mda_infer_element(name)
    let mass = mda_element_mass(element)
    atoms.push(UnivAtom::new(
      name=name,
      resname=resname,
      chain=chain,
      segid="SYS",
      element=element,
      atom_type="",
      resid=resid,
      mass=mass,
      charge=charge,
      is_protein=mda_universe_is_protein(resname),
      is_water=mda_is_water_resname(resname),
      is_backbone=(name == "N" || name == "CA" || name == "C" || name == "O"),
    ))
    coords.push(Vector3::new(x, y, z))
    i = i + 1
  }
  let traj : Array[Array[Vector3]] = Array(capacity=1)
  traj.push(coords)
  Universe::new(atoms, trajectory=traj)
}

// ============================================================================
// PSF + DCD → Universe
// ============================================================================

///|
/// 从 CHARMM PSF 拓扑 + DCD 轨迹构建 Universe。
///
/// 复用 `structure_psf` 的 PsfTopology 和 `dcd` 的 parse_dcd。
pub fn universe_from_psf_dcd(psf_text : String, dcd_text : String) -> Universe {
  let top = try { psf_parse(psf_text) } catch { _ => return Universe::new(Array(capacity=0)) }
  let traj = parse_dcd(dcd_text)
  let natoms = top.atoms.length()
  // PSF 原子（PsfAtom）→ UnivAtom
  let atoms : Array[UnivAtom] = Array(capacity=natoms)
  let mut i = 0
  while i < natoms {
    let pa = top.atoms[i]
    let element = mda_infer_element(pa.atom_name)
    let mass = if pa.mass > 0.0 { pa.mass } else { mda_element_mass(element) }
    let is_prot = mda_universe_is_protein(pa.resname)
    let is_wat = mda_is_water_resname(pa.resname)
    let is_bb = (pa.atom_name == "N" || pa.atom_name == "CA" || pa.atom_name == "C" || pa.atom_name == "O")
    atoms.push(UnivAtom::new(
      name=pa.atom_name,
      resname=pa.resname,
      chain="",
      segid=pa.segid,
      element=element,
      atom_type=pa.atom_type,
      resid=pa.resid,
      mass=mass,
      charge=pa.charge,
      is_protein=is_prot,
      is_water=is_wat,
      is_backbone=is_bb,
    ))
    i = i + 1
  }
  // DCD frames → Array[Array[Vector3]]
  let frames = traj.frames
  let out_frames : Array[Array[Vector3]] = Array(capacity=frames.length())
  let mut f = 0
  while f < frames.length() {
    let frame = frames[f]
    let v3s : Array[Vector3] = Array(capacity=frame.length())
    let mut a = 0
    while a < frame.length() {
      let c = frame[a]
      if c.length() >= 3 {
        v3s.push(Vector3::new(c[0], c[1], c[2]))
      }
      a = a + 1
    }
    out_frames.push(v3s)
    f = f + 1
  }
  Universe::new(atoms, trajectory=out_frames)
}
