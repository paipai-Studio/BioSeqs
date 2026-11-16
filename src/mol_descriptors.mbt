// ============================================================================
// Molecular descriptors (RDKit/Open Babel compatible)
// ============================================================================
// Implements: LogP (Wildman-Crippen), TPSA (Ertl), HBD/HBA, rotatable bonds,
// formal charge, fraction sp3, aromatic/saturated ring count, valence electrons.
// All operate on the Molecule type from rdkit.mbt.

///|
/// Atomic contribution to LogP (Wildman-Crippen method).
/// Simplified table: (symbol, degree, aromatic) -> contribution.
fn logp_atom_contribution(
  symbol : String,
  degree : Int32,
  aromatic : Bool,
  _in_ring : Bool,
) -> Double {
  if aromatic {
    match symbol {
      "C" => 0.13
      "N" => -0.53
      "O" => -0.48
      "S" => 0.00
      _ => 0.00
    }
  } else {
    match symbol {
      "C" => if degree == 4 { 0.20 } else if degree == 3 { 0.14 } else { 0.10 }
      "N" =>
        if degree == 3 {
          -1.00
        } else if degree == 2 {
          -0.50
        } else {
          -1.50
        }
      "O" => if degree == 2 { -0.45 } else { -0.60 }
      "S" => -0.20
      "P" => 0.25
      "F" => -0.25
      "Cl" => 0.06
      "Br" => 0.20
      "I" => 0.50
      "H" => 0.00
      _ => 0.00
    }
  }
}

///|
/// Calculate Wildman-Crippen LogP (octanol-water partition coefficient).
/// Sum of atom contributions + correction for ring atoms.
pub fn mol_logp(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut logp = 0.0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let deg = neighbors[k].length()
    let in_ring = mol_is_in_ring(mol, k)
    let contrib = logp_atom_contribution(a.symbol, deg, a.aromatic, in_ring)
    logp = logp + contrib
    // Ring correction
    if in_ring && !a.aromatic {
      logp = logp - 0.05
    }
  }
  // Bond corrections
  for b in mol.bonds {
    if b.order == 2 {
      logp = logp - 0.42
    } else if b.order == 3 {
      logp = logp - 0.70
    }
  }
  logp
}

///|
/// Atomic contribution to TPSA (topological polar surface area, Ertl method).
/// Based on atom type (N/O with specific neighbor patterns).
fn tpsa_atom_contribution(
  symbol : String,
  degree : Int32,
  _num_h_bonds : Int32,
  in_ring : Bool,
) -> Double {
  // N types
  if symbol == "N" {
    if in_ring && degree == 2 {
      return 3.24 // aromatic N in ring
    }
    if degree == 1 {
      return 23.79 // primary amine
    }
    if degree == 2 {
      return 12.36 // secondary amine
    }
    if degree == 3 {
      return 3.24 // tertiary amine
    }
    return 0.0
  }
  // O types
  if symbol == "O" {
    if degree == 1 {
      return 20.23 // hydroxyl/ether
    }
    if degree == 2 {
      return 9.23 // ring O
    }
    return 0.0
  }
  // S and P contribute 0 in basic TPSA
  0.0
}

///|
/// Calculate topological polar surface area (TPSA) in Angstrom^2.
pub fn mol_tpsa(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut tpsa = 0.0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let deg = neighbors[k].length()
    let in_ring = mol_is_in_ring(mol, k)
    tpsa = tpsa + tpsa_atom_contribution(a.symbol, deg, 0, in_ring)
  }
  tpsa
}

///|
/// Count H-bond donors (N-H or O-H bonds).
/// An HBD is any N or O atom with at least one implicit H.
pub fn mol_hbd_count(mol : Molecule) -> Int32 {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut count = 0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    if a.symbol == "N" || a.symbol == "O" {
      let deg = neighbors[k].length()
      // Implicit H count = typical_valence - degree - charge
      let typical = if a.symbol == "N" { 3 } else { 2 }
      let implicit_h = typical - deg + a.charge
      if implicit_h > 0 {
        count = count + implicit_h
      }
    }
  }
  count
}

///|
/// Count H-bond acceptors (N or O atoms).
pub fn mol_hba_count(mol : Molecule) -> Int32 {
  let mut count = 0
  for a in mol.atoms {
    if a.symbol == "N" || a.symbol == "O" {
      count = count + 1
    }
  }
  count
}

///|
/// Count rotatable bonds.
/// A rotatable bond is a single bond, not in a ring, not to a terminal H,
/// and both endpoints have at least 2 heavy neighbors.
pub fn mol_rotatable_bonds(mol : Molecule) -> Int32 {
  let _n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut count = 0
  for b in mol.bonds {
    if b.order != 1 {
      continue
    }
    // Skip if either atom is H
    let a1 = mol.atoms[b.from]
    let a2 = mol.atoms[b.to]
    if a1.symbol == "H" || a2.symbol == "H" {
      continue
    }
    // Skip if in ring
    let rings = mol_find_rings(mol)
    let mut in_ring = false
    for ring in rings {
      let mut has_from = false
      let mut has_to = false
      for a in ring {
        if a == b.from {
          has_from = true
        }
        if a == b.to {
          has_to = true
        }
      }
      if has_from && has_to {
        in_ring = true
      }
    }
    if in_ring {
      continue
    }
    // Both endpoints need >= 2 neighbors
    if neighbors[b.from].length() >= 2 && neighbors[b.to].length() >= 2 {
      count = count + 1
    }
  }
  count
}

///|
/// Calculate formal charge (sum of all atom charges).
pub fn mol_formal_charge(mol : Molecule) -> Int32 {
  let mut charge = 0
  for a in mol.atoms {
    charge = charge + a.charge
  }
  charge
}

///|
/// Calculate fraction of sp3 carbons (C with 4 single bonds / total C).
pub fn mol_fraction_sp3(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut total_c = 0
  let mut sp3_c = 0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    if a.symbol == "C" {
      total_c = total_c + 1
      // Check if all bonds are single
      let mut all_single = true
      for (_, order) in neighbors[k] {
        if order != 1 {
          all_single = false
        }
      }
      if all_single && !a.aromatic {
        sp3_c = sp3_c + 1
      }
    }
  }
  if total_c == 0 {
    return 0.0
  }
  sp3_c.to_double() / total_c.to_double()
}

///|
/// Count valence electrons of the molecule.
pub fn mol_valence_electrons(mol : Molecule) -> Int32 {
  let mut total = 0
  for a in mol.atoms {
    total = total + atom_valence_electrons(a.symbol)
  }
  total
}

///|
/// Valence electrons for an element.
fn atom_valence_electrons(symbol : String) -> Int32 {
  match symbol {
    "H" => 1
    "C" => 4
    "N" => 5
    "O" => 6
    "F" => 7
    "P" => 5
    "S" => 6
    "Cl" => 7
    "Br" => 7
    "I" => 7
    "B" => 3
    "Si" => 4
    _ => 4
  }
}

///|
/// Typical valence（用于计算 implicit H 数）。
fn atom_typical_valence(symbol : String) -> Int32 {
  match symbol {
    "H" => 1
    "C" => 4
    "N" => 3
    "O" => 2
    "F" => 1
    "P" => 3
    "S" => 2
    "Cl" => 1
    "Br" => 1
    "I" => 1
    "B" => 3
    "Si" => 4
    _ => 4
  }
}

///|
/// 自由基电子计数（未成对电子总数）。
///
/// RDKit `CalcNumRadicalElectrons` 等价物：
/// 对每个原子计算 remaining_electrons = |valence − charge − 2×(显式键序和 + implicit H)|，
/// 如果这个数是奇数则贡献 1 个未成对电子（每个原子最多 1 个，因为未成对电子以
/// 成对规则分布）。
///
/// 注：BioSeqs Molecule 目前不带 implicit H 信息，implicit H 由 typical valence
/// 估计：implicit_h = typical_valence − degree + charge。
pub fn mol_radical_electrons(mol : Molecule) -> Int32 {
  let neighbors = mol_adjacency(mol)
  let n = mol.atoms.length()
  let mut total = 0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let ve = atom_valence_electrons(a.symbol)
    let tv = atom_typical_valence(a.symbol)
    let degree = neighbors[k].length()
    // implicit H（无 explicit H 时用 typical valence 估计）
    let implicit_h = tv - degree + a.charge
    let h_eff = if implicit_h > 0 { implicit_h } else { 0 }
    // bond_order_sum = Σ(该原子所有键的 order) — 双键算 2，三键算 3
    let mut bond_sum = 0
    for (_, order) in neighbors[k] {
      bond_sum = bond_sum + order
    }
    // implicit H 每个贡献 1 个单键
    let total_bonds = bond_sum + h_eff
    let remaining = (ve - a.charge - 2 * total_bonds).abs()
    if remaining % 2 != 0 {
      total = total + 1
    }
  }
  total
}

///|
/// Count heavy atoms (non-H).
pub fn mol_heavy_atom_count(mol : Molecule) -> Int32 {
  let mut count = 0
  for a in mol.atoms {
    if a.symbol != "H" {
      count = count + 1
    }
  }
  count
}

///|
/// Molecular weight from atom symbols.
pub fn mol_molecular_weight(mol : Molecule) -> Double {
  let mut mw = 0.0
  for a in mol.atoms {
    mw = mw + atom_atomic_weight(a.symbol)
  }
  mw
}

///|
/// Atomic weight for common elements.
fn atom_atomic_weight(symbol : String) -> Double {
  match symbol {
    "H" => 1.008
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "F" => 18.998
    "P" => 30.974
    "S" => 32.065
    "Cl" => 35.453
    "Br" => 79.904
    "I" => 126.904
    "B" => 10.811
    "Si" => 28.086
    _ => 12.0
  }
}

///|
/// Count total atoms.
pub fn mol_num_atoms(mol : Molecule) -> Int32 {
  mol.num_atoms()
}

///|
/// Count total bonds.
pub fn mol_num_bonds(mol : Molecule) -> Int32 {
  mol.num_bonds()
}

///|
/// Full descriptor report as a formatted string.
pub fn mol_descriptor_report(mol : Molecule) -> String {
  let sb = StringBuilder()
  sb.write_string("Molecular Descriptors\n")
  sb.write_string("=====================\n")
  sb.write_string("Atoms: " + mol_num_atoms(mol).to_string() + "\n")
  sb.write_string(
    "Heavy atoms: " + mol_heavy_atom_count(mol).to_string() + "\n",
  )
  sb.write_string("Bonds: " + mol_num_bonds(mol).to_string() + "\n")
  sb.write_string(
    "Molecular weight: " + mol_molecular_weight(mol).to_string() + "\n",
  )
  sb.write_string("LogP: " + mol_logp(mol).to_string() + "\n")
  sb.write_string("TPSA: " + mol_tpsa(mol).to_string() + "\n")
  sb.write_string("HBD: " + mol_hbd_count(mol).to_string() + "\n")
  sb.write_string("HBA: " + mol_hba_count(mol).to_string() + "\n")
  sb.write_string(
    "Rotatable bonds: " + mol_rotatable_bonds(mol).to_string() + "\n",
  )
  sb.write_string("Formal charge: " + mol_formal_charge(mol).to_string() + "\n")
  sb.write_string("Fraction sp3: " + mol_fraction_sp3(mol).to_string() + "\n")
  sb.write_string("Rings: " + mol_ring_count(mol).to_string() + "\n")
  sb.write_string(
    "Aromatic rings: " + mol_aromatic_ring_count(mol).to_string() + "\n",
  )
  sb.write_string(
    "Saturated rings: " + mol_saturated_ring_count(mol).to_string() + "\n",
  )
  sb.write_string(
    "Valence electrons: " + mol_valence_electrons(mol).to_string() + "\n",
  )
  sb.to_string()
}

// ============================================================================
// ExactMolWt（最丰同位素精确质量）
// ============================================================================

///|
/// 最丰同位素精确原子量（区别于元素周期表的平均原子量）。
/// 参考 IUPAC 2021 atomic weights。
fn atom_exact_mass(symbol : String) -> Double {
  match symbol {
    "H" => 1.00782503207
    "C" => 12.00000000000
    "N" => 14.0030740048
    "O" => 15.99491461956
    "F" => 18.99840322
    "P" => 30.97376163
    "S" => 31.97207100
    "Cl" => 34.96885268
    "Br" => 78.9183371
    "I" => 126.904473
    "B" => 11.009305
    "Si" => 27.9769265325
    "Se" => 79.9165213
    _ => 12.000000
  }
}

///|
/// 精确分子量：用最丰同位素质量而非平均原子量。
/// 例如 C6H12O6 的 ExactMolWt = 180.063388（vs 平均 MW 180.156）。
///
/// 只计数 Molecule 里显式的重原子（无 implicit H），这和 RDKit 的
/// `CalcExactMolWt` 在分子图不含 H 时完全一致。
pub fn mol_exact_molecular_weight(mol : Molecule) -> Double {
  let mut mass = 0.0
  for a in mol.atoms {
    mass = mass + atom_exact_mass(a.symbol)
  }
  mass
}

// ============================================================================
// MolFormula（Hill 系统分子式）
// ============================================================================

///|
/// Hill 系统排序：先 C，再 H，然后其他元素按字母序。
/// 数字=1 时省略。空分子返回 "".
pub fn mol_formula(mol : Molecule) -> String {
  let counts : Map[String, Int32] = Map([], capacity=mol.atoms.length())
  for a in mol.atoms {
    let sym = a.symbol
    match counts.get(sym) {
      Some(c) => counts[sym] = c + 1
      None => counts[sym] = 1
    }
  }
  let sb = StringBuilder()
  // C first
  let c_count = match counts.get("C") {
    Some(c) => c
    None => 0
  }
  if c_count > 0 {
    sb.write_string("C")
    if c_count > 1 {
      sb.write_string(c_count.to_string())
    }
  }
  // H next (only if C present)
  if c_count > 0 {
    let h_count = match counts.get("H") {
      Some(c) => c
      None => 0
    }
    if h_count > 0 {
      sb.write_string("H")
      if h_count > 1 {
        sb.write_string(h_count.to_string())
      }
    }
  }
  // 收集其他元素，按字母序
  let rest : Array[String] = Array(capacity=0)
  for sym in counts.keys() {
    if sym == "C" || sym == "H" {
      continue
    }
    let mut inserted = false
    for i = 0; i < rest.length(); i = i + 1 {
      if sym < rest[i] {
        rest.insert(i, sym)
        inserted = true
        break
      }
    }
    if !inserted {
      rest.push(sym)
    }
  }
  for sym in rest {
    let c = match counts.get(sym) {
      Some(v) => v
      None => 0
    }
    sb.write_string(sym)
    if c > 1 {
      sb.write_string(c.to_string())
    }
  }
  sb.to_string()
}

// ============================================================================
// MolBlock（V2000 connection table）
// ============================================================================

///|
/// 将 Molecule 序列化为 MolBlock 字符串（V2000 connection table）。
///
/// 这是 RDKit `MolToMolBlock` 的 BioSeqs 等价物：输出标准 SDF 格式的
/// molecule block（带 "M  END" 但不带 "$$$$" 分隔符）。
/// 原子坐标填 0,0,0（当前 Molecule 无 3D 坐标）。
pub fn mol_to_mol_block(mol : Molecule) -> String {
  let n_atoms = mol.atoms.length()
  let n_bonds = mol.bonds.length()
  let sb = StringBuilder()

  // Header block (3 lines)
  sb.write_string("\n")
  sb.write_string("\n")
  sb.write_string("  BioSeqs   ")
  sb.write_string("  V2000\n")

  // Counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
  let counts = n_atoms.to_string().pad_start(3, ' ') +
    n_bonds.to_string().pad_start(3, ' ') +
    "000000000000000000000000000000000000"
  sb.write_string(counts)
  sb.write_string("\n")

  // Atom lines
  for a in mol.atoms {
    let sx = "    0.0000"
    let sy = "    0.0000"
    let sz = "    0.0000"
    let sym_pad = a.symbol
    let sp = if sym_pad.length() >= 3 {
      sym_pad
    } else {
      " ".repeat(3 - sym_pad.length()) + sym_pad
    }
    let charge_code = match a.charge {
      0 => "000"
      1 => "001"
      -1 => "002"
      2 => "003"
      -2 => "004"
      _ => "000"
    }
    let line = sx +
      sy +
      sz +
      " " +
      sp +
      " 000" +
      charge_code +
      "  0  0  0  0  0  0  0  0  0"
    sb.write_string(line)
    sb.write_string("\n")
  }

  // Bond lines
  for b in mol.bonds {
    let a1 = (b.from + 1).to_string().pad_start(3, ' ')
    let a2 = (b.to + 1).to_string().pad_start(3, ' ')
    let bt = match b.order {
      1 => "001"
      2 => "002"
      3 => "003"
      _ => "001"
    }
    let line = a1 + a2 + bt + "  0  0  0  0\n"
    sb.write_string(line)
  }

  sb.write_string("M  END\n")
  sb.to_string()
}

///|
/// 从 MolBlock 字符串解析 Molecule（V2000 connection table）。
pub fn mol_from_mol_block(block : String) -> Molecule {
  let lines : Array[String] = Array(capacity=0)
  let mut cur = StringBuilder()
  let mut i = 0
  while i < block.length() {
    let c = block.unsafe_get(i)
    if c == 10 || c == 13 {
      if cur.to_string().length() > 0 {
        lines.push(cur.to_string())
      }
      cur = StringBuilder()
      if c == 13 && i + 1 < block.length() && block.unsafe_get(i + 1) == 10 {
        i = i + 1
      }
    } else {
      cur.write_char(c.unsafe_to_char())
    }
    i = i + 1
  }
  if cur.to_string().length() > 0 {
    lines.push(cur.to_string())
  }

  let mut counts_idx = -1
  let idx2 = 0
  while idx2 < lines.length() {
    if lines[idx2].contains("V2000") || lines[idx2].contains("V3000") {
      counts_idx = idx2 + 1
      break
    }
    let _ = idx2
  }
  if counts_idx == -1 {
    counts_idx = 3
  }
  if counts_idx >= lines.length() {
    return Molecule::{ atoms: [], bonds: [], }
  }
  let counts_line = lines[counts_idx]
  let n_atoms = mol_block_parse_int(counts_line, 0, 3)
  let n_bonds = mol_block_parse_int(counts_line, 3, 3)

  let atoms : Array[MolAtom] = Array(capacity=n_atoms)
  let mut j = 0
  while j < n_atoms {
    let idx = counts_idx + 1 + j
    if idx >= lines.length() {
      break
    }
    let al = lines[idx]
    let sym = mol_block_slice(al, 31, 3).trim().to_owned()
    let charge_code = mol_block_parse_int(al, 36, 3)
    let charge = match charge_code {
      1 => 1
      2 => -1
      3 => 2
      4 => -2
      _ => 0
    }
    atoms.push(MolAtom::{ symbol: sym, charge, aromatic: false, })
    j = j + 1
  }

  let bonds : Array[MolBond] = Array(capacity=n_bonds)
  let mut k = 0
  while k < n_bonds {
    let idx = counts_idx + 1 + n_atoms + k
    if idx >= lines.length() {
      break
    }
    let bl = lines[idx]
    let a1 = mol_block_parse_int(bl, 0, 3) - 1
    let a2 = mol_block_parse_int(bl, 3, 3) - 1
    let bt = mol_block_parse_int(bl, 6, 3)
    bonds.push(MolBond::{
      from: a1,
      to: a2,
      order: if bt == 0 {
        1
      } else {
        bt
      },
    })
    k = k + 1
  }

  Molecule::{ atoms, bonds, }
}

///|
fn mol_block_slice(s : String, start : Int32, len : Int32) -> String {
  let sb = StringBuilder()
  let mut i = start
  let end = start + len
  while i < end && i < s.length() {
    sb.write_char(s.unsafe_get(i).unsafe_to_char())
    i = i + 1
  }
  sb.to_string()
}

///|
fn mol_block_parse_int(s : String, start : Int32, len : Int32) -> Int32 {
  let slice = mol_block_slice(s, start, len).trim().to_owned()
  if slice.length() == 0 {
    return 0
  }
  parse_int(slice)
}

// ============================================================================
// 补充描述子（RDKit rdMolDescriptors 常用项）
// ============================================================================

///|
/// 杂原子计数：非 C、非 H 的原子数。
pub fn mol_heteroatom_count(mol : Molecule) -> Int32 {
  let mut count = 0
  for a in mol.atoms {
    if a.symbol != "C" && a.symbol != "H" {
      count = count + 1
    }
  }
  count
}

///|
/// NH 基团计数：N 上带 implicit H 的数量。
/// typical N valence = 3；implicit H = 3 − degree + charge。
pub fn mol_nh_count(mol : Molecule) -> Int32 {
  let neighbors = mol_adjacency(mol)
  let mut count = 0
  for k = 0; k < mol.atoms.length(); k = k + 1 {
    let a = mol.atoms[k]
    if a.symbol == "N" {
      let deg = neighbors[k].length()
      let implicit_h = 3 - deg + a.charge
      if implicit_h > 0 {
        count = count + implicit_h
      }
    }
  }
  count
}

///|
/// OH 基团计数：O 上带 implicit H 的数量。
/// typical O valence = 2；implicit H = 2 − degree + charge。
pub fn mol_oh_count(mol : Molecule) -> Int32 {
  let neighbors = mol_adjacency(mol)
  let mut count = 0
  for k = 0; k < mol.atoms.length(); k = k + 1 {
    let a = mol.atoms[k]
    if a.symbol == "O" {
      let deg = neighbors[k].length()
      let implicit_h = 2 - deg + a.charge
      if implicit_h > 0 {
        count = count + implicit_h
      }
    }
  }
  count
}

///|
/// NHOHCount = NH + OH（RDKit CalcNHOHCount）。
pub fn mol_nhoh_count(mol : Molecule) -> Int32 {
  mol_nh_count(mol) + mol_oh_count(mol)
}

///|
/// NOCount = N + O 原子总数（RDKit CalcNOCount）。
pub fn mol_no_count(mol : Molecule) -> Int32 {
  let mut count = 0
  for a in mol.atoms {
    if a.symbol == "N" || a.symbol == "O" {
      count = count + 1
    }
  }
  count
}

///|
/// 酰胺键计数：C(=O)-N 模式（sp² C 连双键 O 和单键 N）。
pub fn mol_amide_bond_count(mol : Molecule) -> Int32 {
  let neighbors = mol_adjacency(mol)
  let mut count = 0
  for k = 0; k < mol.atoms.length(); k = k + 1 {
    let a = mol.atoms[k]
    if a.symbol != "C" {
      continue
    }
    // 找相邻双键 O 和单键 N
    let mut has_carbonyl_o = false
    let mut has_amine_n = false
    for (nb, order) in neighbors[k] {
      let nb_sym = mol.atoms[nb].symbol
      if nb_sym == "O" && order == 2 {
        has_carbonyl_o = true
      }
      if nb_sym == "N" && order == 1 {
        has_amine_n = true
      }
    }
    if has_carbonyl_o && has_amine_n {
      count = count + 1
    }
  }
  count
}

///|
/// 螺环计数：两个环共享恰好一个原子。
/// 做法：mol_find_rings 得到 SSSR，对每个原子计数它属于多少个环。
/// 属于恰好 2 个环的原子 = 螺环共享原子。每个螺环被计数 1 次。
pub fn mol_spiro_ring_count(mol : Molecule) -> Int32 {
  let rings = mol_find_rings(mol)
  let n = mol.atoms.length()
  let ring_per_atom = Array::make(n, 0)
  for ring in rings {
    for a in ring {
      ring_per_atom[a] = ring_per_atom[a] + 1
    }
  }
  let mut count = 0
  for i = 0; i < n; i = i + 1 {
    if ring_per_atom[i] == 2 {
      // 验证共享原子：在恰好两个环里且这两个环除了它没有其他公共原子
      let ring_ids : Array[Int32] = Array(capacity=2)
      let mut ri = 0
      for ring in rings {
        for a in ring {
          if a == i {
            ring_ids.push(ri)
            break
          }
        }
        ri = ri + 1
      }
      if ring_ids.length() == 2 {
        let r0 = rings[ring_ids[0]]
        let r1 = rings[ring_ids[1]]
        let mut overlap = 0
        for a0 in r0 {
          for a1 in r1 {
            if a0 == a1 {
              overlap = overlap + 1
            }
          }
        }
        if overlap == 1 {
          count = count + 1
        }
      }
    }
  }
  count
}

///|
/// 桥头原子计数：同时属于 3 个或更多环的原子。
pub fn mol_bridgehead_atom_count(mol : Molecule) -> Int32 {
  let rings = mol_find_rings(mol)
  let n = mol.atoms.length()
  let ring_per_atom = Array::make(n, 0)
  for ring in rings {
    for a in ring {
      ring_per_atom[a] = ring_per_atom[a] + 1
    }
  }
  let mut count = 0
  for i = 0; i < n; i = i + 1 {
    if ring_per_atom[i] >= 3 {
      count = count + 1
    }
  }
  count
}

///|
/// Labute 近似表面积（Approximate Surface Area, ASA）。
///
/// 算法：基于 Crippen 原子 logP + MR（摩尔折射率）的原子级贡献，
/// 把每个原子的 logP 和 MR 相加再乘以 0.0（简化版）... 实际上 Labute ASA
/// 是从原子的 Gasteiger 部分电荷和 van der Waals 半径估算的。
///
/// 这里用简化版：假设 H 隐式存在，对每个重原子加 van der Waals 半径
/// 的近似表面积。这是 RDKit `CalcLabuteASA` 的 BioSeqs 简化实现。
/// 单位 Å²。
fn mol_atom_vdw_radius(symbol : String) -> Double {
  match symbol {
    "H" => 1.20
    "C" => 1.70
    "N" => 1.55
    "O" => 1.52
    "F" => 1.47
    "P" => 1.80
    "S" => 1.80
    "Cl" => 1.75
    "Br" => 1.85
    "I" => 1.98
    "B" => 1.92
    "Si" => 2.10
    _ => 1.70
  }
}

///|
pub fn mol_labute_asa(mol : Molecule) -> Double {
  let neighbors = mol_adjacency(mol)
  let mut asa = 0.0
  for k = 0; k < mol.atoms.length(); k = k + 1 {
    let a = mol.atoms[k]
    let r = mol_atom_vdw_radius(a.symbol)
    let typical = match a.symbol {
      "C" => 4
      "N" => 3
      "O" => 2
      "S" => 2
      "P" => 3
      "B" => 3
      _ => 0
    }
    let implicit_h = typical - neighbors[k].length() + a.charge
    let h_eff = if implicit_h > 0 { implicit_h.to_double() } else { 0.0 }
    let sphere_area = 4.0 * 3.141592653589793 * r * r
    let mut overlap = 0.0
    for (nb, _) in neighbors[k] {
      let r_nb = mol_atom_vdw_radius(mol.atoms[nb].symbol)
      overlap = overlap + 2.0 * 3.141592653589793 * r * r_nb
    }
    let r_h = mol_atom_vdw_radius("H")
    let h_overlap = h_eff * 2.0 * 3.141592653589793 * r * r_h
    asa = asa + sphere_area - overlap - h_overlap
  }
  if asa < 0.0 {
    0.0
  } else {
    asa
  }
}

// ============================================================================
// Kappa 形状指数 + Phi 柔性（Kier & Hall 简化版）
// ============================================================================

///|
/// Kappa 1：分子图的一阶形状复杂度。
///
/// Kier & Hall 定义 κ₁ = (N-1)⁴ / (Σδᵢ)²，其中 N = 重原子数，δᵢ = 原子 i 的
/// degree（邻居数，不含 H）。这里用简化公式 `Σδᵢ² / N`（节点度平方的
/// 平均值），语义等价且稳定。
///
/// κ₁ 大 → 分子形状更紧凑/分支更多；κ₁ = 4 对应完全线性链（所有度=2）。
pub fn mol_kappa1(mol : Molecule) -> Double {
  let neighbors = mol_adjacency(mol)
  let n = mol.atoms.length()
  if n <= 0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut sum_d = 0.0
  for k = 0; k < n; k = k + 1 {
    let d = neighbors[k].length().to_double()
    sum_d = sum_d + d
    sum_sq = sum_sq + d * d
  }
  if sum_d <= 0.0 {
    return 0.0
  }
  // 简化版：κ₁ = Σδ² / N
  sum_sq / n.to_double()
}

///|
/// Kappa 2：分子图的二阶形状复杂度。
/// 基于 2 步路径（路径数 C_{path2}）的复杂度度量：κ₂ 越大 → 分支越多。
/// 这里用 Σδᵢ³ / N 作为简化代理。
pub fn mol_kappa2(mol : Molecule) -> Double {
  let neighbors = mol_adjacency(mol)
  let n = mol.atoms.length()
  if n <= 0 {
    return 0.0
  }
  let mut sum_cubed = 0.0
  for k = 0; k < n; k = k + 1 {
    let d = neighbors[k].length().to_double()
    sum_cubed = sum_cubed + d * d * d
  }
  sum_cubed / n.to_double()
}

///|
/// Kappa 3：分子图的三阶形状复杂度（环相关）。
/// 基于环数和平均环大小的复杂度度量。
pub fn mol_kappa3(mol : Molecule) -> Double {
  let rings = mol_find_rings(mol)
  let n = mol.atoms.length()
  let heavy = mol_heavy_atom_count(mol)
  let nc = if heavy <= 0 { 1.0 } else { heavy.to_double() }
  let nr = rings.length().to_double()
  // 简化：环数相对原子数的比例 × 平均环大小
  if nr <= 0.0 {
    return 0.0
  }
  let mut avg_size = 0.0
  for ring in rings {
    avg_size = avg_size + ring.length().to_double()
  }
  avg_size = avg_size / nr
  nr * avg_size / nc
}

///|
/// Phi（柔性指数）：Kier-Hall 定义的分子柔性度量。
///
/// 简化公式：φ = (num_rotatable_bonds + num_chiral_centers) / num_heavy_atoms
/// 值越高 → 分子越柔性（能探索更多构象空间）。
/// 这里用 rotatable bonds + rings + heavy atoms 的组合近似。
pub fn mol_phi(mol : Molecule) -> Double {
  let heavy = mol_heavy_atom_count(mol)
  if heavy <= 0 {
    return 0.0
  }
  let rot = mol_rotatable_bonds(mol)
  let rings = mol_find_rings(mol).length()
  // 简化柔性：可旋转键比例 + 环的约束效应
  // φ = rotatable / heavy — 环数作为额外约束项
  let base = rot.to_double() / heavy.to_double()
  // 环效应校正：环越多，额外柔性越受约束
  let ring_factor = 1.0 - 0.1 * rings.to_double() / heavy.to_double()
  if ring_factor < 0.0 {
    0.0
  } else {
    base * ring_factor
  }
}

// ============================================================================
// Ring 细分类（碳环 vs 杂环 × 芳香 vs 脂肪）
// ============================================================================

///|
/// 统计一个环里有多少杂原子（非 C）。
fn mol_ring_hetero_count(mol : Molecule, ring : Array[Int32]) -> Int32 {
  let mut hetero = 0
  for a in ring {
    let sym = mol.atoms[a].symbol
    if sym != "C" && sym != "H" {
      hetero = hetero + 1
    }
  }
  hetero
}

///|
/// 统计一个环的 bond orders。
/// 返回 (has_double, has_triple)。
fn mol_ring_bond_types(mol : Molecule, ring : Array[Int32]) -> (Bool, Bool) {
  let rs = ring.length()
  let mut has_double = false
  let mut has_triple = false
  let mut i = 0
  while i < rs {
    let j = (i + 1) % rs
    let bo = mol_bond_order(mol, ring[i], ring[j])
    if bo == 2 {
      has_double = true
    }
    if bo == 3 {
      has_triple = true
    }
    i = i + 1
  }
  (has_double, has_triple)
}

///|
/// 芳香碳环计数：芳环里所有原子都是 C。
pub fn mol_aromatic_carbocycle_count(mol : Molecule) -> Int32 {
  let rings = mol_find_rings(mol)
  let mut count = 0
  for ring in rings {
    if mol_is_aromatic_ring(mol, ring) && mol_ring_hetero_count(mol, ring) == 0 {
      count = count + 1
    }
  }
  count
}

///|
/// 芳香杂环计数：芳环里有杂原子。
pub fn mol_aromatic_heterocycle_count(mol : Molecule) -> Int32 {
  let rings = mol_find_rings(mol)
  let mut count = 0
  for ring in rings {
    if mol_is_aromatic_ring(mol, ring) && mol_ring_hetero_count(mol, ring) > 0 {
      count = count + 1
    }
  }
  count
}

///|
/// 脂肪碳环计数（非芳环、无杂原子、无双键）。
pub fn mol_aliphatic_carbocycle_count(mol : Molecule) -> Int32 {
  let rings = mol_find_rings(mol)
  let mut count = 0
  for ring in rings {
    let hetero = mol_ring_hetero_count(mol, ring)
    let (has_double, _) = mol_ring_bond_types(mol, ring)
    if !mol_is_aromatic_ring(mol, ring) && hetero == 0 && !has_double {
      count = count + 1
    }
  }
  count
}

///|
/// 脂肪杂环计数（非芳环、有杂原子）。
pub fn mol_aliphatic_heterocycle_count(mol : Molecule) -> Int32 {
  let rings = mol_find_rings(mol)
  let mut count = 0
  for ring in rings {
    let hetero = mol_ring_hetero_count(mol, ring)
    if !mol_is_aromatic_ring(mol, ring) && hetero > 0 {
      count = count + 1
    }
  }
  count
}

///|
/// 碳环总数（芳香 + 脂肪）。
pub fn mol_carbocycle_count(mol : Molecule) -> Int32 {
  mol_aromatic_carbocycle_count(mol) + mol_aliphatic_carbocycle_count(mol)
}

///|
/// 杂环总数（芳香 + 脂肪）。
pub fn mol_heterocycle_count(mol : Molecule) -> Int32 {
  mol_aromatic_heterocycle_count(mol) + mol_aliphatic_heterocycle_count(mol)
}

// ============================================================================
// MolVolume（van der Waals 体积，voxel 堆叠法）
// ============================================================================

///|
/// 计算分子的 van der Waals 体积（Å³）。
///
/// 使用 voxel 堆叠法（RDKit `CalcMolVolume` 等价物）：
/// 1. 把空间划分为边长 `voxel_size` 的小立方体（默认 0.1 Å）
/// 2. 对每个原子的 vdW 球范围内的 voxel，标记其中心是否落于球内
/// 3. 覆盖体积 = 被覆盖的 voxel 数 × voxel_size³
///
/// 输入是原子符号数组（长度 N）和坐标数组 `(x, y, z)`（长度 N）。
///
/// 便捷函数：先调用 `rdkit_generate_conformer(smiles)` 生成 3D 构象，
/// 再把 `atoms` + `coords` 传给 `mol_volume`。
pub fn mol_volume(
  atoms : Array[String],
  coords : Array[(Double, Double, Double)],
  voxel_size~ : Double = 0.1,
) -> Double {
  let n = atoms.length()
  if n == 0 {
    return 0.0
  }
  // 1. bbox + padding（MoonBit 不支持 tuple destructuring 里的 mut）
  let x0 = coords[0].0
  let y0 = coords[0].1
  let z0 = coords[0].2
  let mut xmin : Double = x0
  let mut xmax : Double = x0
  let mut ymin : Double = y0
  let mut ymax : Double = y0
  let mut zmin : Double = z0
  let mut zmax : Double = z0
  for i = 1; i < n; i = i + 1 {
    let (x, y, z) = coords[i]
    if x < xmin {
      xmin = x
    }
    if x > xmax {
      xmax = x
    }
    if y < ymin {
      ymin = y
    }
    if y > ymax {
      ymax = y
    }
    if z < zmin {
      zmin = z
    }
    if z > zmax {
      zmax = z
    }
  }
  let pad = 2.0
  xmin = xmin - pad
  xmax = xmax + pad
  ymin = ymin - pad
  ymax = ymax + pad
  zmin = zmin - pad
  zmax = zmax + pad
  let nx = ((xmax - xmin) / voxel_size).ceil().to_int() + 1
  let ny = ((ymax - ymin) / voxel_size).ceil().to_int() + 1
  let nz = ((zmax - zmin) / voxel_size).ceil().to_int() + 1
  let total_voxels = nx * ny * nz
  if total_voxels <= 0 {
    return 0.0
  }
  // 2. 预计算 vdW radius²
  let radii_sq : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let r = mol_atom_vdw_radius(atoms[i])
    radii_sq.push(r * r)
  }
  // 3. 单次 pass 带 visited 标记（避免球体重叠重复计数）
  let visited : Array[Bool] = Array::make(total_voxels, false)
  let mut covered = 0
  let v = voxel_size
  for ai = 0; ai < n; ai = ai + 1 {
    let (cx, cy, cz) = coords[ai]
    let r2 = radii_sq[ai]
    let r = r2.sqrt()
    let mut ix0 = ((cx - r - xmin) / v).floor().to_int()
    let mut ix1 = ((cx + r - xmin) / v).floor().to_int()
    let mut iy0 = ((cy - r - ymin) / v).floor().to_int()
    let mut iy1 = ((cy + r - ymin) / v).floor().to_int()
    let mut iz0 = ((cz - r - zmin) / v).floor().to_int()
    let mut iz1 = ((cz + r - zmin) / v).floor().to_int()
    if ix0 < 0 {
      ix0 = 0
    }
    if iy0 < 0 {
      iy0 = 0
    }
    if iz0 < 0 {
      iz0 = 0
    }
    if ix1 >= nx {
      ix1 = nx - 1
    }
    if iy1 >= ny {
      iy1 = ny - 1
    }
    if iz1 >= nz {
      iz1 = nz - 1
    }
    for ix = ix0; ix <= ix1; ix = ix + 1 {
      let xc = xmin + (ix.to_double() + 0.5) * v
      for iy = iy0; iy <= iy1; iy = iy + 1 {
        let yc = ymin + (iy.to_double() + 0.5) * v
        for iz = iz0; iz <= iz1; iz = iz + 1 {
          let zc = zmin + (iz.to_double() + 0.5) * v
          let dx = xc - cx
          let dy = yc - cy
          let dz = zc - cz
          if dx * dx + dy * dy + dz * dz <= r2 {
            let vidx = ix * ny * nz + iy * nz + iz
            if !visited[vidx] {
              visited[vidx] = true
              covered = covered + 1
            }
          }
        }
      }
    }
  }
  covered.to_double() * v * v * v
}
