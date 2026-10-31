///|
/// foldx_lite.mbt — FoldX 风格蛋白质稳定性 ΔΔG 与饱和突变扫描
/// (本地确定性实现, Cα 统计势能量差法).
///
/// 定位与边界(诚实声明): FoldX/Rosetta 的力场(侧链 rotamer/
/// 溶剂化/氢键网络)未公开完整且依赖全原子重打包, 本模块采用
/// **Cα-only 统计势能量差**作为筛查级 ΔΔG 代理:
/// - 突变保持 Cα 骨架几何不变, 仅替换残基类型;
/// - ΔΔG = E(mut) - E(wt), 能量为 DOPE 统计势
///   (复用 dope_score.mbt, 210 氨基酸对 x 距离 bin);
/// - ΔΔG < 0 倾向更稳定, > 0 倾向去稳定.
///
/// 功能:
/// - 单点突变 ΔΔG(20 aa 全替换或指定);
/// - 饱和突变扫描: 每位点 x 20 aa 的 ΔΔG 矩阵 + 每位点
///   最稳/最不稳突变 + Boltzmann 适合度(exp(-ddG/RT))与
///   位点 Shannon 熵(可突变性);
/// - 突变文库生成与稳定性排名;
/// - TSV/JSON 导出.
///
/// 数据不出本机.

///|
/// 20 种标准氨基酸三字母码(按字母序).
pub fn foldx_aa3_list() -> Array[String] {
  [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  ]
}

///|
/// 三字母 -> 单字母(20 标准 aa, 其余 '?').
pub fn foldx_aa3_to_aa1(aa3 : String) -> String {
  match aa3 {
    "ALA" => "A"
    "ARG" => "R"
    "ASN" => "N"
    "ASP" => "D"
    "CYS" => "C"
    "GLN" => "Q"
    "GLU" => "E"
    "GLY" => "G"
    "HIS" => "H"
    "ILE" => "I"
    "LEU" => "L"
    "LYS" => "K"
    "MET" => "M"
    "PHE" => "F"
    "PRO" => "P"
    "SER" => "S"
    "THR" => "T"
    "TRP" => "W"
    "TYR" => "Y"
    "VAL" => "V"
    _ => "?"
  }
}

///|
/// 一条单点突变记录.
pub(all) struct FoldxMutation {
  /// 1-based 残基序号(链内).
  pos : Int32
  wt : String
  /// 突变后三字母码.
  mt : String
  /// E(mut) - E(wt), <0 更稳定.
  ddg : Double
  /// Boltzmann 适合度 exp(-ddg / RT), RT=0.593 kcal/mol(298K).
  fitness : Double
}

///|
/// 单个位点的饱和突变扫描结果.
pub(all) struct FoldxScanSite {
  pos : Int32
  wt : String
  /// 20 种替换的 ΔΔG(mt 三字母 -> ΔΔG), 不含 wt 自身.
  ddg : Array[FoldxMutation]
  /// ΔΔG 最小(最稳)的突变.
  best : FoldxMutation
  /// ΔΔG 最大(最不稳)的突变.
  worst : FoldxMutation
  /// 位点 Shannon 熵(基于 fitness 归一分布, 可突变性度量).
  entropy : Double
}

///|
/// 扫描汇总.
pub(all) struct FoldxScanSummary {
  n_sites : Int32
  /// 全部位点最稳突变(整体最优).
  global_best : FoldxMutation
  /// 全部位点最不稳突变.
  global_worst : FoldxMutation
  /// 平均位点熵.
  mean_entropy : Double
}

// ---------------------------------------------------------------------------
// 能量与突变
// ---------------------------------------------------------------------------

///|
/// RT(kcal/mol, 298K).
let foldx_rt : Double = 0.593

///|
/// 链的 DOPE 总能量(残基数 < 2 返回 None).
pub fn foldx_energy(chain : Chain) -> Double? {
  let pot = DopePotential::default()
  match pot.score(chain) {
    Some(s) => Some(s.get_total())
    None => None
  }
}

///|
/// 构造突变链: 替换 pos(1-based)残基类型, Cα 骨架几何保持.
fn foldx_mutate_chain(chain : Chain, pos : Int32, mt : String) -> Chain? {
  let residues = chain.get_residues()
  let idx = pos - 1
  if idx < 0 || idx >= residues.length() {
    return None
  }
  let old = residues[idx]
  // 保留 Cα 原子(统计势仅用 Cα)
  let kept_atoms : Array[Atom] = Array(capacity=0)
  for atom in old.get_atoms() {
    if atom.name == "CA" {
      kept_atoms.push(atom)
    }
  }
  let mutated = Residue::new(
    resname=mt,
    chainid=old.chainid,
    resseq=old.resseq,
    atoms=kept_atoms,
    hetfield=old.hetfield,
  )
  let new_residues : Array[Residue] = Array(capacity=residues.length())
  for i = 0; i < residues.length(); i = i + 1 {
    if i == idx {
      new_residues.push(mutated)
    } else {
      new_residues.push(residues[i])
    }
  }
  Some(Chain::new(id=chain.id, residues=new_residues))
}

///|
/// 单点突变 ΔΔG(E(mut) - E(wt)).
/// pos 为 1-based 链内残基序号; mt 为三字母码.
pub fn foldx_ddg(chain : Chain, pos : Int32, mt : String) -> Double? {
  match foldx_energy(chain) {
    None => None
    Some(wt_e) =>
      match foldx_mutate_chain(chain, pos, mt) {
        None => None
        Some(mut_chain) =>
          match foldx_energy(mut_chain) {
            None => None
            Some(mut_e) => Some(mut_e - wt_e)
          }
      }
  }
}

// ---------------------------------------------------------------------------
// 饱和突变扫描
// ---------------------------------------------------------------------------

///|
fn foldx_make_mutation(
  pos : Int32,
  wt : String,
  mt : String,
  ddg : Double,
) -> FoldxMutation {
  let fitness = @math.exp(-ddg / foldx_rt)
  { pos, wt, mt, ddg, fitness, }
}

///|
/// 饱和突变扫描: 指定位点集(1-based, 默认全部)对 20 种氨基酸替换
/// 逐一计算 ΔΔG; 每位点含最稳/最不稳突变与 Shannon 熵
/// (fitness 归一分布, 度量位点可突变性).
pub fn foldx_saturation_scan(
  chain : Chain,
  positions~ : Array[Int32]? = None,
) -> Array[FoldxScanSite] {
  let residues = chain.get_residues()
  let n = residues.length()
  let sites : Array[Int32] = match positions {
    Some(ps) => {
      let filtered : Array[Int32] = Array(capacity=0)
      for p in ps {
        if p >= 1 && p <= n {
          filtered.push(p)
        }
      }
      filtered.sort()
      filtered
    }
    None => {
      let all : Array[Int32] = Array(capacity=n)
      for i = 1; i <= n; i = i + 1 {
        all.push(i)
      }
      all
    }
  }
  let wt_energy = foldx_energy(chain)
  let wt_e = match wt_energy {
    Some(e) => e
    None => return []
  }
  let out : Array[FoldxScanSite] = Array(capacity=sites.length())
  for p in sites {
    let wt = residues[p - 1].resname
    let muts : Array[FoldxMutation] = Array(capacity=19)
    let mut best : FoldxMutation? = None
    let mut worst : FoldxMutation? = None
    for mt in foldx_aa3_list() {
      if mt == wt {
        continue
      }
      let ddg : Double = match foldx_mutate_chain(chain, p, mt) {
        Some(mc) =>
          match foldx_energy(mc) {
            Some(me) => me - wt_e
            None => continue
          }
        None => continue
      }
      let m = foldx_make_mutation(p, wt, mt, ddg)
      muts.push(m)
      match best {
        Some(b) => if ddg < b.ddg { best = Some(m) }
        None => best = Some(m)
      }
      match worst {
        Some(w) => if ddg > w.ddg { worst = Some(m) }
        None => worst = Some(m)
      }
    }
    match (best, worst) {
      (Some(b), Some(w)) => {
        // Shannon 熵: fitness 归一分布(含 wt fitness=1)
        let mut total_fit = 1.0
        for m in muts {
          total_fit = total_fit + m.fitness
        }
        let mut entropy = 0.0
        let p_wt = 1.0 / total_fit
        if p_wt > 0.0 {
          entropy = entropy - p_wt * @math.ln(p_wt)
        }
        for m in muts {
          let pm = m.fitness / total_fit
          if pm > 0.0 {
            entropy = entropy - pm * @math.ln(pm)
          }
        }
        out.push(FoldxScanSite::{
          pos: p,
          wt,
          ddg: muts,
          best: b,
          worst: w,
          entropy,
        })
      }
      _ => ()
    }
  }
  out
}

///|
/// 扫描汇总(全局最稳/最不稳与平均熵).
pub fn foldx_scan_summary(sites : Array[FoldxScanSite]) -> FoldxScanSummary? {
  if sites.length() == 0 {
    return None
  }
  let mut global_best : FoldxMutation? = None
  let mut global_worst : FoldxMutation? = None
  let mut entropy_sum = 0.0
  for s in sites {
    entropy_sum = entropy_sum + s.entropy
    match global_best {
      Some(b) => if s.best.ddg < b.ddg { global_best = Some(s.best) }
      None => global_best = Some(s.best)
    }
    match global_worst {
      Some(w) => if s.worst.ddg > w.ddg { global_worst = Some(s.worst) }
      None => global_worst = Some(s.worst)
    }
  }
  match (global_best, global_worst) {
    (Some(b), Some(w)) =>
      Some(FoldxScanSummary::{
        n_sites: sites.length(),
        global_best: b,
        global_worst: w,
        mean_entropy: entropy_sum / sites.length().to_double(),
      })
    _ => None
  }
}

// ---------------------------------------------------------------------------
// 突变文库
// ---------------------------------------------------------------------------

///|
/// 文库中的一条突变体.
pub(all) struct FoldxMutant {
  /// "A12G" 风格标签(wt 单字母 + pos + mt 单字母, 多位点逗号分隔).
  label : String
  mutations : Array[FoldxMutation]
  /// 各位点 ΔΔG 之和(加性近似).
  total_ddg : Double
}

///|
/// 生成单点突变文库并按稳定性排名(ΔΔG 升序 = 最稳在前).
/// sites 来自饱和扫描结果.
pub fn foldx_mutant_library(sites : Array[FoldxScanSite]) -> Array[FoldxMutant] {
  let out : Array[FoldxMutant] = Array(capacity=0)
  for s in sites {
    for m in s.ddg {
      let label = foldx_aa3_to_aa1(s.wt) +
        s.pos.to_string() +
        foldx_aa3_to_aa1(m.mt)
      out.push(FoldxMutant::{ label, mutations: [m], total_ddg: m.ddg, })
    }
  }
  out.sort_by(fn(a, b) {
    if a.total_ddg < b.total_ddg {
      -1
    } else if a.total_ddg > b.total_ddg {
      1
    } else {
      0
    }
  })
  out
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 突变记录行 "pos wt mt ddg fitness".
fn foldx_mutation_tsv_row(m : FoldxMutation) -> String {
  "\{m.pos}\t\{m.wt}\t\{m.mt}\t\{m.ddg}\t\{m.fitness}"
}

///|
/// 扫描结果 -> TSV 文本(每替换一行).
pub fn foldx_scan_tsv(sites : Array[FoldxScanSite]) -> String {
  let sb = StringBuilder()
  sb.write_string("pos\twt\tmt\tddg\tfitness\n")
  for s in sites {
    for m in s.ddg {
      sb.write_string(foldx_mutation_tsv_row(m))
      sb.write_char('\n')
    }
  }
  sb.to_string()
}

///|
/// 扫描结果 -> JSON 值.
pub fn foldx_scan_value(sites : Array[FoldxScanSite]) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=sites.length())
  for s in sites {
    let muts : Array[JsonValue] = Array(capacity=s.ddg.length())
    for m in s.ddg {
      muts.push(
        JsonObj(
          Map([
            ("mt", JsonStr(m.mt)),
            ("ddg", JsonNum(m.ddg)),
            ("fitness", JsonNum(m.fitness)),
          ]),
        ),
      )
    }
    arr.push(
      JsonObj(
        Map([
          ("pos", JsonNum(s.pos.to_double())),
          ("wt", JsonStr(s.wt)),
          ("best", JsonStr(s.best.mt)),
          ("best_ddg", JsonNum(s.best.ddg)),
          ("worst", JsonStr(s.worst.mt)),
          ("worst_ddg", JsonNum(s.worst.ddg)),
          ("entropy", JsonNum(s.entropy)),
          ("mutations", JsonArr(muts)),
        ]),
      ),
    )
  }
  JsonArr(arr)
}

///|
/// 扫描结果 -> JSON 文本.
pub fn foldx_scan_json(sites : Array[FoldxScanSite]) -> String {
  json_to_string(foldx_scan_value(sites))
}

///|
/// 汇总摘要行.
pub fn foldx_summary_line(summary : FoldxScanSummary) -> String {
  "foldx-lite: sites=\{summary.n_sites}, best=\{foldx_aa3_to_aa1(summary.global_best.wt)}\{summary.global_best.pos}\{foldx_aa3_to_aa1(summary.global_best.mt)}(\{summary.global_best.ddg}), worst=\{foldx_aa3_to_aa1(summary.global_worst.wt)}\{summary.global_worst.pos}\{foldx_aa3_to_aa1(summary.global_worst.mt)}(\{summary.global_worst.ddg}), mean_entropy=\{summary.mean_entropy}"
}
