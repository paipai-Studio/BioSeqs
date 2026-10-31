///|
/// MHC 多肽结合预测(锚定位置规则框架).
///
/// 定位: 教学与新抗原候选**筛查**级工具. 精度显著低于 netMHCpan/NetMHCpan-4.1
/// 等神经网络方法(其权重未公开, 无法移植); 本模块采用文献教科书级的
/// 经典锚定规则(等位基因特异性 P2/P9 等"锚"位置残基偏好)与半定量打分,
/// 用于缩小候选集后再交给实验或成熟工具验证.
///
/// MHC-I: 呈递 8-11 肽(以 9 肽为主), 锚定规则作用于肽段固定位置;
/// MHC-II: 结合槽开放, 呈递 13-25 肽, 内部 9-mer 核心框架决定结合,
/// 本实现扫描全部 7 个框架取最优.
pub struct MhcAllele {
  name : String
  class : Int32
  anchors : Array[(Int32, String)]
  preferred : Array[(Int32, String)]
  peptide_lengths : Array[Int32]
}

///|
/// 构造等位基因定义. anchors/preferred 为 (1-based 位置, 优选残基串).
pub fn MhcAllele::new(
  name : String,
  class~ : Int32 = 1,
  anchors : Array[(Int32, String)],
  preferred~ : Array[(Int32, String)] = [],
  peptide_lengths~ : Array[Int32] = [9],
) -> MhcAllele {
  { name, class, anchors, preferred, peptide_lengths, }
}

///|
/// 等位基因名称.
pub fn MhcAllele::name(self : MhcAllele) -> String {
  self.name
}

///|
/// MHC 类别(1 或 2).
pub fn MhcAllele::mhc_class(self : MhcAllele) -> Int32 {
  self.class
}

///|
/// 经典等位基因库(锚定规则为教科书级共识, 见各等位基因注释).
pub fn mhc_allele(name : String) -> MhcAllele? {
  let n = name.to_upper()
  if n == "HLA-A*02:01" || n == "A*02:01" || n == "HLA-A0201" {
    // 经典超型原型: P2 容氢phobic/小侧链锚, P9 疏水锚, P6 极性小侧链偏好
    return Some(
      MhcAllele::new(
        "HLA-A*02:01",
        [(2, "ILMVA"), (9, "VLIM")],
        preferred=[(1, "MSAT"), (3, "ILMV"), (6, "NSTV"), (7, "KP")],
        peptide_lengths=[9, 10],
      ),
    )
  }
  if n == "HLA-A*01:01" || n == "A*01:01" {
    return Some(
      MhcAllele::new("HLA-A*01:01", [(2, "TSM"), (9, "YL")], preferred=[
        (3, "DEKRS"),
        (6, "IVL"),
        (7, "PSAT"),
      ]),
    )
  }
  if n == "HLA-A*03:01" || n == "A*03:01" {
    return Some(
      MhcAllele::new("HLA-A*03:01", [(2, "VMILT"), (9, "KY")], preferred=[
        (3, "RK"),
        (7, "KR"),
      ]),
    )
  }
  if n == "HLA-B*07:02" || n == "B*07:02" {
    // B07 超型: 脯氨酸 P2 锚为标志性特征
    return Some(
      MhcAllele::new("HLA-B*07:02", [(2, "P"), (9, "LIMFY")], preferred=[
        (3, "DEKRS"),
        (5, "RS"),
      ]),
    )
  }
  if n == "HLA-B*44:02" || n == "B*44:02" {
    return Some(
      MhcAllele::new("HLA-B*44:02", [(2, "DE"), (9, "FYIL")], preferred=[
        (3, "RK"),
        (5, "RK"),
      ]),
    )
  }
  if n == "H-2KB" || n == "H2KB" {
    // 小鼠 H-2Kb: 经典 F/Y P5 与疏水 P8 锚(8 肽) / P5,P9(9 肽)
    return Some(
      MhcAllele::new(
        "H-2Kb",
        [(5, "FY"), (8, "LIVM")],
        preferred=[(2, "RK"), (3, "ILVY")],
        peptide_lengths=[8, 9],
      ),
    )
  }
  if n == "H-2DB" || n == "H2DB" {
    return Some(
      MhcAllele::new(
        "H-2Db",
        [(5, "N"), (8, "ILVM")],
        preferred=[(2, "RK"), (3, "ILVY")],
        peptide_lengths=[8, 9],
      ),
    )
  }
  if n == "HLA-DRB1*01:01" || n == "DRB1*01:01" {
    // II 类经典: 9-mer 核心框架 P1/P4/P6/P9 锚
    return Some(
      MhcAllele::new(
        "HLA-DRB1*01:01",
        class=2,
        [(1, "FYWILMV"), (4, "LIVMA"), (6, "STPA"), (9, "VLIMA")],
        preferred=[(2, "RK"), (3, "FYWL"), (7, "QKRAE")],
        peptide_lengths=[15],
      ),
    )
  }
  if n == "HLA-DRB1*04:01" || n == "DRB1*04:01" {
    return Some(
      MhcAllele::new(
        "HLA-DRB1*04:01",
        class=2,
        [(1, "FYWLMV"), (4, "LIVMA"), (6, "STPA"), (9, "VLIMA")],
        preferred=[(2, "KR"), (7, "QKRAE")],
        peptide_lengths=[15],
      ),
    )
  }
  None
}

///|
/// 全部内置等位基因名称.
pub fn mhc_allele_names() -> Array[String] {
  [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-B*07:02", "HLA-B*44:02", "H-2Kb",
    "H-2Db", "HLA-DRB1*01:01", "HLA-DRB1*04:01",
  ]
}

///|
/// 一条结合预测命中.
pub struct MhcPeptideHit {
  peptide : String
  start : Int32
  allele : String
  score : Double
  anchor_matches : Int32
  anchor_total : Int32
  peptide_len : Int32
  core_start : Int32
}

///|
/// 命中肽段(核心框架, II 类为 9-mer 核心).
pub fn MhcPeptideHit::peptide(self : MhcPeptideHit) -> String {
  self.peptide
}

///|
/// 蛋白内起始位置(1-based).
pub fn MhcPeptideHit::start(self : MhcPeptideHit) -> Int32 {
  self.start
}

///|
/// 等位基因名.
pub fn MhcPeptideHit::allele(self : MhcPeptideHit) -> String {
  self.allele
}

///|
/// 打分 [0,1]: 锚定匹配占 70%, 偏好位匹配占 30%.
pub fn MhcPeptideHit::score(self : MhcPeptideHit) -> Double {
  self.score
}

///|
/// 匹配的锚数 / 总锚数.
pub fn MhcPeptideHit::anchor_matches(self : MhcPeptideHit) -> Int32 {
  self.anchor_matches
}

///|
/// II 类命中在原始肽段内的核心框架起点(0-based; I 类为 0).
pub fn MhcPeptideHit::core_start(self : MhcPeptideHit) -> Int32 {
  self.core_start
}

///|
/// 强结合子: 全锚匹配且得分 >= 0.8.
pub fn MhcPeptideHit::is_strong(self : MhcPeptideHit) -> Bool {
  self.anchor_matches == self.anchor_total && self.score >= 0.8
}

///|
/// 弱结合子: 全锚匹配(得分 >= 0.6).
pub fn MhcPeptideHit::is_weak(self : MhcPeptideHit) -> Bool {
  self.anchor_matches == self.anchor_total && self.score >= 0.6
}

///|
/// 分级标签("strong"/"weak"/"none").
pub fn MhcPeptideHit::binding_category(self : MhcPeptideHit) -> String {
  if self.is_strong() {
    "strong"
  } else if self.is_weak() {
    "weak"
  } else {
    "none"
  }
}

///|
/// 对单条肽段按等位基因规则打分.
/// 返回 (score, anchor_matches, anchor_total); 位置超界抛 SeqError.
pub fn mhc_score_peptide(
  peptide : String,
  allele : MhcAllele,
) -> (Double, Int32, Int32) raise SeqError {
  let n = peptide.length()
  if n < 3 {
    raise SeqError("peptide too short for MHC scoring")
  }
  let mut anchor_hit = 0
  let mut anchor_total = 0
  for a in allele.anchors {
    let (pos, residues) = a
    if pos > n {
      continue
    }
    anchor_total = anchor_total + 1
    let c = peptide.unsafe_get(pos - 1)
    if residues.contains_code_unit(c) {
      anchor_hit = anchor_hit + 1
    }
  }
  let mut pref_hit = 0
  let mut pref_total = 0
  for p in allele.preferred {
    let (pos, residues) = p
    if pos > n {
      continue
    }
    pref_total = pref_total + 1
    let c = peptide.unsafe_get(pos - 1)
    if residues.contains_code_unit(c) {
      pref_hit = pref_hit + 1
    }
  }
  let anchor_score = if anchor_total > 0 {
    anchor_hit.to_double() / anchor_total.to_double()
  } else {
    0.0
  }
  let pref_score = if pref_total > 0 {
    pref_hit.to_double() / pref_total.to_double()
  } else {
    0.0
  }
  (0.7 * anchor_score + 0.3 * pref_score, anchor_hit, anchor_total)
}

///|
/// MHC-I 预测: 蛋白序列滑窗提取肽段并打分, 按得分降序返回.
/// peptide_len 为提取长度(默认取等位基因支持的第一个长度).
pub fn mhc_predict_i(
  protein : String,
  allele : MhcAllele,
  peptide_len~ : Int32 = 0,
) -> Array[MhcPeptideHit] raise SeqError {
  if allele.class != 1 {
    raise SeqError("mhc_predict_i expects a class I allele, got " + allele.name)
  }
  let len = if peptide_len > 0 {
    peptide_len
  } else {
    allele.peptide_lengths[0]
  }
  let p = protein.to_upper()
  let n = p.length()
  let hits : Array[MhcPeptideHit] = Array(capacity=0)
  if n < len {
    return hits
  }
  for start in 0..<=(n - len) {
    let pep = p[start:start + len].to_owned()
    let (score, am, at) = mhc_score_peptide(pep, allele)
    hits.push({
      peptide: pep,
      start: start + 1,
      allele: allele.name,
      score,
      anchor_matches: am,
      anchor_total: at,
      peptide_len: len,
      core_start: 0,
    })
  }
  // 得分降序
  let mut i = 0
  while i < hits.length() {
    let mut j = i + 1
    while j < hits.length() {
      if hits[j].score > hits[i].score {
        let t = hits[i]
        hits[i] = hits[j]
        hits[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  hits
}

///|
/// MHC-II 预测: 15-mer 肽段内扫描全部 7 个 9-mer 核心框架,
/// 取锚定得分最高的框架作为该肽段结果, 蛋白滑窗步长 1,
/// 按得分降序返回(相邻窗口会覆盖同一强表位, 可用 unique_ 去重).
pub fn mhc_predict_ii(
  protein : String,
  allele : MhcAllele,
  peptide_len~ : Int32 = 15,
) -> Array[MhcPeptideHit] raise SeqError {
  if allele.class != 2 {
    raise SeqError(
      "mhc_predict_ii expects a class II allele, got " + allele.name,
    )
  }
  let p = protein.to_upper()
  let n = p.length()
  let hits : Array[MhcPeptideHit] = Array(capacity=0)
  if n < peptide_len {
    return hits
  }
  for start in 0..<=(n - peptide_len) {
    let pep = p[start:start + peptide_len].to_owned()
    let mut best_score = -1.0
    let mut best_am = 0
    let mut best_at = 0
    let mut best_frame = 0
    for frame in 0..<=(peptide_len - 9) {
      let core = pep[frame:frame + 9].to_owned()
      let (score, am, at) = mhc_score_peptide(core, allele)
      if score > best_score {
        best_score = score
        best_am = am
        best_at = at
        best_frame = frame
      }
    }
    hits.push({
      peptide: pep,
      start: start + 1,
      allele: allele.name,
      score: best_score,
      anchor_matches: best_am,
      anchor_total: best_at,
      peptide_len,
      core_start: best_frame,
    })
  }
  let mut i = 0
  while i < hits.length() {
    let mut j = i + 1
    while j < hits.length() {
      if hits[j].score > hits[i].score {
        let t = hits[i]
        hits[i] = hits[j]
        hits[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  hits
}

///|
/// 统一入口: 按等位基因类别分发 I/II 类预测.
pub fn mhc_predict(
  protein : String,
  allele : MhcAllele,
  peptide_len~ : Int32 = 0,
) -> Array[MhcPeptideHit] raise SeqError {
  if allele.class == 2 {
    if peptide_len > 0 {
      mhc_predict_ii(protein, allele, peptide_len~)
    } else {
      mhc_predict_ii(protein, allele)
    }
  } else if peptide_len > 0 {
    mhc_predict_i(protein, allele, peptide_len~)
  } else {
    mhc_predict_i(protein, allele)
  }
}

///|
/// 只保留弱结合子及以上的命中(过滤 none).
pub fn mhc_binders(
  hits : Array[MhcPeptideHit],
  strong_only~ : Bool = false,
) -> Array[MhcPeptideHit] {
  let out : Array[MhcPeptideHit] = Array(capacity=0)
  for h in hits {
    if strong_only {
      if h.is_strong() {
        out.push(h)
      }
    } else if h.is_weak() {
      out.push(h)
    }
  }
  out
}

///|
/// 按肽段(II 类为核心 9-mer)去重, 保留得分最高的命中.
pub fn mhc_unique_epitopes(hits : Array[MhcPeptideHit]) -> Array[MhcPeptideHit] {
  let best : Map[String, MhcPeptideHit] = Map([], capacity=0)
  let order : Array[String] = Array(capacity=0)
  for h in hits {
    // II 类以核心 9-mer 为表位身份; I 类肽段本身就是核心
    let key = if h.peptide_len > 11 {
      h.peptide[h.core_start:h.core_start + 9].to_owned()
    } else {
      h.peptide
    }
    match best.get(key) {
      Some(prev) => if h.score > prev.score { best[key] = h }
      None => {
        best[key] = h
        order.push(key)
      }
    }
  }
  let out : Array[MhcPeptideHit] = Array(capacity=order.length())
  for k in order {
    out.push(best[k])
  }
  // 得分降序
  let mut i = 0
  while i < out.length() {
    let mut j = i + 1
    while j < out.length() {
      if out[j].score > out[i].score {
        let t = out[i]
        out[i] = out[j]
        out[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  out
}
