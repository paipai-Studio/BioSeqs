///|
/// fba_lite.mbt — 代谢通量平衡分析(cobra-lite, 单纯形 LP 求解).
///
/// 对齐 cobrapy 的核心 FBA 语义(本地实现):
/// - 化学计量模型: metabolites x reactions 矩阵 S, 反应上下界
/// - FBA: 稳态约束 S·v = 0 下最大化 biomass 目标函数
/// - LP 求解: 两阶段单纯形法(Bland 规则防循环, 确定性)
/// - 变量变换: v = base + x+ - x-(base=max(lb,0), 双向非负拆分)
/// - 基因敲除: 反应 bounds 置零重解, 输出表型(生长/非生长)
/// - 模型 IO: JSON 交换格式(契合"内存数据模型+JSON"架构)
///
/// 纯算法无生态依赖. 数据不出本机.
pub suberror FbaError {
  FbaError(String)
}

///|
/// 代谢反应.
pub(all) struct FbaReaction {
  id : String
  name : String
  lower_bound : Double
  upper_bound : Double
  /// 目标系数(biomass 反应为 1)
  objective : Double
  /// 关联基因(GPR 简化: 敲除集合)
  genes : Array[String]
}

///|
/// 代谢模型.
pub(all) struct FbaModel {
  metabolites : Array[String]
  reactions : Array[FbaReaction]
  /// S[metabolite][reaction] 化学计量矩阵
  stoichiometry : Array[Array[Double]]
}

///|
/// FBA 求解结果.
pub(all) struct FbaSolution {
  objective_value : Double
  /// 每反应通量(与模型反应同序)
  fluxes : Array[Double]
  status : String
}

// ---------------------------------------------------------------------------
// 单纯形法(两阶段, 标准型 max c^T x, Ax=b, x>=0)
// ---------------------------------------------------------------------------

///|
priv struct SimplexTableau {
  /// 约束矩阵(含 slack/artificial), m x n_total
  a : Array[Array[Double]]
  /// rhs
  b : Array[Double]
  /// 基变量索引(m 个)
  basis : Array[Int32]
  /// 总变量数(n 原始 + slack + artificial)
  n_total : Int32
}

///|
fn simplex_pivot(t : SimplexTableau, row : Int32, col : Int32) -> Unit {
  let m = t.basis.length()
  let piv = t.a[row][col]
  // 归一化主元行
  for j = 0; j < t.n_total; j = j + 1 {
    t.a[row][j] = t.a[row][j] / piv
  }
  t.b[row] = t.b[row] / piv
  // 消元其它行
  for i = 0; i < m; i = i + 1 {
    if i == row {
      continue
    }
    let factor = t.a[i][col]
    if factor == 0.0 {
      continue
    }
    for j = 0; j < t.n_total; j = j + 1 {
      t.a[i][j] = t.a[i][j] - factor * t.a[row][j]
    }
    t.b[i] = t.b[i] - factor * t.b[row]
  }
  t.basis[row] = col
}

///|
/// 单纯形主循环(Bland 规则): max c^T x.
/// 返回 (最优值, x 解); 不可行/无界返回 None.
fn simplex_solve(
  t : SimplexTableau,
  c : Array[Double],
  artificial_cols : Array[Int32],
) -> (Double, Array[Double])? {
  let m = t.basis.length()
  // 目标行: z = c^T x -> reduced cost r_j = c_j - c_B^T B^-1 A_j
  // 表格法: 维护 z 行等价于每轮计算 reduced cost
  let iter_max = 2000
  let mut iter = 0
  while iter < iter_max {
    iter = iter + 1
    // 基变量的 c 值
    let cb : Array[Double] = Array(capacity=0)
    for i = 0; i < m; i = i + 1 {
      let bv = t.basis[i]
      cb.push(if bv < c.length() { c[bv] } else { 0.0 })
    }
    // reduced cost(Bland: 最小下标正者)
    let mut enter : Int32 = -1
    for j = 0; j < t.n_total; j = j + 1 {
      // 人工变量不入基(第二阶段)
      let is_art = artificial_cols.contains(j)
      if is_art {
        continue
      }
      // 检查是否在基中
      let mut in_basis = false
      for i = 0; i < m; i = i + 1 {
        if t.basis[i] == j {
          in_basis = true
          break
        }
      }
      if in_basis {
        continue
      }
      let mut rj = if j < c.length() { c[j] } else { 0.0 }
      for i = 0; i < m; i = i + 1 {
        rj = rj - cb[i] * t.a[i][j]
      }
      if rj > 1.0e-9 {
        enter = j
        break
      }
    }
    if enter < 0 {
      // 最优
      let mut obj = 0.0
      for i = 0; i < m; i = i + 1 {
        let bv = t.basis[i]
        if bv < c.length() {
          obj = obj + c[bv] * t.b[i]
        }
      }
      let x : Array[Double] = Array(capacity=c.length())
      for j = 0; j < c.length(); j = j + 1 {
        x.push(0.0)
      }
      for i = 0; i < m; i = i + 1 {
        let bv = t.basis[i]
        if bv < c.length() {
          x[bv] = t.b[i]
        }
      }
      return Some((obj, x))
    }
    // 比值检验(最小非负 rhs/主元列; 平局取最小行索引=Bland)
    let mut leave : Int32 = -1
    let mut best_ratio = 1.0e18
    for i = 0; i < m; i = i + 1 {
      let aij = t.a[i][enter]
      if aij > 1.0e-9 {
        let ratio = t.b[i] / aij
        if ratio >= -1.0e-9 && ratio < best_ratio - 1.0e-12 {
          best_ratio = if ratio < 0.0 { 0.0 } else { ratio }
          leave = i
        }
      }
    }
    if leave < 0 {
      // 无界
      return None
    }
    simplex_pivot(t, leave, enter)
  }
  None
}

///|
/// 通用 LP 求解: max c^T x, s.t. A_eq x = b_eq, 0 <= x(分量上界可选).
/// 两阶段: 人工变量求初始基 -> 原目标求解.
pub fn fba_solve_lp(
  a_eq : Array[Array[Double]],
  b_eq : Array[Double],
  c : Array[Double],
  ub : Array[Double],
) -> (Double, Array[Double])? {
  let m = a_eq.length()
  if m == 0 {
    // 无约束: 正系数且上界无穷 -> 无界
    for j = 0; j < c.length(); j = j + 1 {
      if c[j] > 0.0 && ub[j] >= 1.0e15 {
        return None
      }
    }
    let x : Array[Double] = Array(capacity=0)
    let mut obj = 0.0
    for j = 0; j < c.length(); j = j + 1 {
      let v = if c[j] > 0.0 { ub[j] } else { 0.0 }
      x.push(v)
      obj = obj + c[j] * v
    }
    return Some((obj, x))
  }
  let n = c.length()
  // 变量: n 原始 + n 有上界的 slack 上界行 + m 人工
  // 上界 x_j <= ub_j 转为附加行: x_j + s_j = ub_j
  let bounded : Array[Int32] = Array(capacity=0)
  for j = 0; j < n; j = j + 1 {
    if ub[j] < 1.0e15 {
      bounded.push(j)
    }
  }
  let n_slack = bounded.length()
  let n_rows = m + n_slack
  let n_total : Int32 = n + n_slack + n_rows
  // 构造约束矩阵: m 原始行 + n_slack 上界行
  let rows = m + n_slack
  let a : Array[Array[Double]] = Array(capacity=0)
  let b : Array[Double] = Array(capacity=0)
  for i = 0; i < rows; i = i + 1 {
    a.push(Array::make(n_total, 0.0))
  }
  for i = 0; i < m; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      a[i][j] = a_eq[i][j]
    }
    b.push(b_eq[i])
  }
  for k = 0; k < n_slack; k = k + 1 {
    let j = bounded[k]
    a[m + k][j] = 1.0
    a[m + k][n + k] = 1.0
    b.push(ub[j])
  }
  // rhs 负行: 人工变量处理(对 b_i < 0 的行先取负)
  for i = 0; i < rows; i = i + 1 {
    if b[i] < 0.0 {
      for j = 0; j < n_total; j = j + 1 {
        a[i][j] = -a[i][j]
      }
      b[i] = -b[i]
    }
  }
  // 人工变量(所有行, 单位列)
  for i = 0; i < rows; i = i + 1 {
    a[i][n + n_slack + i] = 1.0
  }
  let basis : Array[Int32] = Array(capacity=0)
  for i = 0; i < rows; i = i + 1 {
    basis.push(n + n_slack + i)
  }
  let t : SimplexTableau = { a, b, basis, n_total, }
  // Phase 1: min sum(artificial) = max -sum(artificial)
  let c1 : Array[Double] = Array::make(n_total, 0.0)
  for i = 0; i < rows; i = i + 1 {
    c1[n + n_slack + i] = -1.0
  }
  let art_cols : Array[Int32] = Array(capacity=0)
  for i = 0; i < rows; i = i + 1 {
    art_cols.push(n + n_slack + i)
  }
  let phase1 = simplex_solve(t, c1, art_cols)
  match phase1 {
    None => return None
    Some((obj1, _)) =>
      if obj1 < -1.0e-6 {
        // 人工变量非零 -> 原问题不可行
        return None
      }
  }
  // Phase 2: 原目标(扩展到 n_total, slack/人工为 0)
  let c2 : Array[Double] = Array::make(n_total, 0.0)
  for j = 0; j < n; j = j + 1 {
    c2[j] = c[j]
  }
  match simplex_solve(t, c2, art_cols) {
    None => None
    Some((obj, x_full)) => {
      // 提取原始 n 个分量
      let x : Array[Double] = Array(capacity=n)
      for j = 0; j < n; j = j + 1 {
        x.push(x_full[j])
      }
      Some((obj, x))
    }
  }
}

// ---------------------------------------------------------------------------
// FBA 模型求解
// ---------------------------------------------------------------------------

///|
/// 求解 FBA: 最大化 objective 反应通量和.
pub fn fba_solve(model : FbaModel) -> FbaSolution raise FbaError {
  let n_r = model.reactions.length()
  if n_r == 0 {
    raise FbaError("model has no reactions")
  }
  let m = model.metabolites.length()
  // 变量拆分: v_j = base_j + x+ - x-
  // base_j = max(lb, 0); x+ 上界 ub-base(>0 时), x- 上界 base-lb(lb<0 时)
  let n_var = 2 * n_r
  let c : Array[Double] = Array::make(n_var, 0.0)
  let ub : Array[Double] = Array::make(n_var, 0.0)
  let base : Array[Double] = Array::make(n_r, 0.0)
  for j = 0; j < n_r; j = j + 1 {
    let r = model.reactions[j]
    let b_j = if r.lower_bound > 0.0 { r.lower_bound } else { 0.0 }
    base[j] = b_j
    // x+ 列 j*2, x- 列 j*2+1
    let fwd_ub = r.upper_bound - b_j
    let rev_ub = b_j - r.lower_bound
    ub[j * 2] = if fwd_ub > 0.0 { fwd_ub } else { 0.0 }
    ub[j * 2 + 1] = if rev_ub > 0.0 { rev_ub } else { 0.0 }
    c[j * 2] = r.objective
    c[j * 2 + 1] = -r.objective
  }
  // 稳态约束: S·(base + x+ - x-) = 0 -> S·(x+ - x-) = -S·base
  let a_eq : Array[Array[Double]] = Array(capacity=0)
  let b_eq : Array[Double] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    let row : Array[Double] = Array::make(n_var, 0.0)
    let mut rhs = 0.0
    for j = 0; j < n_r; j = j + 1 {
      let sij = model.stoichiometry[i][j]
      if sij != 0.0 {
        row[j * 2] = sij
        row[j * 2 + 1] = -sij
        rhs = rhs - sij * base[j]
      }
    }
    a_eq.push(row)
    b_eq.push(rhs)
  }
  // 求解
  let result = fba_solve_lp(a_eq, b_eq, c, ub)
  match result {
    None => {
      let fluxes : Array[Double] = Array::make(n_r, 0.0)
      { objective_value: 0.0, fluxes, status: "infeasible_or_unbounded", }
    }
    Some((obj, x)) => {
      let fluxes : Array[Double] = Array(capacity=n_r)
      for j = 0; j < n_r; j = j + 1 {
        fluxes.push(base[j] + x[j * 2] - x[j * 2 + 1])
      }
      { objective_value: obj, fluxes, status: "optimal", }
    }
  }
}

///|
/// 基因敲除: 置零关联反应的 bounds 后重解.
/// 返回 (敲除后模型, 敲除的反应列表).
pub fn fba_knockout(
  model : FbaModel,
  genes : Array[String],
) -> (FbaModel, Array[String]) {
  let reactions : Array[FbaReaction] = Array(capacity=0)
  let knocked : Array[String] = Array(capacity=0)
  for r in model.reactions {
    let hit = {
      let mut h = false
      for g in genes {
        if r.genes.contains(g) {
          h = true
          break
        }
      }
      h
    }
    if hit && r.objective == 0.0 {
      // 敲除非目标反应
      knocked.push(r.id)
      reactions.push({
        id: r.id,
        name: r.name,
        lower_bound: 0.0,
        upper_bound: 0.0,
        objective: r.objective,
        genes: r.genes,
      })
    } else {
      reactions.push(r)
    }
  }
  (
    {
      metabolites: model.metabolites,
      reactions,
      stoichiometry: model.stoichiometry,
    },
    knocked,
  )
}

///|
/// 基因敲除表型模拟: 逐基因敲除返回生长与否.
pub fn fba_knockout_screen(
  model : FbaModel,
  essential_threshold~ : Double = 1.0e-6,
) -> Array[(String, Double, String)] raise FbaError {
  // 野生型基线
  let wt = fba_solve(model)
  let gene_set : Map[String, Bool] = Map([], capacity=0)
  for r in model.reactions {
    for g in r.genes {
      gene_set.set(g, true)
    }
  }
  let out : Array[(String, Double, String)] = Array(capacity=0)
  for g in gene_set.keys() {
    let (ko_model, _) = fba_knockout(model, [g])
    let sol = fba_solve(ko_model)
    let ratio = if wt.objective_value > 1.0e-9 {
      sol.objective_value / wt.objective_value
    } else {
      0.0
    }
    let phenotype = if sol.objective_value < essential_threshold {
      "non-growth"
    } else if ratio < 0.5 {
      "reduced"
    } else {
      "growth"
    }
    out.push((g, sol.objective_value, phenotype))
  }
  out.sort_by(fn(a, b) { compare(a.0, b.0) })
  out
}

// ---------------------------------------------------------------------------
// 模型构造便捷函数与 JSON IO
// ---------------------------------------------------------------------------

///|
/// 空模型.
pub fn fba_model_new(metabolites : Array[String]) -> FbaModel {
  {
    metabolites,
    reactions: Array(capacity=0),
    stoichiometry: Array::make(metabolites.length(), []),
  }
}

///|
/// 添加反应(同时扩展 S 矩阵列).
pub fn FbaModel::add_reaction(
  self : FbaModel,
  reaction : FbaReaction,
  /// 参与计量: (metabolite, 系数)
  participants : Array[(String, Double)],
) -> FbaModel raise FbaError {
  // 校验 metabolite 存在
  for p in participants {
    if !self.metabolites.contains(p.0) {
      raise FbaError("unknown metabolite: \{p.0}")
    }
  }
  let reactions = Array(capacity=0)
  for r in self.reactions {
    reactions.push(r)
  }
  reactions.push(reaction)
  let n_new = reactions.length()
  let stoich : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < self.metabolites.length(); i = i + 1 {
    let row : Array[Double] = Array::make(n_new, 0.0)
    for j = 0; j < n_new - 1; j = j + 1 {
      row[j] = self.stoichiometry[i][j]
    }
    stoich.push(row)
  }
  for p in participants {
    let mi = {
      let mut idx = -1
      for i = 0; i < self.metabolites.length(); i = i + 1 {
        if self.metabolites[i] == p.0 {
          idx = i
          break
        }
      }
      idx
    }
    if mi >= 0 {
      stoich[mi][n_new - 1] = p.1
    }
  }
  { metabolites: self.metabolites, reactions, stoichiometry: stoich, }
}

///|
/// 修改指定反应的上下界(返回新模型).
pub fn fba_set_reaction_bound(
  model : FbaModel,
  reaction_id : String,
  lower_bound : Double,
  upper_bound : Double,
) -> FbaModel raise FbaError {
  let reactions : Array[FbaReaction] = Array(capacity=0)
  let mut found = false
  for r in model.reactions {
    if r.id == reaction_id {
      found = true
      reactions.push({
        id: r.id,
        name: r.name,
        lower_bound,
        upper_bound,
        objective: r.objective,
        genes: r.genes,
      })
    } else {
      reactions.push(r)
    }
  }
  if !found {
    raise FbaError("reaction not found: \{reaction_id}")
  }
  {
    metabolites: model.metabolites,
    reactions,
    stoichiometry: model.stoichiometry,
  }
}

///|
/// 模型 -> JSON 文本(交换格式).
pub fn fba_write_json(model : FbaModel) -> String {
  let rxns : Array[JsonValue] = Array(capacity=0)
  for r in model.reactions {
    let genes : Array[JsonValue] = Array(capacity=0)
    for g in r.genes {
      genes.push(JsonStr(g))
    }
    rxns.push(
      JsonObj(
        Map([
          ("id", JsonStr(r.id)),
          ("name", JsonStr(r.name)),
          ("lower_bound", JsonNum(r.lower_bound)),
          ("upper_bound", JsonNum(r.upper_bound)),
          ("objective", JsonNum(r.objective)),
          ("genes", JsonArr(genes)),
        ]),
      ),
    )
  }
  let mets : Array[JsonValue] = Array(capacity=0)
  for m in model.metabolites {
    mets.push(JsonStr(m))
  }
  let s_rows : Array[JsonValue] = Array(capacity=0)
  for row in model.stoichiometry {
    let vals : Array[JsonValue] = Array(capacity=0)
    for v in row {
      vals.push(JsonNum(v))
    }
    s_rows.push(JsonArr(vals))
  }
  json_to_string(
    JsonObj(
      Map([
        ("metabolites", JsonArr(mets)),
        ("reactions", JsonArr(rxns)),
        ("stoichiometry", JsonArr(s_rows)),
      ]),
    ),
  )
}

///|
/// 解析 JSON 模型文本.
pub fn fba_parse_json(text : String) -> FbaModel? {
  let v = parse_json(text)
  match v {
    Some(JsonObj(fields)) => {
      let mut metabolites : Array[String] = []
      let mut reactions : Array[FbaReaction] = []
      let mut stoich : Array[Array[Double]] = []
      for f in fields {
        if f.0 == "metabolites" {
          match f.1 {
            JsonArr(arr) => {
              metabolites = []
              for x in arr {
                match x {
                  JsonStr(s) => metabolites.push(s)
                  _ => ()
                }
              }
            }
            _ => ()
          }
        } else if f.0 == "reactions" {
          match f.1 {
            JsonArr(arr) => {
              reactions = []
              for x in arr {
                match x {
                  JsonObj(rf) => {
                    let mut id = ""
                    let mut name = ""
                    let mut lb = 0.0
                    let mut ub = 1000.0
                    let mut obj = 0.0
                    let genes : Array[String] = []
                    for kv in rf {
                      match kv.0 {
                        "id" =>
                          match kv.1 {
                            JsonStr(s) => id = s
                            _ => ()
                          }
                        "name" =>
                          match kv.1 {
                            JsonStr(s) => name = s
                            _ => ()
                          }
                        "lower_bound" =>
                          match kv.1 {
                            JsonNum(x) => lb = x
                            _ => ()
                          }
                        "upper_bound" =>
                          match kv.1 {
                            JsonNum(x) => ub = x
                            _ => ()
                          }
                        "objective" =>
                          match kv.1 {
                            JsonNum(x) => obj = x
                            _ => ()
                          }
                        "genes" =>
                          match kv.1 {
                            JsonArr(gs) =>
                              for g in gs {
                                match g {
                                  JsonStr(s) => genes.push(s)
                                  _ => ()
                                }
                              }
                            _ => ()
                          }
                        _ => ()
                      }
                    }
                    reactions.push({
                      id,
                      name,
                      lower_bound: lb,
                      upper_bound: ub,
                      objective: obj,
                      genes,
                    })
                  }
                  _ => ()
                }
              }
            }
            _ => ()
          }
        } else if f.0 == "stoichiometry" {
          match f.1 {
            JsonArr(rows) => {
              stoich = []
              for rv in rows {
                match rv {
                  JsonArr(vals) => {
                    let row : Array[Double] = []
                    for x in vals {
                      match x {
                        JsonNum(d) => row.push(d)
                        _ => ()
                      }
                    }
                    stoich.push(row)
                  }
                  _ => ()
                }
              }
            }
            _ => ()
          }
        }
      }
      if metabolites.length() == 0 && reactions.length() == 0 {
        return None
      }
      Some({ metabolites, reactions, stoichiometry: stoich, })
    }
    _ => None
  }
}

///|
/// 求解摘要.
pub fn FbaSolution::summary(self : FbaSolution, model : FbaModel) -> String {
  let mut n_active = 0
  for v in self.fluxes {
    if v > 1.0e-6 || v < -1.0e-6 {
      n_active = n_active + 1
    }
  }
  "fba-lite: status=\{self.status}, objective=\{self.objective_value}, active fluxes \{n_active}/\{model.reactions.length()}"
}
