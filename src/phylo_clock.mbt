///|
/// phylo_clock.mbt — 分子钟标定与祖先序列重建.
///
/// 对齐 TreeTime 的核心语义(本地确定性实现, 数据不出本机):
/// - root-to-tip 回归: 用叶节点采样日期对 root-to-tip 距离做最小二乘拟合,
///   估计替换速率、根节点日期与 R²
/// - 时间树: 按拟合速率把分支长度换算为时间单位, 为每个节点标注时间
/// - 祖先状态重建: 在 Jukes-Cantor 模型下用 Felsenstein 剪枝算法计算
///   每个节点逐位点的边际后验概率, 输出共识碱基与置信度
/// - JSON 输出(衔接 Web/报告架构)

///|
/// 分子钟处理错误.
pub suberror ClockError {
  ClockError(String)
}

///|
/// 带采样日期的叶节点及其 root-to-tip 距离.
pub(all) struct ClockTip {
  /// 叶节点名称
  name : String
  /// 采样日期(年, 可为小数)
  date : Double
  /// 从根到该叶节点的分支长度之和
  root_to_tip : Double
  /// 回归残差
  residual : Double
}

///|
/// root-to-tip 回归结果.
pub(all) struct RootToTipResult {
  /// 估计的替换速率(分支长度 / 年)
  rate : Double
  /// 回归截距
  intercept : Double
  /// 决定系数 R²
  r_squared : Double
  /// 由回归推算的根节点日期
  root_date : Double
  /// 每个叶节点的拟合明细
  tips : Array[ClockTip]
}

///|
/// 时间树节点: `time` 为该节点在时间轴上的日期.
pub(all) struct DatedNode {
  /// 节点名称(内部无名节点为 None)
  name : String?
  /// 节点日期
  time : Double
  /// 子节点
  children : Array[DatedNode]
}

///|
/// 单个位点的祖先状态重建结果.
pub(all) struct AncestralSite {
  /// 后验概率最高的碱基
  base : Char
  /// 该碱基的后验概率
  probability : Double
  /// A/C/G/T 四态后验概率
  posterior : Array[Double]
}

///|
/// 单个节点的祖先序列.
pub(all) struct AncestralTaxon {
  /// 节点名称
  name : String
  /// 逐位点共识序列
  mut sequence : String
  /// 逐位点后验概率明细
  sites : Array[AncestralSite]
}

///|
/// 祖先序列重建结果.
pub(all) struct AncestralResult {
  /// 整棵树的位点对数似然
  log_likelihood : Double
  /// 每个节点的重建序列(后序遍历)
  taxa : Array[AncestralTaxon]
}

///|
/// 四态顺序: A/C/G/T.
let clock_states : Array[Char] = ['A', 'C', 'G', 'T']

///|
/// 碱基字符 -> 状态下标; 未知碱基(N/-/?)返回 None.
fn clock_state_index(code : Int32) -> Int32? {
  match code {
    65 | 97 => Some(0)
    67 | 99 => Some(1)
    71 | 103 => Some(2)
    84 | 116 => Some(3)
    _ => None
  }
}

///|
/// 节点名称(内部无名节点返回空串).
fn clock_node_name(name : String?) -> String {
  match name {
    Some(n) => n
    None => ""
  }
}

///|
/// 自底向上收集所有叶节点及其 root-to-tip 距离.
fn clock_collect_tips(
  node : SkbioNode,
  distance : Double,
  dates : Map[String, Double],
  out : Array[ClockTip],
) -> Unit raise ClockError {
  let dist = distance + node.length
  if node.children.length() == 0 {
    let name = clock_node_name(node.name)
    if name == "" {
      raise ClockError("tip node has no name")
    }
    match dates.get(name) {
      Some(date) =>
        out.push(ClockTip::{ name, date, root_to_tip: dist, residual: 0.0, })
      None => raise ClockError("missing sampling date for tip: \{name}")
    }
  } else {
    for child in node.children {
      clock_collect_tips(child, dist, dates, out)
    }
  }
}

///|
/// 在采样日期与 root-to-tip 距离之间做最小二乘线性回归.
///
/// 返回 (slope, intercept, r_squared).
fn clock_linear_regression(
  xs : Array[Double],
  ys : Array[Double],
) -> (Double, Double, Double) raise ClockError {
  let n = xs.length()
  if n < 2 {
    raise ClockError("at least two dated tips are required")
  }
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  for i = 0; i < n; i = i + 1 {
    sum_x = sum_x + xs[i]
    sum_y = sum_y + ys[i]
  }
  let mean_x = sum_x / n.to_double()
  let mean_y = sum_y / n.to_double()
  let mut sxx = 0.0
  let mut sxy = 0.0
  let mut syy = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = xs[i] - mean_x
    let dy = ys[i] - mean_y
    sxx = sxx + dx * dx
    sxy = sxy + dx * dy
    syy = syy + dy * dy
  }
  if sxx <= 0.0 {
    raise ClockError("sampling dates must not all be identical")
  }
  let slope = sxy / sxx
  let intercept = mean_y - slope * mean_x
  let mut ss_res = 0.0
  for i = 0; i < n; i = i + 1 {
    let pred = intercept + slope * xs[i]
    let r = ys[i] - pred
    ss_res = ss_res + r * r
  }
  let r_squared = if syy <= 0.0 { 1.0 } else { 1.0 - ss_res / syy }
  (slope, intercept, r_squared)
}

///|
/// 用采样日期标定分子钟: 对 root-to-tip 距离做线性回归.
///
/// 叶节点必须全部在 `dates` 中有记录, 且日期不能完全相同.
pub fn phylo_clock_root_to_tip(
  root : SkbioNode,
  dates : Map[String, Double],
) -> RootToTipResult raise ClockError {
  let raw : Array[ClockTip] = Array(capacity=0)
  clock_collect_tips(root, 0.0, dates, raw)
  let n = raw.length()
  let xs : Array[Double] = Array(capacity=n)
  let ys : Array[Double] = Array(capacity=n)
  for tip in raw {
    xs.push(tip.date)
    ys.push(tip.root_to_tip)
  }
  let (rate, intercept, r_squared) = clock_linear_regression(xs, ys)
  if rate <= 0.0 {
    raise ClockError(
      "root-to-tip distance does not increase with sampling date",
    )
  }
  let root_date = -intercept / rate
  let tips : Array[ClockTip] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let pred = intercept + rate * xs[i]
    tips.push(ClockTip::{
      name: raw[i].name,
      date: raw[i].date,
      root_to_tip: raw[i].root_to_tip,
      residual: raw[i].root_to_tip - pred,
    })
  }
  RootToTipResult::{ rate, intercept, r_squared, root_date, tips, }
}

///|
/// 递归构造时间树.
fn clock_build_dated_tree(
  node : SkbioNode,
  distance : Double,
  fit : RootToTipResult,
) -> DatedNode {
  let dist = distance + node.length
  let children : Array[DatedNode] = Array(capacity=node.children.length())
  for child in node.children {
    children.push(clock_build_dated_tree(child, dist, fit))
  }
  let name = match node.name {
    Option::Some(n) => if n == "" { Option::None } else { Option::Some(n) }
    Option::None => Option::None
  }
  DatedNode::{ name, time: fit.root_date + dist / fit.rate, children, }
}

///|
/// 由 root-to-tip 回归结果生成时间树.
pub fn phylo_clock_dated_tree(
  root : SkbioNode,
  fit : RootToTipResult,
) -> DatedNode {
  clock_build_dated_tree(root, 0.0, fit)
}

///|
/// 扁平化节点, 供逐位点似然计算使用(后序: 子节点在父节点之前).
priv struct ClockFlatNode {
  name : String
  length : Double
  children : Array[Int32]
}

///|
/// 递归扁平化 Newick 树.
fn clock_flatten_tree(node : SkbioNode, nodes : Array[ClockFlatNode]) -> Int32 {
  let child_ids : Array[Int32] = Array(capacity=node.children.length())
  for child in node.children {
    child_ids.push(clock_flatten_tree(child, nodes))
  }
  let id = nodes.length()
  let raw_name = clock_node_name(node.name)
  let name = if raw_name == "" { "node_" + id.to_string() } else { raw_name }
  nodes.push(ClockFlatNode::{ name, length: node.length, children: child_ids, })
  id
}

///|
/// Jukes-Cantor 转换概率: 分支长度为 t 时同态概率.
fn clock_jc_same(t : Double) -> Double {
  let e = @math.exp(-4.0 / 3.0 * t)
  0.25 + 0.75 * e
}

///|
/// Jukes-Cantor 转换概率: 分支长度为 t 时异态概率.
fn clock_jc_diff(t : Double) -> Double {
  let e = @math.exp(-4.0 / 3.0 * t)
  0.25 - 0.25 * e
}

///|
/// 叶节点在指定位点的观测四态似然; 未知碱基四态均设为 1.
fn clock_observed_likelihood(
  sequence : String,
  pos : Int32,
  states : Array[Double],
) -> Unit {
  match clock_state_index(sequence.unsafe_get(pos).to_int()) {
    Some(s) =>
      for i = 0; i < 4; i = i + 1 {
        states[i] = if i == s { 1.0 } else { 0.0 }
      }
    None =>
      for i = 0; i < 4; i = i + 1 {
        states[i] = 1.0
      }
  }
}

///|
/// 计算一棵树的祖先序列.
///
/// 使用 Jukes-Cantor 模型 + Felsenstein 剪枝:
/// 1. 后序遍历计算每个节点在当前位点的子树似然
/// 2. 前序遍历计算每个节点的外部似然
/// 3. 两者相乘得到边际后验概率
pub fn phylo_clock_ancestral(
  root : SkbioNode,
  sequences : Map[String, String],
) -> AncestralResult raise ClockError {
  let nodes : Array[ClockFlatNode] = Array(capacity=0)
  let root_id = clock_flatten_tree(root, nodes)
  let n_nodes = nodes.length()

  let mut seq_len = -1
  for id = 0; id < n_nodes; id = id + 1 {
    let nd = nodes[id]
    if nd.children.length() == 0 {
      match sequences.get(nd.name) {
        Some(seq) =>
          if seq_len < 0 {
            seq_len = seq.length()
          } else if seq.length() != seq_len {
            raise ClockError("all tip sequences must have the same length")
          }
        None => raise ClockError("missing sequence for tip: \{nd.name}")
      }
    }
  }
  if seq_len < 0 {
    raise ClockError("tree has no tips")
  }

  let p_same : Array[Double] = Array(capacity=n_nodes)
  let p_diff : Array[Double] = Array(capacity=n_nodes)
  for id = 0; id < n_nodes; id = id + 1 {
    p_same.push(clock_jc_same(nodes[id].length))
    p_diff.push(clock_jc_diff(nodes[id].length))
  }

  let results : Array[AncestralTaxon] = Array(capacity=n_nodes)
  for id = 0; id < n_nodes; id = id + 1 {
    results.push(AncestralTaxon::{
      name: nodes[id].name,
      sequence: "",
      sites: [],
    })
  }

  let sub : Array[Double] = Array::make(n_nodes * 4, 0.0)
  let outside : Array[Double] = Array::make(n_nodes * 4, 0.0)
  let site_chars : Array[StringBuilder] = Array(capacity=n_nodes)
  for id = 0; id < n_nodes; id = id + 1 {
    site_chars.push(StringBuilder(size_hint=seq_len))
  }

  let mut log_likelihood = 0.0
  for pos = 0; pos < seq_len; pos = pos + 1 {
    // 后序: 子树似然.
    for id = 0; id < n_nodes; id = id + 1 {
      let nd = nodes[id]
      if nd.children.length() == 0 {
        let seq = sequences.get(nd.name).unwrap()
        let states : Array[Double] = Array::make(4, 0.0)
        clock_observed_likelihood(seq, pos, states)
        for s = 0; s < 4; s = s + 1 {
          sub[id * 4 + s] = states[s]
        }
      } else {
        for s = 0; s < 4; s = s + 1 {
          let mut v = 1.0
          for child in nd.children {
            let mut sum = 0.0
            for u = 0; u < 4; u = u + 1 {
              let p = if s == u { p_same[child] } else { p_diff[child] }
              sum = sum + p * sub[child * 4 + u]
            }
            v = v * sum
          }
          sub[id * 4 + s] = v
        }
      }
    }

    // 根节点先验均匀, 外部似然为 1.
    for s = 0; s < 4; s = s + 1 {
      outside[root_id * 4 + s] = 1.0
    }
    // 前序(节点编号从大到小): 由父节点推到子节点.
    for step = 0; step < n_nodes; step = step + 1 {
      let id = n_nodes - 1 - step
      let nd = nodes[id]
      if nd.children.length() == 0 {
        continue
      }
      for c_idx = 0; c_idx < nd.children.length(); c_idx = c_idx + 1 {
        let child = nd.children[c_idx]
        for t = 0; t < 4; t = t + 1 {
          let mut sum = 0.0
          for s = 0; s < 4; s = s + 1 {
            let mut prod = 1.0
            for o_idx = 0; o_idx < nd.children.length(); o_idx = o_idx + 1 {
              if o_idx != c_idx {
                let other = nd.children[o_idx]
                let mut l = 0.0
                for u = 0; u < 4; u = u + 1 {
                  let p = if s == u { p_same[other] } else { p_diff[other] }
                  l = l + p * sub[other * 4 + u]
                }
                prod = prod * l
              }
            }
            let p = if s == t { p_same[child] } else { p_diff[child] }
            sum = sum + outside[id * 4 + s] * p * prod
          }
          outside[child * 4 + t] = sum
        }
      }
    }

    // 根节点似然.
    let mut root_lik = 0.0
    for s = 0; s < 4; s = s + 1 {
      root_lik = root_lik + 0.25 * sub[root_id * 4 + s]
    }
    if root_lik > 0.0 {
      log_likelihood = log_likelihood + @math.ln(root_lik)
    }

    // 逐节点归一化边际后验并记录共识碱基.
    for id = 0; id < n_nodes; id = id + 1 {
      let posterior : Array[Double] = Array::make(4, 0.0)
      let mut total = 0.0
      for s = 0; s < 4; s = s + 1 {
        let v = outside[id * 4 + s] * sub[id * 4 + s]
        posterior[s] = v
        total = total + v
      }
      if total > 0.0 {
        for s = 0; s < 4; s = s + 1 {
          posterior[s] = posterior[s] / total
        }
      } else {
        for s = 0; s < 4; s = s + 1 {
          posterior[s] = 0.25
        }
      }
      let mut best = 0
      for s = 1; s < 4; s = s + 1 {
        if posterior[s] > posterior[best] {
          best = s
        }
      }
      site_chars[id].write_char(clock_states[best])
      results[id].sites.push(AncestralSite::{
        base: clock_states[best],
        probability: posterior[best],
        posterior,
      })
    }
  }

  for id = 0; id < n_nodes; id = id + 1 {
    results[id].sequence = site_chars[id].to_string()
  }
  AncestralResult::{ log_likelihood, taxa: results, }
}

///|
/// root-to-tip 拟合结果 -> JSON 文本.
pub fn phylo_clock_root_to_tip_json(fit : RootToTipResult) -> String {
  let tips : Array[JsonValue] = Array(capacity=fit.tips.length())
  for t in fit.tips {
    tips.push(
      JsonObj(
        Map([
          ("name", JsonStr(t.name)),
          ("date", JsonNum(t.date)),
          ("root_to_tip", JsonNum(t.root_to_tip)),
          ("residual", JsonNum(t.residual)),
        ]),
      ),
    )
  }
  json_to_string(
    JsonObj(
      Map([
        ("rate", JsonNum(fit.rate)),
        ("intercept", JsonNum(fit.intercept)),
        ("r_squared", JsonNum(fit.r_squared)),
        ("root_date", JsonNum(fit.root_date)),
        ("tips", JsonArr(tips)),
      ]),
    ),
  )
}

///|
/// 祖先重建结果 -> JSON 文本.
pub fn phylo_clock_ancestral_json(result : AncestralResult) -> String {
  let taxa : Array[JsonValue] = Array(capacity=result.taxa.length())
  for taxon in result.taxa {
    let sites : Array[JsonValue] = Array(capacity=taxon.sites.length())
    for site in taxon.sites {
      let posterior : Array[JsonValue] = Array(capacity=4)
      for v in site.posterior {
        posterior.push(JsonNum(v))
      }
      sites.push(
        JsonObj(
          Map([
            ("base", JsonStr(site.base.to_string())),
            ("probability", JsonNum(site.probability)),
            ("posterior", JsonArr(posterior)),
          ]),
        ),
      )
    }
    taxa.push(
      JsonObj(
        Map([
          ("name", JsonStr(taxon.name)),
          ("sequence", JsonStr(taxon.sequence)),
          ("sites", JsonArr(sites)),
        ]),
      ),
    )
  }
  json_to_string(
    JsonObj(
      Map([
        ("log_likelihood", JsonNum(result.log_likelihood)),
        ("taxa", JsonArr(taxa)),
      ]),
    ),
  )
}
