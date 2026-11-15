///|
/// Double Principal Coordinates Analysis (DPCoA, Pavoine et al. 2004).
///
/// DPCoA 把系统发育树和 OTU abundance table 结合起来：
///
/// 1. 从带分支长度的 tree 计算所有 tip 两两之间的 phylogenetic distance
///    （叶节点 i,j 之间的路径长度 = depth(i) + depth(j) − 2·depth(LCA(i,j))）
/// 2. 对 OTU × OTU 距离矩阵做 PCoA → 每个 OTU 在欧式 PCoA 空间的坐标
/// 3. 用 sample 的 OTU abundance 作为权重，对 OTU 坐标加权平均
///    → 每个 sample 的主坐标（这就是"Double"的含义：先对 OTU 层面 PCoA，
///    再对 sample 层面 abundance 加权）
///
/// 这种做法等价于 weighted UniFrac distance matrix 直接 PCoA，但能输出
/// **OTU loadings**（哪些 OTU 驱动了样本间差异），这是 DPCoA 相对于
/// UniFrac+PCoA 的独特优势。

///|
/// DPCoA 返回值。
pub struct SkbioDpcoaResult {
  /// 每个 sample 在 PCoA 空间的坐标（sample_coords[i][axis]）。
  sample_coords : Array[Array[Double]]
  /// 每个 OTU 在 PCoA 空间的坐标（otu_coords[otu_idx][axis]）。
  otu_coords : Array[Array[Double]]
  /// 对应每个 axis 的特征值（从大到小，负特征值截断为 0）。
  eigvals : Array[Double]
  /// 每个 axis 的方差解释比例。
  proportion_explained : Array[Double]
  /// OTU 名称（与 otu_coords 行对齐）。
  taxa : Array[String]
  /// OTU 树节点中实际出现在 taxa 里的 depth（到 root 的累计分支长度）。
  otu_depths : Array[Double]
}

// ---------------------------------------------------------------------------
// 内部：从 tree 算 OTU × OTU 距离矩阵
// ---------------------------------------------------------------------------

///|
/// 对每个节点返回累计 depth（从 root 到该节点的累计分支长度），以及所有
/// tip 的 depth + root path。这是一个 O(N) 的 preorder 遍历。
fn skbd_dpcoa_compute_depths(
  node : SkbioNode,
  parent_depth : Double,
  tip_names : Array[String],
  tip_depths : Array[Double],
) -> Unit {
  let my_depth = parent_depth + node.length
  if node.is_tip() {
    match node.name {
      Some(n) => {
        tip_names.push(n)
        tip_depths.push(my_depth)
      }
      None => ()
    }
  } else {
    for child in node.children {
      skbd_dpcoa_compute_depths(child, my_depth, tip_names, tip_depths)
    }
  }
}

///|
/// 在 tree 中收集每个 tip 到 LCA 的 depth。做法：preorder 遍历，对每个 tip
/// 记录从 root 到它的路径上所有祖先的 name → depth 映射。然后对任意两个
/// tip 做路径交集取最深公共节点。
///
/// 为了简单高效，这里用 "每个 tip 到所有祖先的 depth 映射" 方案：
///
/// - 维护一个 `ancestor_depths : Array[(String, Double)]` 栈（从 root 到当前路径）
/// - preorder 时 push(node.name, cumulative_depth) → recurse → pop
/// - 遇到 tip 时保存它的 ancestor_depths 快照
///
/// 最后对每对 tip 从快照交集里找最深的公共祖先。
struct AncestorSnapshot {
  names : Array[String]
  depths : Array[Double]
}

fn skbd_dpcoa_collect_tip_ancestors(
  node : SkbioNode,
  parent_depth : Double,
  stack_names : Array[String],
  stack_depths : Array[Double],
  out_names : Array[String],
  out_snapshots : Array[AncestorSnapshot],
) -> Unit {
  let my_depth = parent_depth + node.length
  // push 当前节点
  match node.name {
    Some(n) => {
      stack_names.push(n)
      stack_depths.push(my_depth)
    }
    None => ()
  }
  if node.is_tip() {
    match node.name {
      Some(n) => {
        out_names.push(n)
        // 复制当前 stack 作为这个 tip 的祖先快照
        let snap_names : Array[String] = Array(capacity=stack_names.length())
        let snap_depths : Array[Double] = Array(capacity=stack_depths.length())
        for s in stack_names { snap_names.push(s) }
        for d in stack_depths { snap_depths.push(d) }
        out_snapshots.push(AncestorSnapshot::{ names: snap_names, depths: snap_depths })
      }
      None => ()
    }
  } else {
    for child in node.children {
      skbd_dpcoa_collect_tip_ancestors(
        child, my_depth, stack_names, stack_depths, out_names, out_snapshots,
      )
    }
  }
  // pop 当前节点（如果 push 了）
  match node.name {
    Some(_) => {
      let _ = stack_names.remove(stack_names.length() - 1)
      let _ = stack_depths.remove(stack_depths.length() - 1)
    }
    None => ()
  }
}

///|
/// 计算两个 snapshot 的 LCA depth（最深的公共祖先 depth）。
fn skbd_dpcoa_lca_depth(a : AncestorSnapshot, b : AncestorSnapshot) -> Double {
  // hash map 查找（这里用线性 scan：树不会太深）
  let mut best = 0.0
  for i = 0; i < a.names.length(); i = i + 1 {
    let na = a.names[i]
    let da = a.depths[i]
    // 在 b 里找
    for j = 0; j < b.names.length(); j = j + 1 {
      if b.names[j] == na && b.depths[j] > best {
        best = b.depths[j]
      }
    }
  }
  best
}

// ---------------------------------------------------------------------------
// DPCoA 主入口
// ---------------------------------------------------------------------------

///|
/// Double Principal Coordinates Analysis.
///
/// - `tree`: 带分支长度的系统发育树（`SkbioNode`）
/// - `samples`: OTU table，每行一个 sample，列顺序与 `taxa` 一致
/// - `taxa`: OTU 名称数组（必须是 tree.tips() 的子集；顺序决定 samples 的列解释）
///
/// 返回 sample 坐标 + OTU 坐标 + 特征值 + 解释比例 + 各 OTU 到 root 的 depth。
pub fn skbio_dpcoa(
  tree : SkbioNode,
  samples : Array[Array[Double]],
  taxa : Array[String],
) -> SkbioDpcoaResult {
  // 1. 从 tree 收集所有 tip + 每个 tip 的祖先快照
  let all_tip_names : Array[String] = Array(capacity=0)
  let all_tip_snaps : Array[AncestorSnapshot] = Array(capacity=0)
  let empty_names : Array[String] = Array(capacity=0)
  let empty_depths : Array[Double] = Array(capacity=0)
  skbd_dpcoa_collect_tip_ancestors(
    tree, 0.0, empty_names, empty_depths,
    all_tip_names, all_tip_snaps,
  )

  // 2. 构建 taxa → tip snapshot 的映射（只保留 taxa 里出现的）
  let n_taxa = taxa.length()
  let taxa_snaps : Array[AncestorSnapshot] = Array(capacity=n_taxa)
  let taxa_depths : Array[Double] = Array(capacity=n_taxa)
  for t in taxa {
    let mut found = false
    for i = 0; i < all_tip_names.length(); i = i + 1 {
      if all_tip_names[i] == t {
        taxa_snaps.push(all_tip_snaps[i])
        // taxa 的 depth = 祖先路径上最后一个（即它自己的 depth）
        let snap = all_tip_snaps[i]
        taxa_depths.push(snap.depths[snap.depths.length() - 1])
        found = true
        break
      }
    }
    if !found {
      // 不在 tree 里 → depth=0, snap 空
      taxa_snaps.push(AncestorSnapshot::{ names: [], depths: [] })
      taxa_depths.push(0.0)
    }
  }

  // 3. 构建 OTU × OTU phylogenetic distance matrix
  let dm : Array[Array[Double]] = Array(capacity=n_taxa)
  for i = 0; i < n_taxa; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_taxa)
    for j = 0; j < n_taxa; j = j + 1 {
      if i == j {
        row.push(0.0)
      } else {
        let lca_d = skbd_dpcoa_lca_depth(taxa_snaps[i], taxa_snaps[j])
        let d = taxa_depths[i] + taxa_depths[j] - 2.0 * lca_d
        row.push(d)
      }
    }
    dm.push(row)
  }

  // 4. PCoA → OTU 在欧式空间的坐标
  let pcoa = skbio_pcoa(dm)
  let otu_coords = pcoa.samples
  let k = if otu_coords.length() > 0 { otu_coords[0].length() } else { 0 }

  // 5. abundance 加权 → sample 坐标
  let n_samples = samples.length()
  let sample_coords : Array[Array[Double]] = Array(capacity=n_samples)
  for i = 0; i < n_samples; i = i + 1 {
    // 归一化 abundance（转相对丰度）
    let row = samples[i]
    let mut total = 0.0
    for v in row { total = total + v }
    let coords_row = Array::make(k, 0.0)
    if total > 0.0 {
      for otu = 0; otu < n_taxa; otu = otu + 1 {
        let w = row[otu] / total
        for ax = 0; ax < k; ax = ax + 1 {
          coords_row[ax] = coords_row[ax] + w * otu_coords[otu][ax]
        }
      }
    }
    sample_coords.push(coords_row)
  }

  SkbioDpcoaResult::{
    sample_coords,
    otu_coords,
    eigvals: pcoa.eigvals,
    proportion_explained: pcoa.proportion_explained,
    taxa,
    otu_depths: taxa_depths,
  }
}
