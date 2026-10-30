///|
/// kraken_lite.mbt — 宏基因组物种丰度定量(Kraken/Bracken 风格).
///
/// 本地确定性实现, 数据不出本机:
/// - 建库: 参考基因组 -> 去重 k-mer -> 所有命中基因组的最小共同祖先(LCA)
/// - 分类: 读段 k-mer 命中后, 沿分类树选择累计支持度最高的路径
/// - 丰度: Kraken 风格累计计数报告 + Bracken 风格属/种级重分摊
/// - JSON / TSV 输出(衔接报告架构)

///|
/// 宏基因组定量错误.
pub suberror KrakenLiteError {
  KrakenLiteError(String)
}

///|
/// 分类树节点.
pub(all) struct KrakenTaxon {
  /// 分类单元编号
  taxid : String
  /// 学名
  name : String
  /// 分类层级(root/kingdom/phylum/class/order/family/genus/species)
  rank : String
  /// 父节点编号(根节点为 None)
  parent : String?
}

///|
/// 构造分类树节点.
pub fn kraken_taxon_new(
  taxid : String,
  name : String,
  rank : String,
  parent : String?,
) -> KrakenTaxon {
  KrakenTaxon::{ taxid, name, rank, parent, }
}

///|
/// 参考基因组条目.
pub(all) struct KrakenReference {
  /// 该基因组归属的分类单元
  taxid : String
  /// 基因组序列
  sequence : String
}

///|
/// 构造参考基因组条目.
pub fn kraken_reference_new(
  taxid : String,
  sequence : String,
) -> KrakenReference {
  KrakenReference::{ taxid, sequence, }
}

///|
/// 构建好的 k-mer 分类数据库.
pub(all) struct KrakenDatabase {
  /// k-mer 长度
  k : Int32
  /// 根节点编号
  root : String
  /// 分类树
  taxa : Map[String, KrakenTaxon]
  /// k-mer -> LCA 编号
  kmer_lca : Map[String, String]
  /// 分类单元 -> 直接分配给该节点的去重 k-mer 数
  direct_kmers : Map[String, Int32]
  /// 分类单元 -> 该节点子树内的去重 k-mer 数
  subtree_kmers : Map[String, Int32]
}

///|
/// 单条读段的分类结果.
pub(all) struct KrakenClassification {
  /// 读段名称
  name : String
  /// 分配到的分类单元(无法分类为 None)
  taxid : String?
  /// 分配到的层级
  rank : String
  /// 最佳路径上的 k-mer 支持数
  hits : Int32
  /// 读段中有效 k-mer 总数
  total_kmers : Int32
  /// 支持比例(hits / total_kmers)
  confidence : Double
}

///|
/// 单个分类单元的丰度.
pub(all) struct KrakenAbundance {
  /// 分类单元编号
  taxid : String
  /// 学名
  name : String
  /// 分类层级
  rank : String
  /// 计数(可为重分摊后的整数或小数, 统一以 Double 表示)
  count : Double
  /// 相对丰度(count / 已分类读段数)
  relative_abundance : Double
}

///|
/// 碱基归一化为大写; 非 A/C/G/T 返回空字符.
fn kraken_normalize_base(c : Char) -> Char {
  match c {
    'A' | 'a' => 'A'
    'C' | 'c' => 'C'
    'G' | 'g' => 'G'
    'T' | 't' => 'T'
    _ => '\u{0}'
  }
}

///|
/// 提取序列中全部只含 A/C/G/T 的 k-mer(统一大写).
fn kraken_extract_kmers(sequence : String, k : Int32) -> Array[String] {
  let out : Array[String] = []
  let n = sequence.length()
  if k <= 0 || n < k {
    return out
  }
  for start = 0; start < n - k + 1; start = start + 1 {
    let sb = StringBuilder(size_hint=k)
    let mut valid = true
    for offset = 0; offset < k; offset = offset + 1 {
      let base = kraken_normalize_base(
        sequence.unsafe_get(start + offset).unsafe_to_char(),
      )
      if base == '\u{0}' {
        valid = false
        break
      }
      sb.write_char(base)
    }
    if valid {
      out.push(sb.to_string())
    }
  }
  out
}

///|
/// 从节点到根的路径, 返回根在前的编号数组.
fn kraken_path(
  taxa : Map[String, KrakenTaxon],
  taxid : String,
) -> Array[String] {
  let upward : Array[String] = []
  let mut current : String? = Some(taxid)
  let mut steps = 0
  while steps <= taxa.length() {
    match current {
      Some(id) => {
        upward.push(id)
        match taxa.get(id) {
          Some(taxon) => current = taxon.parent
          None => current = None
        }
      }
      None => break
    }
    steps = steps + 1
  }
  let path : Array[String] = Array(capacity=upward.length())
  for i = 0; i < upward.length(); i = i + 1 {
    path.push(upward[upward.length() - 1 - i])
  }
  path
}

///|
/// 多个分类单元的最小共同祖先.
fn kraken_lca(
  taxa : Map[String, KrakenTaxon],
  hits : Array[String],
  root : String,
) -> String {
  if hits.length() == 0 {
    return root
  }
  let base = kraken_path(taxa, hits[0])
  let mut depth = base.length()
  for i = 1; i < hits.length(); i = i + 1 {
    let other = kraken_path(taxa, hits[i])
    let mut common = 0
    while common < depth &&
          common < other.length() &&
          base[common] == other[common] {
      common = common + 1
    }
    depth = common
  }
  if depth == 0 {
    root
  } else {
    base[depth - 1]
  }
}

///|
/// 由分类树和参考基因组构建 k-mer LCA 数据库.
pub fn kraken_build_database(
  k : Int32,
  taxa : Array[KrakenTaxon],
  references : Array[KrakenReference],
) -> KrakenDatabase raise KrakenLiteError {
  if k <= 0 {
    raise KrakenLiteError("k-mer length must be positive")
  }
  let taxon_map : Map[String, KrakenTaxon] = Map([], capacity=taxa.length())
  for taxon in taxa {
    if taxon.taxid == "" {
      raise KrakenLiteError("taxon id must not be empty")
    }
    if taxon_map.contains(taxon.taxid) {
      raise KrakenLiteError("duplicate taxon id: \{taxon.taxid}")
    }
    taxon_map.set(taxon.taxid, taxon)
  }
  let mut root = ""
  for taxon in taxa {
    match taxon.parent {
      Some(parent) =>
        if !taxon_map.contains(parent) {
          raise KrakenLiteError(
            "unknown parent \{parent} for taxon \{taxon.taxid}",
          )
        }
      None =>
        if root != "" {
          raise KrakenLiteError("multiple root taxa: \{root}, \{taxon.taxid}")
        } else {
          root = taxon.taxid
        }
    }
  }
  if root == "" {
    raise KrakenLiteError("taxonomy has no root")
  }

  let kmer_hits : Map[String, Array[String]] = Map([], capacity=0)
  for reference in references {
    if !taxon_map.contains(reference.taxid) {
      raise KrakenLiteError("reference has unknown taxon: \{reference.taxid}")
    }
    let seen : Map[String, Int32] = Map([], capacity=0)
    for kmer in kraken_extract_kmers(reference.sequence, k) {
      if !seen.contains(kmer) {
        seen.set(kmer, 1)
        match kmer_hits.get(kmer) {
          Some(hits) => hits.push(reference.taxid)
          None => kmer_hits.set(kmer, [reference.taxid])
        }
      }
    }
  }

  let kmer_lca : Map[String, String] = Map([], capacity=kmer_hits.length())
  let direct_kmers : Map[String, Int32] = Map([], capacity=taxa.length())
  for kmer in kmer_hits.keys() {
    let lca = kraken_lca(taxon_map, kmer_hits.get(kmer).unwrap(), root)
    kmer_lca.set(kmer, lca)
    match direct_kmers.get(lca) {
      Some(count) => direct_kmers.set(lca, count + 1)
      None => direct_kmers.set(lca, 1)
    }
  }
  let subtree_kmers : Map[String, Int32] = Map([], capacity=taxa.length())
  for taxid in direct_kmers.keys() {
    let count = direct_kmers.get(taxid).unwrap()
    for ancestor in kraken_path(taxon_map, taxid) {
      match subtree_kmers.get(ancestor) {
        Some(total) => subtree_kmers.set(ancestor, total + count)
        None => subtree_kmers.set(ancestor, count)
      }
    }
  }

  KrakenDatabase::{
    k,
    root,
    taxa: taxon_map,
    kmer_lca,
    direct_kmers,
    subtree_kmers,
  }
}

///|
/// 查询某 k-mer 的 LCA 分类单元.
pub fn kraken_lookup_kmer(db : KrakenDatabase, kmer : String) -> String? {
  match db.kmer_lca.get(kmer) {
    Some(taxid) => Some(taxid)
    None => None
  }
}

///|
/// 某分类单元子树内的去重 k-mer 数(作为 Bracken 权重).
pub fn kraken_kmer_support(db : KrakenDatabase, taxid : String) -> Int32 {
  match db.subtree_kmers.get(taxid) {
    Some(count) => count
    None => 0
  }
}

///|
/// 分类单元的层级; 未知节点返回空串.
pub fn kraken_taxon_rank(db : KrakenDatabase, taxid : String) -> String {
  match db.taxa.get(taxid) {
    Some(taxon) => taxon.rank
    None => ""
  }
}

///|
/// 子节点编号, 按编号升序以保证确定性.
fn kraken_children(db : KrakenDatabase, taxid : String) -> Array[String] {
  let out : Array[String] = []
  for child in db.taxa.keys() {
    match db.taxa.get(child).unwrap().parent {
      Some(parent) => if parent == taxid { out.push(child) }
      None => ()
    }
  }
  out.sort()
  out
}

///|
/// 分类树中该节点下方指定层级的全部后代编号.
fn kraken_descendants_at_rank(
  db : KrakenDatabase,
  taxid : String,
  level : String,
) -> Array[String] {
  let out : Array[String] = []
  for node in db.taxa.keys() {
    if kraken_taxon_rank(db, node) == level {
      let path = kraken_path(db.taxa, node)
      for ancestor in path {
        if ancestor == taxid {
          out.push(node)
          break
        }
      }
    }
  }
  out.sort()
  out
}

///|
/// 从该节点向上寻找指定层级的祖先(含自身).
fn kraken_ancestor_at_rank(
  db : KrakenDatabase,
  taxid : String,
  level : String,
) -> String? {
  let path = kraken_path(db.taxa, taxid)
  for i = path.length() - 1; i >= 0; i = i - 1 {
    if kraken_taxon_rank(db, path[i]) == level {
      return Some(path[i])
    }
  }
  None
}

///|
/// 用 k-mer LCA 数据库对单条读段分类.
///
/// 先从命中 k-mer 的 LCA 聚合出每个分类单元的子树支持度, 再自根向下
/// 选择支持度最大的子节点; 支持度为 0 时停止. 同分时按分类单元编号
/// 升序选择, 保证结果确定性.
pub fn kraken_classify(
  db : KrakenDatabase,
  name : String,
  sequence : String,
) -> KrakenClassification {
  let kmers = kraken_extract_kmers(sequence, db.k)
  let total_kmers = kmers.length()
  if total_kmers == 0 {
    return KrakenClassification::{
      name,
      taxid: None,
      rank: "unclassified",
      hits: 0,
      total_kmers,
      confidence: 0.0,
    }
  }

  let direct : Map[String, Int32] = Map([], capacity=taxon_length(db))
  let mut valid_hits = 0
  for kmer in kmers {
    match db.kmer_lca.get(kmer) {
      Some(taxid) => {
        valid_hits = valid_hits + 1
        match direct.get(taxid) {
          Some(count) => direct.set(taxid, count + 1)
          None => direct.set(taxid, 1)
        }
      }
      None => ()
    }
  }
  if valid_hits == 0 {
    return KrakenClassification::{
      name,
      taxid: None,
      rank: "unclassified",
      hits: 0,
      total_kmers,
      confidence: 0.0,
    }
  }

  let subtree : Map[String, Int32] = Map([], capacity=taxon_length(db))
  for taxid in direct.keys() {
    let count = direct.get(taxid).unwrap()
    for ancestor in kraken_path(db.taxa, taxid) {
      match subtree.get(ancestor) {
        Some(total) => subtree.set(ancestor, total + count)
        None => subtree.set(ancestor, count)
      }
    }
  }

  let direct_at = fn(taxid : String) -> Int32 {
    match direct.get(taxid) {
      Some(count) => count
      None => 0
    }
  }

  let mut best = db.root
  let mut path_hits = direct_at(best)
  let mut running = true
  while running {
    let children = kraken_children(db, best)
    let mut chosen = ""
    let mut chosen_support = 0
    for child in children {
      let support = match subtree.get(child) {
        Some(value) => value
        None => 0
      }
      if support > chosen_support {
        chosen = child
        chosen_support = support
      }
    }
    if chosen == "" {
      running = false
    } else {
      best = chosen
      path_hits = path_hits + direct_at(best)
    }
  }

  KrakenClassification::{
    name,
    taxid: Some(best),
    rank: kraken_taxon_rank(db, best),
    hits: path_hits,
    total_kmers,
    confidence: path_hits.to_double() / total_kmers.to_double(),
  }
}

///|
/// 批量分类(读段名称 -> 序列).
pub fn kraken_classify_many(
  db : KrakenDatabase,
  reads : Map[String, String],
) -> Array[KrakenClassification] {
  let names : Array[String] = []
  for name in reads.keys() {
    names.push(name)
  }
  names.sort()
  let out : Array[KrakenClassification] = Array(capacity=names.length())
  for name in names {
    out.push(kraken_classify(db, name, reads.get(name).unwrap()))
  }
  out
}

///|
/// 分类树节点总数.
fn taxon_length(db : KrakenDatabase) -> Int32 {
  db.taxa.length()
}

///|
/// 计数排序比较器: 计数降序, 同分按编号升序.
fn kraken_abundance_compare(
  left : KrakenAbundance,
  right : KrakenAbundance,
) -> Int32 {
  if left.count > right.count {
    -1
  } else if left.count < right.count {
    1
  } else if left.taxid < right.taxid {
    -1
  } else if left.taxid > right.taxid {
    1
  } else {
    0
  }
}

///|
/// Kraken 风格分类报告: 统计每个分类单元自身及其后代的读段数.
pub fn kraken_report(
  db : KrakenDatabase,
  classifications : Array[KrakenClassification],
  level : String,
) -> Array[KrakenAbundance] {
  let cumulative : Map[String, Double] = Map([], capacity=taxon_length(db))
  let mut total_classified = 0
  for result in classifications {
    match result.taxid {
      Some(taxid) => {
        total_classified = total_classified + 1
        for ancestor in kraken_path(db.taxa, taxid) {
          add_abundance(cumulative, ancestor, 1.0)
        }
      }
      None => ()
    }
  }
  build_abundance(db, cumulative, total_classified, level)
}

///|
/// Bracken 风格重分摊: 将高层级(如属)的读段计数按数据库 k-mer 支持度
/// 分配到低层级(如种).
pub fn kraken_bracken(
  db : KrakenDatabase,
  classifications : Array[KrakenClassification],
  level : String,
) -> Array[KrakenAbundance] {
  let assigned : Map[String, Double] = Map([], capacity=taxon_length(db))
  let mut total_classified = 0
  for result in classifications {
    match result.taxid {
      Some(taxid) => {
        total_classified = total_classified + 1
        if kraken_taxon_rank(db, taxid) == level {
          add_abundance(assigned, taxid, 1.0)
        } else {
          let descendants = kraken_descendants_at_rank(db, taxid, level)
          if descendants.length() > 0 {
            distribute_to_descendants(db, descendants, assigned)
          } else {
            match kraken_ancestor_at_rank(db, taxid, level) {
              Some(ancestor) => add_abundance(assigned, ancestor, 1.0)
              None => ()
            }
          }
        }
      }
      None => ()
    }
  }
  build_abundance(db, assigned, total_classified, level)
}

///|
/// 将 1 条读段的计数按子树 k-mer 支持度分摊给后代.
fn distribute_to_descendants(
  db : KrakenDatabase,
  descendants : Array[String],
  assigned : Map[String, Double],
) -> Unit {
  let weights : Map[String, Int32] = Map([], capacity=descendants.length())
  let mut total_weight = 0
  for taxid in descendants {
    let support = kraken_kmer_support(db, taxid)
    let weight = if support > 0 { support } else { 1 }
    weights.set(taxid, weight)
    total_weight = total_weight + weight
  }
  if total_weight <= 0 {
    return
  }
  for taxid in descendants {
    let share = weights.get(taxid).unwrap().to_double() /
      total_weight.to_double()
    add_abundance(assigned, taxid, share)
  }
}

///|
/// 累加丰度计数.
fn add_abundance(
  assigned : Map[String, Double],
  taxid : String,
  amount : Double,
) -> Unit {
  match assigned.get(taxid) {
    Some(current) => assigned.set(taxid, current + amount)
    None => assigned.set(taxid, amount)
  }
}

///|
/// 由分类单元计数构建排序后的丰度表.
fn build_abundance(
  db : KrakenDatabase,
  counts : Map[String, Double],
  total_classified : Int32,
  level : String,
) -> Array[KrakenAbundance] {
  let out : Array[KrakenAbundance] = []
  for taxid in counts.keys() {
    let count = counts.get(taxid).unwrap()
    if count > 0.0 && kraken_taxon_rank(db, taxid) == level {
      let taxon = db.taxa.get(taxid).unwrap()
      out.push(KrakenAbundance::{
        taxid,
        name: taxon.name,
        rank: taxon.rank,
        count,
        relative_abundance: if total_classified > 0 {
          count / total_classified.to_double()
        } else {
          0.0
        },
      })
    }
  }
  out.sort_by(kraken_abundance_compare)
  out
}

///|
/// 丰度表 -> JSON 文本.
pub fn kraken_abundance_json(
  db : KrakenDatabase,
  abundances : Array[KrakenAbundance],
) -> String {
  let taxa : Array[JsonValue] = Array(capacity=abundances.length())
  for entry in abundances {
    taxa.push(
      JsonObj(
        Map([
          ("taxid", JsonStr(entry.taxid)),
          ("name", JsonStr(entry.name)),
          ("rank", JsonStr(entry.rank)),
          ("count", JsonNum(entry.count)),
          ("relative_abundance", JsonNum(entry.relative_abundance)),
        ]),
      ),
    )
  }
  json_to_string(
    JsonObj(Map([("k", JsonNum(db.k.to_double())), ("taxa", JsonArr(taxa))])),
  )
}

///|
/// 丰度表 -> Kraken 风格 TSV 文本.
pub fn kraken_abundance_tsv(abundances : Array[KrakenAbundance]) -> String {
  let sb = StringBuilder()
  sb.write_string("taxid\tname\trank\tcount\trelative_abundance\n")
  for entry in abundances {
    sb.write_string(entry.taxid)
    sb.write_char('\t')
    sb.write_string(entry.name)
    sb.write_char('\t')
    sb.write_string(entry.rank)
    sb.write_char('\t')
    sb.write_string(entry.count.to_string())
    sb.write_char('\t')
    sb.write_string(entry.relative_abundance.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}
