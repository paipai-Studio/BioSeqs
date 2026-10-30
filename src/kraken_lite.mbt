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
  /// minimizer 窗口大小; 0 表示精确 k-mer 模式(建库与分类均用全部 k-mer).
  minimizer_w : Int32
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
/// 碱基互补(仅 A/C/G/T).
fn kraken_complement_base(c : Char) -> Char {
  match c {
    'A' => 'T'
    'C' => 'G'
    'G' => 'C'
    'T' => 'A'
    _ => '\u{0}'
  }
}

///|
/// 反向互补(非 A/C/G/T 视为 N).
fn kraken_revcomp(seq : String) -> String {
  let sb = StringBuilder(size_hint=seq.length())
  let mut i = seq.length() - 1
  while i >= 0 {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let comp = kraken_complement_base(c)
    sb.write_char(if comp == '\u{0}' { 'N' } else { comp })
    i = i - 1
  }
  sb.to_string()
}

///|
/// 规范 k-mer: 正链与其反向互补中字典序较小者, 使索引与链方向无关.
fn kraken_canonical_kmer(kmer : String) -> String {
  let rc = kraken_revcomp(kmer)
  if kmer <= rc {
    kmer
  } else {
    rc
  }
}

///|
/// k-mer 的确定性散列(FNV-1a 64 位), 用于 minimizer 选择.
fn kraken_kmer_hash(kmer : String) -> UInt64 {
  let mut h : UInt64 = 0xcbf29ce484222325UL
  for i = 0; i < kmer.length(); i = i + 1 {
    h = h ^ kmer.unsafe_get(i).unsafe_to_char().to_int().to_uint64()
    h = h * 0x100000001b3UL
  }
  h
}

///|
/// 提取全部规范 k-mer(不判重), 仅保留只含 A/C/G/T 的窗口.
fn kraken_extract_canonical_kmers(
  sequence : String,
  k : Int32,
) -> Array[String] {
  let out : Array[String] = []
  for kmer in kraken_extract_kmers(sequence, k) {
    out.push(kraken_canonical_kmer(kmer))
  }
  out
}

///|
/// 提取序列的 minimizer 集合(规范 k-mer 在窗口 w 内取散列最小者, 同散列取字典序小者).
///
/// minimizer 方案只保留每 w 个连续 k-mer 窗口中的一个代表, 显著缩小
/// 数据库体积; 分类时对读段提取同一批 minimizer 回查索引.
fn kraken_extract_minimizers(
  sequence : String,
  k : Int32,
  w : Int32,
) -> Array[String] {
  let canon = kraken_extract_canonical_kmers(sequence, k)
  let n = canon.length()
  let out : Array[String] = []
  let seen : Map[String, Bool] = Map([], capacity=n)
  if n == 0 {
    return out
  }
  let window = if w < 1 { 1 } else { w }
  if n < window {
    for kmer in canon {
      if !seen.contains(kmer) {
        seen.set(kmer, true)
        out.push(kmer)
      }
    }
    return out
  }
  for start = 0; start < n - window + 1; start = start + 1 {
    let mut best = canon[start]
    let mut best_hash = kraken_kmer_hash(best)
    for offset = 1; offset < window; offset = offset + 1 {
      let candidate = canon[start + offset]
      let h = kraken_kmer_hash(candidate)
      if h < best_hash || (h == best_hash && candidate < best) {
        best = candidate
        best_hash = h
      }
    }
    if !seen.contains(best) {
      seen.set(best, true)
      out.push(best)
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
///
/// `w = 0` 使用全部 k-mer; `w > 0` 使用窗口为 w 的 minimizer 建库,
/// 数据库体积显著缩小, 分类时读段同样按 minimizer 回查.
fn kraken_build_inner(
  k : Int32,
  w : Int32,
  taxa : Array[KrakenTaxon],
  references : Array[KrakenReference],
) -> KrakenDatabase raise KrakenLiteError {
  if k <= 0 {
    raise KrakenLiteError("k-mer length must be positive")
  }
  if w < 0 {
    raise KrakenLiteError("minimizer window must be non-negative")
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
    let kmers = if w > 0 {
      kraken_extract_minimizers(reference.sequence, k, w)
    } else {
      kraken_extract_kmers(reference.sequence, k)
    }
    for kmer in kmers {
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
    minimizer_w: w,
  }
}

///|
/// 精确 k-mer 建库(k-mer -> LCA), 保留全部 k-mer.
pub fn kraken_build_database(
  k : Int32,
  taxa : Array[KrakenTaxon],
  references : Array[KrakenReference],
) -> KrakenDatabase raise KrakenLiteError {
  kraken_build_inner(k, 0, taxa, references)
}

///|
/// minimizer 建库: 每个长度 w 的连续 k-mer 窗口只保留散列最小的规范 k-mer.
///
/// 相比精确建库可显著压缩数据库规模(内存开销大的场景), 代价是分类召回率
/// 略降; `w = 1` 等价于精确建库(只保留规范 k-mer).
pub fn kraken_build_database_minimizer(
  k : Int32,
  w : Int32,
  taxa : Array[KrakenTaxon],
  references : Array[KrakenReference],
) -> KrakenDatabase raise KrakenLiteError {
  if w < 1 {
    raise KrakenLiteError("minimizer window must be at least 1")
  }
  kraken_build_inner(k, w, taxa, references)
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
  let kmers = if db.minimizer_w > 0 {
    kraken_extract_minimizers(sequence, db.k, db.minimizer_w)
  } else {
    kraken_extract_kmers(sequence, db.k)
  }
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

// ============================================================================
// 数据库 JSON 序列化 / 反序列化(内存数据模型 + JSON 交换)
// ============================================================================

///|
/// 由成对数组构造 JSON 对象, 保持插入顺序语义(值以 Map 存储, 输出按 key 排序).
fn kraken_json_obj(pairs : Array[(String, JsonValue)]) -> JsonValue {
  let m : Map[String, JsonValue] = Map([], capacity=pairs.length())
  for pair in pairs {
    m.set(pair.0, pair.1)
  }
  JsonObj(m)
}

///|
/// 数据库 -> JSON 值(排序键以保证确定性输出).
pub fn kraken_database_value(db : KrakenDatabase) -> JsonValue {
  let taxids : Array[String] = []
  for id in db.taxa.keys() {
    taxids.push(id)
  }
  taxids.sort()
  let taxa : Array[JsonValue] = Array(capacity=taxids.length())
  for id in taxids {
    let taxon = db.taxa.get(id).unwrap()
    let parent : JsonValue = match taxon.parent {
      Some(p) => JsonStr(p)
      None => JsonNull
    }
    taxa.push(
      kraken_json_obj([
        ("taxid", JsonStr(taxon.taxid)),
        ("name", JsonStr(taxon.name)),
        ("rank", JsonStr(taxon.rank)),
        ("parent", parent),
      ]),
    )
  }
  let kmer_keys : Array[String] = []
  for kmer in db.kmer_lca.keys() {
    kmer_keys.push(kmer)
  }
  kmer_keys.sort()
  let kmer_lca : Array[JsonValue] = Array(capacity=kmer_keys.length())
  for kmer in kmer_keys {
    kmer_lca.push(
      JsonArr([JsonStr(kmer), JsonStr(db.kmer_lca.get(kmer).unwrap())]),
    )
  }
  let direct = kraken_count_pairs(db.direct_kmers)
  let subtree = kraken_count_pairs(db.subtree_kmers)
  kraken_json_obj([
    ("k", JsonNum(db.k.to_double())),
    ("minimizer_w", JsonNum(db.minimizer_w.to_double())),
    ("root", JsonStr(db.root)),
    ("taxa", JsonArr(taxa)),
    ("kmer_lca", JsonArr(kmer_lca)),
    ("direct_kmers", JsonArr(direct)),
    ("subtree_kmers", JsonArr(subtree)),
  ])
}

///|
/// 计数映射 -> 排序后的 [id, count] 数组.
fn kraken_count_pairs(counts : Map[String, Int32]) -> Array[JsonValue] {
  let keys : Array[String] = []
  for id in counts.keys() {
    keys.push(id)
  }
  keys.sort()
  let out : Array[JsonValue] = Array(capacity=keys.length())
  for id in keys {
    out.push(
      JsonArr([JsonStr(id), JsonNum(counts.get(id).unwrap().to_double())]),
    )
  }
  out
}

///|
/// 数据库 -> JSON 文本.
pub fn kraken_database_to_json(db : KrakenDatabase) -> String {
  json_to_string(kraken_database_value(db))
}

///|
/// JSON 数值 -> Int32.
fn kraken_json_to_int(value : JsonValue) -> Int32? {
  match value {
    JsonNum(x) => Some(x.to_int())
    _ => None
  }
}

///|
/// JSON 文本 -> 数据库(roundtrip).
pub fn kraken_database_from_json(text : String) -> KrakenDatabase? {
  let parsed = parse_json(text)
  match parsed {
    Some(JsonObj(fields)) => {
      let mut k = 0
      let mut minimizer_w = 0
      let mut root = ""
      let taxa : Map[String, KrakenTaxon] = Map([], capacity=0)
      let kmer_lca : Map[String, String] = Map([], capacity=0)
      let direct_kmers : Map[String, Int32] = Map([], capacity=0)
      let subtree_kmers : Map[String, Int32] = Map([], capacity=0)
      for field in fields {
        if field.0 == "k" {
          match kraken_json_to_int(field.1) {
            Some(x) => k = x
            None => ()
          }
        } else if field.0 == "minimizer_w" {
          match kraken_json_to_int(field.1) {
            Some(x) => minimizer_w = x
            None => ()
          }
        } else if field.0 == "root" {
          match field.1 {
            JsonStr(s) => root = s
            _ => ()
          }
        } else if field.0 == "taxa" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match kraken_parse_taxon_json(entry) {
                  Some(taxon) => taxa.set(taxon.taxid, taxon)
                  None => ()
                }
              }
            _ => ()
          }
        } else if field.0 == "kmer_lca" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match kraken_parse_str_pair(entry) {
                  Some((key, value)) => kmer_lca.set(key, value)
                  None => ()
                }
              }
            _ => ()
          }
        } else if field.0 == "direct_kmers" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match kraken_parse_int_pair(entry) {
                  Some((key, value)) => direct_kmers.set(key, value)
                  None => ()
                }
              }
            _ => ()
          }
        } else if field.0 == "subtree_kmers" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match kraken_parse_int_pair(entry) {
                  Some((key, value)) => subtree_kmers.set(key, value)
                  None => ()
                }
              }
            _ => ()
          }
        }
      }
      if k <= 0 || root == "" || taxa.length() == 0 {
        return None
      }
      Some(KrakenDatabase::{
        k,
        root,
        taxa,
        kmer_lca,
        direct_kmers,
        subtree_kmers,
        minimizer_w,
      })
    }
    _ => None
  }
}

///|
/// 解析单个分类树节点 JSON.
fn kraken_parse_taxon_json(value : JsonValue) -> KrakenTaxon? {
  match value {
    JsonObj(fields) => {
      let mut taxid = ""
      let mut name = ""
      let mut rank = ""
      let mut parent : String? = None
      for field in fields {
        if field.0 == "taxid" {
          match field.1 {
            JsonStr(s) => taxid = s
            _ => ()
          }
        } else if field.0 == "name" {
          match field.1 {
            JsonStr(s) => name = s
            _ => ()
          }
        } else if field.0 == "rank" {
          match field.1 {
            JsonStr(s) => rank = s
            _ => ()
          }
        } else if field.0 == "parent" {
          match field.1 {
            JsonStr(s) => parent = Some(s)
            JsonNull => parent = None
            _ => ()
          }
        }
      }
      if taxid == "" {
        None
      } else {
        Some(KrakenTaxon::{ taxid, name, rank, parent, })
      }
    }
    _ => None
  }
}

///|
/// 解析 [字符串, 字符串] 对.
fn kraken_parse_str_pair(value : JsonValue) -> (String, String)? {
  match value {
    JsonArr(items) =>
      if items.length() == 2 {
        match (items[0], items[1]) {
          (JsonStr(a), JsonStr(b)) => Some((a, b))
          _ => None
        }
      } else {
        None
      }
    _ => None
  }
}

///|
/// 解析 [字符串, 整数] 对.
fn kraken_parse_int_pair(value : JsonValue) -> (String, Int32)? {
  match value {
    JsonArr(items) =>
      if items.length() == 2 {
        match (items[0], kraken_json_to_int(items[1])) {
          (JsonStr(a), Some(b)) => Some((a, b))
          _ => None
        }
      } else {
        None
      }
    _ => None
  }
}

// ============================================================================
// Kraken 风格 HTML/SVG 报告(衔接 MultiQC 呈现链路)
// ============================================================================

///|
fn kraken_html_escape(s : String) -> String {
  s
  .replace(old="&", new="&amp;")
  .replace(old="<", new="&lt")
  .replace(old=">", new="&gt")
}

///|
/// 单页 HTML 报告: 丰度条形图(SVG) + 明细表.
///
/// 完全本地生成, 不引用外部资源, 保持"数据不出本机".
pub fn kraken_report_html(
  abundances : Array[KrakenAbundance],
  title~ : String = "Kraken-lite Report",
) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<title>",
  )
  sb.write_string(kraken_html_escape(title))
  sb.write_string("</title>\n<style>\n")
  sb.write_string(
    "body{font-family:-apple-system,'Segoe UI',Roboto,sans-serif;margin:0;background:#f5f6fa;color:#2d3436;}\n",
  )
  sb.write_string(
    ".wrap{max-width:1000px;margin:0 auto;padding:24px 20px 40px;}\n",
  )
  sb.write_string("h1{font-size:20px;margin:0 0 16px;}\n")
  sb.write_string(
    "table{border-collapse:collapse;width:100%;background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08);}\n",
  )
  sb.write_string("th,td{padding:8px 12px;font-size:13px;text-align:right;}\n")
  sb.write_string(
    "th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left;}\n",
  )
  sb.write_string("th{background:#dfe6e9;font-weight:600;}\n")
  sb.write_string("</style>\n</head>\n<body>\n<div class=\"wrap\">\n")
  sb.write_string("<h1>")
  sb.write_string(kraken_html_escape(title))
  sb.write_string("</h1>\n")

  let bar_width = 560.0
  let label_width = 240.0
  let row_height = 26.0
  let top_pad = 10.0
  let mut max_rel = 0.0
  for entry in abundances {
    if entry.relative_abundance > max_rel {
      max_rel = entry.relative_abundance
    }
  }
  let chart_h = top_pad * 2.0 + row_height * abundances.length().to_double()
  sb.write_string("<svg width=\"840\" height=\"")
  sb.write_string(chart_h.to_int().to_string())
  sb.write_string("\" viewBox=\"0 0 840 ")
  sb.write_string(chart_h.to_int().to_string())
  sb.write_string("\" xmlns=\"http://www.w3.org/2000/svg\">\n")
  for i, entry in abundances {
    let y = top_pad + row_height * i.to_double()
    let bar_h = 16.0
    let ratio = if max_rel > 0.0 {
      entry.relative_abundance / max_rel
    } else {
      0.0
    }
    let w = bar_width * ratio
    let label = entry.name + " (" + entry.taxid + ")"
    sb.write_string("<text x=\"0\" y=\"")
    sb.write_string((y + 13.0).to_int().to_string())
    sb.write_string("\" font-size=\"12\" fill=\"#2d3436\">")
    sb.write_string(kraken_html_escape(label))
    sb.write_string("</text>\n")
    sb.write_string("<rect x=\"")
    sb.write_string(label_width.to_int().to_string())
    sb.write_string("\" y=\"")
    sb.write_string(y.to_int().to_string())
    sb.write_string("\" width=\"")
    sb.write_string(w.to_int().to_string())
    sb.write_string("\" height=\"")
    sb.write_string(bar_h.to_int().to_string())
    sb.write_string("\" fill=\"#0984e3\" rx=\"2\"/>\n")
    sb.write_string("<text x=\"")
    sb.write_string((label_width + w + 6.0).to_int().to_string())
    sb.write_string("\" y=\"")
    sb.write_string((y + 13.0).to_int().to_string())
    sb.write_string("\" font-size=\"11\" fill=\"#636e72\">")
    sb.write_string(kraken_format_percent(entry.relative_abundance))
    sb.write_string("</text>\n")
  }
  sb.write_string("</svg>\n")

  sb.write_string(
    "<table>\n<tr><th>taxon</th><th>taxid</th><th>rank</th><th>count</th><th>relative</th></tr>\n",
  )
  for entry in abundances {
    sb.write_string("<tr><td>")
    sb.write_string(kraken_html_escape(entry.name))
    sb.write_string("</td><td>")
    sb.write_string(kraken_html_escape(entry.taxid))
    sb.write_string("</td><td>")
    sb.write_string(kraken_html_escape(entry.rank))
    sb.write_string("</td><td>")
    sb.write_string(kraken_format_count(entry.count))
    sb.write_string("</td><td>")
    sb.write_string(kraken_format_percent(entry.relative_abundance))
    sb.write_string("</td></tr>\n")
  }
  sb.write_string("</table>\n")
  sb.write_string(
    "<p style=\"font-size:12px;color:#95a5a6;\">generated by BioSeqs kraken-lite (local, no data leaves this machine)</p>\n",
  )
  sb.write_string("</div>\n</body>\n</html>\n")
  sb.to_string()
}

///|
/// 相对丰度格式化百分比(保留两位小数).
fn kraken_format_percent(value : Double) -> String {
  let rounded = (value * 10000.0).round() / 100.0
  rounded.to_string() + "%"
}

///|
/// 计数格式化(整数显示为整数, 小数保留两位).
fn kraken_format_count(value : Double) -> String {
  if value == value.trunc() && value > -1.0e15 && value < 1.0e15 {
    value.to_int().to_string()
  } else {
    ((value * 100.0).round() / 100.0).to_string()
  }
}
