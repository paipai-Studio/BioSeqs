///|
/// FastQC 风格的 FASTQ 质控指标计算, 输出 JSON 报告.
///
/// 供 `examples/fq_qc_web` 浏览器仪表盘消费; 亦可独立用于流水线.
/// 覆盖 FastQC 核心模块: 基础统计/逐位质量(直方图精确五数)/
/// 序列质量分布/逐位 GC/序列 GC 分布/读长分布/adapter 污染/
/// 高频重复序列.
///
/// 流式逐条处理(复用 fastq_each_record), 逐位质量按 0-50 直方图聚合,
/// 内存 O(位数 x 51 + 高频序列), 与 reads 数无关.
pub let fqqc_max_qual : Int32 = 50

///|
/// 常见 Illumina adapter/污染前缀.
pub let fqqc_adapter_prefixes : Array[String] = [
  "AGATCGGAAGAGC", "CTGTCTCTTATA", "GGGGGGGGGGGG", "TTTTTTTTTTTT", "AATGATACGGCG",
  "CAAGCAGAAGAC",
]

///|
/// 质控统计聚合器.
pub struct FqqcStats {
  mut n_reads : Int32
  mut total_bases : Int32
  lengths : Map[Int32, Int32]
  mean_quals : Array[Int32]
  gc_percents : Array[Int32]
  per_base_hist : Array[Array[Int32]]
  per_base_count : Array[Int32]
  per_base_gc : Array[Int32]
  adapter_hits : Map[String, Int32]
  seq_counts : Map[String, Int32]
}

///|
/// 新建空的统计聚合器.
pub fn fqqc_new_stats() -> FqqcStats {
  {
    n_reads: 0,
    total_bases: 0,
    lengths: Map([], capacity=0),
    mean_quals: [],
    gc_percents: [],
    per_base_hist: [],
    per_base_count: [],
    per_base_gc: [],
    adapter_hits: Map([], capacity=4),
    seq_counts: Map([], capacity=0),
  }
}

///|
/// 流式计算全部质控指标.
pub fn fqqc_compute(fastq : String) -> FqqcStats raise SeqIOError {
  let st = fqqc_new_stats()
  fastq_each_record(fastq, fn(rec) { fqqc_visit(st, rec) })
  st
}

///|
fn fqqc_visit(st : FqqcStats, rec : SeqRecord) -> Unit {
  let s = rec.seq.to_string()
  st.n_reads = st.n_reads + 1
  st.total_bases = st.total_bases + s.length()
  match st.lengths.get(s.length()) {
    Some(v) => st.lengths[s.length()] = v + 1
    None => st.lengths[s.length()] = 1
  }
  match st.seq_counts.get(s) {
    Some(v) => st.seq_counts[s] = v + 1
    None => st.seq_counts[s] = 1
  }
  for prefix in fqqc_adapter_prefixes {
    if s.length() >= prefix.length() && s.contains(prefix) {
      match st.adapter_hits.get(prefix) {
        Some(v) => st.adapter_hits[prefix] = v + 1
        None => st.adapter_hits[prefix] = 1
      }
    }
  }
  let mut q_sum = 0
  let mut gc = 0
  for i in 0..<s.length() {
    let c = s.unsafe_get(i)
    let is_gc = c == 'G' || c == 'C' || c == 'g' || c == 'c'
    if is_gc {
      gc = gc + 1
    }
    // 扩展逐位聚合槽位
    if i >= st.per_base_hist.length() {
      st.per_base_hist.push(Array::make(fqqc_max_qual + 1, 0))
      st.per_base_count.push(0)
      st.per_base_gc.push(0)
    }
    if is_gc {
      st.per_base_gc[i] = st.per_base_gc[i] + 1
    }
  }
  match rec.letter_annotations.get("phred_quality") {
    Some(quals) =>
      for i in 0..<quals.length() {
        let q = quals[i]
        q_sum = q_sum + q
        if i < st.per_base_hist.length() {
          let bucket = if q > fqqc_max_qual {
            fqqc_max_qual
          } else if q < 0 {
            0
          } else {
            q
          }
          st.per_base_hist[i][bucket] = st.per_base_hist[i][bucket] + 1
          st.per_base_count[i] = st.per_base_count[i] + 1
        }
      }
    None => ()
  }
  let mean_q = if s.length() > 0 { q_sum / s.length() } else { 0 }
  st.mean_quals.push(mean_q)
  let gc_pct = if s.length() > 0 { gc * 100 / s.length() } else { 0 }
  st.gc_percents.push(gc_pct)
}

///|
/// 从直方图计算升序第 k 个值(k 0-based).
fn fqqc_hist_kth(hist : Array[Int32], count : Int32, k : Int32) -> Int32 {
  if count == 0 {
    return 0
  }
  let target = if k >= count { count - 1 } else if k < 0 { 0 } else { k }
  let mut acc = 0
  for v, n in hist {
    acc = acc + n
    if acc > target {
      return v
    }
  }
  0
}

///|
/// 直方图分位值(线性插值), p in [0, 100].
fn fqqc_hist_quantile(hist : Array[Int32], count : Int32, p : Double) -> Double {
  if count == 0 {
    return 0.0
  }
  let pos = p / 100.0 * (count - 1).to_double()
  let lo = pos.to_int()
  let hi = if lo + 1 < count { lo + 1 } else { lo }
  let frac = pos - lo.to_double()
  let vlo = fqqc_hist_kth(hist, count, lo).to_double()
  let vhi = fqqc_hist_kth(hist, count, hi).to_double()
  vlo + (vhi - vlo) * frac
}

///|
/// 值域 [0,100] 整数直方图, bins 桶.
fn fqqc_hist(values : Array[Int32], bins : Int32) -> Array[Int32] {
  let hist = Array::make(bins, 0)
  if bins <= 0 {
    return hist
  }
  for v in values {
    let idx = v * bins / 100
    let idx = if idx >= bins { bins - 1 } else { idx }
    if idx >= 0 {
      hist[idx] = hist[idx] + 1
    }
  }
  hist
}

///|
fn fqqc_obj(pairs : Array[(String, JsonValue)]) -> JsonValue {
  let m : Map[String, JsonValue] = Map([], capacity=pairs.length())
  for p in pairs {
    m[p.0] = p.1
  }
  JsonObj(m)
}

///|
fn fqqc_arr_int(arr : Array[Int32]) -> JsonValue {
  let out : Array[JsonValue] = Array(capacity=arr.length())
  for v in arr {
    out.push(JsonNum(v.to_double()))
  }
  JsonArr(out)
}

///|
fn fqqc_arr_num(arr : Array[Double]) -> JsonValue {
  let out : Array[JsonValue] = Array(capacity=arr.length())
  for v in arr {
    out.push(JsonNum(v))
  }
  JsonArr(out)
}

///|
/// 报告树(JsonValue), 便于嵌入或扩展.
pub fn fqqc_report_value(st : FqqcStats) -> JsonValue {
  // 1. 基础统计
  let mean_len = if st.n_reads > 0 {
    st.total_bases.to_double() / st.n_reads.to_double()
  } else {
    0.0
  }
  let mut total_gc = 0
  for g in st.gc_percents {
    total_gc = total_gc + g
  }
  let overall_gc = if st.n_reads > 0 {
    total_gc.to_double() / st.n_reads.to_double()
  } else {
    0.0
  }
  // 2. 逐位质量五数 + 均值(由直方图精确计算)
  let n_pos = st.per_base_count.length()
  let pos_labels : Array[Int32] = Array(capacity=n_pos)
  let lows : Array[Double] = Array(capacity=n_pos)
  let q1s : Array[Double] = Array(capacity=n_pos)
  let meds : Array[Double] = Array(capacity=n_pos)
  let q3s : Array[Double] = Array(capacity=n_pos)
  let highs : Array[Double] = Array(capacity=n_pos)
  let means : Array[Double] = Array(capacity=n_pos)
  let gc_by_pos : Array[Double] = Array(capacity=n_pos)
  for i in 0..<n_pos {
    let hist = st.per_base_hist[i]
    let count = st.per_base_count[i]
    let mut sum = 0
    for v, n in hist {
      sum = sum + v * n
    }
    pos_labels.push(i + 1)
    lows.push(fqqc_hist_quantile(hist, count, 0.0))
    q1s.push(fqqc_hist_quantile(hist, count, 25.0))
    meds.push(fqqc_hist_quantile(hist, count, 50.0))
    q3s.push(fqqc_hist_quantile(hist, count, 75.0))
    highs.push(fqqc_hist_quantile(hist, count, 100.0))
    means.push(
      if count > 0 {
        sum.to_double() / count.to_double()
      } else {
        0.0
      },
    )
    gc_by_pos.push(
      if st.n_reads > 0 {
        st.per_base_gc[i].to_double() / st.n_reads.to_double() * 100.0
      } else {
        0.0
      },
    )
  }
  // 3. 序列质量/GC 分布
  let qual_hist = fqqc_hist(st.mean_quals, 20)
  let gc_hist = fqqc_hist(st.gc_percents, 20)
  // 4. 读长分布(排序键)
  let len_keys : Array[Int32] = []
  for entry in st.lengths {
    len_keys.push(entry.0)
  }
  let mut i = 0
  while i < len_keys.length() {
    let mut j = i + 1
    while j < len_keys.length() {
      if len_keys[j] < len_keys[i] {
        let t = len_keys[i]
        len_keys[i] = len_keys[j]
        len_keys[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  let len_counts : Array[JsonValue] = Array(capacity=len_keys.length())
  for k in len_keys {
    len_counts.push(
      fqqc_obj([
        ("length", JsonNum(k.to_double())),
        ("count", JsonNum(st.lengths[k].to_double())),
      ]),
    )
  }
  // 5. adapter 命中
  let adapter_arr : Array[JsonValue] = Array(capacity=0)
  for prefix in fqqc_adapter_prefixes {
    adapter_arr.push(
      fqqc_obj([
        ("prefix", JsonStr(prefix)),
        (
          "reads",
          JsonNum(
            match st.adapter_hits.get(prefix) {
              Some(v) => v.to_double()
              None => 0.0
            },
          ),
        ),
      ]),
    )
  }
  // 6. 高频序列 top 20
  let pairs : Array[(String, Int32)] = Array(capacity=0)
  for entry in st.seq_counts {
    if entry.1 > 1 {
      pairs.push((entry.0, entry.1))
    }
  }
  let mut i = 0
  while i < pairs.length() {
    let mut j = i + 1
    while j < pairs.length() {
      if pairs[j].1 > pairs[i].1 {
        let t = pairs[i]
        pairs[i] = pairs[j]
        pairs[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  let top_n = if pairs.length() > 20 { 20 } else { pairs.length() }
  let overrep : Array[JsonValue] = Array(capacity=top_n)
  for k in 0..<top_n {
    overrep.push(
      fqqc_obj([
        ("seq", JsonStr(pairs[k].0)),
        ("count", JsonNum(pairs[k].1.to_double())),
      ]),
    )
  }
  fqqc_obj([
    ("module", JsonStr("fqqc")),
    (
      "basic",
      fqqc_obj([
        ("reads", JsonNum(st.n_reads.to_double())),
        ("bases", JsonNum(st.total_bases.to_double())),
        ("mean_length", JsonNum(mean_len)),
        ("mean_gc", JsonNum(overall_gc)),
      ]),
    ),
    (
      "per_base_quality",
      fqqc_obj([
        ("pos", fqqc_arr_int(pos_labels)),
        ("low", fqqc_arr_num(lows)),
        ("q1", fqqc_arr_num(q1s)),
        ("median", fqqc_arr_num(meds)),
        ("q3", fqqc_arr_num(q3s)),
        ("high", fqqc_arr_num(highs)),
        ("mean", fqqc_arr_num(means)),
      ]),
    ),
    ("per_sequence_quality", fqqc_arr_int(qual_hist)),
    ("per_sequence_gc", fqqc_arr_int(gc_hist)),
    ("per_base_gc", fqqc_arr_num(gc_by_pos)),
    ("length_distribution", JsonArr(len_counts)),
    ("adapters", JsonArr(adapter_arr)),
    ("overrepresented", JsonArr(overrep)),
  ])
}

///|
/// 输出 FastQC 风格 JSON 报告文本.
pub fn fqqc_report(fastq : String) -> String raise SeqIOError {
  let st = fqqc_compute(fastq)
  json_to_string(fqqc_report_value(st))
}
