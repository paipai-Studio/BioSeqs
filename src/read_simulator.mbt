///|
/// read_simulator.mbt — 参考序列 -> 测序读段(FASTQ)本地模拟器.
///
/// 对齐 wgsim / ART 的核心语义(本地确定性实现, 数据不出本机):
/// - 单端 / 双端测序
/// - 目标覆盖度、读长、插入片段长度正态分布
/// - 逐碱基 Phred 质量模型(3' 端线性衰减 + 噪声)
/// - 碱基错误、SNP/indel 变异注入
/// - adapter 污染(读段 3' 端追加)
/// - 同种子结果完全可复现, 便于与比对/变异检出模块做闭环测试

///|
/// 读段模拟错误.
pub suberror ReadSimulatorError {
  ReadSimulatorError(String)
}

///|
/// IUPAC 互补表.
let read_sim_complement : Map[String, String] = Map(
  [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")],
  capacity=5,
)

///|
/// 单端/双端模拟配置.
pub(all) struct ReadSimulatorConfig {
  read_length : Int32
  coverage : Double
  paired_end : Bool
  insert_mean : Int32
  insert_sd : Int32
  error_rate : Double
  quality_max : Int32
  quality_min : Int32
  quality_decay : Double
  snp_rate : Double
  indel_rate : Double
  adapter : String
  adapter_rate : Double
  seed : UInt64
}

///|
/// 构造配置并校验参数.
pub fn read_simulator_config(
  read_length~ : Int32 = 150,
  coverage~ : Double = 30.0,
  paired_end~ : Bool = true,
  insert_mean~ : Int32 = 350,
  insert_sd~ : Int32 = 50,
  error_rate~ : Double = 0.001,
  quality_max~ : Int32 = 40,
  quality_min~ : Int32 = 2,
  quality_decay~ : Double = 10.0,
  snp_rate~ : Double = 0.0,
  indel_rate~ : Double = 0.0,
  adapter~ : String = "",
  adapter_rate~ : Double = 0.0,
  seed~ : UInt64 = 42UL,
) -> ReadSimulatorConfig raise ReadSimulatorError {
  if read_length <= 0 {
    raise ReadSimulatorError("read length must be positive")
  }
  if coverage < 0.0 {
    raise ReadSimulatorError("coverage must be non-negative")
  }
  if error_rate < 0.0 || error_rate > 1.0 {
    raise ReadSimulatorError("error rate must be in [0, 1]")
  }
  if snp_rate < 0.0 || snp_rate > 1.0 || indel_rate < 0.0 || indel_rate > 1.0 {
    raise ReadSimulatorError("variant rates must be in [0, 1]")
  }
  if adapter_rate < 0.0 || adapter_rate > 1.0 {
    raise ReadSimulatorError("adapter rate must be in [0, 1]")
  }
  if quality_min < 0 || quality_max < quality_min {
    raise ReadSimulatorError("invalid quality range")
  }
  {
    read_length,
    coverage,
    paired_end,
    insert_mean,
    insert_sd,
    error_rate,
    quality_max,
    quality_min,
    quality_decay,
    snp_rate,
    indel_rate,
    adapter,
    adapter_rate,
    seed,
  }
}

///|
/// 随机选一个与 `base` 不同的碱基.
fn read_sim_other_base(rng : Rng, base : Char) -> Char {
  let bases = ['A', 'C', 'G', 'T']
  let idx = match base {
    'A' => 0
    'C' => 1
    'G' => 2
    'T' => 3
    _ => 0
  }
  let alt = (idx + 1 + rng.randint(3)) % 4
  bases[alt]
}

///|
/// 反向互补.
fn read_sim_revcomp(seq : String) -> String {
  let sb = StringBuilder()
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq.unsafe_get(i).unsafe_to_char().to_ascii_uppercase()
    sb.write_string(read_sim_complement.get(c.to_string()).unwrap_or("N"))
  }
  sb.to_string()
}

///|
/// 向参考序列注入 SNP / indel, 返回 (模拟基因组, 事件数).
fn read_sim_mutate(
  reference : String,
  config : ReadSimulatorConfig,
) -> (String, Int32) {
  let rng = Rng::new(config.seed ^ 0x9E3779B97F4A7C15UL)
  let sb = StringBuilder()
  let mut events = 0
  for i = 0; i < reference.length(); i = i + 1 {
    let mut base = reference.unsafe_get(i).unsafe_to_char().to_ascii_uppercase()
    if rng.next_double() < config.snp_rate {
      base = read_sim_other_base(rng, base)
      events = events + 1
    }
    if rng.next_double() < config.indel_rate {
      events = events + 1
      if rng.next_double() < 0.5 {
        // 缺失: 不写出该碱基
        continue
      } else {
        // 插入: 在当前位置前插入一个随机碱基
        sb.write_char(['A', 'C', 'G', 'T'][rng.randint(4)])
      }
    }
    sb.write_char(base)
  }
  (sb.to_string(), events)
}

///|
/// 生成逐碱基 Phred 质量.
fn read_sim_qualities(
  rng : Rng,
  length : Int32,
  config : ReadSimulatorConfig,
) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=length)
  for i = 0; i < length; i = i + 1 {
    let frac = if length > 1 {
      i.to_double() / (length - 1).to_double()
    } else {
      0.0
    }
    let raw = config.quality_max.to_double() -
      config.quality_decay * frac +
      rng.normal(0.0, 1.5)
    let q = if raw < config.quality_min.to_double() {
      config.quality_min
    } else if raw > config.quality_max.to_double() {
      config.quality_max
    } else {
      raw.to_int()
    }
    out.push(q)
  }
  out
}

///|
/// 按错误率污染读段, 返回 (读段, 错误数).
fn read_sim_errors(
  rng : Rng,
  seq : String,
  config : ReadSimulatorConfig,
) -> (String, Int32) {
  if config.error_rate <= 0.0 {
    return (seq, 0)
  }
  let sb = StringBuilder()
  let mut errors = 0
  for i = 0; i < seq.length(); i = i + 1 {
    let base = seq.unsafe_get(i).unsafe_to_char().to_ascii_uppercase()
    if rng.next_double() < config.error_rate {
      sb.write_char(read_sim_other_base(rng, base))
      errors = errors + 1
    } else {
      sb.write_char(base)
    }
  }
  (sb.to_string(), errors)
}

///|
/// 读取质量编码(Phred+33).
fn read_sim_quality_string(qualities : Array[Int32]) -> String {
  let sb = StringBuilder()
  for q in qualities {
    sb.write_char((q + 33).to_uint16().unsafe_to_char())
  }
  sb.to_string()
}

///|
/// 模拟结果.
pub(all) struct ReadSimulatorResult {
  fastq : String
  n_reads : Int32
  n_pairs : Int32
  mean_read_length : Double
  achieved_coverage : Double
  n_variants : Int32
  n_errors : Int32
}

///|
/// 从参考序列模拟测序读段.
pub fn simulate_reads(
  reference : String,
  config~ : ReadSimulatorConfig,
) -> ReadSimulatorResult raise ReadSimulatorError {
  if reference.length() == 0 {
    raise ReadSimulatorError("reference sequence is empty")
  }
  if reference.length() < config.read_length {
    raise ReadSimulatorError(
      "reference length \{reference.length()} shorter than read length \{config.read_length}",
    )
  }
  let (sample, n_variants) = read_sim_mutate(reference, config)
  if sample.length() < config.read_length {
    raise ReadSimulatorError("mutated reference shorter than read length")
  }
  let rng = Rng::new(config.seed)
  let ref_len = reference.length()
  let read_len = config.read_length
  let frag_min = read_len * 2
  let sb = StringBuilder()
  let mut n_reads = 0
  let mut n_pairs = 0
  let mut total_bases = 0
  let mut n_errors = 0

  if config.paired_end {
    let per_pair = (read_len * 2).to_double()
    let wanted = if per_pair > 0.0 {
      (config.coverage * ref_len.to_double() / per_pair).to_int()
    } else {
      0
    }
    for pi = 0; pi < wanted; pi = pi + 1 {
      let frag = if config.insert_sd > 0 {
        let raw = config.insert_mean.to_double() +
          rng.normal(0.0, config.insert_sd.to_double())
        if raw.to_int() < frag_min {
          frag_min
        } else {
          raw.to_int()
        }
      } else if config.insert_mean < frag_min {
        frag_min
      } else {
        config.insert_mean
      }
      let frag_len = if frag > sample.length() { sample.length() } else { frag }
      let max_start = sample.length() - frag_len
      let start = if max_start > 0 { rng.randint(max_start + 1) } else { 0 }
      let fragment = sample[start:start + frag_len].to_owned()
      let r1 = fragment[0:read_len].to_owned()
      let r2 = read_sim_revcomp(
        fragment[frag_len - read_len:frag_len].to_owned(),
      )
      let (q1, e1) = read_sim_errors(rng, r1, config)
      let (q2, e2) = read_sim_errors(rng, r2, config)
      n_errors = n_errors + e1 + e2
      let features1 = read_sim_qualities(rng, read_len, config)
      let features2 = read_sim_qualities(rng, read_len, config)
      sb.write_string(
        "@read\{pi.to_string()}/1\n\{q1}\n+\n\{read_sim_quality_string(features1)}\n",
      )
      sb.write_string(
        "@read\{pi.to_string()}/2\n\{q2}\n+\n\{read_sim_quality_string(features2)}\n",
      )
      total_bases = total_bases + read_len * 2
      n_reads = n_reads + 2
      n_pairs = n_pairs + 1
    }
  } else {
    let wanted = if read_len > 0 {
      (config.coverage * ref_len.to_double() / read_len.to_double()).to_int()
    } else {
      0
    }
    let max_start = sample.length() - read_len
    for ri = 0; ri < wanted; ri = ri + 1 {
      let start = if max_start > 0 { rng.randint(max_start + 1) } else { 0 }
      let raw_seq = sample[start:start + read_len].to_owned()
      let (seq, errs) = read_sim_errors(rng, raw_seq, config)
      n_errors = n_errors + errs
      let features = read_sim_qualities(rng, read_len, config)
      if config.adapter_rate > 0.0 && rng.next_double() < config.adapter_rate {
        let adapter_qual : Array[Int32] = Array::make(
          config.adapter.length(),
          config.quality_min,
        )
        sb.write_string(
          "@read\{ri.to_string()}\n\{seq}\{config.adapter}\n+\n\{read_sim_quality_string(features)}\{read_sim_quality_string(adapter_qual)}\n",
        )
        total_bases = total_bases + read_len + config.adapter.length()
      } else {
        sb.write_string(
          "@read\{ri.to_string()}\n\{seq}\n+\n\{read_sim_quality_string(features)}\n",
        )
        total_bases = total_bases + read_len
      }
      n_reads = n_reads + 1
    }
  }
  {
    fastq: sb.to_string(),
    n_reads,
    n_pairs,
    mean_read_length: if n_reads > 0 {
      total_bases.to_double() / n_reads.to_double()
    } else {
      0.0
    },
    achieved_coverage: total_bases.to_double() / ref_len.to_double(),
    n_variants,
    n_errors,
  }
}
