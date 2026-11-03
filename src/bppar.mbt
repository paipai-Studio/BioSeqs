///|
/// 确定性并行调度框架(BiocParallel 升级语义).
///
/// MoonBit/wasm 目标无原生线程(核心库无并发原语), 本模块提供并行计算的
/// **语义层**: LPT 贪心负载均衡调度、任务 DAG 依赖调度(关键路径时间线)、
/// fail-fast/carry-on 错误传播、以及并行可复现 RNG 子流(与统一 RNG 协同).
/// 执行为确定性单线程模拟, 但调度决策(分配/时间线/均衡)与真实并行调度器
/// 一致, 输出恒定可复现, 可用于验证并行策略与流水线拓扑.
pub suberror BpSchedError {
  BpSchedError(String)
}

///|
/// 调度策略: BpStatic 轮转分片(基线), BpGreedy LPT 贪心
/// (最重任务优先, 每次分给当前最闲 worker).
pub(all) enum BpStrategy {
  BpStatic
  BpGreedy
}

///|
/// 任务: id + 模拟执行成本(相对单位, 驱动时间线) + 真实 payload +
/// 依赖任务 id 列表(空 = 无依赖) + 失败标记.
pub(all) struct BpTask {
  id : String
  cost : Double
  payload : Array[Double]
  deps : Array[String]
  fails : Bool
}

///|
/// 新建任务(无依赖, 不失败).
pub fn BpTask::new(
  id : String,
  cost : Double,
  payload : Array[Double],
) -> BpTask {
  { id, cost, payload, deps: [], fails: false, }
}

///|
/// 新建带依赖的任务.
pub fn BpTask::with_deps(
  id : String,
  cost : Double,
  payload : Array[Double],
  deps : Array[String],
) -> BpTask {
  { id, cost, payload, deps, fails: false, }
}

///|
/// 新建必定失败的任务(错误传播语义验证).
pub fn BpTask::failing(
  id : String,
  cost : Double,
  payload : Array[Double],
) -> BpTask {
  { id, cost, payload, deps: [], fails: true, }
}

///|
/// 调度器配置.
pub(all) struct BpSchedParam {
  workers : Int32
  strategy : BpStrategy
  stop_on_error : Bool
}

///|
/// 新建调度配置(stop_on_error 为 true 时 fail-fast: 失败后剩余任务跳过).
pub fn BpSchedParam::new(
  workers : Int32,
  strategy : BpStrategy,
  stop_on_error~ : Bool = false,
) -> BpSchedParam {
  { workers, strategy, stop_on_error, }
}

///|
/// 单任务执行结果(value = payload 之和; skipped 任务未执行).
pub(all) struct BpTaskResult {
  id : String
  ok : Bool
  value : Double
  worker : Int32
  start : Double
  finish : Double
  error : String
  skipped : Bool
}

///|
/// 调度报告: 逐任务结果(输入序) + 模拟时间线指标.
/// makespan = 总完工时间; imbalance = (max-min)/max 负载不均衡度;
/// 失败任务按其 cost 计入 worker 负载(执行后才失败).
pub(all) struct BpSchedReport {
  results : Array[BpTaskResult]
  makespan : Double
  worker_loads : Array[Double]
  imbalance : Double
  n_failed : Int32
  n_skipped : Int32
}

///|
/// 无依赖任务集调度(静态轮转 vs LPT 贪心), 输出确定性时间线.
/// workers <= 0 抛 BpSchedError.
pub fn bp_schedule(
  tasks : Array[BpTask],
  param : BpSchedParam,
) -> BpSchedReport raise BpSchedError {
  if param.workers <= 0 {
    raise BpSchedError("workers must be positive")
  }
  let loads : Array[Double] = Array::make(param.workers, 0.0)
  // 调度顺序(产出: 原始下标序)
  let order : Array[Int32] = Array(capacity=tasks.length())
  for i in 0..<tasks.length() {
    order.push(i)
  }
  if param.strategy is BpGreedy {
    bppar_sort_lpt(tasks, order)
  }
  // 时间线: 按 order 依次分配
  let sched : Array[(Int32, Int32, Double)] = Array(capacity=tasks.length())
  let mut seq = 0
  for idx in order {
    let worker : Int32 = if param.strategy is BpStatic {
      seq % param.workers
    } else {
      bppar_min_worker(loads)
    }
    let start = loads[worker]
    loads[worker] = start + tasks[idx].cost
    sched.push((idx, worker, start))
    seq = seq + 1
  }
  // 执行(按调度序): payload 求和; fail-fast 语义
  let n = tasks.length()
  let slots : Array[BpTaskResult?] = Array::make(n, None)
  let mut halted = false
  let mut n_failed = 0
  let mut n_skipped = 0
  for entry in sched {
    let idx = entry.0
    let worker = entry.1
    let start = entry.2
    if halted {
      slots[idx] = Some({
        id: tasks[idx].id,
        ok: false,
        value: 0.0,
        worker: -1,
        start: 0.0,
        finish: 0.0,
        error: "skipped after earlier failure",
        skipped: true,
      })
      n_skipped = n_skipped + 1
      continue
    }
    let finish = start + tasks[idx].cost
    if tasks[idx].fails {
      if param.stop_on_error {
        halted = true
      }
      slots[idx] = Some({
        id: tasks[idx].id,
        ok: false,
        value: 0.0,
        worker,
        start,
        finish,
        error: "task " + tasks[idx].id + " failed",
        skipped: false,
      })
      n_failed = n_failed + 1
    } else {
      let mut value = 0.0
      for v in tasks[idx].payload {
        value = value + v
      }
      slots[idx] = Some({
        id: tasks[idx].id,
        ok: true,
        value,
        worker,
        start,
        finish,
        error: "",
        skipped: false,
      })
    }
  }
  bppar_report(tasks, slots, loads, n_failed, n_skipped)
}

///|
/// DAG 任务调度(依赖拓扑 + 就绪集 LPT: 每步从就绪任务中取 cost 最大者,
/// 分给最早空闲 worker; 任务开始时刻不早于全部依赖的完成时刻).
/// 重复任务 id / 未知依赖 / 循环依赖抛 BpSchedError;
/// 失败任务的下游(直接或间接)标记 skipped.
pub fn bp_dag_schedule(
  tasks : Array[BpTask],
  workers : Int32,
) -> BpSchedReport raise BpSchedError {
  if workers <= 0 {
    raise BpSchedError("workers must be positive")
  }
  let n = tasks.length()
  // id → 下标; 重复 id 检查
  let id_idx : Map[String, Int32] = Map([], capacity=n)
  for i in 0..<n {
    if id_idx.contains(tasks[i].id) {
      raise BpSchedError("duplicate task id: " + tasks[i].id)
    }
    id_idx[tasks[i].id] = i
  }
  // 依赖存在性检查
  for i in 0..<n {
    for d in tasks[i].deps {
      if !id_idx.contains(d) {
        raise BpSchedError("unknown dependency: " + d)
      }
    }
  }
  // 状态: 0 pending / 1 ok / 2 failed / 3 skipped
  let state : Array[Int32] = Array::make(n, 0)
  let start_at : Array[Double] = Array::make(n, 0.0)
  let finish_at : Array[Double] = Array::make(n, 0.0)
  let worker_at : Array[Int32] = Array::make(n, -1)
  let value_at : Array[Double] = Array::make(n, 0.0)
  let loads : Array[Double] = Array::make(workers, 0.0)
  let mut n_failed = 0
  let mut n_skipped = 0
  let mut done = 0
  while done < n {
    // 就绪: pending 且全部依赖已完成(state > 0)
    let mut pick = -1
    for i in 0..<n {
      if state[i] != 0 {
        continue
      }
      let mut ready = true
      for d in tasks[i].deps {
        let di = id_idx[d]
        if state[di] == 0 {
          ready = false
          break
        }
      }
      if ready {
        if pick < 0 || tasks[i].cost > tasks[pick].cost {
          pick = i
        }
      }
    }
    if pick < 0 {
      raise BpSchedError("cyclic dependency in task graph")
    }
    // 上游失败/跳过 → 本任务 skipped(未执行: 不占 worker 负载与时间线)
    let mut upstream_bad = false
    for d in tasks[pick].deps {
      let di = id_idx[d]
      if state[di] == 2 || state[di] == 3 {
        upstream_bad = true
        break
      }
    }
    if upstream_bad {
      state[pick] = 3
      n_skipped = n_skipped + 1
      done = done + 1
    } else {
      // 开始时刻: worker 空闲时刻 与 依赖完成时刻 的最大值
      let worker = bppar_min_worker(loads)
      let mut start = loads[worker]
      for d in tasks[pick].deps {
        let di = id_idx[d]
        if finish_at[di] > start {
          start = finish_at[di]
        }
      }
      let finish = start + tasks[pick].cost
      loads[worker] = finish
      worker_at[pick] = worker
      start_at[pick] = start
      finish_at[pick] = finish
      if tasks[pick].fails {
        state[pick] = 2
        n_failed = n_failed + 1
      } else {
        let mut value = 0.0
        for v in tasks[pick].payload {
          value = value + v
        }
        value_at[pick] = value
        state[pick] = 1
      }
      done = done + 1
    }
  }
  // 汇总(输入序)
  let slots : Array[BpTaskResult?] = Array::make(n, None)
  for i in 0..<n {
    slots[i] = Some({
      id: tasks[i].id,
      ok: state[i] == 1,
      value: value_at[i],
      worker: worker_at[i],
      start: start_at[i],
      finish: finish_at[i],
      error: if state[i] == 2 {
        "task " + tasks[i].id + " failed"
      } else if state[i] == 3 {
        "skipped: upstream failure"
      } else {
        ""
      },
      skipped: state[i] == 3,
    })
  }
  bppar_report(tasks, slots, loads, n_failed, n_skipped)
}

///|
/// 并行可复现 RNG 子流: n 个独立生成器, 种子按 SplitMix64 黄金比例
/// 乘子跳跃, 保证流间低相关; 同一 (seed, n) 恒定重建(并行模拟可复现,
/// 对应 BiocParallel 中 worker 各持独立 RNG 流的最佳实践).
pub fn bp_rng_streams(seed : Int32, n : Int32) -> Array[Rng] raise BpSchedError {
  if n <= 0 {
    raise BpSchedError("rng streams must be positive")
  }
  let base = seed.to_uint64()
  let out : Array[Rng] = Array(capacity=n)
  for i in 0..<n {
    let s = base + (i.to_uint64() + 1UL) * 0x9E3779B97F4A7C15UL
    out.push(Rng::new(s))
  }
  out
}

///|
/// LPT 稳定排序(cost 降序, 同 cost 保持原序).
fn bppar_sort_lpt(tasks : Array[BpTask], order : Array[Int32]) -> Unit {
  let mut i = 0
  while i < order.length() {
    let mut best = i
    let mut j = i + 1
    while j < order.length() {
      let cj = tasks[order[j]].cost
      let cb = tasks[order[best]].cost
      if cj > cb || (cj == cb && order[j] < order[best]) {
        best = j
      }
      j = j + 1
    }
    if best != i {
      let t = order[i]
      order[i] = order[best]
      order[best] = t
    }
    i = i + 1
  }
}

///|
/// 当前累计负载最小的 worker(并列取最小编号).
fn bppar_min_worker(loads : Array[Double]) -> Int32 {
  let mut best = 0
  for i in 1..<loads.length() {
    if loads[i] < loads[best] {
      best = i
    }
  }
  best
}

///|
/// 汇总报告(输入序结果 + makespan/imbalance).
fn bppar_report(
  tasks : Array[BpTask],
  results : Array[BpTaskResult?],
  loads : Array[Double],
  n_failed : Int32,
  n_skipped : Int32,
) -> BpSchedReport {
  let out : Array[BpTaskResult] = Array(capacity=tasks.length())
  for i in 0..<tasks.length() {
    out.push(results[i].unwrap())
  }
  let mut makespan = 0.0
  for r in out {
    if !r.skipped && r.finish > makespan {
      makespan = r.finish
    }
  }
  let mut lo = 0.0
  let mut hi = 0.0
  for l in loads {
    if l > hi {
      hi = l
    }
  }
  let mut first = true
  for l in loads {
    if first || l < lo {
      lo = l
    }
    first = false
  }
  let imbalance = if hi > 0.0 { (hi - lo) / hi } else { 0.0 }
  {
    results: out,
    makespan,
    worker_loads: loads,
    imbalance,
    n_failed,
    n_skipped,
  }
}
