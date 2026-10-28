///|
/// 确定性并行调度演示: 静态轮转 vs LPT 贪心 -> DAG 流水线 -> 失败传播
/// -> 可复现 RNG 子流.
fn main raise {
  // 1. 不均负载: 静态轮转 vs LPT 贪心
  let costs = [10.0, 1.0, 9.0, 2.0, 8.0, 3.0, 7.0, 4.0, 6.0, 5.0]
  let tasks : Array[@src.BpTask] = []
  for i in 0..<costs.length() {
    tasks.push(@src.BpTask::new("task" + i.to_string(), costs[i], [1.0]))
  }
  let static_rep = @src.bp_schedule(
    tasks,
    @src.BpSchedParam::new(3, @src.BpStatic),
  )
  let greedy_rep = @src.bp_schedule(
    tasks,
    @src.BpSchedParam::new(3, @src.BpGreedy),
  )
  println("=== 负载均衡: 10 任务(成本不均) x 3 workers ===")
  println(
    "  静态轮转: makespan " +
    static_rep.makespan.to_string() +
    ", 不均衡度 " +
    static_rep.imbalance.to_string(),
  )
  println(
    "  LPT 贪心: makespan " +
    greedy_rep.makespan.to_string() +
    ", 不均衡度 " +
    greedy_rep.imbalance.to_string(),
  )
  println("  LPT 时间线(任务@worker start~finish):")
  for r in greedy_rep.results {
    println(
      "    " +
      r.id +
      "@w" +
      r.worker.to_string() +
      "  " +
      r.start.to_string() +
      " ~ " +
      r.finish.to_string(),
    )
  }
  // 2. DAG 流水线: 变异 calling(qc -> map -> {snps, indels} -> annotate)
  let pipeline : Array[@src.BpTask] = [
    @src.BpTask::new("qc", 4.0, [1.0]),
    @src.BpTask::with_deps("map", 6.0, [1.0], ["qc"]),
    @src.BpTask::with_deps("call_snps", 5.0, [1.0], ["map"]),
    @src.BpTask::with_deps("call_indels", 3.0, [1.0], ["map"]),
    @src.BpTask::with_deps("annotate", 2.0, [1.0], ["call_snps", "call_indels"]),
  ]
  let dag_rep = @src.bp_dag_schedule(pipeline, 2)
  println("=== DAG 流水线: 5 阶段变异 calling x 2 workers ===")
  for r in dag_rep.results {
    println(
      "    " +
      r.id +
      "@w" +
      r.worker.to_string() +
      "  " +
      r.start.to_string() +
      " ~ " +
      r.finish.to_string() +
      "  " +
      (if r.ok { "ok" } else { "FAIL" }),
    )
  }
  println(
    "  makespan " +
    dag_rep.makespan.to_string() +
    " (关键路径 qc+map+call_snps+annotate)",
  )
  // 3. 失败传播: qc 失败 -> 下游全跳过
  let broken : Array[@src.BpTask] = [
    @src.BpTask::failing("qc", 4.0, []),
    @src.BpTask::with_deps("map", 6.0, [], ["qc"]),
    @src.BpTask::new("index", 2.0, [1.0]),
  ]
  let fail_rep = @src.bp_dag_schedule(broken, 2)
  println("=== 失败传播: qc 失败 ===")
  println(
    "  失败 " +
    fail_rep.n_failed.to_string() +
    ", 跳过 " +
    fail_rep.n_skipped.to_string() +
    " (map 因上游失败未执行; index 独立不受影响)",
  )
  // 4. 可复现 RNG 子流
  let s1 = @src.bp_rng_streams(42, 3)
  let s2 = @src.bp_rng_streams(42, 3)
  let v1 : Array[Int32] = []
  for r in s1 {
    v1.push(r.next_int(1000000))
  }
  let v2 : Array[Int32] = []
  for r in s2 {
    v2.push(r.next_int(1000000))
  }
  println("=== 并行可复现 RNG 子流 (seed=42, 3 流) ===")
  println("  流首值: " + v1.to_string())
  println(
    "  重建一致: " +
    (if v1 == v2 { "true" } else { "false" }) +
    " (worker 各持独立子流, 并行模拟恒定可复现)",
  )
  println("")
  println(
    "注: 本模块为确定性并行调度框架(单线程模拟执行), 调度决策与",
  )
  println(
    "    时间线与真实并行调度器一致; wasm 目标无原生线程.",
  )
}
