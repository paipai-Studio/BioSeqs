///|
test "bppar task construction and defaults" {
  let t = @src.BpTask::new("t1", 2.0, [1.0, 2.0])
  assert_eq(t.id, "t1")
  assert_eq(t.cost, 2.0)
  assert_eq(t.payload.length(), 2)
  assert_true(t.deps.is_empty())
  assert_true(!t.fails)
  let d = @src.BpTask::with_deps("t2", 1.0, [3.0], ["t1"])
  assert_eq(d.deps, ["t1"])
  let f = @src.BpTask::failing("t3", 1.0, [])
  assert_true(f.fails)
}

///|
test "bppar schedule workers must be positive" {
  let tasks : Array[@src.BpTask] = [@src.BpTask::new("a", 1.0, [1.0])]
  let msg = try {
    ignore(@src.bp_schedule(tasks, @src.BpSchedParam::new(0, @src.BpStatic)))
    "NO_ERROR"
  } catch {
    BpSchedError(m) => m
  }
  assert_true(msg.contains("workers must be positive"))
  let msg2 = try {
    ignore(@src.bp_dag_schedule(tasks, -1))
    "NO_ERROR"
  } catch {
    BpSchedError(m) => m
  }
  assert_true(msg2.contains("workers must be positive"))
  let msg3 = try {
    ignore(@src.bp_rng_streams(1, 0))
    "NO_ERROR"
  } catch {
    BpSchedError(m) => m
  }
  assert_true(msg3.contains("streams must be positive"))
}

///|
test "bppar static round robin assignment" {
  let tasks : Array[@src.BpTask] = []
  for i in 0..<5 {
    tasks.push(@src.BpTask::new("t" + i.to_string(), 1.0, [i.to_double()]))
  }
  let rep = @src.bp_schedule(tasks, @src.BpSchedParam::new(2, @src.BpStatic))
  assert_eq(rep.results[0].worker, 0)
  assert_eq(rep.results[1].worker, 1)
  assert_eq(rep.results[2].worker, 0)
  assert_eq(rep.results[3].worker, 1)
  assert_eq(rep.results[4].worker, 0)
  assert_eq(rep.makespan, 3.0)
  assert_eq(rep.worker_loads, [3.0, 2.0])
}

///|
test "bppar greedy perfect balance" {
  let tasks : Array[@src.BpTask] = []
  for i in 0..<3 {
    tasks.push(@src.BpTask::new("t" + i.to_string(), 10.0, [1.0]))
  }
  let rep = @src.bp_schedule(tasks, @src.BpSchedParam::new(3, @src.BpGreedy))
  assert_eq(rep.makespan, 10.0)
  assert_eq(rep.worker_loads, [10.0, 10.0, 10.0])
  assert_eq(rep.imbalance, 0.0)
  assert_eq(rep.n_failed, 0)
  assert_eq(rep.n_skipped, 0)
}

///|
test "bppar greedy beats static on skewed costs" {
  let costs = [10.0, 1.0, 9.0, 2.0, 8.0, 3.0, 7.0, 4.0, 6.0, 5.0]
  let tasks : Array[@src.BpTask] = []
  for i in 0..<costs.length() {
    tasks.push(@src.BpTask::new("t" + i.to_string(), costs[i], [1.0]))
  }
  let static_rep = @src.bp_schedule(
    tasks,
    @src.BpSchedParam::new(3, @src.BpStatic),
  )
  let greedy_rep = @src.bp_schedule(
    tasks,
    @src.BpSchedParam::new(3, @src.BpGreedy),
  )
  // LPT 贪心: makespan 更小且不均衡度更低
  assert_true(greedy_rep.makespan < static_rep.makespan)
  assert_true(greedy_rep.imbalance <= static_rep.imbalance)
  // 负载总量守恒
  let total = 10.0 + 1.0 + 9.0 + 2.0 + 8.0 + 3.0 + 7.0 + 4.0 + 6.0 + 5.0
  let mut s_sum = 0.0
  for l in static_rep.worker_loads {
    s_sum = s_sum + l
  }
  let mut g_sum = 0.0
  for l in greedy_rep.worker_loads {
    g_sum = g_sum + l
  }
  assert_eq(s_sum, total)
  assert_eq(g_sum, total)
}

///|
test "bppar single worker serial timeline" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::new("a", 3.0, [1.0]),
    @src.BpTask::new("b", 5.0, [2.0]),
  ]
  let rep = @src.bp_schedule(tasks, @src.BpSchedParam::new(1, @src.BpStatic))
  assert_eq(rep.results[0].start, 0.0)
  assert_eq(rep.results[0].finish, 3.0)
  assert_eq(rep.results[1].start, 3.0)
  assert_eq(rep.results[1].finish, 8.0)
  assert_eq(rep.makespan, 8.0)
}

///|
test "bppar payload value computed" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::new("a", 1.0, [1.5, 2.5, 3.0]),
    @src.BpTask::new("b", 1.0, []),
  ]
  let rep = @src.bp_schedule(tasks, @src.BpSchedParam::new(2, @src.BpStatic))
  assert_eq(rep.results[0].value, 7.0)
  assert_true(rep.results[0].ok)
  assert_eq(rep.results[1].value, 0.0)
}

///|
test "bppar failure carry on" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::new("a", 1.0, [1.0]),
    @src.BpTask::failing("b", 1.0, [1.0]),
    @src.BpTask::new("c", 1.0, [1.0]),
  ]
  let rep = @src.bp_schedule(tasks, @src.BpSchedParam::new(1, @src.BpStatic))
  assert_eq(rep.n_failed, 1)
  assert_eq(rep.n_skipped, 0)
  assert_true(rep.results[1].error.contains("failed"))
  assert_true(rep.results[2].ok)
}

///|
test "bppar fail fast skips rest" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::new("a", 1.0, [1.0]),
    @src.BpTask::failing("b", 1.0, [1.0]),
    @src.BpTask::new("c", 1.0, [1.0]),
  ]
  let rep = @src.bp_schedule(
    tasks,
    @src.BpSchedParam::new(1, @src.BpStatic, stop_on_error=true),
  )
  assert_eq(rep.n_failed, 1)
  assert_eq(rep.n_skipped, 1)
  assert_true(rep.results[2].skipped)
  assert_eq(rep.results[2].worker, -1)
  assert_true(rep.results[2].error.contains("skipped"))
  // 失败任务自身已执行, 计入负载; skipped 不计入 makespan
  assert_eq(rep.makespan, 2.0)
}

///|
test "bppar dag no deps degenerates" {
  let tasks : Array[@src.BpTask] = []
  for i in 0..<3 {
    tasks.push(@src.BpTask::new("t" + i.to_string(), 2.0, [1.0]))
  }
  let rep = @src.bp_dag_schedule(tasks, 3)
  assert_eq(rep.makespan, 2.0)
  assert_eq(rep.imbalance, 0.0)
  assert_true(rep.results[0].ok)
}

///|
test "bppar dag dependency ordering" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::new("qc", 4.0, [1.0]),
    @src.BpTask::with_deps("map", 3.0, [1.0], ["qc"]),
    @src.BpTask::new("index", 2.0, [1.0]),
  ]
  let rep = @src.bp_dag_schedule(tasks, 2)
  let qc = rep.results[0]
  let m = rep.results[1]
  let idx = rep.results[2]
  // 依赖约束: map 不早于 qc 完成
  assert_true(m.start >= qc.finish)
  // 确定性时间线: qc(4) 先调度, map(3) 次之(就绪集 LPT), index(2) 最后
  assert_eq(qc.start, 0.0)
  assert_eq(qc.finish, 4.0)
  assert_eq(m.start, 4.0)
  assert_eq(m.finish, 7.0)
  assert_eq(idx.start, 4.0)
  assert_eq(idx.finish, 6.0)
  // makespan = 关键路径 qc -> map
  assert_eq(rep.makespan, 7.0)
}

///|
test "bppar dag cyclic raises" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::with_deps("a", 1.0, [], ["b"]),
    @src.BpTask::with_deps("b", 1.0, [], ["a"]),
  ]
  let msg = try {
    ignore(@src.bp_dag_schedule(tasks, 2))
    "NO_ERROR"
  } catch {
    BpSchedError(m) => m
  }
  assert_true(msg.contains("cyclic dependency"))
}

///|
test "bppar dag unknown dependency raises" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::with_deps("a", 1.0, [], ["nope"]),
  ]
  let msg = try {
    ignore(@src.bp_dag_schedule(tasks, 1))
    "NO_ERROR"
  } catch {
    BpSchedError(m) => m
  }
  assert_true(msg.contains("unknown dependency"))
}

///|
test "bppar dag duplicate id raises" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::new("a", 1.0, []),
    @src.BpTask::new("a", 1.0, []),
  ]
  let msg = try {
    ignore(@src.bp_dag_schedule(tasks, 1))
    "NO_ERROR"
  } catch {
    BpSchedError(m) => m
  }
  assert_true(msg.contains("duplicate task id"))
}

///|
test "bppar dag failure propagation" {
  let tasks : Array[@src.BpTask] = [
    @src.BpTask::failing("qc", 1.0, []),
    @src.BpTask::with_deps("map", 1.0, [], ["qc"]),
    @src.BpTask::with_deps("call", 1.0, [], ["map"]),
    @src.BpTask::new("index", 1.0, []),
  ]
  let rep = @src.bp_dag_schedule(tasks, 2)
  assert_eq(rep.n_failed, 1)
  assert_eq(rep.n_skipped, 2)
  assert_true(rep.results[0].error.contains("failed"))
  assert_true(rep.results[1].skipped)
  assert_true(rep.results[2].skipped)
  assert_true(rep.results[3].ok)
  // skipped 任务未执行
  assert_eq(rep.results[1].worker, -1)
}

///|
test "bppar rng streams reproducible and independent" {
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
  assert_eq(v1, v2)
  // 流间独立: 不同流首值应不同
  assert_true(v1[0] != v1[1])
  assert_true(v1[1] != v1[2])
}
