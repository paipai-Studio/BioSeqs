// ===========================================================================
// msprime_lite_test.mbt — msprime 风格群体谱系模拟器测试.
// ===========================================================================

///|
fn ms_test_rng() -> @src.Rng {
  @src.Rng::new(20261031UL)
}

///|
test "msprime: 常数 Ne 合并时间分布" {
  let params = @src.msprime_params()
  let result = @src.msprime_simulate(params, ms_test_rng())
  // 10 样本, 9 次合并
  inspect(result.n_mutations, content="0")
  // 无突变率时 SFS 全 0
  for v in result.joint_sfs {
    assert_eq(v, 0)
  }
}

///|
test "msprime: 时变 Ne 扩张模型" {
  // 最近 1000 代 Ne=1000, 之前 Ne=10000(扩张)
  let epochs = [
    @src.MsprimeEpoch::{ start_time: 0.0, end_time: 1000.0, ne: 1000.0, },
    @src.MsprimeEpoch::{ start_time: 1000.0, end_time: 1.0e30, ne: 10000.0, },
  ]
  let params : @src.MsprimeParams = {
    n_samples: 20,
    n_pops: 1,
    epochs,
    migration: None,
    recombination_rate: 0.0,
    sequence_length: 1000,
    mutation_rate: 1.0e-8,
    n_blocks: 1,
  }
  let result = @src.msprime_simulate(params, ms_test_rng())
  // 有突变率时产生突变
  assert_true(result.n_mutations > 0)
  // SFS 长度 = n-1 = 19
  inspect(result.joint_sfs.length(), content="19")
  // 折叠 SFS 长度 = (n/2)+1 = 11
  inspect(result.folded_sfs.length(), content="11")
  // Watterson theta > 0
  assert_true(result.watterson_theta > 0.0)
}

///|
test "msprime: 重组多块模拟" {
  let epochs = [
    @src.MsprimeEpoch::{ start_time: 0.0, end_time: 1.0e30, ne: 1000.0, },
  ]
  let params : @src.MsprimeParams = {
    n_samples: 10,
    n_pops: 1,
    epochs,
    migration: None,
    recombination_rate: 1.0e-8,
    sequence_length: 10000,
    mutation_rate: 1.0e-8,
    n_blocks: 4,
  }
  let result = @src.msprime_simulate(params, ms_test_rng())
  // 4 个块
  inspect(result.blocks.length(), content="4")
  // 每块起始位置
  assert_eq(result.blocks[0].start, 0)
  assert_eq(result.blocks[1].start, 2500)
  assert_eq(result.blocks[3].start, 7500)
}

///|
test "msprime: demes 风格解析" {
  let demes = @src.msprime_parse_demes(
    "# demes model\n" +
    "- start_time: 0\n  end_time: 500\n  ne: 1000\n" +
    "- start_time: 500\n  end_time: 10000\n  ne: 10000\n",
  )
  inspect(demes.length(), content="2")
  inspect(demes[0].ne, content="1000")
  inspect(demes[1].ne, content="10000")
  inspect(demes[1].start_time, content="500")
}

///|
test "msprime: Tajima's D 与 Watterson theta" {
  // 中等样本数 + 高突变率 -> 可观测统计量
  let epochs = [
    @src.MsprimeEpoch::{ start_time: 0.0, end_time: 1.0e30, ne: 1000.0, },
  ]
  let params : @src.MsprimeParams = {
    n_samples: 30,
    n_pops: 1,
    epochs,
    migration: None,
    recombination_rate: 0.0,
    sequence_length: 5000,
    mutation_rate: 1.0e-7,
    n_blocks: 1,
  }
  let result = @src.msprime_simulate(params, ms_test_rng())
  assert_true(result.n_mutations > 0)
  // Tajima's D 与 Watterson theta 有限(非 NaN)
  assert_eq(result.tajima_d == result.tajima_d, true)
  assert_eq(result.watterson_theta == result.watterson_theta, true)
}

///|
test "msprime: JSON 与 tree sequence 导出" {
  let params = @src.msprime_params()
  let result = @src.msprime_simulate(params, ms_test_rng())
  let json = @src.msprime_json(result)
  assert_eq(json.contains("\"n_mutations\""), true)
  assert_eq(json.contains("\"tajima_d\""), true)
  assert_eq(json.contains("\"watterson_theta\""), true)
  let ts = @src.msprime_tree_sequence(result)
  assert_eq(ts.contains("# position"), true)
  assert_eq(ts.contains("\tcoalescence_times"), true)
  let line = @src.msprime_summary_line(result)
  assert_eq(line.contains("msprime-lite:"), true)
}

///|
test "msprime: 空参数" {
  let params : @src.MsprimeParams = {
    n_samples: 0,
    n_pops: 1,
    epochs: [],
    migration: None,
    recombination_rate: 0.0,
    sequence_length: 1000,
    mutation_rate: 0.0,
    n_blocks: 1,
  }
  let result = @src.msprime_simulate(params, ms_test_rng())
  inspect(result.joint_sfs.length(), content="0")
  inspect(result.n_mutations, content="0")
}
