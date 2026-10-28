///|
test "rng deterministic same seed" {
  let a = @src.Rng::new(42UL)
  let b = @src.Rng::new(42UL)
  let c = @src.Rng::new(43UL)
  let same = @moonbit.assert_close_e9(a.next_double(), b.next_double())
  assert_true(same)
  let diff = @moonbit.assert_close_e9(a.next_double(), c.next_double())
  assert_false(diff)
  let a2 = @src.Rng::new(42UL)
  let b2 = @src.Rng::new(42UL)
  assert_eq(a2.next_u64(), b2.next_u64())
}

///|
test "rng next_double range" {
  let r = @src.Rng::from_int(7)
  let mut i = 0
  while i < 10000 {
    let v = r.next_double()
    assert_true(v >= 0.0 && v < 1.0)
    i = i + 1
  }
}

///|
test "rng next_int range" {
  let r = @src.Rng::from_int(7)
  let mut i = 0
  while i < 10000 {
    let v = r.next_int(6)
    assert_true(v >= 0 && v < 6)
    i = i + 1
  }
  assert_eq(r.next_int(1), 0)
  assert_eq(r.next_int(0), 0)
  assert_eq(r.next_int(-5), 0)
}

///|
test "rng next_int coverage" {
  let r = @src.Rng::from_int(123)
  let seen = [false, false, false, false, false]
  let mut i = 0
  while i < 200 {
    seen[r.next_int(5)] = true
    i = i + 1
  }
  for s in seen {
    assert_true(s)
  }
}

///|
test "rng uniform range" {
  let r = @src.Rng::from_int(99)
  let mut i = 0
  while i < 1000 {
    let v = r.uniform(2.0, 5.0)
    assert_true(v >= 2.0 && v < 5.0)
    i = i + 1
  }
}

///|
test "rng normal mean std" {
  let r = @src.Rng::from_int(2026)
  let n = 20000
  let mut sum = 0.0
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let v = r.normal(0.0, 1.0)
    sum = sum + v
    sum_sq = sum_sq + v * v
    i = i + 1
  }
  let mean = sum / n.to_double()
  let variance = sum_sq / n.to_double() - mean * mean
  // 样本均值标准误 ~ 1/sqrt(n) ≈ 0.007, 方差相对误差阈值放宽
  assert_true(mean.abs() < 0.05)
  assert_true(variance > 0.9 && variance < 1.1)
}

///|
test "rng exp mean" {
  let r = @src.Rng::from_int(555)
  let n = 20000
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + r.exp(2.0)
    i = i + 1
  }
  let mean = sum / n.to_double()
  // 期望 0.5, 标准误 ~ 0.5/sqrt(n) ≈ 0.0035
  assert_true((mean - 0.5).abs() < 0.03)
}

///|
test "rng shuffle preserves multiset" {
  let r = @src.Rng::from_int(77)
  let arr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  r.shuffle(arr)
  assert_eq(arr.length(), 10)
  let sorted = arr.copy()
  let mut i = 0
  while i < sorted.length() {
    let mut j = i + 1
    while j < sorted.length() {
      if sorted[j] < sorted[i] {
        let t = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  assert_eq(sorted, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
}

///|
test "rng shuffle empty and single" {
  let r = @src.Rng::from_int(1)
  let empty : Array[Int32] = []
  r.shuffle(empty)
  assert_eq(empty.length(), 0)
  let one = [5]
  r.shuffle(one)
  assert_eq(one, [5])
}

///|
test "rng sample_int distinct and bounded" {
  let r = @src.Rng::from_int(8888)
  let picked = r.sample_int(100, 10)
  assert_eq(picked.length(), 10)
  let seen : Map[Int32, Bool] = Map([], capacity=20)
  for v in picked {
    assert_true(v >= 0 && v < 100)
    assert_false(seen.contains(v))
    seen[v] = true
  }
}

///|
test "rng sample_int edge cases" {
  let r = @src.Rng::from_int(3)
  assert_eq(r.sample_int(5, 0).length(), 0)
  assert_eq(r.sample_int(5, -1).length(), 0)
  let all = r.sample_int(5, 10)
  assert_eq(all.length(), 5)
}

///|
test "rng next_u32 range" {
  let r = @src.Rng::new(0UL)
  let mut i = 0
  while i < 1000 {
    let v = r.next_u32()
    assert_true(v >= 0UL.to_uint())
    i = i + 1
  }
}
