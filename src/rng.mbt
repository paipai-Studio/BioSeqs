///|
/// 可复现伪随机数生成器(SplitMix64).
///
/// 统一供 bootstrap、置换检验、数据模拟、随机树生成等场景使用, 取代散落
/// 各模块的手写 LCG. 同一 seed 的输出序列恒定, 保证测试与流水线可复现;
/// 生成质量(周期 2^64、均匀性)显著优于 LCG.
///
/// 注意: 既有模块内部的局部 LCG 保持原样以避免统计行为漂移, 新代码应使用本模块.
pub struct Rng {
  mut state : UInt64
}

///|
/// 以 UInt64 种子构造生成器. seed 0 同样有效(SplitMix64 无坏种子).
pub fn Rng::new(seed : UInt64) -> Rng {
  { state: seed, }
}

///|
/// 以 Int 种子构造生成器(常见入参形态, 负数按位模式解释).
pub fn Rng::from_int(seed : Int32) -> Rng {
  { state: seed.to_uint64(), }
}

///|
/// SplitMix64 核心步进, 返回下一个 UInt64.
pub fn Rng::next_u64(self : Rng) -> UInt64 {
  self.state = self.state + 0x9E3779B97F4A7C15UL
  let mut z = self.state
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL
  z ^ (z >> 31)
}

///|
/// 返回下一个 UInt(取高 32 位).
pub fn Rng::next_u32(self : Rng) -> UInt32 {
  (self.next_u64() >> 32).to_uint()
}

///|
/// 均匀采样 [0, 1) 的 Double(53 位精度).
pub fn Rng::next_double(self : Rng) -> Double {
  (self.next_u64() >> 11).to_double() * (1.0 / 9007199254740992.0)
}

///|
/// 均匀采样 [lower, upper) 的 Double.
pub fn Rng::uniform(self : Rng, lower : Double, upper : Double) -> Double {
  lower + self.next_double() * (upper - lower)
}

///|
/// 随机整数 [0, n), n <= 0 时返回 0(兼容既有调用).
pub fn Rng::randint(self : Rng, n : Int32) -> Int32 {
  if n <= 0 {
    return 0
  }
  let r : Int32 = (self.next_u64() % n.to_uint64()).to_int()
  r
}

///|
/// 均匀采样 [0, bound) 的 Int. bound <= 1 时返回 0;
/// 取模偏差在 bound 远小于 2^64 时可忽略.
pub fn Rng::next_int(self : Rng, bound : Int32) -> Int32 {
  if bound <= 1 {
    return 0
  }
  (self.next_u64() % bound.to_uint64()).to_int()
}

///|
/// 正态分布 N(mu, sigma) 抽样(Box-Muller 变换).
pub fn Rng::normal(self : Rng, mu : Double, sigma : Double) -> Double {
  let u1 = self.next_double()
  let u2 = self.next_double()
  // u1 取 (0, 1]: 避开 log(0)
  let a = if u1 <= 0.0 { 1.0 } else { u1 }
  mu + sigma * (-2.0 * @math.ln(a)).sqrt() * @math.cos(2.0 * @math.PI * u2)
}

///|
/// 指数分布(速率 rate)抽样, 均值为 1/rate.
pub fn Rng::exp(self : Rng, rate : Double) -> Double {
  let u = self.next_double()
  // u 取 (0, 1]: 避开 log(0)
  let a = if u <= 0.0 { 1.0 } else { u }
  -@math.ln(a) / rate
}

///|
/// 就地 Fisher-Yates 洗牌.
pub fn[T] Rng::shuffle(self : Rng, arr : Array[T]) -> Unit {
  let n = arr.length()
  let mut i = n - 1
  while i > 0 {
    let j = self.next_int(i + 1)
    let t = arr[i]
    arr[i] = arr[j]
    arr[j] = t
    i = i - 1
  }
}

///|
/// 从 [0, total) 无放回抽取 k 个不重复 Int(部分 Fisher-Yates),
/// 返回顺序随机; k > total 时返回全部, k < 0 时返回空.
pub fn Rng::sample_int(self : Rng, total : Int32, k : Int32) -> Array[Int32] {
  let take = if k > total { total } else if k < 0 { 0 } else { k }
  let pool = Array::makei(total, i => i)
  for i in 0..<take {
    let j = i + self.next_int(total - i)
    let t = pool[i]
    pool[i] = pool[j]
    pool[j] = t
  }
  let out : Array[Int32] = Array(capacity=take)
  for i in 0..<take {
    out.push(pool[i])
  }
  out
}
