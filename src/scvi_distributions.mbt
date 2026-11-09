// scvi_distributions.mbt - Probability distributions for scvi-tools models
// Mirrors scvi.distributions: NegativeBinomial, ZeroInflatedNegativeBinomial,
// Poisson, Bernoulli, NegativeBinomialMixture and MixturePoisson.
//
// All distributions are parameterised element-wise so that they can be used
// directly on the per-gene rate vectors produced by scVI-style decoders.

// ===========================================================================
// Shared numerical helpers
// ===========================================================================

///|
/// Clamp a rate away from zero before taking logs / ratios.
pub fn scvi_dist_clamp(x : Double) -> Double {
  if x < 1.0e-8 {
    1.0e-8
  } else {
    x
  }
}

///|
fn scvi_dist_log_sum_exp(values : Array[Double]) -> Double {
  let mut m = values[0]
  for i in 1..<values.length() {
    if values[i] > m {
      m = values[i]
    }
  }
  let mut s = 0.0
  for i in 0..<values.length() {
    s += @math.exp(values[i] - m)
  }
  m + @math.ln(s)
}

///|
fn scvi_dist_standard_normal(seed : Int32) -> Double {
  let u1 = scvi_dist_uniform(seed)
  let u2 = scvi_dist_uniform(seed + 7919)
  let r = @math.pow(-2.0 * @math.ln(u1 + 1.0e-12), 0.5)
  r * @math.cos(6.283185307179586 * u2)
}

///|
fn scvi_dist_uniform(seed : Int32) -> Double {
  let s = (seed * 1103515245 + 12345).reinterpret_as_uint()
  (s % 100000).to_double() / 100000.0
}

///|
/// Marsaglia-Tsang gamma sampler (shape/scale parameterisation).
fn scvi_dist_gamma_sample(
  shape : Double,
  scale : Double,
  seed : Int32,
) -> Double {
  let safe_shape = if shape < 1.0e-6 { 1.0e-6 } else { shape }
  if safe_shape < 1.0 {
    let u = scvi_dist_uniform(seed)
    let inner = scvi_dist_gamma_sample(safe_shape + 1.0, 1.0, seed + 13)
    return inner * @math.pow(u + 1.0e-12, 1.0 / safe_shape) * scale
  }
  let d = safe_shape - 1.0 / 3.0
  let c = 1.0 / @math.pow(9.0 * d, 0.5)
  let mut result = safe_shape * scale
  for i in 0..<64 {
    let x = scvi_dist_standard_normal(seed + i * 37)
    let v = @math.pow(1.0 + c * x, 3.0)
    if v > 0.0 {
      let u = scvi_dist_uniform(seed + i * 37 + 5)
      if @math.ln(u + 1.0e-12) < 0.5 * x * x + d - d * v + d * @math.ln(v) {
        result = d * v * scale
        break
      }
    }
  }
  result
}

///|
fn scvi_dist_poisson_sample(lam : Double, seed : Int32) -> Double {
  let safe_lam = if lam < 0.0 { 0.0 } else { lam }
  if safe_lam < 30.0 {
    let limit = @math.exp(-safe_lam)
    let mut p = 1.0
    let mut k = 0
    while k < 10000 {
      k = k + 1
      p = p * scvi_dist_uniform(seed + k * 17)
      if p <= limit {
        break
      }
    }
    (k - 1).to_double()
  } else {
    let z = scvi_dist_standard_normal(seed + 101)
    let val = safe_lam + @math.pow(safe_lam, 0.5) * z
    if val < 0.0 {
      0.0
    } else {
      @math.round(val)
    }
  }
}

// ===========================================================================
// Poisson
// ===========================================================================

///|
pub struct ScviPoisson {
  rate : Array[Double]
}

///|
pub fn ScviPoisson::new(rate : Array[Double]) -> ScviPoisson {
  { rate, }
}

///|
pub fn ScviPoisson::log_prob(
  self : ScviPoisson,
  x : Array[Double],
) -> Array[Double] {
  let n = self.rate.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let mu = scvi_dist_clamp(self.rate[i])
    out.push(x[i] * @math.ln(mu) - mu - lgamma(x[i] + 1.0))
  }
  out
}

///|
pub fn ScviPoisson::mean(self : ScviPoisson) -> Array[Double] {
  self.rate
}

///|
pub fn ScviPoisson::variance(self : ScviPoisson) -> Array[Double] {
  self.rate
}

///|
pub fn ScviPoisson::sample(self : ScviPoisson, seed : Int32) -> Array[Double] {
  let n = self.rate.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    out.push(scvi_dist_poisson_sample(self.rate[i], seed + i * 31))
  }
  out
}

// ===========================================================================
// Negative binomial (mean / inverse-dispersion parameterisation)
// ===========================================================================

///|
pub struct ScviNegativeBinomial {
  mu : Array[Double]
  theta : Array[Double]
}

///|
pub fn ScviNegativeBinomial::new(
  mu : Array[Double],
  theta : Array[Double],
) -> ScviNegativeBinomial {
  { mu, theta, }
}

///|
pub fn ScviNegativeBinomial::log_prob(
  self : ScviNegativeBinomial,
  x : Array[Double],
) -> Array[Double] {
  let n = self.mu.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let mu = scvi_dist_clamp(self.mu[i])
    let theta = scvi_dist_clamp(self.theta[i])
    let lg = lgamma(x[i] + theta) - lgamma(theta) - lgamma(x[i] + 1.0)
    let lp = theta * @math.ln(theta / (theta + mu)) +
      x[i] * @math.ln(mu / (theta + mu))
    out.push(lg + lp)
  }
  out
}

///|
pub fn ScviNegativeBinomial::mean(self : ScviNegativeBinomial) -> Array[Double] {
  self.mu
}

///|
pub fn ScviNegativeBinomial::variance(
  self : ScviNegativeBinomial,
) -> Array[Double] {
  let n = self.mu.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let mu = self.mu[i]
    let theta = scvi_dist_clamp(self.theta[i])
    out.push(mu + mu * mu / theta)
  }
  out
}

///|
pub fn ScviNegativeBinomial::sample(
  self : ScviNegativeBinomial,
  seed : Int32,
) -> Array[Double] {
  let n = self.mu.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let mu = scvi_dist_clamp(self.mu[i])
    let theta = scvi_dist_clamp(self.theta[i])
    // Poisson-Gamma mixture: lambda ~ Gamma(theta, mu / theta), x ~ Pois(lambda)
    let lambda = scvi_dist_gamma_sample(theta, mu / theta, seed + i * 53)
    out.push(scvi_dist_poisson_sample(lambda, seed + i * 53 + 11))
  }
  out
}

// ===========================================================================
// Zero-inflated negative binomial
// ===========================================================================

///|
pub struct ScviZeroInflatedNegativeBinomial {
  mu : Array[Double]
  theta : Array[Double]
  zi_logits : Array[Double]
}

///|
pub fn ScviZeroInflatedNegativeBinomial::new(
  mu : Array[Double],
  theta : Array[Double],
  zi_logits : Array[Double],
) -> ScviZeroInflatedNegativeBinomial {
  { mu, theta, zi_logits, }
}

///|
/// Complement logit of the zero-inflation probability.
fn scvi_dist_sigmoid(x : Double) -> Double {
  if x >= 0.0 {
    1.0 / (1.0 + @math.exp(-x))
  } else {
    let e = @math.exp(x)
    e / (1.0 + e)
  }
}

///|
pub fn ScviZeroInflatedNegativeBinomial::log_prob(
  self : ScviZeroInflatedNegativeBinomial,
  x : Array[Double],
) -> Array[Double] {
  let n = self.mu.length()
  let nb = ScviNegativeBinomial::new(self.mu, self.theta)
  let nb_lp = nb.log_prob(x)
  let out = Array(capacity=n)
  for i in 0..<n {
    let pi = scvi_dist_sigmoid(self.zi_logits[i])
    let log_pi = @math.ln(scvi_dist_clamp(pi))
    let log_not_pi = @math.ln(scvi_dist_clamp(1.0 - pi))
    if x[i] == 0.0 {
      let log_zero = log_not_pi + nb_lp[i]
      let a = if log_pi > log_zero { log_pi } else { log_zero }
      let b = if log_pi > log_zero { log_zero } else { log_pi }
      out.push(a + @math.ln(1.0 + @math.exp(b - a)))
    } else {
      out.push(log_not_pi + nb_lp[i])
    }
  }
  out
}

///|
pub fn ScviZeroInflatedNegativeBinomial::mean(
  self : ScviZeroInflatedNegativeBinomial,
) -> Array[Double] {
  let n = self.mu.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    out.push((1.0 - scvi_dist_sigmoid(self.zi_logits[i])) * self.mu[i])
  }
  out
}

///|
pub fn ScviZeroInflatedNegativeBinomial::variance(
  self : ScviZeroInflatedNegativeBinomial,
) -> Array[Double] {
  let n = self.mu.length()
  let nb = ScviNegativeBinomial::new(self.mu, self.theta).variance()
  let out = Array(capacity=n)
  for i in 0..<n {
    let pi = scvi_dist_sigmoid(self.zi_logits[i])
    let mean = (1.0 - pi) * self.mu[i]
    let second = (1.0 - pi) * (nb[i] + self.mu[i] * self.mu[i])
    out.push(second - mean * mean)
  }
  out
}

///|
pub fn ScviZeroInflatedNegativeBinomial::sample(
  self : ScviZeroInflatedNegativeBinomial,
  seed : Int32,
) -> Array[Double] {
  let nb = ScviNegativeBinomial::new(self.mu, self.theta).sample(seed)
  let n = self.mu.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let pi = scvi_dist_sigmoid(self.zi_logits[i])
    let u = scvi_dist_uniform(seed + i * 71 + 3)
    if u < pi {
      out.push(0.0)
    } else {
      out.push(nb[i])
    }
  }
  out
}

// ===========================================================================
// Bernoulli
// ===========================================================================

///|
pub struct ScviBernoulli {
  probs : Array[Double]
}

///|
pub fn ScviBernoulli::new(probs : Array[Double]) -> ScviBernoulli {
  { probs, }
}

///|
pub fn ScviBernoulli::log_prob(
  self : ScviBernoulli,
  x : Array[Double],
) -> Array[Double] {
  let n = self.probs.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let p = self.probs[i]
    let lp = if p < 1.0e-8 {
      1.0e-8
    } else if p > 1.0 - 1.0e-8 {
      1.0 - 1.0e-8
    } else {
      p
    }
    out.push(x[i] * @math.ln(lp) + (1.0 - x[i]) * @math.ln(1.0 - lp))
  }
  out
}

///|
pub fn ScviBernoulli::mean(self : ScviBernoulli) -> Array[Double] {
  self.probs
}

///|
pub fn ScviBernoulli::variance(self : ScviBernoulli) -> Array[Double] {
  let n = self.probs.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    out.push(self.probs[i] * (1.0 - self.probs[i]))
  }
  out
}

///|
pub fn ScviBernoulli::sample(
  self : ScviBernoulli,
  seed : Int32,
) -> Array[Double] {
  let n = self.probs.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let u = scvi_dist_uniform(seed + i * 97)
    if u < self.probs[i] {
      out.push(1.0)
    } else {
      out.push(0.0)
    }
  }
  out
}

// ===========================================================================
// Negative binomial mixture (two components)
// ===========================================================================

///|
pub struct ScviNegativeBinomialMixture {
  mu1 : Array[Double]
  theta1 : Array[Double]
  mu2 : Array[Double]
  theta2 : Array[Double]
  mixture_logits : Array[Double]
}

///|
pub fn ScviNegativeBinomialMixture::new(
  mu1 : Array[Double],
  theta1 : Array[Double],
  mu2 : Array[Double],
  theta2 : Array[Double],
  mixture_logits : Array[Double],
) -> ScviNegativeBinomialMixture {
  { mu1, theta1, mu2, theta2, mixture_logits, }
}

///|
pub fn ScviNegativeBinomialMixture::log_prob(
  self : ScviNegativeBinomialMixture,
  x : Array[Double],
) -> Array[Double] {
  let n = self.mu1.length()
  let lp1 = ScviNegativeBinomial::new(self.mu1, self.theta1).log_prob(x)
  let lp2 = ScviNegativeBinomial::new(self.mu2, self.theta2).log_prob(x)
  let out = Array(capacity=n)
  for i in 0..<n {
    let pi = scvi_dist_sigmoid(self.mixture_logits[i])
    let a = lp1[i] + @math.ln(scvi_dist_clamp(pi))
    let b = lp2[i] + @math.ln(scvi_dist_clamp(1.0 - pi))
    let m = if a > b { a } else { b }
    let l = if a > b { b } else { a }
    out.push(m + @math.ln(1.0 + @math.exp(l - m)))
  }
  out
}

///|
pub fn ScviNegativeBinomialMixture::mean(
  self : ScviNegativeBinomialMixture,
) -> Array[Double] {
  let n = self.mu1.length()
  let out = Array(capacity=n)
  for i in 0..<n {
    let pi = scvi_dist_sigmoid(self.mixture_logits[i])
    out.push(pi * self.mu1[i] + (1.0 - pi) * self.mu2[i])
  }
  out
}

// ===========================================================================
// Mixture of Poisson components
// ===========================================================================

///|
pub struct ScviMixturePoisson {
  rates : Array[Array[Double]] // [component][gene]
  mixture_logits : Array[Array[Double]] // [component][gene]
}

///|
pub fn ScviMixturePoisson::new(
  rates : Array[Array[Double]],
  mixture_logits : Array[Array[Double]],
) -> ScviMixturePoisson {
  { rates, mixture_logits, }
}

///|
pub fn ScviMixturePoisson::log_prob(
  self : ScviMixturePoisson,
  x : Array[Double],
) -> Array[Double] {
  let n_comp = self.rates.length()
  let n_gene = if n_comp == 0 { 0 } else { self.rates[0].length() }
  let out = Array(capacity=n_gene)
  for g in 0..<n_gene {
    let logits = Array(capacity=n_comp)
    for c in 0..<n_comp {
      logits.push(scvi_dist_clamp(self.mixture_logits[c][g]))
    }
    let mut m = logits[0]
    for c in 1..<n_comp {
      if logits[c] > m {
        m = logits[c]
      }
    }
    let mut denom = 0.0
    for c in 0..<n_comp {
      denom += @math.exp(logits[c] - m)
    }
    let log_norm = m + @math.ln(denom)
    let terms = Array(capacity=n_comp)
    for c in 0..<n_comp {
      let mu = scvi_dist_clamp(self.rates[c][g])
      let lp = x[g] * @math.ln(mu) - mu - lgamma(x[g] + 1.0)
      terms.push(lp + logits[c] - log_norm)
    }
    out.push(scvi_dist_log_sum_exp(terms))
  }
  out
}

///|
pub fn ScviMixturePoisson::mean(self : ScviMixturePoisson) -> Array[Double] {
  let n_comp = self.rates.length()
  let n_gene = if n_comp == 0 { 0 } else { self.rates[0].length() }
  let out = Array(capacity=n_gene)
  for g in 0..<n_gene {
    let mut m = self.mixture_logits[0][g]
    for c in 1..<n_comp {
      if self.mixture_logits[c][g] > m {
        m = self.mixture_logits[c][g]
      }
    }
    let mut denom = 0.0
    for c in 0..<n_comp {
      denom += @math.exp(self.mixture_logits[c][g] - m)
    }
    let mut mean = 0.0
    for c in 0..<n_comp {
      let w = @math.exp(self.mixture_logits[c][g] - m) / denom
      mean += w * self.rates[c][g]
    }
    out.push(mean)
  }
  out
}

// ===========================================================================
// Convenience aliases matching scvi-tools constructor names
// ===========================================================================

///|
pub fn scvi_dist_negative_binomial(
  mu : Array[Double],
  theta : Array[Double],
) -> ScviNegativeBinomial {
  ScviNegativeBinomial::new(mu, theta)
}

///|
pub fn scvi_dist_zero_inflated_negative_binomial(
  mu : Array[Double],
  theta : Array[Double],
  zi_logits : Array[Double],
) -> ScviZeroInflatedNegativeBinomial {
  ScviZeroInflatedNegativeBinomial::new(mu, theta, zi_logits)
}

///|
pub fn scvi_dist_poisson(rate : Array[Double]) -> ScviPoisson {
  ScviPoisson::new(rate)
}

///|
pub fn scvi_dist_bernoulli(probs : Array[Double]) -> ScviBernoulli {
  ScviBernoulli::new(probs)
}

///|
pub fn scvi_dist_nb_mixture(
  mu1 : Array[Double],
  theta1 : Array[Double],
  mu2 : Array[Double],
  theta2 : Array[Double],
  mixture_logits : Array[Double],
) -> ScviNegativeBinomialMixture {
  ScviNegativeBinomialMixture::new(mu1, theta1, mu2, theta2, mixture_logits)
}

///|
pub fn scvi_dist_mixture_poisson(
  rates : Array[Array[Double]],
  mixture_logits : Array[Array[Double]],
) -> ScviMixturePoisson {
  ScviMixturePoisson::new(rates, mixture_logits)
}
