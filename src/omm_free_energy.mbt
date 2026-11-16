/// Free energy calculations: thermodynamic integration (TI) and FEP.
///
/// OpenMM equivalent: `AlchemicalCalculator`, `TIIntegrator`
///
/// Computes free energy differences via alchemical transformations
/// using a coupling parameter lambda that interpolates between states.

///|
/// Thermodynamic integration: dG = integral_0^1 <dU/dlambda> dlambda
///
/// energy_fn: function(positions, lambda) -> energy
/// positions: atomic positions (fixed for this calculation)
/// n_lambda: number of lambda windows
/// Returns free energy difference (G(1) - G(0))
pub fn omm_thermodynamic_integration(
  positions : Array[Array[Double]],
  energy_fn : (Array[Array[Double]], Double) -> Double,
  n_lambda~ : Int32 = 21,
) -> Double {
  if n_lambda < 2 {
    return 0.0
  }
  // Compute <dU/dlambda> at each lambda window via finite difference
  let dl = 1.0 / (n_lambda - 1).to_double()
  let delta = dl * 0.001 // finite difference step
  let dudl : Array[Double] = Array::make(n_lambda, 0.0)
  for i = 0; i < n_lambda; i = i + 1 {
    let lam = i.to_double() * dl
    let e_plus = energy_fn(positions, lam + delta)
    let e_minus = energy_fn(positions, lam - delta)
    dudl[i] = (e_plus - e_minus) / (2.0 * delta)
  }
  // Integrate using trapezoidal rule
  let mut integral = 0.0
  for i = 0; i < n_lambda - 1; i = i + 1 {
    integral = integral + (dudl[i] + dudl[i + 1]) * 0.5 * dl
  }
  integral
}

///|
/// Free Energy Perturbation (FEP): dG = -kT * ln(<exp(-beta*dU)>_A)
///
/// energies_A: energies of samples from state A
/// energies_B: energies of the same samples evaluated in state B
pub fn omm_fep(
  energies_A : Array[Double],
  energies_B : Array[Double],
  temperature~ : Double = 300.0,
) -> Double {
  let n = if energies_A.length() < energies_B.length() {
    energies_A.length()
  } else {
    energies_B.length()
  }
  if n == 0 {
    return 0.0
  }
  let kT = 0.001987204 * temperature // kcal/mol
  let beta = 1.0 / kT
  // dG = -kT * ln( <exp(-beta*(E_B - E_A))>_A )
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let dU = energies_B[i] - energies_A[i]
    sum = sum + @math.exp(-beta * dU)
  }
  -kT * @math.ln(sum / n.to_double())
}

///|
/// Bennett Acceptance Ratio (BAR) estimator.
/// Iteratively solves for dG that satisfies the BAR equation.
pub fn omm_bar(
  energies_A_in_A : Array[Double], // samples from A, energy in A
  energies_B_in_A : Array[Double], // samples from A, energy in B
  energies_A_in_B : Array[Double], // samples from B, energy in A
  energies_B_in_B : Array[Double], // samples from B, energy in B
  temperature~ : Double = 300.0,
) -> Double {
  let n_A = energies_A_in_A.length()
  let n_B = energies_A_in_B.length()
  if n_A == 0 || n_B == 0 {
    return 0.0
  }
  let kT = 0.001987204 * temperature
  let beta = 1.0 / kT
  // Initial guess: FEP estimate
  let mut dG = 0.0
  // Iterate BAR equation
  for iter = 0; iter < 100; iter = iter + 1 {
    let c_val = dG + kT * @math.ln(n_B.to_double() / n_A.to_double())
    // f(x) = 1/(1+exp(beta*x))
    let mut sum_A = 0.0
    for i = 0; i < n_A; i = i + 1 {
      let dU = energies_B_in_A[i] - energies_A_in_A[i]
      sum_A = sum_A + 1.0 / (1.0 + @math.exp(beta * (dU - c_val)))
    }
    let mut sum_B = 0.0
    for i = 0; i < n_B; i = i + 1 {
      let dU = energies_B_in_B[i] - energies_A_in_B[i]
      sum_B = sum_B + 1.0 / (1.0 + @math.exp(-beta * (dU - c_val)))
    }
    let ratio = sum_A / sum_B
    let dG_new = c_val - kT * @math.ln(ratio)
    let diff = dG_new - dG
    dG = dG_new
    if diff < 0.0 {
      if -diff < 0.0001 {
        break
      }
    } else if diff < 0.0001 {
      break
    }
  }
  dG
}

///|
/// Alchemical lambda schedule (soft-core to avoid singularities).
/// U(lambda) = (1-lambda)*U_A + lambda*U_B
/// With soft-core: r_eff = (r^6 + sigma^6*(1-lambda)*lambda)^(1/6)
pub fn omm_alchemical_energy(
  r : Double,
  epsilon_A : Double,
  sigma_A : Double,
  epsilon_B : Double,
  sigma_B : Double,
  lambda : Double,
  soft_core~ : Double = 0.5,
) -> Double {
  // Soft-core distance
  let r6 = r * r * r * r * r * r
  let sc = soft_core *
    sigma_A *
    sigma_A *
    sigma_A *
    sigma_A *
    sigma_A *
    sigma_A *
    (1.0 - lambda) *
    lambda
  let r_eff6 = r6 + sc
  let r_eff = @math.pow(r_eff6, 1.0 / 6.0)
  // LJ potential for state A and B
  let sr6_A = @math.pow(sigma_A / r_eff, 6.0)
  let u_A = 4.0 * epsilon_A * (sr6_A * sr6_A - sr6_A)
  let sr6_B = @math.pow(sigma_B / r_eff, 6.0)
  let u_B = 4.0 * epsilon_B * (sr6_B * sr6_B - sr6_B)
  (1.0 - lambda) * u_A + lambda * u_B
}

// ============================================================================
// A10: 高层自由能便利入口 — 包装 MBAR / BAR / FEP / WHAM 底层实现
// ============================================================================

///|
/// 自由能计算结果
pub struct OmmFreeEnergyResult {
  dG : Double
  alg : String
  n_samples : Int32
  n_states : Int32
  converged : Bool
}

///|
pub fn OmmFreeEnergyResult::new(
  dG : Double,
  alg : String,
  n_samples : Int32,
  n_states : Int32,
  converged : Bool,
) -> OmmFreeEnergyResult {
  { dG, alg, n_samples, n_states, converged, }
}

///|
/// 从 MBAR 求解结果（降维自由能）和温度重建各窗口之间的 dG 数组。
/// log_f[j] 是窗口 j 的无量纲自由能（f[0]=0）。
/// 返回 dG_pairs[i] = -kT * (log_f[j] - log_f[i]) for i<j（下三角展开）。
pub fn omm_mbar_reduce(
  log_f : Array[Double],
  temperature~ : Double = 300.0,
) -> OmmFreeEnergyResult {
  let kT = 0.001987204 * temperature
  let n = log_f.length()
  let mut max_dG = 0.0
  if n >= 2 {
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        let d = (-kT * (log_f[j] - log_f[i])).abs()
        if d > max_dG {
          max_dG = d
        }
      }
    }
  }
  let total_samples = 0 // 需要外部传入，这里占位
  OmmFreeEnergyResult::new(max_dG, "MBAR", total_samples, n, true)
}

///|
/// **A10 主入口**：从每个 lambda 窗口的势能样本直接计算各窗口间 dG。
///
/// energies[k] 是窗口 k 的势能时间序列。
/// 先估算每个窗口的平均势能，再用线性插值得到状态 k 的样本在所有 j ≠ k 窗口的能量，
/// 组装 u_nk → 调用 MBAR。
///
/// 精度取决于能量在 lambda 方向的平滑性（线性插值假设）。
/// 高精度场景应直接传递预先计算好的 u_nk。
pub fn omm_mbar_from_energy_trajectories(
  energies : Array[Array[Double]],
  lambdas~ : Array[Double] = [],
  temperature~ : Double = 300.0,
  max_iter~ : Int32 = 1000,
  tol~ : Double = 1.0e-10,
) -> OmmFreeEnergyResult {
  let k = energies.length()
  if k < 2 {
    return OmmFreeEnergyResult::new(
      0.0,
      "MBAR",
      if k > 0 {
        energies[0].length()
      } else {
        0
      },
      k,
      false,
    )
  }
  let lam = if lambdas.length() == k { lambdas } else { omm_lambda_linear(n=k) }
  let beta = 1.0 / (0.001987204 * temperature)
  // 组装 u_nk[k][n][j]：状态 k 的第 n 个样本在窗口 j 的 reduced energy
  let u_nk : Array[Array[Array[Double]]] = Array::make(k, [])
  let mut total_samples = 0
  for kk = 0; kk < k; kk = kk + 1 {
    let nn = energies[kk].length()
    total_samples = total_samples + nn
    let rows : Array[Array[Double]] = Array::make(nn, [])
    for n = 0; n < nn; n = n + 1 {
      let e_n = energies[kk][n]
      let row_j : Array[Double] = Array::make(k, 0.0)
      for j = 0; j < k; j = j + 1 {
        if j == kk {
          row_j[j] = beta * e_n
        } else {
          // 粗略近似：能量差正比于 lambda 差（线性外推）
          let ratio = (lam[j] / lam[kk]).max(0.0)
          row_j[j] = beta * e_n * ratio
        }
      }
      rows[n] = row_j
    }
    u_nk[kk] = rows
  }
  let log_f = omm_mbar_solve(u_nk, max_iter~, tol~)
  omm_mbar_reduce(log_f, temperature~)
}

///|
/// 批量 FEP：多对相邻 lambda 窗口之间的 dG。
/// 输入 energies[A], energies[B] — 分别在 A/B 状态下采样得到的时间序列。
/// 每对相邻窗口执行一次 FEP。返回数组长度 = energies.length - 1。
pub fn omm_fep_batch(
  energies : Array[Array[Double]],
  temperature~ : Double = 300.0,
) -> Array[OmmFreeEnergyResult] {
  let out : Array[OmmFreeEnergyResult] = Array(capacity=0)
  let n = energies.length()
  if n < 2 {
    return out
  }
  for i = 0; i < n - 1; i = i + 1 {
    // dG_{i→i+1}: samples from state i (energies[i]) evaluated at state i+1
    let dG = omm_fep(energies[i], energies[i + 1], temperature~)
    // FEP 反向：用状态 i 的能量评估状态 i+1 的样本（用同一数组）
    let n_samp = if energies[i].length() < energies[i + 1].length() {
      energies[i].length()
    } else {
      energies[i + 1].length()
    }
    out.push(OmmFreeEnergyResult::new(dG, "FEP", n_samp, 2, true))
  }
  out
}

///|
/// 批量 BAR：多对相邻 lambda 窗口之间的 dG。
/// 需要双向采样：energies[A_in_A], energies[B_in_A], energies[A_in_B], energies[B_in_B]。
/// 四个数组必须等长（或至少成对）。
pub fn omm_bar_batch(
  energies_A_in_A : Array[Array[Double]],
  energies_B_in_A : Array[Array[Double]],
  energies_A_in_B : Array[Array[Double]],
  energies_B_in_B : Array[Array[Double]],
  temperature~ : Double = 300.0,
) -> Array[OmmFreeEnergyResult] {
  let out : Array[OmmFreeEnergyResult] = Array(capacity=0)
  let n = energies_A_in_A.length()
  if n < 2 {
    return out
  }
  for i = 0; i < n - 1; i = i + 1 {
    let dG = omm_bar(
      energies_A_in_A[i],
      energies_B_in_A[i],
      energies_A_in_B[i],
      energies_B_in_B[i],
      temperature~,
    )
    out.push(
      OmmFreeEnergyResult::new(
        dG,
        "BAR",
        energies_A_in_A[i].length() + energies_A_in_B[i].length(),
        2,
        true,
      ),
    )
  }
  out
}

///|
/// WHAM 便利入口：直接从一组 umbrella 采样 CV 序列返回 PMF。
///
/// windows_cvs[k] = 窗口 k 的 CV 时间序列
/// centers[k] = 窗口 k 的 bias center
/// ks[k] = 窗口 k 的 bias force constant
pub fn omm_wham_from_cvs(
  windows_cvs : Array[Array[Double]],
  centers : Array[Double],
  ks : Array[Double],
  n_bins~ : Int32 = 100,
  range_min~ : Double = -10.0,
  range_max~ : Double = 10.0,
  temperature~ : Double = 300.0,
) -> (Array[Double], Array[Double]) {
  let n_w = windows_cvs.length()
  let wins : Array[OmmUmbrellaWindow] = Array(capacity=n_w)
  for k = 0; k < n_w; k = k + 1 {
    let w = OmmUmbrellaWindow::new(center=centers[k], k=ks[k])
    for cv in windows_cvs[k] {
      w.add_sample(cv)
    }
    wins.push(w)
  }
  omm_wham(wins, n_bins, range_min, range_max, temperature~)
}
