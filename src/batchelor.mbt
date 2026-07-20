///|
/// Batchelor - Single-cell batch correction methods.
/// Inspired by Bioconductor batchelor package.
/// Provides batch correction methods including rescaleBatches, MNN, and fastMNN.

///|
pub struct BatchCorrectionResult {
  corrected : Array[Array[Double]]
  batch_indices : Array[Int]
  mnn_pairs : Array[(Int, Int)]
  var_explained : Array[Double]
}

///|
pub fn BatchCorrectionResult::new(
  corrected : Array[Array[Double]],
  batch_indices : Array[Int],
  mnn_pairs : Array[(Int, Int)],
  var_explained : Array[Double],
) -> BatchCorrectionResult {
  BatchCorrectionResult::{
    corrected,
    batch_indices,
    mnn_pairs,
    var_explained,
  }
}

///|
pub fn rescale_batches(
  batches : Array[Array[Array[Double]]],
  batch_names : Array[String],
) -> BatchCorrectionResult {
  let n_batches = batches.length()
  if n_batches == 0 {
    return BatchCorrectionResult::new([], [], [], [])
  }

  let n_cells_per_batch : Array[Int] = Array::new()
  let n_genes = batches[0][0].length()

  let mut total_cells = 0
  let mut i = 0
  while i < n_batches {
    n_cells_per_batch.push(batches[i].length())
    total_cells = total_cells + batches[i].length()
    i = i + 1
  }

  let pooled_means : Array[Double] = Array::make(n_genes, 0.0)
  let mut j = 0
  while j < n_batches {
    let batch = batches[j]
    let mut k = 0
    while k < batch.length() {
      let cell = batch[k]
      let mut l = 0
      while l < n_genes {
        pooled_means[l] = pooled_means[l] + cell[l]
        l = l + 1
      }
      k = k + 1
    }
    j = j + 1
  }

  let mut m = 0
  while m < n_genes {
    pooled_means[m] = pooled_means[m] / total_cells.to_double()
    m = m + 1
  }

  let batch_means : Array[Array[Double]] = Array::new()
  let mut n = 0
  while n < n_batches {
    let batch = batches[n]
    let means : Array[Double] = Array::make(n_genes, 0.0)
    let mut o = 0
    while o < batch.length() {
      let cell = batch[o]
      let mut p = 0
      while p < n_genes {
        means[p] = means[p] + cell[p]
        p = p + 1
      }
      o = o + 1
    }
    let mut q = 0
    while q < n_genes {
      means[q] = means[q] / batch.length().to_double()
      q = q + 1
    }
    batch_means.push(means)
    n = n + 1
  }

  let corrected : Array[Array[Double]] = Array::new()
  let batch_indices : Array[Int] = Array::new()
  let mut r = 0
  while r < n_batches {
    let batch = batches[r]
    let means = batch_means[r]
    let mut s = 0
    while s < batch.length() {
      let cell = batch[s]
      let corrected_cell : Array[Double] = Array::new()
      let mut t = 0
      while t < n_genes {
        let correction_factor = if means[t] != 0.0 {
          pooled_means[t] / means[t]
        } else {
          1.0
        }
        corrected_cell.push(cell[t] * correction_factor)
        t = t + 1
      }
      corrected.push(corrected_cell)
      batch_indices.push(r)
      s = s + 1
    }
    r = r + 1
  }

  BatchCorrectionResult::new(corrected, batch_indices, [], [])
}

///|
pub fn find_mutual_nearest_neighbors(
  batch1 : Array[Array[Double]],
  batch2 : Array[Array[Double]],
  k : Int,
) -> Array[(Int, Int)] {
  let n_cells1 = batch1.length()
  let n_cells2 = batch2.length()

  let n1_to_n2 : Array[Array[Int]] = Array::new()
  let mut i = 0
  while i < n_cells1 {
    let cell = batch1[i]
    let distances : Array[(Double, Int)] = Array::new()
    let mut j = 0
    while j < n_cells2 {
      let other = batch2[j]
      let mut dist = 0.0
      let mut kk = 0
      while kk < cell.length() {
        dist = dist + (cell[kk] - other[kk]) * (cell[kk] - other[kk])
        kk = kk + 1
      }
      distances.push((dist.sqrt(), j))
      j = j + 1
    }

    let mut l = 0
    while l < distances.length() - 1 {
      let mut m = l + 1
      while m < distances.length() {
        if distances[l].0 > distances[m].0 {
          let tmp = distances[l]
          distances[l] = distances[m]
          distances[m] = tmp
        }
        m = m + 1
      }
      l = l + 1
    }

    let nn : Array[Int] = Array::new()
    let mut n_neighbors = 0
    while n_neighbors < k && n_neighbors < distances.length() {
      nn.push(distances[n_neighbors].1)
      n_neighbors = n_neighbors + 1
    }
    n1_to_n2.push(nn)
    i = i + 1
  }

  let n2_to_n1 : Array[Array[Int]] = Array::new()
  let mut n = 0
  while n < n_cells2 {
    let cell = batch2[n]
    let distances : Array[(Double, Int)] = Array::new()
    let mut o = 0
    while o < n_cells1 {
      let other = batch1[o]
      let mut dist = 0.0
      let mut kk = 0
      while kk < cell.length() {
        dist = dist + (cell[kk] - other[kk]) * (cell[kk] - other[kk])
        kk = kk + 1
      }
      distances.push((dist.sqrt(), o))
      o = o + 1
    }

    let mut p = 0
    while p < distances.length() - 1 {
      let mut q = p + 1
      while q < distances.length() {
        if distances[p].0 > distances[q].0 {
          let tmp = distances[p]
          distances[p] = distances[q]
          distances[q] = tmp
        }
        q = q + 1
      }
      p = p + 1
    }

    let nn : Array[Int] = Array::new()
    let mut n_neighbors = 0
    while n_neighbors < k && n_neighbors < distances.length() {
      nn.push(distances[n_neighbors].1)
      n_neighbors = n_neighbors + 1
    }
    n2_to_n1.push(nn)
    n = n + 1
  }

  let mnn_pairs : Array[(Int, Int)] = Array::new()
  let mut r = 0
  while r < n_cells1 {
    let neighbors = n1_to_n2[r]
    let mut s = 0
    while s < neighbors.length() {
      let n2_idx = neighbors[s]
      let reverse_neighbors = n2_to_n1[n2_idx]
      let mut t = 0
      while t < reverse_neighbors.length() {
        if reverse_neighbors[t] == r {
          let already_exists = mnn_pairs.any(fn(pair : (Int, Int)) {
            pair.0 == r && pair.1 == n2_idx
          })
          if !already_exists {
            mnn_pairs.push((r, n2_idx))
          }
          break
        }
        t = t + 1
      }
      s = s + 1
    }
    r = r + 1
  }

  mnn_pairs
}

///|
pub fn fast_mnn(
  batches : Array[Array[Array[Double]]],
  batch_names : Array[String],
  k : Int,
  d : Int,
) -> BatchCorrectionResult {
  let n_batches = batches.length()
  if n_batches == 0 {
    return BatchCorrectionResult::new([], [], [], [])
  }

  let n_genes = batches[0][0].length()

  let all_cells : Array[Array[Double]] = Array::new()
  let batch_indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n_batches {
    let batch = batches[i]
    let mut j = 0
    while j < batch.length() {
      all_cells.push(batch[j])
      batch_indices.push(i)
      j = j + 1
    }
    i = i + 1
  }

  let n_cells = all_cells.length()

  let means : Array[Double] = Array::make(n_genes, 0.0)
  let mut kk = 0
  while kk < n_cells {
    let cell = all_cells[kk]
    let mut ll = 0
    while ll < n_genes {
      means[ll] = means[ll] + cell[ll]
      ll = ll + 1
    }
    kk = kk + 1
  }
  let mut mm = 0
  while mm < n_genes {
    means[mm] = means[mm] / n_cells.to_double()
    mm = mm + 1
  }

  let centered : Array[Array[Double]] = Array::new()
  let mut nn = 0
  while nn < n_cells {
    let cell = all_cells[nn]
    let centered_cell : Array[Double] = Array::new()
    let mut oo = 0
    while oo < n_genes {
      centered_cell.push(cell[oo] - means[oo])
      oo = oo + 1
    }
    centered.push(centered_cell)
    nn = nn + 1
  }

  let cov_matrix : Array[Array[Double]] = Array::make(n_genes, Array::make(n_genes, 0.0))
  let mut pp = 0
  while pp < n_cells {
    let cell = centered[pp]
    let mut qq = 0
    while qq < n_genes {
      let mut rr = 0
      while rr < n_genes {
        cov_matrix[qq][rr] = cov_matrix[qq][rr] + cell[qq] * cell[rr]
        rr = rr + 1
      }
      qq = qq + 1
    }
    pp = pp + 1
  }
  let mut ss = 0
  while ss < n_genes {
    let mut tt = 0
    while tt < n_genes {
      cov_matrix[ss][tt] = cov_matrix[ss][tt] / (n_cells - 1).to_double()
      tt = tt + 1
    }
    ss = ss + 1
  }

  let eigenvalues : Array[Double] = Array::make(n_genes, 0.0)
  let eigenvectors : Array[Array[Double]] = Array::make(n_genes, Array::make(n_genes, 0.0))
  let mut uu = 0
  while uu < n_genes {
    eigenvalues[uu] = cov_matrix[uu][uu]
    eigenvectors[uu][uu] = 1.0
    uu = uu + 1
  }

  let mut vv = 0
  while vv < 100 {
    let mut max_diff = 0.0
    let mut ww = 0
    while ww < n_genes - 1 {
      let mut xx = ww + 1
      while xx < n_genes {
        let p = cov_matrix[ww][xx]
        let d_val = eigenvalues[xx] - eigenvalues[ww]
        let c = if d_val.abs() < 0.000001 {
          0.000001
        } else {
          d_val
        }
        let t = p / c
        let cos_val = 1.0 / (1.0 + t * t).sqrt()
        let sin_val = t * cos_val

        let mut yy = 0
        while yy < n_genes {
          let old_wy = eigenvectors[ww][yy]
          let old_xy = eigenvectors[xx][yy]
          eigenvectors[ww][yy] = old_wy * cos_val - old_xy * sin_val
          eigenvectors[xx][yy] = old_wy * sin_val + old_xy * cos_val
          yy = yy + 1
        }

        let mut zz = 0
        while zz < n_genes {
          let old_wz = cov_matrix[ww][zz]
          let old_xz = cov_matrix[xx][zz]
          cov_matrix[ww][zz] = old_wz * cos_val - old_xz * sin_val
          cov_matrix[xx][zz] = old_wz * sin_val + old_xz * cos_val
          zz = zz + 1
        }

        let mut aa = 0
        while aa < n_genes {
          let old_aw = cov_matrix[aa][ww]
          let old_ax = cov_matrix[aa][xx]
          cov_matrix[aa][ww] = old_aw * cos_val - old_ax * sin_val
          cov_matrix[aa][xx] = old_aw * sin_val + old_ax * cos_val
          aa = aa + 1
        }

        eigenvalues[ww] = cov_matrix[ww][ww]
        eigenvalues[xx] = cov_matrix[xx][xx]

        max_diff = max_diff.max(p.abs())
        xx = xx + 1
      }
      ww = ww + 1
    }
    if max_diff < 0.000001 {
      break
    }
    vv = vv + 1
  }

  let mut bb = 0
  while bb < eigenvalues.length() - 1 {
    let mut cc = bb + 1
    while cc < eigenvalues.length() {
      if eigenvalues[bb] < eigenvalues[cc] {
        let tmp = eigenvalues[bb]
        eigenvalues[bb] = eigenvalues[cc]
        eigenvalues[cc] = tmp

        let mut dd = 0
        while dd < n_genes {
          let tmp_vec = eigenvectors[bb][dd]
          eigenvectors[bb][dd] = eigenvectors[cc][dd]
          eigenvectors[cc][dd] = tmp_vec
          dd = dd + 1
        }
      }
      cc = cc + 1
    }
    bb = bb + 1
  }

  let top_components = if d < n_genes { d } else { n_genes }
  let pca_result : Array[Array[Double]] = Array::new()
  let mut ee = 0
  while ee < n_cells {
    let cell = centered[ee]
    let pc : Array[Double] = Array::new()
    let mut ff = 0
    while ff < top_components {
      let mut score = 0.0
      let mut gg = 0
      while gg < n_genes {
        score = score + cell[gg] * eigenvectors[ff][gg]
        gg = gg + 1
      }
      pc.push(score)
      ff = ff + 1
    }
    pca_result.push(pc)
    ee = ee + 1
  }

  let mut total_variance = 0.0
  for val in eigenvalues {
    total_variance = total_variance + val
  }
  let var_explained : Array[Double] = Array::new()
  let mut hh = 0
  while hh < top_components {
    var_explained.push(eigenvalues[hh] / total_variance)
    hh = hh + 1
  }

  let pca_batches : Array[Array[Array[Double]]] = Array::new()
  let mut ii = 0
  while ii < n_batches {
    let batch : Array[Array[Double]] = Array::new()
    let mut jj = 0
    while jj < n_cells {
      if batch_indices[jj] == ii {
        batch.push(pca_result[jj])
      }
      jj = jj + 1
    }
    pca_batches.push(batch)
    ii = ii + 1
  }

  let all_mnn_pairs : Array[(Int, Int)] = Array::new()
  let mut kk = 0
  while kk < n_batches - 1 {
    let ll = kk + 1
    let mnn = find_mutual_nearest_neighbors(pca_batches[kk], pca_batches[ll], k)
    for pair in mnn {
      all_mnn_pairs.push(pair)
    }
    kk = kk + 1
  }

  let correction_vectors : Array[Array[Double]] = Array::make(n_cells, Array::make(top_components, 0.0))

  let mut mm = 0
  while mm < n_batches - 1 {
    let nn = mm + 1
    let batch1 = pca_batches[mm]
    let batch2 = pca_batches[nn]

    let mnn = find_mutual_nearest_neighbors(batch1, batch2, k)

    let batch1_indices : Array[Int] = Array::new()
    let batch2_indices : Array[Int] = Array::new()
    for pair in mnn {
      batch1_indices.push(pair.0)
      batch2_indices.push(pair.1)
    }

    if batch1_indices.length() > 0 {
      let delta : Array[Double] = Array::make(top_components, 0.0)
      let mut oo = 0
      while oo < batch1_indices.length() {
        let idx1 = batch1_indices[oo]
        let idx2 = batch2_indices[oo]
        let mut pp = 0
        while pp < top_components {
          delta[pp] = delta[pp] + batch2[idx2][pp] - batch1[idx1][pp]
          pp = pp + 1
        }
        oo = oo + 1
      }
      let mut qq = 0
      while qq < top_components {
        delta[qq] = delta[qq] / batch1_indices.length().to_double()
        qq = qq + 1
      }

      let batch1_global_indices : Array[Int] = Array::new()
      let mut rr = 0
      let mut count = 0
      while rr < n_cells {
        if batch_indices[rr] == mm {
          if batch1_indices.contains(count) {
            batch1_global_indices.push(rr)
          }
          count = count + 1
        }
        rr = rr + 1
      }

      let sigma = 0.1
      let mut ss = 0
      while ss < n_cells {
        if batch_indices[ss] == mm {
          let mut weight_sum = 0.0
          let weighted_delta : Array[Double] = Array::make(top_components, 0.0)

          let mut tt = 0
          while tt < batch1_global_indices.length() {
            let mnn_idx = batch1_global_indices[tt]
            let mut dist = 0.0
            let mut uu = 0
            while uu < top_components {
              dist = dist + (pca_result[ss][uu] - pca_result[mnn_idx][uu]) * (pca_result[ss][uu] - pca_result[mnn_idx][uu])
              uu = uu + 1
            }
            let weight = exp(-dist / (2.0 * sigma * sigma))
            weight_sum = weight_sum + weight

            let mut vv = 0
            while vv < top_components {
              weighted_delta[vv] = weighted_delta[vv] + weight * delta[vv]
              vv = vv + 1
            }
            tt = tt + 1
          }

          if weight_sum > 0.0 {
            let mut ww = 0
            while ww < top_components {
              correction_vectors[ss][ww] = correction_vectors[ss][ww] + weighted_delta[ww] / weight_sum
              ww = ww + 1
            }
          }
        }
        ss = ss + 1
      }
    }
    mm = mm + 1
  }

  let corrected : Array[Array[Double]] = Array::new()
  let mut xx = 0
  while xx < n_cells {
    let corrected_cell : Array[Double] = Array::new()
    let mut yy = 0
    while yy < top_components {
      corrected_cell.push(pca_result[xx][yy] + correction_vectors[xx][yy])
      yy = yy + 1
    }
    corrected.push(corrected_cell)
    xx = xx + 1
  }

  BatchCorrectionResult::new(corrected, batch_indices, all_mnn_pairs, var_explained)
}

///|
pub fn batchelor_create_example_data(
  n_batches : Int,
  n_cells_per_batch : Int,
  n_genes : Int,
  batch_effect : Double,
) -> (Array[Array[Array[Double]]], Array[String]) {
  let batches : Array[Array[Array[Double]]] = Array::new()
  let batch_names : Array[String] = Array::new()

  let mut seed = 12345
  let mut i = 0
  while i < n_batches {
    let batch : Array[Array[Double]] = Array::new()
    batch_names.push("Batch_" + i.to_string())

    let mut j = 0
    while j < n_cells_per_batch {
      let cell : Array[Double] = Array::new()

      let mut k = 0
      while k < n_genes {
        let mut x = 0.0
        let mut l = 0
        while l < 12 {
          seed = seed * 9301 + 49297
          let r = (seed % 233280).to_double() / 233280.0
          x = x + r
          l = l + 1
        }
        x = x - 6.0

        let base_expression = exp(x)
        let batch_factor = 1.0 + i.to_double() * batch_effect * (if k < n_genes / 2 { 1.0 } else { -1.0 }) * 0.1

        cell.push(base_expression * batch_factor)
        k = k + 1
      }

      batch.push(cell)
      j = j + 1
    }

    batches.push(batch)
    i = i + 1
  }

  (batches, batch_names)
}

///|
pub fn compute_batch_mixing_score(corrected : Array[Array[Double]], batch_indices : Array[Int]) -> Double {
  let n_cells = corrected.length()
  let n_dims = if n_cells > 0 { corrected[0].length() } else { 0 }

  if n_cells < 2 {
    return 0.0
  }

  let unique_batches : Array[Int] = Array::new()
  for idx in batch_indices {
    let mut exists = false
    for u in unique_batches {
      if u == idx {
        exists = true
        break
      }
    }
    if !exists {
      unique_batches.push(idx)
    }
  }
  let n_batches = unique_batches.length()

  let mut score = 0.0
  let mut i = 0
  while i < n_cells {
    let distances : Array[(Double, Int)] = Array::new()
    let mut j = 0
    while j < n_cells {
      if i != j {
        let mut dist = 0.0
        let mut k = 0
        while k < n_dims {
          dist = dist + (corrected[i][k] - corrected[j][k]) * (corrected[i][k] - corrected[j][k])
          k = k + 1
        }
        distances.push((dist.sqrt(), j))
      }
      j = j + 1
    }

    let mut l = 0
    while l < distances.length() - 1 {
      let mut m = l + 1
      while m < distances.length() {
        if distances[l].0 > distances[m].0 {
          let tmp = distances[l]
          distances[l] = distances[m]
          distances[m] = tmp
        }
        m = m + 1
      }
      l = l + 1
    }

    let n_neighbors = (n_cells / 10).max(5)
    let mut same_batch_count = 0
    let mut n = 0
    while n < n_neighbors && n < distances.length() {
      let neighbor_idx = distances[n].1
      if batch_indices[neighbor_idx] == batch_indices[i] {
        same_batch_count = same_batch_count + 1
      }
      n = n + 1
    }

    score = score + (1.0 - same_batch_count.to_double() / n_neighbors.to_double())
    i = i + 1
  }

  score / n_cells.to_double()
}
