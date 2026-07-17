///|
/// K-means Clustering for Biological Data Analysis
pub type Point = Array[Double]

///|
pub struct KMeans {
  k : Int
  centroids : Array[Point]
  labels : Array[Int]
  max_iter : Int
  tol : Double
}

///|
pub fn KMeans::new(k : Int, max_iter : Int, tol : Double) -> KMeans {
  KMeans::{ k, centroids: [], labels: [], max_iter, tol }
}

///|
fn double_max(a : Double, b : Double) -> Double {
  if a > b {
    a
  } else {
    b
  }
}

///|
pub fn euclidean_distance(a : Point, b : Point) -> Double {
  let mut sum = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    let diff = a[i] - b[i]
    sum = sum + diff * diff
  }
  sum.sqrt()
}

///|
pub fn manhattan_distance(a : Point, b : Point) -> Double {
  let mut sum = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    sum = sum + (a[i] - b[i]).abs()
  }
  sum
}

///|
pub fn cosine_distance(a : Point, b : Point) -> Double {
  let mut dot = 0.0
  let mut norm_a = 0.0
  let mut norm_b = 0.0

  for i = 0; i < a.length(); i = i + 1 {
    dot = dot + a[i] * b[i]
    norm_a = norm_a + a[i] * a[i]
    norm_b = norm_b + b[i] * b[i]
  }

  let denom = norm_a.sqrt() * norm_b.sqrt()
  if denom == 0.0 {
    1.0
  } else {
    1.0 - dot / denom
  }
}

///|
pub fn assign_clusters(
  data : Array[Point],
  centroids : Array[Point],
) -> Array[Int] {
  let labels = []

  for point in data {
    let mut min_dist = 1000000000000000000000000000000000.0
    let mut min_idx = 0

    for i = 0; i < centroids.length(); i = i + 1 {
      let dist = euclidean_distance(point, centroids[i])
      if dist < min_dist {
        min_dist = dist
        min_idx = i
      }
    }

    labels.push(min_idx)
  }

  labels
}

///|
pub fn update_centroids(
  data : Array[Point],
  labels : Array[Int],
  k : Int,
) -> Array[Point] {
  let n = data.length()
  let dim = if n > 0 { data[0].length() } else { 0 }

  let centroids = []
  let counts = Array::make(k, 0)

  for idx = 0; idx < k; idx = idx + 1 {
    centroids.push(Array::make(dim, 0.0))
  }

  for i = 0; i < n; i = i + 1 {
    let label = labels[i]
    counts[label] = counts[label] + 1
    for j = 0; j < dim; j = j + 1 {
      centroids[label][j] = centroids[label][j] + data[i][j]
    }
  }

  for i = 0; i < k; i = i + 1 {
    if counts[i] > 0 {
      for j = 0; j < dim; j = j + 1 {
        centroids[i][j] = centroids[i][j] / counts[i].to_double()
      }
    }
  }

  centroids
}

///|
pub fn KMeans::fit(self : KMeans, data : Array[Point]) -> KMeans {
  let n = data.length()
  if n == 0 {
    return KMeans::{ ..self, centroids: [], labels: [] }
  }

  let dim = data[0].length()

  let mut centroids = []
  for i = 0; i < self.k; i = i + 1 {
    let idx = i * (n / self.k)
    centroids.push(data[if idx >= n { n - 1 } else { idx }].copy())
  }

  let mut labels = []

  for iter = 0; iter < self.max_iter; iter = iter + 1 {
    let new_labels = assign_clusters(data, centroids)
    let new_centroids = update_centroids(data, new_labels, self.k)

    let mut max_shift = 0.0
    for i = 0; i < self.k; i = i + 1 {
      let shift = euclidean_distance(centroids[i], new_centroids[i])
      if shift > max_shift {
        max_shift = shift
      }
    }

    labels = new_labels
    centroids = new_centroids

    if max_shift < self.tol {
      break
    }
  }

  KMeans::{ ..self, centroids, labels }
}

///|
pub fn KMeans::predict(self : KMeans, data : Array[Point]) -> Array[Int] {
  assign_clusters(data, self.centroids)
}

///|
pub fn KMeans::predict_single(self : KMeans, point : Point) -> Int {
  let mut min_dist = 1000000000000000000000000000000000.0
  let mut min_idx = 0

  for i = 0; i < self.centroids.length(); i = i + 1 {
    let dist = euclidean_distance(point, self.centroids[i])
    if dist < min_dist {
      min_dist = dist
      min_idx = i
    }
  }

  min_idx
}

///|
pub fn KMeans::get_centroids(self : KMeans) -> Array[Point] {
  self.centroids
}

///|
pub fn KMeans::get_labels(self : KMeans) -> Array[Int] {
  self.labels
}

///|
pub fn KMeans::inertia(self : KMeans, data : Array[Point]) -> Double {
  let mut total = 0.0

  for i = 0; i < data.length(); i = i + 1 {
    let label = self.labels[i]
    let dist = euclidean_distance(data[i], self.centroids[label])
    total = total + dist * dist
  }

  total
}

///|
pub fn silhouette_score(data : Array[Point], labels : Array[Int]) -> Double {
  let n = data.length()
  if n == 0 {
    return 0.0
  }

  let mut score_sum = 0.0
  let mut max_label = 0
  for l in labels {
    if l > max_label {
      max_label = l
    }
  }
  let k = max_label + 1

  for i = 0; i < n; i = i + 1 {
    let label = labels[i]
    let mut a = 0.0
    let mut a_count = 0
    let mut b = 1000000000000000000000000000000000.0

    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }

      let dist = euclidean_distance(data[i], data[j])

      if labels[j] == label {
        a = a + dist
        a_count = a_count + 1
      } else {
        let mut cluster_dist = 0.0
        let mut cluster_count = 0
        for m = 0; m < n; m = m + 1 {
          if labels[m] == labels[j] {
            cluster_dist = cluster_dist + euclidean_distance(data[i], data[m])
            cluster_count = cluster_count + 1
          }
        }
        let avg_dist = if cluster_count > 0 {
          cluster_dist / cluster_count.to_double()
        } else {
          0.0
        }
        if avg_dist < b {
          b = avg_dist
        }
      }
    }

    a = if a_count > 0 { a / a_count.to_double() } else { 0.0 }

    let s = if double_max(a, b) == 0.0 {
      0.0
    } else {
      (b - a) / double_max(a, b)
    }
    score_sum = score_sum + s
  }

  score_sum / n.to_double()
}

///|
pub fn cluster_gene_expression(data : Array[Point], k : Int) -> KMeans {
  let kmeans = KMeans::new(k, 100, 0.000001)
  kmeans.fit(data)
}

///|
pub fn find_optimal_k(data : Array[Point], min_k : Int, max_k : Int) -> Int {
  let mut best_k = min_k
  let mut best_silhouette = -1.0

  for k_val = min_k; k_val <= max_k; k_val = k_val + 1 {
    let kmeans = KMeans::new(k_val, 100, 0.000001).fit(data)
    let score = silhouette_score(data, kmeans.labels)
    if score > best_silhouette {
      best_silhouette = score
      best_k = k_val
    }
  }

  best_k
}
