/// seq_cluster.mbt — sequence clustering (Butina analogue, item B5)
///
/// K-mer multiset Tanimoto distance + the shared Butina leader
/// clustering core (the same implementation used for molecules in
/// rdkit_cluster_butina.mbt).

///|
/// K-mer multiset of a string (list with repetitions).
fn seq_cluster_kmers(s : String, k : Int32) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  let n = s.length()
  if k <= 0 || n < k {
    return out
  }
  let mut i = 0
  while i <= n - k {
    out.push(s[i:i + k].to_owned())
    i = i + 1
  }
  out
}

///|
/// Multiset Tanimoto similarity of two sequences' k-mer bags:
/// min-count intersection over max-count union.
pub fn seq_kmer_tanimoto(a : Seq, b : Seq, k~ : Int32 = 3) -> Double {
  let ka = seq_cluster_kmers(a.to_string(), k)
  let kb = seq_cluster_kmers(b.to_string(), k)
  if ka.length() == 0 && kb.length() == 0 {
    return 0.0
  }
  let ca : Map[String, Int32] = Map([], capacity=16)
  for km in ka {
    ca[km] = match ca.get(km) {
      Some(c) => c + 1
      None => 1
    }
  }
  let cb : Map[String, Int32] = Map([], capacity=16)
  for km in kb {
    cb[km] = match cb.get(km) {
      Some(c) => c + 1
      None => 1
    }
  }
  let mut inter = 0
  for km, cnt in ca {
    match cb.get(km) {
      Some(c2) =>
        if cnt < c2 {
          inter = inter + cnt
        } else {
          inter = inter + c2
        }
      None => ()
    }
  }
  let union = ka.length() + kb.length() - inter
  if union == 0 {
    0.0
  } else {
    inter.to_double() / union.to_double()
  }
}

///|
/// Butina clustering of sequences by k-mer Tanimoto distance
/// (1 - similarity). Returns clusters as arrays of record indices; the
/// first index of each cluster is its leader.
pub fn seq_butina_cluster(
  records : Array[SeqRecord],
  k~ : Int32 = 3,
  threshold~ : Double = 0.3,
) -> Array[Array[Int32]] {
  let n = records.length()
  let dist : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::make(n, 0.0)
    dist.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let d = 1.0 - seq_kmer_tanimoto(records[i].seq, records[j].seq, k~)
      dist[i][j] = d
      dist[j][i] = d
    }
  }
  rdkit_butina_cluster(dist, threshold)
}

///|
/// Largest cluster size for a quick overview of the clustering result.
pub fn seq_butina_largest_cluster(clusters : Array[Array[Int32]]) -> Int32 {
  let mut best = 0
  for c in clusters {
    if c.length() > best {
      best = c.length()
    }
  }
  best
}
