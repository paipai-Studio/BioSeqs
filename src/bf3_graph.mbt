// Bioframe third-round: overlap graph operations (items 23-25)
// overlap_graph, connected_components, overlap_matrix

///|
pub struct Bf3OverlapEdge {
  from_idx : Int32
  to_idx : Int32
  overlap_start : Int32
  overlap_end : Int32
  overlap_length : Int32
}

///|
pub fn bf3_overlap_graph(
  intervals : Array[(String, Int32, Int32)],
) -> Array[Bf3OverlapEdge] {
  // Build overlap graph: edges between overlapping intervals
  let edges : Array[Bf3OverlapEdge] = []
  for i in 0..<intervals.length() {
    let (chrom_i, start_i, end_i) = intervals[i]
    let next_i = i + 1
    for j in next_i..<intervals.length() {
      let (chrom_j, start_j, end_j) = intervals[j]
      if chrom_i == chrom_j {
        let ov_start = if start_i > start_j { start_i } else { start_j }
        let ov_end = if end_i < end_j { end_i } else { end_j }
        if ov_end > ov_start {
          edges.push({
            from_idx: i,
            to_idx: j,
            overlap_start: ov_start,
            overlap_end: ov_end,
            overlap_length: ov_end - ov_start,
          })
        }
      }
    }
  }
  edges
}

///|
pub struct Bf3Component {
  component_id : Int32
  members : Array[Int32] // indices of intervals in this component
  size : Int32
}

///|
pub fn bf3_connected_components(
  intervals : Array[(String, Int32, Int32)],
) -> Array[Bf3Component] {
  // Find connected components of overlapping intervals using union-find
  let n = intervals.length()
  let parent : Array[Int32] = []
  for i in 0..<n {
    parent.push(i)
  }
  fn find(parent : Array[Int32], x : Int32) -> Int32 {
    let mut root = x
    while parent[root] != root {
      root = parent[root]
    }
    // Path compression
    let mut curr = x
    while parent[curr] != root {
      let next = parent[curr]
      parent[curr] = root
      curr = next
    }
    root
  }
  // Build edges from overlaps
  for i in 0..<n {
    let (chrom_i, start_i, end_i) = intervals[i]
    let next_i = i + 1
    for j in next_i..<n {
      let (chrom_j, start_j, end_j) = intervals[j]
      if chrom_i == chrom_j && start_j < end_i && end_j > start_i {
        let ri = find(parent, i)
        let rj = find(parent, j)
        if ri != rj {
          parent[ri] = rj
        }
      }
    }
  }
  // Group by root
  let component_map : Map[Int32, Array[Int32]] = Map([])
  for i in 0..<n {
    let root = find(parent, i)
    let members = component_map.get(root).unwrap_or([])
    members.push(i)
    component_map[root] = members
  }
  let components : Array[Bf3Component] = []
  let mut cid = 0
  component_map.each(fn(_key, members) {
    components.push({ component_id: cid, members, size: members.length(), })
    cid = cid + 1
  })
  components
}

///|
pub fn bf3_overlap_matrix(
  intervals_a : Array[(String, Int32, Int32)],
  intervals_b : Array[(String, Int32, Int32)],
) -> Array[Array[Bool]] {
  // Build boolean overlap matrix: matrix[i][j] = true if A[i] overlaps B[j]
  let na = intervals_a.length()
  let nb = intervals_b.length()
  let matrix : Array[Array[Bool]] = []
  for i in 0..<na {
    let (chrom_i, start_i, end_i) = intervals_a[i]
    let row : Array[Bool] = []
    for j in 0..<nb {
      let (chrom_j, start_j, end_j) = intervals_b[j]
      let overlaps = chrom_i == chrom_j && start_j < end_i && end_j > start_i
      row.push(overlaps)
    }
    matrix.push(row)
  }
  matrix
}

///|
pub fn bf3_overlap_matrix_weighted(
  intervals_a : Array[(String, Int32, Int32)],
  intervals_b : Array[(String, Int32, Int32)],
) -> Array[Array[Int32]] {
  // Build weighted overlap matrix: matrix[i][j] = overlap length (0 if no overlap)
  let na = intervals_a.length()
  let nb = intervals_b.length()
  let matrix : Array[Array[Int32]] = []
  for i in 0..<na {
    let (chrom_i, start_i, end_i) = intervals_a[i]
    let row : Array[Int32] = []
    for j in 0..<nb {
      let (chrom_j, start_j, end_j) = intervals_b[j]
      if chrom_i == chrom_j {
        let ov_start = if start_i > start_j { start_i } else { start_j }
        let ov_end = if end_i < end_j { end_i } else { end_j }
        row.push(if ov_end > ov_start { ov_end - ov_start } else { 0 })
      } else {
        row.push(0)
      }
    }
    matrix.push(row)
  }
  matrix
}
