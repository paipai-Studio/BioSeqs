///|
/// Tree layout and visualization: rectangular and circular layouts.
/// Port of ete3's TreeStyle rectangular and circular (radial) layouts.
/// Generates 2D node positions that can be rendered as ASCII art.

// =========================================================================
// Layout Data Structures
// =========================================================================

///|
/// A node position in a rectangular layout.
pub(all) struct NodePosRect {
  name : String
  x : Double
  y : Double
  is_leaf : Bool
  branch_length : Double
  depth : Int
  confidence : Double?
  children : Array[Int]
}

///|
/// Rectangular tree layout result.
pub(all) struct TreeLayoutRect {
  nodes : Array[NodePosRect]
  width : Double
  height : Double
  n_leaves : Int
}

///|
/// A node position in a circular (radial) layout.
pub(all) struct NodePosCirc {
  name : String
  radius : Double
  angle : Double
  x : Double
  y : Double
  is_leaf : Bool
  branch_length : Double
  depth : Int
  confidence : Double?
  children : Array[Int]
}

///|
/// Circular tree layout result.
pub(all) struct TreeLayoutCirc {
  nodes : Array[NodePosCirc]
  max_radius : Double
  n_leaves : Int
  start_angle : Double
  end_angle : Double
}

// =========================================================================
// Rectangular Layout
// =========================================================================

///|
/// Compute rectangular layout. x = cumulative branch length, y = leaf position.
pub fn phylo_layout_rectangular(tree : Clade) -> TreeLayoutRect {
  let nodes : Array[NodePosRect] = Array(capacity=0)
  let _ = layout_rect_rec(tree, 0.0, 0, nodes)
  let _ = layout_rect_assign_y(nodes, 0, 0)
  let mut max_x = 0.0
  let mut max_y = 0.0
  for n in nodes {
    if n.x > max_x {
      max_x = n.x
    }
    if n.y > max_y {
      max_y = n.y
    }
  }
  TreeLayoutRect::{
    nodes,
    width: max_x,
    height: max_y,
    n_leaves: rect_leaf_count(nodes),
  }
}

///|
fn layout_rect_rec(
  node : Clade,
  x : Double,
  depth : Int,
  out : Array[NodePosRect]
) -> Int {
  let bl = match node.branch_length {
    Some(b) => b
    None => 0.0
  }
  let my_x = x + bl
  let name = match node.name {
    Some(n) => n
    None => ""
  }
  let my_idx = out.length()
  out.push(NodePosRect::{
    name,
    x: my_x,
    y: 0.0,
    is_leaf: node.is_terminal(),
    branch_length: bl,
    depth,
    confidence: node.confidence,
    children: [],
  })
  let child_indices : Array[Int] = Array(capacity=0)
  for ch in node.clades {
    child_indices.push(layout_rect_rec(ch, my_x, depth + 1, out))
  }
  out[my_idx] = NodePosRect::{ ..out[my_idx], children: child_indices }
  my_idx
}

///|
fn layout_rect_assign_y(
  nodes : Array[NodePosRect],
  idx : Int,
  leaf_counter : Int
) -> Int {
  if nodes[idx].is_leaf {
    nodes[idx] = NodePosRect::{ ..nodes[idx], y: leaf_counter.to_double() }
    return leaf_counter + 1
  }
  let mut next_leaf = leaf_counter
  for ci in nodes[idx].children {
    next_leaf = layout_rect_assign_y(nodes, ci, next_leaf)
  }
  let children = nodes[idx].children
  if children.length() > 0 {
    let mut sum_y = 0.0
    for ci in children {
      sum_y = sum_y + nodes[ci].y
    }
    let avg_y = sum_y / children.length().to_double()
    nodes[idx] = NodePosRect::{ ..nodes[idx], y: avg_y }
  }
  next_leaf
}

///|
fn rect_leaf_count(nodes : Array[NodePosRect]) -> Int {
  let mut c = 0
  for n in nodes {
    if n.is_leaf {
      c = c + 1
    }
  }
  c
}

// =========================================================================
// Circular Layout
// =========================================================================

///|
/// Compute circular layout. Leaves evenly distributed around the circle.
pub fn phylo_layout_circular(tree : Clade) -> TreeLayoutCirc {
  phylo_layout_circular_range(tree, 0.0, 2.0 * 3.141592653589793)
}

///|
/// Compute circular layout with custom angle range.
pub fn phylo_layout_circular_range(
  tree : Clade,
  start_angle : Double,
  end_angle : Double
) -> TreeLayoutCirc {
  let nodes : Array[NodePosCirc] = Array(capacity=0)
  let _ = layout_circ_rec(tree, 0.0, 0, nodes)
  let n_leaves = circ_leaf_count(nodes)
  if n_leaves == 0 {
    return TreeLayoutCirc::{
      nodes,
      max_radius: 0.0,
      n_leaves: 0,
      start_angle,
      end_angle,
    }
  }
  let _ = layout_circ_assign_angle(
    nodes, 0, 0, n_leaves, start_angle, end_angle
  )
  let mut max_r = 0.0
  for n in nodes {
    if n.radius > max_r {
      max_r = n.radius
    }
  }
  TreeLayoutCirc::{
    nodes,
    max_radius: max_r,
    n_leaves,
    start_angle,
    end_angle,
  }
}

///|
fn layout_circ_rec(
  node : Clade,
  radius : Double,
  depth : Int,
  out : Array[NodePosCirc]
) -> Int {
  let bl = match node.branch_length {
    Some(b) => b
    None => 0.0
  }
  let my_r = radius + bl
  let name = match node.name {
    Some(n) => n
    None => ""
  }
  let my_idx = out.length()
  out.push(NodePosCirc::{
    name,
    radius: my_r,
    angle: 0.0,
    x: 0.0,
    y: 0.0,
    is_leaf: node.is_terminal(),
    branch_length: bl,
    depth,
    confidence: node.confidence,
    children: [],
  })
  let child_indices : Array[Int] = Array(capacity=0)
  for ch in node.clades {
    child_indices.push(layout_circ_rec(ch, my_r, depth + 1, out))
  }
  out[my_idx] = NodePosCirc::{ ..out[my_idx], children: child_indices }
  my_idx
}

///|
fn layout_circ_assign_angle(
  nodes : Array[NodePosCirc],
  idx : Int,
  leaf_counter : Int,
  n_leaves : Int,
  start_angle : Double,
  end_angle : Double
) -> Int {
  let angle_span = end_angle - start_angle
  if nodes[idx].is_leaf {
    let angle = start_angle +
      leaf_counter.to_double() / n_leaves.to_double() * angle_span
    let x = nodes[idx].radius * @math.cos(angle)
    let y = nodes[idx].radius * @math.sin(angle)
    nodes[idx] = NodePosCirc::{ ..nodes[idx], angle, x, y }
    return leaf_counter + 1
  }
  let mut next_leaf = leaf_counter
  for ci in nodes[idx].children {
    next_leaf = layout_circ_assign_angle(
      nodes, ci, next_leaf, n_leaves, start_angle, end_angle
    )
  }
  let children = nodes[idx].children
  if children.length() > 0 {
    let mut sum_a = 0.0
    for ci in children {
      sum_a = sum_a + nodes[ci].angle
    }
    let avg_a = sum_a / children.length().to_double()
    let x = nodes[idx].radius * @math.cos(avg_a)
    let y = nodes[idx].radius * @math.sin(avg_a)
    nodes[idx] = NodePosCirc::{ ..nodes[idx], angle: avg_a, x, y }
  }
  next_leaf
}

///|
fn circ_leaf_count(nodes : Array[NodePosCirc]) -> Int {
  let mut c = 0
  for n in nodes {
    if n.is_leaf {
      c = c + 1
    }
  }
  c
}


// =========================================================================
// ASCII Rendering
// =========================================================================

///|
/// Render rectangular layout as ASCII art.
pub fn phylo_render_rectangular(
  layout : TreeLayoutRect,
  width : Int,
  height : Int
) -> String {
  if width <= 0 || height <= 0 || layout.nodes.length() == 0 {
    return ""
  }
  let grid = make_grid(width, height)
  let x_scale = if layout.width > 0.0 {
    (width - 1).to_double() / layout.width
  } else {
    1.0
  }
  let y_scale = if layout.height > 0.0 {
    (height - 1).to_double() / layout.height
  } else {
    1.0
  }
  for i in 0..<layout.nodes.length() {
    let node = layout.nodes[i]
    for ci in node.children {
      let child = layout.nodes[ci]
      let y1 = (node.y * y_scale).to_int()
      let y2 = (child.y * y_scale).to_int()
      let x1 = (node.x * x_scale).to_int()
      let x2 = (child.x * x_scale).to_int()
      let xs = if x1 < x2 { x1 } else { x2 }
      let xe = if x1 < x2 { x2 } else { x1 }
      for x in xs..<(xe + 1) {
        if x >= 0 && x < width && y2 >= 0 && y2 < height {
          if grid[y2][x] == ' ' {
            grid[y2][x] = '-'
          }
        }
      }
      let ys = if y1 < y2 { y1 } else { y2 }
      let ye = if y1 < y2 { y2 } else { y1 }
      for y in ys..<(ye + 1) {
        if x1 >= 0 && x1 < width && y >= 0 && y < height {
          if grid[y][x1] == ' ' {
            grid[y][x1] = '|'
          }
        }
      }
    }
  }
  for node in layout.nodes {
    let gx = (node.x * x_scale).to_int()
    let gy = (node.y * y_scale).to_int()
    if gx >= 0 && gx < width && gy >= 0 && gy < height {
      if node.is_leaf && node.name != "" {
        grid[gy][gx] = first_char(node.name)
      } else if not(node.is_leaf) && grid[gy][gx] == ' ' {
        grid[gy][gx] = '+'
      }
    }
  }
  grid_to_string(grid, width, height)
}

///|
/// Render circular layout as ASCII art.
pub fn phylo_render_circular(
  layout : TreeLayoutCirc,
  width : Int,
  height : Int
) -> String {
  if width <= 0 || height <= 0 || layout.nodes.length() == 0 {
    return ""
  }
  let grid = make_grid(width, height)
  let mut max_abs = 0.0
  for n in layout.nodes {
    let ax = n.x.abs()
    let ay = n.y.abs()
    if ax > max_abs {
      max_abs = ax
    }
    if ay > max_abs {
      max_abs = ay
    }
  }
  if max_abs == 0.0 {
    max_abs = 1.0
  }
  let cx = width / 2
  let cy = height / 2
  let scale = {
    let sx = if cx > 0 { (cx - 1).to_double() / max_abs } else { 0.0 }
    let sy = if cy > 0 { (cy - 1).to_double() / max_abs } else { 0.0 }
    if sx < sy {
      sx
    } else {
      sy
    }
  }
  for i in 0..<layout.nodes.length() {
    let node = layout.nodes[i]
    for ci in node.children {
      let child = layout.nodes[ci]
      let x1 = cx + (node.x * scale).to_int()
      let y1 = cy + (node.y * scale).to_int()
      let x2 = cx + (child.x * scale).to_int()
      let y2 = cy + (child.y * scale).to_int()
      draw_line_grid(grid, x1, y1, x2, y2, width, height)
    }
  }
  for node in layout.nodes {
    let gx = cx + (node.x * scale).to_int()
    let gy = cy + (node.y * scale).to_int()
    if gx >= 0 && gx < width && gy >= 0 && gy < height {
      if node.is_leaf && node.name != "" {
        grid[gy][gx] = first_char(node.name)
      } else if not(node.is_leaf) && grid[gy][gx] == ' ' {
        grid[gy][gx] = '+'
      }
    }
  }
  grid_to_string(grid, width, height)
}

// =========================================================================
// Grid helpers
// =========================================================================

///|
fn make_grid(width : Int, height : Int) -> Array[Array[Char]] {
  let grid : Array[Array[Char]] = Array(capacity=height)
  for _ in 0..<height {
    let row : Array[Char] = Array(capacity=width)
    for _ in 0..<width {
      row.push(' ')
    }
    grid.push(row)
  }
  grid
}

///|
fn grid_to_string(grid : Array[Array[Char]], width : Int, height : Int) -> String {
  let sb = StringBuilder()
  for y in 0..<height {
    for x in 0..<width {
      sb.write_char(grid[y][x])
    }
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// Bresenham line algorithm.
fn draw_line_grid(
  grid : Array[Array[Char]],
  x0 : Int,
  y0 : Int,
  x1 : Int,
  y1 : Int,
  width : Int,
  height : Int
) -> Unit {
  let mut dx = x1 - x0
  let dy = y1 - y0
  let mut sx = 1
  if dx < 0 {
    dx = -dx
    sx = -1
  }
  let mut sy = 1
  let abs_dy = if dy < 0 { -dy } else { dy }
  if dy < 0 {
    sy = -1
  }
  let mut err = dx - abs_dy
  let mut x = x0
  let mut y = y0
  let mut done = false
  while not(done) {
    if x >= 0 && x < width && y >= 0 && y < height {
      if grid[y][x] == ' ' {
        grid[y][x] = '/'
      }
    }
    if x == x1 && y == y1 {
      done = true
    }
    let e2 = 2 * err
    if e2 > -abs_dy {
      err = err - abs_dy
      x = x + sx
    }
    if e2 < dx {
      err = err + dx
      y = y + sy
    }
  }
}

///|
fn first_char(s : String) -> Char {
  if s.length() == 0 {
    return '?'
  }
  s.unsafe_get(0).unsafe_to_char()
}
