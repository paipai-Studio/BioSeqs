///|
pub struct Node {
  node_id : String
  label : String
  attributes : Map[String, String]
}

///|
pub fn Node::new(node_id : String, label : String) -> Node {
  Node::{ node_id, label, attributes: Map([], capacity=5) }
}

///|
pub fn Node::to_string(self : Node) -> String {
  self.node_id + ": " + self.label
}

///|
pub struct Edge {
  source : String
  target : String
  weight : Double
  label : String
  attributes : Map[String, String]
}

///|
pub fn Edge::new(source : String, target : String, weight : Double) -> Edge {
  Edge::{ source, target, weight, label: "", attributes: Map([], capacity=5) }
}

///|
pub fn Edge::to_string(self : Edge) -> String {
  self.source + " --[" + self.weight.to_string() + "]--> " + self.target
}

///|
pub enum GraphType {
  Undirected
  Directed
}

///|
pub fn GraphType::new_undirected() -> GraphType {
  Undirected
}

///|
pub fn GraphType::new_directed() -> GraphType {
  Directed
}

///|
pub fn graph_type_is_undirected(gt : GraphType) -> Bool {
  match gt {
    Undirected => true
    Directed => false
  }
}

///|
pub fn GraphType::to_string(self : GraphType) -> String {
  match self {
    Undirected => "Undirected"
    Directed => "Directed"
  }
}

///|
pub struct Graph {
  nodes : Array[Node]
  edges : Array[Edge]
  graph_type : GraphType
  name : String
}

///|
pub fn Graph::new(name : String, graph_type : GraphType) -> Graph {
  Graph::{ nodes: Array::new(), edges: Array::new(), graph_type, name }
}

///|
pub fn graph_add_node(graph : Graph, node : Node) -> Unit {
  graph.nodes.push(node)
}

///|
pub fn graph_add_edge(graph : Graph, edge : Edge) -> Unit {
  graph.edges.push(edge)
}

///|
pub fn graph_add_edge_by_id(
  graph : Graph,
  source_id : String,
  target_id : String,
  weight : Double,
) -> Unit {
  graph.edges.push(Edge::new(source_id, target_id, weight))
}

///|
pub fn graph_get_node_by_id(graph : Graph, node_id : String) -> Node? {
  let mut i = 0
  while i < graph.nodes.length() {
    let node = graph.nodes[i]
    if node.node_id == node_id {
      return Some(node)
    }
    i = i + 1
  }
  None
}

///|
pub fn graph_get_edges_from(graph : Graph, node_id : String) -> Array[Edge] {
  let results : Array[Edge] = Array::new()
  let mut i = 0
  while i < graph.edges.length() {
    let edge = graph.edges[i]
    if edge.source == node_id {
      results.push(edge)
    }
    i = i + 1
  }
  results
}

///|
pub fn graph_get_edges_to(graph : Graph, node_id : String) -> Array[Edge] {
  let results : Array[Edge] = Array::new()
  let mut i = 0
  while i < graph.edges.length() {
    let edge = graph.edges[i]
    if edge.target == node_id {
      results.push(edge)
    }
    i = i + 1
  }
  results
}

///|
pub fn graph_get_neighbors(graph : Graph, node_id : String) -> Array[String] {
  let results : Array[String] = Array::new()
  let mut i = 0
  while i < graph.edges.length() {
    let edge = graph.edges[i]
    if edge.source == node_id {
      results.push(edge.target)
    } else if edge.target == node_id &&
      graph_type_is_undirected(graph.graph_type) {
      results.push(edge.source)
    }
    i = i + 1
  }
  results
}

///|
pub fn graph_get_degree(graph : Graph, node_id : String) -> Int {
  let mut degree = 0
  let mut i = 0
  while i < graph.edges.length() {
    let edge = graph.edges[i]
    if edge.source == node_id {
      degree = degree + 1
    }
    if edge.target == node_id &&
      (graph_type_is_undirected(graph.graph_type) || edge.source != node_id) {
      degree = degree + 1
    }
    i = i + 1
  }
  degree
}

///|
pub fn graph_node_count(graph : Graph) -> Int {
  graph.nodes.length()
}

///|
pub fn graph_edge_count(graph : Graph) -> Int {
  graph.edges.length()
}

///|
pub fn graph_has_edge(
  graph : Graph,
  source_id : String,
  target_id : String,
) -> Bool {
  let mut i = 0
  while i < graph.edges.length() {
    let edge = graph.edges[i]
    if edge.source == source_id && edge.target == target_id {
      return true
    }
    if graph_type_is_undirected(graph.graph_type) &&
      edge.source == target_id &&
      edge.target == source_id {
      return true
    }
    i = i + 1
  }
  false
}

///|
pub fn graph_shortest_path(
  graph : Graph,
  start_id : String,
  end_id : String,
) -> Array[String]? {
  if start_id == end_id {
    return Some([start_id])
  }

  let visited : Map[String, Bool] = Map([], capacity=graph.nodes.length())
  let parent : Map[String, String] = Map([], capacity=graph.nodes.length())
  let queue : Array[String] = Array::new()

  queue.push(start_id)
  visited[start_id] = true

  let mut found = false
  while queue.length() > 0 && !found {
    let current = queue.remove(0)
    let neighbors = graph_get_neighbors(graph, current)

    let mut i = 0
    while i < neighbors.length() && !found {
      let neighbor = neighbors[i]
      let visited_neighbor = visited.get(neighbor)
      if visited_neighbor.is_none() || !visited_neighbor.unwrap() {
        visited[neighbor] = true
        parent[neighbor] = current
        if neighbor == end_id {
          found = true
        } else {
          queue.push(neighbor)
        }
      }
      i = i + 1
    }
  }

  if !found {
    return None
  }

  let path : Array[String] = Array::new()
  let mut current = end_id
  while current != start_id {
    path.push(current)
    let p = parent.get(current)
    if p.is_none() {
      return None
    }
    current = p.unwrap()
  }
  path.push(start_id)

  let reversed : Array[String] = Array::new()
  let mut i = path.length() - 1
  while i >= 0 {
    reversed.push(path[i])
    i = i - 1
  }

  Some(reversed)
}

///|
pub fn graph_connected_components(graph : Graph) -> Array[Array[String]] {
  let visited : Map[String, Bool] = Map([], capacity=graph.nodes.length())
  let components : Array[Array[String]] = Array::new()

  let mut i = 0
  while i < graph.nodes.length() {
    let node = graph.nodes[i]
    let visited_node = visited.get(node.node_id)
    if visited_node.is_none() || !visited_node.unwrap() {
      let component : Array[String] = Array::new()
      let queue : Array[String] = Array::new()

      queue.push(node.node_id)
      visited[node.node_id] = true

      while queue.length() > 0 {
        let current = queue.remove(0)
        component.push(current)

        let neighbors = graph_get_neighbors(graph, current)
        let mut j = 0
        while j < neighbors.length() {
          let neighbor = neighbors[j]
          let visited_neighbor = visited.get(neighbor)
          if visited_neighbor.is_none() || !visited_neighbor.unwrap() {
            visited[neighbor] = true
            queue.push(neighbor)
          }
          j = j + 1
        }
      }

      components.push(component)
    }
    i = i + 1
  }

  components
}

///|
pub fn graph_summarize(graph : Graph) -> String {
  let mut result = "Graph Summary\n"
  result = result + "=============\n"
  result = result + "Name: " + graph.name + "\n"
  result = result + "Type: " + graph.graph_type.to_string() + "\n"
  result = result + "Nodes: " + graph_node_count(graph).to_string() + "\n"
  result = result + "Edges: " + graph_edge_count(graph).to_string() + "\n\n"

  let degrees : Array[Int] = Array::new()
  let mut i = 0
  while i < graph.nodes.length() {
    let node = graph.nodes[i]
    degrees.push(graph_get_degree(graph, node.node_id))
    i = i + 1
  }

  if degrees.length() > 0 {
    let mut min_degree = degrees[0]
    let mut max_degree = degrees[0]
    let mut sum_degree = 0
    i = 0
    while i < degrees.length() {
      sum_degree = sum_degree + degrees[i]
      if degrees[i] < min_degree {
        min_degree = degrees[i]
      }
      if degrees[i] > max_degree {
        max_degree = degrees[i]
      }
      i = i + 1
    }
    let avg_degree = sum_degree.to_double() / degrees.length().to_double()

    result = result + "Degree statistics:\n"
    result = result + "  - Min: " + min_degree.to_string() + "\n"
    result = result + "  - Max: " + max_degree.to_string() + "\n"
    result = result + "  - Avg: " + avg_degree.to_string() + "\n"
  }

  let components = graph_connected_components(graph)
  result = result +
    "\nConnected components: " +
    components.length().to_string() +
    "\n"

  result
}

///|
pub fn graph_to_dot(graph : Graph) -> String {
  let mut result = if graph_type_is_undirected(graph.graph_type) {
    "graph "
  } else {
    "digraph "
  }
  result = result + graph.name + " {\n"

  let mut i = 0
  while i < graph.nodes.length() {
    let node = graph.nodes[i]
    result = result + "  " + node.node_id + " [label=\"" + node.label + "\"];\n"
    i = i + 1
  }

  i = 0
  while i < graph.edges.length() {
    let edge = graph.edges[i]
    let arrow = if graph_type_is_undirected(graph.graph_type) {
      " -- "
    } else {
      " -> "
    }
    let weight_label = if edge.weight != 1.0 {
      " [label=\"" + edge.weight.to_string() + "\"]"
    } else {
      ""
    }
    result = result +
      "  " +
      edge.source +
      arrow +
      edge.target +
      weight_label +
      ";\n"
    i = i + 1
  }

  result = result + "}"
  result
}

///|
pub fn create_example_graph() -> Graph {
  let graph = Graph::new("Protein-Protein Interaction", GraphType::Undirected)

  graph_add_node(graph, Node::new("TP53", "Tumor protein p53"))
  graph_add_node(graph, Node::new("BRCA1", "Breast cancer type 1"))
  graph_add_node(graph, Node::new("BRCA2", "Breast cancer type 2"))
  graph_add_node(graph, Node::new("ATM", "Ataxia telangiectasia mutated"))
  graph_add_node(graph, Node::new("CHEK2", "Checkpoint kinase 2"))
  graph_add_node(graph, Node::new("MDM2", "Mouse double minute 2"))
  graph_add_node(graph, Node::new("BAX", "BCL2-associated X protein"))
  graph_add_node(graph, Node::new("BCL2", "B-cell lymphoma 2"))

  graph_add_edge_by_id(graph, "TP53", "BRCA1", 0.95)
  graph_add_edge_by_id(graph, "TP53", "BRCA2", 0.90)
  graph_add_edge_by_id(graph, "TP53", "ATM", 0.88)
  graph_add_edge_by_id(graph, "TP53", "MDM2", 0.99)
  graph_add_edge_by_id(graph, "TP53", "BAX", 0.92)
  graph_add_edge_by_id(graph, "BRCA1", "BRCA2", 0.85)
  graph_add_edge_by_id(graph, "BRCA1", "ATM", 0.91)
  graph_add_edge_by_id(graph, "BRCA1", "CHEK2", 0.87)
  graph_add_edge_by_id(graph, "BRCA2", "ATM", 0.89)
  graph_add_edge_by_id(graph, "ATM", "CHEK2", 0.93)
  graph_add_edge_by_id(graph, "BAX", "BCL2", 0.94)

  graph
}

///|
pub fn create_pathway_graph() -> Graph {
  let graph = Graph::new("Cell Cycle Pathway", GraphType::Directed)

  graph_add_node(graph, Node::new("G1", "G1 phase"))
  graph_add_node(graph, Node::new("S", "S phase"))
  graph_add_node(graph, Node::new("G2", "G2 phase"))
  graph_add_node(graph, Node::new("M", "M phase"))
  graph_add_node(graph, Node::new("CDK4", "Cyclin-dependent kinase 4"))
  graph_add_node(graph, Node::new("CDK2", "Cyclin-dependent kinase 2"))
  graph_add_node(graph, Node::new("CDK1", "Cyclin-dependent kinase 1"))
  graph_add_node(graph, Node::new("RB", "Retinoblastoma protein"))
  graph_add_node(graph, Node::new("E2F", "E2F transcription factor"))

  graph_add_edge_by_id(graph, "CDK4", "RB", 0.95)
  graph_add_edge_by_id(graph, "RB", "E2F", 0.90)
  graph_add_edge_by_id(graph, "E2F", "S", 0.98)
  graph_add_edge_by_id(graph, "G1", "CDK4", 0.85)
  graph_add_edge_by_id(graph, "CDK2", "S", 0.92)
  graph_add_edge_by_id(graph, "S", "G2", 0.99)
  graph_add_edge_by_id(graph, "G2", "CDK1", 0.93)
  graph_add_edge_by_id(graph, "CDK1", "M", 0.97)

  graph
}
