///|
/// VennDiagram - Venn Diagram Set Operations and Region Calculations.
/// Inspired by Bioconductor's VennDiagram package.
/// Provides set operations, Venn region calculations, and statistics
/// for comparing up to 5 sets (genes, proteins, etc.) commonly used in
/// bioinformatics differential analysis and functional annotation.

///|
/// Region in a Venn diagram identified by a bitmask.
/// For n sets, each region is identified by which sets contain it.
/// The bitmask uses the i-th bit to indicate membership in set i.
pub struct VennRegion {
  /// Unique identifier as a bitmask (e.g., 0b101 means in sets 0 and 2)
  id : Int
  /// Elements in this region
  elements : Array[String]
  /// Count of elements
  count : Int
  /// Description of which sets this region belongs to
  description : String
}

///|
/// Result of a Venn diagram analysis.
pub struct VennDiagramResult {
  /// Number of sets
  n_sets : Int
  /// Set names
  set_names : Array[String]
  /// Set sizes (number of elements in each set)
  set_sizes : Array[Int]
  /// All unique elements across all sets
  universe : Array[String]
  /// Regions of the Venn diagram
  regions : Array[VennRegion]
  /// Total elements in the universe
  total_elements : Int
  /// Euler characteristic (for layout validation)
  euler_char : Double
}

///|
/// Create a Venn diagram result from multiple sets.
/// Each set is an array of element strings (e.g., gene symbols).
///
/// Parameters:
///   - sets: Array of sets (each set is an array of strings)
///   - names: Names for each set
///
/// Returns: VennDiagramResult with region calculations.
pub fn venn_diagram(sets : Array[Array[String]], names? : Array[String] = []) -> VennDiagramResult {
  let n_sets = sets.length()
  let set_names = if names.length() == n_sets {
    names
  } else {
    let default = Array::make(n_sets, "")
    let mut i = 0
    while i < n_sets {
      default[i] = "Set\{i + 1}"
      i = i + 1
    }
    default
  }
  // Compute set sizes
  let set_sizes = Array::make(n_sets, 0)
  let mut i = 0
  while i < n_sets {
    set_sizes[i] = sets[i].length()
    i = i + 1
  }
  // Build universe (all unique elements)
  let universe = build_universe(sets)
  let total_elements = universe.length()
  // For each element, compute its region id (bitmask)
  // region_id is an Int where bit i = 1 if element is in set i
  // Group elements by region
  let n_regions = 1 << n_sets
  let region_elements = Array::make(n_regions, Array::new())
  let mut r = 0
  while r < n_regions {
    region_elements[r] = Array::new()
    r = r + 1
  }
  i = 0
  while i < total_elements {
    let element = universe[i]
    let region_id = compute_region_id(element, sets)
    region_elements[region_id].push(element)
    i = i + 1
  }
  // Build regions
  let regions = Array::make(n_regions, VennRegion::{id: 0, elements: Array::new(), count: 0, description: ""})
  r = 0
  while r < n_regions {
    let count = region_elements[r].length()
    let description = region_description(r, set_names)
    regions[r] = VennRegion::{
      id: r,
      elements: region_elements[r],
      count,
      description,
    }
    r = r + 1
  }
  // Compute Euler characteristic (simplified)
  let euler_char = compute_euler_characteristic(n_sets, regions)
  VennDiagramResult::{
    n_sets,
    set_names,
    set_sizes,
    universe,
    regions,
    total_elements,
    euler_char,
  }
}

///|
/// Build the universe (all unique elements across all sets).
fn build_universe(sets : Array[Array[String]]) -> Array[String] {
  let result = Array::new()
  let n = sets.length()
  let mut i = 0
  while i < n {
    let s = sets[i]
    let mut j = 0
    while j < s.length() {
      let elem = s[j]
      if !array_contains(result, elem) {
        result.push(elem)
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Compute the region ID (bitmask) for an element.
fn compute_region_id(element : String, sets : Array[Array[String]]) -> Int {
  let n = sets.length()
  let mut id = 0
  let mut i = 0
  while i < n {
    if array_contains(sets[i], element) {
      id = id | (1 << i)
    }
    i = i + 1
  }
  id
}

///|
/// Check if an element exists in a string array.
fn array_contains(arr : Array[String], elem : String) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == elem {
      return true
    }
    i = i + 1
  }
  false
}

///|
/// Generate a human-readable description for a region.
fn region_description(region_id : Int, set_names : Array[String]) -> String {
  if region_id == 0 {
    return "Outside all sets"
  }
  let n = set_names.length()
  let parts = Array::new()
  let mut i = 0
  while i < n {
    if ((region_id >> i) & 1) == 1 {
      parts.push(set_names[i])
    }
    i = i + 1
  }
  if parts.length() == 1 {
    "Only in \{parts[0]}"
  } else {
    let mut s = "In "
    let mut j = 0
    while j < parts.length() {
      if j > 0 {
        s = s + " & "
      }
      s = s + parts[j]
      j = j + 1
    }
    s
  }
}

///|
/// Compute Euler characteristic for the Venn diagram.
/// For n sets, chi = V - E + F where:
///   V = number of vertices (intersection points)
///   E = number of edges
///   F = number of faces (regions)
fn compute_euler_characteristic(n_sets : Int, regions : Array[VennRegion]) -> Double {
  // Simplified: just use number of non-empty regions
  let mut non_empty = 0
  let mut i = 0
  while i < regions.length() {
    if regions[i].count > 0 {
      non_empty = non_empty + 1
    }
    i = i + 1
  }
  // For simple Venn diagrams, chi = 1 typically
  if non_empty > 0 { 1.0 } else { 0.0 }
}

///|
/// Get elements in a specific region defined by its bitmask.
/// For example, region 0b101 (5) means elements in set 0 and set 2 but not set 1.
pub fn venn_get_region(result : VennDiagramResult, region_id : Int) -> Array[String] {
  let n = result.regions.length()
  if region_id >= 0 && region_id < n {
    result.regions[region_id].elements
  } else {
    Array::new()
  }
}

///|
/// Get elements only in a specific set (not in any other).
pub fn venn_only_in_set(result : VennDiagramResult, set_index : Int) -> Array[String] {
  let region_id = 1 << set_index
  venn_get_region(result, region_id)
}

///|
/// Get elements in the intersection of specified sets and no others.
/// Sets to include are specified as an array of set indices.
pub fn venn_intersection_only(result : VennDiagramResult, set_indices : Array[Int]) -> Array[String] {
  let mut region_id = 0
  let mut i = 0
  while i < set_indices.length() {
    region_id = region_id | (1 << set_indices[i])
    i = i + 1
  }
  venn_get_region(result, region_id)
}

///|
/// Get elements in the union of all specified sets.
pub fn venn_union(result : VennDiagramResult) -> Array[String] {
  result.universe
}

///|
/// Get elements in the intersection of all sets (regions where all sets overlap).
pub fn venn_all_intersect(result : VennDiagramResult) -> Array[String] {
  let region_id = (1 << result.n_sets) - 1
  venn_get_region(result, region_id)
}

///|
/// Get elements that are unique to each set (only in that set).
pub fn venn_unique_to_each(result : VennDiagramResult) -> Array[Array[String]] {
  let n = result.n_sets
  let result_arr = Array::make(n, Array::new())
  let mut i = 0
  while i < n {
    result_arr[i] = venn_only_in_set(result, i)
    i = i + 1
  }
  result_arr
}

///|
/// Compute pairwise overlap statistics between two Venn sets.
/// Returns overlap count, Jaccard index, and other similarity measures.
pub fn venn_pairwise_overlap(
  result : VennDiagramResult,
  set_i : Int,
  set_j : Int,
) -> Map[String, Double] {
  let elements_i = venn_set_elements(result, set_i)
  let elements_j = venn_set_elements(result, set_j)
  let n_i = elements_i.length()
  let n_j = elements_j.length()
  // Count intersection
  let inter = venn_count_intersection(elements_i, elements_j)
  let union = n_i + n_j - inter
  let jaccard = if union > 0 { inter.to_double() / union.to_double() } else { 0.0 }
  let overlap_coeff = if min_int(n_i, n_j) > 0 {
    inter.to_double() / min_int(n_i, n_j).to_double()
  } else {
    0.0
  }
  let dice = if n_i + n_j > 0 {
    (2.0 * inter.to_double()) / (n_i.to_double() + n_j.to_double())
  } else {
    0.0
  }
  let map = Map::new()
  map["intersection_count"] = inter.to_double()
  map["union_count"] = union.to_double()
  map["jaccard_index"] = jaccard
  map["overlap_coefficient"] = overlap_coeff
  map["dice_coefficient"] = dice
  map["set_i_size"] = n_i.to_double()
  map["set_j_size"] = n_j.to_double()
  map
}

///|
/// Get all elements in a specific set.
fn venn_set_elements(result : VennDiagramResult, set_index : Int) -> Array[String] {
  let elements = Array::new()
  let region_id = 1 << set_index
  // Find all regions that include this set (bit set)
  let mut i = 0
  while i < result.regions.length() {
    if (i & region_id) == region_id {
      // This region is a subset of the set
      let elems = result.regions[i].elements
      let mut j = 0
      while j < elems.length() {
        elements.push(elems[j])
        j = j + 1
      }
    }
    i = i + 1
  }
  elements
}

///|
/// Count intersection of two string arrays.
fn venn_count_intersection(a : Array[String], b : Array[String]) -> Int {
  let mut count = 0
  let mut i = 0
  while i < a.length() {
    if array_contains(b, a[i]) {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Generate a Venn summary report as a string.
pub fn venn_summary(result : VennDiagramResult) -> String {
  let mut s = "=== Venn Diagram Summary ===\n"
  s = s + "Number of sets: \{result.n_sets}\n"
  s = s + "Total elements (universe): \{result.total_elements}\n"
  s = s + "\n"
  s = s + "Set Sizes:\n"
  let mut i = 0
  while i < result.n_sets {
    s = s + "  \{result.set_names[i]}: \{result.set_sizes[i]} elements\n"
    i = i + 1
  }
  s = s + "\n"
  s = s + "Regions:\n"
  let mut r = 0
  while r < result.regions.length() {
    if result.regions[r].count > 0 {
      s = s + "  Region \{r} (0x\{hex_string(r)}): \{result.regions[r].count} elements - \{result.regions[r].description}\n"
    }
    r = r + 1
  }
  s = s + "\n"
  s = s + "Pairwise Overlaps:\n"
  i = 0
  while i < result.n_sets {
    let mut j = i + 1
    while j < result.n_sets {
      if j < result.n_sets {
        let overlap = venn_pairwise_overlap(result, i, j)
        s = s + "  \{result.set_names[i]} & \{result.set_names[j]}: "
        s = s + "intersection=\{overlap["intersection_count"]}, "
        s = s + "Jaccard=\{format_venn_stat(overlap["jaccard_index"])}, "
        s = s + "Overlap=\{format_venn_stat(overlap["overlap_coefficient"])}\n"
      }
      j = j + 1
    }
    i = i + 1
  }
  s
}

///|
/// Compute a Venn diagram with only the counts (no element tracking).
/// More efficient for large datasets where you only need region counts.
pub fn venn_counts_only(sets : Array[Array[String]], names? : Array[String] = []) -> Array[Int] {
  let n_sets = sets.length()
  let n_regions = 1 << n_sets
  let counts = Array::make(n_regions, 0)
  // For each set, build a Map for O(1) lookup
  let set_maps = Array::make(n_sets, Map::new())
  let mut i = 0
  while i < n_sets {
    let m = Map::new()
    let s = sets[i]
    let mut j = 0
    while j < s.length() {
      m[s[j]] = true
      j = j + 1
    }
    set_maps[i] = m
    i = i + 1
  }
  // Build universe
  let universe = build_universe(sets)
  let total = universe.length()
  i = 0
  while i < total {
    let elem = universe[i]
    let mut region_id = 0
    let mut s_idx = 0
    while s_idx < n_sets {
      if set_maps[s_idx].contains(elem) {
        region_id = region_id | (1 << s_idx)
      }
      s_idx = s_idx + 1
    }
    counts[region_id] = counts[region_id] + 1
    i = i + 1
  }
  counts
}

///|
/// Compute Euler diagram layout (regions and their positions).
/// Returns normalized positions for each non-empty region.
pub fn venn_euler_layout(result : VennDiagramResult) -> Array[Map[String, Double]] {
  let n_sets = result.n_sets
  let positions = Array::new()
  let mut r = 0
  while r < result.regions.length() {
    if result.regions[r].count > 0 {
      let pos = compute_region_position(r, n_sets)
      positions.push(pos)
    }
    r = r + 1
  }
  positions
}

///|
/// Compute a approximate 2D position for a Venn region.
fn compute_region_position(region_id : Int, n_sets : Int) -> Map[String, Double] {
  // Place sets in a circle, then compute region positions as weighted average
  let m = Map::new()
  if n_sets == 0 {
    m["x"] = 0.0
    m["y"] = 0.0
    return m
  }
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut count = 0.0
  let mut i = 0
  while i < n_sets {
    let bit = 1 << i
    if (region_id & bit) == bit {
      let angle = 2.0 * 3.141592653589793 * i.to_double() / n_sets.to_double()
      sum_x = sum_x + @math.cos(angle)
      sum_y = sum_y + @math.sin(angle)
      count = count + 1.0
    }
    i = i + 1
  }
  if count > 0.0 {
    m["x"] = sum_x / count
    m["y"] = sum_y / count
  } else {
    m["x"] = 0.0
    m["y"] = 0.0
  }
  m["region_id"] = region_id.to_double()
  m["count"] = region_id.to_double() // placeholder
  m
}

///|
/// Perform a pairwise Venn comparison of two sets only (2-set Venn).
pub fn venn_two_sets(set_a : Array[String], set_b : Array[String], name_a? : String = "A", name_b? : String = "B") -> VennDiagramResult {
  venn_diagram([set_a, set_b], names=[name_a, name_b])
}

///|
/// Compute set intersection of two arrays.
pub fn venn_intersection(a : Array[String], b : Array[String]) -> Array[String] {
  let result = Array::new()
  let mut i = 0
  while i < a.length() {
    if array_contains(b, a[i]) {
      result.push(a[i])
    }
    i = i + 1
  }
  result
}

///|
/// Compute set union of two arrays.
pub fn venn_union_two(a : Array[String], b : Array[String]) -> Array[String] {
  let result = a.copy()
  let mut i = 0
  while i < b.length() {
    if !array_contains(result, b[i]) {
      result.push(b[i])
    }
    i = i + 1
  }
  result
}

///|
/// Compute set difference (elements in a but not in b).
pub fn venn_difference(a : Array[String], b : Array[String]) -> Array[String] {
  let result = Array::new()
  let mut i = 0
  while i < a.length() {
    if !array_contains(b, a[i]) {
      result.push(a[i])
    }
    i = i + 1
  }
  result
}

///|
/// Compute symmetric difference (elements in either a or b but not both).
pub fn venn_symmetric_difference(a : Array[String], b : Array[String]) -> Array[String] {
  let diff_ab = venn_difference(a, b)
  let diff_ba = venn_difference(b, a)
  let result = diff_ab.copy()
  let mut i = 0
  while i < diff_ba.length() {
    result.push(diff_ba[i])
    i = i + 1
  }
  result
}

///|
/// Compute Jaccard index between two sets.
pub fn venn_jaccard(a : Array[String], b : Array[String]) -> Double {
  let inter = venn_intersection(a, b)
  let uni = venn_union_two(a, b)
  if uni.length() > 0 {
    inter.length().to_double() / uni.length().to_double()
  } else {
    0.0
  }
}

///|
/// Compute Sørensen–Dice coefficient between two sets.
pub fn venn_dice(a : Array[String], b : Array[String]) -> Double {
  let inter = venn_intersection(a, b)
  let total = a.length() + b.length()
  if total > 0 {
    (2.0 * inter.length().to_double()) / total.to_double()
  } else {
    0.0
  }
}

///|
/// Compute overlap coefficient between two sets.
pub fn venn_overlap_coefficient(a : Array[String], b : Array[String]) -> Double {
  let inter = venn_intersection(a, b)
  let min_size = if a.length() < b.length() { a.length() } else { b.length() }
  if min_size > 0 {
    inter.length().to_double() / min_size.to_double()
  } else {
    0.0
  }
}

///|
/// Sample data: three gene sets for testing.
pub fn venn_sample_gene_sets() -> Array[Array[String]] {
  [
    ["TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "KRAS", "PTEN", "APC", "RB1", "VHL", "NF1", "MSH2", "MLH1", "MSH6", "PMS2"],
    ["BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2", "BRIP1", "RAD51C", "RAD51D", "BARD1", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG"],
    ["TP53", "EGFR", "MYC", "KRAS", "MTOR", "PIK3CA", "AKT1", "TSC1", "TSC2", "RICTOR", "RAPTOR", "EIF4E", "4EBP1", "S6K1", "PDK1"],
  ]
}

///|
/// Sample gene set names.
pub fn venn_sample_gene_names() -> Array[String] {
  ["Tumor Suppressors", "DNA Repair", "PI3K/MTOR"]
}

///|
/// Format an integer as a hexadecimal string.
fn hex_string(value : Int) -> String {
  if value == 0 {
    return "0"
  }
  let digits = "0123456789ABCDEF"
  let mut result = ""
  let mut v = value
  while v > 0 {
    let d = v % 16
    result = "\{digits[d]}" + result
    v = v / 16
  }
  result
}

///|
/// Format a double for Venn display.
fn format_venn_stat(value : Double) -> String {
  let rounded = (value * 1000.0).round() / 1000.0
  "\{rounded}"
}

///|
/// Create a Venn diagram result with specific set names for testing.
pub fn venn_sample_result() -> VennDiagramResult {
  let sets = venn_sample_gene_sets()
  let names = venn_sample_gene_names()
  venn_diagram(sets, names=names)
}
