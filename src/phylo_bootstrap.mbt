///|
/// Consensus trees (strict, majority-rule, Adam compatibility consensus)
/// and nonparametric bootstrap resampling, ported from Biopython
/// ``Bio.Phylo.Consensus``.  All functions operate on the native
/// [Tree] / [Clade] types defined in this package.
///
/// Following Biopython, clades are identified across trees by the set of
/// their terminal names, encoded as fixed-width bit strings over the terminal
/// order of the first tree (``get_terminals`` preorder).

///|
/// Error raised when consensus or bootstrap input is invalid.
pub suberror PhyloConsensusError {
  PhyloConsensusError(String)
}

///|
/// Fixed-width binary string over terminal taxa.
priv struct BStr {
  bits : Array[Bool]
}

///|
fn bstr_ones(n : Int) -> BStr {
  BStr::{ bits: Array::make(n, true), }
}

///|
/// Serialize to a Biopython-style ``"01101"`` key.
fn bstr_key(bs : BStr) -> String {
  let sb = StringBuilder(size_hint=bs.bits.length())
  for b in bs.bits {
    sb.write_char(if b { '1' } else { '0' })
  }
  sb.to_string()
}

///|
/// Rebuild a bit string from a ``"0110"`` key.
fn bstr_from_key(key : String) -> BStr {
  let bits = Array::make(key.length(), false)
  for i in 0..<key.length() {
    bits[i] = key.unsafe_get(i) == '1'.to_int().to_uint16()
  }
  BStr::{ bits, }
}

///|
fn bstr_and(a : BStr, b : BStr) -> BStr {
  let n = a.bits.length()
  let bits = Array::make(n, false)
  for i in 0..<n {
    bits[i] = a.bits[i] && b.bits[i]
  }
  BStr::{ bits, }
}

///|
fn bstr_xor(a : BStr, b : BStr) -> BStr {
  let n = a.bits.length()
  let bits = Array::make(n, false)
  for i in 0..<n {
    bits[i] = a.bits[i] != b.bits[i]
  }
  BStr::{ bits, }
}

///|
fn bstr_count_ones(bs : BStr) -> Int {
  let mut n = 0
  for b in bs.bits {
    if b {
      n += 1
    }
  }
  n
}

///|
fn bstr_index_one(bs : BStr) -> Array[Int] {
  let idx = Array(capacity=0)
  for i in 0..<bs.bits.length() {
    if bs.bits[i] {
      idx.push(i)
    }
  }
  idx
}

///|
/// Does ``a`` contain ``b`` (ones of ``b`` subset of ones of ``a``)?
fn bstr_contains(a : BStr, b : BStr) -> Bool {
  for i in 0..<a.bits.length() {
    if b.bits[i] && !a.bits[i] {
      return false
    }
  }
  true
}

///|
/// No shared ones.
fn bstr_independent(a : BStr, b : BStr) -> Bool {
  for i in 0..<a.bits.length() {
    if a.bits[i] && b.bits[i] {
      return false
    }
  }
  true
}

///|
/// Containment either way, or independence (mirrors ``_BitString``).
fn bstr_iscompatible(a : BStr, b : BStr) -> Bool {
  bstr_contains(a, b) || bstr_contains(b, a) || bstr_independent(a, b)
}

///|
/// Bit string for a clade given the ordered terminal names.
fn clade_to_bstr(
  clade : Clade,
  term_names : Array[String],
) -> BStr raise PhyloConsensusError {
  let name_set : Map[String, Unit] = Map::new()
  for t in clade.get_terminals() {
    match t.name {
      Some(n) => name_set.set(n, ())
      Option::None =>
        raise PhyloConsensusError("Clade contains a terminal without a name")
    }
  }
  let bits = Array::make(term_names.length(), false)
  for i in 0..<term_names.length() {
    bits[i] = name_set.contains(term_names[i])
  }
  BStr::{ bits, }
}

///|
priv struct CladeCount {
  count : Int
  bl_sum : Double
}

///|
/// Count distinct non-terminal clades across trees (PRIVATE counterpart of
/// Biopython ``_count_clades``).  Returns terminal names (preorder of the
/// first tree), the first tree's terminal clades, and the count table.
fn count_clades(
  trees : Array[Tree],
) -> (Array[String], Array[Clade], Map[String, CladeCount]) raise PhyloConsensusError {
  if trees.is_empty() {
    raise PhyloConsensusError("No trees provided")
  }
  let terms = trees[0].get_terminals()
  let term_names = Array(capacity=0)
  for t in terms {
    match t.name {
      Some(n) => term_names.push(n)
      Option::None =>
        raise PhyloConsensusError("Clade contains a terminal without a name")
    }
  }
  let expected = term_names.copy()
  expected.sort()
  let counts : Map[String, CladeCount] = Map::new()
  for tree in trees {
    let these_names = Array(capacity=0)
    for t in tree.get_terminals() {
      match t.name {
        Some(n) => these_names.push(n)
        Option::None =>
          raise PhyloConsensusError("Clade contains a terminal without a name")
      }
    }
    these_names.sort()
    if these_names.length() != expected.length() {
      raise PhyloConsensusError("Taxons in provided trees should be consistent")
    }
    for i in 0..<expected.length() {
      if these_names[i] != expected[i] {
        raise PhyloConsensusError(
          "Taxons in provided trees should be consistent",
        )
      }
    }
    for clade in tree.find_clades(terminal=Some(false)) {
      let bs = clade_to_bstr(clade, term_names)
      let key = bstr_key(bs)
      let bl = match clade.branch_length {
        Some(v) => v
        Option::None => 0.0
      }
      match counts.get(key) {
        Some(cc) =>
          counts.set(key, CladeCount::{
            count: cc.count + 1,
            bl_sum: cc.bl_sum + bl,
          })
        Option::None => counts.set(key, CladeCount::{ count: 1, bl_sum: bl, })
      }
    }
  }
  (term_names, terms, counts)
}

///|
/// Sorted terminal-name signature of a clade, used as its identity key
/// (clades are considered equal when their terminal sets agree).
fn clade_tip_key(c : Clade) -> String {
  let names = Array(capacity=0)
  for t in c.get_terminals() {
    match t.name {
      Some(n) => names.push(n)
      Option::None => abort("clade_tip_key: unnamed terminal")
    }
  }
  names.sort()
  let sb = StringBuilder()
  for i in 0..<names.length() {
    if i > 0 {
      sb.write_char(',')
    }
    sb.write_string(names[i])
  }
  sb.to_string()
}

///|
/// Same clade within a tree family: identical set of terminal names.
fn same_clade(a : Clade, b : Clade) -> Bool {
  clade_tip_key(a) == clade_tip_key(b)
}

///|
/// Replace the children of a clade in place (the field itself is immutable).
fn set_children(node : Clade, children : Array[Clade]) -> Unit {
  node.clades.clear()
  for c in children {
    node.clades.push(c)
  }
}

///|
/// Strict consensus tree: a clade appears only if present in every tree.
pub fn consensus_strict(trees : Array[Tree]) -> Tree raise PhyloConsensusError {
  let (_, terms, counts) = count_clades(trees)
  let n = terms.length()
  let root_key = bstr_key(bstr_ones(n))

  // strict clades, sorted by decreasing number of terminals
  let strict_keys = Array(capacity=0)
  for kv in counts.iter() {
    if kv.1.count == trees.length() {
      strict_keys.push(kv.0)
    }
  }
  strict_keys.sort_by(fn(a, b) {
    bstr_count_ones(bstr_from_key(b)).compare(bstr_count_ones(bstr_from_key(a)))
  })

  let root = Clade::new()
  if strict_keys.is_empty() || strict_keys[0] != root_key {
    raise PhyloConsensusError("Taxons in provided trees should be consistent")
  }
  for t in terms {
    root.clades.push(t)
  }
  let bitstr_clades : Map[String, Clade] = Map::new()
  bitstr_clades.set(root_key, root)

  for key in strict_keys[1:] {
    let bs = bstr_from_key(key)
    let clade_terms = Array(capacity=0)
    for i in bstr_index_one(bs) {
      clade_terms.push(terms[i])
    }
    let inner = Clade::new()
    for t in clade_terms {
      inner.clades.push(t)
    }
    // find the closest existing clade that contains this one
    let mut found_parent = false
    let entries = bitstr_clades.iter().to_array()
    for pair in entries {
      let parent_key = pair.0
      let parent_clade = pair.1
      if bstr_contains(bstr_from_key(parent_key), bs) {
        bitstr_clades.remove(parent_key)
        let kept = parent_clade.clades.filter(fn(c) {
          !clade_terms.any(fn(t) { same_clade(c, t) })
        })
        kept.push(inner)
        set_children(parent_clade, kept)
        let new_key = bstr_key(bstr_xor(bstr_from_key(parent_key), bs))
        bitstr_clades.set(new_key, parent_clade)
        found_parent = true
        break
      }
    }
    if !found_parent {
      raise PhyloConsensusError("Inconsistent strict consensus state")
    }
    bitstr_clades.set(key, inner)
  }
  Tree::new(root)
}

///|
priv struct MajorityEntry {
  key : String
  count : Int
  bl_sum : Double
}

///|
/// Majority-rule consensus tree with an extended cutoff in [0, 1].
///
/// Clades occurring in at least ``cutoff`` fraction of the trees are kept;
/// each consensus clade is annotated with its percent confidence and the
/// average branch length over its occurrences (as in Biopython).
pub fn consensus_majority(
  trees : Array[Tree],
  cutoff? : Double = 0.0,
) -> Tree raise PhyloConsensusError {
  let (_, terms, counts) = count_clades(trees)
  let n = terms.length()
  let tree_count = trees.length().to_double()

  let entries = Array(capacity=0)
  for kv in counts.iter() {
    let cc = kv.1
    entries.push(MajorityEntry::{
      key: kv.0,
      count: cc.count,
      bl_sum: cc.bl_sum,
    })
  }
  // descending (count, #ones, key)
  entries.sort_by(fn(a, b) {
    if a.count != b.count {
      b.count.compare(a.count)
    } else {
      let oa = bstr_count_ones(bstr_from_key(a.key))
      let ob = bstr_count_ones(bstr_from_key(b.key))
      if oa != ob {
        ob.compare(oa)
      } else {
        b.key.compare(a.key)
      }
    }
  })

  let root = Clade::new()
  if entries.is_empty() || bstr_count_ones(bstr_from_key(entries[0].key)) != n {
    raise PhyloConsensusError("Taxons in provided trees should be consistent")
  }
  for t in terms {
    root.clades.push(t)
  }
  let bitstr_clades : Map[String, Clade] = Map::new()
  bitstr_clades.set(entries[0].key, root)

  for entry in entries[1:] {
    let bitstr = bstr_from_key(entry.key)
    let confidence = 100.0 * entry.count.to_double() / tree_count
    if confidence < cutoff * 100.0 {
      break
    }
    let clade_terms = Array(capacity=0)
    for i in bstr_index_one(bitstr) {
      clade_terms.push(terms[i])
    }
    let inner = Clade::new(
      branch_length=Some(entry.bl_sum / entry.count.to_double()),
      confidence=Some(confidence),
    )
    for t in clade_terms {
      inner.clades.push(t)
    }

    let bsckeys = Array(capacity=0)
    for kv in bitstr_clades.iter() {
      bsckeys.push(kv.0)
    }
    bsckeys.sort_by(fn(a, b) {
      bstr_count_ones(bstr_from_key(b)).compare(
        bstr_count_ones(bstr_from_key(a)),
      )
    })

    let mut compatible = true
    let mut parent_key : String? = Option::None
    let child_keys = Array(capacity=0)
    for k in bsckeys {
      let bs = bstr_from_key(k)
      if !bstr_iscompatible(bs, bitstr) {
        compatible = false
        break
      }
      if bstr_contains(bs, bitstr) {
        parent_key = Some(k)
      }
      if bstr_contains(bitstr, bs) && k != entry.key {
        let mut independent_of_children = true
        for ck in child_keys {
          if !bstr_independent(bstr_from_key(ck), bs) {
            independent_of_children = false
            break
          }
        }
        if independent_of_children {
          child_keys.push(k)
        }
      }
    }
    if !compatible {
      continue
    }

    let mut parent_clade : Clade? = Option::None
    match parent_key {
      Some(pk) => {
        let pc = bitstr_clades.get(pk).unwrap()
        bitstr_clades.remove(pk)
        let kept = pc.clades.filter(fn(c) {
          !clade_terms.any(fn(t) { same_clade(c, t) })
        })
        kept.push(inner)
        set_children(pc, kept)
        bitstr_clades.set(pk, pc)
        parent_clade = Some(pc)
      }
      Option::None => ()
    }

    if !child_keys.is_empty() {
      let remove_idx = Array(capacity=0)
      for ck in child_keys {
        for i in bstr_index_one(bstr_from_key(ck)) {
          remove_idx.push(i)
        }
        let child_clade = bitstr_clades.get(ck).unwrap()
        let pc = parent_clade.unwrap()
        remove_first_identical(pc.clades, child_clade)
        inner.clades.push(child_clade)
      }
      let remove_terms = Array(capacity=0)
      for i in remove_idx {
        remove_terms.push(terms[i])
      }
      let kept_children = inner.clades.filter(fn(c) {
        !remove_terms.any(fn(t) { same_clade(c, t) })
      })
      set_children(inner, kept_children)
    }

    bitstr_clades.set(entry.key, inner)
    let map_size = bitstr_clades.length()
    if map_size == n - 1 || (map_size == n - 2 && root.clades.length() == 3) {
      break
    }
  }
  Tree::new(root)
}

///|
fn remove_first_identical(arr : Array[Clade], target : Clade) -> Unit {
  for i in 0..<arr.length() {
    if same_clade(arr[i], target) {
      let _ = arr.remove(i)
      return
    }
  }
}

///|
/// Adam (compatibility) consensus tree, rooted.
pub fn consensus_adam(trees : Array[Tree]) -> Tree raise PhyloConsensusError {
  if trees.is_empty() {
    raise PhyloConsensusError("No trees provided")
  }
  let roots = Array(capacity=0)
  for t in trees {
    roots.push(t.get_root())
  }
  Tree::new(adam_part(roots), rooted=true)
}

///|
/// Recursive partition step of the Adam algorithm (counterpart of ``_part``).
fn adam_part(clades : Array[Clade]) -> Clade raise PhyloConsensusError {
  let terms = clades[0].get_terminals()
  if terms.length() == 1 || terms.length() == 2 {
    return clades[0]
  }
  let term_names = Array(capacity=0)
  for t in terms {
    match t.name {
      Some(n) => term_names.push(n)
      Option::None =>
        raise PhyloConsensusError("Clade contains a terminal without a name")
    }
  }
  let n = terms.length()
  let mut parts : Array[BStr] = [bstr_ones(n)]
  for clade in clades {
    for child in clade.clades {
      let bitstr = clade_to_bstr(child, term_names)
      let to_remove = Array(capacity=0)
      let to_add = Array(capacity=0)
      for bs in parts.copy() {
        let bkey = bstr_key(bs)
        let ikey = bstr_key(bitstr)
        if bkey == ikey {
          continue
        } else if bstr_contains(bs, bitstr) {
          to_add.push(bitstr)
          to_add.push(bstr_xor(bs, bitstr))
          to_remove.push(bkey)
        } else if bstr_contains(bitstr, bs) {
          to_add.push(bstr_xor(bs, bitstr))
        } else if !bstr_independent(bs, bitstr) {
          let inter = bstr_and(bs, bitstr)
          to_add.push(inter)
          to_add.push(bstr_xor(inter, bitstr))
          to_add.push(bstr_xor(inter, bs))
          to_remove.push(bkey)
        }
      }
      for rk in to_remove {
        parts = parts.filter(fn(bs) { bstr_key(bs) != rk })
      }
      if !to_add.is_empty() {
        let sorted_adds = to_add.copy()
        sorted_adds.sort_by(fn(a, b) {
          bstr_count_ones(a).compare(bstr_count_ones(b))
        })
        for ta in sorted_adds {
          let mut independent = true
          for bs in parts {
            if !bstr_independent(ta, bs) {
              independent = false
              break
            }
          }
          if independent {
            parts.push(ta)
          }
        }
      }
    }
  }
  let new_clade = Clade::new()
  let sorted_parts = parts.copy()
  sorted_parts.sort_by(fn(a, b) { bstr_key(a).compare(bstr_key(b)) })
  for bs in sorted_parts {
    let indices = bstr_index_one(bs)
    if indices.length() == 1 {
      new_clade.clades.push(terms[indices[0]])
    } else if indices.length() == 2 {
      let bifur = Clade::new()
      bifur.clades.push(terms[indices[0]])
      bifur.clades.push(terms[indices[1]])
      new_clade.clades.push(bifur)
    } else if indices.length() > 2 {
      let part_names = Array(capacity=0)
      for i in indices {
        part_names.push(term_names[i])
      }
      let next_clades = Array(capacity=0)
      for clade in clades {
        next_clades.push(adam_sub_clade(clade, part_names))
      }
      new_clade.clades.push(adam_part(next_clades))
    }
  }
  new_clade
}

///|
/// Extract a compatible subclade containing exactly the given terminals
/// (counterpart of ``_sub_clade``).
fn adam_sub_clade(
  clade : Clade,
  term_names : Array[String],
) -> Clade raise PhyloConsensusError {
  let term_clades = Array(capacity=0)
  for n in term_names {
    match clade.find_any(name=Some(n)) {
      Some(tc) => term_clades.push(tc)
      Option::None =>
        raise PhyloConsensusError(
          "Taxons in provided trees should be consistent",
        )
    }
  }
  let mut sub_clade = clade.clade_common_ancestor(term_clades)
  if term_names.length() != sub_clade.count_terminals() {
    let temp_clade = Clade::new()
    for tc in term_clades {
      temp_clade.clades.push(tc)
    }
    let top = clade_top(sub_clade)
    for c in sub_clade.find_clades(terminal=Some(false)) {
      if same_clade(c, top) {
        continue
      }
      let children = Array(capacity=0)
      for t in c.get_terminals() {
        if term_clades.any(fn(x) { same_clade(t, x) }) {
          children.push(t)
        }
      }
      if children.is_empty() {
        continue
      }
      for tc in temp_clade.find_clades(terminal=Some(false)) {
        let tc_new = tc.clades.filter(fn(x) {
          !children.any(fn(ch) { same_clade(x, ch) })
        })
        let all_present = children.all(fn(ch) {
          tc.clades.any(fn(x) { same_clade(ch, x) })
        })
        if all_present && !tc_new.is_empty() {
          set_children(tc, tc_new)
          let child_clade = Clade::new()
          for ch in children {
            child_clade.clades.push(ch)
          }
          tc.clades.push(child_clade)
        }
      }
    }
    sub_clade = temp_clade
  }
  sub_clade
}

///|
/// Most recent common ancestor of terminal clades within ``self``.
fn Clade::clade_common_ancestor(
  self : Clade,
  term_clades : Array[Clade],
) -> Clade {
  let mut common : Array[Clade] = Array(capacity=0)
  for i in 0..<term_clades.length() {
    let t = term_clades[i]
    let name = match t.name {
      Some(n) => n
      Option::None => abort("clade_common_ancestor: unnamed terminal")
    }
    let path = self.get_path(name).unwrap()
    if i == 0 {
      common = path
    } else {
      let mut depth = 0
      while depth < common.length() &&
            depth < path.length() &&
            same_clade(common[depth], path[depth]) {
        depth += 1
      }
      let narrowed = Array(capacity=0)
      for k in 0..<depth {
        narrowed.push(common[k])
      }
      common = narrowed
    }
  }
  common.last().unwrap()
}

///|
/// Walk parent links to the topmost clade.
fn clade_top(start : Clade) -> Clade {
  let mut node = start
  let mut done = false
  while !done {
    match node.parent {
      Some(p) => node = p
      Option::None => done = true
    }
  }
  node
}

///|
/// Annotate each internal clade of ``target`` with its bootstrap support
/// (percent of replicate trees containing the same clade).  Terminal names are
/// ordered alphabetically, as in Biopython ``get_support``.  Clade objects are
/// immutable, so a new annotated tree is returned (terminals are shared).
pub fn tree_get_support(
  target : Tree,
  trees : Array[Tree],
  len_trees? : Int? = None,
) -> Tree raise PhyloConsensusError {
  let term_names = Array(capacity=0)
  for t in target.find_clades(terminal=Some(true)) {
    match t.name {
      Some(n) => term_names.push(n)
      Option::None => ()
    }
  }
  term_names.sort()
  let target_keys : Map[String, Unit] = Map::new()
  for c in target.find_clades(terminal=Some(false)) {
    target_keys.set(bstr_key(clade_to_bstr(c, term_names)), ())
  }
  let seen : Map[String, Int] = Map::new()
  for k in target_keys.keys() {
    seen.set(k, 0)
  }
  let size = match len_trees {
    Some(n) => n
    Option::None => trees.length()
  }
  for tree in trees {
    for c in tree.find_clades(terminal=Some(false)) {
      let key = bstr_key(clade_to_bstr(c, term_names))
      if target_keys.contains(key) {
        seen.set(key, seen.get(key).unwrap() + 1)
      }
    }
  }
  let support : Map[String, Double] = Map::new()
  for kv in seen.iter() {
    support.set(kv.0, kv.1.to_double() * 100.0 / size.to_double())
  }
  let new_root = support_apply_clade(target.get_root(), term_names, support)
  Tree::new(
    new_root,
    rooted=target.rooted,
    name=target.name,
    weight=target.weight,
  )
}

///|
fn support_apply_clade(
  c : Clade,
  term_names : Array[String],
  support : Map[String, Double],
) -> Clade raise PhyloConsensusError {
  if c.is_terminal() {
    return c
  }
  let new_children = Array(capacity=0)
  for child in c.clades {
    new_children.push(support_apply_clade(child, term_names, support))
  }
  let key = bstr_key(clade_to_bstr(c, term_names))
  let confidence = match support.get(key) {
    Some(v) => Some(v)
    Option::None => c.confidence
  }
  Clade::new(
    branch_length=c.branch_length,
    name=c.name,
    confidence~,
    clades=new_children,
  )
}

// ---------------------------------------------------------------------------
// Bootstrap resampling
// ---------------------------------------------------------------------------

///|
/// Park-Miller minimal standard RNG, so column streams are reproducible.
priv struct ParkMiller {
  mut state : Int
}

///|
fn park_miller_next(rng : ParkMiller) -> Int {
  let hi = rng.state / 127773
  let lo = rng.state % 127773
  let mut s = 16807 * lo - 2836 * hi
  if s <= 0 {
    s += 2147483647
  }
  rng.state = s
  s
}

///|
/// Deterministically draw ``times`` bootstrap replicates of ``length``
/// alignment columns (sampling with replacement via a Park-Miller minimal
/// standard PRNG).  Returned column-index arrays can be fed to
/// [bootstrap_trees]; this keeps the random stream explicit and reproducible.
pub fn bootstrap_column_indices(
  length : Int,
  times : Int,
  seed? : Int = 42,
) -> Array[Array[Int]] raise PhyloConsensusError {
  if length <= 0 {
    raise PhyloConsensusError("Alignment length must be positive")
  }
  if times <= 0 {
    raise PhyloConsensusError("Bootstrap times must be positive")
  }
  let seed_state = if seed == 0 { 1 } else if seed < 0 { -seed } else { seed }
  let rng = ParkMiller::{ state: seed_state, }
  let result = Array(capacity=0)
  for _ in 0..<times {
    let cols = Array(capacity=0)
    for _ in 0..<length {
      cols.push((park_miller_next(rng) - 1) % length)
    }
    result.push(cols)
  }
  result
}

///|
/// Build one bootstrap replicate tree per column-index array by resampling
/// alignment columns, recomputing pairwise distances (identity model by
/// default) and reconstructing with the neighbor-joining constructor.
pub fn bootstrap_trees(
  sequences : Array[String],
  names : Array[String],
  columns : Array[Array[Int]],
  model? : String = "identity",
  tree_method? : String = "nj",
) -> Array[Tree] raise PhyloConsensusError {
  validate_alignment(sequences, names)
  if columns.is_empty() {
    raise PhyloConsensusError("No bootstrap columns provided")
  }
  let length = sequences[0].length()
  for cols in columns {
    if cols.length() != length {
      raise PhyloConsensusError(
        "Each bootstrap replicate must sample as many columns as the alignment",
      )
    }
    for j in cols {
      if j < 0 || j >= length {
        raise PhyloConsensusError("Bootstrap column index out of range")
      }
    }
  }
  let calc = DistanceCalculator::new(model~) catch {
    SeqError(msg) =>
      raise PhyloConsensusError("Unsupported distance model: " + msg)
  }
  let tree_constructor = TreeConstructor::new(tree_method~)
  let result = Array(capacity=0)
  for cols in columns {
    let replicate = Array(capacity=0)
    for seq in sequences {
      replicate.push(resample_sequence(seq, cols))
    }
    let ntax = names.length()
    let matrix = Array(capacity=0)
    for _ in 0..<ntax {
      matrix.push(Array::make(ntax, 0.0))
    }
    for i in 0..<ntax {
      for j in (i + 1)..<ntax {
        let d = calc.pairwise(replicate[i], replicate[j]) catch {
          SeqError(msg) =>
            raise PhyloConsensusError("Distance calculation failed: " + msg)
        }
        matrix[i][j] = d
        matrix[j][i] = d
      }
    }
    let root = tree_constructor.build_tree(matrix, names)
    result.push(Tree::new(root))
  }
  result
}

///|
/// Build bootstrap replicate trees and immediately summarize them with a
/// consensus method (``"strict"``, ``"majority"`` or ``"adam"``).
pub fn bootstrap_consensus(
  sequences : Array[String],
  names : Array[String],
  columns : Array[Array[Int]],
  consensus? : String = "majority",
  cutoff? : Double = 0.0,
  model? : String = "identity",
  tree_method? : String = "nj",
) -> Tree raise PhyloConsensusError {
  let trees = bootstrap_trees(sequences, names, columns, model~, tree_method~)
  match consensus {
    "strict" => consensus_strict(trees)
    "majority" => consensus_majority(trees, cutoff~)
    "adam" => consensus_adam(trees)
    other =>
      raise PhyloConsensusError(
        "Unknown consensus method: " +
        other +
        " (expected strict, majority or adam)",
      )
  }
}

///|
fn validate_alignment(
  sequences : Array[String],
  names : Array[String],
) -> Unit raise PhyloConsensusError {
  if sequences.is_empty() {
    raise PhyloConsensusError("No sequences provided")
  }
  if sequences.length() != names.length() {
    raise PhyloConsensusError("Number of sequences and names must match")
  }
  let length = sequences[0].length()
  if length == 0 {
    raise PhyloConsensusError("Sequences must not be empty")
  }
  for s in sequences {
    if s.length() != length {
      raise PhyloConsensusError("All sequences must have equal length")
    }
  }
}

///|
fn resample_sequence(seq : String, columns : Array[Int]) -> String {
  let sb = StringBuilder(size_hint=columns.length())
  for j in columns {
    sb.write_char(Char::from_int(seq.unsafe_get(j).to_int()))
  }
  sb.to_string()
}

// ---------------------------------------------------------------------------
// Newick serialization (supports confidence labels and branch lengths)
// ---------------------------------------------------------------------------

///|
/// Serialize a clade and its descendants to Newick, writing branch lengths
/// when present and using ``confidence`` (or the name) as the internal label.
pub fn Clade::to_newick(self : Clade, lengths? : Bool = true) -> String {
  let sb = StringBuilder()
  write_consensus_newick(self, lengths, sb)
  sb.write_char(';')
  sb.to_string()
}

///|
fn write_consensus_newick(
  clade : Clade,
  lengths : Bool,
  sb : StringBuilder,
) -> Unit {
  if !clade.is_terminal() {
    sb.write_char('(')
    for i in 0..<clade.clades.length() {
      if i > 0 {
        sb.write_char(',')
      }
      write_consensus_newick(clade.clades[i], lengths, sb)
    }
    sb.write_char(')')
  }
  match clade.name {
    Some(n) => sb.write_string(n)
    Option::None =>
      match clade.confidence {
        Some(cf) => sb.write_string(fmt_fixed6(cf))
        Option::None => ()
      }
  }
  if lengths {
    match clade.branch_length {
      Some(bl) => {
        sb.write_char(':')
        sb.write_string(fmt_fixed6(bl))
      }
      Option::None => ()
    }
  }
}

///|
/// Format a double with six decimals, trailing zeros trimmed.
fn fmt_fixed6(v : Double) -> String {
  let neg = v < 0.0
  let av = if neg { -v } else { v }
  let scaled = (av * 1000000.0 + 0.5).to_int()
  let neg_prefix = if neg && scaled > 0 { "-" } else { "" }
  let int_part = (scaled / 1000000).to_string()
  let mut frac = (scaled % 1000000).to_string()
  while frac.length() < 6 {
    frac = "0" + frac
  }
  let mut end = frac.length()
  while end > 0 && frac.get_char(end - 1).unwrap() == '0' {
    end -= 1
  }
  if end == 0 {
    neg_prefix + int_part
  } else {
    neg_prefix + int_part + "." + frac.substring(end~)
  }
}
