///|
/// scirpy-style immune-receptor (IR) object model built on top of the
/// `vdj`/`airr` chain calls.  Chains are paired per cell into the scirpy
/// slot layout (`primary_vj` / `secondary_vj` for TRA/TRG and `primary_vdj` /
/// `secondary_vdj` for TRB/TRD/BCR heavy chains), cells are assigned to
/// clonotypes by the (v, j, cdr3) key of their primary chains, and a clonotype
/// network (`scanpy_ir_clonotype_network`) links clonotypes that share an
/// exact CDR3 or have close CDR3 sequences.

///|
/// One immune-receptor chain call (a minimal AIRR record).
pub(all) struct ScanpyIrChain {
  cell_id : String
  /// locus: "TRA"/"TRG" are VJ chains, everything else VDJ
  locus : String
  v_call : String
  j_call : String
  cdr3_aa : String
  productive : Bool
}

///|
/// A cell with its paired chains in scirpy slot semantics.
pub struct ScanpyIrCell {
  cell_id : String
  primary_vj : ScanpyIrChain?
  secondary_vj : ScanpyIrChain?
  primary_vdj : ScanpyIrChain?
  secondary_vdj : ScanpyIrChain?
  /// more than two chains of one kind (multi-chain cell, excluded from
  /// clonotype assignment like scirpy)
  multichain : Bool
  clonotype_id : Int32
}

///|
/// A clonotype: the set of cells sharing the primary-chain keys.
pub struct ScanpyIrClonotype {
  clonotype_id : Int32
  /// "v:j:cdr3" of the primary VJ chain, or "" when absent
  vj_key : String
  /// "v:j:cdr3" of the primary VDJ chain, or "" when absent
  vdj_key : String
  cell_ids : Array[String]
}

///|
/// Paired IR data: cells with chain slots plus the clonotype table.
pub struct ScanpyIrData {
  cells : Array[ScanpyIrCell]
  clonotypes : Array[ScanpyIrClonotype]
}

///|
/// True when the locus is a VJ chain ("TRA" / "TRG").
pub fn scanpy_ir_is_vj(locus : String) -> Bool {
  locus == "TRA" || locus == "TRG"
}

///|
/// Pair chains into cells and assign clonotypes.  Within each cell the VJ
/// and VDJ chains are ranked productive-first (input order breaks ties);
/// the top chain fills the primary slot, the second the secondary slot and
/// any further chain marks the cell `multichain` (such cells keep
/// clonotype_id -1, scirpy's convention for ambiguous calls).  Clonotypes
/// are the distinct (primary_vj, primary_vdj) key pairs, assigned in
/// first-appearance order.
pub fn scanpy_ir_pair_chains(chains : Array[ScanpyIrChain]) -> ScanpyIrData {
  // group chains by cell, preserving input order
  let order : Array[String] = Array(capacity=0)
  let by_cell : Map[String, (Array[ScanpyIrChain], Array[ScanpyIrChain])] = Map(
    [],
    capacity=16,
  )
  for ch in chains {
    match by_cell.get(ch.cell_id) {
      None => {
        order.push(ch.cell_id)
        let vj : Array[ScanpyIrChain] = Array(capacity=0)
        let vdj : Array[ScanpyIrChain] = Array(capacity=0)
        if scanpy_ir_is_vj(ch.locus) {
          vj.push(ch)
        } else {
          vdj.push(ch)
        }
        by_cell.set(ch.cell_id, (vj, vdj))
      }
      Some((vj, vdj)) =>
        if scanpy_ir_is_vj(ch.locus) {
          vj.push(ch)
        } else {
          vdj.push(ch)
        }
    }
  }
  // chain ranking: productive first, stable otherwise
  fn rank_sort(list : Array[ScanpyIrChain]) -> Array[ScanpyIrChain] {
    let out : Array[ScanpyIrChain] = Array(capacity=list.length())
    for ch in list {
      if ch.productive {
        out.push(ch)
      }
    }
    for ch in list {
      if !ch.productive {
        out.push(ch)
      }
    }
    out
  }
  let cells : Array[ScanpyIrCell] = Array(capacity=order.length())
  let clono_index : Map[String, Int32] = Map([], capacity=16)
  let clonotypes : Array[ScanpyIrClonotype] = Array(capacity=0)
  for cell_id in order {
    let (vj, vdj) = by_cell.get(cell_id).unwrap()
    // VJ chains pair within a single locus (TRA preferred over TRG), like
    // scirpy: a TRA and a TRG chain in the same cell are not cross-paired.
    let tra : Array[ScanpyIrChain] = Array(capacity=0)
    let trg : Array[ScanpyIrChain] = Array(capacity=0)
    for ch in vj {
      if ch.locus == "TRA" {
        tra.push(ch)
      } else {
        trg.push(ch)
      }
    }
    let vj_sorted = if tra.length() > 0 {
      rank_sort(tra)
    } else {
      rank_sort(trg)
    }
    let vdj_sorted = rank_sort(vdj)
    let multichain = tra.length() > 2 ||
      trg.length() > 2 ||
      vdj_sorted.length() > 2
    let primary_vj = if vj_sorted.length() > 0 {
      Some(vj_sorted[0])
    } else {
      None
    }
    let secondary_vj = if vj_sorted.length() > 1 {
      Some(vj_sorted[1])
    } else {
      None
    }
    let primary_vdj = if vdj_sorted.length() > 0 {
      Some(vdj_sorted[0])
    } else {
      None
    }
    let secondary_vdj = if vdj_sorted.length() > 1 {
      Some(vdj_sorted[1])
    } else {
      None
    }
    let mut clonotype_id = -1
    if !multichain {
      let vj_key = match primary_vj {
        Some(c) => c.v_call + ":" + c.j_call + ":" + c.cdr3_aa
        None => ""
      }
      let vdj_key = match primary_vdj {
        Some(c) => c.v_call + ":" + c.j_call + ":" + c.cdr3_aa
        None => ""
      }
      let key = vj_key + "||" + vdj_key
      match clono_index.get(key) {
        None => {
          clonotype_id = clonotypes.length()
          clono_index.set(key, clonotype_id)
          clonotypes.push(ScanpyIrClonotype::{
            clonotype_id,
            vj_key,
            vdj_key,
            cell_ids: [cell_id],
          })
        }
        Some(cid) => {
          clonotype_id = cid
          clonotypes[cid].cell_ids.push(cell_id)
        }
      }
    }
    cells.push(ScanpyIrCell::{
      cell_id,
      primary_vj,
      secondary_vj,
      primary_vdj,
      secondary_vdj,
      multichain,
      clonotype_id,
    })
  }
  ScanpyIrData::{ cells, clonotypes, }
}

///|
/// Number of cells per clonotype, in clonotype-id order.
pub fn scanpy_ir_clonotype_sizes(data : ScanpyIrData) -> Array[(Int32, Int32)] {
  let out : Array[(Int32, Int32)] = Array(capacity=data.clonotypes.length())
  for ct in data.clonotypes {
    out.push((ct.clonotype_id, ct.cell_ids.length()))
  }
  out
}

///|
/// Hamming distance between two equal-length strings; `None` when the lengths
/// differ.
fn scanpy_ir_hamming(a : String, b : String) -> Int32? {
  let n = a.length()
  if n != b.length() {
    return None
  }
  let mut d = 0
  for i = 0; i < n; i = i + 1 {
    if a[i] != b[i] {
      d = d + 1
    }
  }
  Some(d)
}

///|
/// Split a "v:j:cdr3" key into its CDR3 amino-acid part (the segment after
/// the second colon).
fn scanpy_ir_key_cdr3(key : String) -> String {
  let mut colons = 0
  for i = 0; i < key.length(); i = i + 1 {
    if key[i] == ':' {
      colons = colons + 1
      if colons == 2 {
        return key[i + 1:].to_owned()
      }
    }
  }
  ""
}

///|
/// Clonotype network as a `SparseGraph` over clonotypes (self-weight 1):
///
/// - `mode="exact"` (default): an edge joins two clonotypes whenever they
///   share any identical CDR3 sequence among their primary keys (scirpy's
///   clonotype-network distance "clonotype"),
/// - `mode="hamming"`: the primary CDR3s (VDJ preferred, else VJ) are linked
///   when of equal length and Hamming distance <= `max_distance`.
pub fn scanpy_ir_clonotype_network(
  data : ScanpyIrData,
  mode~ : String = "exact",
  max_distance~ : Int32 = 1,
) -> SparseGraph {
  let n = data.clonotypes.length()
  let diag = Array::make(n, 1.0)
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    adj.push([])
  }
  if n == 0 {
    return SparseGraph::{ n, diag, adj, }
  }
  // primary CDR3 per clonotype (vdj first, else vj)
  let primary : Array[String] = Array(capacity=n)
  let all_cdr3s : Array[Array[String]] = Array(capacity=n)
  for ct in data.clonotypes {
    let sets : Array[String] = Array(capacity=2)
    if ct.vj_key != "" {
      sets.push(scanpy_ir_key_cdr3(ct.vj_key))
    }
    if ct.vdj_key != "" {
      sets.push(scanpy_ir_key_cdr3(ct.vdj_key))
    }
    all_cdr3s.push(sets)
    let p = if ct.vdj_key != "" {
      scanpy_ir_key_cdr3(ct.vdj_key)
    } else if ct.vj_key != "" {
      scanpy_ir_key_cdr3(ct.vj_key)
    } else {
      ""
    }
    primary.push(p)
  }
  for a = 0; a < n; a = a + 1 {
    for b = a + 1; b < n; b = b + 1 {
      let linked = match mode {
        "hamming" => {
          let pa = primary[a]
          let pb = primary[b]
          if pa == "" || pb == "" {
            false
          } else {
            match scanpy_ir_hamming(pa, pb) {
              Some(d) => d <= max_distance
              None => false
            }
          }
        }
        _ => {
          // exact: share any identical CDR3
          let mut found = false
          for ca in all_cdr3s[a] {
            if found {
              break
            }
            for cb in all_cdr3s[b] {
              if ca != "" && ca == cb {
                found = true
                break
              }
            }
          }
          found
        }
      }
      if linked {
        adj[a].push((b, 1.0))
        adj[b].push((a, 1.0))
      }
    }
  }
  SparseGraph::{ n, diag, adj, }
}
