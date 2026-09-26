///|
/// Bio.SeqIO.CifIO - extract sequences from mmCIF files (cif-seqres and
/// cif-atom formats), mirroring Biopython's `Bio.SeqIO.CifIO`.
///
///   - cif-seqres: one record per entity from the `entity_poly_seq`
///     category (3-letter codes mapped to 1-letter codes).
///   - cif-atom:   one record per chain from the `atom_site` category,
///     residues aggregated in file order and de-duplicated.

///|
/// One residue of an entity_poly_seq row: residue number and monomer id.
struct CifResidue {
  num : Int
  has_num : Bool
  mon : String
}

///|
/// Fetch a loop_ category's headers and data rows from a block, following
/// the storage convention of parse_mmcif: the first row holds the
/// category.attribute header names, the remaining rows hold data.
fn cif_category(
  block : MmcifDataBlock,
  name : String,
) -> (Array[String], Array[Array[String]]) {
  match block.get_category(name) {
    Some(data) =>
      if data.length() > 1 {
        let headers = data[0]
        let rows : Array[Array[String]] = Array(capacity=0)
        for i = 1; i < data.length(); i = i + 1 {
          rows.push(data[i])
        }
        (headers, rows)
      } else {
        ([], [])
      }
    _ => ([], [])
  }
}

///|
/// Parse cif-seqres content and return one SeqRecord per entity.
///
/// Requires the `entity_poly_seq` category. Residues are grouped by
/// entity_id (in first-appearance order). When residue numbers are
/// present, each entity's sequence is assembled in residue-number order;
/// otherwise file order is kept. Each 3-letter mon_id is mapped to a
/// 1-letter code via `three_to_one_aa` (unknown residues become "X").
/// When the category is missing the function returns an empty array
/// without raising (spec 5.2.3.2).
pub fn parse_cif_seqres(content : String) -> Array[SeqRecord] {
  let block = parse_mmcif(content)
  let records : Array[SeqRecord] = Array(capacity=0)
  let (headers, rows) = cif_category(block, "entity_poly_seq")
  if rows.length() == 0 {
    return records
  }
  let entity_idx = mmcif_find_index(headers, "entity_poly_seq.entity_id")
  let num_idx = mmcif_find_index(headers, "entity_poly_seq.num")
  let mon_idx = mmcif_find_index(headers, "entity_poly_seq.mon_id")
  if entity_idx < 0 || mon_idx < 0 {
    return records
  }
  let entity_ids : Array[String] = Array(capacity=0)
  let entity_residues : Array[Array[CifResidue]] = Array(capacity=0)
  for row in rows {
    let entity = if entity_idx < row.length() { row[entity_idx] } else { "" }
    if entity.length() == 0 {
      continue
    }
    let mon = if mon_idx < row.length() { row[mon_idx] } else { "" }
    let (num, has_num) = if num_idx >= 0 &&
      num_idx < row.length() &&
      row[num_idx].length() > 0 &&
      row[num_idx].unsafe_get(0).to_int() >= '0'.to_int() &&
      row[num_idx].unsafe_get(0).to_int() <= '9'.to_int() {
      (parse_int(row[num_idx]), true)
    } else {
      (0, false)
    }
    let mut slot = -1
    for k = 0; k < entity_ids.length(); k = k + 1 {
      if entity_ids[k] == entity {
        slot = k
        break
      }
    }
    if slot < 0 {
      entity_ids.push(entity)
      entity_residues.push([])
      slot = entity_ids.length() - 1
    }
    entity_residues[slot].push({ num, has_num, mon, })
  }
  for k = 0; k < entity_ids.length(); k = k + 1 {
    let entity = entity_ids[k]
    let residues = entity_residues[k]
    // Determine whether all residues carry usable, distinct numbers
    let mut all_numbered = residues.length() > 0
    for r in residues {
      if !r.has_num {
        all_numbered = false
        break
      }
    }
    let ordered : Array[CifResidue] = if all_numbered {
      let sorted : Array[CifResidue] = Array(capacity=0)
      for r in residues {
        sorted.push(r)
      }
      // insertion sort by residue number (entity_poly_seq rows are few)
      for i = 1; i < sorted.length(); i = i + 1 {
        let key = sorted[i]
        let mut j = i - 1
        while j >= 0 && sorted[j].num > key.num {
          sorted[j + 1] = sorted[j]
          j = j - 1
        }
        sorted[j + 1] = key
      }
      sorted
    } else {
      residues
    }
    let seq_buf = StringBuilder()
    for r in ordered {
      seq_buf.write_string(three_to_one_aa(r.mon))
    }
    records.push(
      SeqRecord::new(
        Seq::new(seq_buf.to_string()),
        id=entity,
        name=entity,
        description="sequence from cif-seqres, entity " + entity,
      ),
    )
  }
  records
}

///|
/// Parse cif-atom content and return one SeqRecord per chain.
///
/// Reads the `atom_site` category, aggregates residues per
/// `label_asym_id` in file order, de-duplicating by
/// (label_seq_id, label_comp_id) so each resolved residue appears once.
/// Water molecules are skipped and non-standard residues become "X".
/// Chains with no atomic coordinates are omitted; a missing `atom_site`
/// category returns an empty array (spec 5.2.3.2).
pub fn parse_cif_atom(content : String) -> Array[SeqRecord] {
  let block = parse_mmcif(content)
  let records : Array[SeqRecord] = Array(capacity=0)
  let (headers, rows) = cif_category(block, "atom_site")
  if rows.length() == 0 {
    return records
  }
  let chain_idx = mmcif_find_index(headers, "atom_site.label_asym_id")
  let comp_idx = mmcif_find_index(headers, "atom_site.label_comp_id")
  let seq_idx = mmcif_find_index(headers, "atom_site.label_seq_id")
  if chain_idx < 0 || comp_idx < 0 {
    return records
  }
  let chain_ids : Array[String] = Array(capacity=0)
  let chain_seqs : Array[StringBuilder] = Array(capacity=0)
  let chain_keys : Array[Array[String]] = Array(capacity=0)
  for row in rows {
    let chain = if chain_idx < row.length() { row[chain_idx] } else { "" }
    let comp = if comp_idx < row.length() { row[comp_idx] } else { "" }
    if chain.length() == 0 || comp.length() == 0 {
      continue
    }
    if comp == "HOH" || comp == "WAT" {
      continue
    }
    let seq_key = if seq_idx >= 0 && seq_idx < row.length() {
      row[seq_idx]
    } else {
      ""
    }
    let dedup_key = seq_key + "|" + comp
    let mut slot = -1
    for k = 0; k < chain_ids.length(); k = k + 1 {
      if chain_ids[k] == chain {
        slot = k
        break
      }
    }
    if slot < 0 {
      chain_ids.push(chain)
      chain_seqs.push(StringBuilder())
      let keys : Array[String] = Array(capacity=0)
      keys.push(dedup_key)
      chain_keys.push(keys)
      chain_seqs[chain_ids.length() - 1].write_string(three_to_one_aa(comp))
      continue
    }
    let keys = chain_keys[slot]
    let mut seen = false
    for key in keys {
      if key == dedup_key {
        seen = true
        break
      }
    }
    if !seen {
      keys.push(dedup_key)
      chain_seqs[slot].write_string(three_to_one_aa(comp))
    }
  }
  for k = 0; k < chain_ids.length(); k = k + 1 {
    let chain = chain_ids[k]
    records.push(
      SeqRecord::new(
        Seq::new(chain_seqs[k].to_string()),
        id=chain,
        name=chain,
        description="sequence from cif-atom, chain " + chain,
      ),
    )
  }
  records
}
