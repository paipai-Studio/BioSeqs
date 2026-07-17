///|
/// mmCIF (Macromolecular Crystallographic Information File) format parsing.
/// Inspired by Biopython Bio.PDB.MMCIFParser.

///|
/// An mmCIF data item.
pub struct MmcifDataItem {
  name : String
  value : String
}

///|
/// Construct an MmcifDataItem.
pub fn MmcifDataItem::new(name : String, value : String) -> MmcifDataItem {
  MmcifDataItem::{ name, value }
}

///|
/// An mmCIF data block containing categories.
pub struct MmcifDataBlock {
  mut name : String
  categories : Map[String, Array[Array[String]]]
}

///|
/// Construct an MmcifDataBlock.
pub fn MmcifDataBlock::new(name : String) -> MmcifDataBlock {
  MmcifDataBlock::{ name, categories: Map([], capacity=32) }
}

///|
/// Get category data by name.
pub fn MmcifDataBlock::get_category(
  self : MmcifDataBlock,
  name : String,
) -> Array[Array[String]]? {
  self.categories.get(name)
}

///|
/// Set category data.
pub fn MmcifDataBlock::set_category(
  self : MmcifDataBlock,
  name : String,
  data : Array[Array[String]],
) -> Unit {
  self.categories[name] = data
}

///|
/// Parse mmCIF content into a data block.
pub fn parse_mmcif(content : String) -> MmcifDataBlock {
  let lines = split_by_char(content, '\n'.to_int())
  let block = MmcifDataBlock::new("default")
  let mut current_category = ""
  let mut current_data : Array[Array[String]] = Array::new()
  let mut headers : Array[String] = Array::new()
  let mut in_loop = false

  let mut i = 0
  while i < lines.length() {
    let line = trim(lines[i])
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' {
      i = i + 1
      continue
    }

    if starts_with(line, 0, "data_") {
      let name = substring(line, 5, line.length() - 5)
      block.name = if name.length() > 0 { name } else { "default" }
      in_loop = false
      current_category = ""
      i = i + 1
      continue
    }

    if starts_with(line, 0, "loop_") {
      in_loop = true
      headers = Array::new()
      current_data = Array::new()
      i = i + 1
      continue
    }

    if starts_with(line, 0, "_") {
      let parts = split_by_char(line, ' '.to_int())
      if parts.length() >= 1 {
        let cat_attr = substring(parts[0], 1, parts[0].length() - 1)
        let dot_idx = {
          let mut idx = -1
          let mut j = 0
          while j < cat_attr.length() {
            if cat_attr.unsafe_get(j).unsafe_to_char() == '.' {
              idx = j
              break
            }
            j = j + 1
          }
          idx
        }
        let cat_name = if dot_idx >= 0 {
          substring(cat_attr, 0, dot_idx)
        } else {
          cat_attr
        }

        if in_loop {
          headers.push(cat_attr)
        } else {
          let value = if parts.length() > 1 {
            let mut v = ""
            let mut j = 1
            while j < parts.length() {
              if j > 1 {
                v = v + " "
              }
              v = v + parts[j]
              j = j + 1
            }
            v
          } else {
            ""
          }
          current_category = cat_name
          let existing : Array[Array[String]] = match
            block.categories.get(cat_name) {
            Some(arr) => arr
            None => Array::new()
          }
          existing.push([cat_attr, value])
          block.categories[cat_name] = existing
        }
      }
      i = i + 1
      continue
    }

    if in_loop && headers.length() > 0 {
      let values = split_by_char(line, ' '.to_int())
      let filtered : Array[String] = Array::new()
      for v in values {
        let tv = trim(v)
        if tv.length() > 0 {
          filtered.push(tv)
        }
      }
      if filtered.length() > 0 {
        current_data.push(filtered)
      }
    } else if current_category.length() > 0 {
      let existing : Array[Array[String]] = match
        block.categories.get(current_category) {
        Some(arr) => arr
        None => Array::new()
      }
      existing.push([current_category + ".value", line])
      block.categories[current_category] = existing
    }
    i = i + 1
  }

  if in_loop && headers.length() > 0 && current_data.length() > 0 {
    let dot_idx = if headers.length() > 0 {
      let h = headers[0]
      let mut idx = -1
      let mut j = 0
      while j < h.length() {
        if h.unsafe_get(j).unsafe_to_char() == '.' {
          idx = j
          break
        }
        j = j + 1
      }
      idx
    } else {
      -1
    }
    let cat_name = if dot_idx >= 0 {
      substring(headers[0], 0, dot_idx)
    } else {
      "unknown"
    }
    let result : Array[Array[String]] = Array::new()
    result.push(headers)
    for row in current_data {
      result.push(row)
    }
    block.categories[cat_name] = result
  }

  block
}

///|
/// Get atom sites from an mmCIF data block.
pub fn MmcifDataBlock::get_atom_sites(
  self : MmcifDataBlock,
) -> Array[Array[String]] {
  match self.categories.get("atom_site") {
    Some(data) => data
    None => Array::new()
  }
}

///|
/// Get entity information from an mmCIF data block.
pub fn MmcifDataBlock::get_entity(
  self : MmcifDataBlock,
) -> Array[Array[String]] {
  match self.categories.get("entity") {
    Some(data) => data
    None => Array::new()
  }
}

///|
/// Get structure information from an mmCIF data block.
pub fn MmcifDataBlock::get_struct(
  self : MmcifDataBlock,
) -> Array[Array[String]] {
  match self.categories.get("struct") {
    Some(data) => data
    None => Array::new()
  }
}

///|
/// Get category headers (first row if it's a loop).
pub fn MmcifDataBlock::get_category_headers(
  self : MmcifDataBlock,
  name : String,
) -> Array[String] {
  match self.categories.get(name) {
    Some(data) =>
      if data.length() > 0 && starts_with(data[0][0], 0, "_") {
        data[0]
      } else {
        Array::new()
      }
    None => Array::new()
  }
}

///|
/// Get category data rows (excluding headers).
pub fn MmcifDataBlock::get_category_rows(
  self : MmcifDataBlock,
  name : String,
) -> Array[Array[String]] {
  match self.categories.get(name) {
    Some(data) =>
      if data.length() > 0 && starts_with(data[0][0], 0, "_") {
        let result : Array[Array[String]] = Array::new()
        let mut j = 1
        while j < data.length() {
          result.push(data[j])
          j = j + 1
        }
        result
      } else {
        data
      }
    None => Array::new()
  }
}

///|
/// Convert mmCIF atom sites to Atom objects (compatible with PDB module).
pub fn MmcifDataBlock::to_atoms(self : MmcifDataBlock) -> Array[Atom] {
  let sites = self.get_atom_sites()
  let atoms : Array[Atom] = Array::new()
  if sites.length() < 2 {
    return atoms
  }

  let headers = sites[0]
  let rows_len = sites.length() - 1
  let rows : Array[Array[String]] = Array::new()
  let mut k = 1
  while k < sites.length() {
    rows.push(sites[k])
    k = k + 1
  }

  let name_idx = mmcif_find_index(headers, "atom_site.label_atom_id")
  let resname_idx = mmcif_find_index(headers, "atom_site.label_comp_id")
  let chainid_idx = mmcif_find_index(headers, "atom_site.label_asym_id")
  let resseq_idx = mmcif_find_index(headers, "atom_site.label_seq_id")
  let x_idx = mmcif_find_index(headers, "atom_site.Cartn_x")
  let y_idx = mmcif_find_index(headers, "atom_site.Cartn_y")
  let z_idx = mmcif_find_index(headers, "atom_site.Cartn_z")
  let occ_idx = mmcif_find_index(headers, "atom_site.occupancy")
  let bfactor_idx = mmcif_find_index(headers, "atom_site.B_iso_or_equiv")
  let element_idx = mmcif_find_index(headers, "atom_site.type_symbol")

  for row in rows {
    if name_idx >= 0 && x_idx >= 0 && y_idx >= 0 && z_idx >= 0 {
      let atom_name = if name_idx < row.length() { row[name_idx] } else { "" }

      let x_val = if x_idx < row.length() { row[x_idx] } else { "0" }
      let x = match parse_double(x_val) {
        Some(v) => v
        None => 0.0
      }

      let y_val = if y_idx < row.length() { row[y_idx] } else { "0" }
      let y = match parse_double(y_val) {
        Some(v) => v
        None => 0.0
      }

      let z_val = if z_idx < row.length() { row[z_idx] } else { "0" }
      let z = match parse_double(z_val) {
        Some(v) => v
        None => 0.0
      }

      let resname = if resname_idx < row.length() {
        row[resname_idx]
      } else {
        "UNK"
      }

      let chainid = if chainid_idx < row.length() &&
        row[chainid_idx].length() > 0 {
        row[chainid_idx].unsafe_get(0).unsafe_to_char()
      } else {
        ' '
      }

      let resseq = if resseq_idx < row.length() {
        parse_int(row[resseq_idx])
      } else {
        1
      }

      let occ_val = if occ_idx < row.length() { row[occ_idx] } else { "1" }
      let occupancy = match parse_double(occ_val) {
        Some(v) => v
        None => 1.0
      }

      let bfactor_val = if bfactor_idx < row.length() {
        row[bfactor_idx]
      } else {
        "0"
      }
      let bfactor = match parse_double(bfactor_val) {
        Some(v) => v
        None => 0.0
      }

      let element = if element_idx < row.length() {
        row[element_idx]
      } else {
        ""
      }

      atoms.push(
        Atom::new(
          name=atom_name,
          coord=Vector3::new(x, y, z),
          resname~,
          chainid~,
          resseq~,
          occupancy~,
          bfactor~,
          element~,
        ),
      )
    }
  }

  atoms
}

///|
/// Helper function: find index of string in array.
fn mmcif_find_index(arr : Array[String], target : String) -> Int {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == target {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Sample mmCIF content for testing.
pub fn sample_mmcif_content() -> String {
  "data_1ABC\n" +
  "#\n" +
  "_entry.id   1ABC\n" +
  "#\n" +
  "loop_\n" +
  "_atom_site.group_PDB\n" +
  "_atom_site.id\n" +
  "_atom_site.type_symbol\n" +
  "_atom_site.label_atom_id\n" +
  "_atom_site.label_comp_id\n" +
  "_atom_site.label_asym_id\n" +
  "_atom_site.label_seq_id\n" +
  "_atom_site.Cartn_x\n" +
  "_atom_site.Cartn_y\n" +
  "_atom_site.Cartn_z\n" +
  "_atom_site.occupancy\n" +
  "_atom_site.B_iso_or_equiv\n" +
  "ATOM   1  N  N   MET A 1  10.0  20.0  30.0  1.0  20.0\n" +
  "ATOM   2  CA C  CA  MET A 1  11.0  21.0  31.0  1.0  21.0\n" +
  "ATOM   3  C  C  C   MET A 1  12.0  22.0  32.0  1.0  22.0\n" +
  "ATOM   4  O  O  O   MET A 1  13.0  23.0  33.0  1.0  23.0\n" +
  "ATOM   5  N  N   ALA A 2  14.0  24.0  34.0  1.0  24.0\n" +
  "ATOM   6  CA C  CA  ALA A 2  15.0  25.0  35.0  1.0  25.0\n" +
  "#\n" +
  "loop_\n" +
  "_entity.id\n" +
  "_entity.type\n" +
  "_entity.src_method\n" +
  "_entity.pdbx_description\n" +
  "1  polymer  man  Alpha-hemoglobin\n" +
  "2  water    nat  Water\n" +
  "#\n" +
  "_struct.title  Crystal structure of human hemoglobin\n"
}

///|
/// Create a sample MmcifDataBlock for testing.
pub fn sample_mmcif_block() -> MmcifDataBlock {
  parse_mmcif(sample_mmcif_content())
}
