///|
/// Nexus format parsing module inspired by Biopython Bio.Nexus.
pub type NexusCell = String

///|
pub struct NexusRow {
  taxon : String
  data : Array[NexusCell]
}

///|
pub fn NexusRow::new(taxon : String, data : Array[NexusCell]) -> NexusRow {
  NexusRow::{ taxon, data }
}

///|
pub enum NexusBlockType {
  Data
  Trees
  Distances
  Characters
  Taxa
  Unknown(String)
}

///|
pub struct NexusBlock {
  block_type : NexusBlockType
  content : String
}

///|
pub fn NexusBlock::new(
  block_type : NexusBlockType,
  content : String,
) -> NexusBlock {
  NexusBlock::{ block_type, content }
}

///|
pub struct NexusFile {
  blocks : Array[NexusBlock]
  taxa : Array[String]
}

///|
pub fn NexusFile::new(
  blocks : Array[NexusBlock],
  taxa : Array[String],
) -> NexusFile {
  NexusFile::{ blocks, taxa }
}

///|
pub fn parse_nexus(content : String) -> NexusFile {
  let blocks : Array[NexusBlock] = Array::new()
  let taxa : Array[String] = Array::new()
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  while i < lines.length() {
    let line_s = lines[i].to_string()
    let line = line_s.trim().to_string()
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if line.has_prefix("#NEXUS") {
      i = i + 1
      continue
    }
    let lower_line = line.to_lower()
    if lower_line.has_prefix("begin ") {
      let raw_name = line[6:].trim().to_owned()
      let block_name : String = if raw_name.has_suffix(";") {
        let without_semicolon = raw_name[:raw_name.length() - 1].to_string()
        without_semicolon.to_upper().to_string()
      } else {
        raw_name.to_upper().to_string()
      }
      let block_type = match block_name {
        "DATA" => NexusBlockType::Data
        "TREES" => NexusBlockType::Trees
        "DISTANCES" => NexusBlockType::Distances
        "CHARACTERS" => NexusBlockType::Characters
        "TAXA" => NexusBlockType::Taxa
        _ => NexusBlockType::Unknown(block_name)
      }
      i = i + 1
      let mut block_content = ""
      while i < lines.length() {
        let inner_line = lines[i].to_string()
        let trimmed_inner = inner_line.trim().to_string()
        if trimmed_inner.to_lower().has_prefix("end;") {
          blocks.push(NexusBlock::new(block_type, block_content))
          if block_type is NexusBlockType::Taxa {
            let t = parse_taxa_from_string(block_content)
            let mut idx = 0
            while idx < t.length() {
              taxa.push(t[idx])
              idx = idx + 1
            }
          }
          i = i + 1
          break
        }
        block_content = block_content + trimmed_inner + "\n"
        i = i + 1
      }
      continue
    }
    i = i + 1
  }
  NexusFile::new(blocks, taxa)
}

///|
fn parse_taxa_from_string(content : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].to_string().trim().to_owned()
    let lower = line.to_lower()
    if lower.has_prefix("taxlabels") {
      let tokens = line.split(" ").collect()
      let mut t_idx = 0
      while t_idx < tokens.length() {
        let token = tokens[t_idx].trim().to_owned()
        if token.length() > 0 &&
          !lower.has_prefix(token.to_lower()) &&
          token != ";" {
          let clean_token : String = if token.has_suffix(";") {
            token[:token.length() - 1].to_string()
          } else {
            token.to_string()
          }
          if clean_token.length() > 0 {
            result.push(clean_token)
          }
        }
        t_idx = t_idx + 1
      }
      break
    }
    i = i + 1
  }
  result
}

///|
pub fn NexusFile::get_data_matrix(self : NexusFile) -> Array[NexusRow] {
  let result : Array[NexusRow] = Array::new()
  let mut block_idx = 0
  while block_idx < self.blocks.length() {
    let block = self.blocks[block_idx]
    if block.block_type is NexusBlockType::Data ||
      block.block_type is NexusBlockType::Characters {
      let lines = split_by_char(block.content, '\n'.to_int())
      let mut in_matrix = false
      let mut i = 0
      while i < lines.length() {
        let line = lines[i].to_string().trim().to_string()
        let lower = line.to_lower()
        if lower.has_prefix("matrix") {
          in_matrix = true
          i = i + 1
          continue
        }
        if line.has_prefix(";") {
          in_matrix = false
          break
        }
        if in_matrix && line.length() > 0 {
          let mut idx = 0
          while idx < line.length() {
            let c = line.unsafe_get(idx).unsafe_to_char()
            if c == ' ' || c == '\t' {
              break
            }
            idx = idx + 1
          }
          if idx > 0 && idx < line.length() {
            let taxon = line[:idx].trim().to_string()
            let data = line[idx:].trim().to_string()
            if taxon.length() > 0 && data.length() > 0 {
              result.push(NexusRow::new(taxon, [data]))
            }
          }
        }
        i = i + 1
      }
    }
    block_idx = block_idx + 1
  }
  result
}

///|
pub fn NexusFile::get_trees(self : NexusFile) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array::new()
  let mut block_idx = 0
  while block_idx < self.blocks.length() {
    let block = self.blocks[block_idx]
    if block.block_type is NexusBlockType::Trees {
      let lines = split_by_char(block.content, '\n'.to_int())
      let mut i = 0
      while i < lines.length() {
        let line = lines[i].to_string().trim().to_string()
        let lower = line.to_lower()
        if lower.has_prefix("tree") {
          let mut j = 0
          while j < line.length() {
            let c = line.unsafe_get(j).unsafe_to_char()
            if c == ' ' || c == '\t' {
              break
            }
            j = j + 1
          }
          let after_tree = line[j:].trim().to_string()
          let mut eq_pos = -1
          let mut k = 0
          while k < after_tree.length() {
            if after_tree.unsafe_get(k).to_int() == '='.to_int() {
              eq_pos = k
              break
            }
            k = k + 1
          }
          if eq_pos >= 0 {
            let name = after_tree[:eq_pos].trim().to_string()
            let newick = after_tree[eq_pos + 1:].trim().to_string()
            result.push((name, newick))
          }
        }
        i = i + 1
      }
    }
    block_idx = block_idx + 1
  }
  result
}

///|
pub fn NexusFile::get_distance_matrix(
  self : NexusFile,
) -> Array[(String, Array[Double])] {
  let _ = self
  []
}

///|
pub fn NexusFile::get_num_taxa(self : NexusFile) -> Int {
  self.taxa.length()
}

///|
pub fn sample_nexus_content() -> String {
  "#NEXUS\n" +
  "Begin Taxa;\n" +
  "  Dimensions NTax=4;\n" +
  "  TaxLabels A B C D;\n" +
  "End;\n" +
  "Begin Characters;\n" +
  "  Dimensions NChar=10;\n" +
  "  Format Datatype=DNA Missing=N Gap=-;\n" +
  "Matrix\n" +
  "  A  ATGCATGCAT\n" +
  "  B  ATGCATGCAT\n" +
  "  C  ATGCATGCAT\n" +
  "  D  ATGCATGCAT\n" +
  ";\n" +
  "End;\n" +
  "Begin Trees;\n" +
  "  Tree tree1 = ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\n" +
  "End;\n"
}

///|
pub fn sample_nexus_file() -> NexusFile {
  parse_nexus(sample_nexus_content())
}
