///|
/// bp3_misc.mbt — Biopython round 3: miscellaneous features
/// #17 motifs.minimal — minimal motif format parser
/// #19 PDBExceptions — PDB-specific exception/warning hierarchy
/// #22 GeneMark — heuristic gene prediction
/// #23 AbstractPropertyMap — PDB property mapping framework
/// #25 ComparativeScatter — comparative genome scatter plot

///|
/// #17 motifs.minimal — minimal motif format parser
pub struct Bp3MinimalMotif {
  id : String
  consensus : String
  instances : Array[String]
  length_ : Int32
}

///|
pub fn bp3_parse_minimal_motif(text : String) -> Array[Bp3MinimalMotif] {
  let results = Array(capacity=0)
  let lines = text.split("\n").to_array()
  let mut current_id = ""
  let mut current_consensus = ""
  let mut current_instances : Array[String] = []
  for i in 0..<lines.length() {
    let line = lines[i].to_owned().trim().to_owned()
    if line.length() == 0 || line.has_prefix("#") {
      continue
    }
    if line.has_prefix("ID") {
      if current_id.length() > 0 {
        results.push({
          id: current_id,
          consensus: current_consensus,
          instances: current_instances,
          length_: current_consensus.length(),
        })
      }
      current_id = line[2:].to_owned().trim().to_owned()
      current_consensus = ""
      current_instances = []
    } else if line.has_prefix("CONSENSUS") {
      current_consensus = line[9:].to_owned().trim().to_owned()
    } else if line.has_prefix("INSTANCE") {
      let inst = line[8:].to_owned().trim().to_owned()
      if inst.length() > 0 {
        current_instances.push(inst)
      }
    }
  }
  if current_id.length() > 0 {
    results.push({
      id: current_id,
      consensus: current_consensus,
      instances: current_instances,
      length_: current_consensus.length(),
    })
  }
  results
}

///|
pub fn bp3_write_minimal_motif(motifs : Array[Bp3MinimalMotif]) -> String {
  let sb = StringBuilder()
  for i in 0..<motifs.length() {
    let m = motifs[i]
    let _ = sb.write("ID " + m.id + "\n")
    let _ = sb.write("CONSENSUS " + m.consensus + "\n")
    for j in 0..<m.instances.length() {
      let _ = sb.write("INSTANCE " + m.instances[j] + "\n")
    }
    let _ = sb.write("//\n")
  }
  sb.to_string()
}

///|
/// #19 PDBExceptions — PDB-specific exception/warning types
pub enum Bp3PDBException {
  PDBConstructionError(String)
  PDBConstructionWarning(String)
  PDBIOError(String)
  PDBParseError(String)
  PDBFormatError(String)
}

///|
pub fn bp3_pdb_exception_to_string(e : Bp3PDBException) -> String {
  match e {
    Bp3PDBException::PDBConstructionError(msg) => "PDBConstructionError: " + msg
    Bp3PDBException::PDBConstructionWarning(msg) =>
      "PDBConstructionWarning: " + msg
    Bp3PDBException::PDBIOError(msg) => "PDBIOError: " + msg
    Bp3PDBException::PDBParseError(msg) => "PDBParseError: " + msg
    Bp3PDBException::PDBFormatError(msg) => "PDBFormatError: " + msg
  }
}

///|
pub fn bp3_pdb_is_warning(e : Bp3PDBException) -> Bool {
  match e {
    Bp3PDBException::PDBConstructionWarning(_) => true
    _ => false
  }
}

///|
pub fn bp3_pdb_is_error(e : Bp3PDBException) -> Bool {
  !bp3_pdb_is_warning(e)
}

///|
/// #22 GeneMark — heuristic gene prediction
pub struct Bp3GeneMarkGene {
  start : Int32
  end : Int32
  strand : Int32
  gene_type : String
  score : Double
}

///|
pub fn bp3_genemark_predict(
  sequence : String,
  min_gene_len~ : Int32 = 30,
  gc_content~ : Double = 0.5,
) -> Array[Bp3GeneMarkGene] {
  let results = Array(capacity=0)
  let n = sequence.length()
  let start_codons = ["ATG", "GTG", "TTG"]
  let stop_codons = ["TAA", "TAG", "TGA"]
  let coding_threshold = 0.6 + (gc_content - 0.5) * 0.4
  let _n3 = n - 3
  for i in 0..<_n3 {
    let codon = sequence[i:i + 3].to_owned()
    let mut is_start = false
    for s in 0..<start_codons.length() {
      if start_codons[s] == codon {
        is_start = true
      }
    }
    if !is_start {
      continue
    }
    let mut j = i + 3
    let mut found_stop = false
    while j + 3 <= n {
      let stop_codon = sequence[j:j + 3].to_owned()
      let mut is_stop = false
      for s in 0..<stop_codons.length() {
        if stop_codons[s] == stop_codon {
          is_stop = true
        }
      }
      if is_stop {
        found_stop = true
        break
      }
      j = j + 3
    }
    if found_stop && j - i >= min_gene_len {
      let gene_seq = sequence[i:j].to_owned()
      let coding_score = bp3_genemark_coding_score(gene_seq, gc_content)
      if coding_score > coding_threshold {
        results.push({
          start: i,
          end: j + 3,
          strand: 1,
          gene_type: "CDS",
          score: coding_score,
        })
      }
    }
  }
  let rc = bp3_genemark_reverse_complement(sequence)
  let _n3 = n - 3
  for i in 0..<_n3 {
    let codon = rc[i:i + 3].to_owned()
    let mut is_start = false
    for s in 0..<start_codons.length() {
      if start_codons[s] == codon {
        is_start = true
      }
    }
    if !is_start {
      continue
    }
    let mut j = i + 3
    let mut found_stop = false
    while j + 3 <= n {
      let stop_codon = rc[j:j + 3].to_owned()
      let mut is_stop = false
      for s in 0..<stop_codons.length() {
        if stop_codons[s] == stop_codon {
          is_stop = true
        }
      }
      if is_stop {
        found_stop = true
        break
      }
      j = j + 3
    }
    if found_stop && j - i >= min_gene_len {
      let gene_seq = rc[i:j].to_owned()
      let coding_score = bp3_genemark_coding_score(gene_seq, gc_content)
      if coding_score > coding_threshold {
        results.push({
          start: n - j - 3,
          end: n - i,
          strand: -1,
          gene_type: "CDS",
          score: coding_score,
        })
      }
    }
  }
  results
}

///|
fn bp3_genemark_coding_score(seq : String, gc_content : Double) -> Double {
  let n = seq.length()
  if n < 6 {
    return 0.0
  }
  let mut gc3 = 0.0
  let mut n_codons = 0
  let mut i = 0
  while i + 3 <= n {
    let c3 = seq[i + 2:i + 3].to_owned()
    if c3 == "G" || c3 == "C" {
      gc3 = gc3 + 1.0
    }
    n_codons = n_codons + 1
    i = i + 3
  }
  if n_codons == 0 {
    return 0.0
  }
  let gc3_ratio = gc3 / n_codons.to_double()
  let expected = gc_content
  let score = 1.0 -
    (if gc3_ratio >= expected {
      gc3_ratio - expected
    } else {
      expected - gc3_ratio
    }) *
    2.0
  if score < 0.0 {
    0.0
  } else {
    score
  }
}

///|
fn bp3_genemark_reverse_complement(seq : String) -> String {
  let comp : Map[Int32, Int32] = Map([
    (65, 84),
    (84, 65),
    (71, 67),
    (67, 71),
    (97, 116),
    (116, 97),
    (103, 99),
    (99, 103),
  ])
  let sb = StringBuilder()
  for i in 0..<seq.length() {
    let ch = seq[seq.length() - 1 - i]
    let v = comp.get(ch.to_int())
    match v {
      Some(c) => {
        let _ = sb.write_char(Int32::unsafe_to_char(c))
      }
      None => {
        let _ = sb.write_char(ch.unsafe_to_char())
      }
    }
  }
  sb.to_string()
}

///|
/// #23 AbstractPropertyMap — PDB property mapping framework
pub struct Bp3PropertyMap {
  keys : Array[(String, Int32, String)]
  values : Array[String]
  property_name : String
}

///|
pub fn Bp3PropertyMap::new(property_name : String) -> Bp3PropertyMap {
  { keys: [], values: [], property_name, }
}

///|
pub fn Bp3PropertyMap::add(
  self : Bp3PropertyMap,
  chain_id : String,
  res_id : Int32,
  res_name : String,
  value : String,
) -> Unit {
  self.keys.push((chain_id, res_id, res_name))
  self.values.push(value)
}

///|
pub fn Bp3PropertyMap::get(
  self : Bp3PropertyMap,
  chain_id : String,
  res_id : Int32,
) -> String? {
  for i in 0..<self.keys.length() {
    let (c, r, _) = self.keys[i]
    if c == chain_id && r == res_id {
      return Some(self.values[i])
    }
  }
  None
}

///|
pub fn Bp3PropertyMap::length(self : Bp3PropertyMap) -> Int32 {
  self.keys.length()
}

///|
pub fn Bp3PropertyMap::to_string(self : Bp3PropertyMap) -> String {
  let sb = StringBuilder()
  let _ = sb.write("# " + self.property_name + "\n")
  for i in 0..<self.keys.length() {
    let (c, r, name) = self.keys[i]
    let _ = sb.write(
      c + "\t" + r.to_string() + "\t" + name + "\t" + self.values[i] + "\n",
    )
  }
  sb.to_string()
}

///|
/// #25 ComparativeScatter — comparative genome scatter plot data
pub struct Bp3ComparativeScatter {
  x_coords : Array[Double]
  y_coords : Array[Double]
  labels : Array[String]
  x_label : String
  y_label : String
  title : String
}

///|
pub fn bp3_comparative_scatter(
  genome1_positions : Array[Int32],
  genome2_positions : Array[Int32],
  labels~ : Array[String] = [],
  title~ : String = "Comparative Genome Scatter",
) -> Bp3ComparativeScatter {
  let n = if genome1_positions.length() < genome2_positions.length() {
    genome1_positions.length()
  } else {
    genome2_positions.length()
  }
  let x = Array(capacity=n)
  let y = Array(capacity=n)
  for i in 0..<n {
    x.push(genome1_positions[i].to_double())
    y.push(genome2_positions[i].to_double())
  }
  let lbl = if labels.length() > 0 {
    labels
  } else {
    let l = Array(capacity=n)
    for i in 0..<n {
      l.push("gene_" + i.to_string())
    }
    l
  }
  {
    x_coords: x,
    y_coords: y,
    labels: lbl,
    x_label: "Genome 1 Position",
    y_label: "Genome 2 Position",
    title,
  }
}

///|
pub fn bp3_comparative_scatter_to_tsv(
  scatter : Bp3ComparativeScatter,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("label\t" + scatter.x_label + "\t" + scatter.y_label + "\n")
  for i in 0..<scatter.x_coords.length() {
    let label = if i < scatter.labels.length() {
      scatter.labels[i]
    } else {
      "point_" + i.to_string()
    }
    let _ = sb.write(
      label +
      "\t" +
      scatter.x_coords[i].to_string() +
      "\t" +
      scatter.y_coords[i].to_string() +
      "\n",
    )
  }
  sb.to_string()
}

///|
pub fn bp3_comparative_scatter_correlation(
  scatter : Bp3ComparativeScatter,
) -> Double {
  let n = scatter.x_coords.length()
  if n < 2 {
    return 0.0
  }
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut sum_y2 = 0.0
  for i in 0..<n {
    let x = scatter.x_coords[i]
    let y = scatter.y_coords[i]
    sum_x = sum_x + x
    sum_y = sum_y + y
    sum_xy = sum_xy + x * y
    sum_x2 = sum_x2 + x * x
    sum_y2 = sum_y2 + y * y
  }
  let nf = n.to_double()
  let num = nf * sum_xy - sum_x * sum_y
  let den = ((nf * sum_x2 - sum_x * sum_x) * (nf * sum_y2 - sum_y * sum_y)).sqrt()
  if den == 0.0 {
    0.0
  } else {
    num / den
  }
}

///|
pub fn _bp3_misc_dummy() -> Unit {
  let _ = Bp3MinimalMotif::{ id: "", consensus: "", instances: [], length_: 0, }
  let _ = Bp3PDBException::PDBConstructionError("")
  let _ = Bp3GeneMarkGene::{
    start: 0,
    end: 0,
    strand: 0,
    gene_type: "",
    score: 0.0,
  }
  let _ = Bp3PropertyMap::{ keys: [], values: [], property_name: "", }
  let _ = Bp3ComparativeScatter::{
    x_coords: [],
    y_coords: [],
    labels: [],
    x_label: "",
    y_label: "",
    title: "",
  }
}

///|
pub fn bp3_misc_init_warnings() -> Unit {
  ignore(PDBConstructionWarning(""))
  ignore(PDBIOError(""))
  ignore(PDBParseError(""))
  ignore(PDBFormatError(""))
}
