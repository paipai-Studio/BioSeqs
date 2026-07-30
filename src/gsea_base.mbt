///|
/// Gene set I/O and collection utilities - inspired by Bioconductor {GSEABase}.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/GSEABase.html
///   - Morgan M, Falcon S, Gentleman R. (2023) GSEABase: Gene set enrichment data
///     structures and methods. R package version 1.64.0.

// ---- helper string ops ----------------------------------------------------

fn gs_str_lower(s : String) -> String {
  let mut result = ""
  let mut i = 0
  let n = s.length()
  while i < n {
    let c = s.unsafe_get(i)
    if c >= 65 && c <= 90 {
      let c2 = (c + 32).to_char()
      match c2 {
        Some(ch) => result = result + ch.to_string()
        None => ()
      }
    } else {
      result = result + s[i:i+1].to_owned()
    }
    i = i + 1
  }
  result
}

fn gs_str_upper(s : String) -> String {
  let mut result = ""
  let mut i = 0
  let n = s.length()
  while i < n {
    let c = s.unsafe_get(i)
    if c >= 97 && c <= 122 {
      let c2 = (c - 32).to_char()
      match c2 {
        Some(ch) => result = result + ch.to_string()
        None => ()
      }
    } else {
      result = result + s[i:i+1].to_owned()
    }
    i = i + 1
  }
  result
}

fn gs_trim(s : String) -> String {
  let mut start = 0
  let n = s.length()
  while start < n {
    let c = s.unsafe_get(start)
    if c == 32 || c == 9 || c == 10 || c == 13 {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let c = s.unsafe_get(end - 1)
    if c == 32 || c == 9 || c == 10 || c == 13 {
      end = end - 1
    } else {
      break
    }
  }
  if start == 0 && end == n { s } else { s[start:end].to_owned() }
}

fn gs_split(s : String, sep_u16 : UInt16) -> Array[String] {
  let out : Array[String] = Array::new()
  let mut buf = StringBuilder::new()
  let n = s.length()
  let mut i = 0
  while i < n {
    let c = s.unsafe_get(i)
    if c == sep_u16 {
      out.push(buf.to_string())
      buf = StringBuilder::new()
    } else {
      let ch = c.to_char()
      match ch {
        Some(cv) => buf.write_char(cv)
        None => ()
      }
    }
    i = i + 1
  }
  out.push(buf.to_string())
  out
}

// ---- types -----------------------------------------------------------------

pub enum GeneSetCollectionType {
  GO_BiologicalProcess
  GO_MolecularFunction
  GO_CellularComponent
  KEGG
  Reactome
  Hallmark
  Canonical
  Positional
  Motif
  Computational
  Other
}

pub fn GeneSetCollectionType::from_string(s : String) -> GeneSetCollectionType {
  let sl = gs_str_lower(s)
  if sl.has_prefix("go_bp") || sl.has_prefix("bp") || sl.contains("biological_process") {
    GO_BiologicalProcess
  } else if sl.has_prefix("go_mf") || sl.has_prefix("mf") || sl.contains("molecular_function") {
    GO_MolecularFunction
  } else if sl.has_prefix("go_cc") || sl.has_prefix("cc") || sl.contains("cellular_component") {
    GO_CellularComponent
  } else if sl.has_prefix("kegg") {
    KEGG
  } else if sl.has_prefix("reactome") {
    Reactome
  } else if sl.has_prefix("hallmark") || sl.has_prefix("h:") {
    Hallmark
  } else if sl.has_prefix("canonical") || sl.has_prefix("c2:") {
    Canonical
  } else if sl.has_prefix("positional") || sl.has_prefix("c1:") {
    Positional
  } else if sl.has_prefix("motif") || sl.has_prefix("c3:") {
    Motif
  } else if sl.has_prefix("computational") || sl.has_prefix("c4:") {
    Computational
  } else {
    Other
  }
}

pub fn GeneSetCollectionType::show(self : GeneSetCollectionType) -> String {
  let s = match self {
    GO_BiologicalProcess => "GO_BiologicalProcess"
    GO_MolecularFunction => "GO_MolecularFunction"
    GO_CellularComponent => "GO_CellularComponent"
    KEGG => "KEGG"
    Reactome => "Reactome"
    Hallmark => "Hallmark"
    Canonical => "Canonical"
    Positional => "Positional"
    Motif => "Motif"
    Computational => "Computational"
    Other => "Other"
  }
  s.to_string()
}

pub struct GmtGeneSet {
  pub name : String
  pub description : String
  pub gene_ids : Array[String]
  pub mut collection_type : GeneSetCollectionType
  pub mut organism : String
  pub mut id : String
}

pub fn GmtGeneSet::new(
  name : String,
  description : String,
  gene_ids : Array[String],
) -> GmtGeneSet {
  let uniq : Array[String] = Array::new()
  let seen : Array[String] = Array::new()
  for x in gene_ids {
    let t = gs_trim(x)
    if t.length() == 0 { continue }
    let mut already = false
    let upper = gs_str_upper(t)
    for s in seen {
      if s == upper { already = true; break }
    }
    if !already {
      seen.push(upper)
      uniq.push(t)
    }
  }
  GmtGeneSet::{
    name: name,
    description: description,
    gene_ids: uniq,
    collection_type: Other,
    organism: "unknown".to_string(),
    id: name,
  }
}

pub fn GmtGeneSet::with_annotation(
  name : String,
  description : String,
  gene_ids : Array[String],
  collection_type : GeneSetCollectionType,
  organism : String,
  id : String,
) -> GmtGeneSet {
  let gs = GmtGeneSet::new(name, description, gene_ids)
  gs.collection_type = collection_type
  gs.organism = organism
  gs.id = id
  gs
}

pub fn GmtGeneSet::size(self : GmtGeneSet) -> Int {
  self.gene_ids.length()
}

pub fn GmtGeneSet::has_gene(self : GmtGeneSet, gene : String) -> Bool {
  let gu = gs_str_upper(gs_trim(gene))
  for g in self.gene_ids {
    if gs_str_upper(g) == gu { return true }
  }
  false
}

// ---- set math --------------------------------------------------------------

pub fn gene_set_overlap(a : GmtGeneSet, b : GmtGeneSet) -> Int {
  let mut cnt = 0
  for g in a.gene_ids {
    if b.has_gene(g) { cnt = cnt + 1 }
  }
  cnt
}

pub fn gene_set_jaccard(a : GmtGeneSet, b : GmtGeneSet) -> Double {
  let ov = gene_set_overlap(a, b).to_double()
  let asz = a.gene_ids.length().to_double()
  let bsz = b.gene_ids.length().to_double()
  let un = asz + bsz - ov
  if un == 0.0 { 0.0 } else { ov / un }
}

pub fn gene_set_overlap_coef(a : GmtGeneSet, b : GmtGeneSet) -> Double {
  let ov = gene_set_overlap(a, b)
  let mut asz = a.gene_ids.length()
  let bsz = b.gene_ids.length()
  if asz > bsz { asz = bsz }
  if asz == 0 { 0.0 } else { ov.to_double() / asz.to_double() }
}

pub fn gene_set_union_arr(a : GmtGeneSet, b : GmtGeneSet) -> Array[String] {
  let out : Array[String] = Array::new()
  let seen : Array[String] = Array::new()
  for g in a.gene_ids {
    let u = gs_str_upper(g)
    let mut f = false
    for s in seen {
      if s == u { f = true; break }
    }
    if !f {
      seen.push(u)
      out.push(g)
    }
  }
  for g in b.gene_ids {
    let u = gs_str_upper(g)
    let mut f = false
    for s in seen {
      if s == u { f = true; break }
    }
    if !f {
      seen.push(u)
      out.push(g)
    }
  }
  out
}

pub fn gene_set_intersect_arr(a : GmtGeneSet, b : GmtGeneSet) -> Array[String] {
  let out : Array[String] = Array::new()
  for g in a.gene_ids {
    if b.has_gene(g) { out.push(g) }
  }
  out
}

pub fn gene_set_setdiff(a : GmtGeneSet, b : GmtGeneSet) -> Array[String] {
  let out : Array[String] = Array::new()
  for g in a.gene_ids {
    if !b.has_gene(g) { out.push(g) }
  }
  out
}

// ---- collection ------------------------------------------------------------

pub struct GmtGeneSetCollection {
  pub sets : Array[GmtGeneSet]
  pub names : Array[String]
  pub mut size : Int
}

pub fn GmtGeneSetCollection::new(sets : Array[GmtGeneSet]) -> GmtGeneSetCollection {
  let names : Array[String] = Array::new()
  for s in sets {
    names.push(s.name)
  }
  GmtGeneSetCollection::{
    sets: sets,
    names: names,
    size: sets.length(),
  }
}

pub fn GmtGeneSetCollection::by_name(
  self : GmtGeneSetCollection,
  name : String,
) -> Option[GmtGeneSet] {
  let n = self.sets.length()
  let mut i = 0
  while i < n {
    if self.sets[i].name == name { return Some(self.sets[i]) }
    i = i + 1
  }
  None
}

pub fn GmtGeneSetCollection::at(
  self : GmtGeneSetCollection,
  i : Int,
) -> Option[GmtGeneSet] {
  if i < 0 || i >= self.sets.length() { None } else { Some(self.sets[i]) }
}

pub fn GmtGeneSetCollection::all_sizes(self : GmtGeneSetCollection) -> Array[Int] {
  let out : Array[Int] = Array::new()
  for s in self.sets { out.push(s.size()) }
  out
}

pub fn GmtGeneSetCollection::push_s(
  self : GmtGeneSetCollection,
  gs : GmtGeneSet,
) -> GmtGeneSetCollection {
  let col = self
  col.sets.push(gs)
  col.names.push(gs.name)
  col.size = col.size + 1
  col
}

pub fn GmtGeneSetCollection::summary(self : GmtGeneSetCollection) -> String {
  let buf = StringBuilder::new()
  buf.write_string("GmtGeneSetCollection: \{self.size} sets\n")
  if self.size == 0 { return buf.to_string() }
  let mut min_sz = self.sets[0].size()
  let mut max_sz = self.sets[0].size()
  let mut sum_sz = 0
  for s in self.sets {
    let sz = s.size()
    if sz < min_sz { min_sz = sz }
    if sz > max_sz { max_sz = sz }
    sum_sz = sum_sz + sz
  }
  let avg = sum_sz.to_double() / self.size.to_double()
  buf.write_string("  size range: [\{min_sz}, \{max_sz}]  avg=\{avg}\n")
  let n_show = if self.size < 5 { self.size } else { 5 }
  let mut i = 0
  while i < n_show {
    let s = self.sets[i]
    buf.write_string("  - \{s.name} (\{s.size()} genes)  \{s.description}\n")
    i = i + 1
  }
  if self.size > 5 {
    buf.write_string("  ... (\{self.size - 5} more)\n")
  }
  buf.to_string()
}

pub fn collection_filter_by_size(
  col : GmtGeneSetCollection,
  min_size : Int,
  max_size : Int,
) -> GmtGeneSetCollection {
  let kept : Array[GmtGeneSet] = Array::new()
  for s in col.sets {
    let sz = s.size()
    if sz >= min_size && sz <= max_size { kept.push(s) }
  }
  GmtGeneSetCollection::new(kept)
}

pub fn collection_filter_by_gene(
  col : GmtGeneSetCollection,
  gene_universe : Array[String],
) -> GmtGeneSetCollection {
  let kept : Array[GmtGeneSet] = Array::new()
  let univ : Array[String] = Array::new()
  for g in gene_universe {
    univ.push(gs_str_upper(gs_trim(g)))
  }
  for s in col.sets {
    let filt : Array[String] = Array::new()
    for g in s.gene_ids {
      let gu = gs_str_upper(g)
      for u in univ {
        if u == gu {
          filt.push(g)
          break
        }
      }
    }
    if filt.length() > 0 {
      let ng = GmtGeneSet::new(s.name, s.description, filt)
      ng.collection_type = s.collection_type
      ng.organism = s.organism
      ng.id = s.id
      kept.push(ng)
    }
  }
  GmtGeneSetCollection::new(kept)
}

pub fn collection_gene_union(col : GmtGeneSetCollection) -> Array[String] {
  let out : Array[String] = Array::new()
  let seen : Array[String] = Array::new()
  for s in col.sets {
    for g in s.gene_ids {
      let u = gs_str_upper(g)
      let mut f = false
      for s2 in seen {
        if s2 == u { f = true; break }
      }
      if !f {
        seen.push(u)
        out.push(g)
      }
    }
  }
  out
}

pub fn collection_pairwise_jaccard(
  col : GmtGeneSetCollection,
) -> Array[Array[Double]] {
  let n = col.sets.length()
  let out : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      let v = if i == j {
        1.0
      } else {
        gene_set_jaccard(col.sets[i], col.sets[j])
      }
      row.push(v)
      j = j + 1
    }
    out.push(row)
    i = i + 1
  }
  out
}

// ---- GMT I/O ---------------------------------------------------------------

pub fn parse_gmt(content : String) -> GmtGeneSetCollection {
  let lines = gs_split(content, '\n'.to_int().to_uint16())
  let result : Array[GmtGeneSet] = Array::new()
  let n = lines.length()
  let mut idx = 0
  while idx < n {
    let raw_line = lines[idx]
    let line = gs_trim(raw_line)
    if line.length() == 0 { idx = idx + 1; continue }
    if line.has_prefix("#") { idx = idx + 1; continue }
    if line.has_prefix("!") { idx = idx + 1; continue }
    let parts = gs_split(line, '\t'.to_int().to_uint16())
    if parts.length() < 3 { idx = idx + 1; continue }
    let name = gs_trim(parts[0])
    let desc = gs_trim(parts[1])
    let genes : Array[String] = Array::new()
    let plen = parts.length()
    let mut i = 2
    while i < plen {
      let g = gs_trim(parts[i])
      if g.length() > 0 { genes.push(g) }
      i = i + 1
    }
    if name.length() > 0 && genes.length() > 0 {
      result.push(GmtGeneSet::new(name, desc, genes))
    }
    idx = idx + 1
  }
  GmtGeneSetCollection::new(result)
}

pub fn write_gmt(col : GmtGeneSetCollection) -> String {
  let buf = StringBuilder::new()
  for gs in col.sets {
    buf.write_string(gs.name)
    buf.write_char('\t')
    buf.write_string(gs.description)
    for g in gs.gene_ids {
      buf.write_char('\t')
      buf.write_string(g)
    }
    buf.write_char('\n')
  }
  buf.to_string()
}

// ---- GMX I/O ---------------------------------------------------------------

pub fn parse_gmx(content : String) -> GmtGeneSetCollection {
  let tab16 = '\t'.to_int().to_uint16()
  let lines = gs_split(content, '\n'.to_int().to_uint16())
  let nonblank : Array[Array[String]] = Array::new()
  for raw in lines {
    let line = gs_trim(raw)
    if line.length() == 0 { continue }
    if line.has_prefix("#") { continue }
    nonblank.push(gs_split(line, tab16))
  }
  if nonblank.length() < 2 { return GmtGeneSetCollection::new([]) }
  let headers = nonblank[0]
  let descs = nonblank[1]
  let n_sets = headers.length()
  let set_genes : Array[Array[String]] = Array::new()
  let mut set_i = 0
  while set_i < n_sets {
    set_genes.push([])
    set_i = set_i + 1
  }
  let mut row_i = 2
  while row_i < nonblank.length() {
    let row = nonblank[row_i]
    let mut j = 0
    while j < n_sets {
      if j >= row.length() { j = j + 1; continue }
      let cell = gs_trim(row[j])
      if cell.length() == 0 { j = j + 1; continue }
      set_genes[j].push(cell)
      j = j + 1
    }
    row_i = row_i + 1
  }
  let result : Array[GmtGeneSet] = Array::new()
  let mut j = 0
  while j < n_sets {
    let nm = gs_trim(headers[j])
    if nm.length() == 0 { j = j + 1; continue }
    let desc = if j < descs.length() { gs_trim(descs[j]) } else { "".to_string() }
    result.push(GmtGeneSet::new(nm, desc, set_genes[j]))
    j = j + 1
  }
  GmtGeneSetCollection::new(result)
}

pub fn write_gmx(col : GmtGeneSetCollection) -> String {
  let buf = StringBuilder::new()
  let n = col.sets.length()
  let mut i = 0
  while i < n {
    if i > 0 { buf.write_char('\t') }
    buf.write_string(col.sets[i].name)
    i = i + 1
  }
  buf.write_char('\n')
  let mut i2 = 0
  while i2 < n {
    if i2 > 0 { buf.write_char('\t') }
    buf.write_string(col.sets[i2].description)
    i2 = i2 + 1
  }
  buf.write_char('\n')
  let mut max_len = 0
  for s in col.sets {
    if s.gene_ids.length() > max_len { max_len = s.gene_ids.length() }
  }
  let mut r = 0
  while r < max_len {
    let mut i3 = 0
    while i3 < n {
      if i3 > 0 { buf.write_char('\t') }
      let genes = col.sets[i3].gene_ids
      if r < genes.length() { buf.write_string(genes[r]) }
      i3 = i3 + 1
    }
    buf.write_char('\n')
    r = r + 1
  }
  buf.to_string()
}

// ---- interop helper: build simple name+genes pair record (plain data) -----

pub fn gsea_as_plain_name_genes(gs : GmtGeneSet) -> (String, Array[String]) {
  (gs.name, gs.gene_ids)
}
