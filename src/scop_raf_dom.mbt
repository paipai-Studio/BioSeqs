///|
/// Bio.SCOP.Raf - ASTRAL RAF (Rapid Access Format) Sequence Maps, and
/// Bio.SCOP.Dom - the (deprecated) SCOP DOMain file.
///
/// The ASTRAL RAF Sequence Maps record the relationship between the PDB
/// SEQRES records (representing the sequence of the molecule used in an
/// experiment) and the ATOM records (representing the atoms experimentally
/// observed). Residues are referenced by residue ID: a PDB residue sequence
/// number (up to 4 digits) plus an optional insertion code, e.g. "1",
/// "10A", "-1". RAF uses an underscore for blank chain IDs.
///
/// DOM records consist of 4 tab-delimited fields:
///   sid, pdbid, residues, hierarchy
/// e.g.  d1sctg_  1sct    g:      1.001.001.001.001.001
///
/// Port of Biopython's Bio.SCOP.Raf (SeqMap / SeqMapIndex / Res) and
/// Bio.SCOP.Dom (Record / parse).

///|
/// Error raised when parsing or querying ASTRAL RAF sequence maps.
pub suberror ScopRafError {
  ScopRafError(String)
}

///|
/// Error raised when parsing SCOP DOM records.
pub suberror ScopDomError {
  ScopDomError(String)
}

///|
fn[T] scop_raf_fail(msg : String) -> T raise ScopRafError {
  raise ScopRafError(msg)
}

///|
fn[T] scop_dom_fail(msg : String) -> T raise ScopDomError {
  raise ScopDomError(msg)
}

///|
/// Convert RAF one-letter amino acid codes into IUPAC standard codes.
/// Letters are uppercased, and "." ("Unknown") is converted to "X".
pub fn scop_raf_normalize_letters(one_letter_code : String) -> String {
  if one_letter_code == "." {
    "X"
  } else {
    one_letter_code.to_upper()
  }
}

///|
/// Three-letter to one-letter amino acid map used by RAF's atom extraction.
/// Mirrors Bio.Data.PDBData.protein_letters_3to1_extended: the 20 standard
/// residues plus the ambiguous IUPAC codes (Asx/Glx/Xle/Xaa, Sec/Pyl) and
/// chemically modified residues commonly found in PDB files.
let scop_raf_protein_3to1 : Map[String, String] = Map([
  ("ALA", "A"), ("ARG", "R"), ("ASN", "N"), ("ASP", "D"), ("CYS", "C"),
  ("GLN", "Q"), ("GLU", "E"), ("GLY", "G"), ("HIS", "H"), ("ILE", "I"),
  ("LEU", "L"), ("LYS", "K"), ("MET", "M"), ("PHE", "F"), ("PRO", "P"),
  ("SER", "S"), ("THR", "T"), ("TRP", "W"), ("TYR", "Y"), ("VAL", "V"),
  ("ASX", "B"), ("GLX", "Z"), ("XLE", "J"), ("XAA", "X"), ("SEC", "U"),
  ("PYL", "O"), ("MSE", "M"), ("SEP", "S"), ("TPO", "T"), ("PTR", "Y"),
  ("MLY", "K"), ("ALY", "K"), ("M3L", "K"), ("KCX", "K"), ("LLP", "K"),
  ("CSD", "C"), ("CSO", "C"), ("CME", "C"), ("CAS", "C"), ("CSX", "C"),
  ("OCS", "C"), ("HYP", "P"), ("ALC", "A"), ("FME", "M"), ("ABA", "A"),
  ("PCA", "E"), ("AYA", "A"), ("TRQ", "W"), ("ASX_", "B"),
])

///|
/// A single residue mapping from a RAF record.
pub struct ScopRafRes {
  chainid : String // single character chain ID
  resid : String // residue ID (sequence number + optional insertion code)
  atom : String // one-letter code from ATOM records
  seqres : String // one-letter code from SEQRES records
}

///|
/// An ASTRAL RAF (Rapid Access Format) Sequence Map.
///
/// Attributes:
///   - pdbid: the PDB 4-character ID
///   - pdb_datestamp: from the PDB file
///   - version: the RAF format version, e.g. "0.01"
///   - flags: RAF flags (see release notes)
///   - res: one ScopRafRes per residue in this sequence map
pub struct ScopSeqMap {
  pdbid : String
  pdb_datestamp : String
  version : String
  mut flags : String
  res : Array[ScopRafRes]
}

///|
/// Parse a single RAF record line into a SeqMap.
pub fn scop_raf_seqmap_parse(line : String) -> ScopSeqMap raise ScopRafError {
  let header_len = 38
  // no trailing whitespace (R-style rstrip; hand-rolled because trailing
  // spaces must not shift the fixed-width column layout)
  let mut text = line
  while text.length() > 0 {
    let last = char_at(text, text.length() - 1)
    if last == " " || last == "\t" || last == "\r" {
      text = text[0:(text.length() - 1)].to_owned()
    } else {
      break
    }
  }
  if text.length() < header_len {
    scop_raf_fail("Incomplete header: " + text)
  }
  let pdbid = text[0:4].to_owned()
  let chainid = text[4:5].to_owned()
  let version = text[6:10].to_owned()
  // RAF format versions 0.01 and 0.02 are identical for practical purposes
  if version != "0.01" && version != "0.02" {
    scop_raf_fail("Incompatible RAF version: " + version)
  }
  let datestamp = text[14:20].to_owned()
  let flags = text[21:27].to_owned()
  let res : Array[ScopRafRes] = Array::new()
  let mut i = header_len
  while i < text.length() {
    let f = text[i:(i + 7)].to_owned()
    if f.length() != 7 {
      scop_raf_fail("Corrupt Field: (" + f + ")")
    }
    let resid = f[0:5].trim().to_owned()
    let atom = scop_raf_normalize_letters(f[5:6].to_owned())
    let seqres = scop_raf_normalize_letters(f[6:7].to_owned())
    res.push(ScopRafRes::{ chainid, resid, atom, seqres })
    i = i + 7
  }
  ScopSeqMap::{ pdbid, pdb_datestamp: datestamp, version, flags, res }
}

///|
/// Iterate over a RAF file, giving a SeqMap object for each line.
/// Empty (all whitespace) lines are skipped.
pub fn scop_raf_parse(text : String) -> Array[ScopSeqMap] raise ScopRafError {
  let result : Array[ScopSeqMap] = Array::new()
  for line in split_by_char(text, '\n'.to_int()) {
    if line.trim().length() == 0 {
      continue
    }
    result.push(scop_raf_seqmap_parse(line))
  }
  result
}

///|
/// Build a RAF file index: PDB ID + chain ID (e.g. "2drcA", "155c_") mapped
/// to the raw record line. Underscore indicates blank chain IDs.
pub fn scop_raf_index(text : String) -> Map[String, String] {
  let index : Map[String, String] = Map([], capacity=16)
  for line in split_by_char(text, '\n'.to_int()) {
    if line.trim().length() == 0 {
      continue
    }
    if line.length() < 5 {
      continue
    }
    let key = line[0:5].to_owned()
    index.set(key, line)
  }
  index
}

///|
/// Number of residues in this sequence map.
pub fn ScopSeqMap::length(self : ScopSeqMap) -> Int {
  self.res.length()
}

///|
/// Get the residue at position i.
pub fn ScopSeqMap::get(self : ScopSeqMap, i : Int) -> ScopRafRes {
  self.res[i]
}

///|
/// All residues (shared array reference).
pub fn ScopSeqMap::residues(self : ScopSeqMap) -> Array[ScopRafRes] {
  self.res
}

///|
/// Return the index of the SeqMap for the given resid and chainid.
pub fn ScopSeqMap::index(
  self : ScopSeqMap,
  resid : String,
  chainid? : String,
) -> Int raise ScopRafError {
  let cid = chainid.unwrap_or("_")
  for i in 0..<self.res.length() {
    if self.res[i].resid == resid && self.res[i].chainid == cid {
      return i
    }
  }
  scop_raf_fail("No such residue " + cid + resid)
}

///|
/// Slice this SeqMap into a fragment [start, end).
pub fn ScopSeqMap::slice(
  self : ScopSeqMap,
  start : Int,
  end_ : Int,
) -> ScopSeqMap {
  let res : Array[ScopRafRes] = Array::new()
  for i in start..<end_ {
    res.push(self.res[i])
  }
  ScopSeqMap::{
    pdbid: self.pdbid,
    pdb_datestamp: self.pdb_datestamp,
    version: self.version,
    flags: self.flags,
    res,
  }
}

///|
/// Concatenate two SeqMaps into a new one. Both must have the same PDB ID,
/// PDB datestamp and RAF version. The RAF flags are erased if they are
/// inconsistent (this may happen when fragments are taken from different
/// chains).
pub fn ScopSeqMap::add(
  self : ScopSeqMap,
  other : ScopSeqMap,
) -> ScopSeqMap raise ScopRafError {
  let out = self.slice(0, self.res.length())
  self.extend(other)
  out
}

///|
/// Append another SeqMap onto the end of self (in place).
pub fn ScopSeqMap::extend(
  self : ScopSeqMap,
  other : ScopSeqMap,
) -> Unit raise ScopRafError {
  if self.pdbid != other.pdbid {
    scop_raf_fail("Cannot add fragments from different proteins")
  }
  if self.version != other.version {
    scop_raf_fail("Incompatible rafs")
  }
  if self.pdb_datestamp != other.pdb_datestamp {
    scop_raf_fail("Different pdb dates!")
  }
  if self.flags != other.flags {
    self.flags = ""
  }
  for r in other.res {
    self.res.push(r)
  }
}

///|
/// Get the sequence map for a collection of residues.
///
/// - index: RAF index as built by scop_raf_index (key = pdbid + chainid)
/// - residues: a ScopResidues instance. If it has no fragments, all
///   residues of the unnamed chain are used.
pub fn scop_raf_get_seqmap(
  index : Map[String, String],
  residues : ScopResidues,
) -> ScopSeqMap raise ScopRafError {
  let pdbid = residues.pdbid
  let mut frags = residues.fragments
  if frags.length() == 0 {
    frags = [ScopFragment::new("_", "", "")]
  }
  let mut seqmap : ScopSeqMap? = None
  for frag in frags {
    let mut chainid = frag.chain
    if chainid == "" || chainid == "-" || chainid == " " || chainid == "_" {
      chainid = "_"
    }
    let id = pdbid + chainid
    let line = match index.get(id) {
      Some(l) => l
      None => scop_raf_fail("RAF record not found for " + id)
    }
    let mut sm = scop_raf_seqmap_parse(line)
    // Cut out fragment of interest
    let start = if frag.start.length() > 0 {
      sm.index(frag.start, chainid~)
    } else {
      0
    }
    let end_ = if frag.end_.length() > 0 {
      sm.index(frag.end_, chainid~) + 1
    } else {
      sm.res.length()
    }
    sm = sm.slice(start, end_)
    match seqmap {
      Some(s) => s.extend(sm)
      None => seqmap = Some(sm)
    }
  }
  match seqmap {
    Some(s) => s
    None => scop_raf_fail("empty residue specification")
  }
}

///|
/// Extract all relevant ATOM and HETATOM records from a PDB file.
///
/// The PDB file is scanned for ATOM and HETATOM records. If the chain ID,
/// residue ID (seqNum and iCode), and residue type match a residue in this
/// sequence map, then the record is echoed to the output. Typically used to
/// find the coordinates of a domain or other residue subset.
pub fn scop_raf_get_atoms(
  seqmap : ScopSeqMap,
  pdb_content : String,
) -> String raise ScopRafError {
  // The set of residues that we have to find records for.
  let res_set : Map[(String, String), ScopRafRes] = Map([], capacity=16)
  for r in seqmap.res {
    if r.atom == "X" {
      // Unknown residue type
      continue
    }
    let mut chainid = r.chainid
    if chainid == "_" {
      chainid = " "
    }
    res_set.set((chainid, r.resid), r)
  }
  let found : Map[(String, String), ScopRafRes] = Map([], capacity=16)
  let buf = StringBuilder::new()
  for line in split_by_char(pdb_content, '\n'.to_int()) {
    if line.length() < 27 {
      continue
    }
    let record = line[0:6].to_owned()
    if record == "ATOM  " || record == "HETATM" {
      let chainid = line[21:22].to_owned()
      let resid = line[22:27].trim().to_owned()
      match res_set.get((chainid, resid)) {
        Some(r) => {
          let res_name = line[17:20].to_owned()
          match scop_raf_protein_3to1.get(res_name.to_upper()) {
            Some(aa) =>
              if aa == r.atom {
                buf.write_string(line)
                buf.write_char('\n')
                found.set((chainid, resid), r)
              }
            None => ()
          }
        }
        None => ()
      }
    }
  }
  if res_set.length() != found.length() {
    // Some residue(s) missing: report the first one for clarity.
    let entries = res_set.iter().collect()
    for entry in entries {
      let (key, _) = entry
      if !found.contains(key) {
        scop_raf_fail(
          "Could not find at least one ATOM or HETATM record for each and every residue in this sequence map. (first missing: chain '\{key.0}' resid '\{key.1}')",
        )
      }
    }
    scop_raf_fail(
      "Could not find at least one ATOM or HETATM record for each and every residue in this sequence map.",
    )
  }
  buf.to_string()
}

// ---------------------------------------------------------------------------
// Bio.SCOP.Dom
// ---------------------------------------------------------------------------

///|
/// Holds information for one SCOP domain.
pub struct DomRecord {
  sid : String // the SCOP ID of the entry, e.g. d1anu1
  residues : ScopResidues // the domain definition
  hierarchy : String // where this domain is in the hierarchy
}

///|
/// Parse a single DOM record line.
///
/// Records consist of 4 tab delimited fields:
/// sid, pdbid, residues, hierarchy
pub fn scop_dom_record_parse(line : String) -> DomRecord raise ScopDomError {
  let mut text = line
  while text.length() > 0 {
    let last = char_at(text, text.length() - 1)
    if last == " " || last == "\t" || last == "\r" || last == "\n" {
      text = text[0:(text.length() - 1)].to_owned()
    } else {
      break
    }
  }
  if text.length() == 0 {
    scop_dom_fail("I don't understand the format of " + line)
  }
  let columns = split_by_char(text, '\t'.to_int())
  if columns.length() != 4 {
    scop_dom_fail("I don't understand the format of " + line)
  }
  let sid = columns[0]
  let pdbid = columns[1]
  let res = columns[2]
  let hierarchy = columns[3]
  let parsed = scop_parse_residues(res)
  // the residues column has no pdbid prefix; set it explicitly
  let residues = ScopResidues::new(pdbid, parsed.fragments)
  DomRecord::{ sid, residues, hierarchy }
}

///|
/// Iterate over a DOM file as a Dom record for each line.
/// Lines starting with "#" are comments and are skipped.
pub fn scop_dom_parse(text : String) -> Array[DomRecord] raise ScopDomError {
  let result : Array[DomRecord] = Array::new()
  for line in split_by_char(text, '\n'.to_int()) {
    if line.length() == 0 {
      continue
    }
    if line[0:1].to_string() == "#" {
      continue
    }
    if line.trim().length() == 0 {
      continue
    }
    result.push(scop_dom_record_parse(line))
  }
  result
}

///|
/// Represent the SCOP domain record as a tab-separated string.
pub fn DomRecord::to_string(self : DomRecord) -> String {
  let res_str = scop_residues_to_string(self.residues)
  let res_str = if res_str == "-" {
    "-"
  } else {
    // replace separating spaces with tabs (sid \t pdbid \t residues)
    let buf = StringBuilder::new()
    for ch in res_str {
      if ch == ' ' {
        buf.write_char('\t')
      } else {
        buf.write_char(ch)
      }
    }
    buf.to_string()
  }
  self.sid + "\t" + res_str + "\t" + self.hierarchy + "\n"
}
