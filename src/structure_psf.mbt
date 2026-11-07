///|
/// CHARMM/XPLOR/NAMD PSF (Protein Structure File) topology reader/writer.
///
/// A PSF file is organized in sections, each introduced by a header line
/// ` N !KEYWORD` where N is the number of entries that follow. Supported
/// keywords: NTITLE, NATOM, NBOND, NTHETA, NPHI, NIMPHI, NDON, NACC, NNB.
/// EXT (extended) mode is assumed: every field is whitespace separated, so
/// variable-width numbers (EXT style) parse correctly with a plain split.

///|
/// One atom record of a PSF topology file (9 NATOM columns minus imove).
pub(all) struct PsfAtom {
  atomid : Int32
  segid : String
  resid : Int32
  resname : String
  atom_name : String
  atom_type : String
  charge : Double
  mass : Double
}

///|
/// Parsed PSF topology: atom records plus connectivity tables.
pub struct PsfTopology {
  atoms : Array[PsfAtom]
  bonds : Array[(Int32, Int32)]
  angles : Array[(Int32, Int32, Int32)]
  dihedrals : Array[(Int32, Int32, Int32, Int32)]
  impropers : Array[(Int32, Int32, Int32, Int32)]
  donors : Array[(Int32, Int32)]
  acceptors : Array[(Int32, Int32)]
  titles : Array[String]
}

///|
/// Error raised for malformed PSF input.
pub suberror PsfError {
  PsfError(String)
}

///|
/// Default constructor for PsfTopology.
pub fn PsfTopology::new(
  atoms : Array[PsfAtom],
  bonds : Array[(Int32, Int32)],
  angles : Array[(Int32, Int32, Int32)],
  dihedrals : Array[(Int32, Int32, Int32, Int32)],
  impropers : Array[(Int32, Int32, Int32, Int32)],
  donors : Array[(Int32, Int32)],
  acceptors : Array[(Int32, Int32)],
  titles : Array[String],
) -> PsfTopology {
  PsfTopology::{
    atoms,
    bonds,
    angles,
    dihedrals,
    impropers,
    donors,
    acceptors,
    titles,
  }
}

///|
/// Number of atoms in the topology.
pub fn PsfTopology::natoms(self : PsfTopology) -> Int32 {
  self.atoms.length()
}

///|
/// Number of bonds.
pub fn PsfTopology::n_bonds(self : PsfTopology) -> Int32 {
  self.bonds.length()
}

///|
/// Number of angles.
pub fn PsfTopology::n_angles(self : PsfTopology) -> Int32 {
  self.angles.length()
}

///|
/// Number of dihedrals.
pub fn PsfTopology::n_dihedrals(self : PsfTopology) -> Int32 {
  self.dihedrals.length()
}

///|
/// Number of improper dihedrals.
pub fn PsfTopology::n_impropers(self : PsfTopology) -> Int32 {
  self.impropers.length()
}

///|
/// Sum of all atomic partial charges.
pub fn PsfTopology::total_charge(self : PsfTopology) -> Double {
  let mut total = 0.0
  for a in self.atoms {
    total = total + a.charge
  }
  total
}

///|
/// Split a PSF line into whitespace-separated tokens.
fn psf_tokens(line : StringView) -> Array[String] {
  let tokens : Array[String] = Array(capacity=0)
  let parts = line.split(" ")
  for p in parts {
    let t = p.trim()
    if t.length() > 0 {
      tokens.push(t.to_owned())
    }
  }
  tokens
}

///|
/// Check whether a token starts with '!' (section-keyword marker).
fn psf_is_bang_token(tok : String) -> Bool {
  tok.length() > 0 && tok.unsafe_get(0).to_int() == 33
}

///|
/// Strip the leading '!' and the trailing ": label" from a section keyword,
/// e.g. "!NBOND: bonds" -> "NBOND".
fn psf_clean_keyword(tok : String) -> String {
  let n = tok.length()
  let mut start = 0
  if n > 0 && tok.unsafe_get(0).to_int() == 33 {
    start = 1
  }
  let mut end = n
  let mut j = start
  while j < n {
    if tok.unsafe_get(j).to_int() == 58 {
      end = j
      break
    }
    j = j + 1
  }
  tok[start:end].to_owned()
}

///|
/// Parse one NATOM row: atomid segid resid resname atomname atomtype charge
/// mass [imove]. Raises when fewer than 8 usable columns are present.
fn psf_parse_atom_line(line : StringView) -> PsfAtom raise PsfError {
  let tokens = psf_tokens(line)
  if tokens.length() < 8 {
    raise PsfError(
      "PSF NATOM entry needs 9 columns, got \{tokens.length()}: \{line.to_owned()}",
    )
  }
  PsfAtom::{
    atomid: parse_int(tokens[0]),
    segid: tokens[1],
    resid: parse_int(tokens[2]),
    resname: tokens[3],
    atom_name: tokens[4],
    atom_type: tokens[5],
    charge: psf_opt_double(tokens[6]),
    mass: psf_opt_double(tokens[7]),
  }
}

///|
/// Parse a double token, defaulting to 0.0 when missing or malformed.
fn psf_opt_double(s : String) -> Double {
  match parse_double(s) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Check whether a token is a plain (optionally signed) integer.
fn psf_is_int_token(s : String) -> Bool {
  let n = s.length()
  if n == 0 {
    return false
  }
  let mut i = 0
  if s.unsafe_get(0).to_int() == 45 && n > 1 {
    i = 1
  }
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < 48 || c > 57 {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Collect up to `total` integer tokens from consecutive all-numeric lines
/// starting at index `start`. Blank lines are skipped; the first line that
/// contains a non-numeric token (e.g. the next section header) stops the
/// stream. Returns the numbers and the index of the first unconsumed line.
fn psf_collect_numbers(
  lines : Array[StringView],
  start : Int32,
  total : Int32,
) -> (Array[Int32], Int32) {
  let nums : Array[Int32] = Array(capacity=total)
  let n = lines.length()
  let mut i = start
  while i < n && nums.length() < total {
    let t = lines[i].trim()
    if t.length() == 0 {
      i = i + 1
      continue
    }
    let tokens = psf_tokens(t)
    let mut line_ok = true
    for tok in tokens {
      if !psf_is_int_token(tok) {
        line_ok = false
        break
      }
    }
    if !line_ok {
      break
    }
    for tok in tokens {
      if nums.length() < total {
        nums.push(parse_int(tok))
      }
    }
    i = i + 1
  }
  (nums, i)
}

///|
/// Skip `count` non-empty entry lines (used for unknown sections and NNB).
fn psf_skip_entries(
  lines : Array[StringView],
  start : Int32,
  count : Int32,
) -> Int32 {
  let n = lines.length()
  let mut i = start
  let mut got = 0
  while got < count && i < n {
    let t = lines[i].trim()
    i = i + 1
    if t.length() == 0 {
      continue
    }
    got = got + 1
  }
  i
}

///|
/// Parse a CHARMM/XPLOR/NAMD PSF topology from its text form.
///
/// Sections are located by their ` N !KEYWORD` header line. Numeric sections
/// (NBOND/NTHETA/NPHI/NIMPHI/NDON/NACC) are parsed as a flat stream of
/// integers split into groups of 2/3/4/4/2/2 atom ids; several groups may
/// share one physical line. Unknown sections (and NNB) skip exactly their
/// declared entry count.
pub fn psf_parse(text : String) -> PsfTopology raise PsfError {
  let atoms : Array[PsfAtom] = Array(capacity=0)
  let bonds : Array[(Int32, Int32)] = Array(capacity=0)
  let angles : Array[(Int32, Int32, Int32)] = Array(capacity=0)
  let dihedrals : Array[(Int32, Int32, Int32, Int32)] = Array(capacity=0)
  let impropers : Array[(Int32, Int32, Int32, Int32)] = Array(capacity=0)
  let donors : Array[(Int32, Int32)] = Array(capacity=0)
  let acceptors : Array[(Int32, Int32)] = Array(capacity=0)
  let titles : Array[String] = Array(capacity=0)
  let lines = text.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  while i < n {
    let trimmed = lines[i].trim()
    i = i + 1
    if trimmed.length() == 0 {
      continue
    }
    let tokens = psf_tokens(trimmed)
    if tokens.length() >= 2 && psf_is_bang_token(tokens[1]) {
      let count = parse_int(tokens[0])
      let kw = psf_clean_keyword(tokens[1])
      match kw {
        "NTITLE" => {
          let mut got = 0
          while got < count && i < n {
            let tl = lines[i].trim()
            i = i + 1
            if tl.length() == 0 {
              continue
            }
            titles.push(tl.to_owned())
            got = got + 1
          }
        }
        "NATOM" => {
          let mut got = 0
          while got < count && i < n {
            let tl = lines[i].trim()
            i = i + 1
            if tl.length() == 0 {
              continue
            }
            atoms.push(psf_parse_atom_line(tl))
            got = got + 1
          }
        }
        "NBOND" => {
          let (nums, ni) = psf_collect_numbers(lines, i, count * 2)
          i = ni
          let mut idx = 0
          while idx + 2 <= nums.length() {
            bonds.push((nums[idx], nums[idx + 1]))
            idx = idx + 2
          }
        }
        "NTHETA" => {
          let (nums, ni) = psf_collect_numbers(lines, i, count * 3)
          i = ni
          let mut idx = 0
          while idx + 3 <= nums.length() {
            angles.push((nums[idx], nums[idx + 1], nums[idx + 2]))
            idx = idx + 3
          }
        }
        "NPHI" => {
          let (nums, ni) = psf_collect_numbers(lines, i, count * 4)
          i = ni
          let mut idx = 0
          while idx + 4 <= nums.length() {
            dihedrals.push(
              (nums[idx], nums[idx + 1], nums[idx + 2], nums[idx + 3]),
            )
            idx = idx + 4
          }
        }
        "NIMPHI" => {
          let (nums, ni) = psf_collect_numbers(lines, i, count * 4)
          i = ni
          let mut idx = 0
          while idx + 4 <= nums.length() {
            impropers.push(
              (nums[idx], nums[idx + 1], nums[idx + 2], nums[idx + 3]),
            )
            idx = idx + 4
          }
        }
        "NDON" => {
          let (nums, ni) = psf_collect_numbers(lines, i, count * 2)
          i = ni
          let mut idx = 0
          while idx + 2 <= nums.length() {
            donors.push((nums[idx], nums[idx + 1]))
            idx = idx + 2
          }
        }
        "NACC" => {
          let (nums, ni) = psf_collect_numbers(lines, i, count * 2)
          i = ni
          let mut idx = 0
          while idx + 2 <= nums.length() {
            acceptors.push((nums[idx], nums[idx + 1]))
            idx = idx + 2
          }
        }
        _ =>
          // Unknown section (including NNB): skip its declared entry count.
          i = psf_skip_entries(lines, i, count)
      }
    }
  }
  // Non-header lines outside any section (e.g. "PSF EXT") are ignored.
  PsfTopology::{
    atoms,
    bonds,
    angles,
    dihedrals,
    impropers,
    donors,
    acceptors,
    titles,
  }
}

///|
/// Format one section header line: " N !KEYWORD".
fn psf_section_header(count : Int32, keyword : String) -> String {
  pad_left(count.to_string(), 10, ' ') + " !" + keyword + "\n"
}

///|
/// Format one NATOM row in the standard fixed-width PSF layout
/// (8d 1x 4s 1x 8d 1x 4s 1x 4s 1x 4s 14.6f 14.4f 8d, imove = 0).
fn psf_atom_line_string(a : PsfAtom) -> String {
  let buf = StringBuilder()
  buf.write_string(pad_left(a.atomid.to_string(), 8, ' '))
  buf.write_char(' ')
  buf.write_string(pad_right(a.segid, 4, ' '))
  buf.write_char(' ')
  buf.write_string(pad_left(a.resid.to_string(), 8, ' '))
  buf.write_char(' ')
  buf.write_string(pad_right(a.resname, 4, ' '))
  buf.write_char(' ')
  buf.write_string(pad_right(a.atom_name, 4, ' '))
  buf.write_char(' ')
  buf.write_string(pad_right(a.atom_type, 4, ' '))
  buf.write_string(format_double(a.charge, 14, 6))
  buf.write_string(format_double(a.mass, 14, 4))
  buf.write_string(pad_left("0", 8, ' '))
  buf.write_char('\n')
  buf.to_string()
}

///|
/// Write 2-id groups, `per_line` groups per physical line.
fn psf_write_pairs(
  buf : StringBuilder,
  pairs : Array[(Int32, Int32)],
  per_line : Int32,
) -> Unit {
  let mut in_line = 0
  for p in pairs {
    buf.write_string(pad_left(p.0.to_string(), 8, ' '))
    buf.write_string(pad_left(p.1.to_string(), 8, ' '))
    in_line = in_line + 1
    if in_line >= per_line {
      buf.write_char('\n')
      in_line = 0
    }
  }
  if in_line > 0 {
    buf.write_char('\n')
  }
}

///|
/// Write 3-id groups, `per_line` groups per physical line.
fn psf_write_triples(
  buf : StringBuilder,
  triples : Array[(Int32, Int32, Int32)],
  per_line : Int32,
) -> Unit {
  let mut in_line = 0
  for p in triples {
    buf.write_string(pad_left(p.0.to_string(), 8, ' '))
    buf.write_string(pad_left(p.1.to_string(), 8, ' '))
    buf.write_string(pad_left(p.2.to_string(), 8, ' '))
    in_line = in_line + 1
    if in_line >= per_line {
      buf.write_char('\n')
      in_line = 0
    }
  }
  if in_line > 0 {
    buf.write_char('\n')
  }
}

///|
/// Write 4-id groups, `per_line` groups per physical line.
fn psf_write_quads(
  buf : StringBuilder,
  quads : Array[(Int32, Int32, Int32, Int32)],
  per_line : Int32,
) -> Unit {
  let mut in_line = 0
  for p in quads {
    buf.write_string(pad_left(p.0.to_string(), 8, ' '))
    buf.write_string(pad_left(p.1.to_string(), 8, ' '))
    buf.write_string(pad_left(p.2.to_string(), 8, ' '))
    buf.write_string(pad_left(p.3.to_string(), 8, ' '))
    in_line = in_line + 1
    if in_line >= per_line {
      buf.write_char('\n')
      in_line = 0
    }
  }
  if in_line > 0 {
    buf.write_char('\n')
  }
}

///|
/// Serialize a topology to standard PSF text (EXT mode, fixed-width NATOM
/// rows, 4/2/1 groups per line for bonds/angles/dihedrals respectively).
/// NNB data is not stored in PsfTopology and therefore not written.
pub fn psf_write(topo : PsfTopology) -> String {
  let buf = StringBuilder()
  buf.write_string("PSF EXT\n\n")
  buf.write_string(psf_section_header(topo.titles.length(), "NTITLE"))
  for t in topo.titles {
    // Titles parsed from a PSF file keep their original " REMARKS ..." text,
    // so a single leading space reproduces the input verbatim.
    buf.write_string(" ")
    buf.write_string(t)
    buf.write_char('\n')
  }
  buf.write_char('\n')
  buf.write_string(psf_section_header(topo.atoms.length(), "NATOM"))
  for a in topo.atoms {
    buf.write_string(psf_atom_line_string(a))
  }
  buf.write_char('\n')
  buf.write_string(psf_section_header(topo.bonds.length(), "NBOND: bonds"))
  psf_write_pairs(buf, topo.bonds, 4)
  buf.write_char('\n')
  buf.write_string(psf_section_header(topo.angles.length(), "NTHETA: angles"))
  psf_write_triples(buf, topo.angles, 2)
  buf.write_char('\n')
  buf.write_string(
    psf_section_header(topo.dihedrals.length(), "NPHI: dihedrals"),
  )
  psf_write_quads(buf, topo.dihedrals, 1)
  buf.write_char('\n')
  buf.write_string(
    psf_section_header(topo.impropers.length(), "NIMPHI: impropers"),
  )
  psf_write_quads(buf, topo.impropers, 1)
  buf.write_char('\n')
  buf.write_string(psf_section_header(topo.donors.length(), "NDON: donors"))
  psf_write_pairs(buf, topo.donors, 4)
  buf.write_char('\n')
  buf.write_string(
    psf_section_header(topo.acceptors.length(), "NACC: acceptors"),
  )
  psf_write_pairs(buf, topo.acceptors, 4)
  buf.write_char('\n')
  buf.to_string()
}

///|
/// Whether a residue name is one of the 20 standard amino acids.
pub fn psf_is_protein_resname(resname : String) -> Bool {
  match resname {
    "ALA"
    | "ARG"
    | "ASN"
    | "ASP"
    | "CYS"
    | "GLN"
    | "GLU"
    | "GLY"
    | "HIS"
    | "ILE"
    | "LEU"
    | "LYS"
    | "MET"
    | "PHE"
    | "PRO"
    | "SER"
    | "THR"
    | "TRP"
    | "TYR"
    | "VAL" => true
    _ => false
  }
}
