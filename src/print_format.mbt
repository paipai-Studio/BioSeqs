///|
/// PrintFormat - print restriction enzyme analysis results.
///
/// Port of Biopython's Bio.Restriction.PrintFormat. Takes the results of a
/// `RestrictionBatch::search` (a map of enzyme name -> 1-based cutting
/// positions) and renders them in three layouts: "list" (alphabetical),
/// "number" (grouped by number of cuts) and "map" (a graphical sequence map
/// with the cut positions and both strands).

///|
pub struct PrintFormatConfig {
  console_width : Int
  name_width : Int
  max_size : Int
  indent : Int
}

///|
/// The default layout, mirroring Biopython's class attributes.
pub fn print_format_config_default() -> PrintFormatConfig {
  PrintFormatConfig::{ console_width: 80, name_width: 10, max_size: 6, indent: 4 }
}

///|
fn pf_linesize(cfg : PrintFormatConfig) -> Int {
  (cfg.console_width - cfg.console_width % cfg.name_width) - cfg.name_width
}

///|
fn pf_repeat(s : String, n : Int) -> String {
  let buf = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    buf.write_string(s)
  }
  buf.to_string()
}

///|
/// Python str.ljust: left-justify to `width` with spaces.
fn pf_ljust(s : String, width : Int) -> String {
  if s.length() >= width {
    s
  } else {
    s + pf_repeat(" ", width - s.length())
  }
}

///|
/// Python str.rjust: right-justify to `width` with spaces.
fn pf_rjust(s : String, width : Int) -> String {
  if s.length() >= width {
    s
  } else {
    pf_repeat(" ", width - s.length()) + s
  }
}

///|
fn pf_join(arr : Array[String], sep : String) -> String {
  let buf = StringBuilder()
  for i = 0; i < arr.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(sep)
    }
    buf.write_string(arr[i])
  }
  buf.to_string()
}

///|
/// The enzyme-name suffix for a cut position: a leading space then the names
/// joined by single spaces (matches Python's " ".join((s, n))).
fn pf_names_suffix(names : Array[String]) -> String {
  let buf = StringBuilder()
  for name in names {
    buf.write_char(' ')
    buf.write_string(name)
  }
  buf.to_string()
}

///|
/// DNA complement (IUPAC) for the reverse strand of the map.
fn pf_complement(seq : String) -> String {
  let buf = StringBuilder()
  for i = 0; i < seq.length(); i = i + 1 {
    match seq.unsafe_get(i).to_int() {
      65 => buf.write_char('T')
      84 => buf.write_char('A')
      67 => buf.write_char('G')
      71 => buf.write_char('C')
      82 => buf.write_char('Y')
      89 => buf.write_char('R')
      83 => buf.write_char('S')
      87 => buf.write_char('W')
      75 => buf.write_char('M')
      77 => buf.write_char('K')
      66 => buf.write_char('V')
      86 => buf.write_char('B')
      68 => buf.write_char('H')
      72 => buf.write_char('D')
      78 => buf.write_char('N')
      97 => buf.write_char('t')
      116 => buf.write_char('a')
      99 => buf.write_char('g')
      103 => buf.write_char('c')
      _ => buf.write_char(seq.unsafe_get(i).unsafe_to_char())
    }
  }
  buf.to_string()
}

///|
/// Split the results into cutting (non-empty) and non-cutting enzymes.
fn pf_split(
  dct : Map[String, Array[Int]],
) -> (Array[(String, Array[Int])], Array[String]) {
  let ls : Array[(String, Array[Int])] = []
  let nc : Array[String] = []
  for name in dct.keys() {
    match dct.get(name) {
      Some(sites) => if sites.length() > 0 { ls.push((name, sites)) } else { nc.push(name) }
      None => ()
    }
  }
  (ls, nc)
}

///|
/// Wrap a long "pos1, pos2, ... ." output into chunks of at most `linesize`
/// characters that end with a comma or period.
fn pf_wrap(output : String, linesize : Int) -> Array[String] {
  let chunks : Array[String] = []
  let n = output.length()
  let mut start = 0
  while start < n {
    let mut end = start + linesize
    if end > n {
      end = n
    }
    let mut cut = -1
    let mut j = end - 1
    while j >= start {
      let c = output.unsafe_get(j).to_int()
      if c == 44 || c == 46 {
        cut = j
        break
      }
      j = j - 1
    }
    if cut < 0 {
      cut = end - 1
    }
    chunks.push(output[start:cut + 1].to_owned())
    start = cut + 1
  }
  chunks
}

///|
/// Render one section of (name, sites) pairs in the "list" style.
fn pf_next_section(
  ls : Array[(String, Array[Int])],
  into : String,
  cfg : PrintFormatConfig,
) -> String {
  let indentation = "\n" + pf_repeat(" ", cfg.name_width + cfg.indent)
  let linesize = pf_linesize(cfg) - cfg.max_size
  let sorted = ls.copy()
  sorted.sort_by(fn(a : (String, Array[Int]), b : (String, Array[Int])) -> Int {
    a.0.compare(b.0)
  })
  let mut out = into
  for (name, sites) in sorted {
    let buf = StringBuilder()
    for j = 0; j < sites.length(); j = j + 1 {
      if j > 0 {
        buf.write_string(", ")
      }
      buf.write_string(sites[j].to_string())
    }
    buf.write_char('.')
    let output = buf.to_string()
    let stringsite = if output.length() > linesize {
      pf_join(pf_wrap(output, linesize), indentation)
    } else {
      output
    }
    out = out + pf_ljust(name, cfg.name_width) + " :  " + stringsite + "\n"
  }
  out
}

///|
/// Render the non-cutting enzymes section.
fn pf_make_nocut_only(
  nc : Array[String],
  s1 : String,
  cfg : PrintFormatConfig,
) -> String {
  if nc.length() == 0 {
    return s1
  }
  let sorted = nc.copy()
  sorted.sort_by(fn(a : String, b : String) -> Int { a.compare(b) })
  let linesize = pf_linesize(cfg)
  let mut stringsite = if s1.length() > 0 {
    s1
  } else {
    "\n   Enzymes which do not cut the sequence.\n\n"
  }
  let mut st = ""
  for key in sorted {
    st = st + pf_ljust(key, cfg.name_width)
    if st.length() > linesize {
      stringsite = stringsite + st + "\n"
      st = ""
    }
  }
  stringsite + st + "\n"
}

///|
fn pf_make_list_only(
  ls : Array[(String, Array[Int])],
  title : String,
  cfg : PrintFormatConfig,
) -> String {
  if ls.length() == 0 {
    title
  } else {
    pf_next_section(ls, title, cfg)
  }
}

///|
fn pf_make_number_only(
  ls : Array[(String, Array[Int])],
  title : String,
  cfg : PrintFormatConfig,
) -> String {
  if ls.length() == 0 {
    return title
  }
  let sorted = ls.copy()
  sorted.sort_by(fn(a : (String, Array[Int]), b : (String, Array[Int])) -> Int {
    a.1.length().compare(b.1.length())
  })
  let mut new_sect : Array[(String, Array[Int])] = []
  let mut cur_len = 1
  let mut out = title
  let mut i = 0
  while i < sorted.length() {
    let (name, sites) = sorted[i]
    let length = sites.length()
    if length > cur_len {
      out = out + "\n\nenzymes which cut \{cur_len} times :\n\n"
      out = pf_next_section(new_sect, out, cfg)
      new_sect = [(name, sites)]
      cur_len = length
      i = i + 1
      continue
    }
    new_sect.push((name, sites))
    i = i + 1
  }
  out = out + "\n\nenzymes which cut \{cur_len} times :\n\n"
  pf_next_section(new_sect, out, cfg)
}

///|
/// Render the graphical sequence map (Port of _make_map_only).
fn pf_make_map_only(
  ls : Array[(String, Array[Int])],
  title : String,
  sequence : String,
  cfg : PrintFormatConfig,
) -> String {
  if ls.length() == 0 {
    return title
  }
  // position -> enzyme names
  let enzymemap : Map[Int, Array[String]] = Map([], capacity=ls.length() + 1)
  for (enzyme, cut) in ls {
    for c in cut {
      match enzymemap.get(c) {
        Some(arr) => {
          arr.push(enzyme)
          enzymemap.set(c, arr)
        }
        None => {
          let arr : Array[String] = []
          arr.push(enzyme)
          enzymemap.set(c, arr)
        }
      }
    }
  }
  let mapping : Array[Int] = []
  for pos in enzymemap.keys() {
    mapping.push(pos)
  }
  mapping.sort_by(fn(a : Int, b : Int) -> Int { a.compare(b) })

  let length = sequence.length()
  // partition positions into 60 bp windows; x ends up as the last window end
  let cutloc : Map[Int, Array[Int]] = Map([], capacity=8)
  let mut remaining = mapping.copy()
  let mut x = 0
  let mut next = 60
  while next < length {
    let counter = next - 60
    let loc : Array[Int] = []
    let rem : Array[Int] = []
    for key in remaining {
      if key <= next {
        loc.push(key)
      } else {
        rem.push(key)
      }
    }
    cutloc.set(counter, loc)
    remaining = rem
    x = next
    next = next + 60
  }
  cutloc.set(x, remaining)

  let revsequence = pf_complement(sequence)
  let a = "|"
  let emptyline = pf_repeat(" ", 60)
  let mut map = title

  // full 60 bp windows
  let mut win = 60
  while win < length {
    let mut line = emptyline
    for key in cutloc.get(win - 60).unwrap_or([]) {
      let s = pf_names_suffix(enzymemap.get(key).unwrap_or([]))
      if key == win {
        let chunk = line[0:59].to_owned()
        map = map + chunk + key.to_string() + s + "\n" + chunk + a + "\n"
        break
      }
      let k = key % 60
      let lineo = line[0:k - 1].to_owned() + key.to_string() + s + "\n"
      line = line[0:k - 1].to_owned() + a + line[k:line.length()].to_owned()
      let line2 = line[0:k - 1].to_owned() + a + line[k:line.length()].to_owned() + "\n"
      map = map + lineo + line2
    }
    let mapunit = pf_join(
      [
        sequence[win - 60:win].to_owned(),
        pf_repeat(a, 60),
        revsequence[win - 60:win].to_owned(),
        pf_ljust((win - 59).to_string(), 15) + pf_repeat(" ", 30) +
          pf_rjust(win.to_string(), 15) + "\n\n",
      ],
      "\n",
    )
    map = map + mapunit
    win = win + 60
  }

  // final (partial) window, from the last window end x
  let mut line = emptyline
  for key in cutloc.get(x).unwrap_or([]) {
    let s = pf_names_suffix(enzymemap.get(key).unwrap_or([]))
    if key == length {
      let chunk = line[0:length - 1].to_owned()
      map = map + chunk + key.to_string() + s + "\n" + chunk + a + "\n"
      break
    }
    let k = key % 60
    let lineo = line[0:k - 1].to_owned() + key.to_string() + s + "\n"
    line = line[0:k - 1].to_owned() + a + line[k:line.length()].to_owned()
    let line2 = line[0:k - 1].to_owned() + a + line[k:line.length()].to_owned() + "\n"
    map = map + lineo + line2
  }
  let mapunit = sequence[x:length].to_owned() + "\n" +
    pf_repeat(a, length - x) + "\n" +
    revsequence[x:length].to_owned() + "\n" +
    pf_ljust((x + 1).to_string(), 15) + pf_repeat(" ", length - x - 30) +
    pf_rjust(length.to_string(), 15) + "\n\n"
  map + mapunit
}

///|
/// Format the results as a list.
pub fn print_format_list(
  dct~ : Map[String, Array[Int]],
  title~ : String,
  s1~ : String,
  config? : PrintFormatConfig,
) -> String {
  let cfg = config.unwrap_or(print_format_config_default())
  let (ls, nc) = pf_split(dct)
  pf_make_list_only(ls, title, cfg) + pf_make_nocut_only(nc, s1, cfg)
}

///|
/// Format the results grouped by number of cutting sites.
pub fn print_format_number(
  dct~ : Map[String, Array[Int]],
  title~ : String,
  s1~ : String,
  config? : PrintFormatConfig,
) -> String {
  let cfg = config.unwrap_or(print_format_config_default())
  let (ls, nc) = pf_split(dct)
  pf_make_number_only(ls, title, cfg) + pf_make_nocut_only(nc, s1, cfg)
}

///|
/// Format the results as a graphical map of `sequence`.
pub fn print_format_map(
  dct~ : Map[String, Array[Int]],
  sequence~ : String,
  title~ : String,
  s1~ : String,
  config? : PrintFormatConfig,
) -> String {
  let cfg = config.unwrap_or(print_format_config_default())
  let (ls, nc) = pf_split(dct)
  pf_make_map_only(ls, title, sequence, cfg) + pf_make_nocut_only(nc, s1, cfg)
}