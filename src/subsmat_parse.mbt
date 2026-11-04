///|
/// Parse a substitution matrix from NCBI/R format.
/// Expected format: header line with amino acids, followed by rows of scores.
pub fn parse_matrix(content : String) -> @data.SubstitutionMatrix? {
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  let mut alphabet : Array[String] = []
  let data = Map([], capacity=20)

  while i < lines.length() {
    let line = lines[i].trim().to_owned()
    if line.length() == 0 ||
      starts_with(line, 0, "#") ||
      starts_with(line, 0, "//") {
      i = i + 1
      continue
    }

    let fields = split_by_char(line, ' '.to_int())
    let tokens : Array[String] = []
    let mut j = 0
    while j < fields.length() {
      let f = fields[j].trim().to_owned()
      if f.length() > 0 {
        tokens.push(f)
      }
      j = j + 1
    }

    if tokens.length() == 0 {
      i = i + 1
      continue
    }

    if alphabet.length() == 0 {
      alphabet = tokens.copy()
    } else {
      if tokens.length() < 2 {
        i = i + 1
        continue
      }
      let aa = tokens[0].to_upper()
      let row = Map([], capacity=alphabet.length())

      j = 1
      while j < tokens.length() && j - 1 < alphabet.length() {
        let score_opt = parse_double(tokens[j])
        if score_opt is Some(_) {
          row[alphabet[j - 1]] = score_opt.unwrap().to_int()
        }
        j = j + 1
      }

      data[aa] = row
    }

    i = i + 1
  }

  if alphabet.length() == 0 || data.length() == 0 {
    return None
  }

  Some(@data.SubstitutionMatrix::new("custom", data, alphabet))
}
