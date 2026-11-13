///|
/// Substitution matrices for sequence alignment, inspired by Biopython Bio.SubsMat.
///
/// Includes BLOSUM62, BLOSUM45, PAM250, and PAM30 matrices commonly used in
/// protein sequence alignment.

///|
/// Substitution matrix type.
pub enum MatrixType {
  BLOSUM62
  BLOSUM45
  PAM250
  PAM30
}

///|
/// Create BLOSUM62 matrix type.
pub fn matrix_blosum62() -> MatrixType {
  MatrixType::BLOSUM62
}

///|
/// Create BLOSUM45 matrix type.
pub fn matrix_blosum45() -> MatrixType {
  MatrixType::BLOSUM45
}

///|
/// Create PAM250 matrix type.
pub fn matrix_pam250() -> MatrixType {
  MatrixType::PAM250
}

///|
/// Create PAM30 matrix type.
pub fn matrix_pam30() -> MatrixType {
  MatrixType::PAM30
}

///|
/// Substitution matrix structure.
pub struct SubstitutionMatrix {
  name : String
  data : Map[String, Map[String, Int32]]
  alphabet : Array[String]
}

///|
/// Create a substitution matrix.
pub fn SubstitutionMatrix::new(
  name : String,
  data : Map[String, Map[String, Int32]],
  alphabet : Array[String],
) -> SubstitutionMatrix {
  SubstitutionMatrix::{ name, data, alphabet, }
}

///|
/// Get matrix name.
pub fn SubstitutionMatrix::name(self : SubstitutionMatrix) -> String {
  self.name
}

///|
/// Get matrix alphabet.
pub fn SubstitutionMatrix::alphabet(self : SubstitutionMatrix) -> Array[String] {
  self.alphabet.copy()
}

///|
/// Get the score for an amino acid pair.
/// Returns 0 if either amino acid is not found in the matrix.
pub fn SubstitutionMatrix::get_score(
  self : SubstitutionMatrix,
  aa1 : String,
  aa2 : String,
) -> Int32 {
  let row = self.data.get(aa1)
  if row is None {
    return 0
  }
  let score = row.unwrap().get(aa2)
  if score is None {
    return 0
  }
  score.unwrap()
}

///|
/// Get the score for an amino acid pair (case-insensitive).
pub fn SubstitutionMatrix::get_score_case_insensitive(
  self : SubstitutionMatrix,
  aa1 : String,
  aa2 : String,
) -> Int32 {
  self.get_score(aa1.to_upper(), aa2.to_upper())
}

///|
/// Check if the matrix contains an amino acid.
pub fn SubstitutionMatrix::contains(
  self : SubstitutionMatrix,
  aa : String,
) -> Bool {
  self.data.contains(aa.to_upper())
}

///|
/// Amino acid alphabet for standard 20 amino acids.
pub fn standard_aa_alphabet() -> Array[String] {
  [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S",
    "T", "W", "Y", "V",
  ]
}

///|
/// Parse a signed integer from a StringView (supports negative values).
fn parse_i32(s : StringView) -> Int32 {
  let str = s.to_owned()
  let n = str.length()
  let mut sign = 1
  let mut i = 0
  if n > 0 && str.unsafe_get(0) == '-' {
    sign = -1
    i = 1
  }
  let mut result = 0
  let zero = '0'.to_int()
  while i < n {
    let c = str.unsafe_get(i).to_int()
    result = result * 10 + (c - zero)
    i = i + 1
  }
  sign * result
}

///|
/// Parse a substitution matrix from a space-separated flat string.
/// Values are in row-major order following standard_aa_alphabet().
fn parse_subs_matrix(
  raw : String,
  alphabet : Array[String],
) -> Map[String, Map[String, Int32]] {
  let tokens = raw.split(" ").to_array()
  let n = alphabet.length()
  let m : Map[String, Map[String, Int32]] = Map([], capacity=n)
  let mut idx = 0
  for i = 0; i < n; i = i + 1 {
    let row_aa = alphabet[i]
    let row : Map[String, Int32] = Map([], capacity=n)
    for j = 0; j < n; j = j + 1 {
      row[alphabet[j]] = parse_i32(tokens[idx])
      idx = idx + 1
    }
    m[row_aa] = row
  }
  m
}

///|
/// BLOSUM62 substitution matrix. Derived from Henikoff & Henikoff (1992).
fn blosum62_data() -> Map[String, Map[String, Int32]] {
  parse_subs_matrix(
    "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0 -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3 -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3 -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -2 -3 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2 -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3 -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3 -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1 -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2 -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1 -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1 -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0 -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3 -2 -2 -2 -2 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4",
    standard_aa_alphabet(),
  )
}

///|
/// BLOSUM45 substitution matrix.
fn blosum45_data() -> Map[String, Map[String, Int32]] {
  parse_subs_matrix(
    "5 -2 -1 -2 -1 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -2 -2 0 -2 7 0 -1 -3 1 0 -2 0 -3 -2 3 -1 -2 -2 -1 -1 -2 -1 -2 -1 0 6 2 -2 0 0 0 1 -2 -3 0 -2 -2 -2 1 0 -4 -2 -3 -2 -1 2 7 -3 0 2 -1 0 -4 -3 0 -3 -4 -1 0 -1 -4 -2 -3 -1 -3 -2 -3 12 -3 -3 -3 -3 -3 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -1 1 0 0 -3 6 2 -2 1 -2 -2 1 0 -4 -1 0 -1 -2 -1 -3 -1 0 0 2 -3 2 6 -2 0 -3 -2 1 -2 -3 0 0 -1 -3 -2 -3 0 -2 0 -1 -3 -2 -2 7 -2 -4 -3 -2 -2 -3 -2 0 -2 -2 -3 -3 -2 0 1 0 -3 1 0 -2 10 -3 -2 -1 0 -2 -2 -1 -2 -3 2 -3 -1 -3 -2 -4 -3 -2 -3 -4 -3 5 2 -3 2 0 -2 -2 -1 -2 0 3 -1 -2 -3 -3 -2 -2 -2 -3 -2 2 5 -3 2 1 -3 -3 -1 -2 0 1 -1 3 0 0 -3 1 1 -2 -1 -3 -3 5 -1 -3 -1 -1 -1 -2 -1 -2 -1 -1 -2 -3 -2 0 -2 -2 0 2 2 -1 6 0 -2 -2 -1 -2 0 1 -2 -2 -2 -4 -2 -4 -3 -3 -2 0 1 -3 0 8 -3 -2 -1 1 3 0 -1 -2 -2 -1 -4 -1 0 -2 -2 -2 -3 -1 -2 -3 9 -1 -1 -3 -3 -3 1 -1 1 0 -1 0 0 0 -1 -2 -3 -1 -2 -2 -1 4 2 -4 -2 -1 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -1 -1 2 5 -3 -1 0 -2 -2 -4 -4 -5 -2 -3 -2 -3 -2 -2 -2 -2 1 -3 -4 -3 15 3 -3 -2 -1 -2 -2 -3 -1 -2 -3 2 0 0 -1 0 3 -3 -2 -1 3 8 -1 0 -2 -3 -3 -1 -3 -3 -3 -3 3 1 -2 1 0 -3 -1 0 -3 -1 5",
    standard_aa_alphabet(),
  )
}

///|
/// PAM250 substitution matrix.
fn pam250_data() -> Map[String, Map[String, Int32]] {
  parse_subs_matrix(
    "2 -2 0 0 -2 0 0 1 -1 -1 -2 -1 -1 -3 1 1 1 -6 -3 0 -2 6 0 -1 -4 1 0 -1 2 -2 -3 3 -1 -3 0 -1 -2 -2 -3 -2 0 0 2 2 -4 1 1 0 2 -2 -3 0 -2 -3 0 1 0 -4 -2 -2 0 -1 2 4 -5 2 3 1 1 -2 -4 0 -2 -3 -1 0 0 -5 -3 -2 -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3 0 -2 -8 0 -2 0 1 1 2 -5 4 2 -1 3 -2 -2 1 0 -3 0 -1 -1 -4 -2 -2 0 0 1 3 -5 2 4 0 1 -2 -3 1 -1 -3 -1 0 0 -4 -2 -2 1 -1 0 1 -3 -1 0 5 -2 -3 -4 -2 -3 -4 0 1 0 -5 -4 -3 -1 2 2 1 -3 3 1 -2 6 -2 -2 0 -1 0 -1 -1 -2 -3 1 -2 -1 -2 -2 -2 -2 -2 -2 -3 -2 5 2 -2 2 0 -2 -2 0 -5 -2 2 -2 -3 -3 -4 -6 -2 -3 -4 -2 2 6 -3 4 0 -3 -3 -2 -4 -3 1 -1 3 0 0 -5 1 1 -2 0 -2 -3 5 0 -3 -1 0 -1 -4 -3 -2 -1 -1 -2 -2 -5 0 -1 -3 -1 2 4 0 6 0 -2 -2 -1 -4 -2 2 -3 -3 -3 -3 -4 -3 -3 -4 0 0 0 -3 0 9 -4 -3 -3 1 3 -1 1 0 0 -1 -3 0 -1 0 -1 -2 -3 -1 -2 -4 6 0 1 -6 -3 -1 1 -1 1 0 0 -1 0 1 -1 -2 -3 0 -2 -3 0 4 1 -5 -2 -2 1 -2 0 0 -2 -1 0 0 -2 0 -2 -1 -1 -3 1 1 5 -5 -3 0 -6 -2 -4 -5 -8 -4 -4 -5 -3 -5 -4 -4 -4 1 -6 -5 -5 17 0 -6 -3 -3 -2 -3 0 -2 -2 -4 1 -2 -3 -3 -2 3 -3 -2 -3 0 7 -1 0 -2 -2 -2 -2 -2 -2 -3 -2 2 1 -2 2 -1 -1 -2 0 -6 -1 5",
    standard_aa_alphabet(),
  )
}

///|
/// PAM30 substitution matrix.
fn pam30_data() -> Map[String, Map[String, Int32]] {
  parse_subs_matrix(
    "9 -3 -2 -2 -2 -2 -2 -1 -3 -2 -3 -2 -2 -4 -1 -1 -1 -6 -4 -1 -3 10 -1 -3 -4 2 -1 -3 0 -4 -4 3 -3 -5 -2 -2 -3 -4 -5 -4 -2 -1 8 2 -4 0 1 -1 0 -4 -4 -1 -3 -5 -2 0 -1 -6 -4 -3 -2 -3 2 9 -5 1 3 0 -2 -4 -5 -2 -4 -5 -2 -1 -2 -7 -4 -4 -2 -4 -4 -5 12 -5 -5 -4 -4 -3 -6 -5 -5 -5 -4 -1 -3 -7 -2 -3 -2 2 0 1 -5 10 3 -2 1 -4 -3 1 -2 -5 -1 -1 -2 -5 -4 -3 -2 -1 1 3 -5 3 9 -2 0 -4 -4 1 -3 -5 -2 -1 -2 -6 -4 -4 -1 -3 -1 0 -4 -2 -2 10 -3 -4 -5 -3 -4 -5 -2 -1 -2 -5 -5 -4 -3 0 0 -2 -4 1 0 -3 9 -4 -4 -1 -3 -3 -2 -2 -3 -5 1 -4 -2 -4 -4 -4 -3 -4 -4 -4 -4 9 2 -3 2 -1 -3 -3 -1 -6 -3 3 -3 -4 -4 -5 -6 -3 -4 -5 -4 2 9 -3 3 -1 -4 -3 -2 -5 -3 2 -2 3 -1 -2 -5 1 1 -3 -1 -3 -3 10 -2 -5 -2 -1 -2 -6 -4 -3 -2 -3 -3 -4 -5 -2 -3 -4 -3 2 3 -2 9 -1 -3 -2 -2 -5 -3 2 -4 -5 -5 -5 -5 -5 -5 -5 -3 -1 -1 -5 -1 10 -5 -4 -3 -1 2 -2 -1 -2 -2 -2 -4 -1 -2 -2 -2 -3 -4 -2 -3 -5 10 -1 -1 -6 -4 -2 -1 -2 0 -1 -1 -1 -1 -1 -2 -3 -3 -1 -2 -4 -1 8 0 -6 -3 -2 -1 -3 -1 -2 -3 -2 -2 -2 -3 -1 -2 -2 -2 -3 -1 0 9 -6 -4 -1 -6 -4 -6 -7 -7 -5 -6 -5 -5 -6 -5 -6 -5 -1 -6 -6 -6 14 -2 -6 -4 -5 -4 -4 -2 -4 -4 -5 1 -3 -3 -4 -3 2 -4 -3 -4 -2 10 -2 -1 -4 -3 -4 -3 -3 -4 -4 -4 3 2 -3 2 -2 -2 -2 -1 -6 -2 9",
    standard_aa_alphabet(),
  )
}

///|
/// Get BLOSUM62 substitution matrix.
pub fn blosum62() -> SubstitutionMatrix {
  SubstitutionMatrix::new("BLOSUM62", blosum62_data(), standard_aa_alphabet())
}

///|
/// Get BLOSUM45 substitution matrix.
pub fn blosum45() -> SubstitutionMatrix {
  SubstitutionMatrix::new("BLOSUM45", blosum45_data(), standard_aa_alphabet())
}

///|
/// Get PAM250 substitution matrix.
pub fn pam250() -> SubstitutionMatrix {
  SubstitutionMatrix::new("PAM250", pam250_data(), standard_aa_alphabet())
}

///|
/// Get PAM30 substitution matrix.
pub fn pam30() -> SubstitutionMatrix {
  SubstitutionMatrix::new("PAM30", pam30_data(), standard_aa_alphabet())
}

///|
/// Get a substitution matrix by type.
pub fn get_matrix(matrix_type : MatrixType) -> SubstitutionMatrix {
  match matrix_type {
    MatrixType::BLOSUM62 => blosum62()
    MatrixType::BLOSUM45 => blosum45()
    MatrixType::PAM250 => pam250()
    MatrixType::PAM30 => pam30()
  }
}

///|
/// Format substitution matrix to string.
pub fn SubstitutionMatrix::to_string(self : SubstitutionMatrix) -> String {
  let sb = StringBuilder()
  sb.write_string(" " + self.alphabet.join(" ") + "\n")

  let mut i = 0
  while i < self.alphabet.length() {
    let aa = self.alphabet[i]
    sb.write_string(aa)
    let mut j = 0
    while j < self.alphabet.length() {
      let score = self.get_score(aa, self.alphabet[j])
      sb.write_string(" " + score.to_string())
      j = j + 1
    }
    sb.write_string("\n")
    i = i + 1
  }

  sb.to_string()
}
