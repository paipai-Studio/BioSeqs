///|
/// Error raised by skbio alignment operations.
pub suberror AlignError {
  AlignError(String)
}

///|
/// DNA sequence type.
pub struct DNA {
  seq : Seq
} derive(Eq)

///|
/// RNA sequence type.
pub struct RNA {
  seq : Seq
} derive(Eq)

///|
/// Protein sequence type.
pub struct Protein {
  seq : Seq
} derive(Eq)

///|
/// Create a DNA sequence from a string.
pub fn DNA::new(s : String) -> DNA {
  DNA::{ seq: Seq::new(s) }
}

///|
/// Create an RNA sequence from a string.
pub fn RNA::new(s : String) -> RNA {
  RNA::{ seq: Seq::new(s) }
}

///|
/// Create a Protein sequence from a string.
pub fn Protein::new(s : String) -> Protein {
  Protein::{ seq: Seq::new(s) }
}

///|
/// Get the underlying sequence.
pub fn DNA::to_seq(self : DNA) -> Seq {
  self.seq
}

///|
/// Get the underlying sequence.
pub fn RNA::to_seq(self : RNA) -> Seq {
  self.seq
}

///|
/// Get the underlying sequence.
pub fn Protein::to_seq(self : Protein) -> Seq {
  self.seq
}

///|
/// Convert to string.
pub fn DNA::to_string(self : DNA) -> String {
  self.seq.to_string()
}

///|
/// Convert to string.
pub fn RNA::to_string(self : RNA) -> String {
  self.seq.to_string()
}

///|
/// Convert to string.
pub fn Protein::to_string(self : Protein) -> String {
  self.seq.to_string()
}

///|
/// Get length.
pub fn DNA::length(self : DNA) -> Int {
  self.seq.length()
}

///|
/// Get length.
pub fn RNA::length(self : RNA) -> Int {
  self.seq.length()
}

///|
/// Get length.
pub fn Protein::length(self : Protein) -> Int {
  self.seq.length()
}

///|
/// Get character at position.
pub fn DNA::get_char(self : DNA, i : Int) -> Char? {
  self.seq.get_char(i)
}

///|
/// Get character at position.
pub fn RNA::get_char(self : RNA, i : Int) -> Char? {
  self.seq.get_char(i)
}

///|
/// Get character at position.
pub fn Protein::get_char(self : Protein, i : Int) -> Char? {
  self.seq.get_char(i)
}

///|
/// TabularMSA stores a multiple sequence alignment in tabular form.
pub struct TabularMSA {
  records : Array[SeqRecord]
  dtype : String
} derive(Eq)

///|
/// Create a new TabularMSA from SeqRecords.
pub fn TabularMSA::new(
  records : Array[SeqRecord],
  dtype : String,
) -> TabularMSA {
  TabularMSA::{ records, dtype }
}

///|
/// Create a TabularMSA from DNA sequences.
pub fn TabularMSA::from_dna(seqs : Array[DNA]) -> TabularMSA {
  let records = seqs.map(fn(s) {
    SeqRecord::new(s.to_seq(), id="<unknown id>", name="<unknown name>")
  })
  TabularMSA::new(records, "DNA")
}

///|
/// Return the number of sequences (rows) in the alignment.
pub fn TabularMSA::num_records(self : TabularMSA) -> Int {
  self.records.length()
}

///|
/// Return the length of the alignment (columns).
pub fn TabularMSA::get_alignment_length(self : TabularMSA) -> Int {
  if self.records.length() == 0 {
    0
  } else {
    self.records[0].length()
  }
}

///|
/// Get a column from the alignment.
pub fn TabularMSA::column(self : TabularMSA, col : Int) -> String? {
  let align_len = self.get_alignment_length()
  if col < 0 || col >= align_len {
    return None
  }
  let result = StringBuilder::new()
  for rec in self.records {
    match rec.seq.get_char(col) {
      Some(c) => result.write_char(c)
      None => result.write_char('-')
    }
  }
  Some(result.to_string())
}

///|
/// Compute consensus sequence.
pub fn TabularMSA::consensus(self : TabularMSA) -> DNA {
  let align_len = self.get_alignment_length()
  let result = StringBuilder::new()
  for col in 0..<align_len {
    let col_str = self.column(col).unwrap()
    let counts = @hashmap.HashMap([])
    let mut max_char = '-'
    let mut max_count = 0
    for i = 0; i < col_str.length(); i = i + 1 {
      let c = col_str.unsafe_get(i).unsafe_to_char()
      let current = counts.get(c.to_string()).unwrap_or(0) + 1
      counts.set(c.to_string(), current)
      if current > max_count {
        max_count = current
        max_char = c
      }
    }
    result.write_char(max_char)
  }
  DNA::new(result.to_string())
}

///|
/// Compute conservation scores per column.
pub fn TabularMSA::conservation(self : TabularMSA) -> Array[Double] {
  let align_len = self.get_alignment_length()
  let n = self.num_records().to_double()
  let result = Array::new()
  for col in 0..<align_len {
    let col_str = self.column(col).unwrap()
    let counts = @hashmap.HashMap([])
    for i = 0; i < col_str.length(); i = i + 1 {
      let c = col_str.unsafe_get(i).unsafe_to_char()
      counts.set(c.to_string(), counts.get(c.to_string()).unwrap_or(0) + 1)
    }
    let values = counts.values().to_array()
    let mut max_count = 0
    for v in values {
      if v > max_count {
        max_count = v
      }
    }
    result.push(max_count.to_double() / n)
  }
  result
}

///|
/// Global pairwise alignment for nucleotide sequences.
pub fn global_pairwise_align_nucleotide(
  seq1 : DNA,
  seq2 : DNA,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()
  let (msa_str, score) = global_align(s1, s2, 1, -1, -1)
  let rec1 = SeqRecord::new(Seq::new(msa_str.0), id="seq1")
  let rec2 = SeqRecord::new(Seq::new(msa_str.1), id="seq2")
  (
    TabularMSA::new([rec1, rec2], "DNA"),
    score,
    (0, seq1.length()),
    (0, seq2.length()),
  )
}

///|
/// Local pairwise alignment for nucleotide sequences.
pub fn local_pairwise_align_nucleotide(
  seq1 : DNA,
  seq2 : DNA,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()
  let (msa_str, score, pos1, pos2) = local_align(s1, s2, 1, -1, -1)
  let rec1 = SeqRecord::new(Seq::new(msa_str.0), id="seq1")
  let rec2 = SeqRecord::new(Seq::new(msa_str.1), id="seq2")
  (TabularMSA::new([rec1, rec2], "DNA"), score, pos1, pos2)
}

///|
/// Global pairwise alignment for protein sequences.
pub fn global_pairwise_align_protein(
  seq1 : Protein,
  seq2 : Protein,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()
  let (msa_str, score) = global_align(s1, s2, 1, -1, -1)
  let rec1 = SeqRecord::new(Seq::new(msa_str.0), id="seq1")
  let rec2 = SeqRecord::new(Seq::new(msa_str.1), id="seq2")
  (
    TabularMSA::new([rec1, rec2], "Protein"),
    score,
    (0, seq1.length()),
    (0, seq2.length()),
  )
}

///|
/// Local pairwise alignment for protein sequences.
pub fn local_pairwise_align_protein(
  seq1 : Protein,
  seq2 : Protein,
) -> (TabularMSA, Double, (Int, Int), (Int, Int)) {
  let s1 = seq1.to_seq().to_string()
  let s2 = seq2.to_seq().to_string()
  let (msa_str, score, pos1, pos2) = local_align(s1, s2, 1, -1, -1)
  let rec1 = SeqRecord::new(Seq::new(msa_str.0), id="seq1")
  let rec2 = SeqRecord::new(Seq::new(msa_str.1), id="seq2")
  (TabularMSA::new([rec1, rec2], "Protein"), score, pos1, pos2)
}

///|
/// Compute alignment score.
pub fn align_score(msa : TabularMSA) -> Double {
  let align_len = msa.get_alignment_length()
  let n = msa.num_records()
  let mut score = 0.0
  for col in 0..<align_len {
    let col_str = msa.column(col).unwrap()
    let counts = @hashmap.HashMap([])
    for i = 0; i < col_str.length(); i = i + 1 {
      let c = col_str.unsafe_get(i).unsafe_to_char()
      counts.set(c.to_string(), counts.get(c.to_string()).unwrap_or(0) + 1)
    }
    let values = counts.values().to_array()
    let mut max_count = 0
    for v in values {
      if v > max_count {
        max_count = v
      }
    }
    if max_count == n {
      score = score + 1.0
    }
  }
  score + n.to_double()
}

///|
/// Reverse a string.
fn reverse_string(s : String) -> String {
  let n = s.length()
  let result = StringBuilder::new()
  let mut i = n - 1
  while i >= 0 {
    result.write_char(s.unsafe_get(i).unsafe_to_char())
    i = i - 1
  }
  result.to_string()
}

///|
/// Global alignment using Needleman-Wunsch algorithm.
fn global_align(
  s1 : String,
  s2 : String,
  match_score : Double,
  mismatch_score : Double,
  gap_score : Double,
) -> ((String, String), Double) {
  let n = s1.length()
  let m = s2.length()
  let dp = Array::new()
  for i = 0; i <= n; i = i + 1 {
    dp.push(FixedArray::make(m + 1, 0.0))
  }

  for i = 0; i <= n; i = i + 1 {
    dp[i][0] = i.to_double() * gap_score
  }
  for j = 0; j <= m; j = j + 1 {
    dp[0][j] = j.to_double() * gap_score
  }

  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
      let match_mismatch = dp[i - 1][j - 1] +
        (if c1 == c2 { match_score } else { mismatch_score })
      let delete = dp[i - 1][j] + gap_score
      let insert = dp[i][j - 1] + gap_score
      dp[i][j] = if match_mismatch > delete && match_mismatch > insert {
        match_mismatch
      } else if delete > insert {
        delete
      } else {
        insert
      }
    }
  }

  let mut i = n
  let mut j = m
  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()

  while i > 0 || j > 0 {
    if i > 0 && j > 0 {
      let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
      let match_mismatch = dp[i - 1][j - 1] +
        (if c1 == c2 { match_score } else { mismatch_score })
      if dp[i][j] == match_mismatch {
        r1.write_char(c1)
        r2.write_char(c2)
        i = i - 1
        j = j - 1
        continue
      }
    }
    if i > 0 && dp[i][j] == dp[i - 1][j] + gap_score {
      r1.write_char(s1.unsafe_get(i - 1).unsafe_to_char())
      r2.write_char('-')
      i = i - 1
    } else if j > 0 {
      r1.write_char('-')
      r2.write_char(s2.unsafe_get(j - 1).unsafe_to_char())
      j = j - 1
    }
  }

  ((reverse_string(r1.to_string()), reverse_string(r2.to_string())), dp[n][m])
}

///|
/// Local alignment using Smith-Waterman algorithm.
fn local_align(
  s1 : String,
  s2 : String,
  match_score : Double,
  mismatch_score : Double,
  gap_score : Double,
) -> ((String, String), Double, (Int, Int), (Int, Int)) {
  let n = s1.length()
  let m = s2.length()
  let dp = Array::new()
  for i = 0; i <= n; i = i + 1 {
    dp.push(FixedArray::make(m + 1, 0.0))
  }
  let mut max_score = 0.0
  let mut max_i = 0
  let mut max_j = 0

  for i = 1; i <= n; i = i + 1 {
    for j = 1; j <= m; j = j + 1 {
      let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
      let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
      let match_mismatch = dp[i - 1][j - 1] +
        (if c1 == c2 { match_score } else { mismatch_score })
      let delete = dp[i - 1][j] + gap_score
      let insert = dp[i][j - 1] + gap_score
      let val = if match_mismatch > delete &&
        match_mismatch > insert &&
        match_mismatch > 0.0 {
        match_mismatch
      } else if delete > insert && delete > 0.0 {
        delete
      } else if insert > 0.0 {
        insert
      } else {
        0.0
      }
      dp[i][j] = val
      if val > max_score {
        max_score = val
        max_i = i
        max_j = j
      }
    }
  }

  let mut i = max_i
  let mut j = max_j
  let r1 = StringBuilder::new()
  let r2 = StringBuilder::new()

  while i > 0 && j > 0 && dp[i][j] > 0.0 {
    let c1 = s1.unsafe_get(i - 1).unsafe_to_char()
    let c2 = s2.unsafe_get(j - 1).unsafe_to_char()
    let match_mismatch = dp[i - 1][j - 1] +
      (if c1 == c2 { match_score } else { mismatch_score })
    if dp[i][j] == match_mismatch {
      r1.write_char(c1)
      r2.write_char(c2)
      i = i - 1
      j = j - 1
    } else if dp[i][j] == dp[i - 1][j] + gap_score {
      r1.write_char(c1)
      r2.write_char('-')
      i = i - 1
    } else {
      r1.write_char('-')
      r2.write_char(c2)
      j = j - 1
    }
  }

  (
    (reverse_string(r1.to_string()), reverse_string(r2.to_string())),
    max_score,
    (i, max_i),
    (j, max_j),
  )
}
