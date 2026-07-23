///|
/// Bio.SeqUtils.CheckSum - Sequence checksums.

///|
pub struct CheckSumResult {
  checksum : String
  method : String
}

///|
pub fn CheckSumResult::new(checksum : String, method : String) -> CheckSumResult {
  CheckSumResult::{ checksum, method }
}

///|
pub fn checksum_gcg(seq : String) -> CheckSumResult {
  let seq_upper = seq.to_upper()
  
  let mut checksum = 0
  let mut idx = 0
  while idx < seq_upper.length() {
    let c = substring(seq_upper, idx, 1)
    if c == "A" { checksum = (checksum + 1) * 11 }
    else if c == "B" { checksum = (checksum + 2) * 11 }
    else if c == "C" { checksum = (checksum + 3) * 11 }
    else if c == "D" { checksum = (checksum + 4) * 11 }
    else if c == "E" { checksum = (checksum + 5) * 11 }
    else if c == "F" { checksum = (checksum + 6) * 11 }
    else if c == "G" { checksum = (checksum + 7) * 11 }
    else if c == "H" { checksum = (checksum + 8) * 11 }
    else if c == "I" { checksum = (checksum + 9) * 11 }
    else if c == "J" { checksum = (checksum + 10) * 11 }
    else if c == "K" { checksum = (checksum + 11) * 11 }
    else if c == "L" { checksum = (checksum + 12) * 11 }
    else if c == "M" { checksum = (checksum + 13) * 11 }
    else if c == "N" { checksum = (checksum + 14) * 11 }
    else if c == "O" { checksum = (checksum + 15) * 11 }
    else if c == "P" { checksum = (checksum + 16) * 11 }
    else if c == "Q" { checksum = (checksum + 17) * 11 }
    else if c == "R" { checksum = (checksum + 18) * 11 }
    else if c == "S" { checksum = (checksum + 19) * 11 }
    else if c == "T" { checksum = (checksum + 20) * 11 }
    else if c == "U" { checksum = (checksum + 21) * 11 }
    else if c == "V" { checksum = (checksum + 22) * 11 }
    else if c == "W" { checksum = (checksum + 23) * 11 }
    else if c == "X" { checksum = (checksum + 24) * 11 }
    else if c == "Y" { checksum = (checksum + 25) * 11 }
    else if c == "Z" { checksum = (checksum + 26) * 11 }
    idx = idx + 1
  }
  
  checksum = checksum % 10000
  
  let checksum_str = if checksum < 10 {
    "000" + checksum.to_string()
  } else if checksum < 100 {
    "00" + checksum.to_string()
  } else if checksum < 1000 {
    "0" + checksum.to_string()
  } else {
    checksum.to_string()
  }
  
  CheckSumResult::new(checksum_str, "gcg")
}

///|
pub fn checksum_seguid(seq : String) -> CheckSumResult {
  let seq_upper = seq.to_upper()
  
  let mut hash_val = 0
  let mut idx = 0
  while idx < seq_upper.length() {
    let c = substring(seq_upper, idx, 1)
    if c == "A" { hash_val = hash_val * 31 + 65 }
    else if c == "T" { hash_val = hash_val * 31 + 84 }
    else if c == "C" { hash_val = hash_val * 31 + 67 }
    else if c == "G" { hash_val = hash_val * 31 + 71 }
    else { hash_val = hash_val * 31 + 78 }
    idx = idx + 1
  }
  
  let table = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  
  let mut result = ""
  let mut h = hash_val
  while h > 0 {
    let rem = h % 64
    let mut j = 0
    while j < table.length() {
      if j == rem {
        result = substring(table, j, 1) + result
        break
      }
      j = j + 1
    }
    h = h / 64
  }
  
  CheckSumResult::new(result, "seguid")
}

///|
pub fn create_example_sequence() -> String {
  "ACDEFGHIKLMNPQRSTVWY"
}

///|
pub fn verify_checksum(seq : String, expected_checksum : String, method : String) -> Bool {
  let result = if method == "gcg" {
    checksum_gcg(seq)
  } else {
    checksum_seguid(seq)
  }
  
  result.checksum == expected_checksum
}