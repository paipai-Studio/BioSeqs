///|
/// Debug BAM parsing.
fn main {
  let test_bam_base64 = "QkFNNgAAAAFodHRwOi8vd3d3LnczLm9yZy8yMDA1LzA1L3htbGEtaW5zdGFuY2UtdjEwMC9yZWdpb24uZmkKQGhkIhZNbjoxLjYJUU86Y29vcmRpbmF0ZQpBU1FMQ0hSTUVEIhJOczpjaHIxIhtMbjoxMDAwCgAAAAAAAABKAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAEAAAACAAAAAAAAAAIAAAAAAAAAAQAAAAAAAAAFAAAAAAAAAAUAAAAAAAAAAgAAAAAAAAACAAAAAAAAAAIAAAAAAAAAAEAAAAAAAAAAQAAAAAAAAABAAAAAAAAAAIAAAAAAAAAAQAAAAAAAAACAAAAAAAAAAQAAAAAAAAADAAAAAAAAAA=="
  
  println("Base64 length: \{test_bam_base64.length()}")
  
  let raw_data = @src.base64_decode(test_bam_base64)
  println("Decoded length: \{raw_data.length()}")
  
  println("Magic bytes:")
  for i = 0; i < 4; i = i + 1 {
    let byte = raw_data.unsafe_get(i).to_int() & 0xFF
    println("  byte \{i}: 0x\{byte.to_string()} (\{byte})")
  }
  
  println("BAM magic: 66, 65, 77, 1")
  
  let is_bam = @src.is_bam_magic(raw_data)
  println("Is BAM magic: \{is_bam}")
  
  let parsed = @src.parse_bam(raw_data) catch { e => 
    println("Parse error: \{e}")
    return
  }
  
  println("Success!")
  println("References: \{parsed.num_references()}")
  println("Records: \{parsed.len()}")
}
