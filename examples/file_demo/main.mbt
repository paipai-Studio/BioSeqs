///|
/// Example: file_demo
/// Demonstrates the Bio.File module (SmartFile with compression support).

fn main {
  // Example 1: Create a SmartFile with automatic compression detection
  println("=== SmartFile Demo ===")
  println("")
  
  // Plain text file
  let plain_file = @bio.SmartFile::new("data.txt", mode="r")
  let plain_path = plain_file.get_path()
  let plain_fmt = plain_file.get_format().to_string()
  let plain_mode = plain_file.get_mode()
  println("Plain file: \{plain_path}")
  println("Compression: \{plain_fmt}")
  println("Mode: \{plain_mode}")
  println("")
  
  // Gzip compressed file
  let gz_file = @bio.SmartFile::new("data.txt.gz", mode="r")
  let gz_path = gz_file.get_path()
  let gz_fmt = gz_file.get_format().to_string()
  println("Gzip file: \{gz_path}")
  println("Compression: \{gz_fmt}")
  println("")
  
  // Bzip2 compressed file
  let bz2_file = @bio.SmartFile::new("data.txt.bz2", mode="r")
  let bz2_path = bz2_file.get_path()
  let bz2_fmt = bz2_file.get_format().to_string()
  println("Bzip2 file: \{bz2_path}")
  println("Compression: \{bz2_fmt}")
  println("")
  
  // Example 2: Reading from a SmartFile
  println("=== Reading from SmartFile ===")
  let content = "FASTA sequence data\nline 2\nline 3"
  let reader = @bio.smart_file_from_string("sequence.fasta", content)
  println("File content:")
  println(reader.read_all())
  println("")
  
  // Example 3: Writing to a SmartFile
  println("=== Writing to SmartFile ===")
  let writer = @bio.SmartFile::new("output.txt", mode="w")
  let w1 = writer.write_line(">seq1 Human gene")
  let w2 = w1.write_line("ACGTACGTACGT")
  let written = w2.get_written_content()
  println("Written: \{written}")
  println("")
  
  // Example 4: Open/Close lifecycle
  println("=== File Lifecycle ===")
  let f = @bio.SmartFile::new("test.txt", mode="r")
  let is_open = f.is_open()
  println("Is open: \{is_open}")
  let opened = f.open()
  let opened_open = opened.is_open()
  println("After open: \{opened_open}")
  let closed = opened.close()
  let closed_open = closed.is_open()
  println("After close: \{closed_open}")
  println("")
  
  // Example 5: Compression detection utilities
  println("=== Compression Utilities ===")
  let has_gz = @bio.has_compression_extension("data.gz")
  println("Has compression extension (data.gz): \{has_gz}")
  let has_txt = @bio.has_compression_extension("data.txt")
  println("Has compression extension (data.txt): \{has_txt}")
  let stripped_gz = @bio.strip_compression_extension("data.txt.gz")
  println("Strip extension (data.txt.gz): \{stripped_gz}")
  let stripped_txt = @bio.strip_compression_extension("data.txt")
  println("Strip extension (data.txt): \{stripped_txt}")
  println("")
  
  // Example 6: Summary
  println("=== Summary ===")
  println(plain_file.summary())
}
