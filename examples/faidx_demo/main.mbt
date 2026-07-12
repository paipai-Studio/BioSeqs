///|
/// FASTA indexing example (pyfaidx style).
///
/// Demonstrates building .fai index, random access to subsequences,
/// and index serialization — compatible with samtools faidx format.
///
/// Run: moon run examples/faidx_demo/main.mbt
fn main {
  // Multi-line FASTA content (like a real genome file)
  let fasta_content = ">chr1\nATGCATGCATGCATGCATGCATGCATGCATGC\nATGCATGCATGCATGCATGCATGCATGCATGC\nATGCATGCATGCATGC\n>chr2\nGGGGAAAATTTTCCCC\nGGGGAAAATTTTCCCC\n>chr3\nACGTACGTACGTACGTACGTACGTACGTACGT"

  // === Build FASTA index ===
  let fasta = @bio.Fasta::from_content(fasta_content)
  println("=== FASTA Index (pyfaidx style) ===")

  let names = fasta.get_names()
  println("  Indexed sequences: \{names}")

  // === Check if a sequence exists ===
  let has_chr1 = fasta.contains("chr1")
  let has_chrX = fasta.contains("chrX")
  println("\n  contains('chr1'): \{has_chr1}")
  println("  contains('chrX'): \{has_chrX}")

  // === Get sequence lengths ===
  println("\n=== Sequence Lengths ===")
  for i = 0; i < names.length(); i = i + 1 {
    let len = fasta.get_length(names[i])
    match len {
      Some(l) => println("  \{names[i]}: \{l} bp")
      None => println("  \{names[i]}: unknown")
    }
  }

  // === Random access: fetch subsequence ===
  println("\n=== Random Access (fetch) ===")

  try {
    // Fetch chr1[0..16]
    let sub1 = fasta.fetch("chr1", 0, 16)
    match sub1 {
      Some(s) => println("  chr1[0:16]:  \{s.to_string()}")
      None => println("  chr1[0:16]:  not found")
    }

    // Fetch chr1[16..32] (crosses line boundary)
    let sub2 = fasta.fetch("chr1", 16, 32)
    match sub2 {
      Some(s) => println("  chr1[16:32]: \{s.to_string()}")
      None => println("  chr1[16:32]: not found")
    }

    // Fetch chr2[0..8]
    let sub3 = fasta.fetch("chr2", 0, 8)
    match sub3 {
      Some(s) => println("  chr2[0:8]:   \{s.to_string()}")
      None => println("  chr2[0:8]:   not found")
    }
  } catch {
    _ => println("  fetch error occurred")
  }

  // === Get full sequence ===
  println("\n=== Full Sequence ===")
  let full = fasta.get_seq("chr3")
  match full {
    Some(s) => println("  chr3: \{s.to_string()}")
    None => println("  chr3: not found")
  }

  // === Build .fai index and serialize ===
  println("\n=== .fai Index (samtools faidx compatible) ===")
  let index = @bio.build_fai(fasta_content)
  let fai_str = @bio.write_fai(index)
  println(fai_str)

  // === Parse .fai index back ===
  try {
    let parsed_index = @bio.parse_fai(fai_str)
    println("=== Parsed .fai Index ===")
    println("  Sequences: \{parsed_index.num_sequences()}")
    let idx_names = parsed_index.get_names()
    for i = 0; i < idx_names.length(); i = i + 1 {
      let rec = parsed_index.get_record(idx_names[i])
      match rec {
        Some(r) =>
          println(
            "  \{idx_names[i]}: len=\{r.length} offset=\{r.offset} line_len=\{r.line_len}",
          )
        None => println("  \{idx_names[i]}: not found")
      }
    }
  } catch {
    _ => println("Failed to parse .fai")
  }
}
