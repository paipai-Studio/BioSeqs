// Performance benchmark for the MoonBit Seq implementation.
//
// Generates a large synthetic DNA sequence and runs each Seq operation
// repeatedly, printing the elapsed time (us) for each. Run with:
//   moon run cmd/bench

///|
/// Generate a deterministic pseudo-random DNA sequence of length `n`
/// using a simple LCG so results are reproducible across runs.
fn gen_dna(n : Int) -> String {
  let buf = StringBuilder::new()
  let bases = ['A', 'C', 'G', 'T']
  let mut state = 12345
  for i = 0; i < n; i = i + 1 {
    state = (state * 1103515245 + 12345) & 0x7fffffff
    let idx = state % 4
    buf.write_char(bases[idx])
  }
  buf.to_string()
}

///|
/// Measure how long `f` takes when called `iters` times; returns total us.
fn bench(label : String, iters : Int, f : () -> Unit) -> Double {
  let ts = @mb.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    f()
  }
  let elapsed_us = @mb.monotonic_clock_end(ts)
  let per_iter = elapsed_us / iters.to_double()
  println(
    "\{label}: \{elapsed_us} us total, \{per_iter} us/iter (\{iters} iters)",
  )
  elapsed_us
}

///|
/// Measure how long `f` takes when called `iters` times; returns total us.
/// Supports functions that may raise errors.
fn bench_err(label : String, iters : Int, f : () -> Unit) -> Double {
  let ts = @mb.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    f()
  }
  let elapsed_us = @mb.monotonic_clock_end(ts)
  let per_iter = elapsed_us / iters.to_double()
  println(
    "\{label}: \{elapsed_us} us total, \{per_iter} us/iter (\{iters} iters)",
  )
  elapsed_us
}

///|

///|
/// Generate a FASTA string with `n_records` records, each of length `seq_len`.
fn gen_fasta(n_records : Int, seq_len : Int) -> String {
  let buf = StringBuilder::new()
  for i = 0; i < n_records; i = i + 1 {
    buf.write_char('>')
    buf.write_string("seq")
    buf.write_string(i.to_string())
    buf.write_string(" Description for sequence ")
    buf.write_string(i.to_string())
    buf.write_char('\n')
    let seq = gen_dna(seq_len)
    buf.write_string(seq)
    buf.write_char('\n')
  }
  buf.to_string()
}

///|
/// Generate a FASTQ string with `n_records` records, each of length `seq_len`.
fn gen_fastq(n_records : Int, seq_len : Int) -> String {
  let buf = StringBuilder::new()
  for i = 0; i < n_records; i = i + 1 {
    buf.write_char('@')
    buf.write_string("seq")
    buf.write_string(i.to_string())
    buf.write_string(" Description ")
    buf.write_string(i.to_string())
    buf.write_char('\n')
    buf.write_string(gen_dna(seq_len))
    buf.write_char('\n')
    buf.write_string("+\n")
    for j = 0; j < seq_len; j = j + 1 {
      buf.write_char('I')
    }
    buf.write_char('\n')
  }
  buf.to_string()
}

///|
fn main {
  let dna = gen_dna(100_000)
  let rna = gen_dna(100_000).replace_all(old="T", new="U")
  let s = @lib.Seq::new(dna)
  let sr = @lib.Seq::new(rna)
  let iters_small = 20
  let iters_find = 200

  println("=== MoonBit Seq benchmark ===")
  println("sequence length: \{dna.char_length()} chars")
  println("iterations (ops): \{iters_small}")
  println("")

  // complement
  ignore(bench("complement", iters_small, fn() { ignore(s.complement()) }))
  ignore(
    bench("complement_rna", iters_small, fn() { ignore(s.complement_rna()) }),
  )
  ignore(
    bench("reverse_complement", iters_small, fn() {
      ignore(s.reverse_complement())
    }),
  )
  ignore(bench("reverse", iters_small, fn() { ignore(s.reverse()) }))

  // transcribe
  ignore(bench("transcribe", iters_small, fn() { ignore(s.transcribe()) }))
  ignore(
    bench("back_transcribe", iters_small, fn() { ignore(sr.back_transcribe()) }),
  )

  // translate
  ignore(
    bench("translate", iters_small, fn() {
      ignore(s.translate()) catch {
        _ => ()
      }
    }),
  )

  // upper / lower
  ignore(bench("upper", iters_small, fn() { ignore(s.upper()) }))
  ignore(bench("lower", iters_small, fn() { ignore(s.lower()) }))

  // count / find
  ignore(bench("count AAA", iters_find, fn() { ignore(s.count("AAA")) }))
  ignore(
    bench("count_overlap AAA", iters_find, fn() {
      ignore(s.count_overlap("AAA"))
    }),
  )
  ignore(bench("find ACGT", iters_find, fn() { ignore(s.find("ACGT")) }))
  ignore(bench("rfind ACGT", iters_find, fn() { ignore(s.rfind("ACGT")) }))

  // startswith / endswith / contains
  ignore(
    bench("startswith ATG", iters_find, fn() { ignore(s.startswith("ATG")) }),
  )
  ignore(bench("endswith TAA", iters_find, fn() { ignore(s.endswith("TAA")) }))
  ignore(bench("contains GGG", iters_find, fn() { ignore(s.contains("GGG")) }))

  // slice
  ignore(
    bench("slice 0..1000", iters_find, fn() {
      ignore(s.slice(start=0, end=1000))
    }),
  )

  // replace
  ignore(
    bench("replace A->N", iters_small, fn() {
      ignore(s.replace(old="A", new="N"))
    }),
  )

  // ---------------------------------------------------------------------------
  // SeqIO benchmarks
  // ---------------------------------------------------------------------------
  println("")
  println("=== MoonBit SeqIO benchmark ===")
  println("")

  let fasta_content = gen_fasta(100, 1000)
  let fastq_content = gen_fastq(100, 1000)

  // FASTA parse
  ignore(
    bench("parse_fasta (100x1000)", iters_small, fn() {
      ignore(@lib.parse_fasta(fasta_content))
    }),
  )

  // FASTA write
  let fasta_records = @lib.parse_fasta(fasta_content)
  ignore(
    bench("write_fasta (100x1000)", iters_small, fn() {
      ignore(@lib.write_fasta(fasta_records))
    }),
  )

  // FASTQ parse
  ignore(
    bench("parse_fastq (100x1000)", iters_small, fn() {
      ignore(@lib.parse_fastq(fastq_content)) catch {
        _ => ()
      }
    }),
  )

  // FASTQ write
  let fq_records = @lib.parse_fastq(fastq_content) catch { _ => Array::new() }
  ignore(
    bench("write_fastq (100x1000)", iters_small, fn() {
      ignore(@lib.write_fastq(fq_records)) catch {
        _ => ()
      }
    }),
  )

  // seqio_parse dispatch
  ignore(
    bench("seqio_parse fasta", iters_small, fn() {
      ignore(@lib.seqio_parse(fasta_content, "fasta")) catch {
        _ => ()
      }
    }),
  )

  ignore(
    bench("seqio_parse fastq", iters_small, fn() {
      ignore(@lib.seqio_parse(fastq_content, "fastq")) catch {
        _ => ()
      }
    }),
  )

  // seqio_to_dict
  ignore(
    bench("seqio_to_dict", iters_small, fn() {
      ignore(@lib.seqio_to_dict(fasta_records))
    }),
  )

  // ---------------------------------------------------------------------------
  // SeqFeature benchmarks
  // ---------------------------------------------------------------------------
  println("")
  println("=== MoonBit SeqFeature benchmark ===")
  println("")

  let feat_qual = @lib.SeqFeature::new(feature_type="gene")
  let pos0 = @lib.Position::new_exact(0)
  let pos1000 = @lib.Position::new_exact(1000)
  let pos2000 = @lib.Position::new_exact(2000)
  let pos3000 = @lib.Position::new_exact(3000)

  // Position creation
  ignore(
    bench("Position::new_exact", iters_find, fn() {
      ignore(@lib.Position::new_exact(100))
    }),
  )

  // FeatureLocation creation using safe constructor
  ignore(
    bench("FeatureLocation::new_safe", iters_find, fn() {
      ignore(@lib.FeatureLocation::new_safe(pos0, pos1000))
    }),
  )

  // SeqFeature creation
  ignore(
    bench("SeqFeature::new", iters_find, fn() {
      ignore(@lib.SeqFeature::new(feature_type="CDS"))
    }),
  )

  // SeqFeature with location (uses pre-created positions)
  ignore(
    bench("SeqFeature::new with location", iters_find, fn() {
      ignore(
        @lib.SeqFeature::new(
          location=Some(@lib.FeatureLocation::new_safe(pos0, pos1000)),
          feature_type="CDS",
        ),
      )
    }),
  )

  // SeqFeature add_qualifier
  ignore(
    bench("SeqFeature::add_qualifier", iters_find, fn() {
      ignore(feat_qual.add_qualifier("gene", "test"))
    }),
  )

  // CompoundLocation
  let loc1 = @lib.FeatureLocation::new_safe(pos0, pos1000)
  let loc2 = @lib.FeatureLocation::new_safe(pos2000, pos3000)
  ignore(
    bench("CompoundLocation::new", iters_find, fn() {
      ignore(@lib.CompoundLocation::new([loc1, loc2], operator="join"))
    }),
  )

  // ---------------------------------------------------------------------------
  // SeqUtils benchmarks
  // ---------------------------------------------------------------------------
  println("")
  println("=== MoonBit SeqUtils benchmark ===")
  println("")

  let dna_seq = gen_dna(10000)

  // CheckSum
  ignore(bench("crc32", iters_find, fn() { ignore(@lib.crc32(dna_seq)) }))
  ignore(bench("gcg", iters_find, fn() { ignore(@lib.gcg(dna_seq)) }))
  ignore(bench("seguid", iters_small, fn() { ignore(@lib.seguid(dna_seq)) }))

  // GC content
  ignore(bench("gc", iters_find, fn() { ignore(@lib.gc(dna_seq)) }))

  // ProteinAnalysis
  let protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV"
  let pa = @lib.ProteinAnalysis::new(protein_seq)
  ignore(
    bench("ProteinAnalysis::count_amino_acids", iters_find, fn() {
      ignore(pa.count_amino_acids())
    }),
  )
  ignore(
    bench("ProteinAnalysis::molecular_weight", iters_find, fn() {
      ignore(pa.molecular_weight())
    }),
  )
  ignore(
    bench("ProteinAnalysis::aromaticity", iters_find, fn() {
      ignore(pa.aromaticity())
    }),
  )
  ignore(
    bench("ProteinAnalysis::gravy", iters_find, fn() { ignore(pa.gravy()) }),
  )

  // MeltingTemp
  ignore(
    bench("tm_wallace", iters_find, fn() { ignore(@lib.tm_wallace(dna_seq)) }),
  )
  ignore(bench("tm_gc", iters_find, fn() { ignore(@lib.tm_gc(dna_seq)) }))

  // ---------------------------------------------------------------------------
  // Bio.Phylo benchmarks
  // ---------------------------------------------------------------------------
  println("")
  println("=== MoonBit Bio.Phylo benchmark ===")
  println("")

  let simple_newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let tree = @lib.parse_newick(simple_newick)

  // Newick parsing
  ignore(
    bench("parse_newick (simple)", iters_find, fn() {
      ignore(@lib.parse_newick(simple_newick))
    }),
  )

  // Newick writing
  ignore(
    bench("write_newick", iters_find, fn() {
      ignore(@lib.write_newick(tree, true))
    }),
  )

  // Count terminals
  ignore(
    bench("Tree::get_leaf_count", iters_find, fn() {
      ignore(@lib.get_leaf_count(tree))
    }),
  )

  // Get terminals
  ignore(
    bench("Tree::get_leaves", iters_find, fn() { ignore(@lib.get_leaves(tree)) }),
  )

  // Get nonterminals
  ignore(
    bench("Tree::get_node_count - get_leaf_count", iters_find, fn() {
      ignore(@lib.get_node_count(tree) - @lib.get_leaf_count(tree))
    }),
  )

  // Find node by name
  ignore(
    bench("Tree::find_node by name", iters_find, fn() {
      ignore(@lib.find_node(tree, "A"))
    }),
  )

  // Get tree depth
  ignore(
    bench("Tree::get_tree_depth", iters_find, fn() {
      ignore(@lib.get_tree_depth(tree))
    }),
  )

  // Total branch length
  ignore(
    bench("Tree::total_branch_length", iters_find, fn() {
      ignore(@lib.total_branch_length(tree))
    }),
  )

  // Get node count
  ignore(
    bench("Tree::get_node_count", iters_find, fn() {
      ignore(@lib.get_node_count(tree))
    }),
  )

  // Draw ASCII
  ignore(
    bench("Tree::write_newick", iters_small, fn() {
      ignore(@lib.write_newick(tree, true))
    }),
  )

  // Generate a larger tree for more meaningful benchmarks
  let large_newick_buf = StringBuilder::new()
  large_newick_buf.write_char('(')
  for i = 0; i < 50; i = i + 1 {
    if i > 0 {
      large_newick_buf.write_char(',')
    }
    large_newick_buf.write_string("X")
    large_newick_buf.write_string(i.to_string())
    large_newick_buf.write_char(':')
    large_newick_buf.write_string(0.1.to_string())
  }
  large_newick_buf.write_string(");")
  let large_newick = large_newick_buf.to_string()
  let large_tree = @lib.parse_newick(large_newick)

  // Large tree operations
  ignore(
    bench("parse_newick (50 leaves)", iters_small, fn() {
      ignore(@lib.parse_newick(large_newick))
    }),
  )

  ignore(
    bench("Tree::get_leaf_count (50 leaves)", iters_small, fn() {
      ignore(@lib.get_leaf_count(large_tree))
    }),
  )

  ignore(
    bench("Tree::get_leaves (50 leaves)", iters_small, fn() {
      ignore(@lib.get_leaves(large_tree))
    }),
  )

  // ---------------------------------------------------------------------------
  // Bio.PDB benchmarks
  // ---------------------------------------------------------------------------
  println("")
  println("=== MoonBit Bio.PDB benchmark ===")
  println("")

  // Generate PDB content
  let pdb_buf = StringBuilder::new()
  pdb_buf.write_string("MODEL     1\n")
  for i = 0; i < 100; i = i + 1 {
    pdb_buf.write_string("ATOM      1  N   ALA A   ")
    pdb_buf.write_string((i + 1).to_string())
    pdb_buf.write_string(
      "      0.000   0.000   0.000  1.00  0.00           N\n",
    )
    pdb_buf.write_string("ATOM      2  CA  ALA A   ")
    pdb_buf.write_string((i + 1).to_string())
    pdb_buf.write_string(
      "      1.500   0.000   0.000  1.00  0.00           C\n",
    )
    pdb_buf.write_string("ATOM      3  C   ALA A   ")
    pdb_buf.write_string((i + 1).to_string())
    pdb_buf.write_string(
      "      3.000   0.000   0.000  1.00  0.00           C\n",
    )
  }
  pdb_buf.write_string("ENDMDL\n")
  pdb_buf.write_string("END\n")
  let pdb_content = pdb_buf.to_string()

  // PDB parsing
  ignore(
    bench("parse_pdb (300 atoms)", iters_small, fn() {
      ignore(@lib.parse_pdb(pdb_content))
    }),
  )

  // PDB writing
  let pdb_struct = @lib.parse_pdb(pdb_content)
  ignore(
    bench("write_pdb (300 atoms)", iters_small, fn() {
      ignore(@lib.write_pdb(pdb_struct))
    }),
  )

  // Structure traversal
  ignore(
    bench("Structure::get_atoms", iters_find, fn() {
      ignore(pdb_struct.get_atoms())
    }),
  )

  ignore(
    bench("Structure::get_residues", iters_find, fn() {
      ignore(pdb_struct.get_residues())
    }),
  )

  ignore(
    bench("Structure::get_chains", iters_find, fn() {
      ignore(pdb_struct.get_chains())
    }),
  )

  // Distance calculation
  let atoms = pdb_struct.get_atoms()
  if atoms.length() >= 2 {
    let atom1 = atoms[0]
    let atom2 = atoms[1]
    ignore(
      bench("Atom::distance", iters_find, fn() { ignore(atom1.distance(atom2)) }),
    )
  }

  // Center of mass
  let residues = pdb_struct.get_residues()
  if residues.length() > 0 {
    let residue = residues[0]
    ignore(
      bench("Residue::center_of_mass", iters_find, fn() {
        ignore(residue.center_of_mass())
      }),
    )
  }

  // Find atom by name
  if residues.length() > 0 {
    let residue = residues[0]
    ignore(
      bench("Residue::find_atom", iters_find, fn() {
        ignore(residue.find_atom("CA"))
      }),
    )
  }

  // Structure distance
  ignore(
    bench("Structure::distance", iters_find, fn() {
      ignore(pdb_struct.distance("CA", "N"))
    }),
  )

  // ---------------------------------------------------------------------------
  // BAM benchmarks
  // ---------------------------------------------------------------------------
  println("")
  println("=== MoonBit BAM benchmark ===")
  println("")

  let test_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAADkAAAAAAAAAZAAAAAAAAAA8AAAGAAEACAAAAP////8AAAAAAAAAAHJlYWQxAIAAAAAAAAAAAAAAAAAAAAA="

  // BAM parse from base64
  ignore(
    bench_err("parse_bam_base64 (1 record)", iters_find, () => {
      ignore(@lib.parse_bam_base64(test_bam_base64)) catch {
        _ => ()
      }
    }),
  )

  // BAM to_sam
  let bam_file_opt = @lib.try_parse_bam_base64(test_bam_base64)
  if bam_file_opt.is_some() {
    let bam_file = bam_file_opt.unwrap()
    ignore(
      bench("BamFile::to_sam", iters_find, fn() { ignore(bam_file.to_sam()) }),
    )

    // BAM get_reference_names
    ignore(
      bench("BamFile::get_reference_names", iters_find, fn() {
        ignore(bam_file.get_reference_names())
      }),
    )

    // BAM num_references
    ignore(
      bench("BamFile::num_references", iters_find, fn() {
        ignore(bam_file.num_references())
      }),
    )

    // BAM len
    ignore(bench("BamFile::len", iters_find, fn() { ignore(bam_file.len()) }))
  }

  // ---------------------------------------------------------------------------
  // CRAM benchmarks
  // ---------------------------------------------------------------------------
  println("")
  println("=== MoonBit CRAM benchmark ===")
  println("")

  // Generate CRAM test data
  let cram_buf = StringBuilder::new()
  cram_buf.write_string("CRAM")
  cram_buf.write_char('\u{0001}')
  cram_buf.write_char('\u{0003}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0001}')
  cram_buf.write_char('\u{0004}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_char('\u{0000}')
  cram_buf.write_string("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n")
  let cram_content = cram_buf.to_string()

  // CRAM magic check
  ignore(
    bench("is_cram_magic", iters_find, fn() {
      ignore(@lib.is_cram_magic(cram_content))
    }),
  )

  // CRAM magic check (negative)
  ignore(
    bench("is_cram_magic (negative)", iters_find, fn() {
      ignore(@lib.is_cram_magic("BAMtest"))
    }),
  )

  println("")
  println("done.")
}
