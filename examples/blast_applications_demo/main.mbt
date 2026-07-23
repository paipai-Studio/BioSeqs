///|
/// Bio.Blast.Applications Demo - BLAST command-line tool wrappers.
/// Inspired by Biopython's Bio.Blast.Applications package.

fn main {
  println("=== Bio.Blast.Applications Demo ===")

  // 1. Create a blastp command line
  println("")
  println("1. Create a blastp command line:")
  let blastp = @bio.blastapp_create_example_blastp()
  println("   Executable: " + blastp.executable)
  println("   Command: " + blastp.blastapp_build_command())

  // 2. Create a blastn command line
  println("")
  println("2. Create a blastn command line:")
  let blastn = @bio.blastapp_create_example_blastn()
  println("   Executable: " + blastn.executable)
  println("   Command: " + blastn.blastapp_build_command())

  // 3. Quick builders
  println("")
  println("3. Quick BLAST command builders:")
  let quick_blastp = @bio.blastapp_quick_blastp("myquery.fasta", "swissprot", "1e-10")
  println("   Quick blastp: " + quick_blastp.blastapp_build_command())
  let quick_blastn = @bio.blastapp_quick_blastn("gene.fa", "nt", "1e-50")
  println("   Quick blastn: " + quick_blastn.blastapp_build_command())
  let quick_blastx = @bio.blastapp_quick_blastx("nt_query.fa", "nr", "1e-5")
  println("   Quick blastx: " + quick_blastx.blastapp_build_command())
  let quick_tblastn = @bio.blastapp_quick_tblastn("protein.fa", "genome_db", "1e-5")
  println("   Quick tblastn: " + quick_tblastn.blastapp_build_command())

  // 4. Makeblastdb command
  println("")
  println("4. Create a makeblastdb command:")
  let mkdb = @bio.blastapp_create_example_makeblastdb()
  println("   Executable: " + mkdb.executable)
  println("   Command: " + mkdb.blastapp_build_command())

  // 5. Inspect available parameters
  println("")
  println("5. Inspect blastp available parameters:")
  let specs = blastp.blastapp_available_params()
  println("   Total parameters: " + specs.length().to_string())
  let mut shown = 0
  let mut i = 0
  while i < specs.length() && shown < 5 {
    let spec = specs[i]
    let val_type = if spec.takes_value { "value" } else { "flag" }
    println("   " + spec.name + " (" + val_type + "): " + spec.description)
    shown = shown + 1
    i = i + 1
  }
  println("   ... and " + (specs.length() - 5).to_string() + " more")

  // 6. Modify parameters
  println("")
  println("6. Modify parameters on a command line:")
  let modified = blastp
    .blastapp_set_parameter("-matrix", "BLOSUM62")
    .blastapp_set_parameter("-gapopen", "11")
    .blastapp_set_parameter("-gapextend", "1")
    .blastapp_set_flag("-html", true)
  println("   Modified command: " + modified.blastapp_build_command())

  // 7. Validate
  println("")
  println("7. Validate command line:")
  let errors = blastp.blastapp_validate()
  println("   Validation errors: " + errors.length().to_string())

  println("")
  println("=== Bio.Blast.Applications Demo Completed ===")
}
