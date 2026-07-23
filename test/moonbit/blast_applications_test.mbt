///|
/// Bio.Blast.Applications module tests - BLAST command-line tool wrappers

///|
test "blastapp_ncbi_blastp_commandline" {
  let cmd = @src.ncbi_blastp_commandline()
  assert_eq(cmd.executable, "blastp")
  // Should have common params
  let specs = cmd.blastapp_available_params()
  assert_true(specs.length() > 20)
  // Should have blastp-specific params
  let mut has_remote = false
  let mut has_subject = false
  let mut i = 0
  while i < specs.length() {
    if specs[i].name == "-remote" {
      has_remote = true
    }
    if specs[i].name == "-subject" {
      has_subject = true
    }
    i = i + 1
  }
  assert_true(has_remote)
  assert_true(has_subject)
}

///|
test "blastapp_ncbi_blastn_commandline" {
  let cmd = @src.ncbi_blastn_commandline()
  assert_eq(cmd.executable, "blastn")
  let specs = cmd.blastapp_available_params()
  let mut has_task = false
  let mut has_strand = false
  let mut i = 0
  while i < specs.length() {
    if specs[i].name == "-task" {
      has_task = true
    }
    if specs[i].name == "-strand" {
      has_strand = true
    }
    i = i + 1
  }
  assert_true(has_task)
  assert_true(has_strand)
}

///|
test "blastapp_ncbi_blastx_commandline" {
  let cmd = @src.ncbi_blastx_commandline()
  assert_eq(cmd.executable, "blastx")
  let specs = cmd.blastapp_available_params()
  let mut has_query_gencode = false
  let mut i = 0
  while i < specs.length() {
    if specs[i].name == "-query_gencode" {
      has_query_gencode = true
    }
    i = i + 1
  }
  assert_true(has_query_gencode)
}

///|
test "blastapp_ncbi_tblastn_commandline" {
  let cmd = @src.ncbi_tblastn_commandline()
  assert_eq(cmd.executable, "tblastn")
  let specs = cmd.blastapp_available_params()
  let mut has_db_gencode = false
  let mut i = 0
  while i < specs.length() {
    if specs[i].name == "-db_gencode" {
      has_db_gencode = true
    }
    i = i + 1
  }
  assert_true(has_db_gencode)
}

///|
test "blastapp_ncbi_tblastx_commandline" {
  let cmd = @src.ncbi_tblastx_commandline()
  assert_eq(cmd.executable, "tblastx")
  let specs = cmd.blastapp_available_params()
  let mut has_strand = false
  let mut has_query_gencode = false
  let mut has_db_gencode = false
  let mut i = 0
  while i < specs.length() {
    if specs[i].name == "-strand" {
      has_strand = true
    }
    if specs[i].name == "-query_gencode" {
      has_query_gencode = true
    }
    if specs[i].name == "-db_gencode" {
      has_db_gencode = true
    }
    i = i + 1
  }
  assert_true(has_strand)
  assert_true(has_query_gencode)
  assert_true(has_db_gencode)
}

///|
test "blastapp_ncbi_psiblast_commandline" {
  let cmd = @src.ncbi_psiblast_commandline()
  assert_eq(cmd.executable, "psiblast")
  let specs = cmd.blastapp_available_params()
  let mut has_num_iterations = false
  let mut has_in_pssm = false
  let mut i = 0
  while i < specs.length() {
    if specs[i].name == "-num_iterations" {
      has_num_iterations = true
    }
    if specs[i].name == "-in_pssm" {
      has_in_pssm = true
    }
    i = i + 1
  }
  assert_true(has_num_iterations)
  assert_true(has_in_pssm)
}

///|
test "blastapp_ncbi_makeblastdb_commandline" {
  let cmd = @src.ncbi_makeblastdb_commandline()
  assert_eq(cmd.executable, "makeblastdb")
  let specs = cmd.blastapp_available_params()
  let mut has_in = false
  let mut has_dbtype = false
  let mut has_parse_seqids = false
  let mut i = 0
  while i < specs.length() {
    if specs[i].name == "-in" {
      has_in = true
    }
    if specs[i].name == "-dbtype" {
      has_dbtype = true
    }
    if specs[i].name == "-parse_seqids" {
      has_parse_seqids = true
    }
    i = i + 1
  }
  assert_true(has_in)
  assert_true(has_dbtype)
  assert_true(has_parse_seqids)
}

///|
test "blastapp_set_parameter" {
  let cmd = @src.ncbi_blastp_commandline()
    .blastapp_set_parameter("-query", "test.fasta")
    .blastapp_set_parameter("-db", "nr")
    .blastapp_set_parameter("-evalue", "0.001")
  let q = cmd.blastapp_get_parameter("-query")
  assert_true(q is Some(_))
  assert_eq(q.unwrap(), "test.fasta")
  let db = cmd.blastapp_get_parameter("-db")
  assert_eq(db.unwrap(), "nr")
  let ev = cmd.blastapp_get_parameter("-evalue")
  assert_eq(ev.unwrap(), "0.001")
}

///|
test "blastapp_set_flag" {
  let cmd = @src.ncbi_blastp_commandline()
    .blastapp_set_flag("-show_gis", true)
    .blastapp_set_flag("-html", false)
  let sg = cmd.blastapp_get_flag("-show_gis")
  assert_true(sg is Some(_))
  assert_true(sg.unwrap())
  let html = cmd.blastapp_get_flag("-html")
  assert_true(html is Some(_))
  assert_true(!html.unwrap())
}

///|
test "blastapp_build_command" {
  let cmd = @src.blastapp_quick_blastp("query.fasta", "nr", "0.001")
  let cmd_str = cmd.blastapp_build_command()
  // Command should start with "blastp"
  assert_true(cmd_str.length() > 6)
  // Should contain the query file
  let has_query = cmd_str.contains("query.fasta")
  assert_true(has_query)
  let has_db = cmd_str.contains("nr")
  assert_true(has_db)
  let has_evalue = cmd_str.contains("0.001")
  assert_true(has_evalue)
}

///|
test "blastapp_quick_blastn" {
  let cmd = @src.blastapp_quick_blastn("gene.fasta", "nt", "1e-50")
  assert_eq(cmd.executable, "blastn")
  let q = cmd.blastapp_get_parameter("-query")
  assert_eq(q.unwrap(), "gene.fasta")
  let db = cmd.blastapp_get_parameter("-db")
  assert_eq(db.unwrap(), "nt")
}

///|
test "blastapp_quick_blastx" {
  let cmd = @src.blastapp_quick_blastx("nt_query.fasta", "swissprot", "1e-10")
  assert_eq(cmd.executable, "blastx")
  let q = cmd.blastapp_get_parameter("-query")
  assert_eq(q.unwrap(), "nt_query.fasta")
}

///|
test "blastapp_quick_tblastn" {
  let cmd = @src.blastapp_quick_tblastn("protein.fasta", "genome_db", "1e-5")
  assert_eq(cmd.executable, "tblastn")
  let db = cmd.blastapp_get_parameter("-db")
  assert_eq(db.unwrap(), "genome_db")
}

///|
test "blastapp_quick_makeblastdb" {
  let cmd = @src.blastapp_quick_makeblastdb("genome.fa", "nucl", "MyGenome")
  assert_eq(cmd.executable, "makeblastdb")
  let infile = cmd.blastapp_get_parameter("-in")
  assert_eq(infile.unwrap(), "genome.fa")
  let dbtype = cmd.blastapp_get_parameter("-dbtype")
  assert_eq(dbtype.unwrap(), "nucl")
  let title = cmd.blastapp_get_parameter("-title")
  assert_eq(title.unwrap(), "MyGenome")
}

///|
test "blastapp_create_example_blastp" {
  let cmd = @src.blastapp_create_example_blastp()
  let cmd_str = cmd.blastapp_build_command()
  assert_true(cmd_str.contains("blastp"))
  assert_true(cmd_str.contains("query.fasta"))
  assert_true(cmd_str.contains("nr"))
  assert_true(cmd_str.contains("0.001"))
  // Flag should appear
  assert_true(cmd_str.contains("-show_gis"))
}

///|
test "blastapp_create_example_blastn" {
  let cmd = @src.blastapp_create_example_blastn()
  let cmd_str = cmd.blastapp_build_command()
  assert_true(cmd_str.contains("blastn"))
  assert_true(cmd_str.contains("gene.fasta"))
  assert_true(cmd_str.contains("nt"))
}

///|
test "blastapp_create_example_makeblastdb" {
  let cmd = @src.blastapp_create_example_makeblastdb()
  let cmd_str = cmd.blastapp_build_command()
  assert_true(cmd_str.contains("makeblastdb"))
  assert_true(cmd_str.contains("genome.fa"))
  assert_true(cmd_str.contains("nucl"))
  assert_true(cmd_str.contains("-parse_seqids"))
  assert_true(cmd_str.contains("-hash_index"))
}

///|
test "blastapp_validate" {
  // No required parameters in our specs, so validation should pass
  let cmd = @src.ncbi_blastp_commandline()
  let errors = cmd.blastapp_validate()
  assert_eq(errors.length(), 0)
}

///|
test "blastapp_param_spec" {
  let spec = @src.BlastParamSpec::new("-evalue", "E-value threshold", true, "10.0", false)
  assert_eq(spec.name, "-evalue")
  assert_eq(spec.description, "E-value threshold")
  assert_true(spec.takes_value)
  assert_eq(spec.default_value, "10.0")
  assert_true(!spec.required)
}
