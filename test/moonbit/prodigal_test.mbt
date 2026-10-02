///|
test "prodigal_default_params" {
  let p = @src.prodigal_default_params()
  assert_eq(p.min_gene_len, 60)
  assert_eq(p.max_overlap, 60)
  assert_eq(p.start_codons.length(), 1)
  assert_eq(p.stop_codons.length(), 3)
  assert_true(p.use_rbs)
  assert_true(p.allow_partial)
}

///|
test "prodigal_predict empty sequence" {
  let result = @src.prodigal_predict("", "empty")
  assert_eq(result.num_genes, 0)
  assert_eq(result.sequence_length, 0)
}

///|
test "prodigal_predict short sequence" {
  let result = @src.prodigal_predict("ATG", "short")
  assert_eq(result.num_genes, 0)
}

///|
test "prodigal_predict single gene" {
  // ATG + 19 codons + TAA = 63 nt, just above min_gene_len
  let seq = "ATG" + "GCT".repeat(19) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  assert_true(result.num_genes >= 1)
  let gene = result.genes[0]
  assert_eq(gene.strand, "+")
  assert_eq(gene.start_codon, "ATG")
  assert_eq(gene.stop_codon, "TAA")
  assert_eq(gene.length, 63)
}

///|
test "prodigal_predict gene on reverse strand" {
  // Build a gene on the reverse strand
  let gene = "ATG" + "GCT".repeat(19) + "TAA"
  let rc_gene = @src.reverse_complement(gene)
  // Add some padding
  let seq = "AAAAAA" + rc_gene + "AAAAAA"
  let result = @src.prodigal_predict(seq, "seq1")
  assert_true(result.num_genes >= 1)
  let mut has_reverse = false
  for g in result.genes {
    if g.strand == "-" {
      has_reverse = true
      assert_eq(g.start_codon, "ATG")
    }
  }
  assert_true(has_reverse)
}

///|
test "prodigal_predict multiple genes" {
  // Two genes on forward strand separated by a gap
  let gene1 = "ATG" + "GCT".repeat(20) + "TAA"
  let gene2 = "ATG" + "GAA".repeat(20) + "TGA"
  let gap = "AAAAAA".repeat(5)
  let seq = gene1 + gap + gene2
  let result = @src.prodigal_predict(seq, "seq1")
  assert_true(result.num_genes >= 2)
}

///|
test "prodigal_predict non-overlapping" {
  // Genes should not overlap (beyond max_overlap)
  let gene1 = "ATG" + "GCT".repeat(20) + "TAA"
  let gene2 = "ATG" + "GAA".repeat(20) + "TGA"
  let gap = "AAAAAA".repeat(5)
  let seq = gene1 + gap + gene2
  let result = @src.prodigal_predict(seq, "seq1")
  for i = 0; i < result.genes.length(); i = i + 1 {
    for j = i + 1; j < result.genes.length(); j = j + 1 {
      let g1 = result.genes[i]
      let g2 = result.genes[j]
      let ov_start = if g1.start > g2.start { g1.start } else { g2.start }
      let ov_end = if g1.end_ < g2.end_ { g1.end_ } else { g2.end_ }
      let ov = if ov_end > ov_start { ov_end - ov_start } else { 0 }
      assert_true(ov <= 60)
    }
  }
}

///|
test "prodigal_predict with params" {
  let seq = "ATG" + "GCT".repeat(10) + "TAA"
  let params = @src.prodigal_params(30, 30)
  let result = @src.prodigal_predict_with_params(seq, "seq1", params)
  assert_true(result.num_genes >= 1)
}

///|
test "prodigal_predict min gene length" {
  // Gene shorter than min_gene_len should not be predicted
  let seq = "ATGAAATAA" // 9 nt, much shorter than 60
  let result = @src.prodigal_predict(seq, "seq1")
  assert_eq(result.num_genes, 0)
}

///|
test "prodigal_to_gff format" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  let gff = @src.prodigal_to_gff(result)
  assert_true(gff.contains("##gff-version 3"))
  assert_true(gff.contains("seq1"))
  assert_true(gff.contains("Prodigal"))
  assert_true(gff.contains("gene"))
  assert_true(gff.contains("CDS"))
}

///|
test "prodigal_to_fasta format" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  let fasta = @src.prodigal_to_fasta(result)
  assert_true(fasta.contains(">"))
  assert_true(fasta.contains("score="))
}

///|
test "prodigal_to_fasta_nt format" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  let fasta = @src.prodigal_to_fasta_nt(result, seq)
  assert_true(fasta.contains(">"))
}

///|
test "prodigal_summary format" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  let summary = @src.prodigal_summary(result)
  assert_true(summary.contains("Prodigal gene prediction summary"))
  assert_true(summary.contains("seq1"))
  assert_true(summary.contains("Genes predicted:"))
}

///|
test "prodigal gene protein translation" {
  // ATG GCT TAA = M A *
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  for gene in result.genes {
    if gene.strand == "+" && gene.start_codon == "ATG" {
      // Protein should start with M (methionine)
      assert_true(gene.protein.length() > 0)
      assert_eq(gene.protein[0].unsafe_to_char(), 'M')
    }
  }
}

///|
test "prodigal gene GC content" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  for gene in result.genes {
    assert_true(gene.gc_content >= 0.0 && gene.gc_content <= 1.0)
  }
}

///|
test "prodigal gene coordinates valid" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA" + "GGG" + "ATG" + "GAA".repeat(20) + "TGA"
  let result = @src.prodigal_predict(seq, "seq1")
  for gene in result.genes {
    assert_true(gene.start >= 1)
    assert_true(gene.end_ <= seq.length())
    assert_true(gene.end_ > gene.start)
    assert_eq(gene.length, gene.end_ - gene.start + 1)
  }
}

///|
test "prodigal genes sorted by position" {
  let gene1 = "ATG" + "GCT".repeat(20) + "TAA"
  let gene2 = "ATG" + "GAA".repeat(20) + "TGA"
  let gap = "AAAAAA".repeat(5)
  let seq = gene1 + gap + gene2
  let result = @src.prodigal_predict(seq, "seq1")
  for i = 0; i < result.genes.length() - 1; i = i + 1 {
    assert_true(result.genes[i].start <= result.genes[i + 1].start)
  }
}

///|
test "prodigal all 6 frames" {
  // Create a sequence with genes in different frames
  let seq = "ATG" + "GCT".repeat(20) + "TAA" +
    "A" + "ATG" + "GAA".repeat(20) + "TAG" +
    "AA" + "ATG" + "GGC".repeat(20) + "TGA"
  let result = @src.prodigal_predict(seq, "seq1")
  assert_true(result.num_genes >= 1)
}

///|
test "prodigal gene with GTG start" {
  // GTG is an alternative start codon in bacteria
  let params = @src.prodigal_params_full(
    60, 60, ["ATG", "GTG", "TTG"], ["TAA", "TAG", "TGA"], true, true, false,
  )
  let seq = "GTG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict_with_params(seq, "seq1", params)
  assert_true(result.num_genes >= 1)
}

///|
test "prodigal mean_gc" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  assert_true(result.mean_gc >= 0.0 && result.mean_gc <= 1.0)
}

///|
test "prodigal result fields" {
  let seq = "ATG" + "GCT".repeat(20) + "TAA"
  let result = @src.prodigal_predict(seq, "seq1")
  assert_eq(result.seqid, "seq1")
  assert_eq(result.method, "prodigal")
  assert_eq(result.sequence_length, seq.length())
  assert_eq(result.genes.length(), result.num_genes)
}

///|
test "prodigal partial gene at boundary" {
  // Gene without stop codon at the end
  let seq = "ATG" + "GCT".repeat(30)
  let result = @src.prodigal_predict(seq, "seq1")
  // Should find at least a partial gene
  assert_true(result.num_genes >= 1)
  let mut has_partial = false
  for gene in result.genes {
    if gene.partial {
      has_partial = true
    }
  }
  assert_true(has_partial)
}
