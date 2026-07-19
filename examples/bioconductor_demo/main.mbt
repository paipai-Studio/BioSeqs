///|
fn main {
  println("=== ChIPseeker Peak Annotation Example ===")
  println("Inspired by Bioconductor ChIPseeker package")
  println("")

  let peaks = [
    ("chr1", 1000, 2000),
    ("chr1", 5000, 6000),
    ("chr2", 10000, 11000),
  ]

  let genes = [
    @bio.Gene::new("gene1", "chr1", 1, 5000, @bio.strand_plus()),
    @bio.Gene::new("gene2", "chr2", 10001, 20000, @bio.strand_minus()),
  ]
  let transcripts : Array[@bio.Transcript] = []
  let exons : Array[@bio.Exon] = []
  let cds : Array[@bio.CDS] = []
  let txdb = @bio.TxDb::new(
    genes, transcripts, exons, cds, "hg38", "Homo sapiens",
  )

  let peak_objs : Array[@bio.Peak] = Array::new()
  let mut i = 0
  while i < peaks.length() {
    let (chr, start, end) = peaks[i]
    peak_objs.push(@bio.Peak::new(chr, start, end, "peak_" + i.to_string()))
    i = i + 1
  }

  let result = @bio.annotate_peaks(peak_objs, txdb)

  println("1. Peak Annotation Results")
  println("   Total annotations: \{result.annotations.length().to_string()}")
  println("")

  let mut i = 0
  while i < result.annotations.length() && i < 3 {
    let ann = result.annotations[i]
    println("   Peak \{i + 1}: \{ann.chr}:\{ann.start}-\{ann.end}")
    println("     - Annotation: \{ann.annotation}")
    println("     - Gene: \{ann.gene_name}")
    println("     - Distance to TSS: \{ann.distance_to_tss.to_string()}")
    println("     - Feature type: \{ann.feature_type}")
    println("")
    i = i + 1
  }

  println("2. Annotation Summary")
  let keys = result.summary.keys().collect()
  let mut j = 0
  while j < keys.length() {
    let key = keys[j]
    println("   \{key}: \{result.summary[key].to_string()}")
    j = j + 1
  }
  println("")

  println("3. Annotation Statistics")
  let stats = @bio.get_annotation_statistics(result)
  let stat_keys = stats.keys().collect()
  let mut k = 0
  while k < stat_keys.length() {
    let stat_key = stat_keys[k]
    println("   \{stat_key}: \{stats[stat_key].to_string()}")
    k = k + 1
  }
  println("")

  println("4. Promoter Annotations")
  let promoters = @bio.get_promoter_annotations(result)
  println("   Promoter peaks: \{promoters.length().to_string()}")
  println("")

  println("5. Annotation Plot")
  println(@bio.plot_annotation_bar(result.summary))
  println("")

  println("\n")
  println("=== DOSE Disease Ontology Enrichment Example ===")
  println("Inspired by Bioconductor DOSE package")
  println("")

  let gene_list = ["TP53", "BRCA1", "BRCA2", "ATM", "CDKN1A"]
  let background_genes = [
    "TP53", "BRCA1", "BRCA2", "ATM", "CDKN1A", "GAPDH", "ACTB", "MYC", "EGFR", "VEGFA",
  ]

  let disease_db = Map([], capacity=4)
  disease_db["DOID:1234"] = ["TP53", "BRCA1", "CDKN1A"]
  disease_db["DOID:5678"] = ["BRCA2", "ATM"]
  disease_db["DOID:9012"] = ["GAPDH", "ACTB"]
  disease_db["DOID:3456"] = ["MYC", "EGFR"]

  let dose_result = @bio.enrich_disease(gene_list, disease_db, background_genes)

  println("1. Enrichment Results")
  println("   Total enriched terms: \{dose_result.terms.length().to_string()}")
  println("   Query genes: \{dose_result.query_count.to_string()}")
  println("   Background genes: \{dose_result.background_count.to_string()}")
  println("")

  println("2. Top 3 Enriched Terms")
  let top_terms = @bio.get_top_terms(dose_result, 3)
  i = 0
  while i < top_terms.length() {
    let term = top_terms[i]
    println("   Term \{i + 1}:")
    println("     - ID: \{term.term_id}")
    println("     - Name: \{term.term_name}")
    println("     - Gene count: \{term.gene_count.to_string()}")
    println("     - P-value: \{term.p_value.to_string()}")
    println("     - Adjusted P-value: \{term.adjusted_p_value.to_string()}")
    println("     - Enrichment score: \{term.enrichment_score.to_string()}")
    println("")
    i = i + 1
  }

  println("3. Filtered Results (adjusted p-value < 0.05)")
  let filtered = @bio.filter_results(dose_result, 0.05)
  println("   Significant terms: \{filtered.terms.length().to_string()}")
  println("")

  println("4. Summary")
  println(@bio.summarize_result(dose_result))
  println("")

  println("\n")
  println("=== ReactomePA Pathway Analysis Example ===")
  println("Inspired by Bioconductor ReactomePA package")
  println("")

  let reactome_gene_list = ["TP53", "BRCA1", "ATM", "CDKN1A"]
  let reactome_bg_genes = [
    "TP53", "BRCA1", "ATM", "CDKN1A", "GAPDH", "ACTB", "MYC",
  ]

  let reactome_db = Map([], capacity=4)
  reactome_db["R-HSA-1234"] = (
    "Cell Cycle",
    ["TP53", "BRCA1", "CDKN1A"],
    "Cellular Processes",
  )
  reactome_db["R-HSA-5678"] = (
    "DNA Repair",
    ["BRCA1", "ATM"],
    "DNA Replication",
  )
  reactome_db["R-HSA-9012"] = (
    "Metabolism",
    ["GAPDH", "ACTB"],
    "Metabolic Processes",
  )
  reactome_db["R-HSA-3456"] = (
    "Signal Transduction",
    ["MYC"],
    "Signal Transduction",
  )

  let reactome_result = @bio.enrich_reactome(
    reactome_gene_list, reactome_db, reactome_bg_genes,
  )

  println("1. Pathway Enrichment Results")
  println(
    "   Total enriched pathways: \{reactome_result.pathways.length().to_string()}",
  )
  println("   Query genes: \{reactome_result.query_count.to_string()}")
  println(
    "   Background genes: \{reactome_result.background_count.to_string()}",
  )
  println("")

  println("2. Top 3 Pathways")
  let top_pathways = @bio.get_top_pathways(reactome_result, 3)
  i = 0
  while i < top_pathways.length() {
    let pathway = top_pathways[i]
    println("   Pathway \{i + 1}:")
    println("     - ID: \{pathway.pathway_id}")
    println("     - Name: \{pathway.pathway_name}")
    println("     - Gene count: \{pathway.gene_count.to_string()}")
    println("     - P-value: \{pathway.p_value.to_string()}")
    println("     - Adjusted P-value: \{pathway.adjusted_p_value.to_string()}")
    println("     - Enrichment score: \{pathway.enrichment_score.to_string()}")
    println("     - Top level term: \{pathway.top_level_term}")
    println("")
    i = i + 1
  }

  println("3. Pathways Grouped by Top Level Term")
  let groups = @bio.group_by_top_level(reactome_result)
  let group_keys = groups.keys().collect()
  j = 0
  while j < group_keys.length() {
    let group_key = group_keys[j]
    println(
      "   \{group_key}: \{groups[group_key].length().to_string()} pathways",
    )
    j = j + 1
  }
  println("")

  println("4. Filtered Results (adjusted p-value < 0.05)")
  let reactome_filtered = @bio.reactome_filter_results(reactome_result, 0.05)
  println(
    "   Significant pathways: \{reactome_filtered.pathways.length().to_string()}",
  )
  println("")

  println("5. Summary")
  println(@bio.reactome_summarize_result(reactome_result))

  println("\n")
  println("=== AnnotationDbi Gene Annotation Database Example ===")
  println("Inspired by Bioconductor AnnotationDbi package")
  println("")

  let db = @bio.create_annotation_db("org.Hs.eg.db")

  let gene_info1 = @bio.GeneInfo::new(
    "ENSG00000130203", "TP53", "Tumor protein p53", "chr17", 7661779, 7687550, "+",
    "protein_coding",
  )
  let gene_info2 = @bio.GeneInfo::new(
    "ENSG00000012048", "BRCA1", "Breast cancer type 1 susceptibility protein", "chr17",
    43044295, 43125370, "-", "protein_coding",
  )

  @bio.add_gene_info(db, gene_info1)
  @bio.add_gene_info(db, gene_info2)

  println("1. Database Summary")
  println(@bio.summarize_db(db))
  println("")

  println("2. Gene Info Query")
  let tp53_info = @bio.annotation_dbi_get_gene_info(db, "ENSG00000130203")
  if tp53_info.is_some() {
    println("   Found TP53: \{tp53_info.unwrap().symbol}")
  } else {
    println("   TP53 not found")
  }
  println("")

  println("3. ID Mapping")
  let id_mappings = Map(
    [("TP53", "ENSG00000130203"), ("BRCA1", "ENSG00000012048")],
    capacity=2,
  )
  @bio.add_id_map(db, "symbol", id_mappings)
  let mapped = @bio.map_ids(db, ["TP53", "BRCA1"], "symbol", "ensembl")
  println("   Mapped IDs: \{mapped.length().to_string()}")
  println("")

  println("4. Select Columns")
  let selected = @bio.select(db, ["ENSG00000130203"], [
    "gene_id", "symbol", "chr",
  ])
  println("   Selected rows: \{selected.length().to_string()}")
  println("")

  println("5. GO and KEGG Annotations")
  let go_terms = [
    @bio.AnnotationGOTerm::new("GO:0005515", "protein binding", "IEA", "MF"),
  ]
  let kegg_terms = [
    @bio.AnnotationKEGGTerm::new("hsa04115", "p53 signaling pathway", "Pathway"),
  ]
  @bio.add_go_terms(db, "ENSG00000130203", go_terms)
  @bio.add_kegg_terms(db, "ENSG00000130203", kegg_terms)
  let go_annotations = @bio.get_go_annotations(db, ["ENSG00000130203"])
  let kegg_annotations = @bio.get_kegg_annotations(db, ["ENSG00000130203"])
  println("   GO annotations: \{go_annotations.length().to_string()}")
  println("   KEGG annotations: \{kegg_annotations.length().to_string()}")
  println("")

  println("\n")
  println("=== ClusterProfiler Generic Enrichment Analysis Example ===")
  println("Inspired by Bioconductor clusterProfiler package")
  println("")

  let cp_gene_list = ["TP53", "BRCA1", "BRCA2", "ATM"]
  let cp_background = ["TP53", "BRCA1", "BRCA2", "ATM", "GAPDH", "ACTB", "MYC"]

  let cp_term_db = Map([], capacity=3)
  cp_term_db["GO:0005515"] = ("protein binding", ["TP53", "BRCA1", "BRCA2"])
  cp_term_db["GO:0007283"] = ("spermatogenesis", ["BRCA2"])
  cp_term_db["GO:0046982"] = ("protein heterodimerization", ["TP53", "BRCA1"])

  let cp_result = @bio.cp_enrich_generic(
    cp_gene_list, cp_term_db, cp_background, "Homo sapiens", "GO",
  )

  println("1. Generic Enrichment Results")
  println("   Total terms: \{cp_result.terms.length().to_string()}")
  println("   Query genes: \{cp_result.query_count.to_string()}")
  println("")

  println("2. Top 2 Terms")
  let cp_top = @bio.cp_get_top_terms(cp_result, 2)
  i = 0
  while i < cp_top.length() {
    let term = cp_top[i]
    println("   Term \{i + 1}:")
    println("     - Name: \{term.term_name}")
    println("     - P-value: \{term.p_value.to_string()}")
    println("     - Gene ratio: \{term.gene_ratio.to_string()}")
    println("")
    i = i + 1
  }

  println("3. Filtered Results")
  let cp_filtered = @bio.cp_filter_results(cp_result, 0.05)
  println("   Significant terms: \{cp_filtered.terms.length().to_string()}")
  println("")

  println("4. Dot Plot")
  println(@bio.cp_dotplot(cp_result, 3))
  println("")

  println("5. Enrichment Map")
  println(@bio.cp_enrich_map(cp_result))
  println("")

  println("\n")
  println("=== WGCNA Weighted Gene Co-expression Network Analysis Example ===")
  println("Inspired by Bioconductor WGCNA package")
  println("")

  let wgcna_gene_names = ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6"]
  let wgcna_sample_names = ["sample1", "sample2", "sample3", "sample4"]
  let wgcna_expression = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [10.0, 11.0, 12.0, 13.0],
    [11.0, 12.0, 13.0, 14.0],
    [12.0, 13.0, 14.0, 15.0],
  ]

  let wgcna_data = @bio.create_wgcna_data(
    wgcna_expression, wgcna_gene_names, wgcna_sample_names,
  )

  println("1. Create WGCNA Data")
  println("   Genes: \{wgcna_data.gene_names.length().to_string()}")
  println("   Samples: \{wgcna_data.sample_names.length().to_string()}")
  println("")

  println("2. Pick Soft Threshold")
  let powers = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let best_power = @bio.pick_soft_threshold(wgcna_data, powers)
  println("   Best power: \{best_power.to_string()}")
  println("")

  println("3. Build Adjacency Matrix")
  let wgcna_with_adj = @bio.adjacency(wgcna_data, best_power)
  println(
    "   Adjacency matrix size: \{wgcna_with_adj.adjacency_matrix.length().to_string()}x\{wgcna_with_adj.adjacency_matrix.length().to_string()}",
  )
  println("")

  println("4. Calculate TOM Similarity")
  let wgcna_with_tom = @bio.tom_similarity(wgcna_with_adj)
  println(
    "   TOM matrix size: \{wgcna_with_tom.tom_matrix.length().to_string()}x\{wgcna_with_tom.tom_matrix.length().to_string()}",
  )
  println("")

  println("5. Detect Modules")
  let wgcna_with_modules = @bio.blockwise_modules(wgcna_with_tom, best_power, 2)
  println(
    "   Number of modules: \{wgcna_with_modules.modules.length().to_string()}",
  )
  println("")

  let module_names = wgcna_with_modules.modules.keys().collect()
  j = 0
  while j < module_names.length() {
    let name = module_names[j]
    let color = if wgcna_with_modules.module_colors.contains(name) {
      wgcna_with_modules.module_colors[name]
    } else {
      "gray"
    }
    println(
      "   \{name} (\{color}): \{wgcna_with_modules.modules[name].length().to_string()} genes",
    )
    j = j + 1
  }
  println("")

  println("6. Module Eigengenes")
  let eigengenes = @bio.module_eigengenes(wgcna_with_modules)
  println("   Number of eigengenes: \{eigengenes.length().to_string()}")
  println("")

  println("7. Summary")
  println(@bio.summarize_wgcna(wgcna_with_modules))
}
