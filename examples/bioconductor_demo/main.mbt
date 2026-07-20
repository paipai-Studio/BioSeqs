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
  println("")

  println("\n")
  println("=== Batchelor Single-cell Batch Correction Example ===")
  println("Inspired by Bioconductor batchelor package")
  println("")

  println("1. Generate simulated single-cell data with batch effects")
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 50, 200, 2.0)
  let mut idx = 0
  while idx < batches.length() {
    println("   \{batch_names[idx]}: \{batches[idx].length().to_string()} cells, \{batches[idx][0].length().to_string()} genes")
    idx = idx + 1
  }
  println("")

  println("2. Evaluate batch mixing before correction")
  let raw_data : Array[Array[Double]] = Array::new()
  for cell in batches[0] {
    raw_data.push(cell)
  }
  for cell in batches[1] {
    raw_data.push(cell)
  }
  let raw_batch_indices : Array[Int] = Array::new()
  let mut bi = 0
  while bi < 50 {
    raw_batch_indices.push(0)
    bi = bi + 1
  }
  bi = 0
  while bi < 50 {
    raw_batch_indices.push(1)
    bi = bi + 1
  }
  let raw_score = @bio.compute_batch_mixing_score(raw_data, raw_batch_indices)
  println("   Batch mixing score (higher = better): \{raw_score.to_string()}")
  println("")

  println("3. Method 1: rescaleBatches")
  let rescale_result = @bio.rescale_batches(batches, batch_names)
  let rescale_score = @bio.compute_batch_mixing_score(rescale_result.corrected, rescale_result.batch_indices)
  println("   Corrected cells: \{rescale_result.corrected.length().to_string()}")
  println("   Batch mixing score: \{rescale_score.to_string()}")
  println("   Improvement: +\{((rescale_score - raw_score) / raw_score * 100.0).to_string()}%")
  println("")

  println("4. Method 2: fastMNN")
  let fastmnn_result = @bio.fast_mnn(batches, batch_names, 5, 20)
  let fastmnn_score = @bio.compute_batch_mixing_score(fastmnn_result.corrected, fastmnn_result.batch_indices)
  println("   Corrected dimensions: \{fastmnn_result.corrected.length().to_string()} cells, \{fastmnn_result.corrected[0].length().to_string()} PCs")
  println("   Number of MNN pairs: \{fastmnn_result.mnn_pairs.length().to_string()}")
  println("   Batch mixing score: \{fastmnn_score.to_string()}")
  println("   Improvement: +\{((fastmnn_score - raw_score) / raw_score * 100.0).to_string()}%")
  println("")

  println("5. Find Mutual Nearest Neighbors")
  let mnn_pairs = @bio.find_mutual_nearest_neighbors(batches[0], batches[1], 5)
  println("   MNN pairs found: \{mnn_pairs.length().to_string()}")
  if mnn_pairs.length() > 0 {
    println("   First 3 pairs:")
    let mut pi = 0
    while pi < 3 && pi < mnn_pairs.length() {
      let pair = mnn_pairs[pi]
      println("     (Batch_0 cell \{pair.0.to_string()}, Batch_1 cell \{pair.1.to_string()})")
      pi = pi + 1
    }
  }
  println("")

  println("6. Comparison")
  println("   Method           | Score")
  println("   -----------------|------")
  println("   Raw data         | \{raw_score.to_string()}")
  println("   rescaleBatches   | \{rescale_score.to_string()}")
  println("   fastMNN          | \{fastmnn_score.to_string()}")
  println("")

  println("\n")
  println("=== Seurat Single-cell RNA-seq Analysis Example ===")
  println("Inspired by Bioconductor Seurat package")
  println("")

  println("1. Create simulated single-cell data")
  let seurat_obj = @bio.seurat_create_example_data()
  println("   Cells: \{seurat_obj.col_names.length().to_string()}")
  println("   Genes: \{seurat_obj.row_names.length().to_string()}")
  println("")

  println("2. Normalize data")
  let seurat_normalized = @bio.normalize_total(seurat_obj)
  println("   Normalization completed")
  println("")

  println("3. Find highly variable genes")
  let seurat_hvg = @bio.find_variable_features(seurat_normalized)
  println("   HVGs found: \{seurat_hvg.var_features.length().to_string()}")
  println("")

  println("4. Scale data")
  let seurat_scaled = @bio.scale_data(seurat_hvg)
  println("   Scaling completed")
  println("")

  println("5. Run PCA")
  let seurat_pca = @bio.run_pca(seurat_scaled)
  println("   Number of PCs: \{seurat_pca.pca[0].length().to_string()}")
  let mut pc_idx = 0
  while pc_idx < 3 && pc_idx < seurat_pca.var_explained.length() {
    println("   PC\{(pc_idx + 1).to_string()}: \{(seurat_pca.var_explained[pc_idx] * 100.0).to_string()}% variance")
    pc_idx = pc_idx + 1
  }
  println("")

  println("6. Find nearest neighbors")
  let seurat_neighbors = @bio.find_neighbors(seurat_pca)
  println("   Neighbors per cell: \{seurat_neighbors.neighbors[0].length().to_string()}")
  println("")

  println("7. Find clusters")
  let seurat_clusters = @bio.find_clusters(seurat_neighbors)
  
  let unique_clusters : Array[Int] = Array::new()
  for c in seurat_clusters.clusters {
    let mut exists = false
    for u in unique_clusters {
      if u == c {
        exists = true
        break
      }
    }
    if !exists {
      unique_clusters.push(c)
    }
  }
  println("   Number of clusters: \{unique_clusters.length().to_string()}")
  println("")

  println("8. Run UMAP")
  let seurat_umap = @bio.run_umap(seurat_pca)
  println("   UMAP completed with 2 dimensions")
  println("")

  println("9. Find marker genes")
  let markers = @bio.find_all_markers(seurat_clusters)
  println("   Markers found: \{markers.length().to_string()}")
  
  if markers.length() > 0 {
    let top_markers = @bio.get_top_markers(markers, 3)
    println("   Top 3 markers:")
    let mut mi = 0
    while mi < top_markers.length() {
      let marker = top_markers[mi]
      println("     \{marker.gene}: log2FC=\{marker.avg_log2_fc.to_string()}, padj=\{marker.adjusted_p_value.to_string()}")
      mi = mi + 1
    }
  }
  println("")

  println("=== All examples completed! ===")
}
