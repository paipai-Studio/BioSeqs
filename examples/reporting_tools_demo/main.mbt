///|
/// ReportingTools module example - Bioinformatics report generation.
///
/// Run: moon run examples/reporting_tools_demo
fn main {
  println("=== ReportingTools Module Demo ===")
  println("")

  // Create a new report
  println("=== Report Creation ===")
  let doc = @bio.ReportDocument::new("RNA-Seq Differential Expression Analysis")
  println("Title: \{doc.get_title()}")
  println("Author: \{doc.get_author()}")
  println("Sections: \{doc.get_n_sections()}")
  println("")

  // Set author
  let doc2 = doc.set_author("BioSeqs Analysis Pipeline")
  println("Updated author: \{doc2.get_author()}")
  println("")

  // Add text sections
  println("=== Adding Text Sections ===")
  let doc3 = doc2.add_text("intro", "Introduction", "This report summarizes the results of an RNA-Seq differential expression analysis comparing two conditions: treatment and control.")
  let doc4 = doc3.add_text("methods", "Methods", "Total RNA was extracted and sequenced using Illumina NovaSeq. Reads were aligned to GRCh38 using HISAT2, quantified with featureCounts, and differential expression analyzed with DESeq2.")
  let doc5 = doc4.add_text("sample_info", "Sample Information", "A total of 20 samples were analyzed (10 treatment, 10 control).")
  println("Sections after adding text: \{doc5.get_n_sections()}")
  println("")

  // Add a results table
  println("=== Adding Results Table ===")
  let table = @bio.ReportTable::new("de_results", "Top Differential Genes",
    ["Gene", "log2FC", "p-value", "padj", "Regulation"],
    [["TP53", "BRCA1", "MYC", "KRAS", "EGFR"],
     ["2.54", "-3.02", "3.21", "-1.48", "0.87"],
     ["0.0012", "0.0035", "0.0008", "0.0125", "0.0620"],
     ["0.012", "0.035", "0.009", "0.085", "0.215"],
     ["Up", "Down", "Up", "Down", "NS"]],
    "Table 1. Top 5 differentially expressed genes"
  )
  let doc6 = doc5.add_table("results_table", table)
  println("Table added. Total tables: \{doc6.get_n_tables()}")
  println("Table rows: \{table.get_n_rows()}")
  println("Table columns: \{table.get_n_columns()}")
  println("")

  // Render table as ASCII
  println("=== Table Rendering ===")
  println(table.to_ascii())
  println("")

  // Add a volcano plot
  println("=== Adding Plot ===")
  let plot = @bio.volcano_plot_default(
    ["TP53", "BRCA1", "MYC", "KRAS", "EGFR", "PTEN", "APC", "RB1"],
    [2.5, -3.0, 3.2, -1.5, 0.8, -2.2, 0.3, -1.8],
    [0.001, 0.005, 0.0005, 0.01, 0.08, 0.002, 0.06, 0.003]
  )
  let plot_ascii = plot.to_ascii(width=50, height=12)
  let doc7 = doc6.add_plot("volcano", "Volcano Plot of DE Results", plot_ascii)
  println("Plot section added.")
  println("")

  // Full report
  println("=== Full Report ===")
  let rendered = doc7.render()
  println(rendered)

  // Summary
  println("=== Report Summary ===")
  println(doc7.summary())

  println("=== Done ===")
}