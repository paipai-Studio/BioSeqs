///|
/// Metadata-aware transcript quantification import, ported from
/// Bioconductor's `tximeta` package.
///
/// The reference package wraps `tximport`: it reads Salmon/Kallisto
/// quantifications, automatically looks up the transcriptome used for
/// quantification, links it to a known source (`linkTxome`) and
/// injects genome/organism/release metadata into the resulting
/// SummarizedExperiment.  This port keeps that pipeline over the
/// existing `TxImportData`-style matrices.

///|
/// Transcriptome provenance metadata attached by `link_txome`.
pub struct TximetaTxome {
  sha : String
  genome : String
  organism : String
  release : String
  source : String
   versioned_tags : Array[String]
} derive(Eq, Show, Debug)
pub extend TximetaTxome with Eq::{equal, not_equal}
pub extend TximetaTxome with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn TximetaTxome::new(
  sha~ : String,
  genome~ : String,
  organism~ : String,
  release~ : String,
  source~ : String,
  versioned_tags~ : Array[String] = [],
) -> TximetaTxome {
  TximetaTxome::{ sha, genome, organism, release, source, versioned_tags }
}

///|
/// Quantification data plus injected transcriptome metadata.
pub struct TximetaData {
  sample_names : Array[String]
  transcript_names : Array[String]
  counts : Array[Array[Int]]
  abundance : Array[Array[Double]]
  coldata : Map[String, Map[String, String]]
  txome : TximetaTxome?
   linked : Bool
} derive(Eq, Show, Debug)
pub extend TximetaData with Eq::{equal, not_equal}
pub extend TximetaData with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn TximetaData::new(
  sample_names~ : Array[String],
  transcript_names~ : Array[String],
  counts~ : Array[Array[Int]],
  abundance~ : Array[Array[Double]],
  coldata~ : Map[String, Map[String, String]] = {
    let m : Map[String, Map[String, String]] = Map::new()
    m
  },
) -> TximetaData {
  TximetaData::{
    sample_names,
    transcript_names,
    counts,
    abundance,
    coldata,
    txome: None,
    linked: false,
  }
}

///|
/// Number of samples.
pub fn TximetaData::n_samples(self : TximetaData) -> Int {
  self.sample_names.length()
}

///|
/// Number of transcripts.
pub fn TximetaData::n_transcripts(self : TximetaData) -> Int {
  self.transcript_names.length()
}

///|
/// Link the quantification to a known transcriptome
/// (`linkTxome`).  Returns a new linked data object.
pub fn link_txome(
  data : TximetaData,
  txome : TximetaTxome,
) -> TximetaData {
  TximetaData::{
    sample_names: data.sample_names,
    transcript_names: data.transcript_names,
    counts: data.counts,
    abundance: data.abundance,
    coldata: data.coldata,
    txome: Some(txome),
    linked: true,
  }
}

///|
/// Validate that the data is linked and carries a non-empty sha.
pub fn tximeta_validate_linked(data : TximetaData) -> Bool {
  match data.txome {
    Some(t) => data.linked && t.sha != ""
    None => false
  }
}

///|
/// Merged metadata: coldata keys plus transcriptome provenance
/// (`tximeta`-style metadata slot).
pub fn tximeta_metadata(data : TximetaData) -> Map[String, String] {
  let meta : Map[String, String] = Map::new()
  for sample in data.coldata.keys() {
    match data.coldata.get(sample) {
      Some(cols) => {
        for k in cols.keys() {
          match cols.get(k) {
            Some(v) => meta.set(sample + "." + k, v)
            None => ()
          }
        }
      }
      None => ()
    }
  }
  match data.txome {
    Some(t) => {
      meta.set("txome.sha", t.sha)
      meta.set("txome.genome", t.genome)
      meta.set("txome.organism", t.organism)
      meta.set("txome.release", t.release)
      meta.set("txome.source", t.source)
      let mut idx = 0
      while idx < t.versioned_tags.length() {
        meta.set(
          "txome.tag." + idx.to_string(),
          t.versioned_tags[idx],
        )
        idx += 1
      }
    }
    None => ()
  }
  meta
}

///|
/// Summarize transcript-level counts and abundance to the gene level
/// (counts summed, abundance summed within gene).
///
/// Returns `(gene_names, gene_counts, gene_abundance)`.
pub fn tximeta_summarize_to_gene(
  data : TximetaData,
  tx2gene~ : Map[String, String],
) -> (Array[String], Array[Array[Int]], Array[Array[Double]]) {
  let gene_order : Array[String] = Array::new()
  let gene_index : Map[String, Int] = Map::new()
  for i in 0..<data.transcript_names.length() {
    let tx = data.transcript_names[i]
    match tx2gene.get(tx) {
      Some(g) => {
        if gene_index.get(g) is None {
          gene_index.set(g, gene_order.length())
          gene_order.push(g)
        }
      }
      None => ()
    }
  }
  let n_genes = gene_order.length()
  let n_samples = data.sample_names.length()
  let gene_counts : Array[Array[Int]] = Array::new()
  let gene_abundance : Array[Array[Double]] = Array::new()
  for _g in 0..<n_genes {
    gene_counts.push(Array::make(n_samples, 0))
    gene_abundance.push(Array::make(n_samples, 0.0))
  }
  for i in 0..<data.transcript_names.length() {
    let tx = data.transcript_names[i]
    match tx2gene.get(tx) {
      Some(g) => {
        match gene_index.get(g) {
          Some(row) => {
            let row_counts = gene_counts[row]
            let row_abund = gene_abundance[row]
            let mut j = 0
            while j < n_samples {
              let new_val = row_counts[j] + data.counts[i][j]
              let new_ab = row_abund[j] + data.abundance[i][j]
              row_counts.set(j, new_val)
              row_abund.set(j, new_ab)
              j += 1
            }
            gene_counts[row] = row_counts
            gene_abundance[row] = row_abund
          }
          None => ()
        }
      }
      None => ()
    }
  }
  (gene_order, gene_counts, gene_abundance)
}

///|
/// Human-readable summary of the imported data.
pub fn tximeta_data_summary(data : TximetaData) -> String {
  let mut out = "TximetaData: " +
  data.n_samples().to_string() +
  " samples x " +
  data.n_transcripts().to_string() +
  " transcripts"
  match data.txome {
    Some(t) => {
      out = out +
      " | txome: " +
      t.organism +
      " " +
      t.genome +
      " (release " +
      t.release +
      ", sha " +
      t.sha +
      ")"
    }
    None => out = out + " | txome: unlinked"
  }
  out
}

///|
/// A small sample dataset with a linked transcriptome.
pub fn tximeta_sample_data() -> TximetaData {
  let samples = ["ctrl1", "ctrl2", "treated1", "treated2"]
  let transcripts = ["ENST1", "ENST2", "ENST3"]
  let counts = [
    [10, 12, 40, 45],
    [100, 110, 95, 90],
    [5, 4, 20, 22],
  ]
  let abundance = [
    [1.0, 1.2, 4.0, 4.5],
    [10.0, 11.0, 9.5, 9.0],
    [0.5, 0.4, 2.0, 2.2],
  ]
  let coldata : Map[String, Map[String, String]] = Map::new()
  let g1 : Map[String, String] = Map::new()
  g1.set("condition", "control")
  coldata.set("ctrl1", g1)
  let g2 : Map[String, String] = Map::new()
  g2.set("condition", "control")
  coldata.set("ctrl2", g2)
  let g3 : Map[String, String] = Map::new()
  g3.set("condition", "treated")
  coldata.set("treated1", g3)
  let g4 : Map[String, String] = Map::new()
  g4.set("condition", "treated")
  coldata.set("treated2", g4)
  let data = TximetaData::new(
    sample_names=samples,
    transcript_names=transcripts,
    counts=counts,
    abundance=abundance,
    coldata=coldata,
  )
  link_txome(
    data,
    TximetaTxome::new(
      sha="a1b2c3d4",
      genome="GRCh38",
      organism="Homo sapiens",
      release="110",
      source="GENCODE",
      versioned_tags=["v38"],
    ),
  )
}

///|
/// Sample tx2gene mapping matching the sample dataset.
pub fn tximeta_sample_tx2gene() -> Map[String, String] {
  let m : Map[String, String] = Map::new()
  m.set("ENST1", "GENE_A")
  m.set("ENST2", "GENE_B")
  m.set("ENST3", "GENE_A")
  m
}