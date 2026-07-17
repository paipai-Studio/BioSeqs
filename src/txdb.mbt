///|
/// TxDb - Transcript Database inspired by Bioconductor GenomicFeatures.

///|
pub struct Gene {
  gene_id : String
  seqname : String
  start : Int
  end : Int
  strand : Strand
}

///|
pub fn Gene::new(
  gene_id : String,
  seqname : String,
  start : Int,
  end : Int,
  strand : Strand,
) -> Gene {
  Gene::{ gene_id, seqname, start, end, strand }
}

///|
pub struct Transcript {
  transcript_id : String
  gene_id : String
  seqname : String
  start : Int
  end : Int
  strand : Strand
}

///|
pub fn Transcript::new(
  transcript_id : String,
  gene_id : String,
  seqname : String,
  start : Int,
  end : Int,
  strand : Strand,
) -> Transcript {
  Transcript::{ transcript_id, gene_id, seqname, start, end, strand }
}

///|
pub struct Exon {
  exon_id : String
  transcript_id : String
  gene_id : String
  seqname : String
  start : Int
  end : Int
  strand : Strand
  exon_rank : Int
}

///|
pub fn Exon::new(
  exon_id : String,
  transcript_id : String,
  gene_id : String,
  seqname : String,
  start : Int,
  end : Int,
  strand : Strand,
  exon_rank : Int,
) -> Exon {
  Exon::{
    exon_id,
    transcript_id,
    gene_id,
    seqname,
    start,
    end,
    strand,
    exon_rank,
  }
}

///|
pub struct CDS {
  cds_id : String
  transcript_id : String
  gene_id : String
  seqname : String
  start : Int
  end : Int
  strand : Strand
  exon_rank : Int
}

///|
pub fn CDS::new(
  cds_id : String,
  transcript_id : String,
  gene_id : String,
  seqname : String,
  start : Int,
  end : Int,
  strand : Strand,
  exon_rank : Int,
) -> CDS {
  CDS::{
    cds_id,
    transcript_id,
    gene_id,
    seqname,
    start,
    end,
    strand,
    exon_rank,
  }
}

///|
pub struct TxDb {
  genes : Array[Gene]
  transcripts : Array[Transcript]
  exons : Array[Exon]
  cds : Array[CDS]
  genome : String
  organism : String
}

///|
pub fn TxDb::new(
  genes : Array[Gene],
  transcripts : Array[Transcript],
  exons : Array[Exon],
  cds : Array[CDS],
  genome : String,
  organism : String,
) -> TxDb {
  TxDb::{ genes, transcripts, exons, cds, genome, organism }
}

///|
pub fn txdb(
  genes : Array[Gene],
  transcripts : Array[Transcript],
  exons : Array[Exon],
  cds : Array[CDS],
  genome : String,
  organism : String,
) -> TxDb {
  TxDb::{ genes, transcripts, exons, cds, genome, organism }
}

///|
pub fn genes(txdb : TxDb) -> GRanges {
  let n = txdb.genes.length()
  let seqnames : Array[String] = Array::make(n, "")
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let strands : Array[Strand] = Array::make(n, Strand::Star)

  let mut i = 0
  while i < n {
    let g = txdb.genes[i]
    seqnames[i] = g.seqname
    ranges[i] = (g.start, g.end)
    strands[i] = g.strand
    i = i + 1
  }

  granges(seqnames, ranges, strands)
}

///|
pub fn transcripts(txdb : TxDb) -> GRanges {
  let n = txdb.transcripts.length()
  let seqnames : Array[String] = Array::make(n, "")
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let strands : Array[Strand] = Array::make(n, Strand::Star)

  let mut i = 0
  while i < n {
    let tx = txdb.transcripts[i]
    seqnames[i] = tx.seqname
    ranges[i] = (tx.start, tx.end)
    strands[i] = tx.strand
    i = i + 1
  }

  granges(seqnames, ranges, strands)
}

///|
pub fn exons(txdb : TxDb) -> GRanges {
  let n = txdb.exons.length()
  let seqnames : Array[String] = Array::make(n, "")
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let strands : Array[Strand] = Array::make(n, Strand::Star)

  let mut i = 0
  while i < n {
    let e = txdb.exons[i]
    seqnames[i] = e.seqname
    ranges[i] = (e.start, e.end)
    strands[i] = e.strand
    i = i + 1
  }

  granges(seqnames, ranges, strands)
}

///|
pub fn cds(txdb : TxDb) -> GRanges {
  let n = txdb.cds.length()
  let seqnames : Array[String] = Array::make(n, "")
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let strands : Array[Strand] = Array::make(n, Strand::Star)

  let mut i = 0
  while i < n {
    let c = txdb.cds[i]
    seqnames[i] = c.seqname
    ranges[i] = (c.start, c.end)
    strands[i] = c.strand
    i = i + 1
  }

  granges(seqnames, ranges, strands)
}

///|
pub fn transcripts_by_gene(txdb : TxDb) -> Map[String, Array[Transcript]] {
  let result = Map([], capacity=32)

  let mut i = 0
  while i < txdb.transcripts.length() {
    let tx = txdb.transcripts[i]
    if !result.contains(tx.gene_id) {
      result[tx.gene_id] = []
    }
    result[tx.gene_id].push(tx)
    i = i + 1
  }

  result
}

///|
pub fn exons_by_transcript(txdb : TxDb) -> Map[String, Array[Exon]] {
  let result = Map([], capacity=32)

  let mut i = 0
  while i < txdb.exons.length() {
    let exon = txdb.exons[i]
    if !result.contains(exon.transcript_id) {
      result[exon.transcript_id] = []
    }
    result[exon.transcript_id].push(exon)
    i = i + 1
  }

  result
}

///|
pub fn exons_by_gene(txdb : TxDb) -> Map[String, Array[Exon]] {
  let result = Map([], capacity=32)

  let mut i = 0
  while i < txdb.exons.length() {
    let exon = txdb.exons[i]
    if !result.contains(exon.gene_id) {
      result[exon.gene_id] = []
    }
    result[exon.gene_id].push(exon)
    i = i + 1
  }

  result
}

///|
pub fn cds_by_transcript(txdb : TxDb) -> Map[String, Array[CDS]] {
  let result = Map([], capacity=32)

  let mut i = 0
  while i < txdb.cds.length() {
    let c = txdb.cds[i]
    if !result.contains(c.transcript_id) {
      result[c.transcript_id] = []
    }
    result[c.transcript_id].push(c)
    i = i + 1
  }

  result
}

///|
pub fn introns_by_transcript(txdb : TxDb) -> Map[String, Array[(Int, Int)]] {
  let result = Map([], capacity=32)
  let exons_map = exons_by_transcript(txdb)
  let keys = exons_map.keys().collect()

  let mut i = 0
  while i < keys.length() {
    let tx_id = keys[i]
    let exons = exons_map[tx_id]
    if exons.length() <= 1 {
      result[tx_id] = []
    } else {
      let n = exons.length()
      let sorted : Array[Exon] = Array::make(n, exons[0])
      let mut j = 0
      while j < n {
        sorted[j] = exons[j]
        j = j + 1
      }

      let mut k = 0
      while k < n - 1 {
        let mut l = k + 1
        while l < n {
          if sorted[k].exon_rank > sorted[l].exon_rank {
            let temp = sorted[k]
            sorted[k] = sorted[l]
            sorted[l] = temp
          }
          l = l + 1
        }
        k = k + 1
      }

      let introns : Array[(Int, Int)] = []
      let mut m = 0
      while m < n - 1 {
        let start = sorted[m].end + 1
        let end = sorted[m + 1].start - 1
        if start <= end {
          introns.push((start, end))
        }
        m = m + 1
      }
      result[tx_id] = introns
    }
    i = i + 1
  }

  result
}

///|
pub fn find_transcript_by_id(txdb : TxDb, tx_id : String) -> Transcript? {
  let mut i = 0
  while i < txdb.transcripts.length() {
    if txdb.transcripts[i].transcript_id == tx_id {
      return Some(txdb.transcripts[i])
    }
    i = i + 1
  }
  None
}

///|
pub fn find_min_cds_start(cds_regions : Array[CDS]) -> Int {
  if cds_regions.length() == 0 {
    return 0
  }
  let mut min_start = cds_regions[0].start
  let mut i = 1
  while i < cds_regions.length() {
    if cds_regions[i].start < min_start {
      min_start = cds_regions[i].start
    }
    i = i + 1
  }
  min_start
}

///|
pub fn find_max_cds_end(cds_regions : Array[CDS]) -> Int {
  if cds_regions.length() == 0 {
    return 0
  }
  let mut max_end = cds_regions[0].end
  let mut i = 1
  while i < cds_regions.length() {
    if cds_regions[i].end > max_end {
      max_end = cds_regions[i].end
    }
    i = i + 1
  }
  max_end
}

///|
pub fn sort_exons_by_rank(exons : Array[Exon]) -> Array[Exon] {
  let n = exons.length()
  let sorted : Array[Exon] = Array::make(n, exons[0])
  let mut i = 0
  while i < n {
    sorted[i] = exons[i]
    i = i + 1
  }

  let mut j = 0
  while j < n - 1 {
    let mut k = j + 1
    while k < n {
      if sorted[j].exon_rank > sorted[k].exon_rank {
        let temp = sorted[j]
        sorted[j] = sorted[k]
        sorted[k] = temp
      }
      k = k + 1
    }
    j = j + 1
  }

  sorted
}

///|
pub fn five_prime_utrs(txdb : TxDb) -> GRanges {
  let results : GRanges = GRanges::{
    seqnames: [],
    starts: [],
    ends: [],
    strands: [],
    widths: [],
  }

  let exons_map = exons_by_transcript(txdb)
  let cds_map = cds_by_transcript(txdb)
  let keys = exons_map.keys().collect()

  let mut i = 0
  while i < keys.length() {
    let tx_id = keys[i]
    let exons = exons_map[tx_id]

    if !cds_map.contains(tx_id) {
      i = i + 1
      continue
    }
    let cds_regions = cds_map[tx_id]
    if cds_regions.length() == 0 {
      i = i + 1
      continue
    }

    let tx_opt = find_transcript_by_id(txdb, tx_id)
    if tx_opt.is_none() {
      i = i + 1
      continue
    }
    let transcript = tx_opt.unwrap()

    if transcript.strand == Strand::Plus {
      let cds_start = find_min_cds_start(cds_regions)
      let sorted = sort_exons_by_rank(exons)

      let mut j = 0
      while j < sorted.length() {
        let e = sorted[j]
        if e.end < cds_start {
          results.seqnames.push(e.seqname)
          results.starts.push(e.start)
          results.ends.push(e.end)
          results.strands.push(e.strand)
          results.widths.push(e.end - e.start + 1)
        } else if e.start < cds_start && e.end >= cds_start {
          results.seqnames.push(e.seqname)
          results.starts.push(e.start)
          results.ends.push(cds_start - 1)
          results.strands.push(e.strand)
          results.widths.push(cds_start - e.start)
        }
        j = j + 1
      }
    } else {
      let cds_end = find_max_cds_end(cds_regions)
      let sorted = sort_exons_by_rank(exons)

      let mut j = 0
      while j < sorted.length() {
        let e = sorted[j]
        if e.start > cds_end {
          results.seqnames.push(e.seqname)
          results.starts.push(e.start)
          results.ends.push(e.end)
          results.strands.push(e.strand)
          results.widths.push(e.end - e.start + 1)
        } else if e.end > cds_end && e.start <= cds_end {
          results.seqnames.push(e.seqname)
          results.starts.push(cds_end + 1)
          results.ends.push(e.end)
          results.strands.push(e.strand)
          results.widths.push(e.end - cds_end)
        }
        j = j + 1
      }
    }

    i = i + 1
  }

  results
}

///|
pub fn three_prime_utrs(txdb : TxDb) -> GRanges {
  let results : GRanges = GRanges::{
    seqnames: [],
    starts: [],
    ends: [],
    strands: [],
    widths: [],
  }

  let exons_map = exons_by_transcript(txdb)
  let cds_map = cds_by_transcript(txdb)
  let keys = exons_map.keys().collect()

  let mut i = 0
  while i < keys.length() {
    let tx_id = keys[i]
    let exons = exons_map[tx_id]

    if !cds_map.contains(tx_id) {
      i = i + 1
      continue
    }
    let cds_regions = cds_map[tx_id]
    if cds_regions.length() == 0 {
      i = i + 1
      continue
    }

    let tx_opt = find_transcript_by_id(txdb, tx_id)
    if tx_opt.is_none() {
      i = i + 1
      continue
    }
    let transcript = tx_opt.unwrap()

    if transcript.strand == Strand::Plus {
      let cds_end = find_max_cds_end(cds_regions)
      let sorted = sort_exons_by_rank(exons)

      let mut j = 0
      while j < sorted.length() {
        let e = sorted[j]
        if e.start > cds_end {
          results.seqnames.push(e.seqname)
          results.starts.push(e.start)
          results.ends.push(e.end)
          results.strands.push(e.strand)
          results.widths.push(e.end - e.start + 1)
        } else if e.end > cds_end && e.start <= cds_end {
          results.seqnames.push(e.seqname)
          results.starts.push(cds_end + 1)
          results.ends.push(e.end)
          results.strands.push(e.strand)
          results.widths.push(e.end - cds_end)
        }
        j = j + 1
      }
    } else {
      let cds_start = find_min_cds_start(cds_regions)
      let sorted = sort_exons_by_rank(exons)

      let mut j = 0
      while j < sorted.length() {
        let e = sorted[j]
        if e.end < cds_start {
          results.seqnames.push(e.seqname)
          results.starts.push(e.start)
          results.ends.push(e.end)
          results.strands.push(e.strand)
          results.widths.push(e.end - e.start + 1)
        } else if e.start < cds_start && e.end >= cds_start {
          results.seqnames.push(e.seqname)
          results.starts.push(e.start)
          results.ends.push(cds_start - 1)
          results.strands.push(e.strand)
          results.widths.push(cds_start - e.start)
        }
        j = j + 1
      }
    }

    i = i + 1
  }

  results
}

///|
pub fn promoters(txdb : TxDb, upstream : Int, downstream : Int) -> GRanges {
  let txs = txdb.transcripts
  let n = txs.length()
  let seqnames : Array[String] = Array::make(n, "")
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let strands : Array[Strand] = Array::make(n, Strand::Star)

  let mut i = 0
  while i < n {
    let tx = txs[i]
    seqnames[i] = tx.seqname
    strands[i] = tx.strand

    if tx.strand == Strand::Plus {
      let p_start = tx.start - upstream
      let p_end = tx.start + downstream - 1
      ranges[i] = (if p_start < 1 { 1 } else { p_start }, p_end)
    } else {
      let p_start = tx.end - downstream + 1
      let p_end = tx.end + upstream
      ranges[i] = (p_start, p_end)
    }

    i = i + 1
  }

  granges(seqnames, ranges, strands)
}

///|
pub fn gene_ids(txdb : TxDb) -> Array[String] {
  let seen = Map([], capacity=32)
  let ids : Array[String] = []

  let mut i = 0
  while i < txdb.genes.length() {
    let gene = txdb.genes[i]
    if !seen.contains(gene.gene_id) {
      seen[gene.gene_id] = ()
      ids.push(gene.gene_id)
    }
    i = i + 1
  }

  ids
}

///|
pub fn transcript_ids(txdb : TxDb) -> Array[String] {
  let seen = Map([], capacity=32)
  let ids : Array[String] = []

  let mut i = 0
  while i < txdb.transcripts.length() {
    let tx = txdb.transcripts[i]
    if !seen.contains(tx.transcript_id) {
      seen[tx.transcript_id] = ()
      ids.push(tx.transcript_id)
    }
    i = i + 1
  }

  ids
}

///|
pub fn exon_ids(txdb : TxDb) -> Array[String] {
  let seen = Map([], capacity=32)
  let ids : Array[String] = []

  let mut i = 0
  while i < txdb.exons.length() {
    let exon = txdb.exons[i]
    if !seen.contains(exon.exon_id) {
      seen[exon.exon_id] = ()
      ids.push(exon.exon_id)
    }
    i = i + 1
  }

  ids
}

///|
pub fn txdb_summary(txdb : TxDb) -> String {
  let gene_count = gene_ids(txdb).length()
  let tx_count = transcript_ids(txdb).length()
  let exon_count = exon_ids(txdb).length()
  let cds_count = txdb.cds.length()

  "TxDb object:" +
  "\n" +
  "  Genome: " +
  txdb.genome +
  "\n" +
  "  Organism: " +
  txdb.organism +
  "\n" +
  "  Genes: " +
  gene_count.to_string() +
  "\n" +
  "  Transcripts: " +
  tx_count.to_string() +
  "\n" +
  "  Exons: " +
  exon_count.to_string() +
  "\n" +
  "  CDS: " +
  cds_count.to_string()
}

///|
pub fn parse_gtf(gtf_content : String) -> TxDb {
  let genes : Array[Gene] = []
  let transcripts : Array[Transcript] = []
  let exons : Array[Exon] = []
  let cds_regions : Array[CDS] = []

  let lines = split_by_char(gtf_content, '\n'.to_int())

  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    let trimmed = line.trim().to_owned()
    if trimmed.length() == 0 || starts_with(trimmed, 0, "#") {
      i = i + 1
      continue
    }

    let fields = split_by_char(trimmed, '\t'.to_int())
    if fields.length() < 9 {
      i = i + 1
      continue
    }

    let seqname = fields[0]
    let feature_type = fields[2]
    let start_opt = parse_double(fields[3])
    if start_opt.is_none() {
      i = i + 1
      continue
    }
    let start = start_opt.unwrap().to_int()

    let end_opt = parse_double(fields[4])
    if end_opt.is_none() {
      i = i + 1
      continue
    }
    let end = end_opt.unwrap().to_int()

    let strand_str = fields[6]
    let attributes = fields[8]

    let strand = match strand_str {
      "+" => Strand::Plus
      "-" => Strand::Minus
      _ => Strand::Star
    }

    let attr_map = parse_gtf_attributes(attributes)

    let gene_id = if attr_map.contains("gene_id") {
      let id = attr_map["gene_id"]
      id.replace(old="\"", new="").replace(old=";", new="")
    } else {
      "unknown"
    }

    let transcript_id = if attr_map.contains("transcript_id") {
      let id = attr_map["transcript_id"]
      id.replace(old="\"", new="").replace(old=";", new="")
    } else {
      gene_id
    }

    if feature_type == "gene" {
      genes.push(Gene::{ gene_id, seqname, start, end, strand })
    } else if feature_type == "transcript" {
      transcripts.push(Transcript::{
        transcript_id,
        gene_id,
        seqname,
        start,
        end,
        strand,
      })
    } else if feature_type == "exon" {
      let exon_rank = if attr_map.contains("exon_number") {
        let r = attr_map["exon_number"]
        let cleaned = r.replace(old="\"", new="").replace(old=";", new="")
        let d_opt = parse_double(cleaned)
        if d_opt.is_some() {
          d_opt.unwrap().to_int()
        } else {
          1
        }
      } else {
        1
      }

      exons.push(Exon::{
        exon_id: gene_id + "." + transcript_id + "." + exon_rank.to_string(),
        transcript_id,
        gene_id,
        seqname,
        start,
        end,
        strand,
        exon_rank,
      })
    } else if feature_type == "CDS" {
      let exon_rank = if attr_map.contains("exon_number") {
        let r = attr_map["exon_number"]
        let cleaned = r.replace(old="\"", new="").replace(old=";", new="")
        let d_opt = parse_double(cleaned)
        if d_opt.is_some() {
          d_opt.unwrap().to_int()
        } else {
          1
        }
      } else {
        1
      }

      cds_regions.push(CDS::{
        cds_id: gene_id + "." + transcript_id + ".cds." + exon_rank.to_string(),
        transcript_id,
        gene_id,
        seqname,
        start,
        end,
        strand,
        exon_rank,
      })
    }

    i = i + 1
  }

  TxDb::{
    genes,
    transcripts,
    exons,
    cds: cds_regions,
    genome: "unknown",
    organism: "unknown",
  }
}

///|
fn parse_gtf_attributes(attr_str : String) -> Map[String, String] {
  let result = Map([], capacity=16)

  let pairs = split_by_char(attr_str, ';'.to_int())
  let mut i = 0
  while i < pairs.length() {
    let pair = pairs[i]
    let trimmed = pair.trim().to_owned()
    if trimmed.length() == 0 {
      i = i + 1
      continue
    }

    let parts = split_by_char(trimmed, ' '.to_int())
    if parts.length() >= 2 {
      let key = parts[0]
      let value = if parts.length() > 2 {
        substring(
          trimmed,
          key.length() + 1,
          trimmed.length() - key.length() - 1,
        )
      } else {
        parts[1]
      }
      result[key] = value
    }

    i = i + 1
  }

  result
}
