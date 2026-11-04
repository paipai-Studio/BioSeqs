/// pysam_bam_transcript.mbt — Transcript-to-genome coordinate mapping
/// Inspired by pysam transcript coordinate handling

///|
/// A transcript exon.
pub struct PysamExon {
  chrom : String
  start : Int32 // 0-based genome start
  end : Int32 // exclusive genome end
  transcript_start : Int32 // 0-based transcript position
}

///|
/// A transcript with exons.
pub struct PysamTranscript {
  name : String
  chrom : String
  strand : Bool // true = forward
  exons : Array[PysamExon]
}

///|
/// Map transcript coordinate to genome coordinate.
pub fn pysam_transcript_to_genome(
  transcript : PysamTranscript,
  transcript_pos : Int32,
) -> Int32? {
  for i = 0; i < transcript.exons.length(); i = i + 1 {
    let exon = transcript.exons[i]
    let exon_len = exon.end - exon.start
    if transcript_pos >= exon.transcript_start &&
      transcript_pos < exon.transcript_start + exon_len {
      if transcript.strand {
        return Some(exon.start + (transcript_pos - exon.transcript_start))
      } else {
        return Some(exon.end - 1 - (transcript_pos - exon.transcript_start))
      }
    }
  }
  None
}

///|
/// Map genome coordinate to transcript coordinate.
pub fn pysam_genome_to_transcript(
  transcript : PysamTranscript,
  genome_pos : Int32,
) -> Int32? {
  for i = 0; i < transcript.exons.length(); i = i + 1 {
    let exon = transcript.exons[i]
    if genome_pos >= exon.start && genome_pos < exon.end {
      if transcript.strand {
        return Some(exon.transcript_start + (genome_pos - exon.start))
      } else {
        return Some(exon.transcript_start + (exon.end - 1 - genome_pos))
      }
    }
  }
  None
}

///|
/// Get transcript length (sum of exon lengths).
pub fn pysam_transcript_length(transcript : PysamTranscript) -> Int32 {
  let mut total = 0
  for i = 0; i < transcript.exons.length(); i = i + 1 {
    total = total + (transcript.exons[i].end - transcript.exons[i].start)
  }
  total
}

///|
/// Build transcript from exon coordinates.
pub fn pysam_build_transcript(
  name : String,
  chrom : String,
  strand : Bool,
  exon_coords : Array[(Int32, Int32)],
) -> PysamTranscript {
  // Sort exons by start position
  let sorted = exon_coords.copy()
  sorted.sort_by(fn(a, b) {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })
  let exons : Array[PysamExon] = Array(capacity=0)
  let mut transcript_pos = 0
  for i = 0; i < sorted.length(); i = i + 1 {
    let (start, end) = sorted[i]
    exons.push({ chrom, start, end, transcript_start: transcript_pos, })
    transcript_pos = transcript_pos + (end - start)
  }
  { name, chrom, strand, exons, }
}

///|
/// Get intron coordinates from transcript.
pub fn pysam_get_introns(transcript : PysamTranscript) -> Array[(Int32, Int32)] {
  let introns : Array[(Int32, Int32)] = Array(capacity=0)
  if transcript.exons.length() < 2 {
    return introns
  }
  for i = 0; i < transcript.exons.length() - 1; i = i + 1 {
    introns.push((transcript.exons[i].end, transcript.exons[i + 1].start))
  }
  introns
}

///|
/// Check if a genome position is within a transcript.
pub fn pysam_in_transcript(
  transcript : PysamTranscript,
  genome_pos : Int32,
) -> Bool {
  for i = 0; i < transcript.exons.length(); i = i + 1 {
    if genome_pos >= transcript.exons[i].start &&
      genome_pos < transcript.exons[i].end {
      return true
    }
  }
  false
}

///|
/// Convert CIGAR from transcript to genome coordinates.
/// Returns a new CIGAR with N operations for introns.
pub fn pysam_cigar_transcript_to_genome(
  transcript : PysamTranscript,
  transcript_pos : Int32,
  cigar : String,
) -> String {
  // For each position in the transcript alignment,
  // map to genome and insert N for introns
  let sb = StringBuilder()
  let mut tpos = transcript_pos
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      // For M operations, check if we cross exon boundaries
      if c == 77 { // M
        let mut remaining = num
        while remaining > 0 {
          // Find current exon
          let mut current_exon = -1
          for e = 0; e < transcript.exons.length(); e = e + 1 {
            let exon = transcript.exons[e]
            let exon_len = exon.end - exon.start
            if tpos >= exon.transcript_start &&
              tpos < exon.transcript_start + exon_len {
              current_exon = e
              break
            }
          }
          if current_exon < 0 {
            // In intron, skip
            tpos = tpos + remaining
            remaining = 0
          } else {
            let exon = transcript.exons[current_exon]
            let exon_len = exon.end - exon.start
            let pos_in_exon = tpos - exon.transcript_start
            let remaining_in_exon = exon_len - pos_in_exon
            let consume = if remaining < remaining_in_exon {
              remaining
            } else {
              remaining_in_exon
            }
            let _ = sb.write(consume.to_string() + "M")
            tpos = tpos + consume
            remaining = remaining - consume
            // If more remaining, we need to cross an intron
            if remaining > 0 && current_exon + 1 < transcript.exons.length() {
              let intron_len = transcript.exons[current_exon + 1].start -
                exon.end
              let _ = sb.write(intron_len.to_string() + "N")
            }
          }
        }
      } else {
        // Pass through other operations
        let _ = sb.write(num.to_string())
        let _ = sb.write_char(c.unsafe_to_char())
        if c == 77 || c == 73 || c == 83 || c == 61 || c == 88 {
          tpos = tpos + num
        }
      }
      num = 0
    }
  }
  sb.to_string()
}
