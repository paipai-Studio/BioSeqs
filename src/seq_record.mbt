// SeqRecord: an annotated biological sequence, mirroring BioPython's
// Bio.SeqRecord.SeqRecord.
//
// A SeqRecord holds:
//   - seq:          the biological sequence (a Seq object)
//   - id:           primary identifier (e.g. accession number)
//   - name:         common name / secondary id
//   - description:  human-readable description (FASTA title line content)
//   - letter_annotations: per-position annotations (e.g. "phred_quality"
//                    stores a list of quality scores, one per base)

///|
/// Error raised by SeqIO operations, mirroring BioPython's ValueError.
pub suberror SeqIOError {
  SeqIOError(String)
}

///|
/// An annotated biological sequence record, equivalent to
/// Bio.SeqRecord.SeqRecord.
pub struct SeqRecord {
  seq : Seq
  id : String
  name : String
  description : String
  letter_annotations : Map[String, Array[Int]]
} derive(Eq)

///|
/// Construct a SeqRecord with the given sequence and metadata.
pub fn SeqRecord::new(
  seq : Seq,
  id? : String = "<unknown id>",
  name? : String = "<unknown name>",
  description? : String = "<unknown description>",
) -> SeqRecord {
  SeqRecord::{
    seq,
    id,
    name,
    description,
    letter_annotations: Map([], capacity=0),
  }
}

///|
/// Construct a SeqRecord and return it.
pub fn SeqRecord::from_seq(seq : Seq) -> SeqRecord {
  SeqRecord::{
    seq,
    id: "<unknown id>",
    name: "<unknown name>",
    description: "<unknown description>",
    letter_annotations: Map([], capacity=0),
  }
}

///|
/// Set the phred quality scores (per-base integer quality values).
/// Returns a new SeqRecord (does not mutate the original).
pub fn SeqRecord::set_phred_quality(
  self : SeqRecord,
  quals : Array[Int],
) -> SeqRecord {
  let new_annotations : Map[String, Array[Int]] = Map(
    [],
    capacity=self.letter_annotations.length() + 1,
  )
  for k, v in self.letter_annotations {
    new_annotations[k] = v
  }
  new_annotations["phred_quality"] = quals
  SeqRecord::{ ..self, letter_annotations: new_annotations }
}

///|
/// Get the phred quality scores, if present.
pub fn SeqRecord::phred_quality(self : SeqRecord) -> Array[Int]? {
  self.letter_annotations.get("phred_quality")
}

///|
/// Return the sequence length.
pub fn SeqRecord::length(self : SeqRecord) -> Int {
  self.seq.length()
}

///|
/// Return a debug representation like BioPython's repr.
pub fn SeqRecord::repr(self : SeqRecord) -> String {
  "SeqRecord(seq=Seq('" +
  self.seq.to_string() +
  "'), id='" +
  self.id +
  "', name='" +
  self.name +
  "', description='" +
  self.description +
  "')"
}
