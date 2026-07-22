///|
fn main {
  println("=== Bio.Variation Demo ===")
  println("Inspired by Biopython Bio.Variation module")
  println("")

  println("1. Analyze SNP")
  let snp = @bio.analyze_snp(100, "A", "T", 50.0, "PASS")
  println("  Position: \{snp.position}")
  println("  Ref allele: \{snp.ref_allele}")
  println("  Alt allele: \{snp.alt_allele}")
  println("  Quality: \{snp.quality}")
  println("")

  println("2. Detect mutations from DNA sequence")
  let dna = "ATGTCTGGCTTA"
  let mutations = @bio.detect_mutations("BRCA1", "BRCA1_HUMAN", dna, 4, "A")
  println("  DNA sequence: \{dna}")
  println("  Mutations found: \{mutations.length()}")
  let mut i = 0
  while i < mutations.length() {
    let m = mutations[i]
    println("    Gene: \{m.gene}, Protein: \{m.protein}")
    println("    Position: \{m.position}")
    println("    Original AA: \{m.original_aa} -> Mutated AA: \{m.mutated_aa}")
    println("    Type: \{m.mutation_type}, Impact: \{m.impact}")
    i = i + 1
  }
  println("")

  println("3. Analyze amino acid substitution")
  let sub = @bio.analyze_amino_acid_substitution("G", "A")
  println("  Original: \{sub.original}")
  println("  Mutated: \{sub.mutated}")
  println("  BLOSUM62 score: \{sub.blosum62_score}")
  println("  Grantham distance: \{sub.grantham_distance}")
  println("")

  println("4. Predict damaging substitutions")
  let sub_damaging = @bio.analyze_amino_acid_substitution("G", "W")
  let is_damaging = @bio.predict_damaging(sub_damaging)
  println("  G -> W: Damaging = \{is_damaging}")
  
  let sub_safe = @bio.analyze_amino_acid_substitution("V", "I")
  let is_safe = @bio.predict_damaging(sub_safe)
  println("  V -> I: Damaging = \{is_safe}")
  println("")

  println("5. Parse VCF line")
  let vcf_line = "17\t7577538\t.\tG\tA\t.\tPASS\t."
  let vcf_result = @bio.parse_vcf_line(vcf_line)
  if vcf_result.is_some() {
    let record = vcf_result.unwrap()
    println("  Chromosome: \{record.chromosome}")
    println("  Position: \{record.position}")
    println("  Ref: \{record.ref_allele}")
    println("  Alt: \{record.alt_allele}")
  }
  println("")

  println("6. Create example variation data")
  let records = @bio.create_example_variation_data()
  println("  Created \{records.length()} variation records")
  i = 0
  while i < records.length() {
    let r = records[i]
    println("    \{r.chromosome}:\{r.position} \{r.ref_allele}->\{r.alt_allele}")
    let mut j = 0
    while j < r.mutations.length() {
      let m = r.mutations[j]
      println("      -> \{m.gene}:\{m.position} \{m.original_aa}->\{m.mutated_aa} (\{m.impact})")
      j = j + 1
    }
    i = i + 1
  }
}