///|
fn show_optional_int(value : Int?) -> String {
  match value {
    Some(number) => number.to_string()
    None => "unknown"
  }
}

///|
fn show_optional_double(value : Double?) -> String {
  match value {
    Some(number) => number.to_string()
    None => "unknown"
  }
}

///|
fn main {
  println("=== Biopython Alignment.counts Demo ===")

  println("\n1. Classify and score an internal deletion")
  let gapped = @src.coordinate_pairwise_alignment(
    "target",
    "GAACT",
    "query",
    "GAT",
    [0, 2, 4, 5],
    [0, 2, 2, 3],
  ) catch {
    _ => abort("failed to build gapped alignment")
  }
  let affine = @src.AlignmentGapScores::affine(-5.0, -1.0) catch {
    _ => abort("failed to configure affine gap scores")
  }
  let scored = @src.AlignmentCountsConfig::scored(
    match_score=2.0,
    mismatch_score=-1.0,
    gap_scores=affine,
  ) catch {
    _ => abort("failed to configure alignment scoring")
  }
  let counts = gapped.alignment_counts(config=scored) catch {
    _ => abort("failed to calculate detailed counts")
  }
  println("  " + counts.summary())
  println(
    "  internal deletions=" +
    counts.internal_deletions().to_string() +
    " (open=" +
    counts.open_internal_deletions.to_string() +
    ", extend=" +
    counts.extend_internal_deletions.to_string() +
    ")",
  )
  println(
    "  substitution=" +
    show_optional_double(counts.substitution_score) +
    ", gap=" +
    show_optional_double(counts.gap_score) +
    ", total=" +
    show_optional_double(counts.score),
  )

  println("\n2. Count positives with the official BLOSUM45 matrix")
  let protein = @src.coordinate_pairwise_alignment(
    "protein_1",
    "EPQSDPSVEPPLSQETFSDLWKLLPE",
    "protein_2",
    "EPSSETGMDPPLSQETFEDLWSLLPD",
    [0, 26],
    [0, 26],
  ) catch {
    _ => abort("failed to build protein alignment")
  }
  let protein_counts = protein.alignment_counts(
    config=@src.AlignmentCountsConfig::with_matrix(@src.blosum45()),
  ) catch {
    _ => abort("failed to score protein alignment")
  }
  println(
    "  identities=" +
    protein_counts.identities.to_string() +
    ", mismatches=" +
    protein_counts.mismatches.to_string() +
    ", positives=" +
    show_optional_int(protein_counts.positives) +
    ", substitution score=" +
    show_optional_double(protein_counts.substitution_score),
  )

  println("\n3. Read residues through a reverse-strand coordinate path")
  let reverse = @src.coordinate_pairwise_alignment(
    "reference",
    "GAACT",
    "reverse_read",
    "ATC",
    [0, 2, 4, 5],
    [3, 1, 1, 0],
  ) catch {
    _ => abort("failed to build reverse-strand alignment")
  }
  let reverse_counts = reverse.alignment_counts() catch {
    _ => abort("failed to count reverse-strand alignment")
  }
  println("  " + reverse_counts.summary())
  println(
    "  reverse internal deletions=" +
    reverse_counts.internal_deletions().to_string(),
  )

  println("\n4. Sum counts and scores over every unordered MSA row pair")
  let msa = @src.coordinate_multiple_alignment(["a", "b", "c"], [
    "A-AA", "AAAA", "A--A",
  ]) catch {
    _ => abort("failed to build multiple alignment")
  }
  let msa_counts = msa.alignment_counts(config=scored) catch {
    _ => abort("failed to count multiple alignment")
  }
  println("  " + msa_counts.summary())
  println(
    "  gap opens=" +
    msa_counts.open_gaps().to_string() +
    ", gap extensions=" +
    msa_counts.extend_gaps().to_string() +
    ", total score=" +
    show_optional_double(msa_counts.score),
  )

  println("\n=== Demo Complete ===")
}
