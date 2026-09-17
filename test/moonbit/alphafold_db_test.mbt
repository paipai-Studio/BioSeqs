///|
/// Tests for the Bio.PDB.alphafold_db port: URL construction, JSON prediction
/// parsing, file-name extraction, and the injected-fetch path.

///|
test "alphafold_predictions_url" {
  assert_eq(
    @src.alphafold_predictions_url("P00520"),
    "https://alphafold.com/api/prediction/P00520",
  )
}

///|
test "alphafold_cif_file_name" {
  assert_eq(
    @src.alphafold_cif_file_name(
      "https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.cif",
    ),
    "AF-P00520-F1-model_v4.cif",
  )
}

///|
test "alphafold_parse_predictions" {
  let preds = @src.alphafold_parse_predictions(@src.alphafold_sample_json()) catch {
    _ => abort("sample JSON should parse")
  }
  assert_eq(preds.length(), 2)

  let p0 = preds[0]
  assert_eq(p0.entry_id, "P00520")
  assert_eq(p0.pdb_id, "AF-P00520-F1")
  assert_eq(p0.uniprot_accession, "P00520")
  assert_eq(p0.uniprot_start, 1)
  assert_eq(p0.uniprot_end, 284)
  assert_eq(p0.latest_version, 4)
  assert_eq(p0.gene, "ABL1")
  assert_true(p0.cif_url.length() > 0)
  assert_true(@src.alphafold_cif_file_name(p0.cif_url) == "AF-P00520-F1-model_v4.cif")

  let p1 = preds[1]
  assert_eq(p1.entry_id, "P69905")
  assert_eq(p1.pdb_id, "AF-P69905-F1")
  assert_eq(p1.gene, "HBA1")
}

///|
test "alphafold_get_predictions_with_injected_fetch" {
  let calls : Array[String] = []
  let preds = @src.alphafold_get_predictions(
    fetch=fn(url : String, _post : String?) -> String {
      calls.push(url)
      @src.alphafold_sample_json()
    },
    qualifier="P00520",
  )
  assert_eq(calls[0], "https://alphafold.com/api/prediction/P00520")
  assert_eq(preds.length(), 2)
  assert_eq(preds[0].entry_id, "P00520")
}

///|
test "alphafold_parse_not_array" {
  let r = try {
    ignore(@src.alphafold_parse_predictions("not json"))
    false
  } catch {
    AlphafoldError(msg) => msg == "Expected a JSON array"
  }
  assert_true(r)
}