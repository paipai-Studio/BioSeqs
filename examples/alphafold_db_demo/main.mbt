///|
/// Bio.PDB.alphafold_db Demo.
///
/// Demonstrates the AlphaFold Protein Structure Database client: query URL
/// construction, parsing the predictions JSON, and file-name extraction. The
/// HTTP transfer is injected, so this runs against a canned response.

///|
fn main raise {
  println("=== 1. Query URL ===")
  println(@src.alphafold_predictions_url("P00520"))

  println("\n=== 2. Parse predictions (canned JSON) ===")
  let preds = @src.alphafold_parse_predictions(@src.alphafold_sample_json())
  for i = 0; i < preds.length(); i = i + 1 {
    let p = preds[i]
    println(
      "\{p.pdb_id}  entry=\{p.entry_id}  \{p.gene}  \{p.uniprot_start}-\{p.uniprot_end}",
    )
    println("  cif = \{p.cif_url}")
  }

  println("\n=== 3. mmCIF file name from a prediction ===")
  println(@src.alphafold_cif_file_name(preds[0].cif_url))

  println("\n=== 4. get_predictions through an injected fetch ===")
  let calls : Array[String] = []
  let fetched = @src.alphafold_get_predictions(
    fetch=fn(url : String, _post : String?) -> String {
      calls.push(url)
      @src.alphafold_sample_json()
    },
    qualifier="P00520",
  )
  println("requested URL: \{calls[0]}")
  println("predictions: \{fetched.length()}")
}