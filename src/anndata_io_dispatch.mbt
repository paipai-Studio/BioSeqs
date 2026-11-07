/// anndata_io_dispatch.mbt — generic `adata.write()` / read dispatch.
///
/// Mirrors `AnnData.write` / `anndata.read_*`: the target format is chosen
/// from the filename extension (.h5ad / .zarr / .loom / .csv) and routed
/// to the matching writer/reader.

///|
/// Format for a filename: "h5ad" | "zarr" | "loom" | "csvs" (unknown
/// extensions fall back to "h5ad", like AnnData.write's default).
pub fn anndata_write_format(filename : String) -> String {
  if filename.has_suffix(".zarr") {
    "zarr"
  } else if filename.has_suffix(".loom") {
    "loom"
  } else if filename.has_suffix(".csv") || filename.has_suffix(".csvs") {
    "csvs"
  } else {
    "h5ad"
  }
}

///|
/// Serialize an AnnData according to the filename extension.
/// Returns the file content (h5ad JSON / zarr text / loom text / csv bundle).
pub fn anndata_write(adata : AnnData, filename : String) -> String {
  match anndata_write_format(filename) {
    "zarr" => anndata_write_zarr(adata)
    "loom" => anndata_write_loom(adata)
    "csvs" => {
      // anndata_write_csvs returns (obs_csv, var_csv, x_csv)
      let (obs_csv, var_csv, x_csv) = anndata_write_csvs(adata)
      x_csv + "\n===OBS===\n" + obs_csv + "\n===VAR===\n" + var_csv
    }
    _ => scanpy_write_h5ad(adata)
  }
}

///|
/// Deserialize an AnnData according to the filename extension.
/// Inverse of `anndata_write` for h5ad / zarr / loom / csvs. The csvs
/// bundle written by `anndata_write` ("===OBS===" / "===VAR===" markers)
/// is split back into its X / obs / var sections before parsing.
pub fn anndata_read(filename : String, content : String) -> AnnData? {
  match anndata_write_format(filename) {
    "zarr" => anndata_read_zarr(content)
    "loom" => anndata_read_loom(content)
    "csvs" => {
      let x_part = dispatch_bundle_section(content, "===OBS===")
      anndata_read_csv(x_part)
    }
    _ => scanpy_read_h5ad(content)
  }
}

///|
/// Extract the section of a csvs bundle that precedes `marker`.
fn dispatch_bundle_section(content : String, marker : String) -> String {
  let lines = content.split("\n").to_array()
  let sb = StringBuilder()
  for line in lines {
    let s = line.to_owned()
    if s == marker || s == "===VAR===" {
      break
    }
    let _ = sb.write(s)
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
/// Round-trip helper: write to the format implied by `filename`, then read
/// it back. Used to validate the dispatch table end to end.
pub fn anndata_write_read_roundtrip(
  adata : AnnData,
  filename : String,
) -> AnnData? {
  anndata_read(filename, anndata_write(adata, filename))
}
