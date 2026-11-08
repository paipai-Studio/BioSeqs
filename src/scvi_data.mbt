// scvi_data.mbt - Data preparation helpers mirroring scvi.data
//
// Provides a lightweight AnnData bundle, field registration, synthetic data
// generation and multiome merging. The real AnnData container lives in
// anndata.mbt; this module focuses on the setup_anndata workflow used by
// scvi-tools models.

// ===========================================================================
// Bundle
// ===========================================================================

///|
/// Minimal in-memory single-cell data set used by the scVI model family.
pub struct ScviAnnDataBundle {
  x : Array[Array[Double]]
  obs_names : Array[String]
  var_names : Array[String]
  batch : Array[String]
  labels : Array[String]
  protein : Array[Array[Double]]
}

///|
pub fn ScviAnnDataBundle::new(
  x : Array[Array[Double]],
  obs_names : Array[String],
  var_names : Array[String],
  batch : Array[String],
  labels : Array[String],
  protein : Array[Array[Double]],
) -> ScviAnnDataBundle {
  { x, obs_names, var_names, batch, labels, protein, }
}

///|
pub fn ScviAnnDataBundle::n_obs(self : ScviAnnDataBundle) -> Int32 {
  self.x.length()
}

///|
pub fn ScviAnnDataBundle::n_vars(self : ScviAnnDataBundle) -> Int32 {
  if self.x.length() == 0 {
    0
  } else {
    self.x[0].length()
  }
}

///|
/// Ensure every bundle has a batch column (scvi-tools uses a single
/// synthetic batch when none is provided).
pub fn scvi_data_ensure_batch(bundle : ScviAnnDataBundle) -> ScviAnnDataBundle {
  let n = bundle.x.length()
  let batch = if bundle.batch.length() == n {
    bundle.batch
  } else {
    let out = Array(capacity=n)
    for _i in 0..<n {
      out.push("batch_0")
    }
    out
  }
  ScviAnnDataBundle::new(
    bundle.x,
    bundle.obs_names,
    bundle.var_names,
    batch,
    bundle.labels,
    bundle.protein,
  )
}

///|
fn scvi_data_unique(values : Array[String]) -> Array[String] {
  let out = Array(capacity=values.length())
  for i in 0..<values.length() {
    let mut seen = false
    for j in 0..<out.length() {
      if out[j] == values[i] {
        seen = true
      }
    }
    if not(seen) {
      out.push(values[i])
    }
  }
  out
}

// ===========================================================================
// Registry (setup_anndata)
// ===========================================================================

///|
/// Bookkeeping returned by `setup_anndata`: registered fields and categories.
pub struct ScviDataRegistry {
  n_obs : Int32
  n_vars : Int32
  batch_categories : Array[String]
  label_categories : Array[String]
  fields : Array[String]
  registered : Bool
}

///|
pub fn scvi_data_setup_anndata(bundle : ScviAnnDataBundle) -> ScviDataRegistry {
  let prepared = scvi_data_ensure_batch(bundle)
  let fields = Array(capacity=0)
  fields.push("X")
  fields.push("batch")
  if prepared.labels.length() > 0 {
    fields.push("labels")
  }
  {
    n_obs: prepared.x.length(),
    n_vars: prepared.n_vars(),
    batch_categories: scvi_data_unique(prepared.batch),
    label_categories: scvi_data_unique(prepared.labels),
    fields,
    registered: true,
  }
}

///|
pub fn scvi_data_setup_anndata_scvi(
  bundle : ScviAnnDataBundle,
) -> ScviDataRegistry {
  scvi_data_setup_anndata(bundle)
}

///|
pub fn scvi_data_setup_anndata_scanvi(
  bundle : ScviAnnDataBundle,
) -> ScviDataRegistry {
  scvi_data_setup_anndata(bundle)
}

///|
pub fn scvi_data_setup_anndata_totalvi(
  bundle : ScviAnnDataBundle,
) -> ScviDataRegistry {
  let registry = scvi_data_setup_anndata(bundle)
  scvi_data_register_field(registry, "protein_expression")
  registry
}

///|
pub fn scvi_data_register_field(
  registry : ScviDataRegistry,
  name : String,
) -> Unit {
  let mut seen = false
  for i in 0..<registry.fields.length() {
    if registry.fields[i] == name {
      seen = true
    }
  }
  if not(seen) {
    registry.fields.push(name)
  }
}

///|
pub fn scvi_data_registry_has_field(
  registry : ScviDataRegistry,
  name : String,
) -> Bool {
  for i in 0..<registry.fields.length() {
    if registry.fields[i] == name {
      return true
    }
  }
  false
}

///|
pub fn scvi_data_registry_summary(registry : ScviDataRegistry) -> String {
  let parts = Array(capacity=registry.fields.length())
  for i in 0..<registry.fields.length() {
    parts.push(registry.fields[i])
  }
  let mut joined = ""
  for i in 0..<parts.length() {
    if i > 0 {
      joined = joined + ","
    }
    joined = joined + parts[i]
  }
  "n_obs=" +
  registry.n_obs.to_string() +
  ";n_vars=" +
  registry.n_vars.to_string() +
  ";fields=" +
  joined
}

///|
/// Copy batch / label columns between two bundles (scvi.data.transfer_fields).
pub fn scvi_data_transfer_fields(
  target : ScviAnnDataBundle,
  source : ScviAnnDataBundle,
) -> ScviAnnDataBundle {
  let n = if target.x.length() < source.x.length() {
    target.x.length()
  } else {
    source.x.length()
  }
  let batch = Array(capacity=n)
  let labels = Array(capacity=n)
  for i in 0..<n {
    if i < source.batch.length() {
      batch.push(source.batch[i])
    } else {
      batch.push("batch_0")
    }
    if i < source.labels.length() {
      labels.push(source.labels[i])
    } else {
      labels.push("unknown")
    }
  }
  ScviAnnDataBundle::new(
    target.x,
    target.obs_names,
    target.var_names,
    batch,
    labels,
    target.protein,
  )
}

// ===========================================================================
// Synthetic data
// ===========================================================================

///|
/// Generate an `synthetic_iid`-style count matrix with batch-specific scaling
/// factors and repeating cell-type labels.
pub fn scvi_data_synthetic_iid(
  n_cells : Int32,
  n_genes : Int32,
  n_batches : Int32,
  seed : Int32,
) -> ScviAnnDataBundle {
  let x = Array(capacity=n_cells)
  let obs_names = Array(capacity=n_cells)
  let batch = Array(capacity=n_cells)
  let labels = Array(capacity=n_cells)
  let var_names = Array(capacity=n_genes)
  for g in 0..<n_genes {
    var_names.push("gene_" + g.to_string())
  }
  let n_b = if n_batches < 1 { 1 } else { n_batches }
  for i in 0..<n_cells {
    let b = i % n_b
    let batch_scale = 0.6 + 0.4 * (b + 1).to_double()
    let row = Array(capacity=n_genes)
    for g in 0..<n_genes {
      let gene_scale = 0.5 + 2.5 * scvi_dist_uniform(seed + g * 31)
      let rate = batch_scale * gene_scale
      row.push(scvi_dist_poisson_sample(rate, seed + i * 131 + g * 7))
    }
    x.push(row)
    obs_names.push("cell_" + i.to_string())
    if b == 0 {
      batch.push("batch_0")
    } else if b == 1 {
      batch.push("batch_1")
    } else if b == 2 {
      batch.push("batch_2")
    } else {
      batch.push("batch_" + b.to_string())
    }
    if i % 3 == 0 {
      labels.push("type_A")
    } else if i % 3 == 1 {
      labels.push("type_B")
    } else {
      labels.push("type_C")
    }
  }
  { x, obs_names, var_names, batch, labels, protein: Array(capacity=0), }
}

// ===========================================================================
// Multiome / multi-modal input
// ===========================================================================

///|
/// Merge RNA and ATAC count matrices column-wise (scvi.data.read_10x_multiome
/// keeps modality-specific feature names).
pub fn scvi_data_read_10x_multiome(
  rna_x : Array[Array[Double]],
  rna_names : Array[String],
  atac_x : Array[Array[Double]],
  atac_names : Array[String],
) -> ScviAnnDataBundle {
  let n = if rna_x.length() < atac_x.length() {
    rna_x.length()
  } else {
    atac_x.length()
  }
  let x = Array(capacity=n)
  let obs_names = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=rna_names.length() + atac_names.length())
    for g in 0..<rna_names.length() {
      if g < rna_x[i].length() {
        row.push(rna_x[i][g])
      } else {
        row.push(0.0)
      }
    }
    for g in 0..<atac_names.length() {
      if g < atac_x[i].length() {
        row.push(atac_x[i][g])
      } else {
        row.push(0.0)
      }
    }
    x.push(row)
    obs_names.push("cell_" + i.to_string())
  }
  let var_names = Array(capacity=rna_names.length() + atac_names.length())
  for g in 0..<rna_names.length() {
    var_names.push(rna_names[g])
  }
  for g in 0..<atac_names.length() {
    var_names.push(atac_names[g])
  }
  {
    x,
    obs_names,
    var_names,
    batch: Array(capacity=0),
    labels: Array(capacity=0),
    protein: Array(capacity=0),
  }
}

///|
/// Stack a protein (CITE-seq) matrix onto a bundle.
pub fn scvi_data_attach_protein(
  bundle : ScviAnnDataBundle,
  protein : Array[Array[Double]],
) -> ScviAnnDataBundle {
  ScviAnnDataBundle::new(
    bundle.x,
    bundle.obs_names,
    bundle.var_names,
    bundle.batch,
    bundle.labels,
    protein,
  )
}
