///|
/// Module layer components: SCANVAE, TOTALVAE, AutoZIVAE, PEAKVAE,
/// MULTIVAE, Classifier, MRDeconv modules.
pub struct ScanvaeModule {
  n_genes : Int32
  n_latent : Int32
  n_labels : Int32
  vae_weights : Array[Array[Double]]
  classifier_weights : Array[Array[Double]]
  classifier_bias : Array[Double]
}

///|
pub fn scanvae_module_init(
  n_genes : Int32,
  n_latent : Int32,
  n_labels : Int32,
) -> ScanvaeModule {
  let make_w = (r : Int32, c : Int32) => {
    let w = Array(capacity=r)
    for _i in 0..<r {
      let row = Array(capacity=c)
      for _j in 0..<c {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  let bias = Array(capacity=n_labels)
  for _i in 0..<n_labels {
    bias.push(0.0)
  }
  {
    n_genes,
    n_latent,
    n_labels,
    vae_weights: make_w(n_latent, n_genes),
    classifier_weights: make_w(n_labels, n_latent),
    classifier_bias: bias,
  }
}

///|
pub fn scanvae_module_inference(
  mdl : ScanvaeModule,
  x : Array[Double],
) -> (Array[Double], Array[Double]) {
  let z = Array(capacity=mdl.n_latent)
  for i in 0..<mdl.n_latent {
    let mut sum = 0.0
    for j in 0..<x.length() {
      if j < mdl.vae_weights[i].length() {
        sum = sum + x[j] * mdl.vae_weights[i][j]
      }
    }
    z.push(sum)
  }
  let logits = Array(capacity=mdl.n_labels)
  for l in 0..<mdl.n_labels {
    let mut score = mdl.classifier_bias[l]
    for j in 0..<z.length() {
      if j < mdl.classifier_weights[l].length() {
        score = score + z[j] * mdl.classifier_weights[l][j]
      }
    }
    logits.push(score)
  }
  (z, logits)
}

///|
pub struct TotalvaeModule {
  n_genes : Int32
  n_proteins : Int32
  n_latent : Int32
  rna_encoder : Array[Array[Double]]
  protein_encoder : Array[Array[Double]]
  rna_decoder : Array[Array[Double]]
  protein_decoder : Array[Array[Double]]
  background_decoder : Array[Array[Double]]
}

///|
pub fn totalvae_module_init(
  n_genes : Int32,
  n_proteins : Int32,
  n_latent : Int32,
) -> TotalvaeModule {
  let make_w = (r : Int32, c : Int32) => {
    let w = Array(capacity=r)
    for _i in 0..<r {
      let row = Array(capacity=c)
      for _j in 0..<c {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  {
    n_genes,
    n_proteins,
    n_latent,
    rna_encoder: make_w(n_latent, n_genes),
    protein_encoder: make_w(n_latent, n_proteins),
    rna_decoder: make_w(n_genes, n_latent),
    protein_decoder: make_w(n_proteins, n_latent),
    background_decoder: make_w(n_proteins, n_latent),
  }
}

///|
pub fn totalvae_module_forward(
  mdl : TotalvaeModule,
  rna : Array[Double],
  protein : Array[Double],
) -> Array[Double] {
  let z = Array(capacity=mdl.n_latent)
  for i in 0..<mdl.n_latent {
    let mut sum = 0.0
    for j in 0..<rna.length() {
      if j < mdl.rna_encoder[i].length() {
        sum = sum + rna[j] * mdl.rna_encoder[i][j]
      }
    }
    for j in 0..<protein.length() {
      if j < mdl.protein_encoder[i].length() {
        sum = sum + protein[j] * mdl.protein_encoder[i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub struct AutozivaeModule {
  n_genes : Int32
  n_latent : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  dropout_logits : Array[Double]
}

///|
pub fn autozivae_module_init(
  n_genes : Int32,
  n_latent : Int32,
) -> AutozivaeModule {
  let make_w = (r : Int32, c : Int32) => {
    let w = Array(capacity=r)
    for _i in 0..<r {
      let row = Array(capacity=c)
      for _j in 0..<c {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  let logits = Array(capacity=n_genes)
  for _i in 0..<n_genes {
    logits.push(0.0)
  }
  {
    n_genes,
    n_latent,
    encoder: make_w(n_latent, n_genes),
    decoder: make_w(n_genes, n_latent),
    dropout_logits: logits,
  }
}

///|
pub fn autozivae_get_dropout_probs(mdl : AutozivaeModule) -> Array[Double] {
  let probs = Array(capacity=mdl.n_genes)
  for g in 0..<mdl.n_genes {
    probs.push(1.0 / (1.0 + @math.exp(-mdl.dropout_logits[g])))
  }
  probs
}

///|
pub struct PeakvaeModule {
  n_regions : Int32
  n_latent : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
}

///|
pub fn peakvae_module_init(
  n_regions : Int32,
  n_latent : Int32,
) -> PeakvaeModule {
  let make_w = (r : Int32, c : Int32) => {
    let w = Array(capacity=r)
    for _i in 0..<r {
      let row = Array(capacity=c)
      for _j in 0..<c {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  {
    n_regions,
    n_latent,
    encoder: make_w(n_latent, n_regions),
    decoder: make_w(n_regions, n_latent),
  }
}

///|
pub fn peakvae_module_loglik(
  mdl : PeakvaeModule,
  x : Array[Double],
  z : Array[Double],
) -> Double {
  let mut ll = 0.0
  for r in 0..<mdl.n_regions {
    let mut logit = 0.0
    for j in 0..<z.length() {
      if j < mdl.decoder[r].length() {
        logit = logit + z[j] * mdl.decoder[r][j]
      }
    }
    let p = 1.0 / (1.0 + @math.exp(-logit))
    if r < x.length() {
      if x[r] > 0.5 && p > 0.0 {
        ll = ll + @math.ln(p)
      } else if p < 1.0 {
        ll = ll + @math.ln(1.0 - p)
      }
    }
  }
  ll
}

///|
pub struct MultivaeModule {
  n_genes : Int32
  n_regions : Int32
  n_latent : Int32
  rna_encoder : Array[Array[Double]]
  atac_encoder : Array[Array[Double]]
  rna_decoder : Array[Array[Double]]
  atac_decoder : Array[Array[Double]]
}

///|
pub fn multivae_module_init(
  n_genes : Int32,
  n_regions : Int32,
  n_latent : Int32,
) -> MultivaeModule {
  let make_w = (r : Int32, c : Int32) => {
    let w = Array(capacity=r)
    for _i in 0..<r {
      let row = Array(capacity=c)
      for _j in 0..<c {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  {
    n_genes,
    n_regions,
    n_latent,
    rna_encoder: make_w(n_latent, n_genes),
    atac_encoder: make_w(n_latent, n_regions),
    rna_decoder: make_w(n_genes, n_latent),
    atac_decoder: make_w(n_regions, n_latent),
  }
}

///|
pub struct ClassifierModule {
  n_input : Int32
  n_labels : Int32
  weights : Array[Array[Double]]
  bias : Array[Double]
}

///|
pub fn classifier_module_init(
  n_input : Int32,
  n_labels : Int32,
) -> ClassifierModule {
  let w = Array(capacity=n_labels)
  for _i in 0..<n_labels {
    let row = Array(capacity=n_input)
    for _j in 0..<n_input {
      row.push(0.01)
    }
    w.push(row)
  }
  let b = Array(capacity=n_labels)
  for _i in 0..<n_labels {
    b.push(0.0)
  }
  { n_input, n_labels, weights: w, bias: b, }
}

///|
pub fn classifier_module_forward(
  cls : ClassifierModule,
  z : Array[Double],
) -> Array[Double] {
  let scores = Array(capacity=cls.n_labels)
  for l in 0..<cls.n_labels {
    let mut score = cls.bias[l]
    for j in 0..<z.length() {
      if j < cls.weights[l].length() {
        score = score + z[j] * cls.weights[l][j]
      }
    }
    scores.push(score)
  }
  scores
}

///|
pub struct MrDeconvModule {
  n_cell_types : Int32
  n_genes : Int32
  weights : Array[Array[Double]]
}

///|
pub fn mr_deconv_module_init(
  n_cell_types : Int32,
  n_genes : Int32,
) -> MrDeconvModule {
  let w = Array(capacity=n_cell_types)
  for _i in 0..<n_cell_types {
    let row = Array(capacity=n_genes)
    for _j in 0..<n_genes {
      row.push(0.01)
    }
    w.push(row)
  }
  { n_cell_types, n_genes, weights: w, }
}

///|
pub fn mr_deconv_module_deconvolve(
  mdl : MrDeconvModule,
  spot : Array[Double],
) -> Array[Double] {
  let props = Array(capacity=mdl.n_cell_types)
  for ct in 0..<mdl.n_cell_types {
    let mut score = 0.0
    for g in 0..<spot.length() {
      if g < mdl.weights[ct].length() {
        score = score + spot[g] * mdl.weights[ct][g]
      }
    }
    props.push(1.0 / (1.0 + @math.exp(-score)))
  }
  let mut total = 0.0
  for i in 0..<props.length() {
    total = total + props[i]
  }
  if total > 0.0 {
    for i in 0..<props.length() {
      props[i] = props[i] / total
    }
  }
  props
}

///|
pub fn _scvi_modules_dummy() -> Unit {
  let _ = ScanvaeModule::{
    n_genes: 0,
    n_latent: 0,
    n_labels: 0,
    vae_weights: [],
    classifier_weights: [],
    classifier_bias: [],
  }
  let _ = TotalvaeModule::{
    n_genes: 0,
    n_proteins: 0,
    n_latent: 0,
    rna_encoder: [],
    protein_encoder: [],
    rna_decoder: [],
    protein_decoder: [],
    background_decoder: [],
  }
  let _ = AutozivaeModule::{
    n_genes: 0,
    n_latent: 0,
    encoder: [],
    decoder: [],
    dropout_logits: [],
  }
  let _ = PeakvaeModule::{
    n_regions: 0,
    n_latent: 0,
    encoder: [],
    decoder: [],
  }
  let _ = MultivaeModule::{
    n_genes: 0,
    n_regions: 0,
    n_latent: 0,
    rna_encoder: [],
    atac_encoder: [],
    rna_decoder: [],
    atac_decoder: [],
  }
  let _ = ClassifierModule::{ n_input: 0, n_labels: 0, weights: [], bias: [], }
  let _ = MrDeconvModule::{ n_cell_types: 0, n_genes: 0, weights: [], }
}
