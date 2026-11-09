// scvi_e2e_demo - end-to-end scvi-tools-style workflow in pure MoonBit.
//
// Synthetic counts -> setup_anndata -> train/validation split -> SCVI
// (autograd + Adam) -> latent representation -> SCANVI annotation ->
// Stereoscope deconvolution -> differential expression.

///|
fn scvi_demo_mean(x : Array[Double]) -> Double {
  if x.length() == 0 {
    return 0.0
  }
  let mut s = 0.0
  for i in 0..<x.length() {
    s += x[i]
  }
  s / x.length().to_double()
}

///|
fn main {
  println("=== scvi-tools style end-to-end demo ===")

  let bundle = @src.scvi_data_synthetic_iid(60, 40, 2, 7)
  let registry = @src.scvi_data_setup_anndata_scanvi(bundle)
  println("\n1. setup_anndata")
  println("   " + @src.scvi_data_registry_summary(registry))

  println("\n2. train / validation / test split")
  let split = @src.scvi_data_splitter_split(
    bundle.x.length(),
    @src.ScviDataSplitterConfig::new(0.8, 0.1, 0.1, true, 42),
  )
  println(
    "   train=" +
    split.train_indices.length().to_string() +
    " val=" +
    split.validation_indices.length().to_string() +
    " test=" +
    split.test_indices.length().to_string(),
  )

  println("\n3. SCVI training (autograd + Adam)")
  let model = @src.scvi_init_model(
    bundle.x[0].length(),
    2,
    @src.ScVIConfig::default(),
  )
  let trained = @src.scvi_train(model, bundle.x, 30)
  let latent = @src.scvi_get_latent_representation(trained, bundle.x)
  println(
    "   latent representation: " +
    latent.length().to_string() +
    " cells x " +
    latent[0].length().to_string() +
    " dimensions",
  )
  println("   mean latent[0] = " + scvi_demo_mean(latent[0]).to_string())

  println("\n4. SCANVI semi-supervised annotation")
  let labels = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  ]
  let x_labeled = Array(capacity=40)
  for i in 0..<40 {
    x_labeled.push(bundle.x[i])
  }
  let scanvi = @src.scanvi_init(
    bundle.x[0].length(),
    ["type_A", "type_B"],
    @src.SCANVIConfig::new(16, 5, 2, 0.005, 30),
  )
  let annotated = @src.scanvi_train(
    scanvi,
    x_labeled,
    labels,
    [
      [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      ],
    ],
    2,
  )
  let predictions = @src.scanvi_predict(annotated, x_labeled)
  println("   predicted labels for first 5 cells: " + predictions[0])

  println("\n5. Stereoscope spatial deconvolution")
  let reference = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
  let stereo = @src.scvi_stereoscope_fit(
    reference,
    [0, 1, 2],
    3,
    @src.ScviStereoscopeConfig::default(),
  )
  let proportions = @src.scvi_stereoscope_deconvolve(
    stereo,
    [10.0, 1.0, 0.0],
    @src.ScviStereoscopeConfig::default(),
  )
  println(
    "   spot proportions (type_0, type_1, type_2): " +
    proportions[0].to_string() +
    ", " +
    proportions[1].to_string() +
    ", " +
    proportions[2].to_string(),
  )

  println("\n6. Differential expression (one vs rest)")
  let de = @src.scvi_differential_expression(
    x_labeled, labels, 0, 1, "change", 0.0,
  )
  println("   DE result length = " + de.length().to_string())
  println("\nDone.")
}
