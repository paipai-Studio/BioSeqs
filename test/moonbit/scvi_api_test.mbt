///|
test "scvi_dataloader_batching" {
  let data = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let loader = @src.scvi_dataloader_init(data, 2, false, 1)
  assert_eq(@src.scvi_dataloader_n_obs(loader), 5)
  assert_eq(@src.scvi_dataloader_n_batches(loader), 3)
  assert_eq(@src.scvi_dataloader_batch(loader, 0).length(), 2)
  assert_eq(@src.scvi_dataloader_batch(loader, 2).length(), 1)
  let idx = @src.scvi_dataloader_batch_indices(loader, 1)
  assert_eq(idx.length(), 2)
  assert_eq(@src.scvi_dataloader_all_batches(loader).length(), 3)
  let shuffled = @src.scvi_dataloader_init(
    [[1.0], [2.0], [3.0], [4.0], [5.0], [6.0], [7.0], [8.0]],
    2,
    true,
    3,
  )
  let first = @src.scvi_dataloader_batch_indices(shuffled, 0)
  @src.scvi_dataloader_shuffle(shuffled)
  let second = @src.scvi_dataloader_batch_indices(shuffled, 0)
  assert_true(first[0] != second[0] || first[1] != second[1])
}

///|
test "scvi_data_splitter_proportions" {
  let config = @src.ScviDataSplitterConfig::new(0.6, 0.2, 0.2, false, 42)
  let split = @src.scvi_data_splitter_split(10, config)
  assert_eq(split.train_indices.length(), 6)
  assert_eq(split.validation_indices.length(), 2)
  assert_eq(split.test_indices.length(), 2)
}

///|
test "scvi_data_splitter_stratified" {
  let labels = [0, 0, 0, 1, 1, 1]
  let config = @src.ScviDataSplitterConfig::new(0.5, 0.25, 0.25, false, 42)
  let split = @src.scvi_data_splitter_stratified(labels, config)
  assert_eq(
    split.train_indices.length() +
    split.validation_indices.length() +
    split.test_indices.length(),
    6,
  )
  assert_true(split.train_indices.length() > 0)
}

///|
test "scvi_semisupervised_split" {
  let labels = [0, -1, 1, -1, 0]
  let split = @src.scvi_semisupervised_split(labels, -1)
  assert_eq(split.labeled_indices.length(), 3)
  assert_eq(split.unlabeled_indices.length(), 2)
  let batch = @src.scvi_semisupervised_batch(
    [[1.0], [2.0], [3.0], [4.0], [5.0]],
    split.labeled_indices,
    2,
    0,
  )
  assert_eq(batch.length(), 2)
}

///|
test "scvi_masked_loader" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let masked = @src.scvi_masked_loader_mask(data, 0.5, 7)
  assert_eq(masked.values.length(), 2)
  assert_eq(masked.mask[0].length(), 3)
  for i in 0..<masked.values.length() {
    for j in 0..<masked.values[i].length() {
      if masked.mask[i][j] == 0.0 {
        assert_true(masked.values[i][j] == 0.0)
      } else {
        assert_true(masked.mask[i][j] == 1.0)
      }
    }
  }
}

///|
test "scvi_concatenate_loader" {
  let loader = @src.scvi_concatenate_loader_init(
    [[[1.0], [2.0]], [[3.0], [4.0], [5.0]]],
    2,
  )
  assert_eq(@src.scvi_concatenate_loader_n_obs(loader), 5)
  assert_eq(@src.scvi_concatenate_loader_n_batches(loader), 3)
  assert_eq(@src.scvi_concatenate_loader_batch(loader, 0).length(), 2)
}

///|
test "scvi_setup_anndata_and_synthetic_data" {
  let bundle = @src.scvi_data_synthetic_iid(12, 8, 2, 3)
  let registry = @src.scvi_data_setup_anndata_scanvi(bundle)
  assert_eq(registry.n_obs, 12)
  assert_eq(registry.n_vars, 8)
  assert_eq(registry.batch_categories.length(), 2)
  assert_true(@src.scvi_data_registry_has_field(registry, "X"))
  assert_true(@src.scvi_data_registry_has_field(registry, "labels"))
  assert_true(@src.scvi_data_registry_summary(registry).contains("n_obs=12"))
  assert_eq(bundle.x.length(), 12)
  assert_eq(bundle.x[0].length(), 8)
}

///|
test "scvi_setup_anndata_totalvi_and_transfer" {
  let bundle = @src.scvi_data_synthetic_iid(6, 4, 1, 5)
  let with_protein = @src.scvi_data_attach_protein(bundle, [
    [1.0, 2.0],
    [3.0, 4.0],
    [1.0, 1.0],
    [2.0, 2.0],
    [0.0, 1.0],
    [1.0, 0.0],
  ])
  let registry = @src.scvi_data_setup_anndata_totalvi(with_protein)
  assert_true(@src.scvi_data_registry_has_field(registry, "protein_expression"))
  let target = @src.scvi_data_synthetic_iid(6, 4, 1, 9)
  let merged = @src.scvi_data_transfer_fields(target, bundle)
  assert_eq(merged.batch.length(), 6)
  assert_eq(merged.labels.length(), 6)
}

///|
test "scvi_data_multiome_merge" {
  let bundle = @src.scvi_data_read_10x_multiome(
    [[1.0, 2.0], [3.0, 4.0]],
    ["g1", "g2"],
    [[5.0], [6.0]],
    ["peak1"],
  )
  assert_eq(bundle.x.length(), 2)
  assert_eq(bundle.x[0].length(), 3)
  assert_eq(bundle.var_names.length(), 3)
}

///|
test "scvi_model_base_lifecycle" {
  let base = @src.scvi_model_base_init(10, 5, 2)
  assert_true(@src.scvi_model_base_validate(base))
  assert_true(!@src.scvi_model_base_is_trained(base))
  let config = @src.ScviDataSplitterConfig::new(0.7, 0.2, 0.1, false, 42)
  @src.scvi_model_base_split(base, config)
  assert_eq(@src.scvi_model_base_n_train(base), 7)
  assert_eq(@src.scvi_model_base_n_validation(base), 2)
  assert_eq(@src.scvi_model_base_n_test(base), 1)
  assert_true(@src.scvi_model_base_validate(base))
  @src.scvi_model_base_mark_trained(base)
  assert_true(@src.scvi_model_base_require_trained(base))
  assert_true(@src.scvi_model_base_summary(base).contains("trained=true"))
}

///|
test "scvi_training_plan_early_stopping" {
  let config = @src.ScviTrainingPlanConfig::new(
    "unsupervised", 50, 0.1, true, 3, 0.0,
  )
  let result = @src.scvi_training_plan_run(config, fn(epoch) {
    if epoch == 0 {
      1.0
    } else {
      0.5
    }
  })
  assert_true(result.stopped_early)
  assert_true(result.epochs_run < 50)
  assert_true(result.best_loss <= result.final_loss + 1.0e-12)
  assert_true(@src.scvi_plan_recent_losses(result, 2).length() <= 2)
  assert_true(@src.scvi_semisupervised_loss(1.0, 0.5, 0.2) == 1.1)
  assert_true(@src.scvi_classifier_loss([1.0, 0.0], 0) < 1.0e-6)
  assert_true(@src.scvi_entropy([1.0, 0.0]) < 1.0e-12)
}

///|
test "scvi_plan_adversarial_loss" {
  let loss = @src.scvi_adversarial_loss(0.9, 1)
  assert_true(loss > 0.0 && loss < 0.2)
}

///|
test "scvi_runtime_settings_and_utils" {
  let settings = @src.scvi_runtime_settings_default()
  assert_eq(settings.seed, 0)
  settings.set_seed(7)
  settings.set_verbosity(3)
  settings.set_progress_bar(false)
  assert_eq(settings.seed, 7)
  assert_true(settings.summary().contains("seed=7"))
  assert_eq(@src.scvi_utils_default_keys().length(), 7)
  let pairs = [("a", "1"), ("b", "2")]
  match @src.scvi_utils_get(pairs, "b") {
    Some(v) => assert_true(v == "2")
    None => fail("missing key")
  }
  assert_eq(@src.scvi_utils_subsample_indices(10, 1.0, 1).length(), 10)
  assert_eq(@src.scvi_utils_subsample_indices(10, 0.0, 1).length(), 0)
  assert_true(@src.scvi_utils_format_duration(3661) == "01:01:01")
  assert_true(@src.scvi_utils_accuracy([1, 0, 1], [1, 0, 0]) > 0.66)
}

///|
test "scvi_autotune_random_and_grid" {
  let ranges = [("lr", 0.001, 0.1)]
  let random = @src.scvi_autotune_random(
    ranges,
    @src.ScviTuneConfig::new(8, 1, "loss", false),
    fn(params) {
      match @src.scvi_autotune_param(params, "lr") {
        Some(v) => (v - 0.05).abs()
        None => 1.0
      }
    },
  )
  assert_eq(random.trials.length(), 8)
  assert_true(random.best_score <= 0.05)
  let grid = @src.scvi_autotune_grid([("lr", [0.01, 0.05, 0.1])], false, fn(
    params,
  ) {
    match @src.scvi_autotune_param(params, "lr") {
      Some(v) => (v - 0.05).abs()
      None => 1.0
    }
  })
  assert_eq(grid.trials.length(), 3)
  assert_true(grid.best_score < 1.0e-12)
}

///|
test "scvi_hub_metadata" {
  let meta = @src.ScviHubMetadata::new(
    "pbmc3k", "1.0.0", "demo model", "https://example.org/pbmc3k", "MIT", "1.1.0",
    "0.10.0", 2700, 2000,
  )
  assert_true(@src.scvi_hub_validate(meta))
  assert_true(@src.scvi_hub_matches(meta, "pbm"))
  assert_true(@src.scvi_hub_to_json(meta).contains("pbmc3k"))
  assert_true(
    @src.scvi_hub_cache_path("/cache", "pbmc3k", "1.0.0") ==
    "/cache/pbmc3k/1.0.0/model.pt",
  )
  assert_true(
    @src.scvi_hub_download_url("https://hub", "pbmc3k", "1.0.0") ==
    "https://hub/pbmc3k/resolve/1.0.0/model.pt",
  )
  assert_true(@src.scvi_hub_summary(meta).contains("2700x2000"))
  let rebuilt = @src.scvi_hub_from_pairs(@src.scvi_hub_model_card(meta))
  assert_true(rebuilt.model_name == "pbmc3k")
  assert_true(rebuilt.version == "1.0.0")
  assert_true(rebuilt.license == "MIT")
  assert_eq(rebuilt.n_obs, 2700)
  assert_eq(rebuilt.n_vars, 2000)
}
