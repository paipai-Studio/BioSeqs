///|
/// MultiAssayExperiment - Multi-omics data coordination inspired by Bioconductor MultiAssayExperiment.
/// Coordinates experiments across multiple modalities (RNA-seq, ChIP-seq, mutations, etc.).

///|
/// A single experiment within the multi-assay container.
pub struct Experiment {
  name : String
  assay_type : String
  data : Array[Array[Double]]
  row_names : Array[String]
  col_names : Array[String]
}

///|
/// Construct Experiment.
pub fn Experiment::new(
  name~ : String,
  assay_type~ : String,
  data~ : Array[Array[Double]],
  row_names~ : Array[String],
  col_names~ : Array[String],
) -> Experiment {
  Experiment::{ name, assay_type, data, row_names, col_names }
}

///|
/// Get number of rows.
pub fn Experiment::nrow(self : Experiment) -> Int {
  self.data.length()
}

///|
/// Get number of columns.
pub fn Experiment::ncol(self : Experiment) -> Int {
  if self.data.length() > 0 {
    self.data[0].length()
  } else {
    0
  }
}

///|
/// Sample mapping between biological samples and assay columns.
pub struct SampleMap {
  assays : Array[String]
  primary_ids : Array[String]
  colnames : Array[String]
}

///|
/// Construct SampleMap.
pub fn SampleMap::new(
  assays~ : Array[String],
  primary_ids~ : Array[String],
  colnames~ : Array[String],
) -> SampleMap {
  SampleMap::{ assays, primary_ids, colnames }
}

///|
/// Get unique primary sample IDs.
pub fn SampleMap::primary_ids(self : SampleMap) -> Array[String] {
  let result : Array[String] = Array::new()
  for id in self.primary_ids {
    let mut found = false
    for r in result {
      if r == id {
        found = true
      }
    }
    if not(found) {
      result.push(id)
    }
  }
  result
}

///|
/// MultiAssayExperiment container.
pub struct MultiAssayExperiment {
  experiments : Array[Experiment]
  sample_map : SampleMap
  col_data : Array[(String, Map[String, String])]
  metadata : Map[String, String]
}

///|
/// Construct MultiAssayExperiment.
pub fn MultiAssayExperiment::new(
  experiments~ : Array[Experiment],
  sample_map~ : SampleMap,
  col_data~ : Array[(String, Map[String, String])],
  metadata? : Map[String, String] = Map([]),
) -> MultiAssayExperiment {
  MultiAssayExperiment::{ experiments, sample_map, col_data, metadata }
}

///|
/// Get experiment by name.
pub fn MultiAssayExperiment::get_experiment(
  self : MultiAssayExperiment,
  name : String,
) -> Experiment? {
  for exp in self.experiments {
    if exp.name == name {
      return Some(exp)
    }
  }
  None
}

///|
/// Get all experiment names.
pub fn MultiAssayExperiment::experiment_names(
  self : MultiAssayExperiment,
) -> Array[String] {
  self.experiments.map(fn(e) { e.name })
}

///|
/// Get number of experiments.
pub fn MultiAssayExperiment::num_experiments(
  self : MultiAssayExperiment,
) -> Int {
  self.experiments.length()
}

///|
/// Get primary sample IDs.
pub fn MultiAssayExperiment::primary_ids(
  self : MultiAssayExperiment,
) -> Array[String] {
  self.sample_map.primary_ids()
}

///|
/// Get number of primary samples.
pub fn MultiAssayExperiment::num_samples(self : MultiAssayExperiment) -> Int {
  self.primary_ids().length()
}

///|
/// Find common samples across all experiments.
pub fn MultiAssayExperiment::common_samples(
  self : MultiAssayExperiment,
) -> Array[String] {
  if self.experiments.length() == 0 {
    return []
  }
  let first_exp = self.experiments[0]
  let common : Array[String] = Array::new()
  for col in first_exp.col_names {
    let mut in_all = true
    for exp in self.experiments {
      let mut found = false
      for c in exp.col_names {
        if c == col {
          found = true
        }
      }
      if not(found) {
        in_all = false
      }
    }
    if in_all {
      common.push(col)
    }
  }
  common
}

///|
/// Subset by experiment names.
pub fn MultiAssayExperiment::subset_by_experiment(
  self : MultiAssayExperiment,
  names : Array[String],
) -> MultiAssayExperiment {
  let filtered = self.experiments.filter(fn(e) {
    let mut found = false
    for n in names {
      if n == e.name {
        found = true
      }
    }
    found
  })
  MultiAssayExperiment::new(
    experiments=filtered,
    sample_map=self.sample_map,
    col_data=self.col_data,
    metadata=self.metadata,
  )
}

///|
/// Subset by primary sample IDs.
pub fn MultiAssayExperiment::subset_by_sample(
  self : MultiAssayExperiment,
  sample_ids : Array[String],
) -> MultiAssayExperiment {
  let new_experiments : Array[Experiment] = Array::new()
  for exp in self.experiments {
    let new_data : Array[Array[Double]] = Array::new()
    let new_cols : Array[String] = Array::new()
    let mut col_idx = 0
    while col_idx < exp.col_names.length() {
      let col_name = exp.col_names[col_idx]
      let mut should_keep = false
      let mut sm_idx = 0
      while sm_idx < self.sample_map.colnames.length() {
        if self.sample_map.colnames[sm_idx] == col_name {
          let pid = self.sample_map.primary_ids[sm_idx]
          for sid in sample_ids {
            if sid == pid {
              should_keep = true
            }
          }
        }
        sm_idx = sm_idx + 1
      }
      if should_keep {
        new_cols.push(col_name)
        let row : Array[Double] = Array::new()
        let mut r = 0
        while r < exp.data.length() {
          row.push(exp.data[r][col_idx])
          r = r + 1
        }
        new_data.push(row)
      }
      col_idx = col_idx + 1
    }
    // Transpose back
    let transposed : Array[Array[Double]] = Array::new()
    if new_data.length() > 0 {
      let n_rows = exp.data.length()
      let mut r = 0
      while r < n_rows {
        let row : Array[Double] = Array::new()
        let mut c = 0
        while c < new_data.length() {
          row.push(new_data[c][r])
          c = c + 1
        }
        transposed.push(row)
        r = r + 1
      }
    }
    new_experiments.push(
      Experiment::new(
        name=exp.name,
        assay_type=exp.assay_type,
        data=transposed,
        row_names=exp.row_names,
        col_names=new_cols,
      ),
    )
  }
  MultiAssayExperiment::new(
    experiments=new_experiments,
    sample_map=self.sample_map,
    col_data=self.col_data,
    metadata=self.metadata,
  )
}

///|
/// Get column metadata for a sample.
pub fn MultiAssayExperiment::sample_col_data(
  self : MultiAssayExperiment,
  sample_id : String,
) -> Map[String, String]? {
  let mut i = 0
  while i < self.col_data.length() {
    let entry = self.col_data[i]
    if entry.0 == sample_id {
      return Some(entry.1)
    }
    i = i + 1
  }
  None
}

///|
/// Get assay data for a specific experiment.
pub fn MultiAssayExperiment::assay(
  self : MultiAssayExperiment,
  experiment_name : String,
) -> Array[Array[Double]]? {
  match self.get_experiment(experiment_name) {
    Some(exp) => Some(exp.data)
    None => None
  }
}

///|
/// Get row names for a specific experiment.
pub fn MultiAssayExperiment::rownames(
  self : MultiAssayExperiment,
  experiment_name : String,
) -> Array[String]? {
  match self.get_experiment(experiment_name) {
    Some(exp) => Some(exp.row_names)
    None => None
  }
}

///|
/// Get col names for a specific experiment.
pub fn MultiAssayExperiment::colnames(
  self : MultiAssayExperiment,
  experiment_name : String,
) -> Array[String]? {
  match self.get_experiment(experiment_name) {
    Some(exp) => Some(exp.col_names)
    None => None
  }
}

///|
/// Get dimensions of all experiments.
pub fn MultiAssayExperiment::dims(
  self : MultiAssayExperiment,
) -> Array[(String, Int, Int)] {
  self.experiments.map(fn(e) { (e.name, e.nrow(), e.ncol()) })
}

///|
/// Add metadata.
pub fn MultiAssayExperiment::add_metadata(
  self : MultiAssayExperiment,
  key : String,
  value : String,
) -> MultiAssayExperiment {
  let new_metadata = self.metadata
  new_metadata[key] = value
  MultiAssayExperiment::new(
    experiments=self.experiments,
    sample_map=self.sample_map,
    col_data=self.col_data,
    metadata=new_metadata,
  )
}

///|
/// Get metadata.
pub fn MultiAssayExperiment::get_metadata(
  self : MultiAssayExperiment,
  key : String,
) -> String? {
  self.metadata.get(key)
}

///|
/// Get summary.
pub fn MultiAssayExperiment::summary(self : MultiAssayExperiment) -> String {
  let n_exp = self.num_experiments()
  let n_samples = self.num_samples()
  let common = self.common_samples().length()
  "MultiAssayExperiment: " +
  n_exp.to_string() +
  " experiments, " +
  n_samples.to_string() +
  " primary samples, " +
  common.to_string() +
  " common samples"
}

///|
/// Create a multi-omics experiment from multiple data sources.
pub fn create_multi_omics(
  rna_counts : Array[Array[Double]],
  rna_genes : Array[String],
  rna_samples : Array[String],
  chip_scores : Array[Array[Double]],
  chip_peaks : Array[String],
  chip_samples : Array[String],
  sample_ids : Array[String],
) -> MultiAssayExperiment {
  let rna_exp = Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=rna_counts,
    row_names=rna_genes,
    col_names=rna_samples,
  )
  let chip_exp = Experiment::new(
    name="ChIP-seq",
    assay_type="scores",
    data=chip_scores,
    row_names=chip_peaks,
    col_names=chip_samples,
  )
  let assays : Array[String] = Array::new()
  let primary : Array[String] = Array::new()
  let colnames : Array[String] = Array::new()
  for s in rna_samples {
    assays.push("RNA-seq")
    colnames.push(s)
  }
  for i in chip_samples {
    assays.push("ChIP-seq")
    colnames.push(i)
  }
  let mut idx = 0
  while idx < rna_samples.length() {
    if idx < sample_ids.length() {
      primary.push(sample_ids[idx])
    }
    idx = idx + 1
  }
  idx = 0
  while idx < chip_samples.length() {
    if idx < sample_ids.length() {
      primary.push(sample_ids[idx])
    }
    idx = idx + 1
  }
  let sm = SampleMap::new(assays~, primary_ids=primary, colnames~)
  let col_data : Array[(String, Map[String, String])] = Array::new()
  for sid in sample_ids {
    let m : Map[String, String] = Map([])
    col_data.push((sid, m))
  }
  MultiAssayExperiment::new(
    experiments=[rna_exp, chip_exp],
    sample_map=sm,
    col_data~,
  )
}
