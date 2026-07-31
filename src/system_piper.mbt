///|
/// SystemPipeR module - Pipeline orchestration for bioinformatics analyses.
/// Inspired by Bioconductor SystemPipeR package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/SystemPipeR.html
///   - Backman TWH, Girke T. (2024) SystemPipeR: NGS workflow and report generation
///     environment. R package version 2.24.0.

///|
/// Pipeline step status.
pub enum StepStatus {
  Pending
  Running
  Completed
  Failed
  Skipped
}

///|
/// Create a Pending status.
pub fn step_status_pending() -> StepStatus {
  StepStatus::Pending
}

///|
/// Create a Running status.
pub fn step_status_running() -> StepStatus {
  StepStatus::Running
}

///|
/// Create a Completed status.
pub fn step_status_completed() -> StepStatus {
  StepStatus::Completed
}

///|
/// Create a Failed status.
pub fn step_status_failed() -> StepStatus {
  StepStatus::Failed
}

///|
/// Create a Skipped status.
pub fn step_status_skipped() -> StepStatus {
  StepStatus::Skipped
}

///|
/// Step status to string.
pub fn StepStatus::to_string(self : StepStatus) -> String {
  match self {
    StepStatus::Pending => "pending"
    StepStatus::Running => "running"
    StepStatus::Completed => "completed"
    StepStatus::Failed => "failed"
    StepStatus::Skipped => "skipped"
  }
}

///|
/// A pipeline step/command.
pub struct PipelineStep {
  step_id : String
  name : String
  description : String
  command : String
  args : Array[String]
  status : StepStatus
  exit_code : Int
  log : String
  duration_seconds : Double
  dependencies : Array[String]
  input_files : Array[String]
  output_files : Array[String]
  parameter_file : String
}

///|
/// A pipeline workflow.
pub struct Pipeline {
  pipeline_id : String
  name : String
  description : String
  steps : Array[PipelineStep]
  global_params : Map[String, String]
}

///|
/// A pipeline parameter configuration.
pub struct PipelineConfig {
  work_dir : String
  input_dir : String
  output_dir : String
  cores : Int
  overwrite : Bool
  cache : Bool
}

///|
/// Create a new pipeline.
pub fn Pipeline::new(pipeline_id : String, name : String) -> Pipeline {
  Pipeline::{
    pipeline_id,
    name,
    description: "",
    steps: Array::new(),
    global_params: Map([], capacity=10),
  }
}

///|
/// Set description.
pub fn Pipeline::set_description(self : Pipeline, desc : String) -> Pipeline {
  Pipeline::{
    pipeline_id: self.pipeline_id,
    name: self.name,
    description: desc,
    steps: self.steps,
    global_params: self.global_params,
  }
}

///|
/// Set a global parameter.
pub fn Pipeline::set_param(self : Pipeline, key : String, value : String) -> Pipeline {
  let params = self.global_params.copy()
  params[key] = value
  Pipeline::{
    pipeline_id: self.pipeline_id,
    name: self.name,
    description: self.description,
    steps: self.steps,
    global_params: params,
  }
}

///|
/// Get a global parameter.
pub fn Pipeline::get_param(self : Pipeline, key : String) -> String {
  match self.global_params.get(key) {
    Some(v) => v
    None => ""
  }
}

///|
/// Add a step to the pipeline.
pub fn Pipeline::add_step(self : Pipeline, step : PipelineStep) -> Pipeline {
  let steps = self.steps.copy()
  steps.push(step)
  Pipeline::{
    pipeline_id: self.pipeline_id,
    name: self.name,
    description: self.description,
    steps,
    global_params: self.global_params,
  }
}

///|
/// Get step by id.
pub fn Pipeline::get_step(self : Pipeline, step_id : String) -> PipelineStep {
  let mut i = 0
  while i < self.steps.length() {
    if self.steps[i].step_id == step_id {
      return self.steps[i]
    }
    i = i + 1
  }
  // Return empty step
  PipelineStep::{
    step_id: "",
    name: "not_found",
    description: "",
    command: "",
    args: Array::new(),
    status: StepStatus::Pending,
    exit_code: -1,
    log: "",
    duration_seconds: 0.0,
    dependencies: Array::new(),
    input_files: Array::new(),
    output_files: Array::new(),
    parameter_file: "",
  }
}

///|
/// Get pipeline id.
pub fn Pipeline::get_id(self : Pipeline) -> String {
  self.pipeline_id
}

///|
/// Get pipeline name.
pub fn Pipeline::get_name(self : Pipeline) -> String {
  self.name
}

///|
/// Get number of steps.
pub fn Pipeline::get_n_steps(self : Pipeline) -> Int {
  self.steps.length()
}

///|
/// Get steps.
pub fn Pipeline::get_steps(self : Pipeline) -> Array[PipelineStep] {
  self.steps
}

///|
/// Get completed count.
pub fn Pipeline::get_completed_count(self : Pipeline) -> Int {
  let mut count = 0
  let mut i = 0
  while i < self.steps.length() {
    if self.steps[i].status is StepStatus::Completed {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Get failed count.
pub fn Pipeline::get_failed_count(self : Pipeline) -> Int {
  let mut count = 0
  let mut i = 0
  while i < self.steps.length() {
    if self.steps[i].status is StepStatus::Failed {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Get pending count.
pub fn Pipeline::get_pending_count(self : Pipeline) -> Int {
  let mut count = 0
  let mut i = 0
  while i < self.steps.length() {
    if self.steps[i].status is StepStatus::Pending {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Get skipped count.
pub fn Pipeline::get_skipped_count(self : Pipeline) -> Int {
  let mut count = 0
  let mut i = 0
  while i < self.steps.length() {
    if self.steps[i].status is StepStatus::Skipped {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Get the overall pipeline progress as a percentage.
pub fn Pipeline::get_progress(self : Pipeline) -> Double {
  let n = self.steps.length()
  if n == 0 { return 0.0 }
  let completed = self.get_completed_count()
  let skipped = self.get_skipped_count()
  ((completed.to_double() + skipped.to_double()) / n.to_double()) * 100.0
}

///|
/// Create a new pipeline step.
pub fn PipelineStep::new(step_id : String, name : String, command : String) -> PipelineStep {
  PipelineStep::{
    step_id,
    name,
    description: "",
    command,
    args: Array::new(),
    status: StepStatus::Pending,
    exit_code: 0,
    log: "",
    duration_seconds: 0.0,
    dependencies: Array::new(),
    input_files: Array::new(),
    output_files: Array::new(),
    parameter_file: "",
  }
}

///|
/// Set step description.
pub fn PipelineStep::set_description(self : PipelineStep, desc : String) -> PipelineStep {
  PipelineStep::{
    step_id: self.step_id,
    name: self.name,
    description: desc,
    command: self.command,
    args: self.args,
    status: self.status,
    exit_code: self.exit_code,
    log: self.log,
    duration_seconds: self.duration_seconds,
    dependencies: self.dependencies,
    input_files: self.input_files,
    output_files: self.output_files,
    parameter_file: self.parameter_file,
  }
}

///|
/// Set step arguments.
pub fn PipelineStep::set_args(self : PipelineStep, args : Array[String]) -> PipelineStep {
  PipelineStep::{
    step_id: self.step_id,
    name: self.name,
    description: self.description,
    command: self.command,
    args,
    status: self.status,
    exit_code: self.exit_code,
    log: self.log,
    duration_seconds: self.duration_seconds,
    dependencies: self.dependencies,
    input_files: self.input_files,
    output_files: self.output_files,
    parameter_file: self.parameter_file,
  }
}

///|
/// Add a dependency.
pub fn PipelineStep::add_dependency(self : PipelineStep, dep : String) -> PipelineStep {
  let deps = self.dependencies.copy()
  deps.push(dep)
  PipelineStep::{
    step_id: self.step_id,
    name: self.name,
    description: self.description,
    command: self.command,
    args: self.args,
    status: self.status,
    exit_code: self.exit_code,
    log: self.log,
    duration_seconds: self.duration_seconds,
    dependencies: deps,
    input_files: self.input_files,
    output_files: self.output_files,
    parameter_file: self.parameter_file,
  }
}

///|
/// Add an input file.
pub fn PipelineStep::add_input(self : PipelineStep, file : String) -> PipelineStep {
  let inputs = self.input_files.copy()
  inputs.push(file)
  PipelineStep::{
    step_id: self.step_id,
    name: self.name,
    description: self.description,
    command: self.command,
    args: self.args,
    status: self.status,
    exit_code: self.exit_code,
    log: self.log,
    duration_seconds: self.duration_seconds,
    dependencies: self.dependencies,
    input_files: inputs,
    output_files: self.output_files,
    parameter_file: self.parameter_file,
  }
}

///|
/// Add an output file.
pub fn PipelineStep::add_output(self : PipelineStep, file : String) -> PipelineStep {
  let outputs = self.output_files.copy()
  outputs.push(file)
  PipelineStep::{
    step_id: self.step_id,
    name: self.name,
    description: self.description,
    command: self.command,
    args: self.args,
    status: self.status,
    exit_code: self.exit_code,
    log: self.log,
    duration_seconds: self.duration_seconds,
    dependencies: self.dependencies,
    input_files: self.input_files,
    output_files: outputs,
    parameter_file: self.parameter_file,
  }
}

///|
/// Get step id.
pub fn PipelineStep::get_id(self : PipelineStep) -> String {
  self.step_id
}

///|
/// Get step name.
pub fn PipelineStep::get_name(self : PipelineStep) -> String {
  self.name
}

///|
/// Get step command.
pub fn PipelineStep::get_command(self : PipelineStep) -> String {
  self.command
}

///|
/// Get step status.
pub fn PipelineStep::get_status(self : PipelineStep) -> StepStatus {
  self.status
}

///|
/// Get step duration.
pub fn PipelineStep::get_duration(self : PipelineStep) -> Double {
  self.duration_seconds
}

///|
/// Get dependencies.
pub fn PipelineStep::get_dependencies(self : PipelineStep) -> Array[String] {
  self.dependencies
}

///|
/// Check if all dependencies of a step are completed.
pub fn Pipeline::can_run_step(self : Pipeline, step : PipelineStep) -> Bool {
  let deps = step.dependencies
  if deps.length() == 0 { return true }

  let mut i = 0
  while i < deps.length() {
    let dep_step = self.get_step(deps[i])
    match dep_step.status {
      StepStatus::Completed => ()
      _ => { return false }
    }
    i = i + 1
  }
  true
}

///|
/// Generate a pipeline configuration.
pub fn PipelineConfig::new(work_dir : String, input_dir : String, output_dir : String) -> PipelineConfig {
  PipelineConfig::{
    work_dir,
    input_dir,
    output_dir,
    cores: 1,
    overwrite: false,
    cache: true,
  }
}

///|
/// Set number of cores.
pub fn PipelineConfig::set_cores(self : PipelineConfig, cores : Int) -> PipelineConfig {
  PipelineConfig::{
    work_dir: self.work_dir,
    input_dir: self.input_dir,
    output_dir: self.output_dir,
    cores,
    overwrite: self.overwrite,
    cache: self.cache,
  }
}

///|
/// Get work directory.
pub fn PipelineConfig::get_work_dir(self : PipelineConfig) -> String {
  self.work_dir
}

///|
/// Get input directory.
pub fn PipelineConfig::get_input_dir(self : PipelineConfig) -> String {
  self.input_dir
}

///|
/// Get output directory.
pub fn PipelineConfig::get_output_dir(self : PipelineConfig) -> String {
  self.output_dir
}

///|
/// Get number of cores.
pub fn PipelineConfig::get_cores(self : PipelineConfig) -> Int {
  self.cores
}

///|
/// Generate a text summary of the pipeline.
pub fn Pipeline::summary(self : Pipeline) -> String {
  let mut s = "=== Pipeline Summary ===\n".to_string()
  s = s + "ID: " + self.pipeline_id + "\n"
  s = s + "Name: " + self.name + "\n"
  if self.description.length() > 0 {
    s = s + "Description: " + self.description + "\n"
  }
  s = s + "Steps: " + self.steps.length().to_string() + "\n"
  s = s + "  Completed: " + self.get_completed_count().to_string() + "\n"
  s = s + "  Failed: " + self.get_failed_count().to_string() + "\n"
  s = s + "  Pending: " + self.get_pending_count().to_string() + "\n"
  s = s + "  Skipped: " + self.get_skipped_count().to_string() + "\n"
  s = s + "Progress: " + self.get_progress().to_string() + "%\n"
  s = s + "\nSteps:\n"

  let mut i = 0
  while i < self.steps.length() {
    let step = self.steps[i]
    s = s + "  [" + step.status.to_string() + "] " + step.step_id + ": " + step.name
    if step.duration_seconds > 0.0 {
      s = s + " (" + step.duration_seconds.to_string() + "s)"
    }
    s = s + "\n"
    i = i + 1
  }

  s
}

///|
/// Generate an ASCII visualization of the pipeline.
pub fn Pipeline::to_ascii(self : Pipeline) -> String {
  let mut result = ""
  result = result + "╔══════════════════════════════════════════╗\n"
  result = result + "║  Pipeline: " + self.name + "\n"
  result = result + "║  ID: " + self.pipeline_id + "\n"
  result = result + "║  Steps: " + self.steps.length().to_string() + "\n"
  result = result + "╚══════════════════════════════════════════╝\n\n"

  let progress = self.get_progress()
  let progress_bar_len = 20
  let filled = ((progress / 100.0) * progress_bar_len.to_double()).to_int()
  let mut bar = "["
  let mut i = 0
  while i < progress_bar_len {
    if i < filled { bar = bar + "█" }
    else { bar = bar + "░" }
    i = i + 1
  }
  bar = bar + "]"
  result = result + "Progress: " + bar + " " + progress.to_string() + "%\n\n"

  let mut j = 0
  while j < self.steps.length() {
    let step = self.steps[j]
    let status_icon = match step.status {
      StepStatus::Completed => "✓"
      StepStatus::Failed => "✗"
      StepStatus::Running => "●"
      StepStatus::Skipped => "⊘"
      StepStatus::Pending => "○"
    }
    result = result + "  " + status_icon + " [" + step.status.to_string() + "] "
    result = result + step.step_id + ": " + step.name + "\n"

    if step.dependencies.length() > 0 {
      result = result + "    depends: "
      let mut d = 0
      while d < step.dependencies.length() {
        result = result + step.dependencies[d]
        if d < step.dependencies.length() - 1 {
          result = result + ", "
        }
        d = d + 1
      }
      result = result + "\n"
    }

    if step.input_files.length() > 0 {
      result = result + "    inputs: "
      let mut i2 = 0
      while i2 < step.input_files.length() {
        result = result + step.input_files[i2]
        if i2 < step.input_files.length() - 1 {
          result = result + ", "
        }
        i2 = i2 + 1
      }
      result = result + "\n"
    }

    if step.output_files.length() > 0 {
      result = result + "    outputs: "
      let mut o = 0
      while o < step.output_files.length() {
        result = result + step.output_files[o]
        if o < step.output_files.length() - 1 {
          result = result + ", "
        }
        o = o + 1
      }
      result = result + "\n"
    }

    result = result + "\n"
    j = j + 1
  }

  result
}

///|
/// Create a sample RNA-seq analysis pipeline.
pub fn pipeline_sample() -> Pipeline {
  let mut pipe = Pipeline::new("rnaseq_pipeline_v1", "RNA-Seq Differential Expression Pipeline")
  pipe = pipe.set_description("Standard RNA-Seq analysis pipeline: QC -> Alignment -> Quantification -> DE Analysis -> Functional Enrichment")

  // Step 1: QC
  let step1 = PipelineStep::new("qc", "Quality Control", "fastqc")
  let step1b = step1.set_args(["-o", "qc_output", "--extract", "raw_data/*.fastq.gz"])
  let step1c = step1b.add_input("raw_data/sample1_R1.fastq.gz")
  let step1d = step1c.add_input("raw_data/sample1_R2.fastq.gz")
  let step1e = step1d.add_output("qc_output/sample1_fastqc.html")
  let step1_final = step1e.set_description("Quality control of raw reads using FastQC")
  pipe = pipe.add_step(step1_final)

  // Step 2: Trimming
  let step2 = PipelineStep::new("trim", "Adapter Trimming", "trimmomatic")
  let step2b = step2.set_args(["PE", "-phred33", "raw_data/sample1_R1.fastq.gz", "raw_data/sample1_R2.fastq.gz", "trimmed/sample1_R1_trimmed.fq.gz", "trimmed/sample1_R2_trimmed.fq.gz", "ILLUMINACLIP:adapters.fa:2:30:10"])
  let step2c = step2b.add_dependency("qc")
  let step2d = step2c.add_output("trimmed/sample1_R1_trimmed.fq.gz")
  let step2e = step2d.add_output("trimmed/sample1_R2_trimmed.fq.gz")
  let step2_final = step2e.set_description("Adapter and quality trimming with Trimmomatic")
  pipe = pipe.add_step(step2_final)

  // Step 3: Alignment
  let step3 = PipelineStep::new("align", "Read Alignment", "hisat2")
  let step3b = step3.set_args(["-x", "genome_index/genome", "-1", "trimmed/sample1_R1_trimmed.fq.gz", "-2", "trimmed/sample1_R2_trimmed.fq.gz", "-S", "alignment/sample1.sam"])
  let step3c = step3b.add_dependency("trim")
  let step3d = step3c.add_output("alignment/sample1.sam")
  let step3e = step3d.add_output("alignment/sample1.bam")
  let step3_final = step3e.set_description("Read alignment to reference genome using HISAT2")
  pipe = pipe.add_step(step3_final)

  // Step 4: Quantification
  let step4 = PipelineStep::new("quantify", "Gene Quantification", "featureCounts")
  let step4b = step4.set_args(["-a", "annotation.gtf", "-o", "counts/gene_counts.txt", "alignment/sample1.bam"])
  let step4c = step4b.add_dependency("align")
  let step4d = step4c.add_output("counts/gene_counts.txt")
  let step4_final = step4d.set_description("Gene-level quantification using featureCounts")
  pipe = pipe.add_step(step4_final)

  // Step 5: Differential Expression
  let step5 = PipelineStep::new("de_analysis", "Differential Expression", "DESeq2")
  let step5b = step5.set_args(["-i", "counts/gene_counts.txt", "-c", "metadata.csv", "-o", "de_results/"])
  let step5c = step5b.add_dependency("quantify")
  let step5d = step5c.add_output("de_results/deseq2_results.csv")
  let step5e = step5d.add_output("de_results/volcano_plot.png")
  let step5_final = step5e.set_description("Differential expression analysis using DESeq2")
  pipe = pipe.add_step(step5_final)

  // Step 6: Enrichment
  let step6 = PipelineStep::new("enrich", "Functional Enrichment", "clusterProfiler")
  let step6b = step6.set_args(["-i", "de_results/deseq2_results.csv", "-o", "enrichment/"])
  let step6c = step6b.add_dependency("de_analysis")
  let step6d = step6c.add_output("enrichment/go_terms.csv")
  let step6e = step6d.add_output("enrichment/pathways.csv")
  let step6_final = step6e.set_description("GO and pathway enrichment analysis")
  pipe = pipe.add_step(step6_final)

  // Add some global params
  pipe = pipe.set_param("reference_genome", "GRCh38.p14")
  pipe = pipe.set_param("species", "Homo sapiens")
  pipe = pipe.set_param("read_length", "150")
  pipe = pipe.set_param("sequencer", "Illumina NovaSeq")

  pipe
}