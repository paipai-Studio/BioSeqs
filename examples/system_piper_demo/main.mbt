///|
/// SystemPipeR module example - Pipeline orchestration for bioinformatics.
///
/// Run: moon run examples/system_piper_demo
fn main {
  println("=== SystemPipeR Module Demo ===")
  println("")

  // Create a pipeline
  println("=== Creating Pipeline ===")
  let pipe = @bio.Pipeline::new("demo_pipeline", "Variant Calling Pipeline")
  println("ID: \{pipe.get_id()}")
  println("Name: \{pipe.get_name()}")
  println("Steps: \{pipe.get_n_steps()}")
  println("")

  // Set description and params
  let pipe2 = pipe.set_description("Variant calling pipeline from aligned BAM files")
    .set_param("reference", "GRCh38")
    .set_param("species", "Homo sapiens")
    .set_param("work_dir", "/data/workspace")
  println("Description updated.")
  println("Reference: \{pipe2.get_param("reference")}")
  println("Species: \{pipe2.get_param("species")}")
  println("")

  // Build pipeline steps
  println("=== Building Pipeline Steps ===")

  // Step 1: Index reference
  let step1 = @bio.PipelineStep::new("index", "Index Reference", "bwa")
    .set_args(["index", "-a", "bwtsw", "ref.fa"])
    .set_description("Build BWA index for reference genome")
    .add_output("ref.fa.bwa")
  let pipe3 = pipe2.add_step(step1)
  println("Step 1 added: Index Reference")

  // Step 2: Alignment
  let step2 = @bio.PipelineStep::new("align", "Read Alignment", "bwa")
    .set_args(["mem", "-t", "8", "-R", "@RG\tID:test\tSM:sample", "ref.fa", "reads_R1.fq.gz", "reads_R2.fq.gz"])
    .set_description("Align reads to reference using BWA-MEM")
    .add_dependency("index")
    .add_input("reads_R1.fq.gz")
    .add_input("reads_R2.fq.gz")
    .add_output("aligned.sam")
  let pipe4 = pipe3.add_step(step2)
  println("Step 2 added: Read Alignment (depends on index)")

  // Step 3: Sort BAM
  let step3 = @bio.PipelineStep::new("sort", "Sort & Index BAM", "samtools")
    .set_args(["sort", "-@", "4", "-o", "sorted.bam", "aligned.sam"])
    .set_description("Sort BAM file by coordinates")
    .add_dependency("align")
    .add_output("sorted.bam")
  let pipe5 = pipe4.add_step(step3)
  println("Step 3 added: Sort & Index BAM (depends on align)")

  // Step 4: Mark duplicates
  let step4 = @bio.PipelineStep::new("dedup", "Mark Duplicates", "picard")
    .set_args(["MarkDuplicates", "I=sorted.bam", "O=dedup.bam", "M=dup_metrics.txt"])
    .set_description("Mark PCR duplicates using Picard")
    .add_dependency("sort")
    .add_output("dedup.bam")
  let pipe6 = pipe5.add_step(step4)
  println("Step 4 added: Mark Duplicates (depends on sort)")

  // Step 5: Variant calling
  let step5 = @bio.PipelineStep::new("call", "Variant Calling", "bcftools")
    .set_args(["mpileup", "-f", "ref.fa", "dedup.bam", "|", "bcftools", "call", "-mv", "-o", "variants.vcf"])
    .set_description("Call variants using bcftools")
    .add_dependency("dedup")
    .add_output("variants.vcf")
  let pipe7 = pipe6.add_step(step5)
  println("Step 5 added: Variant Calling (depends on dedup)")

  // Step 6: Filtering
  let step6 = @bio.PipelineStep::new("filter", "Variant Filtering", "bcftools")
    .set_args(["filter", "-e", "INFO/DP<10 | QUAL<30", "variants.vcf", "-o", "filtered.vcf"])
    .set_description("Filter variants by depth and quality")
    .add_dependency("call")
    .add_output("filtered.vcf")
  let pipe8 = pipe7.add_step(step6)
  println("Step 6 added: Variant Filtering (depends on call)")
  println("")

  // Pipeline progress
  println("=== Pipeline Progress ===")
  println("Steps: \{pipe8.get_n_steps()}")
  println("Completed: \{pipe8.get_completed_count()}")
  println("Failed: \{pipe8.get_failed_count()}")
  println("Pending: \{pipe8.get_pending_count()}")
  println("Progress: @{pipe8.get_progress()}%")
  println("")

  // Check if steps can run
  println("=== Dependency Checking ===")
  let all_steps = pipe8.get_steps()
  let mut i = 0
  while i < all_steps.length() {
    let step = all_steps[i]
    let can_run = pipe8.can_run_step(step)
    println("  \{step.step_id}: can run = \{can_run}")
    i = i + 1
  }
  println("")

  // Pipeline configuration
  println("=== Pipeline Configuration ===")
  let config = @bio.PipelineConfig::new("/work/project", "/data/input", "/data/output")
    .set_cores(16)
  println("Work dir: \{config.get_work_dir()}")
  println("Input dir: \{config.get_input_dir()}")
  println("Output dir: \{config.get_output_dir()}")
  println("Cores: \{config.get_cores()}")
  println("")

  // Pipeline summary
  println("=== Pipeline Summary ===")
  println(pipe8.summary())

  // ASCII visualization
  println("=== Pipeline Visualization ===")
  println(pipe8.to_ascii())

  // Sample pipeline
  println("=== Sample RNA-Seq Pipeline ===")
  let sample = @bio.pipeline_sample()
  println(sample.summary())

  println("=== Done ===")
}