///|
/// Test file for SystemPipeR module.

test "pipeline_new" {
  let pipe = @src.Pipeline::new("test_v1", "Test Pipeline")
  assert_eq(pipe.get_id(), "test_v1")
  assert_eq(pipe.get_name(), "Test Pipeline")
  assert_eq(pipe.get_n_steps(), 0)
}

test "pipeline_set_description" {
  let pipe = @src.Pipeline::new("test", "Test")
  let pipe2 = pipe.set_description("A test pipeline")
  // Check via summary
  let s = pipe2.summary()
  assert_true(s.contains("A test pipeline"))
}

test "pipeline_set_param" {
  let pipe = @src.Pipeline::new("test", "Test")
  let pipe2 = pipe.set_param("key1", "value1")
  assert_eq(pipe2.get_param("key1"), "value1")
}

test "pipeline_get_param_missing" {
  let pipe = @src.Pipeline::new("test", "Test")
  assert_eq(pipe.get_param("nonexistent"), "")
}

test "pipeline_add_step" {
  let pipe = @src.Pipeline::new("test", "Test")
  let step = @src.PipelineStep::new("step1", "Step 1", "echo")
  let pipe2 = pipe.add_step(step)
  assert_eq(pipe2.get_n_steps(), 1)
}

test "pipeline_get_step" {
  let pipe = @src.Pipeline::new("test", "Test")
  let step = @src.PipelineStep::new("step1", "Step 1", "echo")
  let pipe2 = pipe.add_step(step)
  let s = pipe2.get_step("step1")
  assert_eq(s.get_id(), "step1")
  assert_eq(s.get_name(), "Step 1")
}

test "pipeline_get_step_not_found" {
  let pipe = @src.Pipeline::new("test", "Test")
  let s = pipe.get_step("missing")
  assert_eq(s.get_id(), "")
}

test "pipeline_get_steps" {
  let pipe = @src.Pipeline::new("test", "Test")
  let pipe2 = pipe.add_step(@src.PipelineStep::new("s1", "S1", "cmd1"))
  let pipe3 = pipe2.add_step(@src.PipelineStep::new("s2", "S2", "cmd2"))
  let steps = pipe3.get_steps()
  assert_eq(steps.length(), 2)
}

test "pipeline_completed_count" {
  let pipe = @src.Pipeline::new("test", "Test")
  assert_eq(pipe.get_completed_count(), 0)
}

test "pipeline_failed_count" {
  let pipe = @src.Pipeline::new("test", "Test")
  assert_eq(pipe.get_failed_count(), 0)
}

test "pipeline_pending_count" {
  let pipe = @src.Pipeline::new("test", "Test")
  assert_eq(pipe.get_pending_count(), 0)
}

test "pipeline_progress_empty" {
  let pipe = @src.Pipeline::new("test", "Test")
  assert_eq(pipe.get_progress(), 0.0)
}

test "step_new" {
  let step = @src.PipelineStep::new("s1", "My Step", "bwa")
  assert_eq(step.get_id(), "s1")
  assert_eq(step.get_name(), "My Step")
  assert_eq(step.get_command(), "bwa")
}

test "step_set_description" {
  let step = @src.PipelineStep::new("s1", "Step", "cmd")
  let step2 = step.set_description("A step for testing")
  assert_eq(step2.description, "A step for testing")
}

test "step_set_args" {
  let step = @src.PipelineStep::new("s1", "Step", "cmd")
  let step2 = step.set_args(["arg1", "arg2", "-v"])
  assert_eq(step2.args.length(), 3)
}

test "step_add_dependency" {
  let step = @src.PipelineStep::new("s2", "Step 2", "cmd")
  let step2 = step.add_dependency("s1")
  let deps = step2.get_dependencies()
  assert_eq(deps.length(), 1)
  assert_eq(deps[0], "s1")
}

test "step_add_input" {
  let step = @src.PipelineStep::new("s1", "Step", "cmd")
  let step2 = step.add_input("input.fq")
  assert_eq(step2.input_files.length(), 1)
}

test "step_add_output" {
  let step = @src.PipelineStep::new("s1", "Step", "cmd")
  let step2 = step.add_output("output.bam")
  assert_eq(step2.output_files.length(), 1)
}

test "step_status" {
  let step = @src.PipelineStep::new("s1", "Step", "cmd")
  assert_true(step.get_status() is @src.StepStatus::Pending)
}

test "step_duration" {
  let step = @src.PipelineStep::new("s1", "Step", "cmd")
  assert_eq(step.get_duration(), 0.0)
}

test "step_status_to_string" {
  // Test default status via PipelineStep (which has Pending by default)
  let step = @src.PipelineStep::new("s1", "Step", "cmd")
  let status_str = step.get_status().to_string()
  assert_eq(status_str, "pending")
}

test "pipeline_config_new" {
  let config = @src.PipelineConfig::new("/work", "/input", "/output")
  assert_eq(config.get_work_dir(), "/work")
  assert_eq(config.get_input_dir(), "/input")
  assert_eq(config.get_output_dir(), "/output")
  assert_eq(config.get_cores(), 1)
}

test "pipeline_config_set_cores" {
  let config = @src.PipelineConfig::new("/w", "/i", "/o")
  let config2 = config.set_cores(8)
  assert_eq(config2.get_cores(), 8)
}

test "pipeline_summary" {
  let pipe = @src.Pipeline::new("test", "My Pipeline")
  let s = pipe.summary()
  assert_true(s.contains("Pipeline Summary"))
  assert_true(s.contains("My Pipeline"))
}

test "pipeline_to_ascii" {
  let pipe = @src.Pipeline::new("test", "My Pipeline")
  let ascii = pipe.to_ascii()
  assert_true(ascii.contains("My Pipeline"))
  assert_true(ascii.contains("Progress"))
}

test "pipeline_sample" {
  let pipe = @src.pipeline_sample()
  assert_true(pipe.get_n_steps() >= 5)
  assert_eq(pipe.get_id(), "rnaseq_pipeline_v1")
  assert_true(pipe.get_progress() >= 0.0)
}

test "pipeline_sample_steps" {
  let pipe = @src.pipeline_sample()
  let steps = pipe.get_steps()
  let first = steps[0]
  assert_eq(first.get_id(), "qc")
  assert_eq(first.get_name(), "Quality Control")
}

test "pipeline_global_params" {
  let pipe = @src.pipeline_sample()
  assert_eq(pipe.get_param("reference_genome"), "GRCh38.p14")
  assert_eq(pipe.get_param("species"), "Homo sapiens")
}

test "step_can_run_no_deps" {
  let pipe = @src.Pipeline::new("test", "Test")
  let step = @src.PipelineStep::new("s1", "Step 1", "cmd")
  assert_true(pipe.can_run_step(step))
}