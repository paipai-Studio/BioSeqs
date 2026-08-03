///|
/// Survival Analysis module demo - Kaplan-Meier, log-rank, Cox PH model.
///
/// Run: moon run examples/survival_demo
fn main {
  println("=== Survival Analysis Module Demo ===")
  println("Kaplan-Meier estimation, log-rank test, Cox proportional hazards model")
  println("Reference: Therneau T (2023) A Package for Survival Analysis")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample survival data
  // ------------------------------------------------------------------
  println("=== 1. Load sample survival data ===")
  let data = @bio.survival_sample_data()
  println("n_samples: \{data.n_samples()}")
  println("n_groups: \{data.n_groups()}")
  let names = data.group_names()
  let mut i = 0
  while i < names.length() {
    println("  group[\{i}]: \{names[i]}")
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Kaplan-Meier curve for the control group
  // ------------------------------------------------------------------
  println("=== 2. Kaplan-Meier curve for 'control' group ===")
  let km_control = @bio.survival_kaplan_meier(data, "control")
  println(@bio.survival_km_summary(km_control))
  println(@bio.survival_km_to_string(km_control))
  println("")

  // ------------------------------------------------------------------
  // 3. Kaplan-Meier curve for the treatment group
  // ------------------------------------------------------------------
  println("=== 3. Kaplan-Meier curve for 'treatment' group ===")
  let km_treatment = @bio.survival_kaplan_meier(data, "treatment")
  println(@bio.survival_km_summary(km_treatment))
  println(@bio.survival_km_to_string(km_treatment))
  println("")

  // ------------------------------------------------------------------
  // 4. Median survival for both groups
  // ------------------------------------------------------------------
  println("=== 4. Median survival ===")
  match @bio.survival_median_survival(km_control) {
    Some(t) => println("Median survival (control): \{t}")
    None => println("Median survival (control): NA (curve did not cross 0.5)")
  }
  match @bio.survival_median_survival(km_treatment) {
    Some(t) => println("Median survival (treatment): \{t}")
    None => println("Median survival (treatment): NA (curve did not cross 0.5)")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Log-rank test comparing control vs treatment
  // ------------------------------------------------------------------
  println("=== 5. Log-rank test: control vs treatment ===")
  let lr = @bio.survival_log_rank(data, "control", "treatment")
  println(@bio.survival_logrank_summary(lr))
  println("")

  // ------------------------------------------------------------------
  // 6. Fit Cox proportional hazards model
  // ------------------------------------------------------------------
  println("=== 6. Cox proportional hazards model ===")
  let covariates = ["age", "treatment"]
  let cox = @bio.survival_cox_fit(data, covariates)
  println(@bio.survival_cox_summary(cox))
  println("")

  // ------------------------------------------------------------------
  // 7. Hazard ratios and p-values for each covariate
  // ------------------------------------------------------------------
  println("=== 7. Hazard ratios and p-values ===")
  println("Covariate\tHR\tPr(>|z|)")
  let p = cox.covariate_names.length()
  let mut k = 0
  while k < p {
    println(
      "\{cox.covariate_names[k]}\t\{cox.hazard_ratios[k]}\t\{cox.p_values[k]}",
    )
    k = k + 1
  }

  println("=== Demo complete ===")
}
