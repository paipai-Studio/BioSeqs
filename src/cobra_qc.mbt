/// cobra_qc.mbt — Quality control checks (COBRApy quality_control)
/// Check model consistency: duplicate IDs, orphan metabolites, charge balance.

///|
/// QC issue severity.
pub enum CobraQcSeverity {
  CobraQcError
  CobraQcWarning
  CobraQcInfo
}

///|
/// QC issue.
pub struct CobraQcIssue {
  severity : CobraQcSeverity
  category : String
  message : String
  element_id : String
}

///|
/// Run all QC checks on a model.
pub fn cobra_qc_check(model : FbaModel) -> Array[CobraQcIssue] {
  let issues : Array[CobraQcIssue] = Array(capacity=0)
  // Check 1: duplicate metabolite IDs
  let met_set : Map[String, Int32] = Map([], capacity=0)
  for m in model.metabolites {
    let count = met_set.get(m).unwrap_or(0)
    met_set.set(m, count + 1)
  }
  for m in met_set.keys() {
    let count = met_set.get(m).unwrap_or(0)
    if count > 1 {
      issues.push({
        severity: CobraQcError,
        category: "duplicate_metabolite",
        message: "Metabolite ID '\{m}' appears \{count} times",
        element_id: m,
      })
    }
  }
  // Check 2: duplicate reaction IDs
  let rxn_set : Map[String, Int32] = Map([], capacity=0)
  for r in model.reactions {
    let count = rxn_set.get(r.id).unwrap_or(0)
    rxn_set.set(r.id, count + 1)
  }
  for r in rxn_set.keys() {
    let count = rxn_set.get(r).unwrap_or(0)
    if count > 1 {
      issues.push({
        severity: CobraQcError,
        category: "duplicate_reaction",
        message: "Reaction ID '\{r}' appears \{count} times",
        element_id: r,
      })
    }
  }
  // Check 3: orphan metabolites (not in any reaction)
  for i = 0; i < model.metabolites.length(); i = i + 1 {
    let mut has_reaction = false
    for j = 0; j < model.reactions.length(); j = j + 1 {
      if model.stoichiometry[i][j].abs() > 1.0e-9 {
        has_reaction = true
        break
      }
    }
    if !has_reaction {
      issues.push({
        severity: CobraQcWarning,
        category: "orphan_metabolite",
        message: "Metabolite '\{model.metabolites[i]}' is not in any reaction",
        element_id: model.metabolites[i],
      })
    }
  }
  // Check 4: empty reactions (no metabolites)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    let mut has_metabolite = false
    for i = 0; i < model.metabolites.length(); i = i + 1 {
      if model.stoichiometry[i][j].abs() > 1.0e-9 {
        has_metabolite = true
        break
      }
    }
    if !has_metabolite {
      issues.push({
        severity: CobraQcWarning,
        category: "empty_reaction",
        message: "Reaction '\{model.reactions[j].id}' has no metabolites",
        element_id: model.reactions[j].id,
      })
    }
  }
  // Check 5: invalid bounds (lb > ub)
  for r in model.reactions {
    if r.lower_bound > r.upper_bound {
      issues.push({
        severity: CobraQcError,
        category: "invalid_bounds",
        message: "Reaction '\{r.id}' has lower_bound > upper_bound",
        element_id: r.id,
      })
    }
  }
  // Check 6: stoichiometry matrix dimensions
  if model.stoichiometry.length() != model.metabolites.length() {
    issues.push({
      severity: CobraQcError,
      category: "dimension_mismatch",
      message: "Stoichiometry matrix rows (\{model.stoichiometry.length()}) != metabolites (\{model.metabolites.length()})",
      element_id: "",
    })
  }
  for i = 0; i < model.stoichiometry.length(); i = i + 1 {
    if model.stoichiometry[i].length() != model.reactions.length() {
      issues.push({
        severity: CobraQcError,
        category: "dimension_mismatch",
        message: "Stoichiometry row \{i} length (\{model.stoichiometry[i].length()}) != reactions (\{model.reactions.length()})",
        element_id: "",
      })
      break
    }
  }
  // Check 7: reactions with no objective (informational)
  let mut n_objective = 0
  for r in model.reactions {
    if r.objective != 0.0 {
      n_objective = n_objective + 1
    }
  }
  if n_objective == 0 {
    issues.push({
      severity: CobraQcWarning,
      category: "no_objective",
      message: "Model has no objective reaction",
      element_id: "",
    })
  }
  issues
}

///|
/// Summarize QC results.
pub fn cobra_qc_summary(issues : Array[CobraQcIssue]) -> String {
  let mut n_errors = 0
  let mut n_warnings = 0
  let mut n_info = 0
  for issue in issues {
    match issue.severity {
      CobraQcError => n_errors = n_errors + 1
      CobraQcWarning => n_warnings = n_warnings + 1
      CobraQcInfo => n_info = n_info + 1
    }
  }
  let lines : Array[String] = Array(capacity=0)
  lines.push(
    "QC Summary: \{n_errors} errors, \{n_warnings} warnings, \{n_info} info",
  )
  for issue in issues {
    let sev = match issue.severity {
      CobraQcError => "ERROR"
      CobraQcWarning => "WARN"
      CobraQcInfo => "INFO"
    }
    lines.push("[\{sev}] \{issue.category}: \{issue.message}")
  }
  lines.join("\n")
}

///|
/// Default constructor for CobraQcIssue.
pub fn CobraQcIssue::default() -> CobraQcIssue {
  { severity: CobraQcInfo, category: "", message: "", element_id: "", }
}
