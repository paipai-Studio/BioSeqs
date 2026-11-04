/// cobra_double_deletion.mbt — Double gene/reaction deletion (COBRApy flux_analysis.deletion)
/// Computes growth for all pairs of gene/reaction knockouts.

///|
/// Double deletion result.
pub struct CobraDoubleDeletionResult {
  gene_a : String
  gene_b : String
  growth : Double
  status : String
}

///|
/// Double gene deletion: knockout all pairs of genes and measure growth.
pub fn cobra_double_gene_deletion(
  model : FbaModel,
  gene_list~ : Array[String] = [],
  threshold~ : Double = 1.0e-6,
) -> Array[CobraDoubleDeletionResult] raise FbaError {
  // Collect genes
  let genes : Array[String] = Array(capacity=0)
  if gene_list.length() > 0 {
    for g in gene_list {
      genes.push(g)
    }
  } else {
    let gene_set : Map[String, Bool] = Map([], capacity=0)
    for r in model.reactions {
      for g in r.genes {
        gene_set.set(g, true)
      }
    }
    for g in gene_set.keys() {
      genes.push(g)
    }
    genes.sort_by(fn(a, b) { compare(a, b) })
  }
  let results : Array[CobraDoubleDeletionResult] = Array(capacity=0)
  for i = 0; i < genes.length(); i = i + 1 {
    for j = i; j < genes.length(); j = j + 1 {
      let ko_genes : Array[String] = [genes[i], genes[j]]
      let (ko_model, _) = fba_knockout(model, ko_genes)
      let sol = fba_solve(ko_model)
      let status = if sol.objective_value < threshold {
        "lethal"
      } else if sol.objective_value < 0.5 {
        "reduced"
      } else {
        "viable"
      }
      results.push({
        gene_a: genes[i],
        gene_b: genes[j],
        growth: sol.objective_value,
        status,
      })
    }
  }
  results
}

///|
/// Double reaction deletion: knockout all pairs of reactions.
pub fn cobra_double_reaction_deletion(
  model : FbaModel,
  reaction_list~ : Array[String] = [],
  threshold~ : Double = 1.0e-6,
) -> Array[CobraDoubleDeletionResult] raise FbaError {
  let rxns : Array[String] = Array(capacity=0)
  if reaction_list.length() > 0 {
    for r in reaction_list {
      rxns.push(r)
    }
  } else {
    for r in model.reactions {
      if r.objective == 0.0 {
        rxns.push(r.id)
      }
    }
  }
  let results : Array[CobraDoubleDeletionResult] = Array(capacity=0)
  for i = 0; i < rxns.length(); i = i + 1 {
    for j = i; j < rxns.length(); j = j + 1 {
      let ko1 = fba_set_reaction_bound(model, rxns[i], 0.0, 0.0)
      let ko2 = fba_set_reaction_bound(ko1, rxns[j], 0.0, 0.0)
      let sol = fba_solve(ko2)
      let status = if sol.objective_value < threshold {
        "lethal"
      } else if sol.objective_value < 0.5 {
        "reduced"
      } else {
        "viable"
      }
      results.push({
        gene_a: rxns[i],
        gene_b: rxns[j],
        growth: sol.objective_value,
        status,
      })
    }
  }
  results
}
