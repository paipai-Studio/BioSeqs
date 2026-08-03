///|
/// velociraptor module demo - Single-cell RNA velocity analysis.
///
/// Demonstrates steady-state and kinetic (EM) models for estimating
/// transcription/splicing/degradation rates, velocity computation,
/// and embedding projection on simulated single-cell data.
///
/// Reference: Bergen V et al. (2020) scVelo. Nat Biotechnol 38:1408.
/// Run: moon run examples/velociraptor_demo
fn main {
  println("=== velociraptor Module Demo ===")
  println("Single-cell RNA velocity analysis")
  println("Reference: Bioconductor velociraptor / scVelo")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Data ===")
  let gene_data = @bio.velocity_sample_data()
  println("Genes: " + gene_data.length().to_string())
  println("Cells: " + gene_data[0].spliced().length().to_string())
  for gd in gene_data {
    let mean_s = gd.spliced()[0]
    let mean_u = gd.unspliced()[0]
    println(
      "  " +
        gd.gene_name() +
        ": spliced[0]=" +
        mean_s.to_string() +
        " unspliced[0]=" +
        mean_u.to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Steady-state model
  // ------------------------------------------------------------------
  println("=== 2. Steady-State Model ===")
  for gd in gene_data {
    let (ratio, ss_s, ss_u) = @bio.velocity_steady_state(
      gd.spliced(), gd.unspliced()
    )
    println(
      "  " +
        gd.gene_name() +
        ": gamma/beta=" +
        ratio.to_string() +
        " ss_s=" +
        ss_s.to_string() +
        " ss_u=" +
        ss_u.to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Kinetic model (EM)
  // ------------------------------------------------------------------
  println("=== 3. Kinetic Model (EM) ===")
  let kinetics : Array[@bio.GeneKinetics] = Array::new()
  for gd in gene_data {
    let k = @bio.velocity_kinetic_model(
      gd.gene_name(), gd.spliced(), gd.unspliced(), 50
    )
    kinetics.push(k)
    println(
      "  " +
        k.gene_name() +
        ": alpha=" +
        k.alpha().to_string() +
        " beta=" +
        k.beta().to_string() +
        " gamma=" +
        k.gamma().to_string() +
        " likelihood=" +
        k.likelihood().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Compute velocity vectors
  // ------------------------------------------------------------------
  println("=== 4. Velocity Vectors ===")
  let velocities = @bio.velocity_compute(gene_data, kinetics)
  println("Cells with velocity: " + velocities.length().to_string())
  let mut max_speed = 0.0
  let mut max_speed_cell = 0
  for i in 0..<velocities.length() {
    if velocities[i].speed() > max_speed {
      max_speed = velocities[i].speed()
      max_speed_cell = i
    }
  }
  println("Max speed cell: " + max_speed_cell.to_string() + " (speed=" + max_speed.to_string() + ")")
  // Show velocity for first 5 cells.
  for i in 0..<5 {
    let v = velocities[i]
    let mut vel_str = "["
    for j in 0..<v.velocity().length() {
      vel_str = vel_str + v.velocity()[j].to_string()
      if j < v.velocity().length() - 1 {
        vel_str = vel_str + ", "
      }
    }
    vel_str = vel_str + "]"
    println("  Cell " + i.to_string() + ": vel=" + vel_str + " speed=" + v.speed().to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Velocity embedding projection
  // ------------------------------------------------------------------
  println("=== 5. Velocity Embedding ===")
  let embeddings = @bio.velocity_sample_embeddings(velocities.length())
  let emb_vel = @bio.velocity_embedding(velocities, embeddings, 5)
  println("Projected velocities: " + emb_vel.length().to_string())
  for i in 0..<5 {
    let ev = emb_vel[i]
    println(
      "  Cell " +
        i.to_string() +
        ": pos=[" +
        ev.embedding()[0].to_string() +
        ", " +
        ev.embedding()[1].to_string() +
        "] vel=[" +
        ev.velocity_embedding()[0].to_string() +
        ", " +
        ev.velocity_embedding()[1].to_string() +
        "]",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Root cell identification
  // ------------------------------------------------------------------
  println("=== 6. Root Cell Identification ===")
  let roots = @bio.velocity_find_root_cells(velocities, 3)
  println("Top root cells (by speed):")
  for r in roots {
    println(
      "  Cell " +
        r.to_string() +
        " (speed=" +
        velocities[r].speed().to_string() +
        ")",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Transition matrix
  // ------------------------------------------------------------------
  println("=== 7. Transition Matrix ===")
  let matrix = @bio.velocity_transition_matrix(velocities, embeddings, 5, 1.0)
  println("Transition matrix: " + matrix.length().to_string() + " x " + matrix.length().to_string())
  // Show transition probabilities for cell 0.
  let mut row_str = "  Cell 0 -> "
  for j in 0..<5 {
    if matrix[0][j] > 0.0 {
      row_str = row_str + "cell" + j.to_string() + "(" + matrix[0][j].to_string() + ") "
    }
  }
  println(row_str)
  println("")

  println("=== Demo complete ===")
}
