// ============================================================================
// A9: OmmSystem::check_constraints 单元测试
// ============================================================================

test "omm_check_constraints_empty_system" {
  let sys = @src.OmmSystem::new()
  let pos = Array(capacity=0)
  let (max_v, mean_v, violations) = sys.check_constraints(pos)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(mean_v, 0.0, 1.0e-12))
  assert_eq(violations.length(), 0)
}

test "omm_check_constraints_single_constraint_ideal" {
  // 两原子相距 1.0 Å，约束目标 1.0 Å → 违反 = 0
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 1.0))
  let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let (max_v, mean_v, violations) = sys.check_constraints(pos)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(mean_v, 0.0, 1.0e-10))
  assert_eq(violations.length(), 1)
  let v = violations[0]
  assert_eq(v.particle1, 0)
  assert_eq(v.particle2, 1)
  assert_true(@moonbit.assert_num_close(v.target_distance, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v.actual_distance, 1.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v.absolute_violation, 0.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v.relative_violation, 0.0, 1.0e-10))
}

test "omm_check_constraints_single_constraint_violated" {
  // 两原子实际相距 1.5 Å，约束目标 1.0 Å → 违反 = 0.5
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 1.0))
  let pos = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]
  let (max_v, mean_v, violations) = sys.check_constraints(pos)
  assert_true(@moonbit.assert_num_close(max_v, 0.5, 1.0e-10))
  assert_true(@moonbit.assert_num_close(mean_v, 0.5, 1.0e-10))
  assert_eq(violations.length(), 1)
  let v = violations[0]
  assert_true(@moonbit.assert_num_close(v.absolute_violation, 0.5, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v.relative_violation, 0.5, 1.0e-10)) // 0.5/1.0
}

test "omm_check_constraints_3_atoms_water_tip3p" {
  // 两约束都要求 0.9572 Å，用正交坐标保证精确（不涉及浮点角度）
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new())) // O
  ignore(sys.add_particle(particle=@src.OmmParticle::new())) // H1
  ignore(sys.add_particle(particle=@src.OmmParticle::new())) // H2
  ignore(sys.add_constraint(0, 1, 0.9572))
  ignore(sys.add_constraint(0, 2, 0.9572))
  // O 在原点，H1 在 +x, H2 在 +y → x²+0²=0.9572², 0²+y²=0.9572²（精确）
  let pos = [
    [0.0, 0.0, 0.0],
    [0.9572, 0.0, 0.0],
    [0.0, 0.9572, 0.0],
  ]
  let (max_v, mean_v, violations) = sys.check_constraints(pos)
  assert_eq(violations.length(), 2)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(mean_v, 0.0, 1.0e-10))
}

test "omm_check_constraints_opposite_direction" {
  // 原子方向反了也不影响距离计算
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 2.0))
  let pos = [[5.0, 0.0, 0.0], [3.0, 0.0, 0.0]] // 距离 2.0
  let (max_v, _, _) = sys.check_constraints(pos)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-10))
}

test "omm_check_constraints_zero_target_distance" {
  // 约束目标 0 Å（虚拟原子重合），相对违反率应该是 0（避免除零）
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 0.0))
  let pos = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] // 实际 0
  let (max_v, _, violations) = sys.check_constraints(pos)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-10))
  assert_eq(violations.length(), 1)
  assert_true(@moonbit.assert_num_close(violations[0].relative_violation, 0.0, 1.0e-10))
}

test "omm_check_constraints_mean_violation_multiple" {
  // 3 个约束：0 违反、0.1 违反、0.2 违反 → mean = 0.1
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 1.0))  // ideal
  ignore(sys.add_constraint(0, 2, 1.0))  // +0.1
  ignore(sys.add_constraint(0, 3, 1.0))  // +0.2
  let pos = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],    // 约束 0-1: 实际 1.0 违反 0
    [1.1, 0.0, 0.0],    // 约束 0-2: 实际 1.1 违反 0.1
    [1.2, 0.0, 0.0],    // 约束 0-3: 实际 1.2 违反 0.2
  ]
  let (max_v, mean_v, violations) = sys.check_constraints(pos)
  assert_eq(violations.length(), 3)
  assert_true(@moonbit.assert_num_close(max_v, 0.2, 1.0e-10))
  assert_true(@moonbit.assert_num_close(mean_v, 0.1, 1.0e-10))
}

test "omm_check_constraints_diagonal_3d_sqrt_exact" {
  // 对角距离：原子在 (1,2,2) vs (4,6,2) → dx=3, dy=4, dz=0 → r=5.0
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 5.0))
  let pos = [[1.0, 2.0, 2.0], [4.0, 6.0, 2.0]]
  let (max_v, _, violations) = sys.check_constraints(pos)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(violations[0].actual_distance, 5.0, 1.0e-10))
}

test "omm_check_constraints_idx_field_matches_order" {
  // 三条约束，验证 constraint_idx 字段正确标记约束序号
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 1.0))
  ignore(sys.add_constraint(1, 2, 2.0))
  ignore(sys.add_constraint(0, 2, 3.0))
  let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
  let (_, _, violations) = sys.check_constraints(pos)
  assert_eq(violations.length(), 3)
  assert_eq(violations[0].constraint_idx, 0)
  assert_eq(violations[1].constraint_idx, 1)
  assert_eq(violations[2].constraint_idx, 2)
  assert_eq(violations[0].particle1, 0)
  assert_eq(violations[0].particle2, 1)
  assert_eq(violations[2].particle1, 0)
  assert_eq(violations[2].particle2, 2)
}

test "omm_check_constraints_no_constraints_all_zero" {
  // 有 N 个原子但 0 约束 → 返回全零且 violations 空
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  let pos = [[0.0, 0.0, 0.0], [10.0, 10.0, 10.0], [-5.0, -5.0, -5.0]]
  let (max_v, mean_v, violations) = sys.check_constraints(pos)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(mean_v, 0.0, 1.0e-12))
  assert_eq(violations.length(), 0)
}

test "omm_check_constraints_positions_oversized_no_crash" {
  // positions 数组比约束引用的原子多 — 不应越界，多余原子被忽略
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_particle(particle=@src.OmmParticle::new()))
  ignore(sys.add_constraint(0, 1, 1.0))
  // 提供 4 个原子位置，约束只引用 0 和 1
  let pos = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [999.0, 999.0, 999.0],
    [888.0, 888.0, 888.0],
  ]
  let (max_v, mean_v, violations) = sys.check_constraints(pos)
  assert_eq(violations.length(), 1)
  assert_true(@moonbit.assert_num_close(max_v, 0.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(mean_v, 0.0, 1.0e-10))
}
