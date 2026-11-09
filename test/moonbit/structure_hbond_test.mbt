///|
/// Tests for structure_hbond.mbt.
///
/// NOTE: structure_hbond.mbt (24 lines) declares only the read-only HBond
/// record. It has no public constructor and no public functions, so an
/// external test package can neither instantiate HBond nor call any logic
/// from that module. The actual hydrogen-bond detection lives in
/// dssp_full.mbt (dssp_find_hbonds, dssp_count_hbonds) and the multi-model
/// frequency analysis lives in hbond_frequency.mbt (covered by
/// hbond_frequency_test.mbt). The test below documents the limitation by
/// verifying that the HBond type is at least referenceable from here.

///|
test "structure_hbond_type_referenceable" {
  let empty : Array[@src.HBond] = []
  assert_eq(empty.length(), 0)
}
