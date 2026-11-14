// Regression tests for GitHub Issue #109: parse_newick was flattening
// nested topology and dropping internal branch lengths.

///|
test "issue109_nested_topology_two_children" {
  let tree = @src.parse_newick("(A:0.1,(B:0.2,C:0.3):0.4);")
  assert_eq(tree.children.length(), 2)
  // Root has two children: leaf A and internal node (B,C).
  let a = tree.children[0]
  let internal_bc = tree.children[1]
  assert_eq(a.name, "A")
  assert_true((a.branch_length - 0.1).abs() < 0.000000001)
  assert_eq(internal_bc.children.length(), 2)
  assert_true((internal_bc.branch_length - 0.4).abs() < 0.000000001)
  // Verify total branch length: 0.1 + 0.4 + 0.2 + 0.3 = 1.0
  let total = @src.total_branch_length(tree)
  assert_true((total - 1.0).abs() < 0.000000001)
  // Max depth: 0.1 (A) vs 0.4 + 0.3 = 0.7 (C) — so 0.7
  let depth = @src.get_tree_depth(tree)
  assert_true((depth - 0.7).abs() < 0.000000001)
}

///|
test "issue109_roundtrip_write" {
  let tree = @src.parse_newick("(A:0.1,(B:0.2,C:0.3):0.4);")
  let written = @src.write_newick(tree, true)
  // Re-parse the serialized form and check invariants again.
  let reparsed = @src.parse_newick(written + ";")
  assert_eq(reparsed.children.length(), 2)
  let total = @src.total_branch_length(reparsed)
  assert_true((total - 1.0).abs() < 0.000000001)
}

///|
test "issue109_old_tests_still_pass_parseable" {
  let tree = @src.parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
  // Now root has 3 children: A, B, and internal (C,D)
  assert_eq(tree.children.length(), 3)
  assert_eq(@src.get_leaf_count(tree), 4)
  // Nodes: root + A + B + internal(C,D) + C + D = 6
  assert_eq(@src.get_node_count(tree), 6)
  // Leaves depth: A=0.1, B=0.2, C=0.5+0.3=0.8, D=0.5+0.4=0.9 -> max 0.9
  let depth = @src.get_tree_depth(tree)
  assert_true((depth - 0.9).abs() < 0.000000001)
  // Total: 0.1+0.2+0.5+0.3+0.4 = 1.5
  let total = @src.total_branch_length(tree)
  assert_true((total - 1.5).abs() < 0.000000001)
}

///|
test "issue109_simple_two_leaves" {
  let tree = @src.parse_newick("(A,B);")
  assert_eq(tree.children.length(), 2)
  assert_eq(@src.get_leaf_count(tree), 2)
}

///|
test "issue109_unparenthesized_single_root" {
  let tree = @src.parse_newick("F:0.5;")
  assert_eq(tree.children.length(), 0)
  assert_eq(tree.name, "F")
  assert_true((tree.branch_length - 0.5).abs() < 0.000000001)
}

///|
test "issue109_root_label_and_length" {
  let tree = @src.parse_newick("(A,B)C:0.3;")
  assert_eq(tree.children.length(), 2)
  assert_eq(tree.name, "C")
  assert_true((tree.branch_length - 0.3).abs() < 0.000000001)
}

///|
test "issue109_single_quoted_label" {
  let tree = @src.parse_newick("('A B':0.1,C);")
  assert_eq(tree.children[0].name, "A B")
  assert_true((tree.children[0].branch_length - 0.1).abs() < 0.000000001)
}

///|
test "issue109_double_quoted_label" {
  let tree = @src.parse_newick("(\"A,B\":0.1,C);")
  assert_eq(tree.children[0].name, "A,B")
}

///|
test "issue109_scientific_notation_branch_length" {
  // '+' in the exponent must not corrupt the number.
  let tree = @src.parse_newick("(A:1.5e+3,B);")
  assert_true((tree.children[0].branch_length - 1500.0).abs() < 0.000000001)
}

///|
test "issue109_whitespace_around_tokens" {
  let tree = @src.parse_newick("( A : 0.1 , B : 0.2 );")
  assert_eq(tree.children[0].name, "A")
  assert_true((tree.children[0].branch_length - 0.1).abs() < 0.000000001)
  assert_true((tree.children[1].branch_length - 0.2).abs() < 0.000000001)
}

///|
test "issue109_multi_tree_input_keeps_first_tree" {
  // Documented behavior: only the first tree in the string is parsed.
  let tree = @src.parse_newick("(A,B);(C,D);")
  assert_eq(tree.children.length(), 2)
  assert_eq(@src.get_leaf_count(tree), 2)
}

// ============================================================
// Probe group 1: quoted labels (Bug 1 regression surface).
// Bug 1: the shared char_at helper maps ' and " to their decimal
// code-point strings ("39"/"34"), which corrupted quoted labels.
// These probes pin down the exact quote-handling semantics.
// ============================================================

///|
test "issue109_probe_single_quote_no_decimal_leak" {
  // The original bug produced the label "39A" (39 = code point of ').
  let tree = @src.parse_newick("('A B':0.1,C);")
  let name = tree.children[0].name
  assert_eq(name, "A B")
  assert_true(!(name.contains("39") || name.contains("34")))
  assert_true((tree.children[0].branch_length - 0.1).abs() < 0.000000001)
}

///|
test "issue109_probe_double_quote_with_colon" {
  // A colon inside quotes must NOT trigger branch-length parsing;
  // the whole quoted text is the label.
  let tree = @src.parse_newick("(\"A:0.5\",C);")
  assert_eq(tree.children[0].name, "A:0.5")
  assert_eq(tree.children[0].branch_length, 0.0)
  assert_eq(tree.children.length(), 2)
}

///|
test "issue109_probe_double_quote_with_parens" {
  // Parentheses inside quotes must not be interpreted as subtree delimiters.
  let tree = @src.parse_newick("(\"(A)\",C);")
  assert_eq(tree.children[0].name, "(A)")
  assert_eq(tree.children[0].children.length(), 0)
  assert_eq(tree.children.length(), 2)
}

///|
test "issue109_probe_single_quote_wrapping_double" {
  let tree = @src.parse_newick("('\"q\"',C);")
  assert_eq(tree.children[0].name, "\"q\"")
}

///|
test "issue109_probe_double_quote_wrapping_single" {
  let tree = @src.parse_newick("(\"'p'\",C);")
  assert_eq(tree.children[0].name, "'p'")
}

///|
test "issue109_probe_quoted_internal_node_label" {
  // Quoted labels apply to internal nodes too, and their branch length
  // must still be parsed after the closing quote.
  let tree = @src.parse_newick("((A:0.1,B:0.2)'clade one':0.5);")
  let internal = tree.children[0]
  assert_eq(internal.name, "clade one")
  assert_eq(internal.children.length(), 2)
  assert_true((internal.branch_length - 0.5).abs() < 0.000000001)
}

///|
test "issue109_probe_empty_quoted_label" {
  // An empty quoted label yields "" (not code-point leakage).
  let tree = @src.parse_newick("('',C);")
  assert_eq(tree.children[0].name, "")
}

///|
test "issue109_probe_quoted_label_without_length" {
  let tree = @src.parse_newick("(\"A B\",C);")
  assert_eq(tree.children[0].name, "A B")
  assert_eq(tree.children[0].branch_length, 0.0)
}

///|
test "issue109_probe_unclosed_quote_tolerated" {
  // Documented lenient behavior: an unterminated quote consumes to the
  // end of input. The leading "(" makes this an internal node whose sole
  // child swallows the rest of the string; no crash, no hang.
  let tree = @src.parse_newick("('unclosed,C);")
  assert_eq(tree.children.length(), 1)
  assert_eq(tree.name, "")
  assert_eq(tree.children[0].name, "unclosed,C)")
}

///|
test "issue109_probe_non_ascii_quoted_label" {
  // Non-ASCII characters also bypassed char_at's whitelist in the old code;
  // quoted labels with CJK text must survive verbatim.
  let tree = @src.parse_newick("('人类':0.1);")
  assert_eq(tree.children[0].name, "人类")
  assert_true((tree.children[0].branch_length - 0.1).abs() < 0.000000001)
}

// ============================================================
// Probe group 2: branch-length scanning (Bug 2 regression surface).
// Bug 2: '+' was not in char_at's whitelist, so "1.5e+3" became
// "1.5e433" and parsed to a wildly wrong mass. These probes cover
// the full numeric grammar accepted after ':'.
// ============================================================

///|
test "issue109_probe_sci_negative_exponent" {
  let tree = @src.parse_newick("(A:1.5e-3,B);")
  assert_true((tree.children[0].branch_length - 0.0015).abs() < 0.000000001)
}

///|
test "issue109_probe_sci_uppercase_e" {
  let tree = @src.parse_newick("(A:2E2,B);")
  assert_true((tree.children[0].branch_length - 200.0).abs() < 0.000000001)
}

///|
test "issue109_probe_negative_branch_length" {
  // Negative lengths are valid Newick; parsing must preserve the sign
  // (write_newick may drop them, but parse must not lose them).
  let tree = @src.parse_newick("(A:-0.5,B);")
  assert_true((tree.children[0].branch_length - -0.5).abs() < 0.000000001)
}

///|
test "issue109_probe_leading_decimal_point" {
  let tree = @src.parse_newick("(A:.5,B);")
  assert_true((tree.children[0].branch_length - 0.5).abs() < 0.000000001)
}

///|
test "issue109_probe_integer_branch_length" {
  let tree = @src.parse_newick("(A:3,B);")
  assert_true((tree.children[0].branch_length - 3.0).abs() < 0.000000001)
}

///|
test "issue109_probe_invalid_length_falls_back_to_zero" {
  // Non-numeric text after ':' scans to an empty number → 0.0 fallback.
  // Lenient continuation: the scanner stops at 'x', so "xyz" is then
  // treated as an extra unquoted label child, and "B" still parses fine.
  let tree = @src.parse_newick("(A:xyz,B:0.2);")
  assert_eq(tree.children.length(), 3)
  assert_eq(tree.children[0].branch_length, 0.0)
  assert_eq(tree.children[1].name, "xyz")
  assert_true((tree.children[2].branch_length - 0.2).abs() < 0.000000001)
}

///|
test "issue109_probe_multiple_decimal_points_lenient_parse" {
  // NOTE: the shared parse_double is lenient — "1.2.3" accumulates both
  // fractional groups (1.2 + 0.3 = 1.5) rather than failing. This is
  // pre-existing parse_double behavior, not a parser bug; the tree
  // structure itself stays intact.
  let tree = @src.parse_newick("(A:1.2.3,B:0.2);")
  assert_eq(tree.children.length(), 2)
  assert_true((tree.children[0].branch_length - 1.5).abs() < 0.000000001)
  assert_true((tree.children[1].branch_length - 0.2).abs() < 0.000000001)
}

///|
test "issue109_probe_internal_branch_length_sum" {
  // Combined probe: leaf + internal branch lengths all contribute to the
  // total (the exact quantity Issue #109 reported as wrong).
  let tree = @src.parse_newick("((A:0.1,B:0.2):0.5,C:0.3);")
  let total = @src.total_branch_length(tree)
  assert_true((total - 1.1).abs() < 0.000000001)
  let internal = tree.children[0]
  assert_true((internal.branch_length - 0.5).abs() < 0.000000001)
}

// ============================================================
// Probe group: structural edge cases beyond label/length scanning.
// ============================================================

///|
test "issue109_probe_deep_nesting" {
  // Four nested single-child levels: root→n3→n2→n1→A.
  let tree = @src.parse_newick("((((A:0.1):0.2):0.3):0.4);")
  assert_eq(@src.get_node_count(tree), 5)
  assert_eq(@src.get_leaf_count(tree), 1)
  // Cumulative branch lengths: 0.4 + 0.3 + 0.2 + 0.1 = 1.0
  let total = @src.total_branch_length(tree)
  assert_true((total - 1.0).abs() < 0.000000001)
  let depth = @src.get_tree_depth(tree)
  assert_true((depth - 1.0).abs() < 0.000000001)
}

///|
test "issue109_probe_single_child_chain" {
  // A single-child internal node is unusual but legal Newick.
  let tree = @src.parse_newick("((A));")
  assert_eq(tree.children.length(), 1)
  assert_eq(tree.children[0].children.length(), 1)
  assert_eq(tree.children[0].children[0].name, "A")
  assert_eq(@src.get_node_count(tree), 3)
}

///|
test "issue109_probe_nhx_annotations_stripped" {
  // NHX comments must be removed before Newick parsing, leaving plain
  // labels and lengths intact.
  let tree = @src.parse_nhx("(A:0.1[&&NHX:S=human],B:0.2);")
  assert_eq(tree.children[0].name, "A")
  assert_true((tree.children[0].branch_length - 0.1).abs() < 0.000000001)
  assert_true((tree.children[1].branch_length - 0.2).abs() < 0.000000001)
}

///|
test "issue109_probe_multiline_whitespace_tolerated" {
  // Real-world Newick files wrap across lines; the parser skips all
  // whitespace between tokens.
  let tree = @src.parse_newick("(A:0.1,\n\tB:0.2)\n;")
  assert_eq(tree.children.length(), 2)
  assert_true((tree.children[0].branch_length - 0.1).abs() < 0.000000001)
  assert_true((tree.children[1].branch_length - 0.2).abs() < 0.000000001)
}

///|
test "issue109_probe_quoted_semicolon_inside_label" {
  // A ';' inside a quoted label must not terminate parsing (only the
  // trailing semicolon is stripped).
  let tree = @src.parse_newick("(\"A;B\",C);")
  assert_eq(tree.children[0].name, "A;B")
  assert_eq(tree.children.length(), 2)
}

///|
test "issue109_probe_find_and_prune_integration" {
  // End-to-end: parse a named internal clade, locate it, prune to its
  // children, and verify the pruned topology. prune_tree keeps the
  // original root shell (empty name, length 0) with surviving clades
  // as its children.
  let tree = @src.parse_newick("((A:0.1,B:0.2)I:0.5,C:0.3);")
  assert_eq(@src.get_leaf_count(tree), 3)
  let found = @src.find_node(tree, "I")
  assert_true(found is Some(_))
  let pruned = @src.prune_tree(tree, ["A", "B"])
  assert_true(pruned is Some(_))
  let p = match pruned {
    Some(t) => t
    None => tree
  }
  assert_eq(@src.get_leaf_count(p), 2)
  // Root shell is preserved; the clade I survives with its length.
  assert_eq(p.name, "")
  assert_eq(p.branch_length, 0.0)
  assert_eq(p.children.length(), 1)
  assert_eq(p.children[0].name, "I")
  assert_true((p.children[0].branch_length - 0.5).abs() < 0.000000001)
}
