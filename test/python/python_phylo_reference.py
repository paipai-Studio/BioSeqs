#!/usr/bin/env python3
"""Python reference script for Bio.Phylo comparison tests."""

import json
from io import StringIO
from Bio import Phylo

def get_clade_info(clade):
    """Get basic info from a clade."""
    return {
        "name": clade.name,
        "is_terminal": clade.is_terminal(),
        "branch_length": clade.branch_length,
        "confidence": clade.confidence,
        "children": len(clade.clades)
    }

def get_tree_summary(tree):
    """Get summary info from a tree."""
    terminals = tree.get_terminals()
    nonterminals = tree.get_nonterminals()
    
    return {
        "rooted": tree.rooted,
        "weight": tree.weight,
        "name": tree.name,
        "num_terminals": len(terminals),
        "num_nonterminals": len(nonterminals),
        "terminal_names": sorted([t.name for t in terminals if t.name]),
        "is_bifurcating": tree.is_bifurcating()
    }

def test_basic_tree():
    """Test basic tree structure."""
    newick = "(((A,B),(C,D)),(E,F,G));"
    tree = Phylo.read(StringIO(newick), "newick")
    
    summary = get_tree_summary(tree)
    print(json.dumps({"test": "basic_tree", "result": summary}))

def test_tree_with_branch_lengths():
    """Test tree with branch lengths."""
    newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
    tree = Phylo.read(StringIO(newick), "newick")
    
    summary = get_tree_summary(tree)
    
    # Get depths
    depths = tree.depths()
    depth_info = {str(clade.name): d for clade, d in depths.items() if clade.name}
    
    # Get distances
    dist_a_b = tree.distance("A", "B")
    dist_a_c = tree.distance("A", "C")
    
    print(json.dumps({
        "test": "tree_with_branch_lengths",
        "result": {
            **summary,
            "depths": depth_info,
            "distance_A_B": dist_a_b,
            "distance_A_C": dist_a_c
        }
    }))

def test_common_ancestor():
    """Test common ancestor finding."""
    newick = "(A:0.5,(B:0.3,C:0.3):0.5);"
    tree = Phylo.read(StringIO(newick), "newick")
    
    # Find MRCA of B and C
    mrca_bc = tree.common_ancestor({"name": "B"}, {"name": "C"})
    mrca_bc_name = mrca_bc.name if mrca_bc.name else "None"
    
    # Find MRCA of A, B, and C
    mrca_abc = tree.common_ancestor({"name": "A"}, {"name": "B"}, {"name": "C"})
    mrca_abc_name = mrca_abc.name if mrca_abc.name else "None"
    
    print(json.dumps({
        "test": "common_ancestor",
        "result": {
            "mrca_B_C": mrca_bc_name,
            "mrca_A_B_C": mrca_abc_name
        }
    }))

def test_find_clades():
    """Test clade finding."""
    newick = "(A,B,(C,D));"
    tree = Phylo.read(StringIO(newick), "newick")
    
    # Find by name
    found_a = list(tree.find_clades(name="A"))
    found_c = list(tree.find_clades(name="C"))
    
    # Find terminals
    terminals = list(tree.find_clades(terminal=True))
    
    # Find any
    any_a = tree.find_any(name="A")
    any_z = tree.find_any(name="Z")
    
    print(json.dumps({
        "test": "find_clades",
        "result": {
            "found_A_count": len(found_a),
            "found_C_count": len(found_c),
            "terminal_count": len(terminals),
            "any_A_exists": any_a is not None,
            "any_Z_exists": any_z is not None
        }
    }))

def test_newick_roundtrip():
    """Test Newick parsing and writing roundtrip."""
    original = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
    tree = Phylo.read(StringIO(original), "newick")
    
    # Write back to Newick
    output = StringIO()
    Phylo.write(tree, output, "newick")
    written = output.getvalue().strip()
    
    # Parse again
    tree2 = Phylo.read(StringIO(written), "newick")
    summary2 = get_tree_summary(tree2)
    
    print(json.dumps({
        "test": "newick_roundtrip",
        "result": {
            "original": original,
            "written": written,
            "roundtrip_terminal_count": summary2["num_terminals"]
        }
    }))

def test_draw_ascii():
    """Test ASCII drawing."""
    newick = "(A,B,(C,D));"
    tree = Phylo.read(StringIO(newick), "newick")
    
    output = StringIO()
    Phylo.draw_ascii(tree, file=output)
    ascii_art = output.getvalue()
    
    contains_all = all(name in ascii_art for name in ["A", "B", "C", "D"])
    
    print(json.dumps({
        "test": "draw_ascii",
        "result": {
            "contains_all_names": contains_all,
            "ascii_length": len(ascii_art)
        }
    }))

if __name__ == "__main__":
    test_basic_tree()
    test_tree_with_branch_lengths()
    test_common_ancestor()
    test_find_clades()
    test_newick_roundtrip()
    test_draw_ascii()