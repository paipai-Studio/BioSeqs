#!/usr/bin/env python3
"""
Comparison script to verify MoonBit Tree implementation against scikit-bio.
"""

try:
    from skbio import TreeNode
except ImportError:
    print("Error: scikit-bio is not installed. Install with: pip install scikit-bio")
    exit(1)

def test_simple_tree():
    """Test simple tree parsing and traversal."""
    newick = "((A,B)C,D)E;"
    
    tree = TreeNode.read([newick])
    
    print("=== Test: Simple Tree ===")
    print(f"Newick: {newick}")
    print(f"Root name: {tree.name}")
    print(f"Number of children: {len(tree.children)}")
    print(f"Number of tips: {len(list(tree.tips()))}")
    print(f"Total nodes: {len(list(tree.preorder()))}")
    print()
    
    print("Preorder traversal:")
    for node in tree.preorder():
        print(f"  {node.name}")
    print()
    
    print("Postorder traversal:")
    for node in tree.postorder():
        print(f"  {node.name}")
    print()
    
    print("Levelorder traversal:")
    for node in tree.levelorder():
        print(f"  {node.name}")
    print()
    
    print("Tips:")
    for tip in tree.tips():
        print(f"  {tip.name}")
    print()

def test_tree_with_lengths():
    """Test tree with branch lengths."""
    newick = "(A:0.1,B:0.2)root:0.3;"
    
    tree = TreeNode.read([newick])
    
    print("=== Test: Tree with Lengths ===")
    print(f"Newick: {newick}")
    print(f"Total length: {tree.total_length()}")
    print()
    
    print("Node details:")
    for node in tree.preorder():
        print(f"  {node.name}: length={node.length}")
    print()

def test_lca():
    """Test LCA calculation."""
    newick = "((A,B)C,D)E;"
    
    tree = TreeNode.read([newick])
    
    print("=== Test: LCA ===")
    print(f"Newick: {newick}")
    
    lca_ab = tree.lca(["A", "B"])
    print(f"LCA(A, B): {lca_ab.name}")
    
    lca_ac = tree.lca(["A", "C"])
    print(f"LCA(A, C): {lca_ac.name}")
    
    lca_ad = tree.lca(["A", "D"])
    print(f"LCA(A, D): {lca_ad.name}")
    print()

def test_distance():
    """Test distance calculation."""
    newick = "(A:0.1,B:0.2)C:0.3;"
    
    tree = TreeNode.read([newick])
    
    node_a = tree.find("A")
    node_b = tree.find("B")
    node_c = tree.find("C")
    
    print("=== Test: Distance ===")
    print(f"Newick: {newick}")
    print(f"Distance(A, B): {node_a.distance(node_b)}")
    print(f"Distance(A, C): {node_a.distance(node_c)}")
    print()

def test_node_finding():
    """Test node finding."""
    newick = "((A,B)C,D)E;"
    
    tree = TreeNode.read([newick])
    
    print("=== Test: Node Finding ===")
    print(f"Newick: {newick}")
    
    found_c = tree.find("C")
    print(f"Found C: {found_c.name}")
    
    try:
        found_x = tree.find("X")
        print(f"Found X: {found_x}")
    except Exception as e:
        print(f"Found X: None (exception: {type(e).__name__})")
    print()

def test_newick_roundtrip():
    """Test Newick parsing and output roundtrip."""
    original = "((A:0.1,B:0.2)C:0.3,D:0.4)E:0.5;"
    
    tree = TreeNode.read([original])
    output = str(tree)
    
    print("=== Test: Newick Roundtrip ===")
    print(f"Original: {original}")
    print(f"Output: {output}")
    
    # Parse output again
    tree2 = TreeNode.read([output])
    output2 = str(tree2)
    print(f"Output2: {output2}")
    print()

if __name__ == "__main__":
    print("=" * 50)
    print("scikit-bio Tree Module Comparison Tests")
    print("=" * 50)
    print()
    
    test_simple_tree()
    test_tree_with_lengths()
    test_lca()
    test_distance()
    test_node_finding()
    test_newick_roundtrip()
    
    print("=" * 50)
    print("All tests completed!")
    print("=" * 50)
