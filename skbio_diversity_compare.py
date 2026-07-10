#!/usr/bin/env python3
"""
Comparison script: scikit-bio Diversity vs MoonBit implementation.
Compares alpha diversity, beta diversity, and phylogenetic diversity metrics.
"""

import numpy as np
from skbio.diversity import alpha_diversity, beta_diversity
from skbio import TreeNode
import time

# Test data
counts1 = np.array([1, 2, 3, 4, 5, 0, 0, 1, 2, 3])
counts2 = np.array([3, 1, 4, 1, 5, 9, 2, 6, 5, 3])
counts_equal = np.array([1, 1, 1, 1, 1, 1, 1, 1])
counts_single = np.array([10, 0, 0, 0, 0])

print("=" * 60)
print("Alpha Diversity Comparison: scikit-bio vs MoonBit")
print("=" * 60)

# Alpha diversity metrics available in scikit-bio
alpha_metrics = [
    'shannon',
    'simpson',
    'observed_otus',
    'chao1',
    'ace',
    'pielou_e',
    'berger_parker_d',
    'brillouin_d',
    'dominance',
    'goods_coverage',
    'margalef',
    'menhinick',
    'mcintosh_d',
    'mcintosh_e',
    'singles',
    'doubles',
]

print("\n--- Test data 1: mixed counts ---")
print(f"Counts: {counts1}")

for metric in alpha_metrics:
    try:
        result = alpha_diversity(metric, counts1).iloc[0]
        print(f"  {metric:25s}: {result:.6f}")
    except Exception as e:
        print(f"  {metric:25s}: ERROR - {e}")

print("\n--- Test data 2: equal counts (8 species) ---")
print(f"Counts: {counts_equal}")

for metric in alpha_metrics:
    try:
        result = alpha_diversity(metric, counts_equal).iloc[0]
        print(f"  {metric:25s}: {result:.6f}")
    except Exception as e:
        print(f"  {metric:25s}: ERROR - {e}")

print("\n" + "=" * 60)
print("Beta Diversity Comparison")
print("=" * 60)

beta_metrics = [
    'braycurtis',
    'jaccard',
    'euclidean',
    'manhattan',
    'chebyshev',
    'canberra',
    'cosine',
]

print(f"\nSample 1: {counts1}")
print(f"Sample 2: {counts2}")

for metric in beta_metrics:
    try:
        result = beta_diversity(metric, [counts1, counts2])
        dist = result[0, 1]
        print(f"  {metric:20s}: {dist:.6f}")
    except Exception as e:
        print(f"  {metric:20s}: ERROR - {e}")

print("\n" + "=" * 60)
print("Phylogenetic Diversity (Faith's PD)")
print("=" * 60)

# Create a simple tree
newick_tree = "((A:1.0,B:1.0):2.0,C:3.0):0.0;"
tree = TreeNode.read([newick_tree])
taxa = ['A', 'B', 'C']
counts_pd = np.array([1, 1, 1])

print(f"\nTree: {newick_tree}")
print(f"Taxa: {taxa}")
print(f"Counts: {counts_pd}")

try:
    from skbio.diversity import alpha_diversity
    pd_result = alpha_diversity('faith_pd', counts_pd, taxa=taxa, tree=tree).iloc[0]
    print(f"  Faith's PD (all taxa): {pd_result:.6f}")
except Exception as e:
    print(f"  Faith's PD: ERROR - {e}")

# Subset
counts_pd_subset = np.array([1, 0, 1])
try:
    pd_result = alpha_diversity('faith_pd', counts_pd_subset, taxa=taxa, tree=tree).iloc[0]
    print(f"  Faith's PD (A and C):  {pd_result:.6f}")
except Exception as e:
    print(f"  Faith's PD: ERROR - {e}")

print("\n" + "=" * 60)
print("Performance Comparison")
print("=" * 60)

# Large dataset for performance test
np.random.seed(42)
large_counts = np.random.randint(0, 100, size=1000)
print(f"\nLarge dataset: {len(large_counts)} features")

# Performance: alpha diversity
metric = 'shannon'
n_runs = 1000
start = time.time()
for _ in range(n_runs):
    alpha_diversity(metric, large_counts)
elapsed = time.time() - start
print(f"\n  scikit-bio {metric}: {elapsed/n_runs*1e6:.2f} us/call ({n_runs} runs)")

metric = 'simpson'
start = time.time()
for _ in range(n_runs):
    alpha_diversity(metric, large_counts)
elapsed = time.time() - start
print(f"  scikit-bio {metric}: {elapsed/n_runs*1e6:.2f} us/call ({n_runs} runs)")

metric = 'chao1'
start = time.time()
for _ in range(n_runs):
    alpha_diversity(metric, large_counts)
elapsed = time.time() - start
print(f"  scikit-bio {metric}: {elapsed/n_runs*1e6:.2f} us/call ({n_runs} runs)")

# Beta diversity performance
counts_a = np.random.randint(0, 100, size=1000)
counts_b = np.random.randint(0, 100, size=1000)
metric = 'braycurtis'
start = time.time()
for _ in range(n_runs):
    beta_diversity(metric, [counts_a, counts_b])
elapsed = time.time() - start
print(f"\n  scikit-bio {metric}: {elapsed/n_runs*1e6:.2f} us/call ({n_runs} runs)")

print("\n" + "=" * 60)
print("Summary of Expected Values (for MoonBit verification)")
print("=" * 60)

print("\nAlpha diversity (counts = [1, 2, 3, 4, 5, 0, 0, 1, 2, 3]):")
for metric in ['shannon', 'simpson', 'observed_otus', 'chao1', 'berger_parker_d']:
    try:
        result = alpha_diversity(metric, counts1).iloc[0]
        print(f"  {metric:25s} = {result:.10f}")
    except Exception as e:
        print(f"  {metric:25s} = ERROR: {e}")

print("\nBeta diversity (sample1 vs sample2):")
for metric in ['braycurtis', 'jaccard', 'euclidean', 'manhattan']:
    try:
        result = beta_diversity(metric, [counts1, counts2])
        print(f"  {metric:20s} = {result[0,1]:.10f}")
    except:
        pass
