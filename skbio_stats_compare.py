#!/usr/bin/env python3
"""
Comparison script: scikit-bio Stats vs MoonBit implementation.
Compares p-value adjustment, permutation tests, ANOSIM, PERMANOVA, etc.
"""

import numpy as np
from skbio import DistanceMatrix
from skbio.stats.distance import anosim, permanova
from skbio.stats import subsample_counts
from scipy.stats import false_discovery_control
import time

print("=" * 60)
print("Stats Comparison: scikit-bio / scipy vs MoonBit")
print("=" * 60)

# ============================================================
# 1. Subsample counts
# ============================================================
print("\n" + "=" * 60)
print("1. Subsample Counts")
print("=" * 60)

counts = np.array([10, 20, 30, 40], dtype=float)
n = 50
seed = 42

print(f"\nCounts: {counts}")
print(f"Subsample n={n}, with replacement=True, seed={seed}")

# Run with replacement
np.random.seed(seed)
result_with = subsample_counts(counts, n, replace=True)
print(f"  With replacement: {result_with}, sum={result_with.sum()}")

# Run without replacement
np.random.seed(seed)
result_without = subsample_counts(counts, n, replace=False)
print(f"  Without replacement: {result_without}, sum={result_without.sum()}")

# ============================================================
# 2. P-value adjustment
# ============================================================
print("\n" + "=" * 60)
print("2. P-value Adjustment")
print("=" * 60)

pvalues = np.array([0.01, 0.02, 0.03, 0.5, 0.005, 0.001])
print(f"\nOriginal p-values: {pvalues}")

# Bonferroni
bonf = np.minimum(pvalues * len(pvalues), 1.0)
print(f"\nBonferroni:")
print(f"  {bonf}")

# Holm-Bonferroni
def holm_adjust(p):
    n = len(p)
    indexed = sorted(enumerate(p), key=lambda x: x[1])
    adjusted = np.zeros(n)
    prev = 0.0
    for i, (idx, val) in enumerate(indexed):
        adj = val * (n - i)
        if adj < prev:
            adj = prev
        if adj > 1.0:
            adj = 1.0
        adjusted[idx] = adj
        prev = adj
    return adjusted

holm = holm_adjust(pvalues)
print(f"\nHolm-Bonferroni:")
print(f"  {holm}")

# Benjamini-Hochberg (FDR)
bh = false_discovery_control(pvalues, method='bh')
print(f"\nBenjamini-Hochberg (FDR):")
print(f"  {bh}")

# ============================================================
# 3. Permutation test (difference of means)
# ============================================================
print("\n" + "=" * 60)
print("3. Permutation Test (difference of means)")
print("=" * 60)

np.random.seed(seed)
sample1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
sample2 = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
n_permutations = 1000

print(f"\nSample 1: {sample1}")
print(f"Sample 2: {sample2}")
print(f"Permutations: {n_permutations}, seed: {seed}")

obs_diff = np.mean(sample1) - np.mean(sample2)
print(f"Observed difference of means: {obs_diff:.6f}")

# Manual permutation test
combined = np.concatenate([sample1, sample2])
n1 = len(sample1)
n_total = len(combined)

np.random.seed(seed)
count_ge = 0
for _ in range(n_permutations):
    perm = np.random.permutation(combined)
    perm_diff = np.mean(perm[:n1]) - np.mean(perm[n1:])
    if perm_diff >= obs_diff:
        count_ge += 1

p_value = count_ge / n_permutations
print(f"Permutation p-value (one-sided, greater): {p_value:.6f}")

# ============================================================
# 4. ANOSIM
# ============================================================
print("\n" + "=" * 60)
print("4. ANOSIM (Analysis of Similarities)")
print("=" * 60)

# Create a distance matrix with clear group structure
data = np.array([
    [0.0, 1.0, 2.0, 8.0, 9.0, 10.0],
    [1.0, 0.0, 3.0, 7.0, 8.0, 9.0],
    [2.0, 3.0, 0.0, 6.0, 7.0, 8.0],
    [8.0, 7.0, 6.0, 0.0, 1.0, 2.0],
    [9.0, 8.0, 7.0, 1.0, 0.0, 3.0],
    [10.0, 9.0, 8.0, 2.0, 3.0, 0.0],
])
ids = ['A', 'B', 'C', 'D', 'E', 'F']
groups = [0, 0, 0, 1, 1, 1]

dm = DistanceMatrix(data, ids)
print(f"\nDistance matrix ({len(ids)} samples)")
print(f"Groups: {groups}")
print(f"Permutations: 999, seed: {seed}")

result_anosim = anosim(dm, groups, permutations=999, seed=seed)
print(f"\nANOSIM result:")
print(f"  R statistic: {result_anosim['test statistic']:.6f}")
print(f"  p-value: {result_anosim['p-value']:.6f}")
print(f"  Number of groups: {result_anosim['number of groups']}")
print(f"  Sample size: {result_anosim['sample size']}")
print(f"  Permutations: {result_anosim['number of permutations']}")

# ============================================================
# 5. PERMANOVA
# ============================================================
print("\n" + "=" * 60)
print("5. PERMANOVA (Permutational MANOVA)")
print("=" * 60)

result_permanova = permanova(dm, groups, permutations=999, seed=seed)
print(f"\nPERMANOVA result:")
print(f"  F statistic: {result_permanova['test statistic']:.6f}")
print(f"  p-value: {result_permanova['p-value']:.6f}")
print(f"  Number of groups: {result_permanova['number of groups']}")
print(f"  Sample size: {result_permanova['sample size']}")
print(f"  Permutations: {result_permanova['number of permutations']}")

# ============================================================
# 6. Confidence interval (normal approximation)
# ============================================================
print("\n" + "=" * 60)
print("6. Confidence Interval (Normal Approximation)")
print("=" * 60)

from scipy.stats import norm

mean_val = 5.0
std_val = 2.0
n_samples = 100
confidence = 0.95

z = norm.ppf(0.5 + confidence / 2)
se = std_val / np.sqrt(n_samples)
lower = mean_val - z * se
upper = mean_val + z * se

print(f"\nMean: {mean_val}, Std: {std_val}, n: {n_samples}, confidence: {confidence}")
print(f"Z-score: {z:.6f}")
print(f"Standard error: {se:.6f}")
print(f"Confidence interval: ({lower:.6f}, {upper:.6f})")

print("\n" + "=" * 60)
print("Comparison summary - use these values to verify MoonBit")
print("=" * 60)
