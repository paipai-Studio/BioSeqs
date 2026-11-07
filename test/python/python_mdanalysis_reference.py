#!/usr/bin/env python3
"""Generate reference MD-analysis results from MDAnalysis as JSON.

The same hard-coded coordinate system is evaluated by the MoonBit companion
program (cmd/md_compare/main.mbt) and the two JSON outputs are diffed by
mdanalysis_compare.sh with a small tolerance.

Requires MDAnalysis (pip install MDAnalysis). Exits with status 2 when the
package is missing so callers can skip the comparison gracefully.
"""
import json
import sys

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    from MDAnalysis.lib import mdamath
    from MDAnalysis.lib.distances import calc_dihedrals, distance_array
except ImportError:
    print("MDAnalysis is not installed (pip install MDAnalysis)", file=sys.stderr)
    sys.exit(2)


def run_cases():
    cases = []

    def record(name, values):
        cases.append({
            "name": name,
            "result": ",".join("%.6f" % float(v) for v in values),
        })

    # Hard-coded system: 3 atoms (N, C, O), 2 frames — must match
    # cmd/md_compare/main.mbt exactly.
    x0 = np.array(
        [[0.5, 0.2, 0.1], [1.5, 0.7, -0.3], [2.2, -0.4, 0.6]], dtype=np.float64
    )
    x1 = np.array(
        [[0.1, -0.3, 0.9], [1.8, 0.5, 0.4], [2.5, 0.9, -0.7]], dtype=np.float64
    )
    masses = np.array([14.007, 12.011, 15.999], dtype=np.float64)

    # 1. pairwise distance matrix (frame 0)
    dm = distance_array(x0, x0)
    record("distance_matrix", dm.ravel())

    # 2. RMSD of frame 1 onto frame 0 (centered + optimal superposition)
    r = rms.rmsd(x0, x1, center=True, superposition=True)
    record("rmsd_sup", [r])

    # 3. per-atom RMSF over the two frames
    f = rms.rmsf(np.array([x0, x1], dtype=np.float64))
    record("rmsf", f)

    # 4. mass-weighted radius of gyration (frame 0)
    record("radius_of_gyration", [mdamath.radius_of_gyration(x0, masses)])

    # 5. geometric center (frame 0)
    record("center_of_geometry", x0.mean(axis=0))

    # 6. center of mass (frame 0)
    record("center_of_mass", mdamath.weighted_center(x0, weights=masses))

    # 7. dihedral angle of four points (radians, absolute value because the
    #    sign conventions of the two implementations may differ)
    p1 = np.array([[0.0, 0.0, 0.0]])
    p2 = np.array([[1.0, 0.0, 0.0]])
    p3 = np.array([[1.0, 1.0, 0.0]])
    p4 = np.array([[1.0, 1.0, 1.0]])
    d = calc_dihedrals(p1, p2, p3, p4)[0]
    record("dihedral_abs", [abs(d)])

    # 8. single bond distance
    record("bond_distance", [distance_array(x0[:1], x0[1:2])[0][0]])

    return cases


def main():
    cases = run_cases()
    json.dump(cases, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
