#!/usr/bin/env python3
"""
generate_pmp.py — Command-line tool to generate SDPB PMP input files.

Usage
-----
    python3 generate_pmp.py [options]

This script generates a Polynomial Matrix Program (PMP) JSON file that can be
fed directly into SDPB (via pmp2sdp + sdpb) to obtain bounds on the ratio
of two Wilson coefficients W_{obj} / W_{norm}.

Quick start (example from problem statement):
    python3 generate_pmp.py \\
        --obj 1 0 --norm 2 0 \\
        --d 4 --K 4 --max-spin 10 --delta0 40 \\
        --direction upper \\
        --output pmp_upper.json

This bounds W_{1,0} / W_{2,0} from above (d=4, K=4, delta0=40, max spin 10).

Wilson coefficients W_{n-m, m}:
  - Objectives (n >= m >= 0): physically measurable operators.
    * (n=1,m=0): W_{1,0}  — corresponds to g_2 in the 2-to-2 amplitude.
    * (n=1,m=1): W_{0,1}  — corresponds to g_3.
    * (n=2,m=0): W_{2,0}  — corresponds to g_4.
    * (n=2,m=1): W_{1,1}  — higher-order.
  - Null constraints (m > n >= 1): enforced automatically.

For a 2D allowed region plot of (g_3/g_2, g_4/g_2):
  Run this script multiple times to get upper and lower bounds on:
    - W_{1,1}/W_{1,0}  (= g_3/g_2 type ratio)
    - W_{2,0}/W_{1,0}  (= g_4/g_2 type ratio)
  See SDPB_GUIDE.md for the full workflow.
"""

import argparse
import json
import os
import sys

# Allow running from the repo root or from within eft_bounds/
_this_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(_this_dir)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from eft_bounds.pmp_generator import generate_csdr_pmp, write_pmp_json


def main():
    parser = argparse.ArgumentParser(
        description="Generate SDPB PMP JSON file for EFT bounds via CSDR dispersion relations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Upper bound on W_{1,0}/W_{2,0} in d=4, K=4, delta0=40, max_spin=10:
  python3 generate_pmp.py --obj 1 0 --norm 2 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper --output upper_W10_W20.json

  # Lower bound on same ratio:
  python3 generate_pmp.py --obj 1 0 --norm 2 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction lower --output lower_W10_W20.json

  # W_{0,1}/W_{1,0} ratio (involves n=m pair):
  python3 generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper --output upper_W01_W10.json

Decision variable notation:
  W_{n-m, m}  indexed by (n, m) where n >= m >= 0.
  Examples: W_{1,0} = (n=1,m=0),  W_{0,1} = (n=1,m=1),  W_{2,0} = (n=2,m=0).
        """,
    )

    parser.add_argument(
        "--obj", nargs=2, type=int, metavar=("N", "M"), required=True,
        help="Objective Wilson coefficient (n, m) to bound. Must have n >= m >= 0."
    )
    parser.add_argument(
        "--norm", nargs=2, type=int, metavar=("N", "M"), required=True,
        help="Normalization Wilson coefficient (n, m), set equal to 1. Must have n >= m >= 0."
    )
    parser.add_argument(
        "--d", type=int, required=True,
        help="Spacetime dimension (e.g. 4)."
    )
    parser.add_argument(
        "--K", type=int, default=8,
        help="Maximum spectral power 2n+m (higher = more operators, slower). Default: 8."
    )
    parser.add_argument(
        "--max-spin", type=int, default=10,
        help="Maximum even spin included in blocks. Default: 10."
    )
    parser.add_argument(
        "--delta0", type=int, default=1,
        help="IR cutoff: integration starts at s1 = delta0. Default: 1."
    )
    parser.add_argument(
        "--direction", choices=["upper", "lower"], default="upper",
        help="Optimize upper or lower bound. Default: upper."
    )
    parser.add_argument(
        "--precision", type=int, default=200,
        help="Decimal precision for output coefficients. Default: 200."
    )
    parser.add_argument(
        "--output", "-o", type=str, default=None,
        help="Output JSON file path. Default: pmp_{direction}_W{obj}_W{norm}.json"
    )
    parser.add_argument(
        "--list-pairs", action="store_true",
        help="List all available objective and null pairs for the given K, then exit."
    )

    args = parser.parse_args()

    # Import enumeration helpers for --list-pairs
    from eft_bounds.csdr import enumerate_obj_pairs_nm, enumerate_null_pairs

    if args.list_pairs:
        obj_pairs = enumerate_obj_pairs_nm(args.K)
        null_pairs = enumerate_null_pairs(args.K)
        print(f"\nObjective pairs (n >= m >= 0, 2n+m <= {args.K}):")
        for n, m in obj_pairs:
            label = "W_{0," + str(m) + "}" if n == m else f"W_{{{n-m},{m}}}"
            print(f"  (n={n}, m={m})  =>  {label}  [spectral power {2*n+m}]")
        print(f"\nNull constraint pairs (m > n >= 1, 2n+m <= {args.K}):")
        for n, m in null_pairs:
            print(f"  (n={n}, m={m})  =>  W_{{{n-m},{m}}}  [spectral power {2*n+m}]")
        return

    obj_index = tuple(args.obj)
    norm_index = tuple(args.norm)

    if args.output is None:
        obj_str = f"{obj_index[0]}{obj_index[1]}"
        norm_str = f"{norm_index[0]}{norm_index[1]}"
        args.output = f"pmp_{args.direction}_W{obj_str}_over_W{norm_str}.json"

    print(f"Generating PMP:")
    print(f"  Bound:     {args.direction} of W_{{{obj_index[0]-obj_index[1]},{obj_index[1]}}} / "
          f"W_{{{norm_index[0]-norm_index[1]},{norm_index[1]}}}")
    print(f"  d={args.d}, K={args.K}, max_spin={args.max_spin}, delta0={args.delta0}")
    print(f"  Output:    {args.output}")

    pmp = generate_csdr_pmp(
        obj_index=obj_index,
        norm_index=norm_index,
        d=args.d,
        K=args.K,
        max_spin=args.max_spin,
        precision=args.precision,
        bound_direction=args.direction,
        delta0=args.delta0,
    )

    write_pmp_json(pmp, args.output)

    meta = pmp["_metadata"]
    print(f"\nDone.")
    print(f"  Decision variables: {meta['n_decision_variables']} "
          f"({meta['n_obj_pairs']} objectives, {meta['n_null_pairs']} null constraints)")
    print(f"  Spin blocks:        {meta['n_spin_blocks']}")
    print(f"\nNext steps:")
    print(f"  1. Convert to SDPB binary: pmp2sdp --input {args.output} --output sdp/")
    print(f"  2. Run solver:             sdpb --sdpDir sdp/ --outDir out/")
    print(f"  3. Read result from:       out/out.txt  (look for 'primalObjective')")


if __name__ == "__main__":
    main()
