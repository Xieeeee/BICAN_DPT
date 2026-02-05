#!/usr/bin/env python3

import sys
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert values line-by-line after the first two header lines."
    )
    p.add_argument("input", nargs="?", default="-", help="input file (default: stdin)")
    p.add_argument("output", nargs="?", default="-", help="output file (default: stdout)")
    p.add_argument("--cutoff", type=float, default=0.5,
                   help="threshold in (0,1); map [0, cutoff) -> 0 and [cutoff, 1] -> 1 (default: 0.5)")
    return p.parse_args()

def convert_value(s: str, cutoff: float) -> str:
    s = s.strip()
    if s == "-1.0":
        return "2"
    try:
        v = float(s)
    except ValueError:
        return s  # keep non-numeric as-is

    # Only classify within [0,1]; pass through other numbers unchanged
    if 0.0 <= v < cutoff:
        return "0"
    if cutoff <= v <= 1.0:
        return "1"
    return s

def main():
    args = parse_args()

    # validate cutoff
    if not (0.0 <= args.cutoff <= 1.0):
        sys.exit(f"[error] --cutoff must be between 0 and 1 (exclusive). Got: {args.cutoff}")

    fin = sys.stdin if args.input == "-" else open(args.input, "r")
    fout = sys.stdout if args.output == "-" else open(args.output, "w")

    try:
        for i, line in enumerate(fin, 1):
            if i <= 2:
                fout.write(line if line.endswith("\n") else line + "\n")
            else:
                fout.write(convert_value(line, args.cutoff) + "\n")
    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

if __name__ == "__main__":
    main()

