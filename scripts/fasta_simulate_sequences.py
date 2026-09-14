#!/usr/bin/env python
"""Print random substrings of a reference FASTA as FASTA records on stdout.

Each output record is a run of ``--out-length`` bases taken from a random
position of a randomly chosen reference sequence; ``--num-seqs`` records are
produced. Gzipped input is accepted. A test-data helper, not part of the
pipeline.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from random import randint


def read_fasta(path: str) -> dict[str, str]:
    """Return record name -> sequence for a plain or gzipped FASTA file."""
    opener = gzip.open if path.endswith(".gz") else open
    sequences: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with opener(path, "rt", encoding="utf-8") as fo:
        for line in fo:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts)
                name, parts = line[1:], []
            elif name is None:
                sys.exit("not a FASTA file: the first line does not start with '>'")
            else:
                parts.append(line)
    if name is not None:
        sequences[name] = "".join(parts)
    return sequences


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Produce random substrings (FASTA) of an input FASTA on stdout."
    )
    parser.add_argument("input", help="Reference FASTA (optionally gzipped).")
    parser.add_argument(
        "-l", "--out-length", type=int, default=1000, help="Length of each output sequence."
    )
    parser.add_argument(
        "-n", "--num-seqs", type=int, default=1000, help="Number of output sequences."
    )
    parser.add_argument(
        "-b", "--basename", help="Record name prefix (read -> read_0, read_1, ...)."
    )
    args = parser.parse_args()

    references = read_fasta(args.input)
    names = list(references)
    written = 0
    while written < args.num_seqs:
        ref_name = names[randint(0, len(names) - 1)]
        ref_seq = references[ref_name]
        start = randint(0, len(ref_seq))
        seq = ref_seq[start : start + args.out_length]
        if not seq:
            continue
        header = args.basename or ref_name.split()[0]
        print(f">{header}_{written}\n{seq}")
        written += 1


if __name__ == "__main__":
    main()
