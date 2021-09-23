#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-07-21
Purpose: Compute GC content
"""

import argparse
import sys
from typing import NamedTuple, TextIO, List, Tuple
from Bio import SeqIO


class Args(NamedTuple):
    """ Command-line arguments """
    file: TextIO


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Compute GC content',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
                        help='Input sequence file',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        default=sys.stdin,
                        nargs='?')

    args = parser.parse_args()

    return Args(args.file)


# --------------------------------------------------
def get_gc(seq: str) -> float:
    """ Calculate GC content"""

    return (seq.count('G') + seq.count('C')) / len(seq)


# --------------------------------------------------
def test_get_gc() -> None:
    """ Test get_gc()"""

    assert get_gc('ATGC') == 0.5


# --------------------------------------------------
def main() -> None:
    """ Do the stuff """

    args = get_args()

    recs: List[Tuple[float, str]] = []

    for rec in SeqIO.parse(args.file, 'fasta'):
        recs.append((get_gc(rec.seq) * 100, rec.id))

    high = max(recs)

    print(f'{high[1]} {high[0]:.6f}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
