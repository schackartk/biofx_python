#!/usr/bin/env python3
"""
Author : Kenneth Schackart <schackartk1@gmail.com>
Date   : 2021-07-21
Purpose: Compute GC content
"""

import argparse
import sys
from typing import NamedTuple, TextIO, List
from Bio import SeqIO


class Args(NamedTuple):
    """ Command-line arguments """
    file: TextIO


# --------------------------------------------------
class MySeq(NamedTuple):
    """ Sequence """
    gc: float
    name: str


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

    assert get_gc('C') == 1.
    assert get_gc('G') == 1.
    assert get_gc('CGCGCG') == 1.
    assert get_gc('ATATAT') == 0.0
    assert get_gc('ATGC') == 0.5


# --------------------------------------------------
def main() -> None:
    """ Do the stuff """

    args = get_args()

    recs: List[MySeq] = []

    for rec in SeqIO.parse(args.file, 'fasta'):
        recs.append(MySeq(get_gc(rec.seq) * 100, rec.id))

    high = max(recs)

    print(f'{high.name} {high.gc:.6f}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
